// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <cmath>
#include <cstddef>
#include <optional>
#include <vector>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tags/TempTensor.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/Slice.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/AreaElement.hpp"
#include "Domain/ExcisionSphere.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/IndexToSliceAt.hpp"
#include "Domain/Tags.hpp"
#include "Domain/TagsTimeDependent.hpp"
#include "Evolution/Systems/CurvedScalarWave/System.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/ElementActions/ReceiveWorldtubeData.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Inboxes.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/SingletonChare.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/DefiniteIntegral.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/RealSphericalHarmonics.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Spherepack.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Tags {
struct TimeStepId;
}  // namespace Tags
/// \endcond

namespace CurvedScalarWave::Worldtube::Actions {
/*!
 * \brief Projects the regular field \f$\Psi^R\f$ and its time derivative
 * \f$\partial_t \Psi^R\f$ onto real spherical harmonics and sends the result to
 * the worldtube.
 *
 * \details The regular field is obtained by subtracting the singular/puncture
 * field from the numerical DG field.
 * All spherical harmonics are computed for \f$l <= n\f$, where \f$n\f$ is the
 * worldtube expansion order. The projection is done by integrating over the DG
 * grid of the element face using \ref definite_integral with the euclidean area
 * element. The worldtube adds up all integrals from the different elements to
 * obtain the integral over the entire sphere.
 *
 * The projection is done in the co-moving grid frame in which the worldtube
 * is at rest: the time derivative of the regular field is transformed to the
 * grid frame with an advective term coming from the mesh velocity, and the
 * spherical harmonics are evaluated at the grid frame coordinates of the
 * face. The Euclidean area element from `Tags::FaceQuantities` is computed in
 * the inertial frame, which agrees with the grid frame one because the
 * grid-to-inertial map of this scheme is a rigid rotation on the worldtube
 * boundary.
 *
 * DataBox:
 * - Uses:
 *    - `tags_to_slice_on_face`
 *    - `Worldtube::Tags::ExpansionOrder`
 *    - `Worldtube::Tags::FaceCoordinates<Dim, Frame::Grid, true>`
 *    - `Worldtube::Tags::GeodesicPunctureField`
 *    - `Worldtube::Tags::ExcisionSphere`
 *    - `Tags::TimeStepId`
 */
struct SendToWorldtube {
  static constexpr size_t Dim = 3;
  using tags_to_send = tmpl::list<CurvedScalarWave::Tags::Psi,
                                  ::Tags::dt<CurvedScalarWave::Tags::Psi>>;
  using tags_to_slice_to_face =
      tmpl::list<CurvedScalarWave::Tags::Psi, CurvedScalarWave::Tags::Pi,
                 CurvedScalarWave::Tags::Phi<Dim>,
                 gr::Tags::Shift<DataVector, Dim>, gr::Tags::Lapse<DataVector>,
                 domain::Tags::InverseJacobian<Dim, Frame::ElementLogical,
                                               Frame::Inertial>>;
  using simple_tags = tmpl::list<Tags::RegularFieldAdvectiveTerm<Dim>>;

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& box,
      tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    const auto& element_id = db::get<domain::Tags::Element<Dim>>(box).id();
    const auto& excision_sphere = db::get<Tags::ExcisionSphere<Dim>>(box);
    const auto direction = excision_sphere.abutting_direction(element_id);
    if (not direction.has_value()) {
      return {Parallel::AlgorithmExecution::Continue, std::nullopt};
    }
    const auto& mesh = db::get<domain::Tags::Mesh<Dim>>(box);
    const auto face_mesh = mesh.slice_away(direction->dimension());
    const size_t face_size = face_mesh.number_of_grid_points();

    Variables<tmpl::list<
        CurvedScalarWave::Tags::Psi, ::Tags::dt<CurvedScalarWave::Tags::Psi>,
        ::Tags::TempScalar<0>, ::Tags::TempScalar<1>, ::Tags::TempScalar<2>,
        ::Tags::Tempi<3, Dim>, ::Tags::TempI<4, Dim>>>
        temporaries(face_size);
    auto& psi_regular_times_det =
        get(get<CurvedScalarWave::Tags::Psi>(temporaries));
    auto& dt_psi_regular_times_det =
        get(get<::Tags::dt<CurvedScalarWave::Tags::Psi>>(temporaries));
    auto& theta = get(get<::Tags::TempScalar<0>>(temporaries));
    auto& phi = get(get<::Tags::TempScalar<1>>(temporaries));
    auto& spherical_harmonic = get(get<::Tags::TempScalar<2>>(temporaries));
    auto& face_phi = get<::Tags::Tempi<3, Dim>>(temporaries);
    auto& face_mesh_velocity = get<::Tags::TempI<4, Dim>>(temporaries);
    const auto& face_quantities = db::get<Tags::FaceQuantities>(box).value();
    const auto& psi_numerical_face =
        get<CurvedScalarWave::Tags::Psi>(face_quantities);
    const auto& dt_psi_numerical_face =
        get<::Tags::dt<CurvedScalarWave::Tags::Psi>>(face_quantities);
    const auto& area_element =
        get<gr::surfaces::Tags::AreaElement<DataVector>>(face_quantities);
    const auto& puncture_field =
        db::get<Tags::CurrentIteration>(box) > 0
            ? db::get<Tags::IteratedPunctureField<Dim>>(box).value()
            : db::get<Tags::GeodesicPunctureField<Dim>>(box).value();
    const auto& psi_puncture = get<CurvedScalarWave::Tags::Psi>(puncture_field);
    const auto& dt_psi_puncture =
        get<::Tags::dt<CurvedScalarWave::Tags::Psi>>(puncture_field);

    // A spherical-harmonic shell block covers the entire worldtube boundary
    // with a single element. In that case the projection integrals are
    // evaluated with the Gauss quadrature of the spherical-harmonic
    // collocation grid, which is spectrally exact. The Gauss weights include
    // the sin(theta) of the area element, so the euclidean area element is
    // replaced by the quadrature weights times the squared worldtube radius.
    const bool spherical_harmonic_face =
        face_mesh.basis(0) == Spectral::Basis::SphericalHarmonic;
    DataVector s2_integration_weights{};
    if (spherical_harmonic_face) {
      const size_t n_theta = face_mesh.extents(0);
      const size_t n_phi = face_mesh.extents(1);
      std::vector<double> gauss_points(n_theta + 1);
      std::vector<double> gauss_weights(n_theta + 1);
      std::vector<double> gaqd_work(n_theta);
      int gaqd_err = 0;
      gaqd_(static_cast<int>(n_theta), gauss_points.data(),
            gauss_weights.data(), gaqd_work.data(),
            static_cast<int>(gauss_weights.size()), &gaqd_err);
      if (UNLIKELY(gaqd_err != 0)) {
        ERROR("gaqd error " << gaqd_err << " in SendToWorldtube");
      }
      s2_integration_weights = DataVector(face_size);
      const double phi_weight_times_r_squared =
          2. * M_PI / static_cast<double>(n_phi) *
          square(excision_sphere.radius());
      for (size_t j = 0; j < n_phi; ++j) {
        for (size_t i = 0; i < n_theta; ++i) {
          s2_integration_weights[i + j * n_theta] =
              gauss_weights[i] * phi_weight_times_r_squared;
        }
      }
    }
    const DataVector& integration_weight =
        spherical_harmonic_face ? s2_integration_weights : get(area_element);

    psi_regular_times_det =
        (get(psi_numerical_face) - get(psi_puncture)) * integration_weight;

    const auto& mesh_velocity = db::get<domain::Tags::MeshVelocity<Dim>>(box);
    ASSERT(mesh_velocity.has_value(),
           "Expected a moving grid for worldtube evolution.");
    data_on_slice(make_not_null(&face_phi),
                  db::get<CurvedScalarWave::Tags::Phi<Dim>>(box),
                  mesh.extents(), direction.value().dimension(),
                  index_to_slice_at(mesh.extents(), direction.value()));
    data_on_slice(make_not_null(&face_mesh_velocity), mesh_velocity.value(),
                  mesh.extents(), direction.value().dimension(),
                  index_to_slice_at(mesh.extents(), direction.value()));
    // The advective term transforms the time derivative of the regular field
    // into the co-moving grid frame in which the worldtube is at rest. It is
    // saved to the DataBox because it is used again in `ReceiveWorldtubeData`
    // to transform the time derivative of the regular field sent back by the
    // worldtube to the inertial frame.
    db::mutate<Tags::RegularFieldAdvectiveTerm<Dim>>(
        [&face_phi, &face_mesh_velocity,
         &di_psi_puncture =
             get<::Tags::deriv<CurvedScalarWave::Tags::Psi, tmpl::size_t<3>,
                               Frame::Inertial>>(puncture_field)](
            const gsl::not_null<Scalar<DataVector>*> regular_advective_term) {
          tenex::evaluate<>(regular_advective_term,
                            (face_phi(ti::i) - di_psi_puncture(ti::i)) *
                                face_mesh_velocity(ti::I));
        },
        make_not_null(&box));
    dt_psi_regular_times_det =
        (get(dt_psi_numerical_face) - get(dt_psi_puncture) +
         get(db::get<Tags::RegularFieldAdvectiveTerm<Dim>>(box))) *
        integration_weight;
    const auto& centered_face_coords =
        db::get<Tags::FaceCoordinates<Dim, Frame::Grid, true>>(box);
    ASSERT(centered_face_coords.has_value(),
           "Should be an abutting element here, but face coords are not "
           "calculated!");
    const auto& x = get<0>(centered_face_coords.value());
    const auto& y = get<1>(centered_face_coords.value());
    const auto& z = get<2>(centered_face_coords.value());

    const size_t order = db::get<Worldtube::Tags::ExpansionOrder>(box);
    const size_t num_modes = (order + 1) * (order + 1);
    Variables<tags_to_send> Ylm_coefs(num_modes);
    theta = atan2(hypot(x, y), z);
    phi = atan2(y, x);
    size_t index = 0;
    // project onto spherical harmonics
    for (size_t l = 0; l <= order; ++l) {
      // NOLINTNEXTLINE(bugprone-narrowing-conversions,cppcoreguidelines-narrowing-conversions)
      for (int m = -l; m <= static_cast<int>(l); ++m, ++index) {
        spherical_harmonic = ylm::real_spherical_harmonic(theta, phi, l, m);
        if (spherical_harmonic_face) {
          // the quadrature weights are already included in the integrands
          double psi_coef = 0.;
          double dt_psi_coef = 0.;
          for (size_t k = 0; k < face_size; ++k) {
            psi_coef += psi_regular_times_det[k] * spherical_harmonic[k];
            dt_psi_coef += dt_psi_regular_times_det[k] * spherical_harmonic[k];
          }
          get(get<CurvedScalarWave::Tags::Psi>(Ylm_coefs)).at(index) = psi_coef;
          get(get<::Tags::dt<CurvedScalarWave::Tags::Psi>>(Ylm_coefs))
              .at(index) = dt_psi_coef;
        } else {
          get(get<CurvedScalarWave::Tags::Psi>(Ylm_coefs)).at(index) =
              definite_integral(psi_regular_times_det * spherical_harmonic,
                                face_mesh);
          get(get<::Tags::dt<CurvedScalarWave::Tags::Psi>>(Ylm_coefs))
              .at(index) =
              definite_integral(dt_psi_regular_times_det * spherical_harmonic,
                                face_mesh);
        }
      }
    }
    ASSERT(index == num_modes, "Internal indexing error. "
                                   << num_modes
                                   << " modes should have been calculated but "
                                   << index << " modes were computed.");

    auto& worldtube_component = Parallel::get_parallel_component<
        Worldtube::WorldtubeSingleton<Metavariables>>(cache);
    Parallel::receive_data<Worldtube::Tags::SphericalHarmonicsInbox<Dim>>(
        worldtube_component, db::get<::Tags::TimeStepId>(box),
        std::make_pair(element_id, std::move(Ylm_coefs)));
    if (db::get<Tags::CurrentIteration>(box) + 1 <
        db::get<Tags::MaxIterations>(box) ) {
      db::mutate<Tags::CurrentIteration>(
          [](const gsl::not_null<size_t*> current_iteration) {
            *current_iteration += 1;
          },
          make_not_null(&box));
      // still iterating, go to `IteratePunctureField`
      return {Parallel::AlgorithmExecution::Continue, std::nullopt};

    } else {
      db::mutate<Tags::CurrentIteration>(
          [](const gsl::not_null<size_t*> current_iteration) {
            *current_iteration = 0;
          },
          make_not_null(&box));
      // done iterating, get data for BCs
      return {Parallel::AlgorithmExecution::Continue,
              tmpl::index_of<ActionList, ReceiveWorldtubeData>::value};
    }
  }
};
}  // namespace CurvedScalarWave::Worldtube::Actions
