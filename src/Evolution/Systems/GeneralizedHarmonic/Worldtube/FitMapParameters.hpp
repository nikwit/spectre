// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ConstraintDampingTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Matcher.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Tags.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "IO/Observer/ReductionActions.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Spherepack.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/SpherepackCache.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Time/Tags/Time.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::Worldtube::Actions {
/*!
 * \brief Element-local online worldtube matching (runs in the step loop).
 *
 * On the element owning the excision face — with the spherical-harmonic
 * angular basis there is exactly one, self-identified as the block-0 element
 * with an external lower-radial boundary — this action slices the evolved
 * fields to the excision face, fits the 13 first-order affine-map parameters
 * per `gh::Worldtube::fit_map_parameters`, stores the result in
 * `Tags::MapParameters`, and writes a diagnostic row to
 * `/WorldtubeMatcher.dat` in the reductions file.
 *
 * No-op on every other element, when `Tags::Matcher` is `None`, and when the
 * time has not advanced by `FitInterval` since the last fit (this also
 * limits the action to one fit per step under substepping time steppers).
 * Rates are estimated by backward differences of the fit history and held
 * fixed during the fit (the ghost form needs no better, findings §15b).
 */
struct FitMapParameters {
  using const_global_cache_tags = tmpl::list<Tags::Matcher>;

  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box,
      tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*component*/) {
    static constexpr size_t Dim = 3;
    const auto& config_opt = db::get<Tags::Matcher>(box);
    if (not config_opt.has_value()) {
      return {Parallel::AlgorithmExecution::Continue, std::nullopt};
    }
    const auto& element = db::get<domain::Tags::Element<Dim>>(box);
    const auto inner_face = Direction<Dim>::lower_xi();
    if (element.id().block_id() != 0 or
        element.external_boundaries().count(inner_face) == 0) {
      return {Parallel::AlgorithmExecution::Continue, std::nullopt};
    }
    const double time = db::get<::Tags::Time>(box);
    const auto& state = db::get<Tags::MapParameters>(box);
    if (state.valid and
        time < state.last_fit_time + config_opt->fit_interval - 1.0e-12) {
      return {Parallel::AlgorithmExecution::Continue, std::nullopt};
    }

    const auto& mesh = db::get<domain::Tags::Mesh<Dim>>(box);
    ASSERT(mesh.basis(0) != Spectral::Basis::SphericalHarmonic and
               mesh.basis(1) == Spectral::Basis::SphericalHarmonic and
               mesh.basis(2) == Spectral::Basis::SphericalHarmonic,
           "The online worldtube matcher requires the spherical-harmonic "
           "angular basis with the radial direction as dimension 0, but got "
           "mesh "
               << mesh);
    ASSERT(element.id().segment_id(1).refinement_level() == 0 and
               element.id().segment_id(2).refinement_level() == 0,
           "The online worldtube matcher requires a single element covering "
           "the excision sphere (no angular h-refinement), but got element "
               << element.id());
    const size_t n_radial = mesh.extents(0);
    const size_t n_theta = mesh.extents(1);
    const size_t n_phi = mesh.extents(2);
    const size_t l_max = n_theta - 1;
    ASSERT(n_phi == 2 * l_max + 1,
           "Angular extents inconsistent with a Spherepack grid: "
               << n_theta << " x " << n_phi);
    const size_t n_face = n_theta * n_phi;

    // slice a volume component to the inner radial face (radial index is
    // the fastest-varying, matching the Ylm filter's storage convention)
    const auto face_slice = [n_radial, n_face](const DataVector& volume) {
      DataVector face(n_face);
      for (size_t k = 0; k < n_face; ++k) {
        face[k] = volume[k * n_radial];
      }
      return face;
    };
    const auto slice_tensor = [&face_slice](auto face_tensor,
                                            const auto& volume_tensor) {
      for (size_t storage = 0; storage < volume_tensor.size(); ++storage) {
        (*face_tensor)[storage] = face_slice(volume_tensor[storage]);
      }
    };

    tnsr::aa<DataVector, Dim> metric_face{};
    tnsr::aa<DataVector, Dim> pi_face{};
    tnsr::iaa<DataVector, Dim> phi_face{};
    Scalar<DataVector> gamma2_face{};
    tnsr::I<DataVector, Dim> coords_face{};
    slice_tensor(make_not_null(&metric_face),
                 db::get<gr::Tags::SpacetimeMetric<DataVector, Dim>>(box));
    slice_tensor(make_not_null(&pi_face),
                 db::get<gh::Tags::Pi<DataVector, Dim>>(box));
    slice_tensor(make_not_null(&phi_face),
                 db::get<gh::Tags::Phi<DataVector, Dim>>(box));
    slice_tensor(make_not_null(&gamma2_face),
                 db::get<gh::Tags::ConstraintGamma2>(box));
    slice_tensor(
        make_not_null(&coords_face),
        db::get<domain::Tags::Coordinates<Dim, Frame::Inertial>>(box));

    const ylm::Spherepack& ylm_transform = ylm::get_spherepack_cache(l_max);

    // The fit runs with zero drives: feeding the backward-difference rate
    // estimate back into the fit is an unstable loop (the optimal p shifts
    // to compensate a drive error by more than dt, so the estimate grows
    // geometrically — measured x1.33 per fit at FitInterval 0.5). Zero
    // drives in the fit is the offline pass-1, which differs from the
    // converged offline fit by only ~5% of |p| (findings §15b). The rates
    // are computed afterwards, purely as an output for the boundary
    // condition's model Pi.
    const std::array<double, num_map_parameters> zero_rates{};
    std::array<double, num_map_parameters> p_start{};
    if (state.valid) {
      p_start = state.p;
    }

    const FitResult result = fit_map_parameters(
        metric_face, pi_face, phi_face, gamma2_face, coords_face,
        ylm_transform, *config_opt, p_start, zero_rates);

    bool finite = true;
    for (size_t a = 0; a < num_map_parameters; ++a) {
      finite = finite and std::isfinite(gsl::at(result.p, a));
    }

    std::array<double, num_map_parameters> pdot_out{};
    if (finite) {
      db::mutate<Tags::MapParameters>(
          [&result, &time, &pdot_out](
              const gsl::not_null<MapParameterData*> data) {
            if (data->valid and time > data->last_fit_time) {
              const double dt = time - data->last_fit_time;
              for (size_t a = 0; a < num_map_parameters; ++a) {
                gsl::at(pdot_out, a) =
                    (gsl::at(result.p, a) - gsl::at(data->p, a)) / dt;
              }
              data->previous_fit_time = data->last_fit_time;
              data->p_previous = data->p;
            }
            data->last_fit_time = time;
            data->p = result.p;
            data->pdot = pdot_out;
            data->valid = true;
          },
          make_not_null(&box));
    }

    // diagnostic row: time, p (13), pdot (13), residuals, iterations
    // (rows with non-finite entries are recorded but do not update the
    // state, so the warm start stays on the last good fit)
    auto& writer = Parallel::get_parallel_component<
        observers::ObserverWriter<Metavariables>>(cache);
    std::vector<std::string> legend{"Time"};
    static const std::array<std::string, num_map_parameters> names{
        {"qdot0", "b_x", "b_y", "b_z", "v_x", "v_y", "v_z", "s_xx", "s_xy",
         "s_xz", "s_yy", "s_yz", "s_zz"}};
    for (const auto& name : names) {
      legend.push_back(name);
    }
    for (const auto& name : names) {
      legend.push_back("dt_" + name);
    }
    legend.emplace_back("ResidualInitial");
    legend.emplace_back("ResidualFinal");
    legend.emplace_back("Iterations");
    std::vector<double> row;
    row.reserve(3 * num_map_parameters);
    row.push_back(time);
    for (size_t a = 0; a < num_map_parameters; ++a) {
      row.push_back(gsl::at(result.p, a));
    }
    for (size_t a = 0; a < num_map_parameters; ++a) {
      row.push_back(gsl::at(pdot_out, a));
    }
    row.push_back(result.residual_initial);
    row.push_back(result.residual_final);
    row.push_back(static_cast<double>(result.iterations));
    Parallel::threaded_action<
        observers::ThreadedActions::WriteReductionDataRow>(
        writer[0], std::string{"/WorldtubeMatcher"}, std::move(legend),
        std::make_tuple(std::move(row)));

    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
}  // namespace gh::Worldtube::Actions
