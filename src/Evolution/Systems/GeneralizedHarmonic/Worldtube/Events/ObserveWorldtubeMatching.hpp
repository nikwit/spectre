// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Creators/Tags/Domain.hpp"
#include "Domain/Domain.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Tags.hpp"
#include "Domain/TagsTimeDependent.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/KretschmannFaceData.hpp"
#include "IO/Observer/Helpers.hpp"
#include "IO/Observer/ObservationId.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "IO/Observer/ReductionActions.hpp"
#include "IO/Observer/TypeOfObservation.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Options/Auto.hpp"
#include "Options/String.hpp"
#include "Parallel/ArrayComponentId.hpp"
#include "Parallel/ArrayIndex.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/Local.hpp"
#include "Parallel/Reduction.hpp"
#include "Parallel/TypeTraits.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Time/Tags/Time.hpp"
#include "Utilities/Functional.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::worldtube::Events {
/// The reduction data of `ObserveWorldtubeMatching`
using MatchingReductionData = Parallel::ReductionData<
    // Time
    Parallel::ReductionDatum<double, funcl::AssertEqual<>>,
    // Number of excision faces
    Parallel::ReductionDatum<size_t, funcl::Plus<>>,
    // Max |Psi0| of the face curvature
    Parallel::ReductionDatum<double, funcl::Max<>>,
    // Max |E|, |B|
    Parallel::ReductionDatum<double, funcl::Max<>>,
    Parallel::ReductionDatum<double, funcl::Max<>>,
    // Min Re Psi2^K, Max Re Psi2^K, Max |Im Psi2^K|
    Parallel::ReductionDatum<double, funcl::Min<>>,
    Parallel::ReductionDatum<double, funcl::Max<>>,
    Parallel::ReductionDatum<double, funcl::Max<>>,
    // Max type-D residual
    Parallel::ReductionDatum<double, funcl::Max<>>,
    // Max |Psi0| of the TypeD target, max |target - Psi0|
    Parallel::ReductionDatum<double, funcl::Max<>>,
    Parallel::ReductionDatum<double, funcl::Max<>>,
    // Max |Psi0| of the Quadrupole target, max |target - Psi0|
    Parallel::ReductionDatum<double, funcl::Max<>>,
    Parallel::ReductionDatum<double, funcl::Max<>>,
    // Psi4 fit relative residual
    Parallel::ReductionDatum<double, funcl::Max<>>,
    // Max |tanh(rapidity)|, then min and max rapidity
    Parallel::ReductionDatum<double, funcl::Max<>>,
    Parallel::ReductionDatum<double, funcl::Min<>>,
    Parallel::ReductionDatum<double, funcl::Max<>>,
    // Min, max measured radius
    Parallel::ReductionDatum<double, funcl::Min<>>,
    Parallel::ReductionDatum<double, funcl::Max<>>,
    // Max |Psi0| of the order-two target from the relaxed moments
    Parallel::ReductionDatum<double, funcl::Max<>>>;

/*!
 * \brief Diagnostics of the worldtube curvature matching on the excision
 * faces, as a reduction over the elements abutting an excision sphere.
 *
 * \details On each excision face the Weyl curvature is rebuilt from the
 * evolved variables exactly as the `gh::BoundaryConditions::WorldtubeTypeD`
 * boundary condition does, and the type-D (order zero) model is evaluated;
 * with a `Mass`, also the quadrupole (order two) model, using the
 * `gh::worldtube::KretschmannFaceData` of the element updated to the
 * observation time. Writes the columns
 *
 * - Time
 * - NumberOfFaces: the number of excision faces contributing
 * - MaxAbsPsi0: the largest \f$|\Psi_0|\f$ of the face curvature in the
 *   adapted tetrad, i.e. the radiation actually entering the domain
 * - MaxAbsElectric, MaxAbsMagnetic: \f$\max (E_{ij}E^{ij})^{1/2}\f$ and the
 *   same for \f$B\f$
 * - MinReCoulomb, MaxReCoulomb, MaxAbsImCoulomb: the range of the Kinnersley
 *   Coulomb scalar \f$\Psi_2^K = -3J/I\f$, which is \f$-M/r^3\f$ for a hole
 *   at rest
 * - MaxTypeDResidual: the largest pointwise relative misfit of the type-D
 *   solve, \f$\|\Psi^{\rm pred} - \Psi\| / \|\Psi\|\f$ over the five scalars;
 *   the size of the non-type-D content of the face curvature
 * - MaxAbsPsi0TypeD, MaxTypeDMismatch: the largest type-D target
 *   \f$|6 b^2 \Psi_2^K|\f$ and the largest \f$|\Psi_0^{\rm target} -
 *   \Psi_0|\f$, the fixed-point residual the boundary condition drives to zero
 * - MaxAbsPsi0Quadrupole, MaxQuadrupoleMismatch: the same for the order-two
 *   target (NaN without `Mass`)
 * - MaxAbsTanhRapidity: the largest \f$|\tanh\eta|\f$ of the invariant boost
 *   over the face. The boost is physical only below 1; where the Kretschmann
 *   gradient is under-resolved (it is a third derivative of the metric) this
 *   exceeds 1, the order-two model cannot be evaluated, and the remaining
 *   order-two columns are NaN for that face.
 * - Psi4FitRelativeResidual: the relative residual of the tidal fit to the
 *   pulled-back \f$\Psi_4\f$ over the face, the largest over the faces
 * - MinRapidity, MaxRapidity: the range of the invariant boost rapidity
 * - MinMeasuredRadius, MaxMeasuredRadius: the range of
 * - MaxAbsPsi0QuadrupoleImposed: the largest order-two target built from the
 *   relaxed moments the boundary condition imposes
 *   (`KretschmannFaceData::filtered_moments`, NaN when the face carries no
 *   order-two condition), as opposed to the instantaneous fit above
 * \f$(-M/\Psi_2^K)^{1/3}\f$
 *
 * Elements abutting no excision sphere neither register nor contribute.
 */
class ObserveWorldtubeMatching : public Event {
 private:
  using ReductionData = MatchingReductionData;

 public:
  struct SubfileName {
    using type = std::string;
    static constexpr Options::String help = {
        "The name of the subfile inside the HDF5 file without an extension and "
        "without a preceding '/'."};
  };
  struct Mass {
    using type = Options::Auto<double, Options::AutoLabel::None>;
    static constexpr Options::String help = {
        "The mass of the excised hole. With a mass the order-two (Quadrupole) "
        "model is evaluated as well; with None only the type-D model is."};
  };

  /// \cond
  explicit ObserveWorldtubeMatching(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(ObserveWorldtubeMatching);  // NOLINT
  /// \endcond

  using options = tmpl::list<SubfileName, Mass>;
  static constexpr Options::String help =
      "Observe diagnostics of the worldtube curvature matching on the "
      "excision faces: the entering Psi0, the Weyl curvature, the Coulomb "
      "scalar, the type-D residual, the mismatch between the model targets "
      "and the face's Psi0, and with a Mass the Psi4 fit residual, the "
      "rapidity and the measured radius of the order-two model.";

  ObserveWorldtubeMatching() = default;
  explicit ObserveWorldtubeMatching(const std::string& subfile_name,
                                    std::optional<double> mass);

  using observed_reduction_data_tags =
      observers::make_reduction_data_tags<tmpl::list<ReductionData>>;

  using compute_tags_for_observation_box = tmpl::list<>;
  using return_tags = tmpl::list<>;
  using argument_tags = tmpl::list<
      ::Tags::Time, gr::Tags::SpacetimeMetric<DataVector, 3>,
      gh::Tags::Pi<DataVector, 3>, gh::Tags::Phi<DataVector, 3>,
      domain::Tags::Mesh<3>,
      domain::Tags::InverseJacobian<3, Frame::ElementLogical, Frame::Inertial>,
      domain::Tags::Element<3>, domain::Tags::Domain<3>,
      domain::Tags::MeshVelocity<3>, worldtube::Tags::KretschmannFaceData<3>>;

  /// The reduction contribution of one element, or none if the element abuts
  /// no excision sphere
  std::optional<ReductionData> compute_reduction_data(
      double time,
      const tnsr::aa<DataVector, 3, Frame::Inertial>& spacetime_metric,
      const tnsr::aa<DataVector, 3, Frame::Inertial>& pi,
      const tnsr::iaa<DataVector, 3, Frame::Inertial>& phi, const Mesh<3>& mesh,
      const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                            Frame::Inertial>& inverse_jacobian,
      const Element<3>& element, const Domain<3>& domain,
      const std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>&
          mesh_velocity,
      const KretschmannFaceData<3>& face_data) const;

  static std::vector<std::string> legend();

  template <typename ArrayIndex, typename ParallelComponent,
            typename Metavariables>
  void operator()(
      const double time,
      const tnsr::aa<DataVector, 3, Frame::Inertial>& spacetime_metric,
      const tnsr::aa<DataVector, 3, Frame::Inertial>& pi,
      const tnsr::iaa<DataVector, 3, Frame::Inertial>& phi, const Mesh<3>& mesh,
      const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                            Frame::Inertial>& inverse_jacobian,
      const Element<3>& element, const Domain<3>& domain,
      const std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>&
          mesh_velocity,
      const KretschmannFaceData<3>& face_data,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& array_index, const ParallelComponent* const /*meta*/,
      const ObservationValue& observation_value) const {
    auto reduction_data = compute_reduction_data(
        time, spacetime_metric, pi, phi, mesh, inverse_jacobian, element,
        domain, mesh_velocity, face_data);
    if (not reduction_data.has_value()) {
      return;
    }
    auto& local_observer = *Parallel::local_branch(
        Parallel::get_parallel_component<
            tmpl::conditional_t<Parallel::is_nodegroup_v<ParallelComponent>,
                                observers::ObserverWriter<Metavariables>,
                                observers::Observer<Metavariables>>>(cache));
    observers::ObservationId observation_id{observation_value.value,
                                            subfile_path_ + ".dat"};
    Parallel::ArrayComponentId array_component_id{
        std::add_pointer_t<ParallelComponent>{nullptr},
        Parallel::ArrayIndex<ArrayIndex>(array_index)};
    if constexpr (Parallel::is_nodegroup_v<ParallelComponent>) {
      Parallel::threaded_action<
          observers::ThreadedActions::CollectReductionDataOnNode>(
          local_observer, std::move(observation_id),
          std::move(array_component_id), subfile_path_, legend(),
          std::move(*reduction_data));
    } else {
      Parallel::simple_action<observers::Actions::ContributeReductionData>(
          local_observer, std::move(observation_id),
          std::move(array_component_id), subfile_path_, legend(),
          std::move(*reduction_data));
    }
  }

  using observation_registration_tags = tmpl::list<::Tags::DataBox>;

  template <typename DbTagsList>
  std::optional<
      std::pair<observers::TypeOfObservation, observers::ObservationKey>>
  get_observation_type_and_key_for_registration(
      const db::DataBox<DbTagsList>& box) const {
    if (not excision_face_direction(
                db::get<domain::Tags::Domain<3>>(box).excision_spheres(),
                db::get<domain::Tags::Element<3>>(box))
                .has_value()) {
      return std::nullopt;
    }
    return {{observers::TypeOfObservation::Reduction,
             observers::ObservationKey(subfile_path_ + ".dat")}};
  }

  using is_ready_argument_tags = tmpl::list<>;

  template <typename Metavariables, typename ArrayIndex, typename Component>
  bool is_ready(Parallel::GlobalCache<Metavariables>& /*cache*/,
                const ArrayIndex& /*array_index*/,
                const Component* const /*meta*/) const {
    return true;
  }

  bool needs_evolved_variables() const override { return true; }

  const std::string& subfile_path() const { return subfile_path_; }
  const std::optional<double>& mass() const { return mass_; }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) override;

 private:
  std::string subfile_path_;
  std::optional<double> mass_{};
};
}  // namespace gh::worldtube::Events
