// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <cstddef>
#include <optional>
#include <type_traits>

#include "Evolution/Systems/CurvedScalarWave/Worldtube/Tags.hpp"

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Slice.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/ExcisionSphere.hpp"
#include "Domain/Structure/IndexToSliceAt.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/PunctureField.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/KerrSchild.hpp"
#include "PointwiseFunctions/GeneralRelativity/Christoffel.hpp"
#include "PointwiseFunctions/GeneralRelativity/DerivativesOfSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/InverseSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeMetric.hpp"
#include "Time/Tags.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/Gsl.hpp"

namespace CurvedScalarWave::Worldtube::Tags {

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsuggest-attribute=noreturn"
#endif  // defined(__GNUC__) && !defined(__clang__)
template <size_t Dim, typename Frame, bool Centered>
void FaceCoordinatesCompute<Dim, Frame, Centered>::function(
    const gsl::not_null<std::optional<tnsr::I<DataVector, Dim, Frame>>*> result,
    const ::ExcisionSphere<Dim>& excision_sphere, const Element<Dim>& element,
    const tnsr::I<DataVector, Dim, Frame>& coords, const Mesh<Dim>& mesh) {
  const auto direction = excision_sphere.abutting_direction(element.id());
  if (direction.has_value()) {
    ASSERT(
        mesh.quadrature(direction.value().dimension()) ==
            Spectral::Quadrature::GaussLobatto,
        "Expected GaussLobatto quadrature. Other quadratures are disabled "
        "because interpolating the coordinates incurs an unnecessary error.");
    const size_t grid_size =
        mesh.slice_away(direction->dimension()).number_of_grid_points();
    if (result->has_value()) {
      destructive_resize_components(make_not_null(&(result->value())),
                                    grid_size);
    } else {
      result->emplace(grid_size);
    }
    data_on_slice(make_not_null(&(result->value())), coords, mesh.extents(),
                  direction.value().dimension(),
                  index_to_slice_at(mesh.extents(), direction.value()));
    if constexpr (Centered) {
      if constexpr (not std::is_same_v<Frame, ::Frame::Grid>) {
        ERROR("Should be grid frame");
      }
      for (size_t i = 0; i < Dim; ++i) {
        result->value().get(i) -= excision_sphere.center().get(i);
      }
    }
  } else {
    result->reset();
  }
}

template <size_t Dim, typename Frame, bool Centered>
void FaceCoordinatesCompute<Dim, Frame, Centered>::function(
    const gsl::not_null<
        std::optional<tnsr::I<DataVector, Dim, ::Frame::Inertial>>*>
        result,
    const ::ExcisionSphere<Dim>& excision_sphere, const Element<Dim>& element,
    const tnsr::I<DataVector, Dim, ::Frame::Inertial>& inertial_coords,
    const Mesh<Dim>& mesh,
    const std::array<tnsr::I<double, Dim, ::Frame::Inertial>, 2>&
        particle_position) {
  if constexpr (not(Centered and std::is_same_v<Frame, ::Frame::Inertial>)) {
    ERROR("Should be centered in inertial frame");
  }
  const auto direction = excision_sphere.abutting_direction(element.id());
  if (direction.has_value()) {
    ASSERT(
        mesh.quadrature(direction.value().dimension()) ==
            Spectral::Quadrature::GaussLobatto,
        "Expected GaussLobatto quadrature. Other quadratures are disabled "
        "because interpolating the coordinates incurs an unnecessary error.");
    const size_t grid_size =
        mesh.slice_away(direction->dimension()).number_of_grid_points();
    if (result->has_value()) {
      destructive_resize_components(make_not_null(&(result->value())),
                                    grid_size);
    } else {
      result->emplace(grid_size);
    }
    data_on_slice(make_not_null(&(result->value())), inertial_coords,
                  mesh.extents(), direction.value().dimension(),
                  index_to_slice_at(mesh.extents(), direction.value()));
    for (size_t i = 0; i < Dim; ++i) {
      result->value().get(i) -= particle_position[0].get(i);
    }
  } else {
    result->reset();
  }
}
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif  // defined(__GNUC__) && !defined(__clang__)

template <size_t Dim>
void PunctureFieldCompute<Dim>::function(
    const gsl::not_null<return_type*> result,
    const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
        inertial_face_coords_centered,
    const ::ExcisionSphere<Dim>& excision_sphere, const double time,
    const size_t expansion_order,
    const std::array<tnsr::I<double, Dim, ::Frame::Inertial>, 2>&
        particle_position_velocity,
    const tnsr::I<double, Dim>& particle_acceleration, const double charge,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time) {
  if (inertial_face_coords_centered.has_value()) {
    if (not result->has_value()) {
      result->emplace(get<0>(inertial_face_coords_centered.value()).size());
    }
    puncture_field_generic_1(
        make_not_null(&(result->value())),
        inertial_face_coords_centered.value(), particle_position_velocity[0],
        particle_position_velocity[1], particle_acceleration, 1.);
    result->value() *= charge;
  } else {
    result->reset();
  }
}

void ConstraintGamma1Compute::function(
    gsl::not_null<Scalar<DataVector>*> gamma1,
    const tnsr::I<DataVector, 3, Frame::Inertial>& coords,
    const std::array<tnsr::I<double, 3, Frame::Inertial>, 2>& pos_vel) {
  get(*gamma1).destructive_resize(get<0>(coords).size());

  get(*gamma1) = 0.;
}

void ConstraintGamma2Compute::function(
    gsl::not_null<Scalar<DataVector>*> gamma2,
    const tnsr::I<DataVector, 3, Frame::Inertial>& coords,
    const std::array<tnsr::I<double, 3, Frame::Inertial>, 2>& pos_vel) {
  get(*gamma2).destructive_resize(get<0>(coords).size());
  auto centered_coords = coords;
  for (size_t i = 0; i < 3; ++i) {
    centered_coords.get(i) -= pos_vel[0].get(i);
  }
  const double amplitude = 10.;
  const double sigma = 1e-1;
  const double constant = 1e-3;
  const auto radius = magnitude(centered_coords);
  get(*gamma2) = amplitude * exp(-square(sigma * radius.get())) + constant;
  //get(*gamma2) += 30. * exp(-square( 2. * radius.get()));
}

template <size_t Dim>
void ParticlePositionVelocityCompute<Dim>::function(
    gsl::not_null<std::array<tnsr::I<double, Dim, Frame::Inertial>, 2>*>
        position,
    const ::ExcisionSphere<Dim>& excision_sphere,
    const domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, 3>& maps,
    const double time,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time) {
  auto values = maps.coords_frame_velocity_jacobians(excision_sphere.center(),
                                                     time, functions_of_time);
  (*position)[0] = std::move(std::get<0>(values));
  (*position)[1] = std::move(std::get<3>(values));
}

template <size_t Dim>
void ParticleAccelerationCompute<Dim>::function(
    gsl::not_null<tnsr::I<double, Dim, Frame::Inertial>*> acceleration,
    const std::array<tnsr::I<double, Dim, Frame::Inertial>, 2>&
        position_velocity,
    const gr::Solutions::KerrSchild& background_spacetime) {
  const auto& inertial_particle_position = position_velocity[0];
  const auto spacetime_vars = background_spacetime.variables(
      inertial_particle_position, 0.,
      tmpl::list<
          gr::Tags::Lapse<double>,
          gr::Tags::Shift<double, Dim, Frame::Inertial>,
          gr::Tags::SpatialMetric<double, Dim, Frame::Inertial>,
          gr::Tags::InverseSpatialMetric<double, Dim, Frame::Inertial>,
          ::Tags::dt<gr::Tags::Lapse<double>>,
          ::Tags::deriv<gr::Tags::Lapse<double>, tmpl::size_t<Dim>,
                        Frame::Inertial>,
          ::Tags::dt<gr::Tags::Shift<double, Dim, Frame::Inertial>>,
          ::Tags::deriv<gr::Tags::Shift<double, Dim, Frame::Inertial>,
                        tmpl::size_t<Dim>, Frame::Inertial>,
          ::Tags::dt<gr::Tags::SpatialMetric<double, Dim, Frame::Inertial>>,
          ::Tags::deriv<gr::Tags::SpatialMetric<double, Dim, Frame::Inertial>,
                        tmpl::size_t<Dim>, Frame::Inertial>>{});
  const auto inverse_spacetime_metric_inertial = gr::inverse_spacetime_metric(
      get<gr::Tags::Lapse<double>>(spacetime_vars),
      get<gr::Tags::Shift<double, Dim, Frame::Inertial>>(spacetime_vars),
      get<gr::Tags::InverseSpatialMetric<double, Dim, Frame::Inertial>>(
          spacetime_vars));
  const auto spacetime_metric_inertial = gr::spacetime_metric(
      get<gr::Tags::Lapse<double>>(spacetime_vars),
      get<gr::Tags::Shift<double, Dim, Frame::Inertial>>(spacetime_vars),
      get<gr::Tags::SpatialMetric<double, Dim, Frame::Inertial>>(
          spacetime_vars));
  const auto d_spacetime_metric = gr::derivatives_of_spacetime_metric(
      get<gr::Tags::Lapse<double>>(spacetime_vars),
      get<::Tags::dt<gr::Tags::Lapse<double>>>(spacetime_vars),
      get<::Tags::deriv<gr::Tags::Lapse<double>, tmpl::size_t<Dim>,
                        Frame::Inertial>>(spacetime_vars),
      get<gr::Tags::Shift<double, Dim, Frame::Inertial>>(spacetime_vars),
      get<::Tags::dt<gr::Tags::Shift<double, Dim, Frame::Inertial>>>(
          spacetime_vars),
      get<::Tags::deriv<gr::Tags::Shift<double, Dim, Frame::Inertial>,
                        tmpl::size_t<Dim>, Frame::Inertial>>(spacetime_vars),
      get<gr::Tags::SpatialMetric<double, Dim, Frame::Inertial>>(
          spacetime_vars),
      get<::Tags::dt<gr::Tags::SpatialMetric<double, Dim, Frame::Inertial>>>(
          spacetime_vars),
      get<::Tags::deriv<gr::Tags::SpatialMetric<double, Dim, Frame::Inertial>,
                        tmpl::size_t<Dim>, Frame::Inertial>>(spacetime_vars));

  const auto christoffel = gr::christoffel_second_kind(
      d_spacetime_metric, inverse_spacetime_metric_inertial);
  const auto& particle_velocity = position_velocity[1];

  for (size_t i = 0; i < Dim; ++i) {
    acceleration->get(i) = particle_velocity.get(i) * christoffel.get(0, 0, 0) -
                           christoffel.get(i + 1, 0, 0);
    for (size_t j = 0; j < Dim; ++j) {
      acceleration->get(i) +=
          2. * particle_velocity.get(j) *
          (particle_velocity.get(i) * christoffel.get(0, j + 1, 0) -
           christoffel.get(i + 1, j + 1, 0));
      for (size_t k = 0; k < Dim; ++k) {
        acceleration->get(i) +=
            particle_velocity.get(j) * particle_velocity.get(k) *
            (particle_velocity.get(i) * christoffel.get(0, j + 1, k + 1) -
             christoffel.get(i + 1, j + 1, k + 1));
      }
    }
  }
}

template struct ParticlePositionVelocityCompute<3>;
template struct ParticleAccelerationCompute<3>;
template struct PunctureFieldCompute<3>;

template struct FaceCoordinatesCompute<3, Frame::Grid, true>;
template struct FaceCoordinatesCompute<3, Frame::Grid, false>;
template struct FaceCoordinatesCompute<3, Frame::Inertial, true>;
template struct FaceCoordinatesCompute<3, Frame::Inertial, false>;

}  // namespace CurvedScalarWave::Worldtube::Tags
