// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>

#include "ControlSystem/UpdateFunctionOfTime.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "Evolution/Systems/CurvedScalarWave/BackgroundSpacetime.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Inboxes.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Tags.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/GlobalCache.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Tags.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/KerrSchild.hpp"
#include "PointwiseFunctions/GeneralRelativity/Christoffel.hpp"
#include "PointwiseFunctions/GeneralRelativity/DerivativesOfSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/InverseSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeMetric.hpp"
#include "Time/Actions/ChangeSlabSize.hpp"
#include "Time/Tags.hpp"
#include "Time/TimeStepId.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace CurvedScalarWave::Worldtube::Actions {

/*!
 * \brief Waits for the data from all neighboring elements and changes the slab
 * size if a change in the global time step is detected.
 * \details We check the slab size of the time step id sent by the elements. If
 * this is different from the slab size currently used by the worldtube
 * singleton, we assume a global slab size change has occurred in the elements
 * and adjust the worldtube slab size accordingly.
 */
struct UpdateFunctionsOfTime {
  static constexpr size_t Dim = 3;
  // using compute_tags =
  // tmpl::list<Tags::InertialParticlePositionCompute<Dim>>;
  using inbox_tags = tmpl::list<>;

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& box,
      tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& /*array_index*/, ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    const auto& time = db::get<::Tags::Time>(box);
    const auto& time_step = db::get<::Tags::TimeStep>(box);
    const auto& functions_of_time =
        Parallel::get<::domain::Tags::FunctionsOfTime>(cache);
    const auto& excision_sphere = db::get<Tags::ExcisionSphere<Dim>>(box);
    const auto& maps = Parallel::get<Tags::WorldtubeCoordinateMaps>(cache);
    const std::string rot_function_of_time_name = "Rotation";
    const std::string expansion_fot_name = "Expansion";
    const std::string size_a_fot_name = "SizeA";
    const std::string size_b_fot_name = "SizeB";

    const auto& rot_function_of_time =
        functions_of_time.at(rot_function_of_time_name);
    const double current_fot_expiration_time =
        rot_function_of_time->time_bounds()[1];
    if (time > current_fot_expiration_time) {
      const auto& inertial_particle_position =
          maps(excision_sphere.center(), db::get<Tags::PreviousTime>(box),
               functions_of_time);
      const auto& particle_velocity =
          std::get<3>(maps.coords_frame_velocity_jacobians(
              excision_sphere.center(), db::get<Tags::PreviousTime>(box),
              functions_of_time));
      const auto& kerr_schild =
          db::get<CurvedScalarWave::Tags::BackgroundSpacetime<
              gr::Solutions::KerrSchild>>(box);
      const auto spacetime_vars = kerr_schild.variables(
          inertial_particle_position, time,
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
              ::Tags::deriv<
                  gr::Tags::SpatialMetric<double, Dim, Frame::Inertial>,
                  tmpl::size_t<Dim>, Frame::Inertial>>{});
      const auto inverse_spacetime_metric = gr::inverse_spacetime_metric(
          get<gr::Tags::Lapse<double>>(spacetime_vars),
          get<gr::Tags::Shift<double, Dim, Frame::Inertial>>(spacetime_vars),
          get<gr::Tags::InverseSpatialMetric<double, Dim, Frame::Inertial>>(
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
                            tmpl::size_t<Dim>, Frame::Inertial>>(
              spacetime_vars),
          get<gr::Tags::SpatialMetric<double, Dim, Frame::Inertial>>(
              spacetime_vars),
          get<::Tags::dt<
              gr::Tags::SpatialMetric<double, Dim, Frame::Inertial>>>(
              spacetime_vars),
          get<::Tags::deriv<
              gr::Tags::SpatialMetric<double, Dim, Frame::Inertial>,
              tmpl::size_t<Dim>, Frame::Inertial>>(spacetime_vars));

      const auto christoffel = gr::christoffel_second_kind(
          d_spacetime_metric, inverse_spacetime_metric);

      tnsr::I<double, Dim> particle_acceleration{};
      for (size_t i = 0; i < Dim; ++i) {
        particle_acceleration.get(i) =
            particle_velocity.get(i) * christoffel.get(0, 0, 0) -
            christoffel.get(i + 1, 0, 0);
        for (size_t j = 0; j < Dim; ++j) {
          particle_acceleration.get(i) +=
              2. * particle_velocity.get(j) *
              (particle_velocity.get(i) * christoffel.get(0, j + 1, 0) -
               christoffel.get(i + 1, j + 1, 0));
          for (size_t k = 0; k < Dim; ++k) {
            particle_acceleration.get(i) +=
                particle_velocity.get(j) * particle_velocity.get(k) *
                (particle_velocity.get(i) * christoffel.get(0, j + 1, k + 1) -
                 christoffel.get(i + 1, j + 1, k + 1));
          }
        }
      }
      const double& x = get<0>(inertial_particle_position);
      const double& y = get<1>(inertial_particle_position);
      const double& xdot = get<0>(particle_velocity);
      const double& ydot = get<1>(particle_velocity);
      const double& xddot = get<0>(particle_acceleration);
      const double& yddot = get<1>(particle_acceleration);
      const double r = hypot(x, y);
      const double radial_acc =
          (xddot * x + yddot * y + xdot * xdot + ydot * ydot) / r -
          square(xdot * x + ydot * y) / cube(r);
      const double angular_acc = (x * yddot - y * xddot) / square(r) -
                                 2. * (x * xdot + y * ydot) *
                                     (x * ydot - y * xdot) / square(square(r));
      DataVector new_angular_acc(3, 0.);
      new_angular_acc.at(2) = angular_acc;

      const double envelope_radius = 50.;
      const double expansion_acc =
          -10. * radial_acc * sqrt(4. * M_PI) * envelope_radius;

      const auto& expansion_fot = functions_of_time.at(expansion_fot_name);
      const double object_radius = 1.;

      const auto expansion_vals =
          expansion_fot->func_and_deriv(db::get<Tags::PreviousTime>(box));
      const double sqrt_4_pi = sqrt(4. * M_PI);
      const double fac_1 =
          1. / (1. - expansion_vals.at(0)[0] / (sqrt_4_pi * envelope_radius));

      double compression_acc =
          -object_radius / envelope_radius * square(fac_1) *
          (expansion_acc + 2. * square(expansion_vals.at(1)[0]) * fac_1 /
                               sqrt_4_pi / envelope_radius);

      DataVector new_expansion_acc(1, expansion_acc);
      DataVector new_compression_acc_b(1, compression_acc);
      DataVector new_compression_acc_a(1, compression_acc);

      const double new_fot_expiration_time = time + time_step.value() * 0.5;

      Parallel::printf(MakeString{} << "Mutating Time from "
                                    << current_fot_expiration_time << " to "
                                    << new_fot_expiration_time << "\n");
      Parallel::mutate<::domain::Tags::FunctionsOfTime,
                       control_system::UpdateFunctionOfTime>(
          cache, rot_function_of_time_name, current_fot_expiration_time,
          new_angular_acc, new_fot_expiration_time);
      Parallel::mutate<::domain::Tags::FunctionsOfTime,
                       control_system::UpdateFunctionOfTime>(
          cache, expansion_fot_name, current_fot_expiration_time,
          new_expansion_acc, new_fot_expiration_time);
      Parallel::mutate<::domain::Tags::FunctionsOfTime,
                       control_system::UpdateFunctionOfTime>(
          cache, size_a_fot_name, current_fot_expiration_time,
          new_compression_acc_a, new_fot_expiration_time);
      Parallel::mutate<::domain::Tags::FunctionsOfTime,
                       control_system::UpdateFunctionOfTime>(
          cache, size_b_fot_name, current_fot_expiration_time,
          new_compression_acc_b, new_fot_expiration_time);
    }
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
}  // namespace CurvedScalarWave::Worldtube::Actions
