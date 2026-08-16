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
#include "Domain/Domain.hpp"
#include "Domain/ExcisionSphere.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Tags.hpp"
#include "Domain/TagsTimeDependent.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Matcher.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Worldtube/Tags.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "IO/Observer/ReductionActions.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Spherepack.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/SpherepackCache.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/AffineMappedHarmonicSchwarzschild.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ConstraintDampingTags.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ExtrinsicCurvature.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/Ricci.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeNormalVector.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylElectric.hpp"
#include "Time/Tags/Time.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::Worldtube::Actions {
/*!
 * \brief Element-local online worldtube matching (runs in the step loop).
 *
 * On the element owning the excision face — with the spherical-harmonic
 * angular basis there is exactly one, self-identified as the element whose
 * block abuts the configured excision sphere and that owns the abutting
 * external boundary — this action slices the evolved
 * fields to the excision face, fits the 13 first-order affine-map parameters
 * per `gh::Worldtube::fit_map_parameters`, stores the result in
 * `Tags::MapParameters`, and writes a diagnostic row to
 * `/WorldtubeMatcher.dat` in the reductions file.
 *
 * No-op on every other element, when `Tags::Matcher` is `None`, and when the
 * time has not advanced by `FitInterval` since the last fit (this also
 * limits the action to one fit per step under substepping time steppers).
 * In the default algebraic mode the result is a strict slow-time first-order
 * state: coefficient rates are neither fitted nor extrapolated.  The
 * derivative ODE modes are retained as explicitly higher-order experiments.
 */
struct FitMapParameters {
  using const_global_cache_tags = tmpl::list<Tags::Matcher>;

  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*component*/) {
    static constexpr size_t Dim = 3;
    const auto& config_opt = db::get<Tags::Matcher>(box);
    if (not config_opt.has_value()) {
      return {Parallel::AlgorithmExecution::Continue, std::nullopt};
    }
    const auto& element = db::get<domain::Tags::Element<Dim>>(box);
    // The worldtube boundary element is the one whose block abuts the
    // configured excision sphere and that owns the abutting external
    // boundary (with the spherical-harmonic angular basis there is exactly
    // one such element per radial refinement level, and only the innermost
    // owns the external face).
    const auto& domain = db::get<domain::Tags::Domain<Dim>>(box);
    const auto& excision_spheres = domain.excision_spheres();
    const auto excision_sphere_it =
        excision_spheres.find(config_opt->excision_sphere_name);
    if (excision_sphere_it == excision_spheres.end()) {
      ERROR("The worldtube matcher's ExcisionSphereName '"
            << config_opt->excision_sphere_name
            << "' is not an excision sphere of the domain. The domain has "
            << excision_spheres.size() << " excision sphere(s).");
    }
    const std::optional<Direction<Dim>> abutting_direction =
        excision_sphere_it->second.abutting_direction(element.id());
    if (not abutting_direction.has_value() or
        element.external_boundaries().count(*abutting_direction) == 0) {
      return {Parallel::AlgorithmExecution::Continue, std::nullopt};
    }
    ASSERT(*abutting_direction == Direction<Dim>::lower_xi(),
           "The online worldtube matcher requires the excision face to be "
           "the lower boundary of the radial logical direction 0 "
           "(spherical-harmonic shell block convention), but the excision "
           "sphere abuts in direction "
               << *abutting_direction);
    const double time = db::get<::Tags::Time>(box);
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

    // slice a volume component at the configured radial collocation index
    // (0 = the excision face; radial is the fastest-varying index, matching
    // the Ylm filter's storage convention)
    const size_t radial_index = config_opt->fit_radial_index;
    ASSERT(radial_index < n_radial, "FitRadialIndex "
                                        << radial_index << " out of range for "
                                        << n_radial << " radial points");
    const auto shell_slice = [n_radial, n_face](const DataVector& volume,
                                                const size_t shell) {
      DataVector face(n_face);
      for (size_t k = 0; k < n_face; ++k) {
        face[k] = volume[k * n_radial + shell];
      }
      return face;
    };
    const auto face_slice = [&shell_slice,
                             radial_index](const DataVector& volume) {
      return shell_slice(volume, radial_index);
    };
    const auto slice_tensor = [&face_slice](auto face_tensor,
                                            const auto& volume_tensor) {
      for (size_t storage = 0; storage < volume_tensor.size(); ++storage) {
        (*face_tensor)[storage] = face_slice(volume_tensor[storage]);
      }
    };
    const auto slice_tensor_at = [&shell_slice](auto face_tensor,
                                                const auto& volume_tensor,
                                                const size_t shell) {
      for (size_t storage = 0; storage < volume_tensor.size(); ++storage) {
        (*face_tensor)[storage] = shell_slice(volume_tensor[storage], shell);
      }
    };

    tnsr::I<DataVector, Dim> coords_face{};
    slice_tensor(make_not_null(&coords_face),
                 db::get<domain::Tags::Coordinates<Dim, Frame::Inertial>>(box));

    const ylm::Spherepack& ylm_transform = ylm::get_spherepack_cache(l_max);
    const std::array<double, 3> current_worldtube_center =
        detail::worldtube_center(coords_face, ylm_transform);
    db::mutate<Tags::MapParameters>(
        [&current_worldtube_center](
            const gsl::not_null<MapParameterData*> data) {
          data->worldtube_center = current_worldtube_center;
          data->worldtube_center_valid = true;
        },
        make_not_null(&box));

    MatcherConfig fit_config = *config_opt;
    fit_config.center = current_worldtube_center;
    const auto& state = db::get<Tags::MapParameters>(box);
    if (config_opt->stepper_ode and config_opt->uplus_anchor == 0.) {
      // The center update above is still required by the boundary condition;
      // the parameter evolution itself is handled after the RHS computation.
      return {Parallel::AlgorithmExecution::Continue, std::nullopt};
    }
    // Algebraic fits change the ghost prescription.  Applying a new
    // prescription on an interior Runge--Kutta stage makes the stages of one
    // step use different boundary data and produces a strong positive
    // feedback in the clean projected solve.  Fit only at step boundaries;
    // the worldtube center above is still refreshed on every invocation.
    if (db::get<::Tags::TimeStepId>(box).substep() != 0) {
      return {Parallel::AlgorithmExecution::Continue, std::nullopt};
    }
    // FitInterval 0 means fit at every step boundary, with no time comparison
    // at all. The worldtube center above is updated even when this gate skips
    // the expensive field fit.
    if (config_opt->fit_interval > 0.) {
      if (config_opt->stepper_ode) {
        if (state.anchor_valid and
            time < state.anchor_time + config_opt->fit_interval - 1.0e-12) {
          return {Parallel::AlgorithmExecution::Continue, std::nullopt};
        }
      } else if (state.valid and time < state.last_fit_time +
                                            config_opt->fit_interval -
                                            1.0e-12) {
        return {Parallel::AlgorithmExecution::Continue, std::nullopt};
      }
    }

    tnsr::aa<DataVector, Dim> metric_face{};
    tnsr::aa<DataVector, Dim> pi_face{};
    tnsr::iaa<DataVector, Dim> phi_face{};
    Scalar<DataVector> gamma2_face{};
    slice_tensor(make_not_null(&metric_face),
                 db::get<gr::Tags::SpacetimeMetric<DataVector, Dim>>(box));
    slice_tensor(make_not_null(&pi_face),
                 db::get<gh::Tags::Pi<DataVector, Dim>>(box));
    slice_tensor(make_not_null(&phi_face),
                 db::get<gh::Tags::Phi<DataVector, Dim>>(box));
    slice_tensor(make_not_null(&gamma2_face),
                 db::get<gh::Tags::ConstraintGamma2>(box));

    // The value fit is the strict first-order slow-time path. It has no
    // coefficient-rate input or output: D_t p_(1) changes the metric only at
    // O(epsilon^2). The zeroth-order center is different: its O(epsilon)
    // velocity is retained kinematically. The derivative ODE branches below
    // remain separate, explicitly higher-order experiments.
    std::array<double, num_map_parameters> p_out{};
    std::array<double, num_map_parameters> pdot_out{};
    std::array<double, 3> center_offset_out{};
    std::array<double, 3> bulk_velocity_out = state.bulk_velocity;
    double residual_initial = 0.;
    double residual_final = 0.;
    double iterations = 0.;
    std::optional<FitResult> value_diagnostics{};
    std::optional<OrderOneFitResult> order_one_diagnostics{};
    bool order_one_improves_held_out = false;
    std::array<double, 3> track_velocity{};
    bool track_velocity_valid = false;

    // Trace pin: the TraceStrainPin constant, or measured live from the
    // curvature. In vacuum the Gauss-Bonnet scalar equals the Kretschmann
    // scalar; the magnetic contribution is quadratic in the perturbation
    // (B = 0 on the background, |B.B/E.E| ~ 1e-10 here) and is dropped, so
    // K ~ 8 E.E. The invariant harmonic radius rho = (48 M^2/K)^{1/6} - M
    // gives tr sigma / 3 = 1 - <rho>/<R> via the l=0 sphere means
    // (findings 12).
    double trace_pin = config_opt->trace_strain_pin;
    if (config_opt->kretschmann_trace_pin and
        state.trace_pin_time > std::numeric_limits<double>::lowest() and
        time <
            state.trace_pin_time + config_opt->trace_pin_interval - 1.0e-12) {
      trace_pin = state.trace_pin_value;
    } else if (config_opt->kretschmann_trace_pin) {
      const auto& inv_jacobian =
          db::get<domain::Tags::InverseJacobian<Dim, Frame::ElementLogical,
                                                Frame::Inertial>>(box);
      const auto deriv_phi_volume = partial_derivative(
          db::get<gh::Tags::Phi<DataVector, Dim>>(box), mesh, inv_jacobian);
      tnsr::ijaa<DataVector, Dim> deriv_phi_face{};
      slice_tensor(make_not_null(&deriv_phi_face), deriv_phi_volume);

      const auto inv4 = determinant_and_inverse(metric_face).second;
      const DataVector lapse_sq = -1. / get<0, 0>(inv4);
      Scalar<DataVector> lapse{sqrt(lapse_sq)};
      tnsr::I<DataVector, Dim> shift(n_face);
      for (size_t i = 0; i < 3; ++i) {
        shift.get(i) = lapse_sq * inv4.get(0, i + 1);
      }
      tnsr::II<DataVector, Dim> inv_spatial(n_face);
      for (size_t i = 0; i < 3; ++i) {
        for (size_t j = i; j < 3; ++j) {
          inv_spatial.get(i, j) =
              inv4.get(i + 1, j + 1) +
              inv4.get(0, i + 1) * inv4.get(0, j + 1) * lapse_sq;
        }
      }
      const auto normal_vec = gr::spacetime_normal_vector(lapse, shift);
      const auto kij = gh::extrinsic_curvature(normal_vec, pi_face, phi_face);
      const auto ricci =
          gh::spatial_ricci_tensor(phi_face, deriv_phi_face, inv_spatial);
      const auto weyl_e = gr::weyl_electric(ricci, kij, inv_spatial);
      const DataVector ee = get(gr::weyl_electric_scalar(weyl_e, inv_spatial));
      const DataVector gb = 8. * ee;
      const DataVector rho_gb =
          pow(48. * square(config_opt->mass) / gb, 1. / 6.) - config_opt->mass;
      DataVector coord_r(n_face, 0.);
      for (size_t i = 0; i < 3; ++i) {
        coord_r +=
            square(coords_face.get(i) - gsl::at(current_worldtube_center, i));
      }
      coord_r = sqrt(coord_r);
      // l = 0 sphere means via Spherepack (the collocation grid has no
      // uniform quadrature weights, so a plain mean is not an area mean)
      const auto sphere_mean = [&ylm_transform](const DataVector& field) {
        const DataVector spec = ylm_transform.phys_to_spec(field);
        DataVector only(spec.size(), 0.);
        ylm::SpherepackIterator iter(ylm_transform.l_max(),
                                     ylm_transform.m_max());
        iter.set(0, 0);
        only[iter()] = spec[iter()];
        const DataVector back = ylm_transform.spec_to_phys(only);
        double sum = 0.;
        for (size_t k = 0; k < back.size(); ++k) {
          sum += back[k];
        }
        return sum / static_cast<double>(back.size());
      };
      trace_pin = 1. - sphere_mean(rho_gb) / sphere_mean(coord_r);
      if (not std::isfinite(trace_pin)) {
        trace_pin = state.trace_pin_time > std::numeric_limits<double>::lowest()
                        ? state.trace_pin_value
                        : config_opt->trace_strain_pin;
      }
      // l = 1 content of the same field is the centre displacement: with the
      // hole at d relative to the sphere centre, rho_GB = |x - d| ~ r -
      // d.nhat, and nhat_z = cos(theta) = sqrt(4 pi / 3) Y_10, so
      // d_i = -sqrt(3 / 4 pi) c_{1m}. Spherepack stores the m != 0
      // coefficients a factor sqrt(2) small relative to unit-L2 harmonics
      // (findings 15m), which the weight below undoes. Open-loop only.
      std::array<double, 3> dipole{};
      {
        const DataVector spec = ylm_transform.phys_to_spec(rho_gb);
        ylm::SpherepackIterator iter(ylm_transform.l_max(),
                                     ylm_transform.m_max());
        const double norm = -std::sqrt(3. / (4. * M_PI));
        const double m_weight = std::sqrt(2.);
        // (l=1, m=0) -> z; the two m=1 coefficients -> x (cos) and y (sin)
        iter.set(1, 0);
        dipole[2] = norm * spec[iter()];
        iter.set(1, 1);
        dipole[0] = norm * m_weight * spec[iter()];
        iter.set(1, -1);
        dipole[1] = norm * m_weight * spec[iter()];
        for (size_t i = 0; i < 3; ++i) {
          if (not std::isfinite(gsl::at(dipole, i))) {
            gsl::at(dipole, i) = 0.;
          }
        }
      }
      db::mutate<Tags::MapParameters>(
          [&trace_pin, &dipole,
           &time](const gsl::not_null<MapParameterData*> data) {
            data->trace_pin_value = trace_pin;
            data->trace_pin_time = time;
            // backward difference against the previous sample: an
            // Hole--worldtube relative velocity. Add the finite difference of
            // WorldtubeCenter_i to compare with fitted inertial qdot^i.
            if (data->gb_dipole_time > std::numeric_limits<double>::lowest() and
                time > data->gb_dipole_time) {
              const double dt = time - data->gb_dipole_time;
              for (size_t i = 0; i < 3; ++i) {
                gsl::at(data->gb_dipole_velocity, i) =
                    (gsl::at(dipole, i) - gsl::at(data->gb_dipole, i)) / dt;
              }
              data->gb_dipole_previous = data->gb_dipole;
              data->gb_dipole_time_previous = data->gb_dipole_time;
            }
            data->gb_dipole = dipole;
            data->gb_dipole_time = time;
          },
          make_not_null(&box));
    }

    if (config_opt->stepper_ode) {
      // u^+ anchor: value-fit the map to the gauge projection of the
      // OUTGOING characteristic (normal_sign -1) — the channel the ghost
      // BC does not set, i.e. the data the ambient evolution feeds the
      // excision. Consumed by AdvanceMapParameterOde as the weak drive
      // -2 kappa (pdot - pdot_anchor) - kappa^2 (p - p_anchor).
      std::array<double, num_map_parameters> p_start{};
      const std::array<double, 3> center_offset_start{};
      if (state.anchor_valid) {
        p_start = state.anchor_p;
      }
      const FitResult result =
          fit_map_parameters(metric_face, pi_face, phi_face, gamma2_face,
                             coords_face, ylm_transform, fit_config, p_start,
                             center_offset_start, -1.0, trace_pin);
      bool finite = true;
      for (size_t a = 0; a < num_map_parameters; ++a) {
        finite = finite and std::isfinite(gsl::at(result.p, a));
      }
      if (finite) {
        db::mutate<Tags::MapParameters>(
            [&result, &time](const gsl::not_null<MapParameterData*> data) {
              if (data->anchor_valid) {
                data->anchor_p_previous = data->anchor_p;
                data->anchor_time_previous = data->anchor_time;
              }
              data->anchor_p = result.p;
              data->anchor_time = time;
              data->anchor_valid = true;
            },
            make_not_null(&box));
      }
      auto& writer = Parallel::get_parallel_component<
          observers::ObserverWriter<Metavariables>>(cache);
      std::vector<std::string> legend{"Time"};
      static const std::array<std::string, num_map_parameters> names{
          {"qdot0", "b_x", "b_y", "b_z", "v_x", "v_y", "v_z", "s_xx", "s_xy",
           "s_xz", "s_yy", "s_yz", "s_zz"}};
      for (const auto& name : names) {
        legend.push_back("anchor_" + name);
      }
      legend.emplace_back("ResidualInitial");
      legend.emplace_back("ResidualFinal");
      legend.emplace_back("Iterations");
      std::vector<double> row;
      row.reserve(num_map_parameters + 4);
      row.push_back(time);
      for (size_t a = 0; a < num_map_parameters; ++a) {
        row.push_back(gsl::at(result.p, a));
      }
      row.push_back(result.residual_initial);
      row.push_back(result.residual_final);
      row.push_back(static_cast<double>(result.iterations));
      Parallel::threaded_action<
          observers::ThreadedActions::WriteReductionDataRow>(
          writer[0], std::string{"/WorldtubeAnchor"}, std::move(legend),
          std::make_tuple(std::move(row)));
      return {Parallel::AlgorithmExecution::Continue, std::nullopt};
    }

    if (config_opt->second_order_ode) {
      if (config_opt->fit_center_offset or config_opt->rate_ode) {
        ERROR(
            "WorldtubeMatcher: SecondOrderOde cannot be combined with "
            "RateOde or FitCenterOffset.");
      }
      // Second-order mode: fit pddot from the evolution equations (the dt
      // variables hold the GH right-hand sides; at the top of the step they
      // are one step stale, a small fraction of FitInterval) and integrate
      // (p, pdot) by velocity Verlet from (0, 0) at the first fit.
      tnsr::aa<DataVector, Dim> dt_metric_face{};
      tnsr::aa<DataVector, Dim> dt_pi_face{};
      tnsr::iaa<DataVector, Dim> dt_phi_face{};
      slice_tensor(
          make_not_null(&dt_metric_face),
          db::get<::Tags::dt<gr::Tags::SpacetimeMetric<DataVector, Dim>>>(box));
      slice_tensor(make_not_null(&dt_pi_face),
                   db::get<::Tags::dt<gh::Tags::Pi<DataVector, Dim>>>(box));
      slice_tensor(make_not_null(&dt_phi_face),
                   db::get<::Tags::dt<gh::Tags::Phi<DataVector, Dim>>>(box));
      const RateFitResult accel = fit_map_parameter_accelerations(
          metric_face, pi_face, phi_face, dt_metric_face, dt_pi_face,
          dt_phi_face, db::get<Tags::MapParameters>(box).pdot, coords_face,
          ylm_transform, fit_config);
      residual_initial = accel.residual_initial;
      residual_final = accel.residual_final;
      iterations = 1.;
      bool finite = true;
      for (size_t a = 0; a < num_map_parameters; ++a) {
        finite = finite and std::isfinite(gsl::at(accel.pdot, a));
      }
      if (finite) {
        db::mutate<Tags::MapParameters>(
            [&accel, &time, &p_out, &pdot_out, &config_opt,
             &current_worldtube_center](
                const gsl::not_null<MapParameterData*> data) {
              std::array<double, num_map_parameters> acc = accel.pdot;
              if (const double gamma = config_opt->gauge_damping; gamma > 0.) {
                // damped-oscillator gauge fixing of the free parameters,
                // unfolded through the pins (see AdvanceMapParameterOde)
                static constexpr std::array<size_t, 9> free_indices{
                    {0, 1, 2, 3, 7, 8, 9, 10, 11}};
                std::array<double, num_map_parameters> damp{};
                for (const size_t a : free_indices) {
                  gsl::at(damp, a) = -2. * gamma * gsl::at(data->pdot, a) -
                                     gamma * gamma * gsl::at(data->p, a);
                }
                for (size_t i = 0; i < 3; ++i) {
                  gsl::at(damp, 4 + i) =
                      damp[0] * gsl::at(config_opt->center_velocity, i);
                }
                damp[12] = -damp[7] - damp[10];
                for (size_t a = 0; a < num_map_parameters; ++a) {
                  gsl::at(acc, a) += gsl::at(damp, a);
                }
              }
              if (data->valid and time > data->last_fit_time) {
                const double dt = time - data->last_fit_time;
                data->previous_fit_time = data->last_fit_time;
                data->p_previous = data->p;
                for (size_t a = 0; a < num_map_parameters; ++a) {
                  // velocity Verlet: advance with the stored acceleration,
                  // then update the rate with the average of old and new
                  gsl::at(data->p, a) +=
                      dt * gsl::at(data->pdot, a) +
                      0.5 * dt * dt * gsl::at(data->pddot, a);
                  gsl::at(data->pdot, a) +=
                      0.5 * dt * (gsl::at(data->pddot, a) + gsl::at(acc, a));
                }
              }
              data->last_fit_time = time;
              data->worldtube_center_at_last_fit = current_worldtube_center;
              data->pddot = acc;
              data->valid = true;
              p_out = data->p;
              pdot_out = data->pdot;
            },
            make_not_null(&box));
      }
    } else if (config_opt->rate_ode) {
      if (config_opt->fit_center_offset) {
        ERROR(
            "WorldtubeMatcher: RateOde and FitCenterOffset cannot be "
            "combined.");
      }
      // Rate mode: fit pdot linearly from the Pi channel and integrate the
      // first-order ODE dp/dt = pdot by the trapezoid rule; p starts at
      // zero (Schwarzschild) at the first fit.
      const RateFitResult rate =
          fit_map_parameter_rates(metric_face, pi_face, phi_face, coords_face,
                                  ylm_transform, fit_config);
      residual_initial = rate.residual_initial;
      residual_final = rate.residual_final;
      iterations = 1.;
      bool finite = true;
      for (size_t a = 0; a < num_map_parameters; ++a) {
        finite = finite and std::isfinite(gsl::at(rate.pdot, a));
      }
      if (finite) {
        db::mutate<Tags::MapParameters>(
            [&rate, &time, &p_out, &current_worldtube_center](
                const gsl::not_null<MapParameterData*> data) {
              if (data->valid and time > data->last_fit_time) {
                const double dt = time - data->last_fit_time;
                data->previous_fit_time = data->last_fit_time;
                data->p_previous = data->p;
                for (size_t a = 0; a < num_map_parameters; ++a) {
                  gsl::at(data->p, a) +=
                      0.5 * dt *
                      (gsl::at(data->pdot, a) + gsl::at(rate.pdot, a));
                }
              }
              data->last_fit_time = time;
              data->worldtube_center_at_last_fit = current_worldtube_center;
              data->pdot = rate.pdot;
              data->valid = true;
              p_out = data->p;
            },
            make_not_null(&box));
      }
      pdot_out = rate.pdot;
    } else {
      // Strict first-order value mode: determine the instantaneous center and
      // p_(1) algebraically. Hold p_(1) fixed between fits and set its rates
      // to zero; both D_t p_(1) and linear extrapolation of epsilon*p_(1)
      // contribute first at O(epsilon^2). In contrast, q^i is a zeroth-order
      // placement and must be predicted with the separate finite bulk
      // velocity when active, or otherwise dq^i/dt = p_(1)[qdot^i].
      std::array<double, num_map_parameters> p_start{};
      std::array<double, 3> center_offset_start{};
      std::array<double, 3> bulk_velocity_start{};
      if (state.valid) {
        p_start = state.p;
        center_offset_start =
            detail::current_center_offset(fit_config, state, time);
        bulk_velocity_start = state.bulk_velocity;
      }

      // Track velocity for the V_c consistency diagnostic of the exact-frame
      // mode: the worldtube-center motion (the domain map, GB-tracked) plus
      // the hole--worldtube relative motion from the GB dipole when measured.
      if (config_opt->fit_exact_frame and state.valid and
          time > state.last_fit_time and state.worldtube_center_valid) {
        const double elapsed = time - state.last_fit_time;
        for (size_t i = 0; i < 3; ++i) {
          gsl::at(track_velocity, i) =
              (gsl::at(current_worldtube_center, i) -
               gsl::at(state.worldtube_center_at_last_fit, i)) /
                  elapsed +
              gsl::at(state.gb_dipole_velocity, i);
        }
        track_velocity_valid = true;
      }

      FitResult result{};
      if (config_opt->fit_exact_frame) {
        std::array<double, num_map_parameters> theta_start{};
        if (state.exact_frame_valid) {
          theta_start = state.exact_frame_theta;
        } else {
          // Cold start at the physical priors (spec Eq. Z7): the boost from
          // the configured trajectory velocity, s0 = -pin, s_ij = pin
          // delta_ij, sigma = 0. Priors initialize the solve only; they are
          // never constraints.
          const std::array<double, 3> prior_rapidity =
              gh::Solutions::exact_frame::rapidity_from_velocity(
                  config_opt->center_velocity);
          for (size_t i = 0; i < 3; ++i) {
            gsl::at(theta_start, i) = gsl::at(prior_rapidity, i);
          }
          theta_start[3] = -trace_pin;
          theta_start[7] = trace_pin;
          theta_start[10] = trace_pin;
          theta_start[12] = trace_pin;
        }
        std::optional<RadialDerivativeStencil> radial_stencil{};
        if (fit_config.fit_radial_derivative) {
          radial_stencil.emplace();
          radial_stencil->fit_shell = radial_index;
          // The adopted time channel follows the extraction sphere. Its
          // velocity is the l=0 content of the grid-to-inertial mesh velocity;
          // static domains have no mesh-velocity value and therefore zero
          // comoving correction.
          const auto& mesh_velocity =
              db::get<domain::Tags::MeshVelocity<Dim, Frame::Inertial>>(box);
          if (mesh_velocity.has_value()) {
            tnsr::I<DataVector, Dim> mesh_velocity_face{};
            slice_tensor(make_not_null(&mesh_velocity_face), *mesh_velocity);
            radial_stencil->center_velocity =
                detail::worldtube_center(mesh_velocity_face, ylm_transform);
          } else if (state.valid and time > state.last_fit_time) {
            const double elapsed = time - state.last_fit_time;
            for (size_t i = 0; i < 3; ++i) {
              radial_stencil->center_velocity[i] =
                  (current_worldtube_center[i] -
                   state.worldtube_center_at_last_fit[i]) /
                  elapsed;
            }
          }
          const auto& metric_volume =
              db::get<gr::Tags::SpacetimeMetric<DataVector, Dim>>(box);
          const auto& pi_volume = db::get<gh::Tags::Pi<DataVector, Dim>>(box);
          const auto& phi_volume = db::get<gh::Tags::Phi<DataVector, Dim>>(box);
          const auto& gamma2_volume = db::get<gh::Tags::ConstraintGamma2>(box);
          const auto& coords_volume =
              db::get<domain::Tags::Coordinates<Dim, Frame::Inertial>>(box);
          const auto& dt_metric_volume =
              db::get<::Tags::dt<gr::Tags::SpacetimeMetric<DataVector, Dim>>>(
                  box);
          const auto& dt_pi_volume =
              db::get<::Tags::dt<gh::Tags::Pi<DataVector, Dim>>>(box);
          const auto& dt_phi_volume =
              db::get<::Tags::dt<gh::Tags::Phi<DataVector, Dim>>>(box);
          const auto& inverse_jacobian =
              db::get<domain::Tags::InverseJacobian<Dim, Frame::ElementLogical,
                                                    Frame::Inertial>>(box);
          const auto deriv_pi_volume =
              partial_derivative(pi_volume, mesh, inverse_jacobian);
          const auto deriv_phi_volume =
              partial_derivative(phi_volume, mesh, inverse_jacobian);
          for (size_t shell = 0; shell < n_radial; ++shell) {
            tnsr::aa<DataVector, Dim> metric_shell{};
            tnsr::aa<DataVector, Dim> pi_shell{};
            tnsr::iaa<DataVector, Dim> phi_shell{};
            tnsr::aa<DataVector, Dim> dt_metric_shell{};
            tnsr::aa<DataVector, Dim> dt_pi_shell{};
            tnsr::iaa<DataVector, Dim> dt_phi_shell{};
            Scalar<DataVector> gamma2_shell{};
            tnsr::I<DataVector, Dim> coords_shell{};
            slice_tensor_at(make_not_null(&metric_shell), metric_volume, shell);
            slice_tensor_at(make_not_null(&pi_shell), pi_volume, shell);
            slice_tensor_at(make_not_null(&phi_shell), phi_volume, shell);
            slice_tensor_at(make_not_null(&dt_metric_shell), dt_metric_volume,
                            shell);
            slice_tensor_at(make_not_null(&dt_pi_shell), dt_pi_volume, shell);
            slice_tensor_at(make_not_null(&dt_phi_shell), dt_phi_volume, shell);
            tnsr::iaa<DataVector, Dim> deriv_pi_shell{};
            tnsr::ijaa<DataVector, Dim> deriv_phi_shell{};
            slice_tensor_at(make_not_null(&deriv_pi_shell), deriv_pi_volume,
                            shell);
            slice_tensor_at(make_not_null(&deriv_phi_shell), deriv_phi_volume,
                            shell);
            for (size_t a = 0; a < Dim + 1; ++a) {
              for (size_t b = a; b < Dim + 1; ++b) {
                for (size_t i = 0; i < Dim; ++i) {
                  dt_metric_shell.get(a, b) +=
                      radial_stencil->center_velocity[i] *
                      phi_shell.get(i, a, b);
                  dt_pi_shell.get(a, b) += radial_stencil->center_velocity[i] *
                                           deriv_pi_shell.get(i, a, b);
                  for (size_t j = 0; j < Dim; ++j) {
                    dt_phi_shell.get(j, a, b) +=
                        radial_stencil->center_velocity[i] *
                        deriv_phi_shell.get(i, j, a, b);
                  }
                }
              }
            }
            slice_tensor_at(make_not_null(&gamma2_shell), gamma2_volume, shell);
            slice_tensor_at(make_not_null(&coords_shell), coords_volume, shell);
            radial_stencil->metric.push_back(std::move(metric_shell));
            radial_stencil->pi.push_back(std::move(pi_shell));
            radial_stencil->phi.push_back(std::move(phi_shell));
            radial_stencil->dt_metric.push_back(std::move(dt_metric_shell));
            radial_stencil->dt_pi.push_back(std::move(dt_pi_shell));
            radial_stencil->dt_phi.push_back(std::move(dt_phi_shell));
            radial_stencil->gamma2.push_back(std::move(gamma2_shell));
            radial_stencil->coords.push_back(std::move(coords_shell));
          }
        }
        if (config_opt->order_one != OrderOneMode::Off) {
          ASSERT(radial_stencil.has_value(),
                 "The adopted order-one fit requires the element-local "
                 "radial stencil.");
          IteratedOrderZeroOneFitResult iterated = fit_iterated_order_zero_one(
              metric_face, pi_face, phi_face, gamma2_face, coords_face,
              ylm_transform, fit_config, theta_start, center_offset_start,
              *radial_stencil);
          result = std::move(iterated.order_zero);
          order_one_diagnostics = std::move(iterated.order_one);
          order_one_improves_held_out =
              order_one_diagnostics->valid and
              std::isfinite(order_one_diagnostics->minus_residual_initial) and
              std::isfinite(order_one_diagnostics->minus_residual_final) and
              order_one_diagnostics->minus_residual_final <=
                  order_one_diagnostics->minus_residual_initial * (1. + 1.e-12);
        } else {
          result = fit_exact_frame_parameters(
              metric_face, pi_face, phi_face, gamma2_face, coords_face,
              ylm_transform, fit_config, theta_start, center_offset_start,
              radial_stencil);
        }
      } else {
        result = fit_map_parameters(
            metric_face, pi_face, phi_face, gamma2_face, coords_face,
            ylm_transform, fit_config, p_start, center_offset_start,
            config_opt->fit_uplus ? -1.0 : 1.0, trace_pin, bulk_velocity_start);
      }
      residual_initial = result.residual_initial;
      residual_final = result.residual_final;
      iterations = static_cast<double>(result.iterations);
      value_diagnostics = result;

      bool finite = true;
      for (size_t a = 0; a < num_map_parameters; ++a) {
        finite = finite and std::isfinite(gsl::at(result.p, a));
        finite = finite and std::isfinite(gsl::at(result.exact_frame_theta, a));
      }
      for (size_t i = 0; i < 3; ++i) {
        finite = finite and std::isfinite(gsl::at(result.center_offset, i));
        finite = finite and std::isfinite(gsl::at(result.bulk_velocity, i));
        finite = finite and
                 std::isfinite(gsl::at(result.exact_frame_center_velocity, i));
      }
      if (finite) {
        const bool exact_frame_mode = config_opt->fit_exact_frame;
        db::mutate<Tags::MapParameters>(
            [&result, &time, &current_worldtube_center, &exact_frame_mode,
             &order_one_diagnostics](
                const gsl::not_null<MapParameterData*> data) {
              if (data->valid and time > data->last_fit_time) {
                data->previous_fit_time = data->last_fit_time;
                data->p_previous = data->p;
              }
              data->last_fit_time = time;
              data->worldtube_center_at_last_fit = current_worldtube_center;
              data->p = result.p;
              data->pdot.fill(0.);
              data->pddot.fill(0.);
              data->center_offset = result.center_offset;
              data->bulk_velocity = result.bulk_velocity;
              if (exact_frame_mode) {
                data->exact_frame_theta = result.exact_frame_theta;
                data->exact_frame_center_velocity =
                    result.exact_frame_center_velocity;
                data->exact_frame_valid = true;
              }
              data->order_one_valid = order_one_diagnostics.has_value() and
                                      order_one_diagnostics->valid;
              if (data->order_one_valid) {
                data->order_one_rates = order_one_diagnostics->rates;
              } else {
                data->order_one_rates.fill(0.);
              }
              data->valid = true;
            },
            make_not_null(&box));
      }
      p_out = result.p;
      center_offset_out = result.center_offset;
      bulk_velocity_out = result.bulk_velocity;
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
    legend.emplace_back("q_x");
    legend.emplace_back("q_y");
    legend.emplace_back("q_z");
    legend.emplace_back("BulkVelocity_x");
    legend.emplace_back("BulkVelocity_y");
    legend.emplace_back("BulkVelocity_z");
    legend.emplace_back("WorldtubeCenter_x");
    legend.emplace_back("WorldtubeCenter_y");
    legend.emplace_back("WorldtubeCenter_z");
    legend.emplace_back("ResidualInitial");
    legend.emplace_back("ResidualFinal");
    legend.emplace_back("Iterations");
    legend.emplace_back("TracePin");
    // Open-loop hole--worldtube offset from the l = 1 Kretschmann radius and
    // its relative velocity. Together with WorldtubeCenter_i this remains an
    // independent comparison for the fitted inertial velocity.
    for (const std::string& c : {"x", "y", "z"}) {
      legend.push_back("GbDipole_" + c);
    }
    for (const std::string& c : {"x", "y", "z"}) {
      legend.push_back("GbDipoleRelativeVel_" + c);
    }
    if (value_diagnostics.has_value()) {
      legend.emplace_back("BaselineResidual");
      legend.emplace_back("HeldOutMinusBaselineResidual");
      legend.emplace_back("HeldOutMinusResidualFinal");
      legend.emplace_back("WeightedDesignConditionNumber");
      legend.emplace_back("MetricBaselineResidual");
      legend.emplace_back("MetricResidualFinal");
      legend.emplace_back("PhiBaselineResidual");
      legend.emplace_back("PhiResidualFinal");
      legend.emplace_back("BulkMetricResidualInitial");
      legend.emplace_back("BulkMetricResidualFinal");
      for (const std::string& prefix :
           {"Closure_", "TargetRms_", "ResidualRms_", "HeldOutMinusClosure_",
            "HeldOutMinusTargetRms_", "HeldOutMinusResidualRms_"}) {
        for (const std::string& field : {"A", "C", "V"}) {
          for (size_t l = 0; l <= 4; ++l) {
            legend.push_back(prefix + field + "_l" + std::to_string(l));
          }
        }
      }
    }
    if (config_opt->fit_exact_frame) {
      for (const std::string& c : {"x", "y", "z"}) {
        legend.push_back("ExactFrameV_" + c);
      }
      for (const std::string& c : {"x", "y", "z"}) {
        legend.push_back("ExactFrameVc_" + c);
      }
      legend.emplace_back("ExactFrameS0");
      for (const std::string& c : {"x", "y", "z"}) {
        legend.push_back("ExactFrameSigma_" + c);
      }
      for (const std::string& c : {"xx", "xy", "xz", "yy", "yz", "zz"}) {
        legend.push_back("ExactFrameStrain_" + c);
      }
      legend.emplace_back("ExactFrameVcTrackError");
    }
    std::vector<double> row;
    row.reserve(3 * num_map_parameters + 100);
    row.push_back(time);
    for (size_t a = 0; a < num_map_parameters; ++a) {
      row.push_back(gsl::at(p_out, a));
    }
    for (size_t a = 0; a < num_map_parameters; ++a) {
      row.push_back(gsl::at(pdot_out, a));
    }
    for (size_t i = 0; i < 3; ++i) {
      row.push_back(gsl::at(center_offset_out, i));
    }
    for (size_t i = 0; i < 3; ++i) {
      row.push_back(gsl::at(bulk_velocity_out, i));
    }
    for (size_t i = 0; i < 3; ++i) {
      row.push_back(gsl::at(current_worldtube_center, i));
    }
    row.push_back(residual_initial);
    row.push_back(residual_final);
    row.push_back(iterations);
    row.push_back(trace_pin);
    {
      const auto& fresh = db::get<Tags::MapParameters>(box);
      for (size_t i = 0; i < 3; ++i) {
        row.push_back(gsl::at(fresh.gb_dipole, i));
      }
      for (size_t i = 0; i < 3; ++i) {
        row.push_back(gsl::at(fresh.gb_dipole_velocity, i));
      }
    }
    if (value_diagnostics.has_value()) {
      const auto& diagnostic = value_diagnostics.value();
      row.push_back(diagnostic.baseline_residual);
      row.push_back(diagnostic.minus_baseline_residual);
      row.push_back(diagnostic.minus_residual_final);
      row.push_back(diagnostic.condition_number);
      row.push_back(diagnostic.metric_baseline_residual);
      row.push_back(diagnostic.metric_residual_final);
      row.push_back(diagnostic.phi_baseline_residual);
      row.push_back(diagnostic.phi_residual_final);
      row.push_back(diagnostic.bulk_residual_initial);
      row.push_back(diagnostic.bulk_residual_final);
      for (const auto* values :
           {&diagnostic.block_closure, &diagnostic.block_target_rms,
            &diagnostic.block_residual_rms, &diagnostic.block_minus_closure,
            &diagnostic.block_minus_target_rms,
            &diagnostic.block_minus_residual_rms}) {
        for (const double value : *values) {
          row.push_back(value);
        }
      }
    }
    if (config_opt->fit_exact_frame) {
      const auto& diagnostic = value_diagnostics.value();
      for (size_t i = 0; i < 3; ++i) {
        row.push_back(gsl::at(diagnostic.exact_frame_velocity, i));
      }
      for (size_t i = 0; i < 3; ++i) {
        row.push_back(gsl::at(diagnostic.exact_frame_center_velocity, i));
      }
      row.push_back(diagnostic.exact_frame_theta[3]);
      for (size_t i = 0; i < 3; ++i) {
        row.push_back(gsl::at(diagnostic.exact_frame_theta, 4 + i));
      }
      for (size_t pair = 0; pair < 6; ++pair) {
        row.push_back(gsl::at(diagnostic.exact_frame_theta, 7 + pair));
      }
      double vc_track_error = 0.;
      if (track_velocity_valid) {
        for (size_t i = 0; i < 3; ++i) {
          vc_track_error +=
              square(gsl::at(diagnostic.exact_frame_center_velocity, i) -
                     gsl::at(track_velocity, i));
        }
        vc_track_error = std::sqrt(vc_track_error);
      }
      row.push_back(vc_track_error);
    }
    Parallel::threaded_action<
        observers::ThreadedActions::WriteReductionDataRow>(
        writer[0], std::string{"/WorldtubeMatcher"}, std::move(legend),
        std::make_tuple(std::move(row)));

    if (order_one_diagnostics.has_value()) {
      static const std::array<std::string, num_order_one_rates> rate_names{{
          "ClockAcceleration",
          "TimeGradientRate_x",
          "TimeGradientRate_y",
          "TimeGradientRate_z",
          "CenterAcceleration_x",
          "CenterAcceleration_y",
          "CenterAcceleration_z",
          "StrainRate_xx",
          "StrainRate_xy",
          "StrainRate_xz",
          "StrainRate_yy",
          "StrainRate_yz",
          "StrainRate_zz",
          "RotationRate_xy",
          "RotationRate_xz",
          "RotationRate_yz",
      }};
      std::vector<std::string> order_one_legend{"Time"};
      for (const auto& name : rate_names) {
        order_one_legend.push_back(name);
      }
      order_one_legend.emplace_back("DtUPlusEll0ResidualInitial");
      order_one_legend.emplace_back("DtUPlusEll0ResidualFinal");
      order_one_legend.emplace_back("DrDtUPlusEll0ResidualInitial");
      order_one_legend.emplace_back("DrDtUPlusEll0ResidualFinal");
      order_one_legend.emplace_back("HeldOutMinusResidualInitial");
      order_one_legend.emplace_back("HeldOutMinusResidualFinal");
      order_one_legend.emplace_back("WeightedDesignConditionNumber");
      order_one_legend.emplace_back("GatedAlternations");
      order_one_legend.emplace_back("FinalFrameStepNorm");
      order_one_legend.emplace_back("Valid");
      order_one_legend.emplace_back("ImprovesHeldOutMinus");
      const auto& diagnostic = *order_one_diagnostics;
      std::vector<double> order_one_row{};
      order_one_row.reserve(1 + num_order_one_rates + 11);
      order_one_row.push_back(time);
      for (const double rate : diagnostic.rates) {
        order_one_row.push_back(rate);
      }
      order_one_row.push_back(diagnostic.time_residual_initial);
      order_one_row.push_back(diagnostic.time_residual_final);
      order_one_row.push_back(diagnostic.radial_time_residual_initial);
      order_one_row.push_back(diagnostic.radial_time_residual_final);
      order_one_row.push_back(diagnostic.minus_residual_initial);
      order_one_row.push_back(diagnostic.minus_residual_final);
      order_one_row.push_back(diagnostic.condition_number);
      order_one_row.push_back(static_cast<double>(diagnostic.alternations));
      order_one_row.push_back(diagnostic.final_frame_step_norm);
      order_one_row.push_back(diagnostic.valid ? 1. : 0.);
      order_one_row.push_back(order_one_improves_held_out ? 1. : 0.);
      Parallel::threaded_action<
          observers::ThreadedActions::WriteReductionDataRow>(
          writer[0], std::string{"/WorldtubeMatcherOrderOne"},
          std::move(order_one_legend),
          std::make_tuple(std::move(order_one_row)));
    }

    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
}  // namespace gh::Worldtube::Actions
