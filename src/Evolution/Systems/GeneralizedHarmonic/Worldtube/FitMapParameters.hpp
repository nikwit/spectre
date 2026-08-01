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
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ExtrinsicCurvature.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/Ricci.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeNormalVector.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylElectric.hpp"
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
#include "Utilities/ErrorHandling/Error.hpp"
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
    if (config_opt->stepper_ode and config_opt->uplus_anchor == 0.) {
      // handled by Actions::AdvanceMapParameterOde after the RHS computation
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
    // FitInterval 0 means fit at every invocation, with no time
    // comparison at all: substep times need not be monotone (a stepper
    // may have equal or decreasing stage times), and a high-water-mark
    // gate would silently skip fits and hold a p fitted at a LATER time
    // than the fields it is applied to.
    if (config_opt->fit_interval > 0.) {
      if (config_opt->stepper_ode) {
        if (state.anchor_valid and
            time < state.anchor_time + config_opt->fit_interval - 1.0e-12) {
          return {Parallel::AlgorithmExecution::Continue, std::nullopt};
        }
      } else if (state.valid and
                 time < state.last_fit_time + config_opt->fit_interval -
                            1.0e-12) {
        return {Parallel::AlgorithmExecution::Continue, std::nullopt};
      }
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

    // slice a volume component at the configured radial collocation index
    // (0 = the excision face; radial is the fastest-varying index, matching
    // the Ylm filter's storage convention)
    const size_t radial_index = config_opt->fit_radial_index;
    ASSERT(radial_index < n_radial,
           "FitRadialIndex " << radial_index << " out of range for "
                             << n_radial << " radial points");
    const auto face_slice = [n_radial, n_face,
                             radial_index](const DataVector& volume) {
      DataVector face(n_face);
      for (size_t k = 0; k < n_face; ++k) {
        face[k] = volume[k * n_radial + radial_index];
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
    std::array<double, num_map_parameters> p_out{};
    std::array<double, num_map_parameters> pdot_out{};
    std::array<double, 3> center_offset_out{};
    double residual_initial = 0.;
    double residual_final = 0.;
    double iterations = 0.;

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
        time < state.trace_pin_time + config_opt->trace_pin_interval -
                   1.0e-12) {
      trace_pin = state.trace_pin_value;
    } else if (config_opt->kretschmann_trace_pin) {
      const auto& inv_jacobian = db::get<domain::Tags::InverseJacobian<
          Dim, Frame::ElementLogical, Frame::Inertial>>(box);
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
          pow(48. * square(config_opt->mass) / gb, 1. / 6.) -
          config_opt->mass;
      DataVector coord_r(n_face, 0.);
      for (size_t i = 0; i < 3; ++i) {
        coord_r += square(coords_face.get(i) -
                          gsl::at(config_opt->center, i));
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
        trace_pin = state.trace_pin_time >
                            std::numeric_limits<double>::lowest()
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
            // independent velocity estimate, for comparison with the
            // fitted qdot^i (findings 15w)
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
      const std::array<double, num_map_parameters> zero_rates{};
      std::array<double, num_map_parameters> p_start{};
      const std::array<double, 3> center_offset_start{};
      if (state.anchor_valid) {
        p_start = state.anchor_p;
      }
      const FitResult result = fit_map_parameters(
          metric_face, pi_face, phi_face, gamma2_face, coords_face,
          ylm_transform, *config_opt, p_start, center_offset_start,
          zero_rates, -1.0, trace_pin);
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
          {"qdot0", "b_x", "b_y", "b_z", "v_x", "v_y", "v_z", "s_xx",
           "s_xy", "s_xz", "s_yy", "s_yz", "s_zz"}};
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
          db::get<::Tags::dt<gr::Tags::SpacetimeMetric<DataVector, Dim>>>(
              box));
      slice_tensor(make_not_null(&dt_pi_face),
                   db::get<::Tags::dt<gh::Tags::Pi<DataVector, Dim>>>(box));
      slice_tensor(make_not_null(&dt_phi_face),
                   db::get<::Tags::dt<gh::Tags::Phi<DataVector, Dim>>>(box));
      const RateFitResult accel = fit_map_parameter_accelerations(
          metric_face, pi_face, phi_face, dt_metric_face, dt_pi_face,
          dt_phi_face, db::get<Tags::MapParameters>(box).pdot, coords_face,
          ylm_transform, *config_opt);
      residual_initial = accel.residual_initial;
      residual_final = accel.residual_final;
      iterations = 1.;
      bool finite = true;
      for (size_t a = 0; a < num_map_parameters; ++a) {
        finite = finite and std::isfinite(gsl::at(accel.pdot, a));
      }
      if (finite) {
        db::mutate<Tags::MapParameters>(
            [&accel, &time, &p_out, &pdot_out, &config_opt](
                const gsl::not_null<MapParameterData*> data) {
              std::array<double, num_map_parameters> acc = accel.pdot;
              if (const double gamma = config_opt->gauge_damping;
                  gamma > 0.) {
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
                  gsl::at(data->p, a) += dt * gsl::at(data->pdot, a) +
                                         0.5 * dt * dt *
                                             gsl::at(data->pddot, a);
                  gsl::at(data->pdot, a) +=
                      0.5 * dt *
                      (gsl::at(data->pddot, a) + gsl::at(acc, a));
                }
              }
              data->last_fit_time = time;
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
      const RateFitResult rate = fit_map_parameter_rates(
          metric_face, pi_face, phi_face, coords_face, ylm_transform,
          *config_opt);
      residual_initial = rate.residual_initial;
      residual_final = rate.residual_final;
      iterations = 1.;
      bool finite = true;
      for (size_t a = 0; a < num_map_parameters; ++a) {
        finite = finite and std::isfinite(gsl::at(rate.pdot, a));
      }
      if (finite) {
        db::mutate<Tags::MapParameters>(
            [&rate, &time, &p_out](
                const gsl::not_null<MapParameterData*> data) {
              if (data->valid and time > data->last_fit_time) {
                const double dt = time - data->last_fit_time;
                data->previous_fit_time = data->last_fit_time;
                data->p_previous = data->p;
                for (size_t a = 0; a < num_map_parameters; ++a) {
                  gsl::at(data->p, a) += 0.5 * dt *
                                         (gsl::at(data->pdot, a) +
                                          gsl::at(rate.pdot, a));
                }
              }
              data->last_fit_time = time;
              data->pdot = rate.pdot;
              data->valid = true;
              p_out = data->p;
            },
            make_not_null(&box));
      }
      pdot_out = rate.pdot;
    } else {
      // Value mode: the fit runs with zero drives — feeding the
      // backward-difference rate estimate back into the fit is an unstable
      // loop (the optimal p shifts to compensate a drive error by more than
      // dt, so the estimate grows geometrically; measured x1.33 per fit at
      // FitInterval 0.5). Zero drives in the fit is the offline pass-1,
      // which differs from the converged offline fit by only ~5% of |p|
      // (findings §15b). The rates are computed afterwards, purely as an
      // output for the boundary condition's model Pi.
      const std::array<double, num_map_parameters> zero_rates{};
      std::array<double, num_map_parameters> p_start{};
      std::array<double, 3> center_offset_start{};
      if (state.valid) {
        p_start = state.p;
        center_offset_start = state.center_offset;
      }

      const FitResult result = fit_map_parameters(
          metric_face, pi_face, phi_face, gamma2_face, coords_face,
          ylm_transform, *config_opt, p_start, center_offset_start,
          zero_rates, config_opt->fit_uplus ? -1.0 : 1.0, trace_pin);
      residual_initial = result.residual_initial;
      residual_final = result.residual_final;
      iterations = static_cast<double>(result.iterations);

      bool finite = true;
      for (size_t a = 0; a < num_map_parameters; ++a) {
        finite = finite and std::isfinite(gsl::at(result.p, a));
      }
      for (size_t i = 0; i < 3; ++i) {
        finite = finite and std::isfinite(gsl::at(result.center_offset, i));
      }

      if (finite) {
        // the rates fed to the boundary condition's model Pi come from the
        // evolution equations (dt g = beta Phi - alpha Pi projected on the
        // rate directions), not from backward differences of the fits:
        // smooth, and independent of the fit history (feeding fit
        // differences back was the measured x1.33-per-fit instability)
        const RateFitResult rate = fit_map_parameter_rates(
            metric_face, pi_face, phi_face, coords_face, ylm_transform,
            *config_opt);
        bool rate_finite = true;
        for (size_t a = 0; a < num_map_parameters; ++a) {
          rate_finite =
              rate_finite and std::isfinite(gsl::at(rate.pdot, a));
        }
        if (rate_finite) {
          pdot_out = rate.pdot;
        }
        db::mutate<Tags::MapParameters>(
            [&result, &time, &pdot_out](
                const gsl::not_null<MapParameterData*> data) {
              if (data->valid and time > data->last_fit_time) {
                data->previous_fit_time = data->last_fit_time;
                data->p_previous = data->p;
              }
              data->last_fit_time = time;
              data->p = result.p;
              data->pdot = pdot_out;
              data->center_offset = result.center_offset;
              data->valid = true;
            },
            make_not_null(&box));
      }
      p_out = result.p;
      center_offset_out = result.center_offset;
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
    legend.emplace_back("ResidualInitial");
    legend.emplace_back("ResidualFinal");
    legend.emplace_back("Iterations");
    legend.emplace_back("TracePin");
    // open-loop centre measurement from the l = 1 Kretschmann radius and its
    // backward difference: the independent comparison for a fitted velocity
    for (const std::string& c : {"x", "y", "z"}) {
      legend.push_back("GbDipole_" + c);
    }
    for (const std::string& c : {"x", "y", "z"}) {
      legend.push_back("GbDipoleVel_" + c);
    }
    std::vector<double> row;
    row.reserve(3 * num_map_parameters);
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
    Parallel::threaded_action<
        observers::ThreadedActions::WriteReductionDataRow>(
        writer[0], std::string{"/WorldtubeMatcher"}, std::move(legend),
        std::make_tuple(std::move(row)));

    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
}  // namespace gh::Worldtube::Actions
