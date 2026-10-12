// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/DampedHarmonic.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/SpecAmr.hpp"

namespace gh::spec_amr {
namespace {
pybind11::dict spectrum_dict(const Spectrum& spectrum) {
  pybind11::dict result{};
  result["error"] = spectrum.error;
  result["next_error"] = spectrum.next_error;
  result["convergence"] = spectrum.convergence;
  result["pile_up"] = spectrum.pile_up;
  result["degenerate"] = spectrum.degenerate;
  return result;
}
}  // namespace
void bind_spec_amr(pybind11::module_& module) {
  namespace py = pybind11;
  module.def("spec_amr_spectrum", [](const DataVector& power) {
    return spectrum_dict(diagnose_power(power));
  });
  module.def(
      "spec_amr_damped_harmonic_diagnostics",
      [](const tnsr::aa<DataVector, 3>& metric,
         const tnsr::aa<DataVector, 3>& pi, const tnsr::iaa<DataVector, 3>& phi,
         const Mesh<3>& mesh,
         const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                               Frame::Inertial>& inverse,
         const tnsr::I<DataVector, 3>& coordinates, const double time,
         const Scalar<DataVector>& gamma2, const double width,
         const std::array<double, 3>& amplitudes,
         const std::array<int, 3>& exponents) {
        const gauges::DampedHarmonic gauge{width, amplitudes, exponents};
        const auto diagnostics = diagnose(metric, pi, phi, mesh, inverse,
                                          coordinates, time, gauge, gamma2);
        py::dict result{};
        py::list spectra{};
        for (const auto& group : diagnostics.spectra) {
          py::list directions{};
          for (const auto& spectrum : group) {
            directions.append(spectrum_dict(spectrum));
          }
          spectra.append(directions);
        }
        result["spectra"] = spectra;
        result["constraints"] = diagnostics.constraints;
        return result;
      },
      py::arg("metric"), py::arg("pi"), py::arg("phi"), py::arg("mesh"),
      py::arg("inverse_jacobian"), py::arg("coordinates"), py::arg("time"),
      py::arg("gamma2"), py::arg("width"), py::arg("amplitudes"),
      py::arg("exponents"));
}
}  // namespace gh::spec_amr
