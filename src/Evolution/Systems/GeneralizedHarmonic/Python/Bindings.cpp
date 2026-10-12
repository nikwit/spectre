// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <pybind11/pybind11.h>

#include "Evolution/Systems/GeneralizedHarmonic/Python/FilledSpherePowerMonitor.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Python/SphericalShellPowerMonitor.hpp"
#include "Utilities/ErrorHandling/SegfaultHandler.hpp"

namespace py = pybind11;

namespace gh::spec_amr {
void bind_spec_amr(py::module_& module);
}

PYBIND11_MODULE(_Pybindings, m) {  // NOLINT
  enable_segfault_handler();
  py::module_::import("spectre.DataStructures");
  py::module_::import("spectre.DataStructures.Tensor");
  py::module_::import("spectre.Domain");
  py::module_::import("spectre.Spectral");

  gh::power_monitor::py_bindings::bind_filled_sphere_power_monitor(m);
  gh::power_monitor::py_bindings::bind_spherical_shell_power_monitor(m);
  gh::spec_amr::bind_spec_amr(m);
}
