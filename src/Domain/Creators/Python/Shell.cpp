// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <array>
#include <cstddef>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>

#include "Domain/Creators/DomainCreator.hpp"
#include "Domain/Creators/Shell.hpp"

namespace py = pybind11;

namespace domain::creators::py_bindings {

void bind_shell(py::module& m) {  // NOLINT
  py::class_<Shell, DomainCreator<3>>(m, "Shell")
      .def(py::init<double, double, std::array<size_t, 2>,
                    std::array<size_t, 2>, bool, double>(),
           py::arg("inner_radius"), py::arg("outer_radius"),
           py::arg("initial_refinement"),
           py::arg("initial_number_of_grid_points"),
           py::arg("use_equiangular_map") = true, py::arg("aspect_ratio") = 1.);
}
}  // namespace domain::creators::py_bindings
