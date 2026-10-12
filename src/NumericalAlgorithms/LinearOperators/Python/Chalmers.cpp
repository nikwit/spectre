// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "NumericalAlgorithms/LinearOperators/Python/Chalmers.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "DataStructures/DataVector.hpp"
#include "NumericalAlgorithms/LinearOperators/Chalmers.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"

namespace py_bindings {
namespace {
template <size_t Dim>
pybind11::list to_python(
    const std::array<Chalmers::DirectionEstimate, Dim>& estimates) {
  pybind11::list result;
  for (const auto& direction : estimates) {
    pybind11::dict item;
    item["error"] = direction.error;
    item["sigma"] = direction.sigma;
    item["max_log_fit_rms"] = direction.max_log_fit_rms;
    item["insufficient_sequences"] = direction.insufficient_sequences;
    item["ignored_sequences"] = direction.ignored_sequences;
    result.append(item);
  }
  return result;
}

template <size_t Dim>
void bind_impl(pybind11::module& m) {
  m.def("chalmers_p_coarsening_error", &Chalmers::p_coarsening_error<Dim>,
        pybind11::arg("values"), pybind11::arg("mesh"),
        pybind11::arg("absolute_variation_budget") = 0.0,
        pybind11::arg("ignored_highest_modes") = 0,
        pybind11::arg("separate_parity") = false);
  m.def(
      "chalmers_estimate",
      [](const DataVector& values, const Mesh<Dim>& mesh,
         const double stripe_budget, const size_t ignored_highest_modes) {
        return to_python(Chalmers::estimate(values, mesh, stripe_budget,
                                            ignored_highest_modes));
      },
      pybind11::arg("values"), pybind11::arg("mesh"),
      pybind11::arg("absolute_stripe_budget") = 0.0,
      pybind11::arg("ignored_highest_modes") = 0);
  m.def(
      "chalmers_energy_estimate",
      [](const DataVector& values, const Mesh<Dim>& mesh,
         const double variation_budget, const size_t ignored_highest_modes,
         const bool separate_parity) {
        return to_python(Chalmers::estimate_directional_energy(
            values, mesh, variation_budget, ignored_highest_modes,
            separate_parity));
      },
      pybind11::arg("values"), pybind11::arg("mesh"),
      pybind11::arg("absolute_variation_budget") = 0.0,
      pybind11::arg("ignored_highest_modes") = 0,
      pybind11::arg("separate_parity") = false);
}
}  // namespace

void bind_chalmers(pybind11::module& m) {
  m.def(
      "chalmers_parity_estimate",
      [](const DataVector& coefficients, const double budget) {
        return pybind11::dict(to_python(std::array{
            Chalmers::estimate_sequence_parity(coefficients, budget)})[0]);
      },
      pybind11::arg("coefficients"), pybind11::arg("variation_budget") = 0.0);
  bind_impl<1>(m);
  bind_impl<2>(m);
  bind_impl<3>(m);
}
}  // namespace py_bindings
