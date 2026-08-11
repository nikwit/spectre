
// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/DynamicBuffer.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "Utilities/Gsl.hpp"

/// \cond
namespace CurvedScalarWave::Worldtube::detail {

// Implementation stages of the auto-generated second-order puncture field,
// split into separate translation units to keep per-function compile time
// bounded. The stages form a pipeline over a shared buffer of temporaries
// and must be called in order.

using Order2Vars =
    Variables<tmpl::list<CurvedScalarWave::Tags::Psi,
                         ::Tags::dt<CurvedScalarWave::Tags::Psi>,
                         ::Tags::deriv<CurvedScalarWave::Tags::Psi,
                                       tmpl::size_t<3>, Frame::Inertial>>>;

inline constexpr size_t order2_n_doubles = 2920;

void puncture_field_2_part_0(const std::array<double, order2_n_doubles>& d,
                             const DataVector& Dx, const DataVector& Dy,
                             const DataVector& z,
                             DynamicBuffer<DataVector>& temps);
void puncture_field_2_part_1(const std::array<double, order2_n_doubles>& d,
                             const DataVector& Dx, const DataVector& Dy,
                             DynamicBuffer<DataVector>& temps);
void puncture_field_2_part_2(const std::array<double, order2_n_doubles>& d,
                             const DataVector& Dx, const DataVector& Dy,
                             DynamicBuffer<DataVector>& temps);
void puncture_field_2_part_3(const std::array<double, order2_n_doubles>& d,
                             const DataVector& Dx, const DataVector& Dy,
                             DynamicBuffer<DataVector>& temps);
void puncture_field_2_part_4(const std::array<double, order2_n_doubles>& d,
                             const DataVector& Dx, const DataVector& Dy,
                             DynamicBuffer<DataVector>& temps);
void puncture_field_2_part_5(const std::array<double, order2_n_doubles>& d,
                             const DataVector& Dx, const DataVector& Dy,
                             DynamicBuffer<DataVector>& temps,
                             const gsl::not_null<Order2Vars*> result);
void puncture_field_2_part_6(const std::array<double, order2_n_doubles>& d,
                             const DataVector& Dx, const DataVector& Dy,
                             DynamicBuffer<DataVector>& temps);
void puncture_field_2_part_7(const std::array<double, order2_n_doubles>& d,
                             const DataVector& Dx, const DataVector& Dy,
                             DynamicBuffer<DataVector>& temps);
void puncture_field_2_part_8(const std::array<double, order2_n_doubles>& d,
                             const DataVector& Dx, const DataVector& Dy,
                             DynamicBuffer<DataVector>& temps);
void puncture_field_2_part_9(const std::array<double, order2_n_doubles>& d,
                             const DataVector& Dx, const DataVector& Dy,
                             DynamicBuffer<DataVector>& temps);
void puncture_field_2_part_10(const std::array<double, order2_n_doubles>& d,
                              const DataVector& Dx, const DataVector& Dy,
                              DynamicBuffer<DataVector>& temps,
                              const gsl::not_null<Order2Vars*> result);
void puncture_field_2_part_11(const std::array<double, order2_n_doubles>& d,
                              const DataVector& Dx, const DataVector& Dy,
                              const DataVector& z,
                              DynamicBuffer<DataVector>& temps,
                              const gsl::not_null<Order2Vars*> result);
}  // namespace CurvedScalarWave::Worldtube::detail
/// \endcond
