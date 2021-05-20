// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/CurvedScalarWave/BoundaryConditions/DirichletAnalytic.hpp"

#include <cstddef>
#include <memory>
#include <pup.h>

#include "Utilities/GenerateInstantiations.hpp"

namespace CurvedScalarWave::BoundaryConditions {
template <size_t Dim>
DirichletAnalytic<Dim>::DirichletAnalytic(CkMigrateMessage* const msg) noexcept
    : BoundaryCondition<Dim>(msg) {}

template <size_t Dim>
std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
DirichletAnalytic<Dim>::get_clone() const noexcept {
  return std::make_unique<DirichletAnalytic>(*this);
}

template <size_t Dim>
void DirichletAnalytic<Dim>::pup(PUP::er& p) {
  BoundaryCondition<Dim>::pup(p);
}

template <size_t Dim>
// NOLINTNEXTLINE
PUP::able::PUP_ID DirichletAnalytic<Dim>::my_PUP_ID = 0;

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATION(r, data) template class DirichletAnalytic<DIM(data)>;

GENERATE_INSTANTIATIONS(INSTANTIATION, (1, 2, 3))

#undef INSTANTIATION
#undef DIM
}  // namespace CurvedScalarWave::BoundaryConditions
