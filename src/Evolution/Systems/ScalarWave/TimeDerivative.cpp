// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/ScalarWave/TimeDerivative.hpp"

#include "DataStructures/DataVector.hpp"
#include "Utilities/GenerateInstantiations.hpp"

namespace ScalarWave {
template class TimeDerivative<1>;
template class TimeDerivative<2>;
template class TimeDerivative<3>;

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)
#define DTYPE(data) BOOST_PP_TUPLE_ELEM(1, data)

#define INSTANTIATION(_, data)                                                \
  template evolution::dg::TimeDerivativeDecisions<DIM(data)>                  \
  TimeDerivative<DIM(data)>::apply(                                           \
      const gsl::not_null<Scalar<DTYPE(data)>*>,                              \
      const gsl::not_null<Scalar<DTYPE(data)>*>,                              \
      const gsl::not_null<tnsr::i<DTYPE(data), DIM(data), Frame::Inertial>*>, \
      const gsl::not_null<Scalar<DTYPE(data)>*>,                              \
      const tnsr::i<DTYPE(data), DIM(data), Frame::Inertial>&,                \
      const tnsr::i<DTYPE(data), DIM(data), Frame::Inertial>&,                \
      const tnsr::ij<DTYPE(data), DIM(data), Frame::Inertial>&,               \
      const Scalar<DTYPE(data)>&,                                             \
      const tnsr::i<DTYPE(data), DIM(data), Frame::Inertial>&,                \
      const Scalar<DTYPE(data)>&);

GENERATE_INSTANTIATIONS(INSTANTIATION, (1, 2, 3), (double, DataVector))

#undef INSTANTIATION
#undef DIM
#undef DTYPE
}  // namespace ScalarWave
