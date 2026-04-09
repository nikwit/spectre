// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <limits>
#include <random>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"

namespace {
template <typename TensorType, typename DataType>
void set_identity_delta(const gsl::not_null<TensorType*> delta,
                        const size_t dimension, const DataType& used_for_size) {
  for (size_t i = 0; i < dimension; ++i) {
    for (size_t j = 0; j < dimension; ++j) {
      delta->get(i, j) =
          make_with_value<DataType>(used_for_size, i == j ? 1.0 : 0.0);
    }
  }
}

template <typename DataType>
void test_lowering_raising_index(const DataType& used_for_size) {
  MAKE_GENERATOR(generator);
  std::uniform_real_distribution<> distribution(-1.0, 1.0);

  auto spatial_delta_Ij =
      make_with_value<tnsr::Ij<DataType, 3, Frame::Grid>>(used_for_size, 0.0);
  auto spatial_delta_iJ =
      make_with_value<tnsr::iJ<DataType, 3, Frame::Grid>>(used_for_size, 0.0);
  set_identity_delta(make_not_null(&spatial_delta_Ij), 3, used_for_size);
  set_identity_delta(make_not_null(&spatial_delta_iJ), 3, used_for_size);

  auto spacetime_delta_Ab =
      make_with_value<tnsr::Ab<DataType, 3, Frame::Grid>>(used_for_size, 0.0);
  auto spacetime_delta_aB =
      make_with_value<tnsr::aB<DataType, 3, Frame::Grid>>(used_for_size, 0.0);
  set_identity_delta(make_not_null(&spacetime_delta_Ab), 4, used_for_size);
  set_identity_delta(make_not_null(&spacetime_delta_aB), 4, used_for_size);

  const auto R_I = make_with_random_values<tnsr::I<DataType, 3, Frame::Grid>>(
      make_not_null(&generator), distribution, used_for_size);
  const auto R_i = make_with_random_values<tnsr::i<DataType, 3, Frame::Grid>>(
      make_not_null(&generator), distribution, used_for_size);
  const auto R_IJ = make_with_random_values<tnsr::IJ<DataType, 3, Frame::Grid>>(
      make_not_null(&generator), distribution, used_for_size);
  const auto R_ij = make_with_random_values<tnsr::ij<DataType, 3, Frame::Grid>>(
      make_not_null(&generator), distribution, used_for_size);
  const auto R_iJ = make_with_random_values<tnsr::iJ<DataType, 3, Frame::Grid>>(
      make_not_null(&generator), distribution, used_for_size);
  const auto R_Ij = make_with_random_values<tnsr::Ij<DataType, 3, Frame::Grid>>(
      make_not_null(&generator), distribution, used_for_size);

  const auto R_A = make_with_random_values<tnsr::A<DataType, 3, Frame::Grid>>(
      make_not_null(&generator), distribution, used_for_size);
  const auto R_a = make_with_random_values<tnsr::a<DataType, 3, Frame::Grid>>(
      make_not_null(&generator), distribution, used_for_size);
  const auto R_AB = make_with_random_values<tnsr::AB<DataType, 3, Frame::Grid>>(
      make_not_null(&generator), distribution, used_for_size);
  const auto R_ab = make_with_random_values<tnsr::ab<DataType, 3, Frame::Grid>>(
      make_not_null(&generator), distribution, used_for_size);
  const auto R_aB = make_with_random_values<tnsr::aB<DataType, 3, Frame::Grid>>(
      make_not_null(&generator), distribution, used_for_size);
  const auto R_Ab = make_with_random_values<tnsr::Ab<DataType, 3, Frame::Grid>>(
      make_not_null(&generator), distribution, used_for_size);

  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::J>(spatial_delta_iJ(ti::i, ti::J) * R_I(ti::I))),
      (tenex::evaluate<ti::J>(
          tenex::spatial_kronecker_delta<3, Frame::Grid>(ti::i, ti::J) *
          R_I(ti::I))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::I>(spatial_delta_Ij(ti::I, ti::j) * R_I(ti::J))),
      (tenex::evaluate<ti::I>(
          tenex::spatial_kronecker_delta<3, Frame::Grid>(ti::I, ti::j) *
          R_I(ti::J))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::j>(spatial_delta_Ij(ti::I, ti::j) * R_i(ti::i))),
      (tenex::evaluate<ti::j>(
          tenex::spatial_kronecker_delta<3, Frame::Grid>(ti::I, ti::j) *
          R_i(ti::i))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::i>(spatial_delta_iJ(ti::i, ti::J) * R_i(ti::j))),
      (tenex::evaluate<ti::i>(
          tenex::spatial_kronecker_delta<3, Frame::Grid>(ti::i, ti::J) *
          R_i(ti::j))));

  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::J, ti::K>(spatial_delta_iJ(ti::i, ti::J) *
                                     R_IJ(ti::I, ti::K))),
      (tenex::evaluate<ti::J, ti::K>(
          tenex::spatial_kronecker_delta<3, Frame::Grid>(ti::i, ti::J) *
          R_IJ(ti::I, ti::K))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::I, ti::K>(spatial_delta_Ij(ti::I, ti::j) *
                                     R_IJ(ti::K, ti::J))),
      (tenex::evaluate<ti::I, ti::K>(
          tenex::spatial_kronecker_delta<3, Frame::Grid>(ti::I, ti::j) *
          R_IJ(ti::K, ti::J))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::j, ti::k>(spatial_delta_Ij(ti::I, ti::j) *
                                     R_ij(ti::i, ti::k))),
      (tenex::evaluate<ti::j, ti::k>(
          tenex::spatial_kronecker_delta<3, Frame::Grid>(ti::I, ti::j) *
          R_ij(ti::i, ti::k))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::i, ti::k>(spatial_delta_iJ(ti::i, ti::J) *
                                     R_ij(ti::k, ti::j))),
      (tenex::evaluate<ti::i, ti::k>(
          tenex::spatial_kronecker_delta<3, Frame::Grid>(ti::i, ti::J) *
          R_ij(ti::k, ti::j))));

  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::j, ti::K>(spatial_delta_Ij(ti::I, ti::j) *
                                     R_iJ(ti::i, ti::K))),
      (tenex::evaluate<ti::j, ti::K>(
          tenex::spatial_kronecker_delta<3, Frame::Grid>(ti::I, ti::j) *
          R_iJ(ti::i, ti::K))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::J, ti::k>(spatial_delta_iJ(ti::i, ti::J) *
                                     R_iJ(ti::k, ti::I))),
      (tenex::evaluate<ti::J, ti::k>(
          tenex::spatial_kronecker_delta<3, Frame::Grid>(ti::i, ti::J) *
          R_iJ(ti::k, ti::I))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::i, ti::K>(spatial_delta_iJ(ti::i, ti::J) *
                                     R_iJ(ti::j, ti::K))),
      (tenex::evaluate<ti::i, ti::K>(
          tenex::spatial_kronecker_delta<3, Frame::Grid>(ti::i, ti::J) *
          R_iJ(ti::j, ti::K))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::I, ti::k>(spatial_delta_Ij(ti::I, ti::j) *
                                     R_iJ(ti::k, ti::J))),
      (tenex::evaluate<ti::I, ti::k>(
          tenex::spatial_kronecker_delta<3, Frame::Grid>(ti::I, ti::j) *
          R_iJ(ti::k, ti::J))));

  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::J, ti::k>(spatial_delta_iJ(ti::i, ti::J) *
                                     R_Ij(ti::I, ti::k))),
      (tenex::evaluate<ti::J, ti::k>(
          tenex::spatial_kronecker_delta<3, Frame::Grid>(ti::i, ti::J) *
          R_Ij(ti::I, ti::k))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::j, ti::K>(spatial_delta_Ij(ti::I, ti::j) *
                                     R_Ij(ti::K, ti::i))),
      (tenex::evaluate<ti::j, ti::K>(
          tenex::spatial_kronecker_delta<3, Frame::Grid>(ti::I, ti::j) *
          R_Ij(ti::K, ti::i))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::i, ti::K>(spatial_delta_iJ(ti::i, ti::J) *
                                     R_Ij(ti::K, ti::j))),
      (tenex::evaluate<ti::i, ti::K>(
          tenex::spatial_kronecker_delta<3, Frame::Grid>(ti::i, ti::J) *
          R_Ij(ti::K, ti::j))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::I, ti::k>(spatial_delta_Ij(ti::I, ti::j) *
                                     R_Ij(ti::J, ti::k))),
      (tenex::evaluate<ti::I, ti::k>(
          tenex::spatial_kronecker_delta<3, Frame::Grid>(ti::I, ti::j) *
          R_Ij(ti::J, ti::k))));

  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::B>(spacetime_delta_aB(ti::a, ti::B) * R_A(ti::A))),
      (tenex::evaluate<ti::B>(
          tenex::spacetime_kronecker_delta<3, Frame::Grid>(ti::a, ti::B) *
          R_A(ti::A))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::A>(spacetime_delta_Ab(ti::A, ti::b) * R_A(ti::B))),
      (tenex::evaluate<ti::A>(
          tenex::spacetime_kronecker_delta<3, Frame::Grid>(ti::A, ti::b) *
          R_A(ti::B))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::b>(spacetime_delta_Ab(ti::A, ti::b) * R_a(ti::a))),
      (tenex::evaluate<ti::b>(
          tenex::spacetime_kronecker_delta<3, Frame::Grid>(ti::A, ti::b) *
          R_a(ti::a))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::a>(spacetime_delta_aB(ti::a, ti::B) * R_a(ti::b))),
      (tenex::evaluate<ti::a>(
          tenex::spacetime_kronecker_delta<3, Frame::Grid>(ti::a, ti::B) *
          R_a(ti::b))));

  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::B, ti::C>(spacetime_delta_aB(ti::a, ti::B) *
                                     R_AB(ti::A, ti::C))),
      (tenex::evaluate<ti::B, ti::C>(
          tenex::spacetime_kronecker_delta<3, Frame::Grid>(ti::a, ti::B) *
          R_AB(ti::A, ti::C))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::A, ti::C>(spacetime_delta_Ab(ti::A, ti::b) *
                                     R_AB(ti::C, ti::B))),
      (tenex::evaluate<ti::A, ti::C>(
          tenex::spacetime_kronecker_delta<3, Frame::Grid>(ti::A, ti::b) *
          R_AB(ti::C, ti::B))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::b, ti::c>(spacetime_delta_Ab(ti::A, ti::b) *
                                     R_ab(ti::a, ti::c))),
      (tenex::evaluate<ti::b, ti::c>(
          tenex::spacetime_kronecker_delta<3, Frame::Grid>(ti::A, ti::b) *
          R_ab(ti::a, ti::c))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::a, ti::c>(spacetime_delta_aB(ti::a, ti::B) *
                                     R_ab(ti::c, ti::b))),
      (tenex::evaluate<ti::a, ti::c>(
          tenex::spacetime_kronecker_delta<3, Frame::Grid>(ti::a, ti::B) *
          R_ab(ti::c, ti::b))));

  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::b, ti::C>(spacetime_delta_Ab(ti::A, ti::b) *
                                     R_aB(ti::a, ti::C))),
      (tenex::evaluate<ti::b, ti::C>(
          tenex::spacetime_kronecker_delta<3, Frame::Grid>(ti::A, ti::b) *
          R_aB(ti::a, ti::C))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::B, ti::c>(spacetime_delta_aB(ti::a, ti::B) *
                                     R_aB(ti::c, ti::A))),
      (tenex::evaluate<ti::B, ti::c>(
          tenex::spacetime_kronecker_delta<3, Frame::Grid>(ti::a, ti::B) *
          R_aB(ti::c, ti::A))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::a, ti::C>(spacetime_delta_aB(ti::a, ti::B) *
                                     R_aB(ti::b, ti::C))),
      (tenex::evaluate<ti::a, ti::C>(
          tenex::spacetime_kronecker_delta<3, Frame::Grid>(ti::a, ti::B) *
          R_aB(ti::b, ti::C))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::A, ti::c>(spacetime_delta_Ab(ti::A, ti::b) *
                                     R_aB(ti::c, ti::B))),
      (tenex::evaluate<ti::A, ti::c>(
          tenex::spacetime_kronecker_delta<3, Frame::Grid>(ti::A, ti::b) *
          R_aB(ti::c, ti::B))));

  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::B, ti::c>(spacetime_delta_aB(ti::a, ti::B) *
                                     R_Ab(ti::A, ti::c))),
      (tenex::evaluate<ti::B, ti::c>(
          tenex::spacetime_kronecker_delta<3, Frame::Grid>(ti::a, ti::B) *
          R_Ab(ti::A, ti::c))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::b, ti::C>(spacetime_delta_Ab(ti::A, ti::b) *
                                     R_Ab(ti::C, ti::a))),
      (tenex::evaluate<ti::b, ti::C>(
          tenex::spacetime_kronecker_delta<3, Frame::Grid>(ti::A, ti::b) *
          R_Ab(ti::C, ti::a))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::a, ti::C>(spacetime_delta_aB(ti::a, ti::B) *
                                     R_Ab(ti::C, ti::b))),
      (tenex::evaluate<ti::a, ti::C>(
          tenex::spacetime_kronecker_delta<3, Frame::Grid>(ti::a, ti::B) *
          R_Ab(ti::C, ti::b))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::A, ti::c>(spacetime_delta_Ab(ti::A, ti::b) *
                                     R_Ab(ti::B, ti::c))),
      (tenex::evaluate<ti::A, ti::c>(
          tenex::spacetime_kronecker_delta<3, Frame::Grid>(ti::A, ti::b) *
          R_Ab(ti::B, ti::c))));
}

template <typename DataType>
void test_trace(const DataType& used_for_size) {
  MAKE_GENERATOR(generator);
  std::uniform_real_distribution<> distribution(-1.0, 1.0);

  auto spatial_delta_Ij =
      make_with_value<tnsr::Ij<DataType, 3, Frame::Grid>>(used_for_size, 0.0);
  auto spatial_delta_iJ =
      make_with_value<tnsr::iJ<DataType, 3, Frame::Grid>>(used_for_size, 0.0);
  set_identity_delta(make_not_null(&spatial_delta_Ij), 3, used_for_size);
  set_identity_delta(make_not_null(&spatial_delta_iJ), 3, used_for_size);

  auto spacetime_delta_Ab =
      make_with_value<tnsr::Ab<DataType, 3, Frame::Grid>>(used_for_size, 0.0);
  auto spacetime_delta_aB =
      make_with_value<tnsr::aB<DataType, 3, Frame::Grid>>(used_for_size, 0.0);
  set_identity_delta(make_not_null(&spacetime_delta_Ab), 4, used_for_size);
  set_identity_delta(make_not_null(&spacetime_delta_aB), 4, used_for_size);

  const auto R_Ij = make_with_random_values<tnsr::Ij<DataType, 3, Frame::Grid>>(
      make_not_null(&generator), distribution, used_for_size);
  const auto R_iJ = make_with_random_values<tnsr::iJ<DataType, 3, Frame::Grid>>(
      make_not_null(&generator), distribution, used_for_size);
  const auto R_Ab = make_with_random_values<tnsr::Ab<DataType, 3, Frame::Grid>>(
      make_not_null(&generator), distribution, used_for_size);
  const auto R_aB = make_with_random_values<tnsr::aB<DataType, 3, Frame::Grid>>(
      make_not_null(&generator), distribution, used_for_size);

  CHECK_ITERABLE_APPROX(
      (tenex::evaluate(spatial_delta_iJ(ti::i, ti::J) * R_Ij(ti::I, ti::j))),
      (tenex::evaluate(
          tenex::spatial_kronecker_delta<3, Frame::Grid>(ti::i, ti::J) *
          R_Ij(ti::I, ti::j))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate(spatial_delta_Ij(ti::I, ti::j) * R_iJ(ti::i, ti::J))),
      (tenex::evaluate(
          tenex::spatial_kronecker_delta<3, Frame::Grid>(ti::I, ti::j) *
          R_iJ(ti::i, ti::J))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate(spacetime_delta_aB(ti::a, ti::B) * R_Ab(ti::A, ti::b))),
      (tenex::evaluate(
          tenex::spacetime_kronecker_delta<3, Frame::Grid>(ti::a, ti::B) *
          R_Ab(ti::A, ti::b))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate(spacetime_delta_Ab(ti::A, ti::b) * R_aB(ti::a, ti::B))),
      (tenex::evaluate(
          tenex::spacetime_kronecker_delta<3, Frame::Grid>(ti::A, ti::b) *
          R_aB(ti::a, ti::B))));
}

template <typename DataType>
void test_multiplying_from_right(const DataType& used_for_size) {
  MAKE_GENERATOR(generator);
  std::uniform_real_distribution<> distribution(-1.0, 1.0);

  auto spatial_delta_iJ =
      make_with_value<tnsr::iJ<DataType, 3, Frame::Grid>>(used_for_size, 0.0);
  auto spatial_delta_Ij =
      make_with_value<tnsr::Ij<DataType, 3, Frame::Grid>>(used_for_size, 0.0);
  set_identity_delta(make_not_null(&spatial_delta_iJ), 3, used_for_size);
  set_identity_delta(make_not_null(&spatial_delta_Ij), 3, used_for_size);

  auto spacetime_delta_aB =
      make_with_value<tnsr::aB<DataType, 3, Frame::Grid>>(used_for_size, 0.0);
  auto spacetime_delta_Ab =
      make_with_value<tnsr::Ab<DataType, 3, Frame::Grid>>(used_for_size, 0.0);
  set_identity_delta(make_not_null(&spacetime_delta_aB), 4, used_for_size);
  set_identity_delta(make_not_null(&spacetime_delta_Ab), 4, used_for_size);

  const auto R_I = make_with_random_values<tnsr::I<DataType, 3, Frame::Grid>>(
      make_not_null(&generator), distribution, used_for_size);
  const auto S_i = make_with_random_values<tnsr::i<DataType, 3, Frame::Grid>>(
      make_not_null(&generator), distribution, used_for_size);
  const auto T_iJ = make_with_random_values<tnsr::iJ<DataType, 3, Frame::Grid>>(
      make_not_null(&generator), distribution, used_for_size);

  const auto R_A = make_with_random_values<tnsr::A<DataType, 3, Frame::Grid>>(
      make_not_null(&generator), distribution, used_for_size);
  const auto S_a = make_with_random_values<tnsr::a<DataType, 3, Frame::Grid>>(
      make_not_null(&generator), distribution, used_for_size);
  const auto T_aB = make_with_random_values<tnsr::aB<DataType, 3, Frame::Grid>>(
      make_not_null(&generator), distribution, used_for_size);

  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::J>(R_I(ti::I) * spatial_delta_iJ(ti::i, ti::J))),
      (tenex::evaluate<ti::J>(
          R_I(ti::I) *
          tenex::spatial_kronecker_delta<3, Frame::Grid>(ti::i, ti::J))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate(S_i(ti::j) * R_I(ti::I) *
                       spatial_delta_iJ(ti::i, ti::J))),
      (tenex::evaluate(
          S_i(ti::j) * R_I(ti::I) *
          tenex::spatial_kronecker_delta<3, Frame::Grid>(ti::i, ti::J))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::j, ti::K>(T_iJ(ti::i, ti::K) *
                                     spatial_delta_Ij(ti::I, ti::j))),
      (tenex::evaluate<ti::j, ti::K>(
          T_iJ(ti::i, ti::K) *
          tenex::spatial_kronecker_delta<3, Frame::Grid>(ti::I, ti::j))));

  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::B>(R_A(ti::A) * spacetime_delta_aB(ti::a, ti::B))),
      (tenex::evaluate<ti::B>(
          R_A(ti::A) *
          tenex::spacetime_kronecker_delta<3, Frame::Grid>(ti::a, ti::B))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate(S_a(ti::b) * R_A(ti::A) *
                       spacetime_delta_aB(ti::a, ti::B))),
      (tenex::evaluate(
          S_a(ti::b) * R_A(ti::A) *
          tenex::spacetime_kronecker_delta<3, Frame::Grid>(ti::a, ti::B))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::b, ti::C>(T_aB(ti::a, ti::C) *
                                     spacetime_delta_Ab(ti::A, ti::b))),
      (tenex::evaluate<ti::b, ti::C>(
          T_aB(ti::a, ti::C) *
          tenex::spacetime_kronecker_delta<3, Frame::Grid>(ti::A, ti::b))));
}

void test_multifactor_delta_positions() {
  MAKE_GENERATOR(generator);
  std::uniform_real_distribution<> distribution(-1.0, 1.0);
  const double used_for_size = std::numeric_limits<double>::signaling_NaN();

  auto spatial_delta_iJ =
      make_with_value<tnsr::iJ<double, 3, Frame::Grid>>(used_for_size, 0.0);
  set_identity_delta(make_not_null(&spatial_delta_iJ), 3, used_for_size);

  const auto R_i = make_with_random_values<tnsr::i<double, 3, Frame::Grid>>(
      make_not_null(&generator), distribution, used_for_size);
  const auto S_I = make_with_random_values<tnsr::I<double, 3, Frame::Grid>>(
      make_not_null(&generator), distribution, used_for_size);
  const auto T_I = make_with_random_values<tnsr::I<double, 3, Frame::Grid>>(
      make_not_null(&generator), distribution, used_for_size);
  const auto U_i = make_with_random_values<tnsr::i<double, 3, Frame::Grid>>(
      make_not_null(&generator), distribution, used_for_size);

  CHECK_ITERABLE_APPROX(
      (tenex::evaluate(spatial_delta_iJ(ti::i, ti::J) * S_I(ti::I) *
                       R_i(ti::j) * T_I(ti::K) * U_i(ti::k))),
      (tenex::evaluate(
          tenex::spatial_kronecker_delta<3, Frame::Grid>(ti::i, ti::J) *
          S_I(ti::I) * R_i(ti::j) * T_I(ti::K) * U_i(ti::k))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate(S_I(ti::I) * spatial_delta_iJ(ti::i, ti::J) *
                       R_i(ti::j) * T_I(ti::K) * U_i(ti::k))),
      (tenex::evaluate(
          S_I(ti::I) *
          tenex::spatial_kronecker_delta<3, Frame::Grid>(ti::i, ti::J) *
          R_i(ti::j) * T_I(ti::K) * U_i(ti::k))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate(S_I(ti::I) * R_i(ti::j) *
                       spatial_delta_iJ(ti::i, ti::J) * T_I(ti::K) *
                       U_i(ti::k))),
      (tenex::evaluate(
          S_I(ti::I) * R_i(ti::j) *
          tenex::spatial_kronecker_delta<3, Frame::Grid>(ti::i, ti::J) *
          T_I(ti::K) * U_i(ti::k))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate(S_I(ti::I) * R_i(ti::j) * T_I(ti::K) *
                       spatial_delta_iJ(ti::i, ti::J) * U_i(ti::k))),
      (tenex::evaluate(
          S_I(ti::I) * R_i(ti::j) * T_I(ti::K) *
          tenex::spatial_kronecker_delta<3, Frame::Grid>(ti::i, ti::J) *
          U_i(ti::k))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate(S_I(ti::I) * R_i(ti::j) * T_I(ti::K) * U_i(ti::k) *
                       spatial_delta_iJ(ti::i, ti::J))),
      (tenex::evaluate(
          S_I(ti::I) * R_i(ti::j) * T_I(ti::K) * U_i(ti::k) *
          tenex::spatial_kronecker_delta<3, Frame::Grid>(ti::i, ti::J))));
}

template <typename DataType>
void test_two_deltas_with_tensor(const DataType& used_for_size) {
  MAKE_GENERATOR(generator);
  std::uniform_real_distribution<> distribution(-1.0, 1.0);

  auto spacetime_delta_Ab_1 =
      make_with_value<tnsr::Ab<DataType, 3, Frame::Inertial>>(used_for_size,
                                                              0.0);
  auto spacetime_delta_Ab_2 =
      make_with_value<tnsr::Ab<DataType, 3, Frame::Inertial>>(used_for_size,
                                                              0.0);
  set_identity_delta(make_not_null(&spacetime_delta_Ab_1), 4, used_for_size);
  set_identity_delta(make_not_null(&spacetime_delta_Ab_2), 4, used_for_size);

  const auto spacetime_vector =
      make_with_random_values<tnsr::A<DataType, 3, Frame::Inertial>>(
          make_not_null(&generator), distribution, used_for_size);

  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::A>(spacetime_delta_Ab_1(ti::A, ti::b) *
                              spacetime_delta_Ab_2(ti::B, ti::c) *
                              spacetime_vector(ti::C))),
      (tenex::evaluate<ti::A>(
          tenex::spacetime_kronecker_delta<3, Frame::Inertial>(ti::A, ti::b) *
          tenex::spacetime_kronecker_delta<3, Frame::Inertial>(ti::B, ti::c) *
          spacetime_vector(ti::C))));
}

void test_addition_with_delta() {
  MAKE_GENERATOR(generator);
  std::uniform_real_distribution<> distribution(-1.0, 1.0);
  const double used_for_size = std::numeric_limits<double>::signaling_NaN();

  auto spatial_delta_iJ =
      make_with_value<tnsr::iJ<double, 3, Frame::Grid>>(used_for_size, 0.0);
  set_identity_delta(make_not_null(&spatial_delta_iJ), 3, used_for_size);

  const auto R_iJ = make_with_random_values<tnsr::iJ<double, 3, Frame::Grid>>(
      make_not_null(&generator), distribution, used_for_size);

  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::i, ti::J>(R_iJ(ti::i, ti::J) +
                                     spatial_delta_iJ(ti::i, ti::J))),
      (tenex::evaluate<ti::i, ti::J>(
          R_iJ(ti::i, ti::J) +
          tenex::spatial_kronecker_delta<3, Frame::Grid>(ti::i, ti::J))));
  CHECK_ITERABLE_APPROX(
      (tenex::evaluate<ti::i, ti::J>(spatial_delta_iJ(ti::i, ti::J) +
                                     R_iJ(ti::i, ti::J))),
      (tenex::evaluate<ti::i, ti::J>(
          tenex::spatial_kronecker_delta<3, Frame::Grid>(ti::i, ti::J) +
          R_iJ(ti::i, ti::J))));
}
}  // namespace

SPECTRE_TEST_CASE("Unit.DataStructures.Tensor.Expressions.KroneckerDelta",
                  "[DataStructures][Unit]") {
  test_lowering_raising_index(std::numeric_limits<double>::signaling_NaN());
  test_trace(std::numeric_limits<double>::signaling_NaN());
  test_multiplying_from_right(std::numeric_limits<double>::signaling_NaN());
  test_two_deltas_with_tensor(std::numeric_limits<double>::signaling_NaN());
  test_multifactor_delta_positions();
  test_addition_with_delta();

  test_lowering_raising_index(DataVector(5));
  test_trace(DataVector(5));
  test_multiplying_from_right(DataVector(5));
  test_two_deltas_with_tensor(DataVector(5));
}
