// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "DataStructures/DynamicBuffer.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"

void check_results(const DynamicBuffer& dynamic_buffer,
                   const std::vector<DataVector>& expected) {
  CHECK(dynamic_buffer.size() == expected.size());

  for (size_t i = 0; i < dynamic_buffer.size(); ++i) {
    CHECK_ITERABLE_APPROX(dynamic_buffer[i], expected[i]);
    CHECK_ITERABLE_APPROX(dynamic_buffer.at(i), expected[i]);
  }

  size_t i = 0;
  for (auto it = dynamic_buffer.begin(); it != dynamic_buffer.end();
       ++it, ++i) {
    CHECK_ITERABLE_APPROX(*it, expected[i]);
  }

  i = 0;
  for (const auto& val : dynamic_buffer) {
    CHECK_ITERABLE_APPROX(val, expected[i]);
    ++i;
  }
}

std::vector<DataVector> create_random_data_vectors(
    const gsl::not_null<std::mt19937*> gen, size_t outer_size,
    size_t inner_size) {
  auto value_distribution = std::uniform_real_distribution(
      std::numeric_limits<double>::min(), std::numeric_limits<double>::max());

  std::vector<DataVector> res(outer_size);
  for (size_t i = 0; i < outer_size; ++i) {
    res[i] = make_with_random_values<DataVector>(gen, value_distribution,
                                                 inner_size);
  }
  return res;
}

void check_constructors() {
  MAKE_GENERATOR(gen);
  auto size_distribution = std::uniform_int_distribution(1, 10);
  const auto inner_size = static_cast<size_t>(size_distribution(gen));
  const auto outer_size = static_cast<size_t>(size_distribution(gen));

  CAPTURE(inner_size);
  CAPTURE(outer_size);

  const auto expected =
      create_random_data_vectors(make_not_null(&gen), outer_size, inner_size);

  auto dynamic_buffer = std::make_unique<DynamicBuffer>(outer_size, inner_size);

  for (size_t i = 0; i < outer_size; ++i) {
    dynamic_buffer->at(i) = expected.at(i);
  }

  check_results(*dynamic_buffer, expected);

  DynamicBuffer dynamic_buffer_copied(*dynamic_buffer);
  DynamicBuffer dynamic_buffer_assigned = *dynamic_buffer;
  DynamicBuffer dynamic_buffer_moved = std::move(*dynamic_buffer);

  dynamic_buffer.reset();  // calls destructor

  check_results(dynamic_buffer_copied, expected);
  check_results(dynamic_buffer_assigned, expected);
  check_results(dynamic_buffer_moved, expected);
}

SPECTRE_TEST_CASE("Unit.DataStructures.DynamicBuffer",
                  "[DataStructures][Unit]") {
  for (size_t i = 0; i < 20; ++i) {
    check_constructors();
  }
}

// [[OutputRegex, Must copy into same size]]
[[noreturn]] SPECTRE_TEST_CASE("Unit.DataStructures.DynamicBuffer.wrong_size",
                               "[DataStructures][Unit]") {
  ASSERTION_TEST();
#ifdef SPECTRE_DEBUG
  DynamicBuffer dynamic_buffer(1, 3);
  dynamic_buffer[0] = DataVector{1., 2.};
  ERROR("Failed to trigger ASSERT in an assertion test");
#endif
}
