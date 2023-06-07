// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/FunctionsOfTime/IntegratedFoT.hpp"

#include <algorithm>
#include <iterator>
#include <memory>
#include <ostream>
#include <pup.h>
#include <pup_stl.h>
#include <utility>  // IWYU pragma: keep

#include "DataStructures/DataVector.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Literals.hpp"
#include "Utilities/MakeArray.hpp"

namespace domain::FunctionsOfTime {
IntegratedFoT::IntegratedFoT(const double t, value_type initial_func_and_derivs,
                             const double expiration_time)
    : deriv_info_at_update_times_{{t, std::move(initial_func_and_derivs)}},
      expiration_time_(expiration_time) {}

std::unique_ptr<FunctionOfTime> IntegratedFoT::get_clone() const {
  return std::make_unique<IntegratedFoT>(*this);
}

template <size_t MaxDerivReturned>
std::array<DataVector, MaxDerivReturned + 1> IntegratedFoT::func_and_derivs(
    const double t) const {
  static_assert(MaxDerivReturned <= 2, "can only return 2 derivatives");
  if (t > expiration_time_) {
    ERROR("Attempt to evaluate IntegratedFoT at a time "
          << t << " that is after the expiration time " << expiration_time_
          << ". The difference between times is " << t - expiration_time_
          << ".");
  }
  const auto& deriv_info_at_t =
      stored_info_from_upper_bound(t, deriv_info_at_update_times_);
  std::array<DataVector, MaxDerivReturned + 1> result{};
  for (size_t i = 0; i <= MaxDerivReturned; ++i) {
    result.at(i) = deriv_info_at_t.stored_quantities.at(i);
  }

  return result;
}

void IntegratedFoT::update(
    // Clang-tidy says to use 'const DataVector& updated_max_deriv'.
    // However, updated_max_deriv is std::moved out of inside this function.
    // NOLINTNEXTLINE(performance-unnecessary-value-param)
    const double time_of_update, DataVector updated_max_deriv,
    const double next_expiration_time) {
  if (time_of_update <= deriv_info_at_update_times_.back().time) {
    ERROR("t must be increasing from call to call. "
          << "Attempted to update at time " << time_of_update
          << ", which precedes the previous update time of "
          << deriv_info_at_update_times_.back().time << ".");
  }
  if (next_expiration_time < expiration_time_) {
    ERROR("expiration_time must be nondecreasing from call to call. "
          << "Attempted to change expiration time to " << next_expiration_time
          << ", which precedes the previous expiration time of "
          << expiration_time_ << ".");
  }
  if (time_of_update < expiration_time_) {
    ERROR("Attempt to update IntegratedFoT at a time "
          << time_of_update
          << " that is earlier than the previous expiration time of "
          << expiration_time_
          << ". This is bad because some asynchronous process may have already "
             "used IntegratedFoT at a time later than the current time "
          << time_of_update << ".");
  }
  if (time_of_update > next_expiration_time) {
    ERROR(
        "Attempt to set the expiration time of IntegratedFoT "
        "to a value "
        << next_expiration_time << " that is earlier than the current time "
        << time_of_update << ".");
  }

  // Normally, func_and_derivs(t) throws an error if t>expiration_time_.
  // But here, we want to allow time_of_update to
  // be greater than the *previous* expiration time, so we need to
  // reset expiration_time_ before the call to func_and_derivs.
  expiration_time_ = next_expiration_time;

  // get the current values, before updating the `MaxDeriv'th deriv
  value_type func = func_and_derivs(time_of_update);
  ASSERT(updated_max_deriv.size() == 3,
         "The size of the DataVector should be 3: position, velocity and "
         "acceleration");
  for (size_t i = 0; i < 3; ++i) {
    func[i][0] = updated_max_deriv[i];
  }
  deriv_info_at_update_times_.emplace_back(time_of_update, std::move(func));
}

void IntegratedFoT::reset_expiration_time(const double next_expiration_time) {
  FunctionOfTimeHelpers::reset_expiration_time(make_not_null(&expiration_time_),
                                               next_expiration_time);
}

void IntegratedFoT::pup(PUP::er& p) {
  FunctionOfTime::pup(p);
  size_t version = 0;
  p | version;
  // Remember to increment the version number when making changes to this
  // function. Retain support for unpacking data written by previous versions
  // whenever possible. See `Domain` docs for details.
  if (version >= 0) {
    p | deriv_info_at_update_times_;
    p | expiration_time_;
  }
}

bool operator==(const IntegratedFoT& lhs, const IntegratedFoT& rhs) {
  return lhs.deriv_info_at_update_times_ == rhs.deriv_info_at_update_times_ and
         lhs.expiration_time_ == rhs.expiration_time_;
}

bool operator!=(const IntegratedFoT& lhs, const IntegratedFoT& rhs) {
  return not(lhs == rhs);
}

std::ostream& operator<<(std::ostream& os,
                         const IntegratedFoT& integrated_fot) {
  os << integrated_fot.deriv_info_at_update_times_.back();
  os << "\n";

  return os;
}
PUP::able::PUP_ID IntegratedFoT::my_PUP_ID = 0;
#define DIMRETURNED(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data)                             \
  template std::array<DataVector, DIMRETURNED(data) + 1> \
  IntegratedFoT::func_and_derivs<DIMRETURNED(data)>(const double) const;

GENERATE_INSTANTIATIONS(INSTANTIATE, (0, 1, 2))

#undef DIMRETURNED
#undef INSTANTIATE
}  // namespace domain::FunctionsOfTime
