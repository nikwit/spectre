// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Framework/TestingFramework.hpp"

#include <pup.h>

#include "Evolution/EventsAndDenseTriggers/DenseTrigger.hpp"
#include "Options/Options.hpp"
#include "Parallel/CharmPupable.hpp"
#include "Utilities/Registration.hpp"
#include "Utilities/TMPL.hpp"

namespace TestHelpers::DenseTriggers {
template <typename RegistrarList>
class TestTrigger;

namespace Registrars {
using TestTrigger =
    Registration::Registrar<::TestHelpers::DenseTriggers::TestTrigger>;
}  // namespace Registrars

template <typename RegistrarList = tmpl::list<Registrars::TestTrigger>>
class TestTrigger : public DenseTrigger<RegistrarList> {
 public:
  /// \cond
  TestTrigger() = default;
  explicit TestTrigger(CkMigrateMessage* const msg) noexcept
      : DenseTrigger<RegistrarList>(msg) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(TestTrigger);  // NOLINT
  /// \endcond

  using Result = typename DenseTrigger<RegistrarList>::Result;

  struct IsReady {
    using type = bool;
    constexpr static Options::String help = "IsReady";
  };

  struct IsTriggered {
    using type = bool;
    constexpr static Options::String help = "IsTriggered";
  };

  struct NextCheck {
    using type = double;
    constexpr static Options::String help = "NextCheck";
  };

  using options = tmpl::list<IsReady, IsTriggered, NextCheck>;
  constexpr static Options::String help = "help";

  TestTrigger(const bool is_ready, const bool is_triggered,
              const double next_check) noexcept
      : is_ready_(is_ready),
        is_triggered_(is_triggered),
        next_check_(next_check) {}

  using is_triggered_argument_tags = tmpl::list<>;
  Result is_triggered() const noexcept {
    CHECK(is_ready_);
    return {is_triggered_, next_check_};
  }

  using is_ready_argument_tags = tmpl::list<>;
  bool is_ready() const noexcept { return is_ready_; }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) noexcept override {
    DenseTrigger<RegistrarList>::pup(p);
    p | is_ready_;
    p | is_triggered_;
    p | next_check_;
  }

 private:
  bool is_ready_;
  bool is_triggered_;
  double next_check_;
};

/// \cond
template <typename RegistrarList>
PUP::able::PUP_ID TestTrigger<RegistrarList>::my_PUP_ID = 0;  // NOLINT
/// \endcond
}  // namespace TestHelpers::DenseTriggers
