// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <memory>
#include <sstream>
#include <string>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Evolution/EventsAndDenseTriggers/DenseTrigger.hpp"
#include "Evolution/EventsAndDenseTriggers/DenseTriggers/Filter.hpp"
#include "Evolution/EventsAndDenseTriggers/DenseTriggers/Times.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/Evolution/EventsAndDenseTriggers/DenseTriggers/TestTrigger.hpp"
#include "Options/Options.hpp"
#include "Parallel/RegisterDerivedClassesWithCharm.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Trigger.hpp"
#include "Time/Tags.hpp"
#include "Time/TimeStepId.hpp"
#include "Time/Triggers/TimeCompares.hpp"
#include "Utilities/Registration.hpp"
#include "Utilities/TMPL.hpp"

namespace {
template <typename RegistrarList>
class TestTrigger : public Trigger<RegistrarList> {
 public:
  TestTrigger() = default;
  explicit TestTrigger(CkMigrateMessage* /*unused*/) noexcept {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(TestTrigger);  // NOLINT

  struct Result {
    using type = bool;
    constexpr static Options::String help = "Result";
  };

  using options = tmpl::list<Result>;
  constexpr static Options::String help = "help";

  explicit TestTrigger(const bool result) noexcept : result_(result) {}

  using argument_tags = tmpl::list<>;
  bool operator()() const noexcept { return result_; }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) noexcept override {
    Trigger<RegistrarList>::pup(p);
    p | result_;
  }

 private:
  bool result_{};
};

template <typename RegistrarList>
PUP::able::PUP_ID TestTrigger<RegistrarList>::my_PUP_ID = 0;  // NOLINT

using trigger_registrars = tmpl::list<Registration::Registrar<TestTrigger>>;
using TriggerBase = Trigger<trigger_registrars>;
using DenseTriggerBase = DenseTrigger<
    tmpl::list<DenseTriggers::Registrars::Filter<trigger_registrars>,
               TestHelpers::DenseTriggers::Registrars::TestTrigger>>;

void check(const bool expected_is_ready,
           const bool expected_is_triggered,
           const double expected_next_check,
           const bool dense_is_ready,
           const bool dense_is_triggered,
           const bool non_dense_is_triggered) noexcept {
  std::stringstream creation_string;
  creation_string << std::boolalpha;
  creation_string << "Filter:\n"
                  << "  Trigger:\n"
                  << "    TestTrigger:\n"
                  << "      IsReady: " << dense_is_ready << "\n"
                  << "      IsTriggered: " << dense_is_triggered << "\n"
                  << "      NextCheck: " << expected_next_check << "\n"
                  << "  Filter:\n"
                  << "    TestTrigger:\n"
                  << "      Result: " << non_dense_is_triggered;
  CAPTURE(creation_string.str());
  const auto box = db::create<db::AddSimpleTags<>>();
  const auto trigger = serialize_and_deserialize(
      TestHelpers::test_creation<std::unique_ptr<DenseTriggerBase>>(
          creation_string.str()));
  CHECK(trigger->is_ready(box) == expected_is_ready);
  if (not expected_is_ready) {
    return;
  }
  const auto result = trigger->is_triggered(box);
  CHECK(result.is_triggered == expected_is_triggered);
  CHECK(result.next_check == expected_next_check);
}

void example() noexcept {
  const std::string input = R"(
# [example]
Filter:
  Trigger:
    Times:
      EvenlySpaced:
        Offset: 0
        Interval: 10
  Filter:
    TimeCompares:
      Comparison: GreaterThanOrEqualTo
      Value: 100
# [example]
)";
  using triggers = tmpl::list<Triggers::Registrars::TimeCompares>;
  using dense_triggers = tmpl::list<DenseTriggers::Registrars::Filter<triggers>,
                                    DenseTriggers::Registrars::Times>;
  const auto trigger =
      TestHelpers::test_creation<std::unique_ptr<DenseTrigger<dense_triggers>>>(
          input);
  const auto run_trigger = [&trigger](const double time) noexcept {
    const auto box =
        db::create<db::AddSimpleTags<::Tags::Time, ::Tags::TimeStepId>>(
            time, TimeStepId(true, 0, Slab(0., 1.).start()));
    return trigger->is_triggered(box).is_triggered;
  };
  CHECK(not run_trigger(90.));
  CHECK(not run_trigger(95.));
  CHECK(run_trigger(100.));
  CHECK(not run_trigger(105.));
  CHECK(run_trigger(110.));
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Evolution.EventsAndDenseTriggers.DenseTriggers.Filter",
                  "[Unit][Evolution]") {
  Parallel::register_derived_classes_with_charm<DenseTriggerBase>();
  Parallel::register_derived_classes_with_charm<TriggerBase>();

  check(false, false, 3.5, false, false, false);
  check(false, false, 3.5, false, false,  true);
  check(false, false, 3.5, false,  true, false);
  check(false, false, 3.5, false,  true,  true);
  check( true, false, 3.5,  true, false, false);
  check( true, false, 3.5,  true, false,  true);
  check( true, false, 3.5,  true,  true, false);
  check( true,  true, 3.5,  true,  true,  true);

  example();
}
