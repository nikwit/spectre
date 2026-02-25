// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <tuple>
#include <vector>

#include <pup.h>

#include "DataStructures/VariablesKokkos.hpp"
#include "Evolution/Executables/ScalarWave/Batched/Tags.hpp"
#include "Evolution/Systems/ScalarWave/System.hpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "IO/Observer/ReductionActions.hpp"
#include "Options/String.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Local.hpp"
#include "Parallel/ParallelComponentHelpers.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarWave::Events {

class ObserveNormsBatched : public Event {
 public:
  static constexpr size_t volume_dim = 3;
  using system = ScalarWave::System<volume_dim>;
  using host_variables_type = typename system::variables_tag::type;

  struct SubfileName {
    using type = std::string;
    static constexpr Options::String help = {
        "Name of reduction subfile (without leading slash)."};
  };
  struct ObserveTensor {
    static constexpr Options::String help = {
        "Tensor name and norm settings for batched reduction output."};
    struct Name {
      using type = std::string;
      static constexpr Options::String help = {"One of Psi, Pi, Phi."};
    };
    struct NormType {
      using type = std::string;
      static constexpr Options::String help = {"One of Max, Min, L2Norm."};
    };
    struct Components {
      using type = std::string;
      static constexpr Options::String help = {"Individual or Sum."};
    };
    using options = tmpl::list<Name, NormType, Components>;

    std::string name{};
    std::string norm_type{};
    std::string components{};

    // NOLINTNEXTLINE(google-runtime-references)
    void pup(PUP::er& p) {
      p | name;
      p | norm_type;
      p | components;
    }
  };
  struct TensorsToObserve {
    using type = std::vector<ObserveTensor>;
    static constexpr Options::String help = {
        "List of tensors and requested norms."};
  };

  explicit ObserveNormsBatched(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(ObserveNormsBatched);  // NOLINT

  using options = tmpl::list<SubfileName, TensorsToObserve>;
  static constexpr Options::String help =
      "Observe norms from packed batched ScalarWave host-synced variables.";

  ObserveNormsBatched() = default;
  ObserveNormsBatched(const std::string& subfile_name,
                      const std::vector<ObserveTensor>& tensors_to_observe)
      : subfile_path_("/" + subfile_name), tensors_to_observe_(tensors_to_observe) {}

  using compute_tags_for_observation_box = tmpl::list<>;
  using return_tags = tmpl::list<>;
  using argument_tags = tmpl::list<ScalarWave::Batched::Tags::DeviceData>;

  template <typename ArrayIndex, typename ParallelComponent,
            typename Metavariables>
  void operator()(const ScalarWave::Batched::DeviceData& device_data,
                  Parallel::GlobalCache<Metavariables>& cache,
                  const ArrayIndex& /*array_index*/,
                  const ParallelComponent* const /*meta*/,
                  const ObservationValue& observation_value) const {
    host_variables_type host_vars{device_data.total_points()};
    copy_to_host(make_not_null(&host_vars), device_data.device_variables());
    const auto& psi = get(get<ScalarWave::Tags::Psi>(host_vars));
    const auto& pi = get(get<ScalarWave::Tags::Pi>(host_vars));

    const auto l2_norm = [](const DataVector& u) {
      if (u.size() == 0) {
        return 0.0;
      }
      double sum_sq = 0.0;
      for (const double value : u) {
        sum_sq += value * value;
      }
      return std::sqrt(sum_sq / static_cast<double>(u.size()));
    };

    const auto max_abs = [](const DataVector& u) {
      double max_value = 0.0;
      for (const double value : u) {
        max_value = std::max(max_value, std::abs(value));
      }
      return max_value;
    };
    const auto max_value = [](const DataVector& u) {
      double result = -std::numeric_limits<double>::infinity();
      for (const double value : u) {
        result = std::max(result, value);
      }
      return result;
    };
    const auto min_value = [](const DataVector& u) {
      double result = std::numeric_limits<double>::infinity();
      for (const double value : u) {
        result = std::min(result, value);
      }
      return result;
    };
    const auto reduce_norm = [&l2_norm, &max_value, &min_value](
                                 const DataVector& u,
                                 const std::string& norm_type) {
      if (norm_type == "Max") {
        return max_value(u);
      }
      if (norm_type == "Min") {
        return min_value(u);
      }
      if (norm_type == "L2Norm") {
        return l2_norm(u);
      }
      ERROR("ObserveNormsBatched only supports Max, Min, and L2Norm. Got '"
            << norm_type << "'.");
    };

    std::vector<std::string> legend{observation_value.name};
    std::vector<double> values{observation_value.value};
    const auto& phi = get<ScalarWave::Tags::Phi<volume_dim>>(host_vars);
    auto append_scalar_norm = [&legend, &values, &reduce_norm](
                                  const std::string& name,
                                  const DataVector& data,
                                  const std::string& norm_type) {
      legend.push_back(name + "_" + norm_type);
      values.push_back(reduce_norm(data, norm_type));
    };
    for (const auto& observe_tensor : tensors_to_observe_) {
      const auto& name = observe_tensor.name;
      const auto& norm_type = observe_tensor.norm_type;
      const auto& components = observe_tensor.components;
      if (name == "Psi") {
        append_scalar_norm("Psi", psi, norm_type);
      } else if (name == "Pi") {
        append_scalar_norm("Pi", pi, norm_type);
      } else if (name == "Phi") {
        if (components == "Individual") {
          append_scalar_norm("Phi_x", phi.get(0), norm_type);
          append_scalar_norm("Phi_y", phi.get(1), norm_type);
          append_scalar_norm("Phi_z", phi.get(2), norm_type);
        } else if (components == "Sum") {
          if (norm_type == "L2Norm") {
            double sum_sq = 0.0;
            for (size_t s = 0; s < phi.get(0).size(); ++s) {
              sum_sq += phi.get(0)[s] * phi.get(0)[s] +
                        phi.get(1)[s] * phi.get(1)[s] +
                        phi.get(2)[s] * phi.get(2)[s];
            }
            const double l2_phi =
                phi.get(0).size() == 0
                    ? 0.0
                    : std::sqrt(sum_sq / static_cast<double>(phi.get(0).size()));
            legend.push_back("Phi_" + norm_type);
            values.push_back(l2_phi);
          } else if (norm_type == "Max") {
            legend.push_back("Phi_" + norm_type);
            values.push_back(
                std::max({max_abs(phi.get(0)), max_abs(phi.get(1)),
                          max_abs(phi.get(2))}));
          } else if (norm_type == "Min") {
            double min_phi = std::numeric_limits<double>::infinity();
            for (size_t s = 0; s < phi.get(0).size(); ++s) {
              min_phi =
                  std::min(min_phi, std::min({phi.get(0)[s], phi.get(1)[s], phi.get(2)[s]}));
            }
            if (phi.get(0).size() == 0) {
              min_phi = 0.0;
            }
            legend.push_back("Phi_" + norm_type);
            values.push_back(min_phi);
          } else {
            ERROR("ObserveNormsBatched only supports Max, Min, and L2Norm.");
          }
        } else {
          ERROR("ObserveNormsBatched components must be Individual or Sum.");
        }
      } else {
        ERROR("ObserveNormsBatched tensor '" << name
                                             << "' is unsupported. Use Psi, Pi, Phi.");
      }
    }

    auto& reduction_writer = Parallel::get_parallel_component<
        observers::ObserverWriter<Metavariables>>(cache);
    Parallel::threaded_action<observers::ThreadedActions::WriteReductionDataRow>(
        reduction_writer[0], subfile_path_, std::move(legend),
        std::make_tuple(std::move(values)));
  }

  using is_ready_argument_tags = tmpl::list<>;
  template <typename Metavariables, typename ArrayIndex, typename Component>
  bool is_ready(Parallel::GlobalCache<Metavariables>& /*cache*/,
                const ArrayIndex& /*array_index*/,
                const Component* const /*meta*/) const {
    return true;
  }

  bool needs_evolved_variables() const override { return false; }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) override {
    Event::pup(p);
    p | subfile_path_;
    p | tensors_to_observe_;
  }

 private:
  std::string subfile_path_ = "/BatchedScalarWave";
  std::vector<ObserveTensor> tensors_to_observe_{};
};

class CompletionBatched : public Event {
 public:
  explicit CompletionBatched(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(CompletionBatched);  // NOLINT

  using compute_tags_for_observation_box = tmpl::list<>;
  using options = tmpl::list<>;
  static constexpr Options::String help = {
      "Set terminate on the local batched singleton driver."};

  CompletionBatched() = default;

  using return_tags = tmpl::list<>;
  using argument_tags = tmpl::list<>;

  template <typename Metavariables, typename ArrayIndex, typename Component>
  void operator()(Parallel::GlobalCache<Metavariables>& cache,
                  const ArrayIndex& /*array_index*/,
                  const Component* const /*meta*/,
                  const ObservationValue& /*observation_value*/) const {
    auto* const local_component =
        Parallel::local(Parallel::get_parallel_component<Component>(cache));
    ASSERT(local_component != nullptr,
           "Failed to find local singleton component for completion.");
    local_component->set_terminate(true);
  }

  using is_ready_argument_tags = tmpl::list<>;
  template <typename Metavariables, typename ArrayIndex, typename Component>
  bool is_ready(Parallel::GlobalCache<Metavariables>& /*cache*/,
                const ArrayIndex& /*array_index*/,
                const Component* const /*meta*/) const {
    return true;
  }

  bool needs_evolved_variables() const override { return false; }
};

}  // namespace ScalarWave::Events
