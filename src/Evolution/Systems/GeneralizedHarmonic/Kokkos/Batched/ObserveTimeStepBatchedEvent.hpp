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
#include "Evolution/Kokkos/PackedTags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "IO/Observer/ReductionActions.hpp"
#include "Options/String.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Local.hpp"
#include "Parallel/ParallelComponentHelpers.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::Events {

class ObserveNormsBatched : public Event {
 public:
  static constexpr size_t volume_dim = 3;
  using system = gh::System<volume_dim>;
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
      static constexpr Options::String help = {
          "One of SpacetimeMetric, Pi, Phi."};
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
      "Observe norms from packed batched GH host-synced variables.";

  ObserveNormsBatched() = default;
  ObserveNormsBatched(const std::string& subfile_name,
                      const std::vector<ObserveTensor>& tensors_to_observe)
      : subfile_path_("/" + subfile_name),
        tensors_to_observe_(tensors_to_observe) {}

  using compute_tags_for_observation_box = tmpl::list<>;
  using return_tags = tmpl::list<>;
  using argument_tags =
      tmpl::list<evolution::Kokkos::Tags::PackedTopology<gh::System<3>>,
                 evolution::Kokkos::Tags::PackedEvolutionState<gh::System<3>>>;

  template <typename ArrayIndex, typename ParallelComponent,
            typename Metavariables>
  void operator()(
      const evolution::Kokkos::PackedTopology<gh::System<3>>& packed_topology,
      const evolution::Kokkos::PackedEvolutionState<gh::System<3>>&
          packed_evolution_state,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& /*array_index*/,
      const ParallelComponent* const /*meta*/,
      const ObservationValue& observation_value) const {
    host_variables_type host_vars{packed_topology.total_points};
    copy_to_host(make_not_null(&host_vars),
                 packed_evolution_state.device_variables);
    const auto& spacetime_metric =
        get<gr::Tags::SpacetimeMetric<DataVector, volume_dim, Frame::Inertial>>(
            host_vars);
    const auto& pi =
        get<gh::Tags::Pi<DataVector, volume_dim, Frame::Inertial>>(host_vars);
    const auto& phi =
        get<gh::Tags::Phi<DataVector, volume_dim, Frame::Inertial>>(host_vars);

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
    const auto max_value = [](const DataVector& u) {
      double result = -std::numeric_limits<double>::infinity();
      for (const double value : u) {
        result = std::max(result, value);
      }
      return u.size() == 0 ? 0.0 : result;
    };
    const auto min_value = [](const DataVector& u) {
      double result = std::numeric_limits<double>::infinity();
      for (const double value : u) {
        result = std::min(result, value);
      }
      return u.size() == 0 ? 0.0 : result;
    };
    const auto max_abs = [](const DataVector& u) {
      double max_result = 0.0;
      for (const double value : u) {
        max_result = std::max(max_result, std::abs(value));
      }
      return max_result;
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
    const auto reduce_norm_sum =
        [&max_abs](const std::vector<const DataVector*>& components,
                   const std::string& norm_type) {
          if (components.empty() or components[0]->size() == 0) {
            return 0.0;
          }
          if (norm_type == "L2Norm") {
            double sum_sq = 0.0;
            const size_t num_points = components[0]->size();
            for (const DataVector* comp : components) {
              for (size_t s = 0; s < num_points; ++s) {
                const double value = (*comp)[s];
                sum_sq += value * value;
              }
            }
            return std::sqrt(sum_sq / static_cast<double>(num_points));
          }
          if (norm_type == "Max") {
            double max_result = 0.0;
            for (const DataVector* comp : components) {
              max_result = std::max(max_result, max_abs(*comp));
            }
            return max_result;
          }
          if (norm_type == "Min") {
            double min_result = std::numeric_limits<double>::infinity();
            for (const DataVector* comp : components) {
              for (const double value : *comp) {
                min_result = std::min(min_result, value);
              }
            }
            return std::isinf(min_result) ? 0.0 : min_result;
          }
          ERROR("ObserveNormsBatched only supports Max, Min, and L2Norm.");
        };

    std::vector<std::string> legend{observation_value.name};
    std::vector<double> values{observation_value.value};
    const auto append_individual = [&legend, &values, &reduce_norm](
                                       const std::string& name,
                                       const DataVector& data,
                                       const std::string& norm_type) {
      legend.push_back(name + "_" + norm_type);
      values.push_back(reduce_norm(data, norm_type));
    };
    const auto append_sum =
        [&legend, &values, &reduce_norm_sum](
            const std::string& name,
            const std::vector<const DataVector*>& components,
            const std::string& norm_type) {
          legend.push_back(name + "_" + norm_type);
          values.push_back(reduce_norm_sum(components, norm_type));
        };

    for (const auto& observe_tensor : tensors_to_observe_) {
      const auto& name = observe_tensor.name;
      const auto& norm_type = observe_tensor.norm_type;
      const auto& components = observe_tensor.components;
      if (name == "SpacetimeMetric") {
        if (components == "Individual") {
          for (size_t a = 0; a < volume_dim + 1; ++a) {
            for (size_t b = a; b < volume_dim + 1; ++b) {
              append_individual("SpacetimeMetric_" + std::to_string(a) + "_" +
                                    std::to_string(b),
                                spacetime_metric.get(a, b), norm_type);
            }
          }
        } else if (components == "Sum") {
          std::vector<const DataVector*> tensor_components{};
          for (size_t a = 0; a < volume_dim + 1; ++a) {
            for (size_t b = a; b < volume_dim + 1; ++b) {
              tensor_components.push_back(&spacetime_metric.get(a, b));
            }
          }
          append_sum("SpacetimeMetric", tensor_components, norm_type);
        } else {
          ERROR("ObserveNormsBatched components must be Individual or Sum.");
        }
      } else if (name == "Pi") {
        if (components == "Individual") {
          for (size_t a = 0; a < volume_dim + 1; ++a) {
            for (size_t b = a; b < volume_dim + 1; ++b) {
              append_individual(
                  "Pi_" + std::to_string(a) + "_" + std::to_string(b),
                  pi.get(a, b), norm_type);
            }
          }
        } else if (components == "Sum") {
          std::vector<const DataVector*> tensor_components{};
          for (size_t a = 0; a < volume_dim + 1; ++a) {
            for (size_t b = a; b < volume_dim + 1; ++b) {
              tensor_components.push_back(&pi.get(a, b));
            }
          }
          append_sum("Pi", tensor_components, norm_type);
        } else {
          ERROR("ObserveNormsBatched components must be Individual or Sum.");
        }
      } else if (name == "Phi") {
        if (components == "Individual") {
          for (size_t d = 0; d < volume_dim; ++d) {
            for (size_t a = 0; a < volume_dim + 1; ++a) {
              for (size_t b = a; b < volume_dim + 1; ++b) {
                append_individual("Phi_" + std::to_string(d) + "_" +
                                      std::to_string(a) + "_" +
                                      std::to_string(b),
                                  phi.get(d, a, b), norm_type);
              }
            }
          }
        } else if (components == "Sum") {
          std::vector<const DataVector*> tensor_components{};
          for (size_t d = 0; d < volume_dim; ++d) {
            for (size_t a = 0; a < volume_dim + 1; ++a) {
              for (size_t b = a; b < volume_dim + 1; ++b) {
                tensor_components.push_back(&phi.get(d, a, b));
              }
            }
          }
          append_sum("Phi", tensor_components, norm_type);
        } else {
          ERROR("ObserveNormsBatched components must be Individual or Sum.");
        }
      } else {
        ERROR("ObserveNormsBatched tensor '" << name
                                             << "' is unsupported. Use "
                                                "SpacetimeMetric, Pi, or Phi.");
      }
    }

    auto& reduction_writer = Parallel::get_parallel_component<
        observers::ObserverWriter<Metavariables>>(cache);
    Parallel::threaded_action<
        observers::ThreadedActions::WriteReductionDataRow>(
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
  std::string subfile_path_ = "/BatchedGeneralizedHarmonic";
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

}  // namespace gh::Events
