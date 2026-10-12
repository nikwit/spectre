// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Domain/Tags.hpp"
#include "Evolution/Systems/CurvedScalarWave/BackgroundSpacetime.hpp"
#include "ParallelAlgorithms/Amr/Protocols/Projector.hpp"
#include "Time/Tags/Time.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace CurvedScalarWave {
namespace detail {
template <typename System, typename Background, typename Tags>
struct ProjectBackground;

template <typename System, typename Background, typename... Tags>
struct ProjectBackground<System, Background, tmpl::list<Tags...>>
    : tt::ConformsTo<amr::protocols::Projector> {
  using return_tags = tmpl::list<Tags...>;
  using argument_tags =
      tmpl::list<domain::Tags::Coordinates<System::volume_dim, Frame::Inertial>,
                 ::Tags::Time,
                 CurvedScalarWave::Tags::BackgroundSpacetime<Background>>;

  template <typename AmrData>
  static void apply(const gsl::not_null<typename Tags::type*>... fields,
                    const tnsr::I<DataVector, System::volume_dim>& coordinates,
                    const double time, const Background& background,
                    const AmrData& /*old_mesh_or_parent_or_children*/) {
    const auto values = background.variables(coordinates, time, return_tags{});
    ((*fields = get<Tags>(values)), ...);
  }
};
}  // namespace detail

/// Reevaluate the analytic background after p refinement, splitting or joining.
/// Must run after ProjectDomain; interpolating the old metric is not
/// sufficient.
template <typename System, typename Background>
using ProjectBackground =
    detail::ProjectBackground<System, Background,
                              typename System::spacetime_tag_list>;
}  // namespace CurvedScalarWave
