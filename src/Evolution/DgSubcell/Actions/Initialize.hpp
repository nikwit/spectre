// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <tuple>
#include <type_traits>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/DgSubcell/ActiveGrid.hpp"
#include "Evolution/DgSubcell/Mesh.hpp"
#include "Evolution/DgSubcell/NeighborData.hpp"
#include "Evolution/DgSubcell/Projection.hpp"
#include "Evolution/DgSubcell/Tags/ActiveGrid.hpp"
#include "Evolution/DgSubcell/Tags/Coordinates.hpp"
#include "Evolution/DgSubcell/Tags/DidRollback.hpp"
#include "Evolution/DgSubcell/Tags/Inactive.hpp"
#include "Evolution/DgSubcell/Tags/Jacobians.hpp"
#include "Evolution/DgSubcell/Tags/Mesh.hpp"
#include "Evolution/DgSubcell/Tags/NeighborData.hpp"
#include "Evolution/DgSubcell/Tags/SubcellOptions.hpp"
#include "Evolution/DgSubcell/Tags/TciGridHistory.hpp"
#include "Evolution/DgSubcell/Tags/TciStatus.hpp"
#include "Evolution/DgSubcell/TciStatus.hpp"
#include "Evolution/Initialization/SetVariables.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace evolution::dg::subcell::Actions {
/*!
 * \brief Initialize the subcell grid, and switch from DG to subcell if the DG
 * solution is inadmissible.
 *
 * Interior cells are marked as troubled if
 * `subcell_options.always_use_subcells()` is `true` or if
 * `Metavariables::SubcellOptions::DgInitialDataTci::apply` reports that the
 * initial data is not well represented on the DG grid for that cell. Exterior
 * cells are never marked as troubled, because subcell doesn't yet support
 * boundary conditions.
 *
 * If the cell is troubled then `Tags::ActiveGrid` is set to
 * `subcell::ActiveGrid::Subcell`, the `System::variables_tag` become the
 * variables on the subcell grid set by calling
 * `evolution::Initialization::Actions::SetVariables`, and
 * ` Tags::Inactive<System::variables_tag>` become the (inadmissible) DG
 * solution.
 *
 * \details `Metavariables::SubcellOptions::DgInitialDataTci::apply` is called
 * with the evolved variables on the DG grid, the projected evolved variables,
 * the DG mesh, the initial RDMP parameters \f$\delta_0\f$ and \f$\epsilon\f$,
 * and the Persson TCI parameter \f$\alpha\f$. The apply function must return a
 * `bool` that is `true` if the cell is troubled.
 *
 * GlobalCache:
 * - Uses:
 *   - `subcell::Tags::SubcellOptions`
 *
 * DataBox:
 * - Uses:
 *   - `domain::Tags::Mesh<Dim>`
 *   - `domain::Tags::Element<Dim>`
 *   - `System::variables_tag`
 * - Adds:
 *   - `subcell::Tags::Mesh<Dim>`
 *   - `subcell::Tags::ActiveGrid`
 *   - `subcell::Tags::DidRollback`
 *   - `subcell::Tags::Inactive<System::variables_tag>`
 *   - `subcell::Tags::TciGridHistory`
 *   - `subcell::Tags::NeighborDataForReconstructionAndRdmpTci<Dim>`
 *   - `subcell::fd::Tags::InverseJacobianLogicalToGrid<Dim>`
 *   - `subcell::fd::Tags::DetInverseJacobianLogicalToGrid`
 *   - `subcell::Tags::LogicalCoordinates<Dim>`
 *   - `subcell::Tags::Corodinates<Dim, Frame::Grid>` (as compute tag)
 *   - `subcell::Tags::TciStatusCompute<Dim>`
 * - Removes: nothing
 * - Modifies:
 *   - `System::variables_tag` if the cell is troubled
 */
template <typename Metavariables>
struct Initialize {
  using const_global_cache_tags = tmpl::list<Tags::SubcellOptions>;

  using simple_tags = tmpl::list<
      Tags::Mesh<Metavariables::volume_dim>, Tags::ActiveGrid,
      Tags::DidRollback,
      Tags::Inactive<typename Metavariables::system::variables_tag>,
      Tags::TciGridHistory,
      Tags::NeighborDataForReconstructionAndRdmpTci<Metavariables::volume_dim>,
      fd::Tags::InverseJacobianLogicalToGrid<Metavariables::volume_dim>,
      fd::Tags::DetInverseJacobianLogicalToGrid>;
  using compute_tags = tmpl::list<
      Tags::LogicalCoordinatesCompute<Metavariables::volume_dim>,
      ::domain::Tags::MappedCoordinates<
          ::domain::Tags::ElementMap<Metavariables::volume_dim, Frame::Grid>,
          subcell::Tags::Coordinates<Metavariables::volume_dim, Frame::Logical>,
          subcell::Tags::Coordinates>,
      Tags::TciStatusCompute<Metavariables::volume_dim>>;

  template <typename DbTagsList, typename... InboxTags, typename ArrayIndex,
            typename ActionList, typename ParallelComponent>
  static std::tuple<db::DataBox<DbTagsList>&&> apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& inboxes,
      const Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& array_index, ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) noexcept {
    const SubcellOptions& subcell_options = db::get<Tags::SubcellOptions>(box);
    const Mesh<Metavariables::volume_dim>& dg_mesh =
        db::get<::domain::Tags::Mesh<Metavariables::volume_dim>>(box);
    const Mesh<Metavariables::volume_dim> subcell_mesh = fd::mesh(dg_mesh);
    // Note: we currently cannot do subcell at boundaries, so only set on
    // interior elements.
    bool cell_is_troubled =
        subcell_options.always_use_subcells() and
        db::get<::domain::Tags::Element<Metavariables::volume_dim>>(box)
            .external_boundaries()
            .empty();
    db::mutate<subcell::Tags::Mesh<Metavariables::volume_dim>, Tags::ActiveGrid,
               Tags::DidRollback,
               Tags::Inactive<typename Metavariables::system::variables_tag>,
               typename Metavariables::system::variables_tag>(
        make_not_null(&box),
        [&cell_is_troubled, &dg_mesh, &subcell_mesh, &subcell_options](
            const gsl::not_null<Mesh<Metavariables::volume_dim>*>
                subcell_mesh_ptr,
            const gsl::not_null<ActiveGrid*> active_grid_ptr,
            const gsl::not_null<bool*> did_rollback_ptr,
            const auto inactive_vars_ptr, const auto active_vars_ptr) noexcept {
          // We don't consider setting the initial grid to subcell as rolling
          // back. Since no time step is undone, we just continue on the
          // subcells as a normal solve.
          *did_rollback_ptr = false;

          *subcell_mesh_ptr = subcell_mesh;
          *active_grid_ptr = ActiveGrid::Dg;
          fd::project(inactive_vars_ptr, *active_vars_ptr, dg_mesh,
                      subcell_mesh.extents());
          // Now check if the DG solution is admissible
          cell_is_troubled |=
              Metavariables::SubcellOptions::DgInitialDataTci::apply(
                  *active_vars_ptr, *inactive_vars_ptr, dg_mesh,
                  subcell_options.initial_data_rdmp_delta0(),
                  subcell_options.initial_data_rdmp_epsilon(),
                  subcell_options.initial_data_persson_exponent());
          if (cell_is_troubled) {
            // Swap grid
            *active_grid_ptr = ActiveGrid::Subcell;
            using std::swap;
            swap(*active_vars_ptr, *inactive_vars_ptr);
          }
        });
    if (cell_is_troubled) {
      // Set variables on subcells.
      evolution::Initialization::Actions::SetVariables<
          Tags::Coordinates<Metavariables::volume_dim, Frame::Logical>>::
          apply(box, inboxes, cache, array_index, ActionList{},
                std::add_pointer_t<ParallelComponent>{nullptr});
    }
    return std::forward_as_tuple(std::move(box));
  }
};
}  // namespace evolution::dg::subcell::Actions
