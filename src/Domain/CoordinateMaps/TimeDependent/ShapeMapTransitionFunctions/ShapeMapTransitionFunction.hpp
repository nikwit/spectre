// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <optional>
#include <pup.h>

#include "DataStructures/DataVector.hpp"
#include "Parallel/CharmPupable.hpp"

namespace domain::CoordinateMaps {

class SphereTransition;
/**
 * \brief Abstract base class for the transition functions used by the shape
 *map.
 *
 * \details This base class defines the required methods of a transition
 * function used by the shape map. Different domains require the shape map to
 * fall off towards the boundary in different ways. This behavior is controlled
 * by the transition function. It is also needed to find the inverse of the
 * shape map. Since the shape map preserves angles, the problem of finding its
 * inverse reduces to the 1-dimensional problem of finding the original radius
 * from the mapped radius. The mapped radius \f$\tilde{r}\f$ is related to the
 * original \f$r\f$ radius by:
 * \f{equation}{\label{eq:shape_map_radius}
 *  \tilde{r} = r (1 - f(r,\theta,\phi) \sum_{lm} \lambda_{lm}(t)Y_{lm}(\theta,
 *  \phi)),
 *  \f}
 * where \f$f(r,\theta,\phi)\f$ is the transition function. Depending
 * on the format of the transition function, it should be possible to
 * analytically derive this map's inverse because it preserves angles and shifts
 * only the radius of each point. Otherwise the inverse has to be computed
 * numerically.
 *
 * The transition function must also be able to compute the gradient and the
 * value of the function divided by the radius. Care must be taken that this
 * does not divide by zero.
 *
 * All member functions with the exception of `original_radius_over_radius`
 * exist as overloads for types `double` and `DataVector` so that they work with
 * the templated shape map methods calling them. To avoid code duplication these
 * can be forwarded to templated implementation methods held by the derived
 * classes only.
 *
 * For an example, see SphereTransition.
 **/
class ShapeMapTransitionFunction : public PUP::able {
 public:
  ShapeMapTransitionFunction() = default;
  using creatable_classes = tmpl::list<SphereTransition>;

  /** Evaluate the transition function at the Cartesian coordinates
   *`source_coords`.
   **/
  virtual double operator()(
      const std::array<double, 3>& source_coords) const = 0;
  virtual DataVector operator()(
      const std::array<DataVector, 3>& source_coords) const = 0;

  /** Given the mapped coordinates `target_coords` and the corresponding
   * spherical harmonic expansion \f$\sum_{lm} \lambda_{lm}(t)Y_{lm}\f$,
   * `distorted_radius`, this method evaluates the original radius from the
   * mapped radius by inverting the shape map \f$\ref{eq:shape_map_radius}\f$.
   * It also divides by the mapped radius to simplify calculations in the shape
   * map.
   **/
  virtual std::optional<double> original_radius_over_radius(
      const std::array<double, 3>& target_coords,
      double distorted_radius) const = 0;

  /** Evaluate the transition function at the Cartesian coordinates divided by
   * the radius. Care must be taken not to divide by zero. **/
  virtual double map_over_radius(
      const std::array<double, 3>& source_coords) const = 0;
  virtual DataVector map_over_radius(
      const std::array<DataVector, 3>& source_coords) const = 0;

  /** Evaluate the gradient of the transition function at the Cartesian
   * coordinates `source_coords`. **/
  virtual std::array<double, 3> gradient(
      const std::array<double, 3>& source_coords) const = 0;
  virtual std::array<DataVector, 3> gradient(
      const std::array<DataVector, 3>& source_coords) const = 0;

  virtual bool operator==(const ShapeMapTransitionFunction& other) const = 0;
  virtual bool operator!=(const ShapeMapTransitionFunction& other) const = 0;

  WRAPPED_PUPable_abstract(ShapeMapTransitionFunction);
  explicit ShapeMapTransitionFunction(CkMigrateMessage* m) : PUP::able(m) {}
  void pup(PUP::er& p) override { PUP::able::pup(p); };
};
}  // namespace domain::CoordinateMaps
