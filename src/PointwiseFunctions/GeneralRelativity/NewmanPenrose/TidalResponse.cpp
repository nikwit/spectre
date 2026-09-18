// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/TidalResponse.hpp"

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "PointwiseFunctions/GeneralRelativity/NewmanPenrose/Tetrad.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"

namespace gr::np {
QuadrupoleProfiles quadrupole_profiles(const Scalar<DataVector>& radius,
                                       const double mass) {
  const DataVector compactness = mass / get(radius);
  const DataVector f = 1. - 2. * compactness;
  for (size_t p = 0; p < f.size(); ++p) {
    if (f[p] <= 0.) {
      ERROR("The profile radius "
            << get(radius)[p] << " lies inside the horizon of mass " << mass);
    }
  }
  const DataVector root_f = sqrt(f);
  return QuadrupoleProfiles{DataVector(f.size(), 1.),
                            root_f * (1. + 2. * compactness),
                            f,
                            DataVector(f.size(), 1.),
                            root_f,
                            f};
}

namespace {
template <typename T>
tnsr::ii<T, 3, Frame::Inertial> stf_from_components_impl(
    const std::array<T, 5>& components) {
  tnsr::ii<T, 3, Frame::Inertial> tensor{};
  get<0, 0>(tensor) = components[0];
  get<1, 1>(tensor) = components[1];
  get<2, 2>(tensor) = -components[0] - components[1];
  get<0, 1>(tensor) = components[2];
  get<0, 2>(tensor) = components[3];
  get<1, 2>(tensor) = components[4];
  return tensor;
}

// The r_hat-irreducible pieces of eq. tide-irreducible for one constant real
// STF tensor: scalar = T:dd, vector = P T d, transverse = P T P + scalar P/2
struct IrreducibleParts {
  DataVector scalar;
  std::array<DataVector, 3> vector;
  std::array<std::array<DataVector, 3>, 3> transverse;
};

IrreducibleParts irreducible_parts(
    const tnsr::ii<double, 3, Frame::Inertial>& tensor,
    const TriadVector& direction) {
  const size_t num_points = get_size(get<0>(direction));
  IrreducibleParts parts{};
  std::array<DataVector, 3> contracted{};
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(contracted, i) = DataVector(num_points, 0.);
    for (size_t j = 0; j < 3; ++j) {
      gsl::at(contracted, i) += tensor.get(i, j) * direction.get(j);
    }
  }
  parts.scalar = DataVector(num_points, 0.);
  for (size_t i = 0; i < 3; ++i) {
    parts.scalar += direction.get(i) * gsl::at(contracted, i);
  }
  for (size_t i = 0; i < 3; ++i) {
    gsl::at(parts.vector, i) =
        gsl::at(contracted, i) - parts.scalar * direction.get(i);
  }
  std::array<std::array<DataVector, 3>, 3> projector{};
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      gsl::at(gsl::at(projector, i), j) =
          (i == j ? 1. : 0.) - direction.get(i) * direction.get(j);
    }
  }
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      DataVector& entry = gsl::at(gsl::at(parts.transverse, i), j);
      entry = 0.5 * parts.scalar * gsl::at(gsl::at(projector, i), j);
      for (size_t k = 0; k < 3; ++k) {
        for (size_t l = 0; l < 3; ++l) {
          entry += gsl::at(gsl::at(projector, i), k) * tensor.get(k, l) *
                   gsl::at(gsl::at(projector, l), j);
        }
      }
    }
  }
  return parts;
}

// eq. tidal-tensor for one parity sector
std::array<std::array<DataVector, 3>, 3> assemble(const IrreducibleParts& parts,
                                                  const TriadVector& direction,
                                                  const DataVector& profile_l,
                                                  const DataVector& profile_v,
                                                  const DataVector& profile_t,
                                                  const bool transverse_only) {
  std::array<std::array<DataVector, 3>, 3> result{};
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      DataVector& entry = gsl::at(gsl::at(result, i), j);
      entry = profile_t * gsl::at(gsl::at(parts.transverse, i), j);
      if (not transverse_only) {
        entry +=
            1.5 * profile_l * parts.scalar *
            (direction.get(i) * direction.get(j) - (i == j ? 1. : 0.) / 3.);
        entry += profile_v * (direction.get(i) * gsl::at(parts.vector, j) +
                              gsl::at(parts.vector, i) * direction.get(j));
      }
    }
  }
  return result;
}
}  // namespace

tnsr::ii<double, 3, Frame::Inertial> stf_from_components(
    const std::array<double, 5>& components) {
  return stf_from_components_impl(components);
}

tnsr::ii<std::complex<double>, 3, Frame::Inertial> stf_from_components(
    const std::array<std::complex<double>, 5>& components) {
  return stf_from_components_impl(components);
}

ComplexMatrix quadrupole_tide_tensor(
    const tnsr::ii<double, 3, Frame::Inertial>& electric,
    const tnsr::ii<double, 3, Frame::Inertial>& magnetic,
    const TriadVector& direction, const Scalar<DataVector>& radius,
    const double mass, const bool transverse_only) {
  const QuadrupoleProfiles profiles = quadrupole_profiles(radius, mass);
  const auto electric_part =
      assemble(irreducible_parts(electric, direction), direction, profiles.e_l,
               profiles.e_v, profiles.e_t, transverse_only);
  const auto magnetic_part =
      assemble(irreducible_parts(magnetic, direction), direction, profiles.b_l,
               profiles.b_v, profiles.b_t, transverse_only);
  const std::complex<double> imaginary_unit{0., 1.};
  ComplexMatrix result(get_size(get(radius)), std::complex<double>{0., 0.});
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      result.get(i, j) = gsl::at(gsl::at(electric_part, i), j) +
                         imaginary_unit * gsl::at(gsl::at(magnetic_part, i), j);
    }
  }
  return result;
}

ComplexMatrix self_dual_boost(const TriadVector& axis,
                              const Scalar<DataVector>& rapidity) {
  const size_t num_points = get_size(get(rapidity));
  const std::complex<double> imaginary_unit{0., 1.};
  const DataVector cosh_eta = cosh(get(rapidity));
  const DataVector sinh_eta = sinh(get(rapidity));
  ComplexMatrix result(num_points, std::complex<double>{0., 0.});
  // [a]_x
  std::array<std::array<DataVector, 3>, 3> cross{};
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      gsl::at(gsl::at(cross, i), j) = DataVector(num_points, 0.);
    }
  }
  cross[0][1] = -axis.get(2);
  cross[0][2] = axis.get(1);
  cross[1][0] = axis.get(2);
  cross[1][2] = -axis.get(0);
  cross[2][0] = -axis.get(1);
  cross[2][1] = axis.get(0);
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      const DataVector parallel = axis.get(i) * axis.get(j);
      result.get(i, j) =
          parallel + cosh_eta * ((i == j ? 1. : 0.) - parallel) +
          imaginary_unit * sinh_eta * gsl::at(gsl::at(cross, i), j);
    }
  }
  return result;
}

ComplexMatrix rest_to_slice_map(const TriadVector& radial_direction,
                                const TriadVector& transverse_velocity,
                                const Scalar<DataVector>& rapidity) {
  const size_t num_points = get_size(get(rapidity));
  const DataVector speed = sqrt(square(get<0>(transverse_velocity)) +
                                square(get<1>(transverse_velocity)) +
                                square(get<2>(transverse_velocity)));
  // At a zero of the transverse speed the boost is the identity and any
  // regular axis represents it
  TriadVector axis(num_points, 0.);
  for (size_t p = 0; p < num_points; ++p) {
    for (size_t i = 0; i < 3; ++i) {
      axis.get(i)[p] = speed[p] > 1.e-14
                           ? transverse_velocity.get(i)[p] / speed[p]
                           : radial_direction.get(i)[p];
    }
  }
  const ComplexMatrix transverse_map =
      self_dual_boost(axis, Scalar<DataVector>{atanh(speed)});
  const ComplexMatrix radial_map = self_dual_boost(radial_direction, rapidity);
  ComplexMatrix result(num_points, std::complex<double>{0., 0.});
  for (size_t i = 0; i < 3; ++i) {
    for (size_t k = 0; k < 3; ++k) {
      for (size_t j = 0; j < 3; ++j) {
        result.get(i, k) += transverse_map.get(i, j) * radial_map.get(j, k);
      }
    }
  }
  return result;
}

std::array<WeylScalars, 5> direct_tide_scalar_columns(
    const TriadVector& radial_direction, const TriadVector& transverse_velocity,
    const Scalar<DataVector>& rapidity, const Scalar<DataVector>& radius,
    const RealMatrix& adapted_rotation, const double mass) {
  const size_t num_points = get_size(get(rapidity));
  const ComplexMatrix boost =
      rest_to_slice_map(radial_direction, transverse_velocity, rapidity);
  const tnsr::ii<double, 3, Frame::Inertial> zero{0.};
  std::array<WeylScalars, 5> columns{};
  for (size_t component = 0; component < 5; ++component) {
    std::array<double, 5> unit{};
    gsl::at(unit, component) = 1.;
    const ComplexMatrix tide_rest = quadrupole_tide_tensor(
        stf_from_components(unit), zero, radial_direction, radius, mass, true);
    // Boost to the slice: B T B^T (no conjugation, the action is complex
    // linear)
    ComplexMatrix tide_slice(num_points, std::complex<double>{0., 0.});
    for (size_t i = 0; i < 3; ++i) {
      for (size_t j = 0; j < 3; ++j) {
        for (size_t k = 0; k < 3; ++k) {
          for (size_t l = 0; l < 3; ++l) {
            tide_slice.get(i, j) +=
                boost.get(i, k) * tide_rest.get(k, l) * boost.get(j, l);
          }
        }
      }
    }
    gsl::at(columns, component) = weyl_scalars_from_tidal_tensor(
        rotate_symmetric(tide_slice, adapted_rotation));
  }
  return columns;
}
}  // namespace gr::np
