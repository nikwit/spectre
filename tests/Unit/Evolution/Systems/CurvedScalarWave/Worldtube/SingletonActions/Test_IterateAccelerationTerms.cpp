// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <optional>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/SingletonActions/IterateAccelerationTerms.hpp"
#include "Evolution/Systems/CurvedScalarWave/Worldtube/Tags.hpp"
#include "NumericalAlgorithms/Strahlkorper/Tags.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/KerrSchild.hpp"
#include "Time/Tags/Time.hpp"
#include "Utilities/TMPL.hpp"

namespace CurvedScalarWave::Worldtube {
namespace {

// The grid-frame worldtube scheme does not support the scalar self-force, so
// `IterateAccelerationTerms` is expected to error when applied.
SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.CSW.Worldtube.IterateAccelerationTerms",
    "[Unit][Evolution]") {
  static constexpr size_t Dim = 3;
  const gr::Solutions::KerrSchild kerr_schild(1.4, {{0.1, 0.2, 0.3}},
                                              {{0., 0., 0.}});
  const tnsr::I<double, Dim, Frame::Inertial> pos{{5., 0., 0.}};
  const tnsr::I<double, Dim, Frame::Inertial> vel{{0., 0.1, 0.}};
  const Scalar<double> psi_monopole{0.1};
  const Scalar<double> dt_psi_monopole{0.2};
  const tnsr::i<double, Dim, Frame::Grid> psi_dipole{{0.1, 0.2, 0.3}};
  const tnsr::i<double, Dim, Frame::Grid> dt_psi_dipole{{0.4, 0.5, 0.6}};
  const double charge = 0.1;
  const double mass = 0.1;
  const double time = 10.;
  const double turn_on_time = 20.;
  const double turn_on_interval = 1.;
  const size_t current_iteration = 1;
  auto box = db::create<
      db::AddSimpleTags<
          Tags::AccelerationTerms, Tags::ParticlePositionVelocity<Dim>,
          CurvedScalarWave::Tags::BackgroundSpacetime<
              gr::Solutions::KerrSchild>,
          Stf::Tags::StfTensor<Tags::PsiWorldtube, 0, Dim, Frame::Grid>,
          Stf::Tags::StfTensor<::Tags::dt<Tags::PsiWorldtube>, 0, Dim,
                               Frame::Grid>,
          Stf::Tags::StfTensor<Tags::PsiWorldtube, 1, Dim, Frame::Grid>,
          Stf::Tags::StfTensor<::Tags::dt<Tags::PsiWorldtube>, 1, Dim,
                               Frame::Grid>,
          Tags::Charge, Tags::Mass, ::Tags::Time, Tags::SelfForceTurnOnTime,
          Tags::SelfForceTurnOnInterval, Tags::CurrentIteration>,
      db::AddComputeTags<Tags::BackgroundQuantitiesCompute<Dim>,
                         Tags::GeodesicAccelerationCompute<Dim>>>(
      Scalar<DataVector>{}, std::array<tnsr::I<double, Dim>, 2>{pos, vel},
      kerr_schild, psi_monopole, dt_psi_monopole, psi_dipole, dt_psi_dipole,
      charge, std::make_optional(mass), time, std::make_optional(turn_on_time),
      std::make_optional(turn_on_interval), current_iteration);
  CHECK_THROWS_WITH(
      db::mutate_apply<IterateAccelerationTerms>(make_not_null(&box)),
      Catch::Matchers::ContainsSubstring(
          "The scalar self-force is not supported"));
}
}  // namespace
}  // namespace CurvedScalarWave::Worldtube
