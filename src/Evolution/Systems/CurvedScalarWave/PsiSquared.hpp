// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <string>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"

#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"

#include "Evolution/Systems/CurvedScalarWave/TagsDeclarations.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/Gsl.hpp"

namespace CurvedScalarWave::Tags {
/*!
 * \brief The square of the scalar field \f$\Psi\f$.
 */
struct PsiSquared : db::SimpleTag {
  using type = Scalar<DataVector>;
};

/*!
 * \brief Compute tag that calculates the square of the scalar field \f$\Psi\f$.
 */
struct PsiSquaredCompute : PsiSquared, db::ComputeTag {
  using base = PsiSquared;
  using return_type = Scalar<DataVector>;
  using argument_tags = tmpl::list<CurvedScalarWave::Tags::Psi>;
  static void function(const gsl::not_null<Scalar<DataVector>*> psi_squared,
                       const Scalar<DataVector>& psi) {
    get(*psi_squared) = get(psi) * get(psi);
  }
};

template <size_t Dim, typename Frame>
struct StressEnergyFlux : db::SimpleTag {
  using type = tnsr::I<DataVector, Dim, Frame>;
};

template <size_t Dim, typename Frame>
struct StressEnergyFluxCompute : StressEnergyFlux<Dim, Frame>, db::ComputeTag {
  using base = StressEnergyFlux<Dim, Frame>;
  using return_type = Scalar<DataVector>;
  using argument_tags =
      tmpl::list<CurvedScalarWave::Tags::Psi, CurvedScalarWave::Tags::Pi,
                 CurvedScalarWave::Tags::Phi<Dim>,
                 gr::Tags::SpacetimeMetric<DataVector, Dim, Frame>,
                 gr::Tags::InverseSpacetimeMetric<DataVector, Dim, Frame>,
                 gr::Tags::Lapse<DataVector>,
                 gr::Tags::Shift<DataVector, Dim, Frame>>;

  static void function(
      const gsl::not_null<tnsr::I<DataVector, Dim, Frame>*> stress_energy_flux,
      const Scalar<DataVector>& psi, const Scalar<DataVector>& pi,
      const tnsr::i<DataVector, Dim, Frame>& phi,
      const tnsr::aa<DataVector, Dim, Frame>& spacetime_metric,
      const tnsr::AA<DataVector, Dim, Frame>& inverse_spacetime_metric,
      const Scalar<DataVector>& lapse,
      const tnsr::I<DataVector, Dim, Frame>& shift) {
    tnsr::a<DataVector, Dim, Frame> dmu_psi(get(psi).size());
    get<0>(dmu_psi) = -get(lapse) * get(pi);
    for (size_t i = 0; i < Dim; ++i) {
      get<0>(dmu_psi) += shift.get(i) * phi.get(i);
      dmu_psi.get(i + 1) = phi.get(i);
    }
    tnsr::A<DataVector, Dim, Frame> timelike_killing_vector(get(psi).size(),
                                                            0.);
    get<0>(timelike_killing_vector) = 1.;
    const auto stress_energy_tensor = tenex::evaluate<ti::a, ti::b>(
        0.25 * M_1_PI *
        (dmu_psi(ti::a) * dmu_psi(ti::b) -
         0.5 * spacetime_metric(ti::a, ti::b) *
             inverse_spacetime_metric(ti::C, ti::D) * dmu_psi(ti::c) *
             dmu_psi(ti::d)));

    const auto spacetime_stress_energy_flux = tenex::evaluate<ti::A>(
        inverse_spacetime_metric(ti::A, ti::B) *
        stress_energy_tensor(ti::b, ti::c) * timelike_killing_vector(ti::C));
    for (size_t i = 0; i < Dim; ++i) {
      stress_energy_flux->get(i) = spacetime_stress_energy_flux.get(i + 1);
    }
  }
};
}  // namespace CurvedScalarWave::Tags
