# Distributed under the MIT License.
# See LICENSE.txt for details.

import unittest

import numpy as np
from numpy.polynomial.legendre import legval, legval2d, legval3d

from spectre.DataStructures import DataVector
from spectre.Domain import ElementId, ElementMap, Side
from spectre.Domain.Creators import Sphere
from spectre.NumericalAlgorithms.LinearOperators import (
    chalmers_energy_estimate,
    chalmers_estimate,
    chalmers_p_coarsening_error,
    chalmers_parity_estimate,
)
from spectre.Spectral import Basis, Mesh, Quadrature
from spectre.tools.Amr.Chalmers import (
    estimate_sequence,
    estimate_tensor,
    legendre_grid,
    refinement_choice,
    sample_modal,
    to_modal,
)
from spectre.tools.Amr.ChalmersCube import (
    evaluate,
    polynomial_difference,
    quadrature,
)
from spectre.tools.Amr.ChalmersKerrAnalysis import PatchWaveforms, harmonic
from spectre.tools.Amr.ChalmersSingleHole import MULTIPLICITY, PAIRS
from spectre.tools.Amr.ChalmersSingleHole import exact as schwarzschild
from spectre.tools.Amr.ChalmersSphere import curved_difference
from spectre.tools.Amr.ChalmersSphere import exact as spherical_wave
from spectre.tools.Amr.ChalmersSphere import geometry, physical_gradient
from spectre.tools.Amr.ChalmersStudy import approximation_error, compare
from spectre.tools.Amr.ChalmersSymmetry import paired, separate_parity
from spectre.tools.Amr.ChalmersWedges import NORMAL
from spectre.tools.Amr.ChalmersWedges import exact as oblique_wave
from spectre.tools.Amr.ChalmersWedges import interface_fields


class TestChalmers(unittest.TestCase):
    def test_complex_patch_waveforms(self):
        domain = Sphere(1.0, 4.0, 0, 8, True, excise=True).create_domain()
        modes = (
            (0, 0),
            (2, 0),
            (2, 1),
            (2, -1),
            (2, 2),
            (3, 1),
            (4, 1),
            (4, 3),
        )
        extractor = PatchWaveforms(domain, radius=2.5, modes=modes)
        fields = {}
        for block in range(6):
            # Include angular seams and an extraction sphere on a radial seam.
            for side in (Side.Lower, Side.Upper):
                eid = ElementId[3](block).id_of_child(block % 3, side)
                xyz, _, _ = geometry(
                    ElementMap(eid, domain),
                    [legendre_grid(n)[0] for n in (20, 20, 4)],
                )
                x, y, z = xyz
                r2 = np.sum(xyz**2, axis=0)
                # Closed-form harmonics independently check scipy conventions.
                y21 = -np.sqrt(15 / (8 * np.pi)) * (x + 1j * y) * z / r2
                y20 = np.sqrt(5 / (16 * np.pi)) * (3 * z * z / r2 - 1)
                np.testing.assert_allclose(harmonic(2, 1, xyz), y21, atol=1e-14)
                np.testing.assert_allclose(
                    harmonic(2, -1, xyz), -y21.conj(), atol=1e-14
                )
                fields[str(eid)] = [
                    to_modal(0.7 + 0.3 * y20 + 1.2 * y21.real - 0.4 * y21.imag)
                ]
        expected = [
            0.7 * np.sqrt(4 * np.pi),
            0.3,
            0.6 + 0.2j,
            -0.6 + 0.2j,
            0,
            0,
            0,
            0,
        ]
        for order in (24, 32):
            np.testing.assert_allclose(
                extractor(fields, order), expected, atol=5e-10, rtol=0
            )
        # Existing m=0 extraction remains real and uses the same normalization.
        axisymmetric = PatchWaveforms(domain, radius=2.5)(fields)
        self.assertFalse(np.iscomplexobj(axisymmetric))
        np.testing.assert_allclose(
            axisymmetric, [expected[0], 0.3, 0], atol=5e-10, rtol=0
        )
        with self.assertRaises(ValueError):
            PatchWaveforms(domain, modes=((0, 1),))
        with self.assertRaises(ValueError):
            harmonic(1, 2, np.ones((3, 1)))
        with self.assertRaises(ValueError):
            extractor({})

    def test_separate_parity(self):
        for estimator in (paired, separate_parity):
            for bad in ([1.0], [1.0, np.nan], [[1.0, 2.0]]):
                with self.assertRaises(ValueError):
                    estimator(bad)
            with self.assertRaises(ValueError):
                estimator([1.0, 2.0], np.nan)
        # Independent Python (hypergeometric) and C++ tail evaluations.
        rng = np.random.default_rng(842)
        for n in (6, 7, 12, 13, 20):
            for slopes in ((0.005, 0.1), (0.1, 0.10001), (1.2, 0.15)):
                k = np.arange(n)
                a = np.where(k % 2, 1e-3, 1.0) * np.exp(
                    -np.where(k % 2, slopes[1], slopes[0]) * k
                )
                for scale in (1e-100, 1.0, 1e100):
                    expected = separate_parity(scale * a, scale * 1e-5)
                    native = chalmers_parity_estimate(
                        DataVector(scale * a), scale * 1e-5
                    )
                    np.testing.assert_allclose(
                        native["error"], expected.error, rtol=1e-10
                    )
                    np.testing.assert_allclose(
                        native["sigma"], expected.sigma, rtol=1e-10
                    )
        # Both alternatives withstand roundoff perturbations of symmetry zeros.
        a = np.zeros(11)
        a[1::2] = np.exp(-0.8 * np.arange(1, 11, 2))
        for _ in range(20):
            perturbed = a + 1e-13 * rng.normal(size=11)
            for estimator in (paired, separate_parity):
                self.assertAlmostEqual(
                    estimator(perturbed, 3e-5).sigma, 0.8, places=7
                )
            self.assertAlmostEqual(
                chalmers_parity_estimate(DataVector(perturbed), 3e-5)["sigma"],
                0.8,
                places=7,
            )
        # A significant high odd mode cannot be hidden by low-order zeros.
        a = np.zeros(12)
        a[9] = 0.01
        self.assertTrue(
            np.isinf(chalmers_parity_estimate(DataVector(a), 1e-5)["error"])
        )
        # Filter-aware p coarsening must retain the cost of its newly filtered
        # slab even when the retained polynomial has no tail.
        a = np.zeros(9)
        a[0], a[7] = 1.0, 0.2
        mesh = Mesh[1]([9], Basis.Legendre, Quadrature.GaussLobatto)
        values = DataVector(legval(legendre_grid(9)[0], a))
        error = chalmers_p_coarsening_error(values, mesh, 1e-8, 1, True)[0]
        self.assertAlmostEqual(error, 0.2 * np.sqrt(2 / 15), places=12)

    def test_stationary_schwarzschild(self):
        xyz = np.array([[1.7, 2.0, -1.2], [0.5, -0.8, 1.8], [1.0, 0.2, -0.6]])
        state = schwarzschild(xyz)
        metric = np.zeros((3, 4, 4))
        for k, (a, b) in enumerate(PAIRS):
            metric[:, a, b] = metric[:, b, a] = state[k]
        np.testing.assert_allclose(np.linalg.det(metric), -1, atol=1e-14)
        inverse = np.linalg.inv(metric)
        lapse = 1 / np.sqrt(-inverse[:, 0, 0])
        shift = (lapse[:, None] ** 2 * inverse[:, 0, 1:]).T
        phi = state[20:].reshape(3, 10, 3)
        np.testing.assert_allclose(
            -lapse * state[10:20] + np.einsum("ik,ijk->jk", shift, phi),
            0,
            atol=1e-14,
        )
        for d in range(3):
            offset = np.zeros_like(xyz)
            offset[d] = 1e-5
            derivative = (
                schwarzschild(xyz + offset)[:10]
                - schwarzschild(xyz - offset)[:10]
            ) / 2e-5
            np.testing.assert_allclose(derivative, phi[d], atol=1e-9)

    def test_spherical_wave(self):
        self.check_wave_equation(spherical_wave)

    def test_oblique_wave(self):
        self.check_wave_equation(oblique_wave)

    def check_wave_equation(self, solution):
        xyz = np.array([[1.2, 2.0, -1.2], [0.5, -0.8, 1.8], [1.0, 0.2, -0.6]])
        time, delta = 0.3, 1e-5
        state = solution(xyz, time)
        dt = (solution(xyz, time + delta) - solution(xyz, time - delta)) / (
            2 * delta
        )
        dx = []
        for d in range(3):
            offset = np.zeros_like(xyz)
            offset[d] = delta
            dx.append(
                (solution(xyz + offset, time) - solution(xyz - offset, time))
                / (2 * delta)
            )
        np.testing.assert_allclose(dt[0], -state[1], atol=1e-8)
        np.testing.assert_allclose([v[0] for v in dx], state[2:], atol=1e-8)
        np.testing.assert_allclose(
            dt[1] + sum(dx[d][d + 2] for d in range(3)), 0, atol=1e-8
        )
        np.testing.assert_allclose(
            dt[2:] + np.array([v[1] for v in dx]), 0, atol=1e-8
        )

    def test_wedge_interface(self):
        domain = Sphere(1.0, 2.0, 0, 8, True, excise=True).create_domain()
        modes = {}
        # Different angular/radial partitions and unequal extents on the two
        # sides exercise sampling on a nonconforming, oriented block interface.
        for block, dimension in [(3, 1), (4, 2)]:
            for side in [Side.Lower, Side.Upper]:
                element_id = ElementId[3](block).id_of_child(dimension, side)
                axes = [legendre_grid(n)[0] for n in (4, 5, 6)]
                xyz, _, _ = geometry(ElementMap(element_id, domain), axes)
                lo, hi = element_id.segment_ids[1].endpoints
                angle = (lo + hi) / 2 + (hi - lo) / 2 * axes[1]
                psi = np.sum(xyz**2, axis=0) * angle[None, :, None]
                fields = [
                    psi,
                    np.ones_like(psi),
                    *(n * np.ones_like(psi) for n in NORMAL),
                ]
                modes[str(element_id)] = [to_modal(v) for v in fields]
        xyz, weights, traces = interface_fields(modes, domain, 24)
        expected = (
            np.sum(xyz**2, axis=0) * 4 / np.pi * np.arctan2(xyz[2], xyz[0])
        )
        for trace in traces:
            np.testing.assert_allclose(trace[0], expected, atol=1e-12)
            np.testing.assert_allclose(trace[1], 1.0, atol=1e-12)
        np.testing.assert_allclose(traces[0], traces[1], atol=1e-12)
        self.assertAlmostEqual(
            np.sum(weights), 3 * np.arctan(1 / np.sqrt(2)), places=12
        )
        # A discontinuous offset on the lower radial half checks that the
        # integral is split at element boundaries, not on one global face grid.
        lower = str(ElementId[3](4).id_of_child(2, Side.Lower))
        modes[lower][0][0, 0, 0] += 0.01
        _, weights, traces = interface_fields(modes, domain, 24)
        self.assertAlmostEqual(
            np.sqrt(np.sum(weights * (traces[0] - traces[1]) ** 2)),
            0.01 * np.sqrt(1.25 * np.arctan(1 / np.sqrt(2))),
            places=12,
        )
        del modes[next(iter(modes))]
        with self.assertRaises(ValueError):
            interface_fields(modes, domain, 24)

    def test_curved_geometry(self):
        domain = Sphere(1.0, 4.0, 0, 8, True, excise=True).create_domain()
        xi, _, weights = quadrature(24)
        volume = 0.0
        for block in range(6):
            element_map = ElementMap(ElementId[3](block), domain)
            xyz, determinant, inverse = geometry(element_map, [xi] * 3)
            volume += np.sum(weights * determinant)
            # r^2 is quadratic in the radial logical coordinate, even though
            # its physical gradient points in all three Cartesian directions.
            nodal_xyz, _, _ = geometry(
                element_map, [legendre_grid(n)[0] for n in (3, 4, 5)]
            )
            modes = to_modal(np.sum(nodal_xyz**2, axis=0))
            np.testing.assert_allclose(
                physical_gradient(modes, inverse, 24), 2 * xyz, atol=1e-12
            )
        self.assertAlmostEqual(volume, 4 * np.pi / 3 * (4**3 - 1), places=10)

    def test_curved_common_partition(self):
        domain = Sphere(1.0, 4.0, 0, 8, True, excise=True).create_domain()
        parent = ElementId[3](0)

        def state(element_id, extents, offset=0.0):
            xyz, _, _ = geometry(
                ElementMap(element_id, domain),
                [legendre_grid(n)[0] for n in extents],
            )
            modes = to_modal(np.sum(xyz**2, axis=0) + offset)
            return [modes] + [np.zeros_like(modes) for _ in range(4)]

        a = {str(parent): state(parent, (3, 4, 5))}
        children = [
            parent.id_of_child(2, side) for side in [Side.Lower, Side.Upper]
        ]
        b = {str(child): state(child, (4, 5, 6), 0.01) for child in children}
        expected = [0.01 * np.sqrt(4 * np.pi / 3 * (4**3 - 1) / 6), 0, 0, 0, 0]
        np.testing.assert_allclose(
            curved_difference(a, b, domain, 24), expected, atol=1e-12
        )
        with self.assertRaises(ValueError):
            curved_difference(
                a, {str(children[0]): b[str(children[0])]}, domain, 24
            )
        # GH has 50 independent components, representing five symmetric 4x4
        # tensors. A constant offset in each gives 80 full tensor entries.
        gh_a = {name: [fields[0]] * 50 for name, fields in a.items()}
        gh_b = {name: [fields[0]] * 50 for name, fields in b.items()}
        differences = curved_difference(gh_a, gh_b, domain, 24)
        np.testing.assert_allclose(differences, expected[0], atol=1e-12)
        self.assertAlmostEqual(
            np.sqrt(np.dot(MULTIPLICITY, differences**2)),
            np.sqrt(80) * expected[0],
            places=11,
        )
        with self.assertRaisesRegex(ValueError, "component counts"):
            curved_difference(a, gh_b, domain, 24)

    def test_cube_quadrature(self):
        # An independent tensor polynomial with unequal extents checks axis
        # ordering and the physical normalization used by the 3D diagnostics.
        modes = np.zeros((4, 5, 6))
        modes[2, 1, 0] = 1.0
        x, _, weights = quadrature(12)
        expected = (0.5 * (3 * x * x - 1))[:, None, None] * x[None, :, None]
        np.testing.assert_allclose(
            evaluate(modes, 12),
            np.broadcast_to(expected, (12, 12, 12)),
            atol=1e-14,
        )
        jacobian = 2.0 * 3.0 * 4.0
        self.assertAlmostEqual(
            jacobian * np.sum(weights * evaluate(modes, 12) ** 2),
            jacobian * (2.0 / 5) * (2.0 / 3) * 2.0,
        )

    def test_cube_common_partition(self):
        # Compare x on two different anisotropic partitions of the same cube.
        def cell(bounds, offset=0.0):
            c = np.zeros((3, 4, 5))
            c[0, 0, 0] = sum(bounds[0]) / 2 + offset
            c[1, 0, 0] = (bounds[0][1] - bounds[0][0]) / 2
            return bounds, [c] * 5

        a = {
            "left": cell([(-1, 0), (-1, 1), (-1, 1)]),
            "right": cell([(0, 1), (-1, 1), (-1, 1)]),
        }
        b = {
            "bottom": cell([(-1, 1), (-1, 0), (-1, 1)]),
            "top": cell([(-1, 1), (0, 1), (-1, 1)]),
        }
        np.testing.assert_allclose(polynomial_difference(a, b), 0, atol=1e-14)
        b["top"] = cell([(-1, 1), (0, 1), (-1, 1)], offset=0.1)
        np.testing.assert_allclose(polynomial_difference(a, b), 0.2, atol=1e-14)
        del b["bottom"]
        with self.assertRaisesRegex(ValueError, "same volume"):
            polynomial_difference(a, b)

    def test_native_estimator(self):
        # Compare two independent implementations, including tensor ordering,
        # parity zeros, noise, and the two infinite-tail summation algorithms.
        rng = np.random.default_rng(123)
        for shape in [(10,), (7, 8), (6, 7, 8)]:
            nodes = np.meshgrid(
                *(legendre_grid(n)[0] for n in shape), indexing="ij"
            )
            mesh = Mesh[len(shape)](
                shape, Basis.Legendre, Quadrature.GaussLobatto
            )
            for values in [
                np.zeros(shape),
                np.ones(shape),
                np.exp(-sum((x / 0.6) ** 2 for x in nodes)),
                rng.uniform(-1, 1, shape),
            ]:
                native = chalmers_estimate(
                    DataVector(values.ravel(order="F")), mesh
                )
                prototype = estimate_tensor(to_modal(values))
                for actual, expected in zip(native, prototype):
                    for key in ["error", "sigma"]:
                        np.testing.assert_allclose(
                            actual[key],
                            getattr(expected, key),
                            rtol=1e-8,
                            atol=1e-10,
                        )
                    self.assertEqual(
                        actual["insufficient_sequences"],
                        sum(
                            s.status == "insufficient_modes"
                            for s in expected.sequences
                        ),
                    )
                    np.testing.assert_allclose(
                        actual["max_log_fit_rms"],
                        max(
                            (
                                s.log_fit_rms
                                for s in expected.sequences
                                if np.isfinite(s.log_fit_rms)
                            ),
                            default=0,
                        ),
                        rtol=1e-8,
                        atol=1e-10,
                    )

    def test_directional_energy_and_transverse_basis(self):
        # A weak slow stripe must not independently label the dominant,
        # exponentially convergent error nonsmooth. Rotating the transverse
        # orthonormal basis leaves the energy and resulting decision invariant.
        c = np.zeros((8, 7))
        c[:, 0] = 2 * np.exp(-1.2 * np.arange(8))
        c[:, 1] = 1e-4 * np.exp(-0.1 * np.arange(8))
        x, y = np.meshgrid(
            legendre_grid(8)[0], legendre_grid(7)[0], indexing="ij"
        )
        mesh = Mesh[2]([8, 7], Basis.Legendre, Quadrature.GaussLobatto)

        def estimate(modes, energy=True):
            return (chalmers_energy_estimate if energy else chalmers_estimate)(
                DataVector(legval2d(x, y, modes).ravel(order="F")), mesh, 1e-4
            )[0]

        original = estimate(c)
        self.assertGreater(original["error"], 1e-4)
        self.assertGreater(original["sigma"], 1)
        self.assertGreater(estimate(c, energy=False)["error"], 1e-4)
        self.assertLess(estimate(c, energy=False)["sigma"], 1)
        norms = np.sqrt(2 / (2 * np.arange(7) + 1))
        rng = np.random.default_rng(321)
        rotation, _ = np.linalg.qr(rng.normal(size=(7, 7)))
        rotated = ((c * norms) @ rotation) / norms
        for key in ["error", "sigma", "max_log_fit_rms"]:
            np.testing.assert_allclose(
                estimate(rotated)[key], original[key], rtol=1e-9, atol=1e-12
            )
        spectrum = np.linalg.norm(c * norms, axis=1)
        expected = estimate_sequence(spectrum)
        np.testing.assert_allclose(
            original["error"], expected.error, rtol=1e-10
        )
        np.testing.assert_allclose(
            original["sigma"], expected.sigma, rtol=1e-10
        )
        # Energy includes a strong transverse contribution rather than hiding
        # it with the guard or allowing cancellation between stripes.
        c[:, 1] = np.exp(-0.1 * np.arange(8))
        self.assertLess(estimate(c)["sigma"], 1)
        self.assertGreater(estimate(c)["error"], original["error"])

    def test_directional_energy_constant_and_parity(self):
        mesh = Mesh[2]([3, 3], Basis.Legendre, Quadrature.GaussLobatto)
        x, y = np.meshgrid(
            legendre_grid(3)[0], legendre_grid(3)[0], indexing="ij"
        )
        values = DataVector((1 + 1e-10 * (3 * y**2 - 1) / 2).ravel(order="F"))
        fit = chalmers_energy_estimate(values, mesh, 1e-7)
        self.assertEqual(fit[1]["error"], 0)
        self.assertEqual(fit[1]["ignored_sequences"], 1)
        fit = chalmers_energy_estimate(values, mesh, 1e-12)
        self.assertTrue(np.isinf(fit[1]["error"]))
        self.assertEqual(fit[1]["insufficient_sequences"], 1)
        for parity in [0, 1]:
            c = np.exp(-0.4 * np.arange(14))
            c[parity::2] = 0
            one_mesh = Mesh[1]([14], Basis.Legendre, Quadrature.GaussLobatto)
            fit = chalmers_energy_estimate(
                DataVector(legval(legendre_grid(14)[0], c)), one_mesh
            )[0]
            np.testing.assert_allclose(fit["sigma"], 0.4, rtol=1e-11)

    def test_exponential_and_infinite_tail(self):
        # Independently sum the squared Legendre norms, including the last
        # retained term, to check the hypergeometric expression.
        for sigma in [0.005, 0.1, 0.5, 1.0, 2.0]:
            coeffs = 2.3 * np.exp(-sigma * np.arange(10))
            estimate = estimate_sequence(coeffs)
            self.assertAlmostEqual(estimate.sigma, sigma)
            self.assertLess(estimate.log_fit_rms, 1e-13)
            modes = np.arange(9, 10000)
            expected = np.sqrt(
                np.sum(
                    2 * 2.3**2 * np.exp(-2 * sigma * modes) / (2 * modes + 1)
                )
            )
            np.testing.assert_allclose(estimate.error, expected, rtol=1e-12)
            self.assertEqual(
                refinement_choice(estimate, 0), "p" if sigma > 1 else "h"
            )

    def test_zero_polynomials_and_roundoff(self):
        for coeffs in [np.zeros(8), [1, 2, 3, 0, 0, 0, 0, 0]]:
            estimate = estimate_sequence(coeffs)
            self.assertEqual(estimate.error, 0)
            self.assertEqual(refinement_choice(estimate, 1e-8), "hold")
        coeffs = np.array([1, 2, 0, 0, 0, 1e-16, 2e-16, 1e-16])
        self.assertEqual(estimate_sequence(coeffs).status, "terminated_tail")
        # A full-degree polynomial is not assumed resolved from its samples.
        estimate = estimate_sequence([0, 0, 0, 1])
        self.assertEqual(estimate.status, "insufficient_modes")
        self.assertEqual(refinement_choice(estimate, 1e-8), "undetermined")

    def test_parity_uses_original_indices(self):
        for parity in [0, 1]:
            coeffs = np.exp(-0.4 * np.arange(14))
            coeffs[parity::2] = 0
            estimate = estimate_sequence(coeffs)
            self.assertAlmostEqual(estimate.sigma, 0.4)
            self.assertEqual(estimate.used_modes, 7)
            self.assertLess(estimate.log_fit_rms, 1e-14)

    def test_plateau_and_nondecaying(self):
        coeffs = np.exp(-1.5 * np.arange(16))
        coeffs[8:] = coeffs[8]
        estimate = estimate_sequence(coeffs)
        self.assertGreater(estimate.log_fit_rms, 1)
        self.assertLess(estimate.sigma, 1)
        self.assertEqual(refinement_choice(estimate, 0), "h")
        for coeffs in [np.ones(6), np.exp(np.arange(6))]:
            estimate = estimate_sequence(coeffs)
            self.assertTrue(np.isinf(estimate.error))
            self.assertEqual(estimate.status, "nondecaying")

    def test_scale_invariance(self):
        coeffs = np.exp(-0.3 * np.arange(10))
        reference = estimate_sequence(coeffs)
        for factor in [1e-100, 1e-9, 1e9, 1e100]:
            result = estimate_sequence(factor * coeffs)
            self.assertAlmostEqual(result.sigma, reference.sigma)
            np.testing.assert_allclose(
                result.error / factor, reference.error, rtol=1e-13
            )

    def test_tensor_directions(self):
        sigmas = [0.3, 1.2, 0.7]
        sequences = [
            np.exp(-sigma * np.arange(n)) for sigma, n in zip(sigmas, [6, 7, 8])
        ]
        coeffs = np.einsum("i,j,k->ijk", *sequences)
        estimates = estimate_tensor(coeffs)
        np.testing.assert_allclose([e.sigma for e in estimates], sigmas)
        np.testing.assert_allclose(
            [e.error for e in estimates],
            [estimate_sequence(c).error for c in sequences],
        )
        self.assertEqual(
            [refinement_choice(e, 0) for e in estimates], ["h", "p", "h"]
        )

    def test_transverse_roundoff_and_insufficient_modes(self):
        coeffs = np.full((8, 8), 1e-17)
        coeffs[:, 0] = np.exp(-1.2 * np.arange(8))
        result = estimate_tensor(coeffs)
        self.assertAlmostEqual(result[0].sigma, 1.2)
        self.assertEqual(result[1].error, 0)
        coeffs[-1, 1] = 0.1
        result = estimate_tensor(coeffs)
        self.assertTrue(np.isnan(result[0].sigma))
        self.assertTrue(np.isinf(result[0].error))

    def test_spectre_normalization(self):
        x = legendre_grid(7)[0]
        np.testing.assert_allclose(
            to_modal(3 + x + x**2), [10 / 3, 1, 2 / 3, 0, 0, 0, 0], atol=2e-14
        )
        coeffs = np.array([1.2, -0.7, 2.5, 0.01, 0.02, 0.001, 0.0003])
        np.testing.assert_allclose(
            to_modal(legval(x, coeffs)), coeffs, atol=2e-14
        )

    def test_tensor_transform_and_mapped_box(self):
        rng = np.random.default_rng(1043)
        coeffs = rng.uniform(-1, 1, (4, 5, 6))
        nodes = np.meshgrid(
            *(legendre_grid(n)[0] for n in coeffs.shape), indexing="ij"
        )
        np.testing.assert_allclose(
            to_modal(legval3d(*nodes, coeffs)), coeffs, atol=2e-14
        )
        result = sample_modal(lambda x: x**2, (5,), [(2, 4)])
        # x = 3 + xi, so x^2 = 28/3 P0 + 6 P1 + 2/3 P2.
        np.testing.assert_allclose(result, [28 / 3, 6, 2 / 3, 0, 0], atol=2e-14)

    def test_invalid_inputs(self):
        for coeffs in [[1], [[1, 2]], [np.nan, 1], [np.inf, 2]]:
            with self.assertRaises(ValueError):
                estimate_sequence(coeffs)
        for floor in [-1, np.inf, np.nan]:
            with self.assertRaises(ValueError):
                estimate_sequence([1, 2, 3], coefficient_floor=floor)
        with self.assertRaises(ValueError):
            estimate_tensor(np.ones((2, 2, 2, 2)))
        for tolerance in [-1, np.inf, np.nan]:
            with self.assertRaises(ValueError):
                refinement_choice(estimate_sequence([1, 2, 3]), tolerance)

    def test_independent_error_and_child_jacobians(self):
        # Linear LGL interpolation of x^2 is 1 on [-1,1]. Its exact L2
        # error is 4/sqrt(15). Splitting into two linear elements reduces
        # that error by four; quadratic interpolation is exact.
        rows = compare(
            lambda x: x**2,
            (2,),
            order=8,
            tolerance=1e-6,
            metadata={"family": "quadratic"},
        )
        self.assertEqual(len(rows), 1)
        np.testing.assert_allclose(rows[0]["error_base"], 4 / np.sqrt(15))
        np.testing.assert_allclose(rows[0]["error_h"], 1 / np.sqrt(15))
        self.assertLess(rows[0]["error_p"], 1e-14)
        self.assertEqual(rows[0]["points_p"], 3)
        self.assertEqual(rows[0]["points_h"], 4)
        error = approximation_error(
            lambda x, y, z: x**2,
            (2, 3, 4),
            [(-1, 1), (-1, 1), (-1, 1)],
            8,
        )
        np.testing.assert_allclose(error, 8 / np.sqrt(15))

    def test_sampling_blind_spot(self):
        # This analytic sinusoid aliases to a constant on three LGL nodes.
        # A coefficient-only estimator cannot see the unsampled oscillation.
        function = lambda x: np.sin(2 * np.pi * x + 0.37)
        result = estimate_tensor(sample_modal(function, (3,)))[0]
        self.assertEqual(refinement_choice(result, 1e-6), "hold")
        self.assertGreater(
            approximation_error(function, (3,), [(-1, 1)], 32), 1.0
        )


if __name__ == "__main__":
    unittest.main()
