# Distributed under the MIT License.
# See LICENSE.txt for details.

import csv
import json
import os
import shutil
import unittest

import numpy as np
import numpy.testing as npt
from click.testing import CliRunner

import spectre.IO.H5 as spectre_h5
from spectre.Domain import ElementId, serialize_domain
from spectre.Domain.Creators import Brick
from spectre.Informer import unit_test_build_path
from spectre.IO.Exporter import ModalSpacetimeInterpolator, Verbosity
from spectre.IO.Exporter.ValidateModalSpacetimeInterpolator import (
    _interior_collapsed_indices,
    _summarize,
    validate_modal_spacetime_interpolator_command,
)
from spectre.IO.H5 import ElementVolumeData, TensorComponent
from spectre.Spectral import Basis, Mesh, Quadrature


class TestModalSpacetimeInterpolator(unittest.TestCase):
    def setUp(self):
        self.test_dir = os.path.join(
            unit_test_build_path(), "IO/Exporter/ModalSpacetimeInterpolator"
        )
        self.h5_filename = os.path.join(self.test_dir, "VolumeData.h5")
        os.makedirs(self.test_dir, exist_ok=True)

        domain = Brick(
            lower_bounds=[0.0, 0.0, 0.0],
            upper_bounds=[1.0, 1.0, 1.0],
            initial_refinement_levels=[0, 0, 0],
            initial_num_points=[4, 4, 4],
            is_periodic=[False, False, False],
        ).create_domain()
        serialized_domain = serialize_domain(domain)
        element_id = ElementId[3]("[B0,(L0I0,L0I0,L0I0)]")

        def field(time):
            return 1.0 + 0.2 * time - 0.03 * time**2 + 0.005 * time**3

        self.field = field
        subfiles = [
            (
                "VeryCoarse",
                Mesh[3](2, Basis.Legendre, Quadrature.GaussLobatto),
                np.linspace(0.0, 4.0, 17),
            ),
            (
                "Full",
                Mesh[3](4, Basis.Legendre, Quadrature.GaussLobatto),
                np.linspace(0.0, 4.0, 5),
            ),
            (
                "Validation",
                Mesh[3](4, Basis.Legendre, Quadrature.GaussLobatto),
                [0.375, 2.25],
            ),
        ]
        with spectre_h5.H5File(self.h5_filename, "w") as h5file:
            for subfile_name, mesh, times in subfiles:
                volfile = h5file.insert_vol("/" + subfile_name, version=0)
                for observation_id, observation_time in enumerate(times):
                    tensor_components = [
                        TensorComponent(
                            "Psi",
                            np.full(
                                mesh.number_of_grid_points(),
                                field(observation_time),
                            ),
                        )
                    ]
                    if subfile_name == "Validation":
                        tensor_components.extend(
                            TensorComponent(
                                f"InertialCoordinates_{axis}",
                                np.full(mesh.number_of_grid_points(), value),
                            )
                            for axis, value in zip("xyz", [0.25, 0.4, 0.6])
                        )
                    volfile.write_volume_data(
                        observation_id=observation_id,
                        observation_value=observation_time,
                        elements=[
                            ElementVolumeData(
                                element_id,
                                tensor_components,
                                mesh,
                            )
                        ],
                        serialized_domain=serialized_domain,
                    )
                h5file.close_current_object()

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_binding(self):
        interpolator = ModalSpacetimeInterpolator[3](
            self.h5_filename,
            subfiles_in_priority_order=["VeryCoarse", "Full"],
            tensor_components=["Psi"],
            verbosity=Verbosity.Silent,
        )
        npt.assert_equal(interpolator.time_bounds(), [0.0, 4.0])
        self.assertEqual(interpolator.tensor_components(), ["Psi"])
        for time in [0.0, 0.375, 2.25, 4.0]:
            (result,) = interpolator.interpolate_to_point(
                np.array([0.25, 0.4, 0.6]), time=time
            )
            self.assertAlmostEqual(result, self.field(time), places=11)

    def test_validation_command(self):
        output_filename = os.path.join(self.test_dir, "validation.csv")
        result = CliRunner().invoke(
            validate_modal_spacetime_interpolator_command,
            [
                self.h5_filename,
                "--training-subfile",
                "VeryCoarse",
                "--training-subfile",
                "Full",
                "--validation-subfile",
                "Validation",
                "--var",
                "Psi",
                "--samples-per-element",
                "2",
                "--output",
                output_filename,
            ],
        )
        self.assertEqual(result.exit_code, 0, msg=result.output)
        with open(output_filename, newline="") as output_file:
            rows = list(csv.DictReader(output_file))
        self.assertEqual(len(rows), 2)
        for row in rows:
            self.assertEqual(row["component"], "Psi")
            self.assertEqual(int(row["sample_count"]), 2)
            self.assertLess(float(row["max_abs_error"]), 1.0e-11)
        with open(output_filename.removesuffix(".csv") + ".json") as metadata:
            self.assertEqual(
                json.load(metadata)["validation_times"], [0.375, 2.25]
            )

    def test_validation_helpers(self):
        indices = _interior_collapsed_indices((5, 6, 7), 12, "seed")
        npt.assert_equal(
            indices,
            _interior_collapsed_indices((5, 6, 7), 12, "seed"),
        )
        multi_indices = np.unravel_index(indices, (5, 6, 7), order="F")
        for index, extent in zip(multi_indices, (5, 6, 7)):
            self.assertTrue(np.all(index > 0))
            self.assertTrue(np.all(index < extent - 1))

        summary = _summarize([1.0, 2.0], [1.5, 1.0])
        self.assertEqual(summary["sample_count"], 2)
        self.assertEqual(summary["sum_squared_error"], 1.25)
        self.assertEqual(summary["sum_squared_reference"], 5.0)
        self.assertAlmostEqual(summary["relative_l2_error"], 0.5)
        self.assertEqual(summary["max_abs_error"], 1.0)
        self.assertEqual(summary["worst_index"], 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
