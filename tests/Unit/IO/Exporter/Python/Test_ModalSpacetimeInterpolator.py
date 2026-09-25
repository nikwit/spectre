# Distributed under the MIT License.
# See LICENSE.txt for details.

import csv
import json
import os
import shutil
import subprocess
import unittest
from pathlib import Path
from unittest.mock import patch

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

    def test_serialization(self):
        interpolator = ModalSpacetimeInterpolator[3](
            self.h5_filename,
            ["VeryCoarse", "Full"],
            ["Psi"],
            verbosity=Verbosity.Silent,
        )
        saved = os.path.join(self.test_dir, "interpolator.bin")
        interpolator.save(saved)
        self.assertFalse(os.path.exists(saved + ".partial"))
        # The saved object is independent of all source H5 files.
        os.remove(self.h5_filename)
        # A fresh interpreter must register the domain's polymorphic classes
        # on load, without first constructing another interpolator.
        python_spectre = (
            Path(unit_test_build_path()).parents[1] / "bin" / "python-spectre"
        )
        cold_load = subprocess.run(
            [
                str(python_spectre),
                "-c",
                (
                    "import sys; import numpy as np; from spectre.IO.Exporter"
                    " import ModalSpacetimeInterpolator; obj ="
                    " ModalSpacetimeInterpolator[3].load(sys.argv[1]);"
                    " np.testing.assert_allclose(obj.interpolate_to_point("
                    "np.array([0.25, 0.4, 0.6]), 2.25), [1 + .2*2.25 -"
                    " .03*2.25**2 +"
                    " .005*2.25**3], rtol=0, atol=1.e-11)"
                ),
                saved,
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(cold_load.returncode, 0, msg=cold_load.stderr)
        restored = ModalSpacetimeInterpolator[3].load(saved)
        npt.assert_equal(restored.time_bounds(), interpolator.time_bounds())
        self.assertEqual(restored.tensor_components(), ["Psi"])
        for time in [0.0, 0.375, 2.25, 4.0]:
            npt.assert_equal(
                restored.interpolate_to_point(np.array([0.25, 0.4, 0.6]), time),
                interpolator.interpolate_to_point(
                    np.array([0.25, 0.4, 0.6]), time
                ),
            )
        with self.assertRaisesRegex(RuntimeError, "Refusing to overwrite"):
            interpolator.save(saved)
        with self.assertRaisesRegex(RuntimeError, "header or dimension"):
            ModalSpacetimeInterpolator[2].load(saved)
        truncated = os.path.join(self.test_dir, "truncated.bin")
        Path(truncated).write_bytes(Path(saved).read_bytes()[:-1])
        with self.assertRaisesRegex(RuntimeError, "possibly truncated"):
            ModalSpacetimeInterpolator[3].load(truncated)
        Path(truncated).write_bytes(b"bad header")
        with self.assertRaisesRegex(RuntimeError, "header or dimension"):
            ModalSpacetimeInterpolator[3].load(truncated)

    def test_build_only_and_reload(self):
        saved = os.path.join(self.test_dir, "interpolator.bin")
        output = os.path.join(self.test_dir, "validation.csv")
        runner = CliRunner()
        result = runner.invoke(
            validate_modal_spacetime_interpolator_command,
            [
                self.h5_filename,
                "--training-subfile",
                "VeryCoarse",
                "--training-subfile",
                "Full",
                "--save-interpolator",
                saved,
                "--build-only",
                "--output",
                output,
            ],
        )
        self.assertEqual(result.exit_code, 0, msg=result.output)
        self.assertTrue(os.path.isfile(saved))
        self.assertTrue(os.path.isfile(saved + ".json"))
        self.assertFalse(os.path.exists(output))
        cls = ModalSpacetimeInterpolator[3]
        # Fail if validation tries to construct from the training data again.
        with patch(
            "spectre.IO.Exporter.ValidateModalSpacetimeInterpolator."
            "ModalSpacetimeInterpolator"
        ) as bindings:
            bindings.__getitem__.return_value.load = cls.load
            bindings.__getitem__.return_value.side_effect = AssertionError(
                "Unexpected reconstruction"
            )
            result = runner.invoke(
                validate_modal_spacetime_interpolator_command,
                [
                    self.h5_filename,
                    "--load-interpolator",
                    saved,
                    "--training-subfile",
                    "DoesNotExist",
                    "--output",
                    output,
                ],
            )
        self.assertEqual(result.exit_code, 0, msg=result.output)
        with open(output) as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 2)
        self.assertTrue(all(float(r["max_abs_error"]) < 1.0e-11 for r in rows))
        metadata = json.loads(Path(output).with_suffix(".json").read_text())
        self.assertEqual(metadata["construction_seconds"], 0.0)
        self.assertEqual(metadata["loaded_interpolator"], saved)

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
                "--observation-batch-size",
                "4",
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
            run_metadata = json.load(metadata)
            self.assertEqual(run_metadata["validation_times"], [0.375, 2.25])
            self.assertEqual(run_metadata["observation_batch_size"], 4)

    def test_validation_multiple_files(self):
        domain = Brick(
            lower_bounds=[0.0, 0.0, 0.0],
            upper_bounds=[1.0, 1.0, 1.0],
            initial_refinement_levels=[2, 0, 0],
            initial_num_points=[4, 4, 4],
            is_periodic=[False, False, False],
        ).create_domain()
        files = []
        for shard in range(4):
            filename = os.path.join(self.test_dir, f"BbhVolume{shard}.h5")
            files.append(filename)
            with spectre_h5.H5File(filename, "w") as h5file:
                for name, extent, times in [
                    ("VeryCoarse", 2, np.linspace(0.0, 4.0, 17)),
                    ("Coarse", 3, np.linspace(0.0, 4.0, 9)),
                    ("Full", 4, np.linspace(0.0, 4.0, 5)),
                    ("Validation", 4, [0.375, 2.25]),
                ]:
                    volfile = h5file.insert_vol(name, version=0)
                    for obs_id, time in enumerate(times):
                        value = self.field(time) + shard
                        components = [
                            TensorComponent("Psi", np.full(extent**3, value)),
                            TensorComponent(
                                "Chi", np.full(extent**3, 2 * value)
                            ),
                        ]
                        if name == "Validation":
                            components.extend(
                                TensorComponent(
                                    f"InertialCoordinates_{axis}",
                                    np.full(extent**3, coordinate),
                                )
                                for axis, coordinate in zip(
                                    "xyz", [(shard + 0.5) / 4, 0.4, 0.6]
                                )
                            )
                        volfile.write_volume_data(
                            observation_id=obs_id,
                            observation_value=time,
                            elements=[
                                ElementVolumeData(
                                    ElementId[3](
                                        f"[B0,(L2I{shard},L0I0,L0I0)]"
                                    ),
                                    components,
                                    Mesh[3](
                                        extent,
                                        Basis.Legendre,
                                        Quadrature.GaussLobatto,
                                    ),
                                )
                            ],
                            serialized_domain=serialize_domain(domain),
                        )
                    h5file.close_current_object()
        output = os.path.join(self.test_dir, "multiple-files.csv")
        result = CliRunner().invoke(
            validate_modal_spacetime_interpolator_command,
            [*files, "--samples-per-element", "8", "--output", output],
        )
        self.assertEqual(result.exit_code, 0, msg=result.output)
        with open(output, newline="") as output_file:
            rows = list(csv.DictReader(output_file))
        self.assertEqual(len(rows), 4)
        self.assertEqual({row["component"] for row in rows}, {"Psi", "Chi"})
        for row in rows:
            self.assertEqual(int(row["sample_count"]), 32)
            self.assertLess(float(row["max_abs_error"]), 1.0e-11)
        self.assertIn("Reading validation file 4/4", result.output)
        self.assertIn("relative L2=", result.output)

    def test_validation_outside_time_bounds(self):
        with spectre_h5.H5File(self.h5_filename, "a") as h5file:
            volfile = h5file.get_vol("Validation")
            volfile.write_volume_data(
                observation_id=2,
                observation_value=5.0,
                elements=[
                    ElementVolumeData(
                        ElementId[3]("[B0,(L0I0,L0I0,L0I0)]"),
                        [TensorComponent("Psi", np.ones(4**3))],
                        Mesh[3](4, Basis.Legendre, Quadrature.GaussLobatto),
                    )
                ],
            )
        output = os.path.join(self.test_dir, "outside.csv")
        result = CliRunner().invoke(
            validate_modal_spacetime_interpolator_command,
            [
                self.h5_filename,
                "--training-subfile",
                "VeryCoarse",
                "--training-subfile",
                "Full",
                "--var",
                "Psi",
                "--save-interpolator",
                os.path.join(self.test_dir, "before_validation.bin"),
                "--output",
                output,
            ],
        )
        self.assertNotEqual(result.exit_code, 0)
        self.assertIn("Validation times [5.0] are outside", result.output)
        self.assertFalse(os.path.exists(output))

        self.assertTrue(
            os.path.isfile(os.path.join(self.test_dir, "before_validation.bin"))
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
