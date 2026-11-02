# Distributed under the MIT License.
# See LICENSE.txt for details.

import logging
import os
from pathlib import Path
import shutil
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
import yaml

import spectre.IO.H5 as spectre_h5
from spectre.Informer import unit_test_build_path
from spectre.support.Logging import configure_logging
from spectre.tools.Status.ExecutableStatus import match_executable_status
from spectre.tools.Status.Status import (
    AVAILABLE_COLUMNS,
    fetch_status,
    get_executable_name,
    get_input_file,
    get_work_dir,
)


class TestExecutableStatus(unittest.TestCase):
    def setUp(self):
        self.input_file = {
            "Observers": {"ReductionFileName": "Reductions"},
            "EventsAndTriggersAtSlabs": [
                {
                    "Trigger": "Always",
                    "Events": [
                        {"ObserveTimeStep": {"SubfileName": "TimeSteps"}}
                    ],
                }
            ],
        }
        self.work_dir = os.path.join(
            unit_test_build_path(), "tools/ExecutableStatus"
        )
        shutil.rmtree(self.work_dir, ignore_errors=True)
        os.makedirs(self.work_dir, exist_ok=True)
        with spectre_h5.H5File(
            os.path.join(self.work_dir, "Reductions.h5"), "w"
        ) as open_h5_file:
            # Time data
            time_subfile = open_h5_file.insert_dat(
                "/TimeSteps",
                legend=[
                    "Time",
                    "Slab size",
                    "Minimum Walltime",
                    "Maximum Walltime",
                ],
                version=0,
            )
            time_subfile.append([0.0, 0.5, 0.5, 1.0])
            time_subfile.append([1.0, 0.5, 2.0, 3.0])
            time_subfile.append([2.0, 0.5, 3.0, 4.0])
            open_h5_file.close_current_object()
            # Control systems data
            rotation_subfile = open_h5_file.insert_dat(
                "/ControlSystems/Rotation/z",
                legend=["FunctionOfTime"],
                version=0,
            )
            rotation_subfile.append([0.0])
            rotation_subfile.append([np.pi])
            open_h5_file.close_current_object()
            # AH data
            aha_subfile = open_h5_file.insert_dat(
                "/ApparentHorizons/ControlSystemAhA_Centers",
                legend=[
                    "InertialCenter_x",
                    "InertialCenter_y",
                    "InertialCenter_z",
                ],
                version=0,
            )
            aha_subfile.append([1.0, 0.0, 0.0])
            open_h5_file.close_current_object()
            ahb_subfile = open_h5_file.insert_dat(
                "/ApparentHorizons/ControlSystemAhB_Centers",
                legend=[
                    "InertialCenter_x",
                    "InertialCenter_y",
                    "InertialCenter_z",
                ],
                version=0,
            )
            ahb_subfile.append([-1.0, 0.0, 0.0])
            open_h5_file.close_current_object()
            # Constraints
            constraints_subfile = open_h5_file.insert_dat(
                "/Norms",
                legend=[
                    "L2Norm(ConstraintEnergy)",
                    "L2Norm(PointwiseL2Norm(ThreeIndexConstraint))",
                ],
                version=0,
            )
            constraints_subfile.append([1.0e-3, 1.0e-4])
            open_h5_file.close_current_object()
            # Elliptic solver residuals
            residuals_subfile = open_h5_file.insert_dat(
                "/NewtonRaphsonResiduals",
                legend=["Iteration", "Residual"],
                version=0,
            )
            residuals_subfile.append([0.0, 1.0e-1])
            residuals_subfile.append([1.0, 1.0e-4])
            open_h5_file.close_current_object()
            residuals_subfile = open_h5_file.insert_dat(
                "/GmresResiduals",
                legend=["Iteration", "Residual"],
                version=0,
            )
            residuals_subfile.append([0.0, 1.0e-1])
            residuals_subfile.append([1.0, 1.0e-4])
            residuals_subfile.append([0.0, 1.0e-2])
            residuals_subfile.append([1.0, 1.0e-5])

    def tearDown(self):
        shutil.rmtree(self.work_dir, ignore_errors=True)

    def test_evolution_status(self):
        executable_status = match_executable_status("EvolveSomething")
        status = executable_status.status(self.input_file, self.work_dir)
        self.assertEqual(status, {"Time": 2.0, "Speed": 2400.0})
        self.assertEqual(executable_status.format("Time", 1.5), "1.5")
        self.assertEqual(executable_status.format("Speed", 1.2), "1.2")

    def test_evolve_bbh_status(self):
        executable_status = match_executable_status("EvolveGhBinaryBlackHole")
        status = executable_status.status(self.input_file, self.work_dir)
        self.assertEqual(status["Time"], 2.0)
        self.assertEqual(status["Speed"], 2400.0)
        self.assertEqual(status["Orbits"], 0.5)
        self.assertEqual(status["Separation"], 2.0)
        self.assertEqual(status["3-Index Constraint"], 1.0e-4)

    def test_evolve_single_bh_status(self):
        executable_status = match_executable_status("EvolveGhSingleBlackHole")
        status = executable_status.status(self.input_file, self.work_dir)
        self.assertEqual(status["Time"], 2.0)
        self.assertEqual(status["Speed"], 2400.0)
        self.assertEqual(status["Constraint Energy"], 1.0e-3)

    def test_elliptic_status(self):
        executable_status = match_executable_status("SolveSomething")
        status = executable_status.status(self.input_file, self.work_dir)
        self.assertEqual(status, {"Iteration": 1.0, "Residual": 1.0e-4})
        self.assertEqual(executable_status.format("Iteration", 1.0), "1")
        self.assertEqual(
            executable_status.format("Residual", 1.0e-4), "1.0e-04"
        )

    def test_xcts_status(self):
        executable_status = match_executable_status("SolveXcts")
        status = executable_status.status(self.input_file, self.work_dir)
        self.assertEqual(
            status,
            {
                "Nonlinear iteration": 1,
                "Nonlinear residual": 1.0e-4,
                "Linear iteration": 2,
                "Linear residual": 1.0e-5,
            },
        )
        self.assertEqual(
            executable_status.format("Nonlinear iteration", 1.0), "1"
        )
        self.assertEqual(executable_status.format("Linear iteration", 1.0), "1")
        self.assertEqual(
            executable_status.format("Nonlinear residual", 1.0e-4), "1.0e-04"
        )
        self.assertEqual(
            executable_status.format("Linear residual", 1.0e-4), "1.0e-04"
        )


class TestStatus(unittest.TestCase):
    def setUp(self):
        self.slurm_comment = (
            "SPECTRE_INPUT_FILE=path/to/input/file\n"
            "SPECTRE_EXECUTABLE=path/to/executable\n"
        )
        self.work_dir = os.path.join(unit_test_build_path(), "tools/Status")
        self.input_file_path = os.path.join(self.work_dir, "InputFile.yaml")
        os.makedirs(self.work_dir, exist_ok=True)
        with open(self.input_file_path, "w") as open_input_file:
            yaml.safe_dump_all([{"Executable": "MyExec"}, {}], open_input_file)

    def tearDown(self):
        shutil.rmtree(self.work_dir, ignore_errors=True)

    def test_get_input_file(self):
        self.assertEqual(
            get_input_file(self.slurm_comment, self.work_dir),
            os.path.join(self.work_dir, "path/to/input/file"),
        )
        self.assertIsNone(
            get_input_file("", os.path.join(self.work_dir, "nonexistent"))
        )
        self.assertEqual(
            get_input_file("", self.work_dir), self.input_file_path
        )

    def test_get_work_dir(self):
        submit_dir = str(Path(self.work_dir) / "submit")
        output_path = str(Path(self.work_dir) / "spectre.out")
        self.assertEqual(get_work_dir(submit_dir, ""), submit_dir)
        self.assertEqual(get_work_dir(submit_dir, output_path), submit_dir)
        context_path = Path(self.work_dir) / "SchedulerContext.yaml"
        with context_path.open("w") as context_file:
            yaml.safe_dump(
                {"run_dir": self.work_dir, "out_file": output_path},
                context_file,
            )
        self.assertEqual(get_work_dir(submit_dir, output_path), self.work_dir)
        self.assertEqual(get_work_dir(self.work_dir, output_path), self.work_dir)
        # A separate logs directory must not become the run directory just
        # because it happens to contain a context from another submission.
        for context in (
            {"run_dir": submit_dir, "out_file": output_path},
            {"run_dir": self.work_dir, "out_file": "another.out"},
            None,
        ):
            with context_path.open("w") as context_file:
                yaml.safe_dump(context, context_file)
            self.assertEqual(get_work_dir(submit_dir, output_path), submit_dir)

    def test_jobs_submitted_from_same_directory(self):
        jobs = []
        run_dirs = []
        for i in range(3):
            run_dir = Path(self.work_dir) / str(i) / "segments/0000_Inspiral"
            run_dir.mkdir(parents=True)
            run_dirs.append(str(run_dir))
            output_path = str(run_dir / "spectre.out")
            with (run_dir / "SchedulerContext.yaml").open("w") as context_file:
                yaml.safe_dump(
                    {
                        "run_dir": str(run_dir),
                        "out_file": output_path,
                        "input_file": "InputFile.yaml",
                    },
                    context_file,
                )
            shutil.copyfile(self.input_file_path, run_dir / "InputFile.yaml")
            job = {key: "" for key in AVAILABLE_COLUMNS if key != "SegmentId"}
            job.update(
                JobID=str(i),
                State="RUNNING",
                WorkDir=self.work_dir,
                StdOut=output_path,
                StdErr=output_path,
            )
            jobs.append(job)
        with patch(
            "spectre.tools.Status.Status.fetch_job_data",
            return_value=pd.DataFrame(jobs),
        ):
            status = fetch_status()
        self.assertEqual(set(status["WorkDir"]), set(run_dirs))
        self.assertEqual(list(status["ExecutableName"]), ["MyExec"] * 3)
        self.assertEqual(len(set(status["SegmentsDir"])), 3)

    def test_get_executable_name(self):
        self.assertEqual(
            get_executable_name(self.slurm_comment, self.input_file_path),
            "executable",
        )
        self.assertIsNone(get_executable_name("", None))
        self.assertEqual(
            get_executable_name("", self.input_file_path), "MyExec"
        )


if __name__ == "__main__":
    configure_logging(log_level=logging.DEBUG)
    unittest.main(verbosity=2)
