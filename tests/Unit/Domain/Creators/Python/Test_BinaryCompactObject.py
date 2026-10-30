# Distributed under the MIT License.
# See LICENSE.txt for details.

import unittest

from spectre.Domain.Creators import BinaryCompactObject, DomainCreator3D
from spectre.Domain.Creators.TimeDependentOptions import (
    BinaryCompactObjectTimeDependentOptions,
    RotationMapOptions,
)


class TestBinaryCompactObject(unittest.TestCase):
    def test_construction(self):
        binary_compact_object_no_time_dependence = BinaryCompactObject(
            inner_radius_a=0.5,
            outer_radius_a=2.0,
            x_coord_a=5.0,
            excise_a=True,
            shell_log_map_strength_a=1.0,
            cube_log_map_strength_a=None,
            inner_radius_b=0.5,
            outer_radius_b=2.0,
            x_coord_b=-5.0,
            excise_b=True,
            shell_log_map_strength_b=1.0,
            cube_log_map_strength_b=1.5,
            center_of_mass_offset=[0.1, 0.2],
            envelope_radius=50.0,
            outer_radius=600.0,
            cube_scale=1.2,
            initial_refinement=1,
            initial_number_of_grid_points=5,
            use_equiangular_map=True,
            radial_partitioning_outer_shell=[],
            opening_angle_in_degrees=120.0,
            time_dependent_options=None,
        )
        rotation_map = RotationMapOptions([[0.0, 0.0, 0.0, 1.0]], 100.0)
        bco_time_dependent_options = BinaryCompactObjectTimeDependentOptions(
            0.0,
            None,
            rotation_map,
            None,
            None,
            None,
            None,
            None,
        )
        binary_compact_object_time_dependence = BinaryCompactObject(
            inner_radius_a=0.5,
            outer_radius_a=2.0,
            x_coord_a=5.0,
            excise_a=True,
            shell_log_map_strength_a=1.0,
            cube_log_map_strength_a=None,
            inner_radius_b=0.5,
            outer_radius_b=2.0,
            x_coord_b=-5.0,
            excise_b=True,
            shell_log_map_strength_b=1.0,
            cube_log_map_strength_b=1.5,
            center_of_mass_offset=[0.1, 0.2],
            envelope_radius=50.0,
            outer_radius=600.0,
            cube_scale=1.2,
            initial_refinement=1,
            initial_number_of_grid_points=5,
            use_equiangular_map=True,
            radial_partitioning_outer_shell=[],
            opening_angle_in_degrees=120.0,
            time_dependent_options=bco_time_dependent_options,
        )
        # A single spherical-harmonic shell block around each object. On these
        # blocks the scalar number of grid points is interpreted as
        # [n_radial, L_max] and only radial h-refinement is supported.
        binary_compact_object_spherical_harmonics = BinaryCompactObject(
            inner_radius_a=0.5,
            outer_radius_a=2.0,
            x_coord_a=5.0,
            excise_a=True,
            shell_log_map_strength_a=1.0,
            cube_log_map_strength_a=None,
            inner_radius_b=0.5,
            outer_radius_b=2.0,
            x_coord_b=-5.0,
            excise_b=True,
            shell_log_map_strength_b=None,
            cube_log_map_strength_b=None,
            center_of_mass_offset=[0.1, 0.2],
            envelope_radius=50.0,
            outer_radius=600.0,
            cube_scale=1.0,
            initial_refinement=0,
            initial_number_of_grid_points=7,
            use_equiangular_map=True,
            radial_partitioning_outer_shell=[],
            opening_angle_in_degrees=120.0,
            time_dependent_options=None,
            use_spherical_harmonics_a=True,
            use_spherical_harmonics_b=True,
        )
        self.assertIsInstance(
            binary_compact_object_no_time_dependence, DomainCreator3D
        )
        self.assertIsInstance(
            binary_compact_object_time_dependence, DomainCreator3D
        )
        self.assertNotEqual(
            binary_compact_object_time_dependence,
            binary_compact_object_no_time_dependence,
        )
        # 1 shell + 6 cubes per object, 10 envelope and 10 outer shell blocks
        block_names = binary_compact_object_spherical_harmonics.block_names()
        self.assertEqual(len(block_names), 34)
        self.assertEqual(block_names[0], "ObjectAShell0")
        self.assertEqual(block_names[7], "ObjectBShell0")
        self.assertEqual(
            binary_compact_object_spherical_harmonics.initial_extents()[0],
            [7, 8, 15],
        )
        domain = binary_compact_object_spherical_harmonics.create_domain()
        self.assertEqual(len(domain.blocks), 34)


if __name__ == "__main__":
    unittest.main(verbosity=2)
