import importlib.util
import math
import unittest
from pathlib import Path

import numpy as np


SCRIPT = Path(__file__).parents[1] / "build_o6u_global_torsion_design_matrix_v1.py"
SPEC = importlib.util.spec_from_file_location("global_design", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class GlobalDesignMatrixTests(unittest.TestCase):
    def test_dihedral_angle(self):
        coordinates = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 1.0],
        ])
        angle = MODULE.torsion_angle_deg(coordinates, (1, 2, 3, 4))
        self.assertTrue(math.isclose(abs(angle), 90.0, abs_tol=1e-10))

    def test_basis_is_reference_subtracted(self):
        value = MODULE.periodic_basis(60.0, 60.0, 2, 180.0)
        self.assertTrue(math.isclose(value, 0.0, abs_tol=1e-12))

    def test_parameter_line_parser(self):
        term = MODULE.parse_parameter_line("X CG2R51 NG2S1 X 3.87 1 180.0 ! penalty=28")
        self.assertEqual(term["periodicity"], 1)
        self.assertEqual(term["phase_deg"], 180.0)


if __name__ == "__main__":
    unittest.main()
