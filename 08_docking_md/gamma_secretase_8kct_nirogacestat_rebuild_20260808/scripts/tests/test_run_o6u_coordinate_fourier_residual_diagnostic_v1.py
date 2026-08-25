import importlib.util
import math
import unittest
from pathlib import Path

import numpy as np


SCRIPT = Path(__file__).parents[1] / "run_o6u_coordinate_fourier_residual_diagnostic_v1.py"
SPEC = importlib.util.spec_from_file_location("fourier_diagnostic", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class FourierResidualDiagnosticTests(unittest.TestCase):
    def test_design_is_zero_at_reference(self):
        matrix = MODULE.design_matrix(np.array([30.0]), 30.0, 2)
        np.testing.assert_allclose(matrix, np.zeros((1, 4)), atol=1e-12)

    def test_recovers_single_harmonic_residual(self):
        angles = np.array([-30.0, -15.0, 15.0, 30.0])
        reference = 0.0
        matrix = MODULE.design_matrix(angles, reference, 1)
        expected_coefficients = np.array([2.0, -0.5])
        residual = matrix @ expected_coefficients
        fit = MODULE.fit_order(angles, residual, reference, 1)
        np.testing.assert_allclose(fit["coefficients"], expected_coefficients, atol=1e-10)
        self.assertLess(fit["training_rmse_kcal_mol"], 1e-10)

    def test_phase_amplitude_round_trip(self):
        item = MODULE.coefficient_terms([3.0, 4.0])[0]
        self.assertTrue(math.isclose(item["amplitude_kcal_mol"], 5.0))
        self.assertTrue(math.isclose(item["phase_deg"], math.degrees(math.atan2(4.0, 3.0))))


if __name__ == "__main__":
    unittest.main()
