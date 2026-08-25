import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "run_o6u_adaptive_scan_initial_cgenff_fitting_input_v1.py"
SPEC = importlib.util.spec_from_file_location("adaptive_fitting_input", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class AdaptiveFittingInputMetadataTests(unittest.TestCase):
    def test_reframes_numeric_result_without_changing_comparison(self):
        source = {
            "status": "pass_tier_b_initial_cgenff_mismatch_prescreen",
            "scope": {"rotor_id": "ROT_C15_N05", "signed_step_index": -2},
            "comparison": {
                "qm_delta_kcal_mol": 1.25,
                "initial_cgenff_mm_delta_kcal_mol": 2.5,
                "signed_mm_minus_qm_delta_kcal_mol": 1.25,
                "absolute_mm_minus_qm_delta_kcal_mol": 1.25,
            },
        }
        result = MODULE.reframe_report(
            source,
            rotor_id="ROT_C15_N05",
            signed_step_index=-2,
            engine_path=Path("engine.py"),
            engine_sha="a" * 64,
            engine_report_sha="b" * 64,
        )
        self.assertEqual(result["comparison"], source["comparison"])
        self.assertEqual(result["status"], "pass_adaptive_scan_initial_cgenff_fitting_input")
        self.assertFalse(result["method"]["parameter_mutation"])
        self.assertIn("not a fitted parameter set", result["interpretation_boundary"])

    def test_rejects_rotor_outside_closed_adaptive_dataset(self):
        with self.assertRaises(ValueError):
            MODULE.validate_scope("ROT_C20_C18", 1)

    def test_rejects_reference_step(self):
        with self.assertRaises(ValueError):
            MODULE.validate_scope("ROT_C09_N04", 0)


if __name__ == "__main__":
    unittest.main()
