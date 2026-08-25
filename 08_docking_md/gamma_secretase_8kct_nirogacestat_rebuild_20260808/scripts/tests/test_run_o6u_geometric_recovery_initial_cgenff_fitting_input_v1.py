import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "run_o6u_geometric_recovery_initial_cgenff_fitting_input_v1.py"
SPEC = importlib.util.spec_from_file_location("recovery_fitting_input", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class RecoveryAdapterTests(unittest.TestCase):
    def test_changes_only_status(self):
        source = {
            "status": "pass_relaxed_mp2_torsion_geometric_recovery",
            "rotor_id": "ROT_C17_C15",
            "signed_step_index": 1,
            "final_energy_hartree": -100.0,
        }
        adapted = MODULE.adapt_recovery_report(source)
        self.assertEqual(adapted["status"], "pass_relaxed_mp2_torsion_scan_point")
        restored = dict(adapted)
        restored["status"] = source["status"]
        self.assertEqual(restored, source)

    def test_rejects_wrong_recovery_scope(self):
        with self.assertRaises(ValueError):
            MODULE.adapt_recovery_report({
                "status": "pass_relaxed_mp2_torsion_geometric_recovery",
                "rotor_id": "ROT_C15_N05",
                "signed_step_index": 1,
            })

    def test_adds_manifest_bound_step_when_recovery_report_omits_it(self):
        source = {
            "status": "pass_relaxed_mp2_torsion_geometric_recovery",
            "rotor_id": "ROT_C17_C15",
        }
        adapted = MODULE.adapt_recovery_report(source)
        self.assertEqual(adapted["signed_step_index"], 1)


if __name__ == "__main__":
    unittest.main()
