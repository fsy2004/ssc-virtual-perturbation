import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[2]
SCRIPT = ROOT / "scripts" / "seal_secondary_endpoint_all_three_gate.py"


class AllThreeGateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = None
        if SCRIPT.is_file():
            spec = importlib.util.spec_from_file_location("all_three_gate", SCRIPT)
            cls.module = importlib.util.module_from_spec(spec)
            assert spec.loader is not None
            spec.loader.exec_module(cls.module)

    def test_implementation_exists(self):
        self.assertTrue(SCRIPT.is_file())

    def test_completion_contract_requires_exact_500ns(self):
        if self.module is None:
            self.skipTest("implementation not created yet")
        payload = {
            "status": "pass", "report_type": "production_500ns_completion",
            "replica": "rep01", "final_step": 125000000, "final_time_ps": 500000.0,
            "production_release_sha256": self.module.PRODUCTION_RELEASE_SHA256,
        }
        self.assertEqual(self.module.validate_completion("rep01", payload)["status"], "pass")
        payload["final_time_ps"] = 499999.0
        with self.assertRaisesRegex(ValueError, "500"):
            self.module.validate_completion("rep01", payload)

    def test_analysis_evidence_supports_status_or_three_membrane_statuses(self):
        if self.module is None:
            self.skipTest("implementation not created yet")
        self.assertEqual(self.module.validate_analysis_evidence(
            "pbc", {"status": "pass", "maximum_absolute_difference_nm": 0.001}
        )["status"], "pass")
        membrane = {"technical_status": "pass", "sampling_status": "pass", "preproduction_status": "pass"}
        self.assertEqual(self.module.validate_analysis_evidence("membrane", membrane)["status"], "pass")
        with self.assertRaisesRegex(ValueError, "energy"):
            self.module.validate_analysis_evidence("energy", {"status": "fail"})


if __name__ == "__main__":
    unittest.main()
