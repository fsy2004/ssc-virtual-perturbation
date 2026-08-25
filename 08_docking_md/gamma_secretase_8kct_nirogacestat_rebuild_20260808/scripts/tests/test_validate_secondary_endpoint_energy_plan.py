import copy
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[2]
SCRIPT = ROOT / "scripts" / "validate_secondary_endpoint_energy_plan.py"
PLAN = ROOT / "config" / "secondary_endpoint_energy_plan_v1.json"


class SecondaryEndpointPlanTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = None
        cls.valid = json.loads(PLAN.read_text(encoding="utf-8"))
        if not SCRIPT.is_file():
            return
        spec = importlib.util.spec_from_file_location("endpoint_plan", SCRIPT)
        cls.module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(cls.module)

    def setUp(self):
        if self._testMethodName != "test_implementation_exists" and self.module is None:
            self.skipTest("implementation not created yet")

    def test_implementation_exists(self):
        self.assertTrue(SCRIPT.is_file(), f"missing implementation: {SCRIPT.name}")

    def write_plan(self, payload):
        handle = tempfile.NamedTemporaryFile(
            mode="w", suffix=".json", encoding="utf-8", delete=False
        )
        with handle:
            json.dump(payload, handle)
        self.addCleanup(Path(handle.name).unlink, missing_ok=True)
        return Path(handle.name)

    def validate_mutation(self, mutate):
        payload = copy.deepcopy(self.valid)
        mutate(payload)
        return self.module.validate_plan(self.write_plan(payload), bind_frozen_hash=False)

    def test_accepts_the_frozen_plan_and_reports_hash(self):
        result = self.module.validate_plan(PLAN)
        self.assertEqual(result["status"], "pass")
        self.assertEqual(result["analysis_id"], "o6u_secondary_endpoint_energy_v2_20260822")
        self.assertEqual(len(result["plan_sha256"]), 64)
        self.assertEqual(result["models"], ["PB_membrane_indi4"])
        self.assertEqual(result["withdrawn_sensitivity_models_not_run"], ["PB_membrane_indi1", "GB_OBC2", "GB_Neck2"])

    def test_rejects_release_hash_drift(self):
        with self.assertRaisesRegex(ValueError, "release_archive_sha256"):
            self.validate_mutation(lambda p: p.__setitem__("release_archive_sha256", "0" * 64))

    def test_rejects_missing_replica(self):
        with self.assertRaisesRegex(ValueError, "realizations"):
            self.validate_mutation(lambda p: p.__setitem__("realizations", ["rep01", "rep02"]))

    def test_rejects_changed_frame_count(self):
        with self.assertRaisesRegex(ValueError, "frames_per_realization"):
            self.validate_mutation(lambda p: p["sampling"].__setitem__("frames_per_realization", 299))

    def test_rejects_enabled_entropy(self):
        with self.assertRaisesRegex(ValueError, "entropy"):
            self.validate_mutation(lambda p: p["entropy"].__setitem__("normal_mode", True))

    def test_rejects_unapproved_decomposition_selection(self):
        with self.assertRaisesRegex(ValueError, "decomposition.selection"):
            self.validate_mutation(lambda p: p["decomposition"].__setitem__("selection", "top_20_by_energy"))

    def test_rejects_affinity_claim(self):
        with self.assertRaisesRegex(ValueError, "affinity_potency_efficacy_claim"):
            self.validate_mutation(lambda p: p["inference"].__setitem__("affinity_potency_efficacy_claim", True))

    def test_rejects_formal_execution_during_production(self):
        with self.assertRaisesRegex(ValueError, "concurrent_with_production"):
            self.validate_mutation(lambda p: p["execution_gate"].__setitem__("concurrent_with_production", True))


if __name__ == "__main__":
    unittest.main()
