import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[2]
SCRIPT = ROOT / "scripts" / "validate_secondary_endpoint_execution_defaults.py"
DEFAULTS = ROOT / "config" / "secondary_endpoint_energy_execution_defaults_v1.json"
SUPPLEMENT = ROOT / "SECONDARY_ENDPOINT_ENERGY_EXECUTION_SUPPLEMENT_20260820.md"
PLIP = ROOT / "docking_native_redock" / "figures" / "native_8kct_o6u" / "8KCT_O6U_native_contacts.interactions.normalized.json"


class EndpointExecutionDefaultsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = None
        if SCRIPT.is_file():
            spec = importlib.util.spec_from_file_location("execution_defaults", SCRIPT)
            cls.module = importlib.util.module_from_spec(spec)
            assert spec.loader is not None
            spec.loader.exec_module(cls.module)

    def setUp(self):
        if self._testMethodName != "test_implementation_exists" and self.module is None:
            self.skipTest("implementation not created yet")

    def test_implementation_exists(self):
        self.assertTrue(SCRIPT.is_file())

    def test_frozen_defaults_and_sources_pass(self):
        report = self.module.validate(DEFAULTS, SUPPLEMENT, PLIP)
        self.assertEqual(report["status"], "pass")
        self.assertEqual(report["decomposition_print_res"], "B/261,268,272,282,287,380-381,431-432,502")

    def test_defaults_hash_binding_rejects_modified_copy(self):
        payload = DEFAULTS.read_text(encoding="utf-8").replace('"idecomp": 2', '"idecomp": 1')
        temporary = DEFAULTS.parent / "_temporary_modified_defaults.json"
        temporary.write_text(payload, encoding="utf-8")
        try:
            with self.assertRaisesRegex(ValueError, "sha256"):
                self.module.validate(temporary, SUPPLEMENT, PLIP)
        finally:
            temporary.unlink()


if __name__ == "__main__":
    unittest.main()
