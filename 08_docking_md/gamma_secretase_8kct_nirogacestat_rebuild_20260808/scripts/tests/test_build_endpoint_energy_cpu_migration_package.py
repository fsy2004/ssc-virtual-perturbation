import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[2]
SCRIPT = ROOT / "scripts" / "build_endpoint_energy_cpu_migration_package.py"


class CpuMigrationPackageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = None
        if SCRIPT.is_file():
            spec = importlib.util.spec_from_file_location("migration_package", SCRIPT)
            cls.module = importlib.util.module_from_spec(spec)
            assert spec.loader is not None
            spec.loader.exec_module(cls.module)

    def setUp(self):
        if self._testMethodName != "test_implementation_exists" and self.module is None:
            self.skipTest("implementation not created yet")

    def test_implementation_exists(self):
        self.assertTrue(SCRIPT.is_file(), f"missing implementation: {SCRIPT.name}")

    def test_replica_mapping_requires_all_three_once(self):
        result = self.module.parse_replica_mapping([
            "rep01=/a", "rep02=/b", "rep03=/c"
        ])
        self.assertEqual(set(result), {"rep01", "rep02", "rep03"})
        with self.assertRaisesRegex(ValueError, "rep03"):
            self.module.parse_replica_mapping(["rep01=/a", "rep02=/b"])
        with self.assertRaisesRegex(ValueError, "duplicate"):
            self.module.parse_replica_mapping(["rep01=/a", "rep01=/b", "rep02=/c", "rep03=/d"])

    def test_topology_relative_path_rejects_escape(self):
        self.assertEqual(self.module.safe_relative_path("topology/topol.top"), Path("topology/topol.top"))
        with self.assertRaisesRegex(ValueError, "relative"):
            self.module.safe_relative_path("../topol.top")
        with self.assertRaisesRegex(ValueError, "relative"):
            self.module.safe_relative_path("/root/topol.top")

    def test_manifest_records_are_sorted_and_hashed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "b.txt").write_text("b", encoding="ascii")
            (root / "a.txt").write_text("a", encoding="ascii")
            records = self.module.manifest_file_records(root)
            self.assertEqual([row["path"] for row in records], ["a.txt", "b.txt"])
            self.assertTrue(all(len(row["sha256"]) == 64 for row in records))

    def test_completion_gate_requires_runners_stopped(self):
        passing = {
            "status": "pass",
            "all_three_500ns_complete": True,
            "all_required_gates_passed": True,
            "production_runners_active": False,
        }
        self.assertEqual(self.module.validate_completion_gate(passing)["status"], "pass")
        passing["production_runners_active"] = True
        with self.assertRaisesRegex(ValueError, "runner"):
            self.module.validate_completion_gate(passing)


if __name__ == "__main__":
    unittest.main()
