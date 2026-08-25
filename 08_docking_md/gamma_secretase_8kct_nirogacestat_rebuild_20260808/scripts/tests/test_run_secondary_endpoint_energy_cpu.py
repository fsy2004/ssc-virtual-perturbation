import importlib.util
import unittest
from pathlib import Path, PurePosixPath


ROOT = Path(__file__).parents[2]
SCRIPT = ROOT / "scripts" / "run_secondary_endpoint_energy_cpu.py"


class FormalEndpointEnergyCpuRunnerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = None
        if SCRIPT.is_file():
            spec = importlib.util.spec_from_file_location("formal_cpu_runner", SCRIPT)
            cls.module = importlib.util.module_from_spec(spec)
            assert spec.loader is not None
            spec.loader.exec_module(cls.module)

    def setUp(self):
        if self._testMethodName != "test_implementation_exists" and self.module is None:
            self.skipTest("implementation not created yet")

    def test_implementation_exists(self):
        self.assertTrue(SCRIPT.is_file(), f"missing implementation: {SCRIPT.name}")

    def test_command_uses_cpu_partition_mpi_and_all_300_frames(self):
        command = self.module.build_formal_command(
            prefix=PurePosixPath("/opt/o6u"),
            cpu_set="2-69",
            mpi_ranks=68,
            input_file=Path("model.in"),
            structure=Path("system.gro"),
            trajectory=Path("trajectory.xtc"),
            index=Path("groups.ndx"),
            topology=Path("topol.top"),
            reference=Path("complex.pdb"),
            final_text=Path("result.dat"),
            final_csv=Path("result.csv"),
            decomp_text=Path("decomp.dat"),
            decomp_csv=Path("decomp.csv"),
        )
        self.assertEqual(command[:3], ["taskset", "-c", "2-69"])
        self.assertEqual(command[3:6], ["/opt/o6u/bin/mpirun", "--bind-to", "none"])
        self.assertEqual(command[command.index("-np") + 1], "68")
        self.assertEqual(command[command.index("-cg") + 1:command.index("-cg") + 3], ["0", "1"])
        self.assertEqual(command[command.index("-do") + 1], "decomp.dat")
        self.assertEqual(command[command.index("-deo") + 1], "decomp.csv")

    def test_resource_plan_requires_three_disjoint_cpu_sets(self):
        model = {
            "concurrent_jobs": 3,
            "mpi_ranks_per_job": 2,
            "total_mpi_ranks": 6,
            "replica_batch": ["rep01", "rep02", "rep03"],
            "cpu_sets": {"rep01": "2-3", "rep02": "4-5", "rep03": "6-7"},
        }
        self.assertEqual(self.module.validate_model_resource_plan(model)["status"], "pass")
        model["cpu_sets"]["rep03"] = "3-4"
        with self.assertRaisesRegex(ValueError, "overlap"):
            self.module.validate_model_resource_plan(model)

    def test_migration_manifest_validator_rejects_hash_drift(self):
        with self.assertRaisesRegex(ValueError, "hash"):
            self.module.validate_manifest_records(
                Path("/not/used"),
                [{"path": "missing.txt", "bytes": 1, "sha256": "a" * 64}],
            )


if __name__ == "__main__":
    unittest.main()
