import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[2]
SCRIPT = ROOT / "scripts" / "plan_secondary_endpoint_resources.py"
GIB = 1024 ** 3
PEAKS = {
    "PB_membrane_indi4": 2 * GIB,
}


class ResourcePlanTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = None
        if SCRIPT.is_file():
            spec = importlib.util.spec_from_file_location("resource_plan", SCRIPT)
            cls.module = importlib.util.module_from_spec(spec)
            assert spec.loader is not None
            spec.loader.exec_module(cls.module)

    def setUp(self):
        if self._testMethodName != "test_implementation_exists" and self.module is None:
            self.skipTest("implementation not created yet")

    def test_implementation_exists(self):
        self.assertTrue(SCRIPT.is_file(), f"missing implementation: {SCRIPT.name}")

    def test_uses_nearly_all_208_cpus_in_three_replica_batches(self):
        result = self.module.build_resource_plan(
            logical_cpus=208,
            mem_total_bytes=754 * GIB,
            disk_free_bytes=270 * GIB,
            canary_peak_rss_bytes=PEAKS,
            projected_disk_bytes=40 * GIB,
        )
        self.assertEqual(result["cpu"]["reserved"], 2)
        self.assertEqual(result["models"]["PB_membrane_indi4"]["concurrent_jobs"], 3)
        self.assertEqual(result["models"]["PB_membrane_indi4"]["mpi_ranks_per_job"], 68)
        self.assertEqual(result["models"]["PB_membrane_indi4"]["total_mpi_ranks"], 204)
        self.assertEqual(
            result["models"]["PB_membrane_indi4"]["cpu_sets"],
            {"rep01": "2-69", "rep02": "70-137", "rep03": "138-205"},
        )
        self.assertEqual(result["environment"]["OMP_NUM_THREADS"], "1")
        self.assertFalse(result["gpu_required"])

    def test_memory_peak_limits_total_pb_ranks(self):
        peaks = dict(PEAKS)
        peaks["PB_membrane_indi4"] = 10 * GIB
        result = self.module.build_resource_plan(
            logical_cpus=208,
            mem_total_bytes=100 * GIB,
            disk_free_bytes=270 * GIB,
            canary_peak_rss_bytes=peaks,
            projected_disk_bytes=40 * GIB,
        )
        pb = result["models"]["PB_membrane_indi4"]
        self.assertEqual(pb["concurrent_jobs"], 3)
        self.assertEqual(pb["mpi_ranks_per_job"], 2)
        self.assertEqual(pb["total_mpi_ranks"], 6)
        self.assertEqual(pb["limiting_resource"], "memory")

    def test_rejects_missing_canary_peak(self):
        peaks = dict(PEAKS)
        peaks.pop("PB_membrane_indi4")
        with self.assertRaisesRegex(ValueError, "PB_membrane_indi4"):
            self.module.build_resource_plan(208, 754 * GIB, 270 * GIB, peaks, 40 * GIB)

    def test_rejects_insufficient_disk_reserve(self):
        with self.assertRaisesRegex(ValueError, "disk"):
            self.module.build_resource_plan(208, 754 * GIB, 80 * GIB, PEAKS, 40 * GIB)

    def test_rejects_resource_plan_with_fewer_than_three_pb_ranks(self):
        peaks = dict(PEAKS)
        peaks["PB_membrane_indi4"] = 50 * GIB
        with self.assertRaisesRegex(ValueError, "PB_membrane_indi4"):
            self.module.build_resource_plan(208, 100 * GIB, 270 * GIB, peaks, 40 * GIB)


if __name__ == "__main__":
    unittest.main()
