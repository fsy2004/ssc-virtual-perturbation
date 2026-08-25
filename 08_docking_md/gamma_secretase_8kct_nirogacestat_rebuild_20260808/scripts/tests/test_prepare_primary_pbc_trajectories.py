import hashlib
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
SCRIPT = ROOT / "scripts" / "prepare_primary_pbc_trajectories.py"


class PreparePrimaryPbcTrajectoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = None
        if SCRIPT.is_file():
            spec = importlib.util.spec_from_file_location("prepare_primary_pbc", SCRIPT)
            cls.module = importlib.util.module_from_spec(spec)
            assert spec.loader is not None
            spec.loader.exec_module(cls.module)

    def setUp(self):
        if self._testMethodName != "test_implementation_exists" and self.module is None:
            self.skipTest("implementation not created yet")

    def test_implementation_exists(self):
        self.assertTrue(SCRIPT.is_file(), f"missing implementation: {SCRIPT.name}")

    def test_completion_gate_binds_exact_endpoint_release_and_live_sizes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            work = root / "rep01" / "work"
            work.mkdir(parents=True)
            artifacts = {}
            for name, data in {
                "production.tpr": b"tpr",
                "production.xtc": b"xtc",
                "production.edr": b"edr",
                "production.log": b"log",
                "production.gro": b"gro",
                "production.cpt": b"cpt",
            }.items():
                path = work / name
                path.write_bytes(data)
                artifacts[name] = {
                    "bytes": len(data),
                    "sha256": hashlib.sha256(data).hexdigest(),
                }
            report = {
                "schema_version": "1.0",
                "report_type": "production_500ns_completion",
                "status": "pass",
                "replica": "rep01",
                "final_step": 125000000,
                "final_time_ps": 500000.0,
                "production_release_sha256": self.module.PRODUCTION_RELEASE_SHA256,
                "production_tpr_sha256": artifacts["production.tpr"]["sha256"],
                "checks": {"gmx_readability": {key: "pass" for key in ("cpt", "edr", "gro", "xtc")}},
                "artifacts": artifacts,
            }
            report_path = root / "rep01" / "PRODUCTION_COMPLETION_500NS.json"
            report_path.write_text(json.dumps(report), encoding="utf-8")
            self.module.TPR_SHA256["rep01"] = artifacts["production.tpr"]["sha256"]
            result = self.module.validate_completion_gate(root, "rep01", report_path)
            self.assertEqual(result["status"], "pass")

            (work / "production.xtc").write_bytes(b"changed-size")
            with self.assertRaisesRegex(ValueError, "byte count"):
                self.module.validate_completion_gate(root, "rep01", report_path)

    def test_pipeline_order_is_fixed_and_intermediates_are_explicit(self):
        plan = self.module.pipeline_contract()
        self.assertEqual(
            [item["mode"] for item in plan],
            [
                "whole",
                "cluster_complex_if_required",
                "processed_first_frame_reference",
                "nojump",
                "center_and_rebox",
                "fit_analysis_selection",
                "fixed_window_extract",
            ],
        )
        self.assertEqual(plan[4]["retain"], True)
        self.assertEqual(plan[5]["retain"], True)
        self.assertEqual(plan[6]["retain"], True)
        self.assertEqual([item["mode"] for item in plan if not item["retain"]], ["whole", "cluster_complex_if_required", "nojump"])


if __name__ == "__main__":
    unittest.main()
