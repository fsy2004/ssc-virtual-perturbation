from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "seal_primary_postprocessing_manifest.py"


def load_module():
    spec = importlib.util.spec_from_file_location("seal_primary_manifest_under_test", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class SealPrimaryPostprocessingManifestTests(unittest.TestCase):
    def setUp(self):
        self.module = load_module()
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        (self.root / "config").mkdir()
        self.draft = self.root / "config" / "primary_postprocessing_manifest.json"
        self.output = self.root / "config" / "primary_postprocessing_manifest.approved.json"
        self.protocol = self.root / "config" / "production_protocol.json"
        self.draft.write_text(json.dumps({
            "schema_version": "1.0",
            "approval_status": "draft_not_for_execution",
            "acceptance_gates": {"approval_status": "draft_not_for_execution"},
            "realizations": [
                {"realization_id": rid, "velocity_seed": "TODO", "topology": {},
                 "centered_system_trajectory": {}, "energy_edr": {}, "production_log": {}}
                for rid in ("rep01", "rep02", "rep03")
            ],
        }), encoding="utf-8")
        self.protocol.write_text(json.dumps({
            "realizations": [
                {"id": "rep01", "velocity_seed": 26081601},
                {"id": "rep02", "velocity_seed": 26081602},
                {"id": "rep03", "velocity_seed": 26081603},
            ]
        }), encoding="utf-8")
        for rid in ("rep01", "rep02", "rep03"):
            work = self.root / rid / "work"
            work.mkdir(parents=True)
            for name in ("production.tpr", "production.edr", "production.log"):
                (work / name).write_bytes(f"{rid}:{name}".encode())
            completion = {
                "status": "pass",
                "artifacts": {
                    name: {
                        "path": str(work / name),
                        "bytes": (work / name).stat().st_size,
                        "sha256": self.module.sha256(work / name),
                    }
                    for name in ("production.tpr", "production.edr", "production.log")
                },
            }
            (self.root / rid / "PRODUCTION_COMPLETION_500NS.json").write_text(
                json.dumps(completion), encoding="utf-8"
            )
            trajectory_dir = self.root / "analysis" / "trajectories" / "8kct_nirogacestat_native" / rid
            trajectory_dir.mkdir(parents=True)
            centered = trajectory_dir / "05_centered_reboxed.xtc"
            centered.write_bytes(f"{rid}:centered".encode())
            provenance = {
                "status": "pass_pending_scientific_qc_seal",
                "realization_id": rid,
                "production_tpr_sha256": completion["artifacts"]["production.tpr"]["sha256"],
                "retained_outputs": {
                    "center_and_rebox": {
                        "path": str(centered),
                        "bytes": centered.stat().st_size,
                        "sha256": self.module.sha256(centered),
                    }
                },
            }
            (trajectory_dir / "trajectory_provenance.pre_qc.json").write_text(
                json.dumps(provenance), encoding="utf-8"
            )

    def tearDown(self):
        self.temporary.cleanup()

    def test_binds_all_three_realizations_from_hash_verified_records(self):
        payload = self.module.bind_manifest(self.draft, self.root, self.protocol, self.output)
        self.assertEqual(payload["approval_status"], "approved_for_server_execution")
        self.assertEqual(
            payload["acceptance_gates"]["approval_status"],
            "approved_and_frozen_before_production",
        )
        self.assertEqual(
            [item["velocity_seed"] for item in payload["realizations"]],
            [26081601, 26081602, 26081603],
        )
        self.assertEqual(
            payload["realizations"][0]["centered_system_trajectory"]["path"],
            "analysis/trajectories/8kct_nirogacestat_native/rep01/05_centered_reboxed.xtc",
        )
        self.assertTrue(self.output.with_suffix(".json.sha256").is_file())

    def test_rejects_tpr_hash_disagreement_between_completion_and_provenance(self):
        path = self.root / "analysis" / "trajectories" / "8kct_nirogacestat_native" / "rep02" / "trajectory_provenance.pre_qc.json"
        payload = json.loads(path.read_text(encoding="utf-8"))
        payload["production_tpr_sha256"] = "0" * 64
        path.write_text(json.dumps(payload), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "TPR hash"):
            self.module.bind_manifest(self.draft, self.root, self.protocol, self.output)


if __name__ == "__main__":
    unittest.main()
