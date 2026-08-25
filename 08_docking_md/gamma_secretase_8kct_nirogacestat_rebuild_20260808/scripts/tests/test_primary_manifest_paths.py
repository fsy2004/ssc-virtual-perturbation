from __future__ import annotations

import tempfile
import unittest
import importlib.util
import json
from pathlib import Path

from scripts import primary_postprocessing_common as common

BUILDER = Path(__file__).resolve().parents[1] / "build_primary_manifest.py"


class PrimaryManifestPathTests(unittest.TestCase):
    def test_config_manifest_resolves_records_from_package_root(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest = root / "config" / "primary_postprocessing_manifest.json"
            manifest.parent.mkdir()
            self.assertEqual(common.manifest_package_root(manifest), root.resolve())

    def test_non_config_manifest_uses_its_own_directory(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest = root / "synthetic_manifest.json"
            self.assertEqual(common.manifest_package_root(manifest), root.resolve())

    def test_builder_binds_the_same_protonated_reference_used_by_mapping(self):
        spec = importlib.util.spec_from_file_location("primary_manifest_builder_under_test", BUILDER)
        module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as tmpdir:
            module.OUT = Path(tmpdir) / "primary.json"
            self.assertEqual(module.main(), 0)
            payload = json.loads(module.OUT.read_text(encoding="utf-8"))
        expected_path = "docking_native_redock/plip_native/8KCT_protonated.pdb"
        self.assertEqual(payload["reference"]["topology"]["path"], expected_path)
        self.assertEqual(payload["reference"]["coordinates"]["path"], expected_path)
        self.assertEqual(
            payload["reference"]["topology"]["sha256"],
            common.sha256_file(module.REFERENCE),
        )
        self.assertEqual(
            payload["acceptance_gates"]["source_record"]["sha256"],
            common.sha256_file(module.FREEZE),
        )


if __name__ == "__main__":
    unittest.main()
