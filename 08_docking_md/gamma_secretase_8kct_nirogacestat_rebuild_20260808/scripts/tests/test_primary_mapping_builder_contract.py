from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

import MDAnalysis as mda


SCRIPTS = Path(__file__).resolve().parents[1]
ROOT = SCRIPTS.parent
sys.path.insert(0, str(SCRIPTS))

from analyze_primary_structure_mdanalysis import (  # noqa: E402
    _validate_mapping_record,
    topology_identity_sha256,
)


class PrimaryMappingBuilderContractTests(unittest.TestCase):
    def test_builder_output_matches_structural_analyzer_schema_and_identity(self):
        builder_path = SCRIPTS / "build_primary_mapping_records.py"
        spec = importlib.util.spec_from_file_location("primary_mapping_builder_under_test", builder_path)
        builder = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(builder)
        with tempfile.TemporaryDirectory() as tmpdir:
            builder.OUT = Path(tmpdir) / "mapping.json"
            self.assertEqual(builder.main(), 0)
            mapping = json.loads(builder.OUT.read_text(encoding="utf-8"))
        reference = mda.Universe(str(builder.REFERENCE))
        trajectory = mda.Universe(str(builder.STEP5))
        self.assertEqual(mapping["trajectory_atom_identity_sha256"], topology_identity_sha256(trajectory))
        mapping["approval_status"] = "approved"
        validated = _validate_mapping_record(mapping, reference, trajectory)
        self.assertGreaterEqual(len(validated["native_contacts"]), 1)
        self.assertGreaterEqual(len(validated["hydrogen_bonds"]), 1)


if __name__ == "__main__":
    unittest.main()
