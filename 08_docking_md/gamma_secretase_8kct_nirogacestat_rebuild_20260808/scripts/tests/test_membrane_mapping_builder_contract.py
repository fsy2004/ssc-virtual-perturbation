import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
BUILDER = ROOT / "scripts" / "build_membrane_mapping.py"
ANALYZER = ROOT / "scripts" / "analyze_membrane_qc_mdanalysis.py"
STEP5 = ROOT / "analysis_config_work" / "step5_input.pdb"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class MembraneMappingBuilderContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.builder = load_module("build_membrane_mapping_contract", BUILDER)
        cls.analyzer = load_module("analyze_membrane_qc_contract", ANALYZER)

    def test_draft_without_external_tools_is_internally_valid_and_placeholder_free(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "membrane_qc_mapping.json"
            old_output = self.builder.OUT
            self.builder.OUT = output
            try:
                self.assertEqual(self.builder.main(STEP5), 0)
            finally:
                self.builder.OUT = old_output

            mapping = json.loads(output.read_text(encoding="utf-8"))
            mapping["approval_status"] = "synthetic_self_test"
            universe = self.analyzer.mda.Universe(str(STEP5))
            groups = self.analyzer._validate_mapping(mapping, universe)

            self.assertGreater(len(groups["upper"]), 1)
            self.assertGreater(len(groups["lower"]), 1)
            self.assertFalse(self.analyzer.has_placeholder(mapping))
            self.assertIsNone(mapping["qc_gates"]["maximum_absolute_scd_adjacent_block_change"])
            self.assertIsNone(mapping["qc_gates"]["maximum_absolute_scd_first_last_change"])


if __name__ == "__main__":
    unittest.main()
