import io
import importlib.util
import unittest
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).parents[2]
SCRIPT = ROOT / "scripts" / "build_analysis_ndx.py"


class BuildAnalysisIndexTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location("build_analysis_ndx", SCRIPT)
        cls.module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(cls.module)

    def test_gromacs_index_file_is_one_based(self):
        handle = io.StringIO()
        self.module.write_group(handle, "System", [0, 1, 2])
        self.assertEqual(handle.getvalue(), "[ System ]\n1 2 3\n\n")

    def test_cluster_group_closes_every_seed_fragment(self):
        seed = SimpleNamespace(
            fragments=[
                SimpleNamespace(indices=[0, 1, 2, 3]),
                SimpleNamespace(indices=[8, 9]),
            ]
        )
        self.assertEqual(self.module.complete_fragment_indices(seed), [0, 1, 2, 3, 8, 9])


if __name__ == "__main__":
    unittest.main()
