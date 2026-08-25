import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[2]
SCRIPT = ROOT / "scripts" / "prepare_secondary_endpoint_energy_inputs.py"


class PrepareEndpointInputsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = None
        if SCRIPT.is_file():
            spec = importlib.util.spec_from_file_location("prepare_endpoint", SCRIPT)
            cls.module = importlib.util.module_from_spec(spec)
            assert spec.loader is not None
            spec.loader.exec_module(cls.module)

    def setUp(self):
        if self._testMethodName != "test_implementation_exists" and self.module is None:
            self.skipTest("implementation not created yet")

    def test_implementation_exists(self):
        self.assertTrue(SCRIPT.is_file(), f"missing implementation: {SCRIPT.name}")

    def test_selects_exactly_300_midpoint_frames(self):
        times = [index / 10.0 for index in range(5001)]
        selected = self.module.select_midpoint_frames(times)
        self.assertEqual(len(selected), 300)
        self.assertEqual(selected[0]["target_time_ns"], 200.5)
        self.assertEqual(selected[0]["source_time_ns"], 200.5)
        self.assertEqual(selected[-1]["target_time_ns"], 499.5)
        self.assertEqual(selected[-1]["source_time_ns"], 499.5)
        self.assertEqual([sum(row["block_index"] == i for row in selected) for i in range(5)], [60] * 5)

    def test_midpoint_tie_selects_earlier_frame(self):
        times = [200.4, 200.6] + [201.5 + i for i in range(299)]
        selected = self.module.select_midpoint_frames(times)
        self.assertEqual(selected[0]["source_index_zero_based"], 0)
        self.assertEqual(selected[0]["source_time_ns"], 200.4)

    def test_rejects_sparse_times_that_reuse_a_frame(self):
        with self.assertRaisesRegex(ValueError, "unique"):
            self.module.select_midpoint_frames([200.0, 500.0])

    def test_frame_index_is_zero_based_and_hashable(self):
        selected = self.module.select_midpoint_frames([index / 10.0 for index in range(5001)])
        text = self.module.render_frame_index(selected)
        self.assertTrue(text.startswith("[ endpoint_midpoint_frames_zero_based ]\n"))
        self.assertIn("2005", text.splitlines()[1:])
        self.assertEqual(len(self.module.text_sha256(text)), 64)

    def test_output_names_separate_full_system_gro_from_complex_reference_pdb(self):
        names = self.module.preparation_output_names("rep01")
        self.assertEqual(names["structure"], "rep01_endpoint_structure.gro")
        self.assertEqual(names["reference"], "rep01_endpoint_complex_reference.pdb")
        self.assertEqual(names["trajectory"], "rep01_endpoint_300frames_midplane0.xtc")
        self.assertEqual(names["canary_trajectory"], "rep01_endpoint_canary_3frames_midplane0.xtc")

    def test_validates_frozen_receptor_ligand_groups(self):
        receptor = list(range(1, 1001))
        ligand = list(range(1001, 1077))
        result = self.module.validate_endpoint_groups(
            {"Receptor": receptor, "Ligand_O6U": ligand, "Complex": receptor + ligand},
            {"O6U": 1, "NAG": 18, "BMA": 3, "CLR": 3, "PC1": 2, "DSPC": 0},
        )
        self.assertEqual(result["ligand_atom_count"], 76)
        self.assertEqual(result["structural_lipid_residue_count"], 2)

    def test_rejects_wrong_o6u_atom_count(self):
        with self.assertRaisesRegex(ValueError, "76"):
            self.module.validate_endpoint_groups(
                {"Receptor": [1, 2], "Ligand_O6U": list(range(3, 78)), "Complex": list(range(1, 78))},
                {"O6U": 1, "NAG": 18, "BMA": 3, "CLR": 3, "PC1": 2, "DSPC": 0},
            )

    def test_rejects_component_count_drift(self):
        with self.assertRaisesRegex(ValueError, "CLR"):
            self.module.validate_endpoint_groups(
                {"Receptor": [1, 2], "Ligand_O6U": list(range(3, 79)), "Complex": list(range(1, 79))},
                {"O6U": 1, "NAG": 18, "BMA": 3, "CLR": 2, "PC1": 2, "DSPC": 0},
            )

    def test_selected_distance_invariance_requires_matching_times(self):
        raw = [(200.5, 0.25), (201.5, 0.30)]
        derived = [(200.5, 0.251), (201.5, 0.299)]
        result = self.module.compare_selected_distances(raw, derived, 0.01)
        self.assertEqual(result["status"], "pass")
        with self.assertRaisesRegex(ValueError, "times"):
            self.module.compare_selected_distances(raw, derived[::-1], 0.01)

    def test_refuses_to_overwrite_manifest(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "manifest.json"
            path.write_text("{}", encoding="utf-8")
            with self.assertRaises(FileExistsError):
                self.module.write_new_json(path, {"status": "pass"})


if __name__ == "__main__":
    unittest.main()
