import hashlib
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
SCRIPT = ROOT / "scripts" / "resume_primary_pbc_trajectories.py"


class ResumePrimaryPbcTrajectoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = None
        if SCRIPT.is_file():
            spec = importlib.util.spec_from_file_location("resume_primary_pbc", SCRIPT)
            cls.module = importlib.util.module_from_spec(spec)
            assert spec.loader is not None
            spec.loader.exec_module(cls.module)

    def setUp(self):
        if self._testMethodName != "test_implementation_exists" and self.module is None:
            self.skipTest("implementation not created yet")

    def test_implementation_exists(self):
        self.assertTrue(SCRIPT.is_file(), f"missing implementation: {SCRIPT.name}")

    def test_archives_only_partial_mindist_outputs_with_hashes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "analysis" / "rep01"
            output.mkdir(parents=True)
            retained = output / "07_fixed_200_500ns.xtc"
            retained.write_bytes(b"retained")
            raw = output / "09_raw_minimum_image_protein_O6U_heavy.xvg"
            processed = output / "10_processed_minimum_image_protein_O6U_heavy.xvg"
            command = output / "mindist_raw.command.json"
            raw.write_bytes(b"partial-raw")
            processed.write_bytes(b"partial-processed")
            command.write_bytes(b"partial-command")

            records = self.module.archive_partial_mindist_outputs(
                output, root / "audit" / "resume-001"
            )

            self.assertEqual(
                [record["name"] for record in records],
                [
                    "09_raw_minimum_image_protein_O6U_heavy.xvg",
                    "10_processed_minimum_image_protein_O6U_heavy.xvg",
                    "mindist_raw.command.json",
                ],
            )
            self.assertFalse(raw.exists())
            self.assertFalse(processed.exists())
            self.assertFalse(command.exists())
            self.assertEqual(retained.read_bytes(), b"retained")
            self.assertEqual(
                records[0]["sha256"], hashlib.sha256(b"partial-raw").hexdigest()
            )
            self.assertEqual(
                (root / "audit" / "resume-001" / raw.name).read_bytes(),
                b"partial-raw",
            )

    def test_refuses_to_archive_when_final_pbc_report_exists(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "analysis" / "rep01"
            output.mkdir(parents=True)
            (output / "11_pbc_distance_invariance.json").write_text("{}", encoding="utf-8")
            (output / "09_raw_minimum_image_protein_O6U_heavy.xvg").write_bytes(b"raw")

            with self.assertRaisesRegex(FileExistsError, "final PBC output"):
                self.module.archive_partial_mindist_outputs(output, root / "audit")

    def test_retained_inputs_require_exact_expected_files_and_no_final_provenance(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            expected = {
                "03_processed_first_frame.gro": b"gro",
                "05_centered_reboxed.xtc": b"centered",
                "06_fitted_analysis.xtc": b"fitted",
                "07_fixed_200_500ns.xtc": b"fixed",
            }
            for name, content in expected.items():
                (output / name).write_bytes(content)

            records = self.module.validate_retained_inputs(output)
            self.assertEqual(set(records), set(expected))
            self.assertEqual(records["05_centered_reboxed.xtc"]["bytes"], 8)

            (output / "trajectory_provenance.pre_qc.json").write_text("{}", encoding="utf-8")
            with self.assertRaisesRegex(FileExistsError, "final provenance"):
                self.module.validate_retained_inputs(output)


if __name__ == "__main__":
    unittest.main()
