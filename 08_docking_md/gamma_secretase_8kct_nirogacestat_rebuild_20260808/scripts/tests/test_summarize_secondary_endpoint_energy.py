import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[2]
SCRIPT = ROOT / "scripts" / "summarize_secondary_endpoint_energy.py"


class EndpointEnergySummaryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = None
        if SCRIPT.is_file():
            spec = importlib.util.spec_from_file_location("endpoint_summary", SCRIPT)
            cls.module = importlib.util.module_from_spec(spec)
            assert spec.loader is not None
            spec.loader.exec_module(cls.module)

    def setUp(self):
        if self._testMethodName != "test_implementation_exists" and self.module is None:
            self.skipTest("implementation not created yet")

    def test_implementation_exists(self):
        self.assertTrue(SCRIPT.is_file(), f"missing implementation: {SCRIPT.name}")

    def test_parses_delta_section_only(self):
        text = (
            ",Complex\n,Frame #,VDWAALS,TOTAL\n,1,-1,9\n"
            ",Delta (Complex - Receptor - Ligand)\n"
            ",Frame #,VDWAALS,TOTAL\n,1,-10,-8\n,2,-12,-9\n"
        )
        rows = self.module.parse_delta_energy_csv(text)
        self.assertEqual(rows, [
            {"frame": 1, "VDWAALS": -10.0, "TOTAL": -8.0},
            {"frame": 2, "VDWAALS": -12.0, "TOTAL": -9.0},
        ])

    def test_fixed_blocks_are_five_by_sixty(self):
        rows = [{"frame": index + 1, "TOTAL": float(index)} for index in range(300)]
        blocks = self.module.fixed_block_means(rows)
        self.assertEqual(len(blocks), 5)
        self.assertEqual([row["frame_count"] for row in blocks], [60] * 5)
        self.assertEqual(blocks[0]["start_ns"], 200.0)
        self.assertEqual(blocks[-1]["end_ns"], 500.0)

    def test_hierarchical_bootstrap_is_seeded_and_has_no_p_value(self):
        replica_blocks = {
            "rep01": [float(i) for i in range(5)],
            "rep02": [float(i + 1) for i in range(5)],
            "rep03": [float(i + 2) for i in range(5)],
        }
        first = self.module.hierarchical_block_bootstrap(replica_blocks, seed=20260818, draws=1000)
        second = self.module.hierarchical_block_bootstrap(replica_blocks, seed=20260818, draws=1000)
        self.assertEqual(first, second)
        self.assertEqual(set(first), {"estimate", "ci95_low", "ci95_high", "draws", "seed"})
        self.assertLessEqual(first["ci95_low"], first["estimate"])
        self.assertGreaterEqual(first["ci95_high"], first["estimate"])

    def test_three_replica_descriptive_summary_reports_mean_and_sample_sd(self):
        result = self.module.three_replica_descriptive(
            {"rep01": 1.0, "rep02": 2.0, "rep03": 3.0}
        )
        self.assertEqual(result, {
            "n_realizations": 3,
            "mean": 2.0,
            "sample_sd": 1.0,
        })
        self.assertNotIn("p_value", result)

    def test_rejects_nonfinite_or_non_300_frame_formal_series(self):
        with self.assertRaisesRegex(ValueError, "300"):
            self.module.fixed_block_means([{"frame": 1, "TOTAL": 1.0}])

    def test_parses_fixed_per_residue_decomposition_without_ranking(self):
        header = (
            "DELTAS:\nTotal Decomposition Contribution (TDC)\n"
            "Frame #,Residue,Internal,van der Waals,Electrostatic,Polar Solvation,Non-Polar Solv.,TOTAL\n"
        )
        body = "\n".join(
            f"{frame},{residue},0,-1,-2,1,0,-2"
            for frame in range(1, 301)
            for residue in self.module.FROZEN_DECOMP_RESIDUES
        )
        rows = self.module.parse_decomposition_csv(header + body + "\n")
        self.assertEqual(len(rows), 300 * 10)
        self.assertEqual(
            self.module.validate_fixed_decomposition_rows(rows, frame_count=300)["status"],
            "pass",
        )
        self.assertEqual(
            list(self.module.fixed_decomposition_summary(rows)),
            list(self.module.FROZEN_DECOMP_RESIDUES),
        )
        first = self.module.fixed_decomposition_summary(rows)[self.module.FROZEN_DECOMP_RESIDUES[0]]
        self.assertEqual(first["frame_count"], 300)
        self.assertEqual(len(first["blocks"]), 5)
        self.assertNotIn("rank", first)
        self.assertNotIn("p_value", first)

        with self.assertRaisesRegex(ValueError, "frames"):
            self.module.validate_fixed_decomposition_rows(rows[:-1], frame_count=300)

    def test_decomposition_rejects_unfrozen_or_pairwise_residues(self):
        pairwise = (
            "DELTAS:\nTotal Decomposition Contribution (TDC)\n"
            "Frame #,Resid 1,Resid 2,Internal,van der Waals,Electrostatic,Polar Solvation,Non-Polar Solv.,TOTAL\n"
            "1,R:B:VAL:261,L:B:O6U:502,0,-1,-2,1,0,-2\n"
        )
        with self.assertRaisesRegex(ValueError, "per-residue"):
            self.module.parse_decomposition_csv(pairwise)

        unfrozen = (
            "DELTAS:\nTotal Decomposition Contribution (TDC)\n"
            "Frame #,Residue,Internal,van der Waals,Electrostatic,Polar Solvation,Non-Polar Solv.,TOTAL\n"
            "1,R:B:TYR:999,0,-1,-2,1,0,-2\n"
        )
        with self.assertRaisesRegex(ValueError, "frozen"):
            self.module.parse_decomposition_csv(unfrozen)

    def test_end_to_end_summary_writes_hash_bound_energy_and_decomposition_outputs(self):
        energy = (
            ",Delta (Complex - Receptor - Ligand)\n,Frame #,VDWAALS,TOTAL\n"
            + "\n".join(f",{frame},-1,-2" for frame in range(1, 301))
            + "\n"
        )
        decomposition = (
            "DELTAS:\nTotal Decomposition Contribution (TDC)\n"
            "Frame #,Residue,Internal,van der Waals,Electrostatic,Polar Solvation,Non-Polar Solv.,TOTAL\n"
            + "\n".join(
                f"{frame},{residue},0,-1,-2,1,0,-2"
                for frame in range(1, 301)
                for residue in self.module.FROZEN_DECOMP_RESIDUES
            )
            + "\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run_dir = root / "run"
            for model in self.module.MODELS:
                for replica in self.module.REPLICAS:
                    work = run_dir / model / replica
                    work.mkdir(parents=True)
                    (work / "FINAL_RESULTS_MMPBSA.csv").write_text(energy, encoding="utf-8")
                    (work / "FINAL_DECOMP_MMPBSA.csv").write_text(decomposition, encoding="utf-8")
            completion = root / "completion.json"
            completion.write_text(json.dumps({"status": "pass"}), encoding="utf-8")
            output = root / "summary"
            result = self.module.summarize(run_dir, completion, output)
            self.assertEqual(result["status"], "pass")
            manifest = json.loads((output / "SUMMARY_MANIFEST.json").read_text(encoding="utf-8"))
            self.assertEqual(len(manifest["outputs"]), 3)
            summary = json.loads(
                (output / "SECONDARY_ENDPOINT_ENERGY_SUMMARY.json").read_text(encoding="utf-8")
            )
            self.assertFalse(summary["decomposition_data_driven_ranking"])
            self.assertEqual(
                summary["inference"]["PB_membrane_indi4"]["TOTAL"]
                ["all_three_descriptive"]["sample_sd"],
                0.0,
            )


if __name__ == "__main__":
    unittest.main()
