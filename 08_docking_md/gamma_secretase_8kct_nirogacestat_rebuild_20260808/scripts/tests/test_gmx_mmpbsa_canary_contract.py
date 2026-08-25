import importlib.util
import math
import unittest
from pathlib import Path, PurePosixPath


ROOT = Path(__file__).parents[2]
SCRIPT = ROOT / "scripts" / "run_gmx_mmpbsa_canary.py"


class GmxMmpbsaCanaryContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = None
        if SCRIPT.is_file():
            spec = importlib.util.spec_from_file_location("endpoint_canary", SCRIPT)
            cls.module = importlib.util.module_from_spec(spec)
            assert spec.loader is not None
            spec.loader.exec_module(cls.module)

    def setUp(self):
        if self._testMethodName != "test_implementation_exists" and self.module is None:
            self.skipTest("implementation not created yet")

    def test_implementation_exists(self):
        self.assertTrue(SCRIPT.is_file(), f"missing implementation: {SCRIPT.name}")

    def test_canary_frames_are_prespecified_and_span_the_window(self):
        selected = self.module.prespecified_canary_frames()
        self.assertEqual([row["output_index_zero_based"] for row in selected], [0, 150, 299])
        self.assertEqual([row["target_time_ns"] for row in selected], [200.5, 350.5, 499.5])

    def test_model_inputs_match_frozen_scientific_settings(self):
        models = self.module.render_model_inputs(mthick_angstrom=38.75, frame_count=3)
        self.assertEqual(set(models), {"PB_membrane_indi4"})
        self.assertIn("memopt=1", models["PB_membrane_indi4"])
        self.assertIn("emem=7.0", models["PB_membrane_indi4"])
        self.assertIn("indi=4.0", models["PB_membrane_indi4"])
        self.assertIn("mthick=38.750000", models["PB_membrane_indi4"])
        self.assertIn("mctrdz=0.0", models["PB_membrane_indi4"])
        self.assertIn("fillratio=1.25", models["PB_membrane_indi4"])
        self.assertIn("solvopt=2", models["PB_membrane_indi4"])
        self.assertIn("eneopt=1", models["PB_membrane_indi4"])
        self.assertIn("cutnb=99.0", models["PB_membrane_indi4"])
        for text in models.values():
            self.assertIn("startframe=1", text)
            self.assertIn("endframe=3", text)
            self.assertIn("idecomp=2", text)
            self.assertIn('print_res="B/261,268,272,282,287,380-381,431-432,502"', text)
            self.assertNotIn("entropy", text.lower())
            self.assertNotIn("nmode", text.lower())

    def test_toolchain_contract_is_exact_and_cpu_only(self):
        record = {
            "gmx_mmpbsa": "1.6.5",
            "gmx_mmpbsa_git_commit": self.module.GMX_MMPBSA_COMMIT,
            "python": "3.11.8",
            "ambertools": "23.3",
            "gromacs": "2023.4",
            "openmpi": "4.1.6",
            "mpi4py": "4.0.1",
            "numpy": "1.26.4",
            "pandas": "1.5.3",
            "matplotlib": "3.7.3",
            "seaborn": "0.11.2",
            "scipy": "1.14.1",
            "tqdm": "4.67.1",
            "parmed": "4.3.0",
            "gpu_required": False,
            "executables": {name: {"sha256": "a" * 64} for name in self.module.REQUIRED_EXECUTABLES},
        }
        self.assertEqual(self.module.validate_toolchain_record(record)["status"], "pass")
        record["gromacs"] = "2025.2"
        with self.assertRaisesRegex(ValueError, "gromacs"):
            self.module.validate_toolchain_record(record)

    def test_canary_report_requires_all_models_finite_and_charge_bound(self):
        models = {}
        for name in self.module.MODEL_NAMES:
            models[name] = {
                "frame_count": 3,
                "finite_components": True,
                "component_values": [-10.0, 2.0, 4.0],
                "generated_topology_sha256": {
                    "complex": "a" * 64,
                    "receptor": "b" * 64,
                    "ligand": "c" * 64,
                },
                "decomposition_status": "pass",
                "decomposition_output_sha256": "d" * 64,
                "decomposition_residues_in_frozen_order": list(
                    self.module.FROZEN_DECOMP_RESIDUES
                ),
            }
        report = {
            "models": models,
            "complex_atom_count": 5000,
            "receptor_atom_count": 4924,
            "ligand_atom_count": 76,
            "ligand_total_charge": 0.0,
            "ligand_partial_charges_all_zero": False,
            "toolchain_status": "pass",
        }
        self.assertEqual(self.module.validate_canary_report(report)["status"], "pass")
        report["models"]["PB_membrane_indi4"]["component_values"][0] = math.nan
        with self.assertRaisesRegex(ValueError, "finite"):
            self.module.validate_canary_report(report)

    def test_canary_rejects_missing_decomposition(self):
        report = self.module.example_passing_canary_report()
        report["models"]["PB_membrane_indi4"]["decomposition_status"] = "missing"
        with self.assertRaisesRegex(ValueError, "decomposition"):
            self.module.validate_canary_report(report)

    def test_canary_rejects_decomposition_mapping_drift(self):
        report = self.module.example_passing_canary_report()
        report["models"]["PB_membrane_indi4"][
            "decomposition_residues_in_frozen_order"
        ][0] = "R:B:TYR:999"
        with self.assertRaisesRegex(ValueError, "residue mapping"):
            self.module.validate_canary_report(report)

    def test_canary_rejects_zero_partial_charges(self):
        report = self.module.example_passing_canary_report()
        report["ligand_partial_charges_all_zero"] = True
        with self.assertRaisesRegex(ValueError, "partial charges"):
            self.module.validate_canary_report(report)

    def test_canary_command_uses_three_mpi_ranks_and_frozen_groups(self):
        command = self.module.build_gmx_mmpbsa_command(
            prefix=PurePosixPath("/opt/o6u-gmxmmpbsa"),
            input_file=Path("model.in"),
            structure=Path("system.gro"),
            trajectory=Path("canary.xtc"),
            index=Path("groups.ndx"),
            topology=Path("topol.top"),
            reference=Path("complex.pdb"),
            final_text=Path("result.dat"),
            final_csv=Path("result.csv"),
            decomp_text=Path("decomp.dat"),
            decomp_csv=Path("decomp.csv"),
        )
        self.assertEqual(command[:3], ["/opt/o6u-gmxmmpbsa/bin/mpirun", "-np", "3"])
        self.assertIn("MPI", command)
        self.assertEqual(command[command.index("-cg") + 1:command.index("-cg") + 3], ["0", "1"])
        self.assertEqual(command[command.index("-do") + 1], "decomp.dat")
        self.assertEqual(command[command.index("-deo") + 1], "decomp.csv")
        self.assertNotIn("shell=True", command)

    def test_result_parser_requires_finite_components(self):
        values = self.module.parse_result_components("VDWAALS  -12.5\nDELTA TOTAL  3.25\n")
        self.assertEqual(values, [-12.5, 3.25])
        with self.assertRaisesRegex(ValueError, "finite"):
            self.module.parse_result_components("DELTA TOTAL NaN\n")

    def test_generated_topology_classifier_requires_one_each(self):
        paths = [
            Path("_GMXMMPBSA_COM.prmtop"),
            Path("_GMXMMPBSA_REC.prmtop"),
            Path("_GMXMMPBSA_LIG.prmtop"),
        ]
        result = self.module.classify_generated_topologies(paths)
        self.assertEqual(result["complex"].name, "_GMXMMPBSA_COM.prmtop")
        with self.assertRaisesRegex(ValueError, "ligand"):
            self.module.classify_generated_topologies(paths[:2])

    def test_gnu_time_peak_rss_is_converted_from_kib_to_bytes_per_rank(self):
        text = "Maximum resident set size (kbytes): 123456\n"
        self.assertEqual(self.module.parse_peak_rss_bytes(text, mpi_ranks=3), 123456 * 1024)
        with self.assertRaisesRegex(ValueError, "resident"):
            self.module.parse_peak_rss_bytes("Elapsed: 1.0", mpi_ranks=3)


if __name__ == "__main__":
    unittest.main()
