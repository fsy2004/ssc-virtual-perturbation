import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[2]
CAPTURE = ROOT / "scripts" / "capture_gmx_mmpbsa_toolchain.py"
CPU_INSTALLER = ROOT / "scripts" / "install_gmx_mmpbsa_1_6_5_cpu.sh"
PREP_INSTALLER = ROOT / "scripts" / "install_endpoint_preprocess_env.sh"


class EndpointEnergyEnvironmentContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = None
        if CAPTURE.is_file():
            spec = importlib.util.spec_from_file_location("capture_toolchain", CAPTURE)
            cls.module = importlib.util.module_from_spec(spec)
            assert spec.loader is not None
            spec.loader.exec_module(cls.module)

    def test_required_environment_files_exist(self):
        self.assertTrue(CAPTURE.is_file())
        self.assertTrue(CPU_INSTALLER.is_file())
        self.assertTrue(PREP_INSTALLER.is_file())

    def test_cpu_installer_is_isolated_and_commit_pinned(self):
        text = CPU_INSTALLER.read_text(encoding="utf-8")
        self.assertIn("3.11.8", text)
        self.assertIn("ambertools=23.3", text)
        self.assertIn("gromacs=2023.4", text)
        self.assertIn("openmpi=4.1.6 c-compiler", text)
        self.assertIn("mpi4py==4.0.1", text)
        self.assertIn("tqdm==4.67.1", text)
        self.assertIn("64e994c71aaff315f3c82dd0852919aecb1ab62e", text)
        self.assertIn("refusing to overwrite", text)
        self.assertNotIn("cuda", text.lower())

    def test_cpu_installer_uses_retryable_exact_commit_archive_and_working_mirror(self):
        text = CPU_INSTALLER.read_text(encoding="utf-8")
        self.assertIn("--override-channels", text)
        self.assertIn("mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge", text)
        self.assertIn("--retry 10", text)
        self.assertIn("codeload.github.com/Valdes-Tresanco-MS/gmx_MMPBSA/tar.gz/${GMX_MMPBSA_COMMIT}", text)
        self.assertNotIn("git+https://", text)

    def test_preprocess_installer_pins_mdanalysis_without_gmx_mmpbsa(self):
        text = PREP_INSTALLER.read_text(encoding="utf-8")
        self.assertIn("python=3.11.8", text)
        self.assertIn("MDAnalysis==2.10.0", text)
        self.assertIn("numpy==1.26.4", text)
        self.assertNotIn("gmx_MMPBSA", text)
        self.assertIn("refusing to overwrite", text)

    def test_conda_package_index_rejects_duplicate_names(self):
        if self.module is None:
            self.skipTest("capture implementation not created yet")
        with self.assertRaisesRegex(ValueError, "duplicate"):
            self.module.index_conda_packages([
                {"name": "python", "version": "3.11.8"},
                {"name": "python", "version": "3.11.8"},
            ])

    def test_conda_package_index_preserves_build_and_channel(self):
        if self.module is None:
            self.skipTest("capture implementation not created yet")
        result = self.module.index_conda_packages([
            {"name": "gromacs", "version": "2023.4", "build_string": "nompi_0", "channel": "conda-forge"}
        ])
        self.assertEqual(result["gromacs"]["version"], "2023.4")
        self.assertEqual(result["gromacs"]["build_string"], "nompi_0")
        self.assertEqual(result["gromacs"]["channel"], "conda-forge")


if __name__ == "__main__":
    unittest.main()
