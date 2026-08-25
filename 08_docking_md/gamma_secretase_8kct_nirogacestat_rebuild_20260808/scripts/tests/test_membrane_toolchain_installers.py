import unittest
from pathlib import Path


ROOT = Path(__file__).parents[2]
GORDER = ROOT / "scripts" / "install_gorder_1_5_0.sh"
FATSLIM = ROOT / "scripts" / "install_fatslim_0_2_2.sh"


class MembraneToolchainInstallerTests(unittest.TestCase):
    def test_gorder_download_is_retryable_and_bound_to_exact_commit_archive(self):
        text = GORDER.read_text(encoding="utf-8")
        self.assertIn("1beece37dc58a819be0a20b3ec691ef6cade365d", text)
        self.assertIn("RUST_VERSION=1.87.0", text)
        self.assertIn("--retry 10", text)
        self.assertIn("codeload.github.com/VachaLab/gorder/tar.gz/${GORDER_COMMIT}", text)
        self.assertNotIn("git clone", text)

    def test_fatslim_uses_working_mirror_and_exact_commit_archive(self):
        text = FATSLIM.read_text(encoding="utf-8")
        self.assertIn("ad79df027b62f10edf8e7d65298b13088d46f151", text)
        self.assertIn("--override-channels", text)
        self.assertIn("mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge", text)
        self.assertIn("--retry 10", text)
        self.assertIn("codeload.github.com/FATSLiM/fatslim/tar.gz/${FATSLIM_COMMIT}", text)
        self.assertIn("pip==23.3.2 numpy==1.21.6 cython==0.29.36 pytest==7.4.4", text)
        self.assertNotIn("git clone", text)


if __name__ == "__main__":
    unittest.main()
