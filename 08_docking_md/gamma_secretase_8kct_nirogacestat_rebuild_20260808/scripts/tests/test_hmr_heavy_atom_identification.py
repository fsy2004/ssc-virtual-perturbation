import importlib.util
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(ROOT / "scripts"))


def load_module(name: str, filename: str):
    path = ROOT / "scripts" / filename
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class HmrHeavyAtomIdentificationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.structure = load_module("structure_hmr_contract", "analyze_primary_structure_mdanalysis.py")
        cls.membrane = load_module("membrane_hmr_contract", "analyze_membrane_qc_mdanalysis.py")

    def test_hydrogen_identity_does_not_depend_on_hmr_mass(self):
        hydrogen = SimpleNamespace(name="H17", element="H", mass=3.024)
        carbon = SimpleNamespace(name="C1", element="C", mass=7.979)
        for module in (self.structure, self.membrane):
            self.assertTrue(module.atom_is_hydrogen(hydrogen))
            self.assertFalse(module.atom_is_hydrogen(carbon))

    def test_numeric_prefixed_hydrogen_name_is_recognized_without_element(self):
        hydrogen = SimpleNamespace(name="1H2", mass=3.024)
        for module in (self.structure, self.membrane):
            self.assertTrue(module.atom_is_hydrogen(hydrogen))


if __name__ == "__main__":
    unittest.main()
