import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "map_o6u_rotor_to_charmm_dihedrals_v1.py"
SPEC = importlib.util.spec_from_file_location("rotor_map", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


PSF = """PSF

       1 !NTITLE
 REMARKS test
       4 !NATOM
       1 LIG  1 O6U A1   CT1  0.0 12.0
       2 LIG  1 O6U C15  CG2R 0.0 12.0
       3 LIG  1 O6U N05  NG2S 0.0 14.0
       4 LIG  1 O6U A4   HGA  0.0  1.0
       1 !NPHI: dihedrals
       1       2       3       4
"""


class RotorMappingTests(unittest.TestCase):
    def test_maps_central_bond_and_parameter(self):
        atoms, dihedrals = MODULE.parse_psf_text(PSF)
        mapped = MODULE.map_rotor(atoms, dihedrals, "C15", "N05")
        self.assertEqual(mapped[0]["atom_types"], ["CT1", "CG2R", "NG2S", "HGA"])
        matches = MODULE.match_parameter_lines(
            mapped[0]["atom_types"],
            ["X  CG2R  NG2S  X  2.000  2  180.0 ! wildcard"],
        )
        self.assertEqual(len(matches), 1)

    def test_rejects_ambiguous_atom_name(self):
        atoms, dihedrals = MODULE.parse_psf_text(PSF.replace("A1", "C15"))
        with self.assertRaises(ValueError):
            MODULE.map_rotor(atoms, dihedrals, "C15", "N05")


if __name__ == "__main__":
    unittest.main()
