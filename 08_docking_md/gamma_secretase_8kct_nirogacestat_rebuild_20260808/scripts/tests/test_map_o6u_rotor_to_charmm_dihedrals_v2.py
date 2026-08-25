import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "map_o6u_rotor_to_charmm_dihedrals_v2.py"
SPEC = importlib.util.spec_from_file_location("rotor_map_v2", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class CorrespondenceTests(unittest.TestCase):
    def test_maps_ccd_to_cgenff_names(self):
        text = (
            "ordinal\tccd_atom_id\tcgenff_atom_name\n"
            "1\tC15\tC3\n"
            "2\tN05\tN2\n"
        )
        mapping = MODULE.parse_correspondence_text(text)
        self.assertEqual(mapping["C15"], "C3")
        self.assertEqual(mapping["N05"], "N2")

    def test_rejects_duplicate_ccd_name(self):
        text = (
            "ordinal\tccd_atom_id\tcgenff_atom_name\n"
            "1\tC15\tC3\n"
            "2\tC15\tC4\n"
        )
        with self.assertRaises(ValueError):
            MODULE.parse_correspondence_text(text)


if __name__ == "__main__":
    unittest.main()
