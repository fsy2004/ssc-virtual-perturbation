import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "analyze_o6u_torsion_type_conflicts_v1.py"
SPEC = importlib.util.spec_from_file_location("torsion_conflicts", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class TorsionConflictTests(unittest.TestCase):
    def test_signature_is_orientation_invariant(self):
        forward = ["A", "B", "C", "D"]
        reverse = ["D", "C", "B", "A"]
        self.assertEqual(MODULE.canonical_signature(forward), MODULE.canonical_signature(reverse))

    def test_detects_nonlocal_central_bond_collision(self):
        occurrences = [
            {"central_bond": [2, 3]},
            {"central_bond": [5, 6]},
        ]
        self.assertTrue(MODULE.has_nonlocal_collision(occurrences, {frozenset((2, 3))}))
        self.assertFalse(MODULE.has_nonlocal_collision(occurrences[:1], {frozenset((2, 3))}))


if __name__ == "__main__":
    unittest.main()
