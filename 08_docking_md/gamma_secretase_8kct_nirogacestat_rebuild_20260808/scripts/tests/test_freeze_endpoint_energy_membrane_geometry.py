import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[2]
SCRIPT = ROOT / "scripts" / "freeze_endpoint_energy_membrane_geometry.py"


def rows(thickness, flip=0):
    return [
        {
            "time_ns": 200.0 + index,
            "phosphate_peak_thickness_nm": thickness,
            "upper_phosphate_peak_z_relative_nm": thickness / 2.0,
            "lower_phosphate_peak_z_relative_nm": -thickness / 2.0,
            "upper_leaflet_mismatch_count": 0,
            "lower_leaflet_mismatch_count": 0,
            "cumulative_leaflet_flip_events": flip,
        }
        for index in range(301)
    ]


class MembraneGeometryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = None
        if SCRIPT.is_file():
            spec = importlib.util.spec_from_file_location("membrane_geometry", SCRIPT)
            cls.module = importlib.util.module_from_spec(spec)
            assert spec.loader is not None
            spec.loader.exec_module(cls.module)

    def setUp(self):
        if self._testMethodName != "test_implementation_exists" and self.module is None:
            self.skipTest("implementation not created yet")

    def test_implementation_exists(self):
        self.assertTrue(SCRIPT.is_file(), f"missing implementation: {SCRIPT.name}")

    def test_freezes_across_replica_median_in_angstrom(self):
        result = self.module.freeze_geometry({
            "rep01": rows(4.0),
            "rep02": rows(4.2),
            "rep03": rows(4.4),
        })
        self.assertEqual(result["replica_median_thickness_nm"], {"rep01": 4.0, "rep02": 4.2, "rep03": 4.4})
        self.assertAlmostEqual(result["mthick_angstrom"], 42.0)
        self.assertEqual(result["mctrdz_angstrom"], 0.0)

    def test_rejects_missing_replica(self):
        with self.assertRaisesRegex(ValueError, "rep03"):
            self.module.freeze_geometry({"rep01": rows(4.0), "rep02": rows(4.2)})

    def test_rejects_leaflet_flip_confounded_input(self):
        with self.assertRaisesRegex(ValueError, "flip"):
            self.module.freeze_geometry({"rep01": rows(4.0), "rep02": rows(4.2, flip=1), "rep03": rows(4.4)})

    def test_rejects_inverted_leaflets(self):
        invalid = rows(4.2)
        invalid[0]["upper_phosphate_peak_z_relative_nm"] = -2.1
        with self.assertRaisesRegex(ValueError, "leaflet"):
            self.module.freeze_geometry({"rep01": rows(4.0), "rep02": invalid, "rep03": rows(4.4)})

    def test_requires_prepared_midplane_zero_and_z_normal(self):
        manifests = {
            rep: {"status": "pass", "membrane_normal_axis": "z", "membrane_midplane_z_angstrom": 0.0}
            for rep in ("rep01", "rep02", "rep03")
        }
        self.module.validate_preparation_manifests(manifests)
        manifests["rep03"]["membrane_normal_axis"] = "x"
        with self.assertRaisesRegex(ValueError, "rep03"):
            self.module.validate_preparation_manifests(manifests)


if __name__ == "__main__":
    unittest.main()
