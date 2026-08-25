import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[2]
SCRIPT = ROOT / "scripts" / "collect_endpoint_cpu_inventory.py"


class EndpointCpuInventoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = None
        if SCRIPT.is_file():
            spec = importlib.util.spec_from_file_location("cpu_inventory", SCRIPT)
            cls.module = importlib.util.module_from_spec(spec)
            assert spec.loader is not None
            spec.loader.exec_module(cls.module)

    def test_implementation_exists(self):
        self.assertTrue(SCRIPT.is_file())

    def test_meminfo_parser_converts_kib(self):
        if self.module is None:
            self.skipTest("implementation not created yet")
        self.assertEqual(self.module.parse_mem_total_bytes("MemTotal:       12345 kB\n"), 12345 * 1024)
        with self.assertRaisesRegex(ValueError, "MemTotal"):
            self.module.parse_mem_total_bytes("MemFree: 1 kB\n")


if __name__ == "__main__":
    unittest.main()
