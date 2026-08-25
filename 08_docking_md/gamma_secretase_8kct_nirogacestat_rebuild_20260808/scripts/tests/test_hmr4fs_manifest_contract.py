import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
BUILDER = ROOT / "scripts" / "build_study_manifest.py"
CONTRACT = ROOT / "scripts" / "md_contract.py"
PREFLIGHT = ROOT / "scripts" / "validate_preflight.py"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def contains_placeholder(value) -> bool:
    if isinstance(value, str):
        upper = value.upper()
        return any(token in upper for token in ("TODO", "TBD", "PENDING", "PLACEHOLDER"))
    if isinstance(value, dict):
        return any(contains_placeholder(item) for item in value.values())
    if isinstance(value, list):
        return any(contains_placeholder(item) for item in value)
    return False


class Hmr4fsManifestContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.builder = load_module("build_study_manifest_hmr4fs", BUILDER)
        cls.contract = load_module("md_contract_hmr4fs", CONTRACT)
        cls.preflight = load_module("validate_preflight_hmr4fs", PREFLIGHT)

    def test_builder_uses_frozen_hmr4fs_protocol_and_does_not_false_release(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.builder.OUT_MANIFEST = root / "study_manifest.json"
            self.builder.OUT_PLAN = root / "analysis_plan.json"
            self.assertEqual(self.builder.main(), 0)
            manifest = json.loads(self.builder.OUT_MANIFEST.read_text(encoding="utf-8"))

        self.assertEqual(manifest["global_model"]["temperature_k"], 303.15)
        self.assertEqual(manifest["simulation"]["time_step_ps"], 0.004)
        self.assertIs(manifest["simulation"]["hydrogen_mass_repartitioning"], True)
        mdp = manifest["simulation"]["production_mdp_contract"]
        self.assertEqual(mdp["thermostat"], "v-rescale")
        self.assertEqual(mdp["thermostat_groups"], ["SOLU", "MEMB", "SOLV"])
        self.assertEqual(mdp["barostat"], "C-rescale")
        self.assertEqual(mdp["com_removal_groups"], ["SOLU_MEMB", "SOLV"])
        self.assertEqual(mdp["output_cadence_steps"]["nstxout_compressed"], 5000)
        seeds = {
            item["id"]: (item["velocity_seed"], item["thermostat_seed"])
            for item in manifest["systems"][0]["realizations"]
        }
        self.assertEqual(
            seeds,
            {
                "rep01": (26081601, 126081601),
                "rep02": (26081602, 126081602),
                "rep03": (26081603, 126081603),
            },
        )
        if manifest["manifest_status"] == "ready_for_trajectory_processing":
            self.assertFalse(contains_placeholder(manifest))

    def test_validator_accepts_actual_frozen_hmr4fs_production_mdp(self):
        mdp_text = """\
integrator = md
dt = 0.004
nsteps = 125000000
continuation = yes
gen-vel = no
pbc = xyz
periodic-molecules = no
nstxout-compressed = 5000
compressed-x-precision = 1000
nstxout = 0
nstvout = 0
nstfout = 0
nstcalcenergy = 100
nstenergy = 5000
nstlog = 5000
cutoff-scheme = Verlet
nstlist = 20
rlist = 1.2
vdwtype = Cut-off
vdw-modifier = Force-switch
rvdw-switch = 1.0
rvdw = 1.2
DispCorr = no
coulombtype = PME
rcoulomb = 1.2
pme-order = 4
fourierspacing = 0.12
tcoupl = v-rescale
tc-grps = SOLU MEMB SOLV
tau-t = 1.0 1.0 1.0
ref-t = 303.15 303.15 303.15
pcoupl = C-rescale
pcoupltype = semiisotropic
tau-p = 5.0
compressibility = 4.5e-5 4.5e-5
ref-p = 1.0 1.0
constraints = h-bonds
constraint-algorithm = LINCS
nstcomm = 100
comm-mode = linear
comm-grps = SOLU_MEMB SOLV
"""
        manifest = {
            "global_model": {"temperature_k": 303.15, "pressure_bar": 1.0},
            "simulation": {
                "production_ns": 500.0,
                "time_step_ps": 0.004,
                "hydrogen_mass_repartitioning": True,
                "constraints": "h-bonds",
                "pressure_coupling": "semiisotropic",
                "production_mdp_contract": {
                    "integrator": "md",
                    "thermostat": "v-rescale",
                    "thermostat_groups": ["SOLU", "MEMB", "SOLV"],
                    "tau_t_ps": 1.0,
                    "barostat": "C-rescale",
                    "barostat_tau_p_ps": 5.0,
                    "compressibility_bar_inverse": [4.5e-5, 4.5e-5],
                    "cutoff_scheme": "verlet",
                    "neighbor_list_update_steps": 20,
                    "rlist_nm": 1.2,
                    "rcoulomb_nm": 1.2,
                    "vdw_type": "cut-off",
                    "vdw_modifier": "force-switch",
                    "rvdw_switch_nm": 1.0,
                    "rvdw_nm": 1.2,
                    "dispersion_correction": "no",
                    "pme_order": 4,
                    "fourier_spacing_nm": 0.12,
                    "constraint_algorithm": "lincs",
                    "com_removal_mode": "linear",
                    "com_removal_groups": ["SOLU_MEMB", "SOLV"],
                    "com_removal_interval_steps": 100,
                    "output_cadence_steps": {
                        "nstxout": 0,
                        "nstvout": 0,
                        "nstfout": 0,
                        "nstxout_compressed": 5000,
                        "nstcalcenergy": 100,
                        "nstenergy": 5000,
                        "nstlog": 5000,
                        "compressed_x_precision": 1000,
                    },
                },
            },
        }
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "production_500ns.mdp"
            path.write_text(mdp_text, encoding="utf-8")
            result = self.contract.validate_production_mdp(path, manifest)

        self.assertEqual(result["duration_ns"], 500.0)
        self.assertEqual(result["temperature_k"], 303.15)
        self.assertIs(result["hydrogen_mass_repartitioning"], True)

    def test_preflight_protocol_field_gate_accepts_hmr4fs_and_rejects_legacy_values(self):
        valid = {
            "global_model": {"temperature_k": 303.15},
            "simulation": {
                "time_step_ps": 0.004,
                "hydrogen_mass_repartitioning": True,
                "production_ns": 500.0,
            },
        }
        audit = self.preflight.Audit()
        self.preflight.validate_frozen_production_protocol_fields(valid, audit)
        self.assertEqual(audit.errors, [])

        legacy = json.loads(json.dumps(valid))
        legacy["global_model"]["temperature_k"] = 310.15
        legacy["simulation"]["time_step_ps"] = 0.002
        legacy["simulation"]["hydrogen_mass_repartitioning"] = False
        audit = self.preflight.Audit()
        self.preflight.validate_frozen_production_protocol_fields(legacy, audit)
        self.assertEqual(len(audit.errors), 3)


if __name__ == "__main__":
    unittest.main()
