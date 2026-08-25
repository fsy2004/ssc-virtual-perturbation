#!/usr/bin/env python3
"""Contract tests for the three-realization HMR production release builder."""

from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).with_name("build_three_replica_release.py")


def load_builder():
    spec = importlib.util.spec_from_file_location("build_three_replica_release", MODULE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load builder module")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ReleaseContractTests(unittest.TestCase):
    def test_step61_replaces_only_random_seeds(self) -> None:
        builder = load_builder()
        source = (
            "integrator = md\n"
            "dt = 0.001\n"
            "nsteps = 125000\n"
            "tcoupl = v-rescale\n"
            "ref_t = 303.15 303.15 303.15\n"
            "constraints = h-bonds\n"
            "gen-vel = yes\n"
            "gen-temp = 303.15\n"
            "gen-seed = -1\n"
        )
        rendered = builder.render_step61(source, velocity_seed=26081601, thermostat_seed=126081601)
        self.assertIn("gen-seed                = 26081601", rendered)
        self.assertIn("ld-seed                 = 126081601", rendered)
        self.assertNotIn("gen-seed = -1", rendered)
        self.assertEqual(rendered.count("gen-vel"), 1)
        self.assertIn("gen-vel = yes", rendered)

    def test_production_contract_rejects_velocity_regeneration(self) -> None:
        builder = load_builder()
        bad = """
integrator = md
dt = 0.004
nsteps = 125000000
tcoupl = v-rescale
ref_t = 303.15 303.15 303.15
pcoupl = C-rescale
pcoupltype = semiisotropic
constraints = h-bonds
constraint-algorithm = LINCS
continuation = yes
gen-vel = yes
"""
        with self.assertRaisesRegex(ValueError, "gen-vel must be no"):
            builder.validate_production_mdp(bad)

    def test_production_contract_accepts_frozen_500ns_scheme(self) -> None:
        builder = load_builder()
        good = """
integrator = md
dt = 0.004
nsteps = 125000000
tcoupl = v-rescale
tc-grps = SOLU MEMB SOLV
ref-t = 303.15 303.15 303.15
pcoupl = C-rescale
pcoupltype = semiisotropic
constraints = h-bonds
constraint-algorithm = lincs
continuation = yes
gen-vel = no
"""
        parsed = builder.validate_production_mdp(good)
        self.assertEqual(parsed["nsteps"], "125000000")
        self.assertEqual(parsed["dt"], "0.004")

    def test_manifest_has_exact_three_distinct_realizations(self) -> None:
        builder = load_builder()
        protocol = {
            "realizations": [
                {"id": "rep01", "velocity_seed": 26081601, "thermostat_seed": 126081601},
                {"id": "rep02", "velocity_seed": 26081602, "thermostat_seed": 126081602},
                {"id": "rep03", "velocity_seed": 26081603, "thermostat_seed": 126081603},
            ]
        }
        builder.validate_realizations(protocol)
        protocol["realizations"][2]["velocity_seed"] = 26081602
        with self.assertRaisesRegex(ValueError, "velocity seeds must be distinct"):
            builder.validate_realizations(protocol)

    def test_runner_does_not_treat_normal_infinite_dielectric_as_failure(self) -> None:
        builder = load_builder()
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary)
            builder._write_runner(destination)
            runner = (destination / "run_replica.sh").read_text(encoding="utf-8")
            self.assertNotIn("(NaN|Inf)", runner)
            self.assertIn("NaN([^A-Za-z]|$)", runner)
            self.assertIn('read -r -a mdrun_args', runner)
            self.assertIn('mdrun "${mdrun_args[@]}"', runner)

    def test_runner_appends_only_when_a_checkpoint_exists(self) -> None:
        builder = load_builder()
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary)
            builder._write_runner(destination)
            runner = (destination / "run_replica.sh").read_text(encoding="utf-8")
            self.assertIn('if [[ -s "$work/${prefix}.cpt" ]]', runner)
            self.assertIn('restart_args=(-cpi "$work/${prefix}.cpt" -append)', runner)
            self.assertIn('restart_args=()', runner)
            self.assertIn('"${restart_args[@]}"', runner)

    def test_write_release_is_fail_closed_on_existing_destination(self) -> None:
        builder = load_builder()
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary) / "release"
            destination.mkdir()
            with self.assertRaisesRegex(FileExistsError, "release destination already exists"):
                builder.ensure_new_destination(destination)


if __name__ == "__main__":
    unittest.main()
