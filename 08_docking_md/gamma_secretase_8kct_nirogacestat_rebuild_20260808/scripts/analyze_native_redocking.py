#!/usr/bin/env python3
"""Export and symmetry-correct the 8KCT/O6U self-redocking results."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolAlign


SCORE_RE = re.compile(r"REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-root", type=Path, required=True)
    ap.add_argument("--native-sdf", type=Path, required=True)
    ap.add_argument("--mk-export", type=Path, required=True)
    args = ap.parse_args()

    manifest_path = args.run_root / "redocking_run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    native = Chem.SDMolSupplier(str(args.native_sdf), removeHs=True)[0]
    if native is None or native.GetNumHeavyAtoms() != 35:
        raise SystemExit("NO-GO: native O6U SDF is unreadable or not 35 heavy atoms")

    rows = []
    seed_summary = []
    for run in manifest["runs"]:
        seed = int(run["seed"])
        pdbqt = Path(run["output"])
        if sha256(pdbqt) != run["output_sha256"]:
            raise SystemExit(f"NO-GO: seed {seed} output hash mismatch")
        sdf = pdbqt.with_suffix(".sdf")
        proc = subprocess.run(
            [str(args.mk_export), str(pdbqt), "-s", str(sdf)],
            capture_output=True,
            text=True,
        )
        (pdbqt.parent / "mk_export.stdout.log").write_text(proc.stdout, encoding="utf-8")
        (pdbqt.parent / "mk_export.stderr.log").write_text(proc.stderr, encoding="utf-8")
        if proc.returncode != 0 or not sdf.is_file():
            raise SystemExit(f"NO-GO: Meeko export failed for seed {seed}")

        scores = [float(x) for x in SCORE_RE.findall(pdbqt.read_text(encoding="utf-8"))]
        poses = [m for m in Chem.SDMolSupplier(str(sdf), removeHs=True) if m is not None]
        if len(scores) != len(poses) or not poses:
            raise SystemExit(
                f"NO-GO: seed {seed} has {len(scores)} scores but {len(poses)} readable poses"
            )
        for rank, (score, pose) in enumerate(zip(scores, poses), start=1):
            if pose.GetNumHeavyAtoms() != native.GetNumHeavyAtoms():
                raise SystemExit(f"NO-GO: seed {seed} rank {rank} atom-count mismatch")
            try:
                rmsd = float(rdMolAlign.GetBestRMS(native, pose, maxMatches=1_000_000))
            except RuntimeError as exc:
                raise SystemExit(f"NO-GO: seed {seed} rank {rank} graph mismatch: {exc}")
            rows.append(
                {
                    "seed": seed,
                    "rank": rank,
                    "vina_score_kcal_mol": score,
                    "symmetry_corrected_heavy_atom_rmsd_A": round(rmsd, 6),
                    "within_preregistered_2A_threshold": rmsd <= 2.0,
                }
            )
        this = [r for r in rows if r["seed"] == seed]
        top = this[0]
        best = min(this, key=lambda r: r["symmetry_corrected_heavy_atom_rmsd_A"])
        seed_summary.append(
            {
                "seed": seed,
                "top_ranked_rmsd_A": top["symmetry_corrected_heavy_atom_rmsd_A"],
                "top_ranked_score_kcal_mol": top["vina_score_kcal_mol"],
                "scoring_success_top_ranked_le_2A": top["within_preregistered_2A_threshold"],
                "best_sampled_rank": best["rank"],
                "best_sampled_rmsd_A": best["symmetry_corrected_heavy_atom_rmsd_A"],
                "sampling_success_any_pose_le_2A": best["within_preregistered_2A_threshold"],
            }
        )

    csv_path = args.run_root / "redocking_all_poses.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    report = {
        "schema_version": 1,
        "rmsd_definition": "RDKit GetBestRMS symmetry-corrected heavy-atom RMSD",
        "preregistered_success_threshold_A": 2.0,
        "native_sdf": str(args.native_sdf.resolve()),
        "native_sdf_sha256": sha256(args.native_sdf),
        "run_manifest_sha256": sha256(manifest_path),
        "seed_results": seed_summary,
        "all_seeds_scoring_success": all(x["scoring_success_top_ranked_le_2A"] for x in seed_summary),
        "all_seeds_sampling_success": all(x["sampling_success_any_pose_le_2A"] for x in seed_summary),
        "claim_boundary": "protocol QA only; no affinity, stability, potency, or mechanism inference",
        "all_pose_csv": str(csv_path.resolve()),
        "all_pose_csv_sha256": sha256(csv_path),
    }
    report_path = args.run_root / "redocking_analysis_report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
