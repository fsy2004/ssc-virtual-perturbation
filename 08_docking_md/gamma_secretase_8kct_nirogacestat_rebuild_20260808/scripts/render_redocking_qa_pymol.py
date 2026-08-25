#!/usr/bin/env python3
"""Render a deterministic redocking-QA overlay, never a binding figure."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import struct
import subprocess
import sys
import time
import zlib
from pathlib import Path


WIDTH = 3780
HEIGHT = 2835


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def pml_path(path: Path) -> str:
    return path.resolve().as_posix().replace("'", "\\'")


def validate_png(path: Path) -> None:
    data = path.read_bytes()
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError("not a PNG")
    length = struct.unpack(">I", data[8:12])[0]
    if data[12:16] != b"IHDR" or length != 13:
        raise ValueError("missing canonical IHDR")
    width, height = struct.unpack(">II", data[16:24])
    if (width, height) != (WIDTH, HEIGHT):
        raise ValueError(f"wrong PNG dimensions: {(width, height)}")
    expected_crc = struct.unpack(">I", data[29:33])[0]
    if zlib.crc32(data[12:29]) & 0xFFFFFFFF != expected_crc:
        raise ValueError("IHDR CRC mismatch")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pymol", type=Path, required=True)
    ap.add_argument("--receptor-pdb", type=Path, required=True)
    ap.add_argument("--native-sdf", type=Path, required=True)
    ap.add_argument("--run-root", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--staging-dir", type=Path, required=True,
                    help="dedicated short Windows path outside the project for PyMOL I/O")
    args = ap.parse_args()

    seed_files = {
        seed: args.run_root / f"seed{seed}" / "O6U_redocked.sdf"
        for seed in (11111, 22222, 33333)
    }
    inputs = [args.pymol, args.receptor_pdb, args.native_sdf, *seed_files.values()]
    for path in inputs:
        if not path.is_file():
            ap.error(f"missing input: {path}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.staging_dir.mkdir(parents=True, exist_ok=True)
    declared = {
        "redocking_qa.pml", "redocking_qa.pse", "redocking_qa_600dpi.png",
        "pymol_version.txt", "camera.txt", "pymol.stdout.log",
        "pymol.stderr.log", "redocking_qa_figure_manifest.json",
    }
    extras = [p.name for p in args.output_dir.iterdir() if p.name not in declared]
    if extras:
        raise SystemExit(f"NO-GO: undeclared files already exist in output directory: {extras}")

    stage_names = {
        "receptor.pdb", "native.sdf", "seed11111.sdf", "seed22222.sdf",
        "seed33333.sdf", "redocking_qa.pml", "redocking_qa.pse",
        "redocking_qa_600dpi.png", "pymol_version.txt", "camera.txt",
    }
    stage_extras = [p.name for p in args.staging_dir.iterdir() if p.name not in stage_names]
    if stage_extras:
        raise SystemExit(f"NO-GO: staging directory is not dedicated: {stage_extras}")
    for name in stage_names:
        target = args.staging_dir / name
        if target.is_file():
            target.unlink()
    staged_receptor = args.staging_dir / "receptor.pdb"
    staged_native = args.staging_dir / "native.sdf"
    shutil.copy2(args.receptor_pdb, staged_receptor)
    shutil.copy2(args.native_sdf, staged_native)
    staged_seeds = {}
    for seed, source in seed_files.items():
        target = args.staging_dir / f"seed{seed}.sdf"
        shutil.copy2(source, target)
        staged_seeds[seed] = target

    pml = args.output_dir / "redocking_qa.pml"
    pse = args.output_dir / "redocking_qa.pse"
    png = args.output_dir / "redocking_qa_600dpi.png"
    version = args.output_dir / "pymol_version.txt"
    camera = args.output_dir / "camera.txt"
    stage_pml = args.staging_dir / pml.name
    stage_pse = args.staging_dir / pse.name
    stage_png = args.staging_dir / png.name
    stage_version = args.staging_dir / version.name
    stage_camera = args.staging_dir / camera.name
    lines = [
        "reinitialize",
        "set retain_order, 1",
        "set quiet, 1",
        "bg_color white",
        "set orthoscopic, 1",
        "set depth_cue, 0",
        "set fog, 0",
        "set antialias, 2",
        "set ray_trace_mode, 0",
        "set ray_shadows, 0",
        "set ambient, 0.48",
        "set direct, 0.72",
        "set specular, 0.12",
        "set shininess, 10",
        "set_color protein_grey, [0.70,0.70,0.70]",
        "set_color native_magenta, [0.80,0.47,0.65]",
        "set_color seed_blue, [0.00,0.45,0.70]",
        "set_color seed_green, [0.00,0.62,0.45]",
        "set_color seed_orange, [0.90,0.62,0.00]",
        "set_color atom_oxygen, [0.84,0.15,0.16]",
        "set_color atom_nitrogen, [0.12,0.33,0.78]",
        "set_color atom_fluorine, [0.45,0.78,0.46]",
        f"load {pml_path(staged_receptor)}, receptor",
        f"load {pml_path(staged_native)}, native_O6U",
    ]
    for seed, path in staged_seeds.items():
        lines.extend([
            f"load {pml_path(path)}, seed{seed}_all",
            f"create seed{seed}_top1, seed{seed}_all, 1, 1",
            f"delete seed{seed}_all",
        ])
    lines.extend([
        "hide everything, all",
        "select literature_pocket, receptor and chain B and resi 77+261+268+271+272+282+287+379+380+381+425+431+432",
        "show sticks, literature_pocket",
        "color protein_grey, literature_pocket and elem C",
        "show sticks, native_O6U or seed11111_top1 or seed22222_top1 or seed33333_top1",
        "set stick_radius, 0.18, literature_pocket",
        "set stick_radius, 0.26, native_O6U or seed11111_top1 or seed22222_top1 or seed33333_top1",
        "color native_magenta, native_O6U and elem C",
        "color seed_blue, seed11111_top1 and elem C",
        "color seed_green, seed22222_top1 and elem C",
        "color seed_orange, seed33333_top1 and elem C",
        "color atom_oxygen, (native_O6U or seed11111_top1 or seed22222_top1 or seed33333_top1) and elem O",
        "color atom_nitrogen, (native_O6U or seed11111_top1 or seed22222_top1 or seed33333_top1) and elem N",
        "color atom_fluorine, (native_O6U or seed11111_top1 or seed22222_top1 or seed33333_top1) and elem F",
        "orient native_O6U",
        "zoom native_O6U or seed11111_top1 or seed22222_top1 or seed33333_top1, 7.0",
        "turn x, -12",
        "turn y, 18",
        "python",
        "import json",
        "from pymol import cmd",
        "names = ['receptor','native_O6U','seed11111_top1','seed22222_top1','seed33333_top1']",
        "counts = {name: {'all': cmd.count_atoms(name), 'heavy': cmd.count_atoms(f'({name}) and not hydro')} for name in names}",
        "assert counts['receptor']['all'] > 10000 and all(counts[name]['heavy'] == 35 for name in names if name != 'receptor'), counts",
        f"open(r'{pml_path(stage_version)}','w',encoding='utf-8').write(repr(cmd.get_version())+'\\n')",
        f"open(r'{pml_path(stage_camera)}','w',encoding='utf-8').write(json.dumps({{'view': list(cmd.get_view()), 'atom_counts': counts}}, indent=2)+'\\n')",
        f"cmd.save(r'{pml_path(stage_pse)}')",
        f"cmd.png(r'{pml_path(stage_png)}', width={WIDTH}, height={HEIGHT}, dpi=600, ray=1, quiet=0)",
        "cmd.sync()",
        "python end",
        "quit",
    ])
    pml.write_text("\n".join(lines) + "\n", encoding="utf-8")
    stage_pml.write_text(pml.read_text(encoding="utf-8"), encoding="utf-8")
    proc = subprocess.run([str(args.pymol), "-cq", str(stage_pml)], capture_output=True, text=True)
    stdout_log = args.output_dir / "pymol.stdout.log"
    stderr_log = args.output_dir / "pymol.stderr.log"
    stdout_log.write_text(proc.stdout, encoding="utf-8")
    stderr_log.write_text(proc.stderr, encoding="utf-8")
    if proc.returncode != 0:
        raise SystemExit(f"NO-GO: PyMOL exited {proc.returncode}")
    # The Windows GUI launcher can return before its rendering child finishes.
    expected_stage = (stage_pse, stage_png, stage_version, stage_camera)
    deadline = time.monotonic() + 300.0
    prior_sizes = None
    stable_polls = 0
    while time.monotonic() < deadline:
        if all(p.is_file() and p.stat().st_size > 0 for p in expected_stage):
            sizes = tuple(p.stat().st_size for p in expected_stage)
            if sizes == prior_sizes:
                stable_polls += 1
                if stable_polls >= 4:
                    break
            else:
                stable_polls = 0
                prior_sizes = sizes
        time.sleep(0.25)
    else:
        raise SystemExit("NO-GO: timed out waiting for PyMOL outputs")
    for source, target in ((stage_pse, pse), (stage_png, png),
                           (stage_version, version), (stage_camera, camera)):
        if source.is_file() and source.stat().st_size:
            shutil.copy2(source, target)
    for path in (pse, png, version, camera):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"NO-GO: missing/empty PyMOL output {path}")
    validate_png(png)
    manifest = {
        "schema_version": 1,
        "figure_type": "failed-self-redocking QA overlay",
        "scientific_status": "not manuscript evidence",
        "legend": {
            "native_O6U": "#CC79A7",
            "seed11111_top_ranked": "#0072B2",
            "seed22222_top_ranked": "#009E73",
            "seed33333_top_ranked": "#E69F00",
            "protein": "#B3B3B3",
        },
        "interaction_segments": 0,
        "reason_no_interaction_segments": "PLIP endpoint geometry is reserved for the native-structure figure; docking QA failed its all-seed top-rank criterion",
        "inputs": {str(p.resolve()): sha256(p) for p in inputs},
        "outputs": {
            str(p.resolve()): sha256(p)
            for p in (pml, pse, png, version, camera, stdout_log, stderr_log)
        },
        "png": {"width_px": WIDTH, "height_px": HEIGHT, "requested_dpi": 600},
    }
    out = args.output_dir / "redocking_qa_figure_manifest.json"
    out.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
