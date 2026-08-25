#!/usr/bin/env python3
"""Create or validate a complete immutable SHA-256 manifest for a directory tree."""

from __future__ import annotations

import argparse
import json
import tempfile
from datetime import datetime, timezone
from pathlib import Path

from md_contract import canonical_json_sha256, load_json, sha256


def inventory(root: Path, excluded: Path | None = None) -> list[dict[str, object]]:
    root = root.resolve()
    if not root.is_dir():
        raise ValueError(f"Tree root is not a directory: {root}")
    rows: list[dict[str, object]] = []
    for path in sorted((item for item in root.rglob("*") if item.is_file()), key=lambda item: item.as_posix()):
        if excluded is not None and path.resolve() == excluded.resolve():
            continue
        relative = path.relative_to(root).as_posix()
        if relative.startswith("../") or relative == "..":
            raise ValueError(f"Tree entry escapes root: {path}")
        rows.append({"path": relative, "bytes": path.stat().st_size, "sha256": sha256(path)})
    if not rows:
        raise ValueError("Refusing to manifest an empty directory tree")
    return rows


def validate_report(report: dict[str, object], root: Path, report_path: Path) -> list[str]:
    errors: list[str] = []
    if report.get("schema_version") != "1.0" or report.get("report_type") != "complete_directory_sha256_manifest":
        errors.append("tree-manifest schema/report type is invalid")
    try:
        observed = inventory(root, report_path if report_path.parent.resolve() == root.resolve() else None)
    except ValueError as exc:
        return [str(exc)]
    expected = report.get("files")
    if expected != observed:
        errors.append("directory contents, byte sizes, or file hashes differ from tree manifest")
    if report.get("file_count") != len(observed):
        errors.append("tree-manifest file_count is stale")
    if report.get("tree_payload_sha256") != canonical_json_sha256(observed):
        errors.append("tree-manifest payload SHA-256 is stale")
    return errors


def self_test() -> None:
    with tempfile.TemporaryDirectory(prefix="tree_manifest_") as temporary:
        root = Path(temporary) / "tree"
        root.mkdir()
        (root / "a.txt").write_text("a\n", encoding="utf-8")
        nested = root / "nested"
        nested.mkdir()
        (nested / "b.txt").write_text("b\n", encoding="utf-8")
        report_path = Path(temporary) / "tree.json"
        files = inventory(root)
        report = {
            "schema_version": "1.0", "report_type": "complete_directory_sha256_manifest",
            "root": str(root.resolve()), "file_count": len(files), "files": files,
            "tree_payload_sha256": canonical_json_sha256(files),
        }
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        if validate_report(load_json(report_path), root, report_path):
            raise RuntimeError("valid tree manifest failed")
        (nested / "b.txt").write_text("tampered\n", encoding="utf-8")
        if not validate_report(load_json(report_path), root, report_path):
            raise RuntimeError("tampered tree passed")
    print("SELF-TEST PASS: complete tree manifest detects changed content.")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path, help="Validate this existing manifest instead of writing")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    if args.root is None:
        parser.error("--root is required")
    root = args.root.resolve()
    if args.validate:
        report_path = args.validate.resolve()
        errors = validate_report(load_json(report_path), root, report_path)
        for error in errors:
            print(f"ERROR: {error}")
        return 1 if errors else 0
    if args.output is None:
        parser.error("--output is required when not using --validate")
    output = args.output.resolve()
    if output.exists():
        raise SystemExit(f"Refusing to overwrite tree manifest: {output}")
    files = inventory(root, output if output.parent == root else None)
    report = {
        "schema_version": "1.0",
        "report_type": "complete_directory_sha256_manifest",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "root": str(root),
        "file_count": len(files),
        "files": files,
        "tree_payload_sha256": canonical_json_sha256(files),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")
    print(json.dumps({"status": "pass", "output": str(output), "sha256": sha256(output)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
