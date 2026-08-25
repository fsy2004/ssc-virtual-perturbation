#!/usr/bin/env python3
"""Seal a reviewed release JSON with an immutable canonical payload checksum."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from md_contract import load_json, report_payload_sha256, sha256


TODO_RE = re.compile(r"\b(?:TODO|TBD|PENDING|PLACEHOLDER)\b", re.IGNORECASE)


def has_todo(value: object) -> bool:
    if isinstance(value, str):
        return bool(TODO_RE.search(value))
    if isinstance(value, dict):
        return any(has_todo(item) for item in value.values())
    if isinstance(value, list):
        return any(has_todo(item) for item in value)
    return False


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    source = args.input.resolve(); output = args.output.resolve()
    if output.exists():
        raise SystemExit(f"Refusing to overwrite sealed report: {output}")
    report = load_json(source)
    integrity = report.get("integrity")
    if not isinstance(integrity, dict) or "payload_sha256" not in integrity:
        raise SystemExit("Report must contain integrity.payload_sha256")
    report["integrity"]["payload_sha256"] = "UNSEALED"
    if has_todo(report):
        raise SystemExit("Refusing to seal a report containing TODO/TBD values")
    if report.get("status") != "pass" or report.get("approval_status") != "approved":
        raise SystemExit("Refusing to seal a release report that is not approved PASS")
    report["integrity"]["payload_sha256"] = report_payload_sha256(report, ("integrity", "payload_sha256"))
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")
    print(json.dumps({"status": "sealed", "output": str(output), "sha256": sha256(output)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
