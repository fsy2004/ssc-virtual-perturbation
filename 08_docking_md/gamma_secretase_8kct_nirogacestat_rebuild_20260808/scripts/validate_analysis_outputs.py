#!/usr/bin/env python3
"""Hard prohibition validator for absent PCA/FEL/3D-map outputs."""

from __future__ import annotations

import argparse
import sys


PROHIBITION = (
    "PCA/FEL/3D-map outputs are outside the frozen 3x500 ns protocol and "
    "cannot be validated or used for a publication claim."
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true", help="Confirm that no PCA/FEL validation path exists")
    args = parser.parse_args()
    if args.self_test:
        print(f"SELF-TEST PASS: {PROHIBITION}")
        return 0
    print(f"PROHIBITED: {PROHIBITION}", file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
