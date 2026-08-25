#!/usr/bin/env python3
"""Hard prohibition for PCA/FEL/3D population maps in the 3x500 ns study."""

from __future__ import annotations

import argparse
import sys


PROHIBITION = (
    "PCA, occupancy-derived -kBT ln(P) maps, FEL labels, and 3D population "
    "surfaces are prohibited for this fixed 3x500 ns / 200-500 ns protocol. "
    "This executable intentionally has no analysis path."
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true", help="Confirm that no PCA/FEL execution path exists")
    args = parser.parse_args()
    if args.self_test:
        print(f"SELF-TEST PASS: {PROHIBITION}")
        return 0
    print(f"PROHIBITED: {PROHIBITION}", file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
