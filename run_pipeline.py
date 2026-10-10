#!/usr/bin/env python3
"""
Simple command-line runner for the calculation pipeline.

Works whether this file lives in scripts/ or in the project root.
Looks for the oxide table at:
    data/master/xenolith_geochem.csv
    data/xenolith_geochem.csv
"""

import sys
from pathlib import Path


def find_root() -> Path:
    here = Path(__file__).resolve().parent
    candidates = [here, here.parent]
    for root in candidates:
        if (root / "calculations").is_dir() and (
            (root / "data" / "xenolith_geochem.csv").is_file()
            or (root / "data" / "master" / "xenolith_geochem.csv").is_file()
        ):
            return root
        if (root / "calculations").is_dir():
            return root
    return here


ROOT = find_root()
sys.path.insert(0, str(ROOT))

from calculations import run_full_pipeline


def find_input() -> Path:
    for path in (
        ROOT / "data" / "master" / "xenolith_geochem.csv",
        ROOT / "data" / "xenolith_geochem.csv",
    ):
        if path.is_file():
            return path
    raise FileNotFoundError(
        "Cannot find xenolith_geochem.csv under "
        f"{ROOT / 'data'} or {ROOT / 'data' / 'master'}"
    )


DATA_FILE = find_input()
DERIVED = ROOT / "data" / "derived"
OUTPUT_FILE = DERIVED / "p-T-fO2.csv"
FULL_FILE = DERIVED / "screened_PT.csv"

if __name__ == "__main__":
    print(f"Project root: {ROOT}")
    print(f"Input:        {DATA_FILE}")
    result, stats = run_full_pipeline(
        DATA_FILE,
        verbose=True,
        write_csv=OUTPUT_FILE,
        write_full_csv=FULL_FILE,
    )

    print("\n--- Summary ---")
    for k, v in stats.items():
        print(f"  {k}: {v}")
    print(f"\nWrote {OUTPUT_FILE}")
    print(f"Wrote {FULL_FILE}")
