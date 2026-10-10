#!/usr/bin/env python3
"""
Generate the Fig. 12 reference buffers (EMOG/D, Fe–Ni precipitation,
BSE trend, Nimis 2015 mantle band) and overlay them on P–ΔFMQ.

Usage (from the project root):
    python scripts/make_fig12_buffers.py

Writes
    data/derived/fig12_buffers.csv
    Figures/Fig12_P_fO2.png / .pdf
and prints a literature-anchor check.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent if (HERE.parent / "calculations").is_dir() else HERE
sys.path.insert(0, str(ROOT))

from calculations.pt_fo2.buffers import (
    buffer_table,
    nimis2015_band,
    verify_against_literature,
)

FULL = ROOT / "data" / "derived" / "screened_PT.csv"
OUT_CSV = ROOT / "data" / "derived" / "fig12_buffers.csv"
OUT_FIG = ROOT / "Figures"


CRATON_STYLE = {
    "Slave Craton": dict(marker="o", color="#c44e52", label="Slave"),
    "Siberian Craton": dict(marker="s", color="#4c72b0", label="Siberian"),
    "Kaapvaal Craton": dict(marker="^", color="#55a868", label="Kaapvaal"),
}


def main():
    print("=== Fig. 12 buffer check ===")
    verify_against_literature(verbose=True)

    df = None
    if FULL.exists():
        df = pd.read_csv(FULL, encoding="utf-8-sig")
        print(f"Loaded {len(df)} rows from {FULL.name}")
    else:
        print(f"No {FULL} — drawing buffers only.")

    P = np.linspace(20, 80, 121)
    buf = buffer_table(P, df=df)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    buf.to_csv(OUT_CSV, index=False)
    print(f"Wrote {OUT_CSV}")

    fig, ax = plt.subplots(figsize=(4.8, 5.6))
    lo, hi = nimis2015_band(P)
    ax.fill_betweenx(P, lo, hi, color="0.85", lw=0, zorder=0,
                     label="typical mantle (Nimis 2015)")
    ax.plot(buf["dEMOD"], P, color="#d4a017", lw=1.4, ls="--",
            label="EMOG/D")
    ax.plot(buf["dFeNi"], P, color="#4c72b0", lw=1.4, ls="-.",
            label="Fe–Ni precipitation")
    ax.plot(buf["dBSE"], P, color="k", lw=1.2, label="BSE")
    ax.axvline(0, color="0.4", lw=0.6)

    if df is not None and "delta_logfO2_FMQ" in df.columns:
        use = df
        if "screen_pass" in df.columns:
            use = df[df["screen_pass"] == True]
        for craton, sty in CRATON_STYLE.items():
            sub = use[use["Craton"].fillna("").astype(str) == craton]
            ax.scatter(sub["delta_logfO2_FMQ"], sub["P_CA91"],
                       s=28, edgecolors="k", linewidths=0.35,
                       marker=sty["marker"], c=sty["color"],
                       label=sty["label"], zorder=3)
        other = use[~use["Craton"].fillna("").astype(str).isin(CRATON_STYLE)]
        if not other.empty:
            ax.scatter(other["delta_logfO2_FMQ"], other["P_CA91"],
                       s=16, edgecolors="0.5", facecolors="none",
                       linewidths=0.5, marker="o", label="Other", zorder=2)

    ax.set_xlabel(r"log $f$O$_2$ ($\Delta$FMQ)")
    ax.set_ylabel("P (kbar)")
    ax.set_xlim(-6, 1)
    ax.set_ylim(80, 20)
    ax.tick_params(direction="in", top=True, right=True)
    ax.legend(frameon=False, fontsize=7, loc="lower left")
    ax.set_title("Fig. 12  P–fO₂ with buffers", loc="left", fontsize=10)
    OUT_FIG.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(OUT_FIG / "Fig12_P_fO2.png", dpi=300)
    fig.savefig(OUT_FIG / "Fig12_P_fO2.pdf")
    plt.close(fig)
    print(f"Wrote {OUT_FIG / 'Fig12_P_fO2.png'}")


if __name__ == "__main__":
    main()
