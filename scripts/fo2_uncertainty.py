#!/usr/bin/env python3
"""
Residual-based fO2 uncertainty for comparative-geothermometry estimates.

Usage (from the project root, after ``python run_pipeline.py``):

    python scripts/fo2_uncertainty.py

For the 57 calibration + validation samples (30 screened Slave and 27
screened Siberian + Kaapvaal samples with non-zero measured garnet
Fe3+/ΣFe), ΔFMQ is calculated (Stagno et al. 2013, same P, T and mineral
chemistry as the pipeline) from the estimated and from the measured
garnet Fe3+/ΣFe.  The fO2 1σ is the standard deviation of
ΔFMQ(estimated) − ΔFMQ(measured).  Error groups are retained (the
conservative case).  The samples are split at their median estimated
ΔFMQ into reduced and oxidised halves; these two 1σ values are the
Fig. 10 bars, and the per-dataset values are the Table 1 log column.

Outputs:
  data/derived/fO2_residuals.csv      per-sample values
  data/derived/fO2_uncertainty.json   overall, reduced and oxidised 1σ
  Figures/diagnostics/fO2_residual_vs_Fe3.png
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from calculations.pt_fo2 import calibration as cal  # noqa: E402

DERIVED = ROOT / "data" / "derived"
DIAG = ROOT / "Figures" / "diagnostics"


def _plot(res, unc, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    g = res["error_group"].astype(bool)
    slave = res["dataset"] == "Slave calibration"
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.9))

    def sc(a, x, y):
        a.scatter(x[~g & slave], y[~g & slave], c="k", marker="*", s=60, label="Slave calibration")
        a.scatter(x[~g & ~slave], y[~g & ~slave], facecolors="none", edgecolors="k",
                  marker="*", s=60, label="Validation")
        a.scatter(x[g], y[g], c="#c0392b", s=30, label="Error groups (Fig. 5)")

    sc(ax[0], res["XFe3_meas"], res["res_log"].abs())
    ax[0].set_xlabel(r"Measured garnet Fe$^{3+}$/$\Sigma$Fe")
    ax[0].set_ylabel(r"|$\Delta$FMQ(est) − $\Delta$FMQ(meas)| (log units)")
    ax[0].legend(fontsize=7)
    sc(ax[1], res["XFe3_meas"], res["sens_per_002"])
    ax[1].set_xlabel(r"Measured garnet Fe$^{3+}$/$\Sigma$Fe")
    ax[1].set_ylabel(r"log units per 0.02 Fe$^{3+}$/$\Sigma$Fe")
    sc(ax[2], res["dFMQ_est"], res["res_log"])
    ax[2].axhline(0, color="0.5", lw=0.6)
    ax[2].axvline(unc["split"], color="0.4", ls="--")
    ax[2].set_xlabel(r"Estimated $\Delta$FMQ")
    ax[2].set_ylabel("Residual (log units)")
    ax[2].text(0.03, 0.97,
               f"split {unc['split']:.2f}\n"
               f"reduced 1σ {unc['reduced']['sigma_log']:.2f} (n={unc['reduced']['n']})\n"
               f"oxidised 1σ {unc['oxidised']['sigma_log']:.2f} (n={unc['oxidised']['n']})\n"
               f"all 1σ {unc['all']['sigma_log']:.2f} (n={unc['all']['n']})",
               transform=ax[2].transAxes, va="top", fontsize=7)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main():
    df = pd.read_csv(DERIVED / "screened_PT.csv", encoding="utf-8-sig", low_memory=False)
    df = df.rename(columns={df.columns[0]: "SampleID"})
    res = cal.fo2_residuals(df)
    unc = cal.fo2_uncertainty(res)
    no_g = cal.fo2_uncertainty(res[~res["error_group"].astype(bool)], split=unc["split"])
    res.round(4).to_csv(DERIVED / "fO2_residuals.csv", index=False)
    out = dict(unc, error_groups_removed=dict(all=no_g["all"], reduced=no_g["reduced"],
                                              oxidised=no_g["oxidised"]))
    rnd = lambda d: {k: (round(v, 4) if isinstance(v, float) else rnd(v) if isinstance(v, dict) else v)
                     for k, v in d.items()}
    (DERIVED / "fO2_uncertainty.json").write_text(json.dumps(rnd(out), indent=2) + "\n")
    _plot(res, unc, DIAG / "fO2_residual_vs_Fe3.png")
    print(f"n = {unc['all']['n']}: 1σ {unc['all']['sigma_log']:.2f}, "
          f"RMS {unc['all']['rms_log']:.2f}, bias {unc['all']['bias_log']:+.2f} log units")
    print(f"split at median estimated ΔFMQ {unc['split']:.2f}")
    for k in ("reduced", "oxidised"):
        print(f"  {k:8s} n = {unc[k]['n']:2d}  1σ {unc[k]['sigma_log']:.2f}  "
              f"bias {unc[k]['bias_log']:+.2f}")
    print(f"error groups removed: n = {no_g['all']['n']}, 1σ {no_g['all']['sigma_log']:.2f}")
    print(f"median sensitivity at measured Fe3+/ΣFe: {np.median(res['sens_per_002']):.2f} "
          "log units per 0.02")


if __name__ == "__main__":
    main()
