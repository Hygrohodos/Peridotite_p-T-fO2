#!/usr/bin/env python3
"""Fig. 4 (revision v1): weighted single-regression figure with background
circles whose AREA is proportional to each run's normalised weight, a size
key, n_eff, and a side table of the weights.

    python3 scripts/make_fig04_weightcircles.py                       # current data
    python3 scripts/make_fig04_weightcircles.py --input runs.csv --outdir out/

Method: weights and fit come only from scripts/weighted_thermometer_regression.py
(w = 1/(sigma_a^2 + excess^2), WLS centred at 1350 C); the figure adds nothing.

Runs from a clean checkout of the repository. Needs: numpy, pandas,
matplotlib, and the repo itself (scripts/make_figures.py for the journal
style, which imports calculations/pt_fo2). Optional cross-checks, each
skipped with a notice when it cannot apply (--no-crosscheck skips both):
  * pipeline: rebuilds the 13 runs from the raw experimental data through
    calculations.pt_fo2.exp_calibration and compares weights, coefficients
    and the SAME-convention covariance (both scaled by max(1, chi2_red)).
    Only when the input's SampleIDs equal the pipeline's usable runs.
  * item 1: data/reference/item1_experiment_weights.csv (in-repo copy of
    ../item1_weighted_calibration/results/experiment_weights.csv; also supplies the sensitivity-only quality_factor
    column of the weights CSV, else NaN).

Output (--outdir; default $REV1_OUT, else ../figures_revision_v1 next to the
repo if it exists, else Figures/): Fig04_TA98_WZ07_measured_Fe3.png/.pdf and
Fig04_weights.csv.

Colour (11 Nov, t1615u): palette taken from the paper's other figures
(scripts/make_figures.py, plot_exp_ta98_wz07.py), no new colours:
  Matjuschkin 2014  #5e3c99 purple  (accent of Figs 2 and 6)
  Rzehak 2020       #4c72b0 blue    (plot_exp_ta98_wz07.STYLE; Fig. 3a)
  other references  FALLBACK_COLOURS, also from make_figures.py, in order
  fit / 95 % band   black / grey 0.85 (identical to Fig. X b, same fit)
Fonts/line weights: make_figures._apply_journal_style (Arial 8 pt).
"""
from __future__ import annotations
import argparse
import os
import re
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(REPO))
import make_figures as mf
import weighted_thermometer_regression as wtr

EXP_STYLE = {  # colours from make_figures.py / plot_exp_ta98_wz07.py
    "Matjuschkin et al. 2014": dict(marker="o", color="#5e3c99", label="Matjuschkin 2014"),
    "Rzehak et al. 2020": dict(marker="s", color="#4c72b0", label="Rzehak 2020"),
}
# for references not in EXP_STYLE (all used elsewhere in make_figures.py)
FALLBACK_COLOURS = ("#b35806", "#55a868", "#c44e52", "#d4a017", "#8c5a2b", "#6b4c9a", "#8c8c8c")
FALLBACK_MARKERS = ("^", "D", "v", "P", "X", "h", "<")

SIZE_PER_PCT = 25.0   # matplotlib marker size s (pt^2, = diameter^2) per 1 % of total weight;
                      # drawn circle area = (pi/4) * s, so area is proportional to weight
MARK_S = 16           # fixed marker size s (pt^2)
LABEL_OFF = {  # manual label offsets (points) to avoid overlaps; others default to (5, 4)
    "F1100o": (5, -9), "F1100r": (6, 4), "F1200o": (7, 5), "F1200r": (7, -3),
    "F1300o": (9, -3), "F1300r": (-31, 4), "F14002o": (7, -8), "F14002r": (6, 4),
    "MaLm6": (6, 4), "MaLm11": (-40, 9), "MaLm1": (9, 7), "MaLm9": (-30, -3),
    "MaLm12": (7, -3),
}
LEADER = ("MaLm11", "MaLm1")   # close pair: short leader lines
XLIM, YLIM = (1120, 1660), (-330, 230)   # expanded automatically if data fall outside
ITEM1 = REPO / "data" / "reference" / "item1_experiment_weights.csv"   # in-repo copy of ../item1_weighted_calibration/results/experiment_weights.csv (rev. v3)
STEM = "Fig04_TA98_WZ07_measured_Fe3"


def default_outdir():
    if os.environ.get("REV1_OUT"):
        return Path(os.environ["REV1_OUT"])
    rev = REPO.parent / "figures_revision_v1"
    return rev if rev.is_dir() else REPO / "Figures"


def reference_styles(refs):
    """Style per reference: EXP_STYLE first, then the fallback cycle in order of appearance."""
    styles, k = {}, 0
    for r in list(EXP_STYLE) + [r for r in dict.fromkeys(refs) if r not in EXP_STYLE]:
        if r not in set(refs):
            continue
        if r in EXP_STYLE:
            styles[r] = EXP_STYLE[r]
        else:
            lab = re.sub(r"\s+et al\.?", "", r).strip() or "Unspecified source"
            styles[r] = dict(marker=FALLBACK_MARKERS[k % len(FALLBACK_MARKERS)],
                             color=FALLBACK_COLOURS[k % len(FALLBACK_COLOURS)], label=lab)
            k += 1
    return styles


def crosscheck_pipeline(runs, wt, fr):
    """Compare with calculations.pt_fo2.exp_calibration on the raw-data rebuild."""
    import plot_exp_ta98_wz07 as pe
    from calculations.pt_fo2 import exp_calibration as xc
    ex = pe.temperatures(pe.assemble())
    ex = ex[ex["usable"]].copy()
    if set(ex["SampleID"].astype(str)) != set(runs["SampleID"]) or not np.isclose(fr.t_ref, xc.T_REF) \
            or fr.thresholds != xc.SCREEN_THR or "sigma_analytical_C" in runs:
        print("notice: pipeline cross-check skipped (input runs/options differ from the pipeline's 13-run set)")
        return
    ex = ex.set_index(ex["SampleID"].astype(str)).loc[runs["SampleID"]].reset_index(drop=True)
    w_x = xc.weights(ex, "invvar", fr.sigma0)
    f = xc.wls(ex["T_TA98"].to_numpy(float), ex["dT_Fe2"].to_numpy(float), w_x, 1)
    scale = max(1.0, f["chi2_red"])   # pipeline returns the UNSCALED cov; band_cov applies this
    assert np.allclose(w_x, wt["w_raw_per_C2"], rtol=1e-12)
    assert np.allclose(f["beta"], [fr.b0, fr.b1], rtol=1e-10)
    assert np.allclose(f["cov"] * scale, np.asarray(fr.cov), rtol=1e-10)   # scaled vs scaled
    print("cross-check vs exp_calibration: OK")


def crosscheck_item1(wt):
    if not ITEM1.is_file():
        print(f"notice: item-1 cross-check skipped ({ITEM1} not found); quality_factor column = NaN")
        return np.full(len(wt), np.nan)
    ref = pd.read_csv(ITEM1).set_index("SampleID")
    if not set(wt["SampleID"]) <= set(ref.index.astype(str)):
        print("notice: item-1 cross-check skipped (input has runs not in item-1 results); quality_factor = NaN")
        return np.full(len(wt), np.nan)
    ref = ref.loc[wt["SampleID"]]
    assert np.allclose(ref["w_rel_primary"].to_numpy(), wt["weight_frac"], atol=6e-5)
    print("cross-check vs item-1 weights: OK")
    return ref["quality_factor"].to_numpy()


def _lims(v, default, pad):
    return (min(default[0], np.min(v) - pad), max(default[1], np.max(v) + pad))


def make_figure(runs, wt, fr, outdir):
    mf._apply_journal_style()
    e = runs.assign(dT_Fe2=runs["y_TA98_minus_WZ07_C"])
    pct = wt["weight_pct"].to_numpy(float)
    excess = wt["screen_excess_C"].to_numpy(float)
    sigma = wt["sigma_i_C"].to_numpy(float)
    styles = reference_styles(e["Reference"].tolist())

    fig = plt.figure(figsize=(174 / 25.4, 95 / 25.4))
    ax = fig.add_axes([0.085, 0.115, 0.575, 0.80])
    tx = fig.add_axes([0.675, 0.115, 0.32, 0.80]); tx.axis("off")

    ax.axhline(0, color="0.6", lw=0.7, zorder=0)
    # fit drawn over the data range, rounded to 0.01 C as params.wz07_offset_tmin/tmax_c
    tt = np.linspace(round(e["T_TA98"].min(), 2), round(e["T_TA98"].max(), 2), 100)
    lo, hi = fr.band(tt)
    ax.fill_between(tt, lo, hi, color="0.85", lw=0, zorder=1, label="95 % confidence band")
    ax.plot(tt, fr.predict(tt), color="k", lw=1.4, zorder=4, label=f"Weighted fit ({fr.n} experiments)")
    for ref_, sty in styles.items():
        s = (e["Reference"] == ref_).to_numpy()
        col = sty["color"]
        ax.scatter(e.loc[s, "T_TA98"], e.loc[s, "dT_Fe2"], s=SIZE_PER_PCT * pct[s],
                   marker="o", facecolors=col, alpha=0.20, edgecolors="none", zorder=2)
        ax.scatter(e.loc[s, "T_TA98"], e.loc[s, "dT_Fe2"], s=SIZE_PER_PCT * pct[s],
                   marker="o", facecolors="none", edgecolors=col, linewidths=0.7, alpha=0.85, zorder=2)
        ax.scatter(e.loc[s, "T_TA98"], e.loc[s, "dT_Fe2"], s=MARK_S, marker=sty["marker"],
                   color=col, edgecolors="k", linewidths=0.4, zorder=5, label=sty["label"])
    for i, r in e.iterrows():
        lead = r["SampleID"] in LEADER
        ax.annotate(r["SampleID"], (r["T_TA98"], r["dT_Fe2"]), fontsize=6,
                    xytext=LABEL_OFF.get(r["SampleID"], (5, 4)), textcoords="offset points", zorder=6,
                    arrowprops=dict(arrowstyle="-", lw=0.5, color="0.25", shrinkA=0.5, shrinkB=2.5) if lead else None)
    ax.set_title(rf"$T_{{\mathrm{{TA98}}}}-T_{{\mathrm{{WZ07}}}}$ = {fr.b0:.0f} ($\pm${fr.se_b0:.0f}) "
                 rf"{fr.b1:+.3f} ($\pm${fr.se_b1:.3f})($T_{{\mathrm{{TA98}}}}$ − {fr.t_ref:.0f})", fontsize=7.5)
    ax.set_xlim(*_lims(e["T_TA98"], XLIM, 40)); ax.set_ylim(*_lims(e["dT_Fe2"], YLIM, 50))
    mf._style_ax(ax, "TA98 temperature (°C)", r"$T_{\mathrm{TA98}}-T_{\mathrm{WZ07}}$ (°C), Fe$^{2+}$ only")
    h, l = ax.get_legend_handles_labels()
    leg1 = ax.legend(h, l, frameon=True, facecolor="white", edgecolor="0.55", framealpha=0.92,
                     fontsize=6.5, loc="upper right")
    ax.add_artist(leg1)
    # size key (axes coordinates): circle area proportional to weight
    keys = (5, 10, 20)
    ax.add_patch(FancyBboxPatch((0.02, 0.025), 0.42, 0.20, boxstyle="round,pad=0.005,rounding_size=0.01",
                                transform=ax.transAxes, fc="white", ec="0.55", lw=0.6, alpha=0.92, zorder=7))
    ax.text(0.23, 0.205, "Weight (circle area, % of total)", transform=ax.transAxes, fontsize=6.3,
            ha="center", va="top", zorder=8)
    AXW_PT = 0.575 * 174 / 25.4 * 72  # axes width in points
    for k, xk in zip(keys, (0.065, 0.19, 0.335)):
        ax.scatter([xk], [0.10], s=SIZE_PER_PCT * k, transform=ax.transAxes, facecolors="0.5",
                   alpha=0.25, edgecolors="none", zorder=8, clip_on=False)
        ax.scatter([xk], [0.10], s=SIZE_PER_PCT * k, transform=ax.transAxes, facecolors="none",
                   edgecolors="0.3", linewidths=0.6, zorder=8, clip_on=False)
        r_ax = np.sqrt(SIZE_PER_PCT * k) / 2 / AXW_PT   # marker radius = sqrt(s)/2 points
        ax.text(xk + r_ax + 0.012, 0.10, f"{k} %", transform=ax.transAxes, fontsize=6.3,
                va="center", ha="left", zorder=8)
    ax.text(0.02, 0.97, rf"$n_{{\mathrm{{eff}}}}$ = ($\Sigma w$)$^2$/$\Sigma w^2$ = {fr.n_eff:.1f} of {fr.n}",
            transform=ax.transAxes, fontsize=6.5, va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.55", lw=0.6))

    # side table, sorted by weight
    order = np.argsort(-pct)
    rows = [[e.loc[i, "SampleID"], f"{excess[i]:.0f}", f"{sigma[i]:.0f}", f"{pct[i]:.1f}"] for i in order]
    tb = tx.table(cellText=rows, colLabels=["Run", "Excess\n(°C)", r"$\sigma_i$" + "\n(°C)", "Weight\n(%)"],
                  loc="upper center", cellLoc="center", colLoc="center",
                  colWidths=[0.30, 0.22, 0.22, 0.24], bbox=[0.0, 0.17, 1.0, 0.83])
    tb.auto_set_font_size(False); tb.set_fontsize(6.2 if fr.n <= 13 else max(4.0, 6.2 * 13 / fr.n))
    for (rr, cc), cell in tb.get_celld().items():
        cell.set_linewidth(0.3); cell.set_edgecolor("0.6")
        if rr == 0:
            cell.set_facecolor("0.92"); cell.set_text_props(fontweight="bold")
        else:
            cell.set_facecolor("white")
            if cc == 0:
                refn = e.loc[order[rr - 1], "Reference"]
                cell.get_text().set_color(styles[refn]["color"])
                cell.get_text().set_fontweight("bold")
    thr = fr.thresholds
    n_fail = int((wt["screens_failed"] != "").sum())
    fail_txt = (f"All {fr.n} runs fail at least one screen." if n_fail == fr.n
                else f"{n_fail} of {fr.n} runs fail at least one screen.")
    if wt["sigma_analytical_is_default"].all():
        sig_txt = rf"$\sigma_i^2$ = {fr.sigma0:g}$^2$ + excess$^2$;"
    else:
        sig_txt = rf"$\sigma_i^2$ = $\sigma_{{a,i}}^2$ + excess$^2$ ($\sigma_a$ = {fr.sigma0:g} °C unless given);"
    tx.text(0.0, 0.155,
            "Excess = quadrature sum of the amounts by which\n"
            f"the equilibrium screens are exceeded (opx–cpx {thr['dT_opx_cpx']:g},\n"
            f"pyx–gt {thr['dT_pyx_gt']:g}, gt–ol {thr['dT_gt_ol']:g} °C); " + sig_txt + "\n"
            r"$w_i$ = 1/$\sigma_i^2$, shown as % of $\Sigma w$." "\n" + fail_txt,
            fontsize=5.8, va="top", ha="left", transform=tx.transAxes, linespacing=1.3)
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"{STEM}.png", dpi=600)
    fig.savefig(outdir / f"{STEM}.pdf")
    plt.close(fig)
    print(f"  wrote {outdir / STEM}.png / .pdf")
    return order


def write_weights_csv(runs, wt, fr, order, qf, outdir):
    out = pd.DataFrame({
        "SampleID": wt["SampleID"], "Reference": wt["Reference"],
        "T_TA98_C": runs["T_TA98"].round(1), "dT_TA98_minus_WZ07_Fe2only_C": runs["y_TA98_minus_WZ07_C"].round(1),
        "dT_opx_cpx_C": runs["dT_opx_cpx"].round(1), "dT_pyx_gt_C": runs["dT_pyx_gt"].round(1),
        "dT_gt_ol_C": runs["dT_gt_ol"].round(1), "screen_excess_C": wt["screen_excess_C"].round(1),
        "sigma_i_C": wt["sigma_i_C"].round(1), "w_raw_per_C2": wt["w_raw_per_C2"],
        "weight_pct": wt["weight_pct"].round(2),
        "marker_size_s_pt2": (SIZE_PER_PCT * wt["weight_pct"]).round(1),
        "quality_factor_sensitivity_only": qf,
    }).iloc[order]
    path = outdir / "Fig04_weights.csv"
    out.to_csv(path, index=False)
    n_dflt = int(wt["sigma_analytical_is_default"].sum())
    with open(path, "a") as fh:
        fh.write(f"# n_eff=(sum w)^2/sum w^2={fr.n_eff:.3f}; sigma0={fr.sigma0:g} C "
                 f"(default for {n_dflt} of {fr.n} runs); weights = primary item-1 inverse-variance "
                 f"(disequilibrium only); marker_size_s_pt2 is matplotlib s (diameter^2, pt^2), drawn circle "
                 f"area = pi/4*s; quality_factor is the oxide-total placeholder used only in the sensitivity run\n")
    print(f"  wrote {path}")


def main(argv=None):
    p = argparse.ArgumentParser(description="Fig. 4 weighted-regression figure (deg C input).")
    p.add_argument("--input", default=str(wtr.DEFAULT_INPUT), help="run CSV, schema as weighted_thermometer_regression.py")
    p.add_argument("--outdir", default=None, help="output folder (default: $REV1_OUT, ../figures_revision_v1, or Figures/)")
    p.add_argument("--col", action="append", default=[], metavar="NAME=CSVCOL")
    p.add_argument("--no-crosscheck", action="store_true", help="skip the pipeline and item-1 cross-checks")
    a = p.parse_args(argv)
    outdir = Path(a.outdir) if a.outdir else default_outdir()
    try:
        runs = wtr.load_runs(a.input, wtr.parse_colmap(a.col))
        wt, fr = wtr.fit(runs)
    except wtr.InputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2)
    if a.no_crosscheck:
        print("notice: cross-checks skipped (--no-crosscheck)")
        qf = np.full(len(wt), np.nan)
    else:
        crosscheck_pipeline(runs, wt, fr)
        qf = crosscheck_item1(wt)
    order = make_figure(runs, wt, fr, outdir)
    write_weights_csv(runs, wt, fr, order, qf, outdir)
    pct = wt["weight_pct"]
    print("n_eff", fr.n_eff, "pct range", pct.min(), pct.max())


if __name__ == "__main__":
    main()
