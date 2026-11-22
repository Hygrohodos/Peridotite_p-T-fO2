#!/usr/bin/env python3
"""Fig. 5: the 13 experiments and their weighted single regression, with background circles whose AREA is
proportional to each run's normalised weight, a size key and n_eff; and the Online Resource 1 table of the
runs' screen excess, sigma_i and weight (formerly drawn beside the plot; moved to the supplement, 14 Nov 2026).

    python3 scripts/make_fig05_weightcircles.py --pivot 1207.0                  # current data
    python3 scripts/make_fig05_weightcircles.py --input runs.csv --outdir out/  # pivot from a pipeline run

Method: weights and fit come only from scripts/weighted_thermometer_regression.py
(w = 1/(sigma_a^2 + excess^2), WLS centred at 1350 C); the figure adds nothing. The title equation is the same
fit written about the pivot of Eq. 4 (median TA98 of the screened Slave calibration samples, 1207.0 C):
refit_about() refits with t_ref = pivot and asserts the same slope, SE, line and 95 % band.

Style: scripts/paper_figures/unified.py (fonts, ticks, line widths, frameless key), as Fig. 4.
Runs from a clean checkout of the repository. Needs: numpy, pandas, matplotlib, and the repo itself.
Optional cross-checks, each skipped with a notice when it cannot apply (--no-crosscheck skips both):
  * pipeline: rebuilds the 13 runs from the raw experimental data through
    calculations.pt_fo2.exp_calibration and compares weights, coefficients
    and the SAME-convention covariance (both scaled by max(1, chi2_red)).
    Only when the input's SampleIDs equal the pipeline's usable runs.
  * item 1: data/reference/item1_experiment_weights.csv (in-repo copy of
    ../item1_weighted_calibration/results/experiment_weights.csv; also supplies the sensitivity-only quality_factor
    column of the weights CSV, else NaN).

Output (--outdir; default $REV1_OUT, else ../figures_revision_v1 next to the repo if it exists, else Figures/):
Fig05_TA98_WZ07_measured_Fe3.png/.pdf, Fig05_TA98_WZ07_measured_Fe3_equation.json (coefficients at 1350 C and at
the pivot), Fig05_weights.csv, and Online_Resource_1/Online_Resource_1_experiment_weights.csv/.md.

Colour (11 Nov, t1615u): palette taken from the paper's other figures, no new colours:
  Matjuschkin 2014  #5e3c99 purple  (accent colour of the original-submission figures)
  Rzehak 2020       #009E73 green, diamond (15 Nov, t1832u: was #4c72b0 blue square, confusable with the Siberia craton symbol)
  other references  FALLBACK_COLOURS, in order
  fit / 95 % band   black / grey 0.85
"""
from __future__ import annotations
import argparse
import json
import os
import re
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(REPO))
import weighted_thermometer_regression as wtr

EXP_STYLE = {  # colours from the original-submission figures / plot_exp_ta98_wz07.py
    "Matjuschkin et al. 2014": dict(marker="o", color="#5e3c99", label="Matjuschkin 2014"),
    "Rzehak et al. 2020": dict(marker="D", color="#009E73", label="Rzehak 2020"),
}
# for references not in EXP_STYLE (all used elsewhere in make_figures.py)
FALLBACK_COLOURS = ("#b35806", "#55a868", "#c44e52", "#d4a017", "#8c5a2b", "#6b4c9a", "#8c8c8c")
FALLBACK_MARKERS = ("^", "D", "v", "P", "X", "h", "<")

SIZE_PER_PCT = 14.0   # matplotlib marker size s (pt^2, = diameter^2) per 1 % of total weight;
                      # drawn circle area = (pi/4) * s, so area is proportional to weight
MARK_S = 16           # fixed marker size s (pt^2)
LABEL_OFF = {  # manual label offsets (points) to avoid overlaps; others default to (5, 4)
    "F1100o": (-12, -14), "F1100r": (6, 4), "F1200o": (7, 5), "F1200r": (7, -3),
    "F1300o": (9, -3), "F1300r": (-36, -11), "F14002o": (7, -8), "F14002r": (-8, 8),
    "MaLm6": (6, 4), "MaLm11": (-40, 9), "MaLm1": (-12, 17), "MaLm9": (-30, -3),
    "MaLm12": (7, -3),
}
LEADER = ("MaLm11", "MaLm1")   # close pair: short leader lines
XLIM, YLIM = (1120, 1660), (-330, 230)   # expanded automatically if data fall outside
ITEM1 = REPO / "data" / "reference" / "item1_experiment_weights.csv"   # in-repo copy of ../item1_weighted_calibration/results/experiment_weights.csv (rev. v3)
STEM = "Fig05_TA98_WZ07_measured_Fe3"


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


def refit_about(runs, fr, pivot):
    """The same weighted fit written about the pivot: T_TA98 - T_WZ07 = b0' + b1 (T_TA98 - pivot).

    Refitted with t_ref = pivot and checked against fr: same slope and SE, b0' = b0 + b1 (pivot - t_ref),
    and the same line and 95 % band at every temperature (a change of origin only)."""
    fp = wtr.fit(runs, t_ref=pivot)[1]
    tt = np.linspace(1000.0, 1700.0, 71)
    assert np.isclose(fp.b1, fr.b1, rtol=1e-10) and np.isclose(fp.se_b1, fr.se_b1, rtol=1e-10)
    assert np.isclose(fp.b0, fr.b0 + fr.b1 * (pivot - fr.t_ref), rtol=1e-10, atol=1e-9)
    assert np.allclose(fp.predict(tt), fr.predict(tt), atol=1e-9)
    assert np.allclose(np.array(fp.band(tt)), np.array(fr.band(tt)), atol=1e-9)
    return fp


def make_figure(runs, wt, fr, outdir, pivot):
    """Fig. 5 in the shared paper-figure style (scripts/paper_figures/unified.py, as Fig. 4); the run/weight
    table is not drawn (it is Online Resource 1, write_supplement_table). Returns the run order by weight."""
    sys.path.insert(0, str(HERE / "paper_figures"))
    import unified as uf
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MultipleLocator, FuncFormatter
    fp = refit_about(runs, fr, pivot)
    e = runs.assign(dT_Fe2=runs["y_TA98_minus_WZ07_C"])
    pct = wt["weight_pct"].to_numpy(float)
    styles = reference_styles(e["Reference"].tolist())

    W, H = uf.COL1, 97
    AX = (15.5, 11.0, 66.0, 66.0)   # axes left, bottom, width, height (mm)
    fig = uf.newfig(W, H)
    ax = uf.axmm(fig, W, H, *AX)
    ax.axhline(0, color="0.6", lw=uf.LW_THIN, zorder=0)
    # fit drawn over the data range, rounded to 0.01 C as params.wz07_offset_tmin/tmax_c
    tt = np.linspace(round(e["T_TA98"].min(), 2), round(e["T_TA98"].max(), 2), 100)
    lo, hi = fp.band(tt)
    ax.fill_between(tt, lo, hi, color="0.85", lw=0, zorder=1)
    ax.plot(tt, fp.predict(tt), color="k", lw=uf.LW_MAIN, zorder=4)
    for ref_, sty in styles.items():
        s = (e["Reference"] == ref_).to_numpy()
        col = sty["color"]
        ax.scatter(e.loc[s, "T_TA98"], e.loc[s, "dT_Fe2"], s=SIZE_PER_PCT * pct[s],
                   marker="o", facecolors=col, alpha=0.18, edgecolors="none", zorder=2)
        ax.scatter(e.loc[s, "T_TA98"], e.loc[s, "dT_Fe2"], s=SIZE_PER_PCT * pct[s],
                   marker="o", facecolors="none", edgecolors=col, linewidths=0.6, alpha=0.85, zorder=2)
        ax.scatter(e.loc[s, "T_TA98"], e.loc[s, "dT_Fe2"], s=uf.S_STD, marker=sty["marker"],
                   color=col, edgecolors="k", linewidths=uf.EW_STD, zorder=5)
    for i, r in e.iterrows():
        lead = r["SampleID"] in LEADER
        ax.annotate(r["SampleID"], (r["T_TA98"], r["dT_Fe2"]), fontsize=uf.FS_SMALL,
                    xytext=LABEL_OFF.get(r["SampleID"], (5, 4)), textcoords="offset points", zorder=6,
                    arrowprops=dict(arrowstyle="-", lw=uf.LW_THIN, color="0.25", shrinkA=0.5, shrinkB=2.5) if lead else None)
    eq = uf.M(rf"$T_\mathrm{{TA98}}-T_\mathrm{{WZ07}}$ = {fp.b0:.0f} (±{fp.se_b0:.0f})" "\n"
              rf"{fp.b1:+.3f} (±{fp.se_b1:.3f})($T_\mathrm{{TA98}}$ - {pivot:.1f})")   # two lines to fit one column
    fig.text((AX[0] + AX[2] / 2) / W, 1 - 1.5 / H, eq, ha="center", va="top", fontsize=uf.FS_LABEL, linespacing=1.3)
    ax.set_xlim(*_lims(e["T_TA98"], XLIM, 40)); ax.set_ylim(*_lims(e["dT_Fe2"], YLIM, 50))
    ax.xaxis.set_major_locator(MultipleLocator(100)); ax.xaxis.set_minor_locator(MultipleLocator(25))
    ax.yaxis.set_major_locator(MultipleLocator(100)); ax.yaxis.set_minor_locator(MultipleLocator(25))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: uf.M(f"{v:.0f}")))
    ax.set_xlabel(uf.XT_TA98)
    ax.set_ylabel(uf.DT_LABEL + ", Fe$^{2{+}}$ only")
    key = [Line2D([], [], color="k", lw=uf.LW_MAIN, label=f"weighted fit (n = {fr.n})"),
           plt.Rectangle((0, 0), 1, 1, fc="0.85", ec="none", label="95 % confidence band")]
    key += [Line2D([], [], marker=sty["marker"], ls="", mfc=sty["color"], mec="k", mew=uf.EW_STD, ms=uf.MS_KEY,
                   label=sty["label"]) for sty in styles.values()]
    # key above the axes, two columns: fit and band | references
    ax.legend(handles=key, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=uf.FS_KEY,
              borderaxespad=0.5, columnspacing=1.2)
    ax.text(0.02, 0.975, uf.M(rf"$n_\mathrm{{eff}}$ = ($\Sigma w$)$^2$/$\Sigma w^2$ = {fr.n_eff:.1f} of {fr.n}"),
            transform=ax.transAxes, fontsize=uf.FS_KEY, va="top", ha="left")
    # size key (axes coordinates): circle area proportional to weight
    ax.text(0.015, 0.155, "Weight (circle area, % of total)", transform=ax.transAxes, fontsize=uf.FS_SMALL,
            ha="left", va="top", zorder=8)
    axw_pt = AX[2] / 25.4 * 72  # axes width in points
    for k, xk in zip((5, 10, 20), np.array((2.5, 11.5, 23.0)) / AX[2]):   # circle centres, mm from the axes' left edge
        ax.scatter([xk], [0.07], s=SIZE_PER_PCT * k, transform=ax.transAxes, facecolors="0.5",
                   alpha=0.25, edgecolors="none", zorder=8, clip_on=False)
        ax.scatter([xk], [0.07], s=SIZE_PER_PCT * k, transform=ax.transAxes, facecolors="none",
                   edgecolors="0.3", linewidths=0.6, zorder=8, clip_on=False)
        r_ax = np.sqrt(SIZE_PER_PCT * k) / 2 / axw_pt   # marker radius = sqrt(s)/2 points
        ax.text(xk + r_ax + 0.6 / AX[2], 0.07, f"{k} %", transform=ax.transAxes, fontsize=uf.FS_SMALL,
                va="center", ha="left", zorder=8)
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"{STEM}.png", dpi=600)
    fig.savefig(outdir / f"{STEM}.pdf")
    plt.close(fig)
    print(f"  wrote {outdir / STEM}.png / .pdf")
    stats = dict(n=fr.n, n_eff=fr.n_eff, slope=fr.b1, se_slope=fr.se_b1, t_ref_fit=fr.t_ref, b0_at_t_ref=fr.b0,
                 se_b0_at_t_ref=fr.se_b0, pivot=pivot, b0_at_pivot=fp.b0, se_b0_at_pivot=fp.se_b0)
    (outdir / f"{STEM}_equation.json").write_text(json.dumps(stats, indent=1, default=float))
    return np.argsort(-pct)


def write_supplement_table(runs, wt, fr, order, outdir):
    """Online Resource 1 table: the run/weight table formerly drawn beside Fig. 5, with its notes."""
    styles = reference_styles(runs["Reference"].tolist())
    thr = fr.thresholds
    t = pd.DataFrame({
        "Run": wt["SampleID"], "Reference": wt["Reference"],
        "TA98 (°C)": runs["T_TA98"].round(0).astype(int),
        "TA98 − WZ07, Fe2+ only (°C)": runs["y_TA98_minus_WZ07_C"].round(0).astype(int),
        "Excess (°C)": wt["screen_excess_C"].round(0).astype(int),
        "σi (°C)": wt["sigma_i_C"].round(0).astype(int),
        "Weight (%)": wt["weight_pct"].map(lambda v: f"{v:.1f}"),
    }).iloc[order]
    d = outdir / "Online_Resource_1"
    d.mkdir(parents=True, exist_ok=True)
    t.to_csv(d / "Online_Resource_1_experiment_weights.csv", index=False, encoding="utf-8")
    n_fail = int((wt["screens_failed"] != "").sum())
    fail = (f"All {fr.n} runs fail at least one screen." if n_fail == fr.n
            else f"{n_fail} of {fr.n} runs fail at least one screen.")
    lines = ["# Online Resource 1: weights of the experimental runs in the fit of Fig. 5", "",
             f"**Table.** The {fr.n} experimental run products with measured garnet Fe3+/ΣFe used in the weighted fit "
             "of Fig. 5 (Section 3.3), in order of weight. TA98 − WZ07 is calculated with the measured garnet Fe3+. "
             "Excess is the quadrature sum of the amounts by which the equilibrium screens of Section 3.2 are exceeded "
             f"(opx–cpx {thr['dT_opx_cpx']:g} °C, pyx–gt {thr['dT_pyx_gt']:g} °C, gt–ol {thr['dT_gt_ol']:g} °C). "
             f"σi² = {fr.sigma0:g}² + excess² (°C²); wi = 1/σi², shown as a percentage of Σw. {fail} "
             f"Effective number of runs n_eff = (Σw)²/Σw² = {fr.n_eff:.1f} of {fr.n}.", "",
             "| " + " | ".join(t.columns) + " |", "|" + "---|" * len(t.columns)]
    lines += ["| " + " | ".join(str(v) for v in r) + " |" for r in t.itertuples(index=False)]
    lines += ["", "Generated by `scripts/make_fig05_weightcircles.py` (called by "
              "`scripts/paper_figures/make_paper_figures.py`); full-precision weights are in `Fig5_weights.csv`.", ""]
    (d / "Online_Resource_1_experiment_weights.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"  wrote {d}")


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
    path = outdir / "Fig05_weights.csv"
    out.to_csv(path, index=False)
    n_dflt = int(wt["sigma_analytical_is_default"].sum())
    with open(path, "a") as fh:
        fh.write(f"# n_eff=(sum w)^2/sum w^2={fr.n_eff:.3f}; sigma0={fr.sigma0:g} C "
                 f"(default for {n_dflt} of {fr.n} runs); weights = primary item-1 inverse-variance "
                 f"(disequilibrium only); marker_size_s_pt2 is matplotlib s (diameter^2, pt^2), drawn circle "
                 f"area = pi/4*s; quality_factor is the oxide-total placeholder used only in the sensitivity run\n")
    print(f"  wrote {path}")


def slave_pivot():
    """Median TA98 of the screened Slave calibration samples (the pivot of Eq. 4), from a pipeline run."""
    sys.path.insert(0, str(HERE / "paper_figures"))
    import make_paper_figures as mp
    return mp.slave_match_level()[1]


def main(argv=None):
    p = argparse.ArgumentParser(description="Fig. 5 weighted-regression figure (deg C input).")
    p.add_argument("--input", default=str(wtr.DEFAULT_INPUT), help="run CSV, schema as weighted_thermometer_regression.py")
    p.add_argument("--outdir", default=None, help="output folder (default: $REV1_OUT, ../figures_revision_v1, or Figures/)")
    p.add_argument("--col", action="append", default=[], metavar="NAME=CSVCOL")
    p.add_argument("--no-crosscheck", action="store_true", help="skip the pipeline and item-1 cross-checks")
    p.add_argument("--pivot", type=float, default=None,
                   help="pivot (deg C) about which the title equation is written; default: the median TA98 of the "
                        "screened Slave calibration samples, from a pipeline run (as make_paper_figures.py)")
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
    pivot = a.pivot if a.pivot is not None else slave_pivot()
    order = make_figure(runs, wt, fr, outdir, pivot)
    write_supplement_table(runs, wt, fr, order, outdir)
    write_weights_csv(runs, wt, fr, order, qf, outdir)
    pct = wt["weight_pct"]
    print("n_eff", fr.n_eff, "pct range", pct.min(), pct.max())


if __name__ == "__main__":
    main()
