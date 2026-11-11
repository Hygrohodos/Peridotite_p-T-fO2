"""Fig. 4 (revision v1, Stephanie 11 Nov 2026): weighted single-regression
figure with background circles whose AREA is proportional to each run's
normalised weight, a size key, n_eff, and a side table of the weights.

Weights and fit come from scripts/weighted_thermometer_regression.py (the
stand-alone module; w = 1/(sigma0^2 + excess^2), WLS centred at 1350 C).
They are cross-checked here against calculations.pt_fo2.exp_calibration
(the pipeline's item-1 fit) so the figure and the pipeline cannot drift.

Colour (11 Nov, t1615u): palette taken from the paper's other figures,
scripts/make_figures.py -- no new colours:
  Matjuschkin 2014  #5e3c99 purple  (accent of Figs 2 and 6)
  Rzehak 2020       #4c72b0 blue    (plot_exp_ta98_wz07.STYLE; Fig. 3a)
  fit / 95 % band   black / grey 0.85 (identical to Fig. X b, same fit)
Fonts/line weights: make_figures._apply_journal_style (Arial 8 pt).
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent))
import make_figures_rev1 as r1
import make_figures as mf
import plot_exp_ta98_wz07 as pe
import weighted_thermometer_regression as wtr
from calculations.pt_fo2 import exp_calibration as xc

EXP_STYLE = {  # colours from make_figures.py / plot_exp_ta98_wz07.py
    "Matjuschkin et al. 2014": dict(marker="o", color="#5e3c99", label="Matjuschkin 2014"),
    "Rzehak et al. 2020": dict(marker="s", color="#4c72b0", label="Rzehak 2020"),
}

AREA_PER_PCT = 25.0   # pt^2 of circle area per 1 % of total weight
MARK_S = 16           # fixed marker size (pt^2)
LABEL_OFF = {  # manual label offsets (points) to avoid overlaps
    "F1100o": (5, -9), "F1100r": (6, 4), "F1200o": (7, 5), "F1200r": (7, -3),
    "F1300o": (9, -3), "F1300r": (-31, 4), "F14002o": (7, -8), "F14002r": (6, 4),
    "MaLm6": (6, 4), "MaLm11": (-40, 9), "MaLm1": (9, 7), "MaLm9": (-30, -3),
    "MaLm12": (7, -3),
}


def main():
    mf._apply_journal_style()
    runs = wtr.load_runs()                       # data/derived/exp_TA98_WZ07_Fe3.csv
    wt, fr = wtr.fit(runs)                       # module = single source of the method
    ex, w_x, f = r1.exp_fit()                    # pipeline fit, for the cross-check only
    ex = ex.set_index("SampleID").loc[runs["SampleID"]].reset_index()
    assert np.allclose(xc.weights(ex, "invvar", 50.0), wt["w_raw_per_C2"], rtol=1e-12)
    assert np.allclose(f["beta"], [fr.b0, fr.b1], rtol=1e-10)
    assert np.allclose(np.sqrt(np.diag(f["cov"])), [fr.se_b0, fr.se_b1], rtol=1e-10)
    e = runs.assign(dT_Fe2=runs["y_TA98_minus_WZ07_C"])
    w = wt["w_raw_per_C2"].to_numpy(float)
    pct = wt["weight_pct"].to_numpy(float)
    n_eff = fr.n_eff
    excess = wt["screen_excess_C"].to_numpy(float)
    sigma = wt["sigma_i_C"].to_numpy(float)
    ref = pd.read_csv(r1.HERE.parent.parent / "item1_weighted_calibration/results/experiment_weights.csv")
    ref = ref.set_index("SampleID").loc[e["SampleID"]]
    assert np.allclose(ref["w_rel_primary"].to_numpy(), wt["weight_frac"], atol=6e-5)

    fig = plt.figure(figsize=(174 / 25.4, 95 / 25.4))
    ax = fig.add_axes([0.085, 0.115, 0.575, 0.80])
    tx = fig.add_axes([0.675, 0.115, 0.32, 0.80]); tx.axis("off")

    ax.axhline(0, color="0.6", lw=0.7, zorder=0)
    tt = np.linspace(r1.TLO, r1.THI, 100)
    lo, hi = fr.band(tt)
    ax.fill_between(tt, lo, hi, color="0.85", lw=0, zorder=1, label="95 % confidence band")
    ax.plot(tt, fr.predict(tt), color="k", lw=1.4, zorder=4, label=f"Weighted fit ({fr.n} experiments)")
    for ref_, sty in EXP_STYLE.items():
        s = (e["Reference"] == ref_).to_numpy()
        if not s.any():
            continue
        col = sty["color"]
        ax.scatter(e.loc[s, "T_TA98"], e.loc[s, "dT_Fe2"], s=AREA_PER_PCT * pct[s],
                   marker="o", facecolors=col, alpha=0.20, edgecolors="none", zorder=2)
        ax.scatter(e.loc[s, "T_TA98"], e.loc[s, "dT_Fe2"], s=AREA_PER_PCT * pct[s],
                   marker="o", facecolors="none", edgecolors=col, linewidths=0.7, alpha=0.85, zorder=2)
        ax.scatter(e.loc[s, "T_TA98"], e.loc[s, "dT_Fe2"], s=MARK_S, marker=sty["marker"],
                   color=col, edgecolors="k", linewidths=0.4, zorder=5, label=sty["label"])
    for i, r in e.iterrows():
        lead = r["SampleID"] in ("MaLm11", "MaLm1")  # close pair: short leader lines
        ax.annotate(r["SampleID"], (r["T_TA98"], r["dT_Fe2"]), fontsize=6,
                    xytext=LABEL_OFF.get(r["SampleID"], (5, 4)), textcoords="offset points", zorder=6,
                    arrowprops=dict(arrowstyle="-", lw=0.5, color="0.25", shrinkA=0.5, shrinkB=2.5) if lead else None)
    b0, b1 = fr.b0, fr.b1; se = (fr.se_b0, fr.se_b1)
    ax.set_title(rf"$T_{{\mathrm{{TA98}}}}-T_{{\mathrm{{WZ07}}}}$ = {b0:.0f} ($\pm${se[0]:.0f}) "
                 rf"{b1:+.3f} ($\pm${se[1]:.3f})($T_{{\mathrm{{TA98}}}}$ − 1350)", fontsize=7.5)
    ax.set_xlim(1120, 1660); ax.set_ylim(-330, 230)
    mf._style_ax(ax, "TA98 temperature (°C)", r"$T_{\mathrm{TA98}}-T_{\mathrm{WZ07}}$ (°C), Fe$^{2+}$ only")
    h, l = ax.get_legend_handles_labels()
    leg1 = ax.legend(h, l, frameon=True, facecolor="white", edgecolor="0.55", framealpha=0.92,
                     fontsize=6.5, loc="upper right")
    ax.add_artist(leg1)
    # size key (drawn in axes coordinates): circle area proportional to weight
    from matplotlib.patches import FancyBboxPatch
    keys = (5, 10, 20)
    ax.add_patch(FancyBboxPatch((0.02, 0.025), 0.42, 0.20, boxstyle="round,pad=0.005,rounding_size=0.01",
                                transform=ax.transAxes, fc="white", ec="0.55", lw=0.6, alpha=0.92, zorder=7))
    ax.text(0.23, 0.205, "Weight (circle area, % of total)", transform=ax.transAxes, fontsize=6.3,
            ha="center", va="top", zorder=8)
    AXW_PT = 0.575 * 174 / 25.4 * 72  # axes width in points
    for k, xk in zip(keys, (0.065, 0.19, 0.335)):
        ax.scatter([xk], [0.10], s=AREA_PER_PCT * k, transform=ax.transAxes, facecolors="0.5",
                   alpha=0.25, edgecolors="none", zorder=8, clip_on=False)
        ax.scatter([xk], [0.10], s=AREA_PER_PCT * k, transform=ax.transAxes, facecolors="none",
                   edgecolors="0.3", linewidths=0.6, zorder=8, clip_on=False)
        r_ax = np.sqrt(AREA_PER_PCT * k / np.pi) / AXW_PT
        ax.text(xk + r_ax + 0.012, 0.10, f"{k} %", transform=ax.transAxes, fontsize=6.3,
                va="center", ha="left", zorder=8)
    ax.text(0.02, 0.97, rf"$n_{{\mathrm{{eff}}}}$ = ($\Sigma w$)$^2$/$\Sigma w^2$ = {n_eff:.1f} of {fr.n}",
            transform=ax.transAxes, fontsize=6.5, va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.55", lw=0.6))

    # side table, sorted by weight
    order = np.argsort(-pct)
    rows = [[e.loc[i, "SampleID"], f"{excess[i]:.0f}", f"{sigma[i]:.0f}", f"{pct[i]:.1f}"] for i in order]
    tb = tx.table(cellText=rows, colLabels=["Run", "Excess\n(°C)", r"$\sigma_i$" + "\n(°C)", "Weight\n(%)"],
                  loc="upper center", cellLoc="center", colLoc="center",
                  colWidths=[0.30, 0.22, 0.22, 0.24], bbox=[0.0, 0.17, 1.0, 0.83])
    tb.auto_set_font_size(False); tb.set_fontsize(6.2)
    for (rr, cc), cell in tb.get_celld().items():
        cell.set_linewidth(0.3); cell.set_edgecolor("0.6")
        if rr == 0:
            cell.set_facecolor("0.92"); cell.set_text_props(fontweight="bold")
        else:
            cell.set_facecolor("white")
            if cc == 0:
                refn = e.loc[order[rr - 1], "Reference"]
                cell.get_text().set_color(EXP_STYLE[refn]["color"])
                cell.get_text().set_fontweight("bold")
    tx.text(0.0, 0.155,
            "Excess = quadrature sum of the amounts by which\n"
            "the equilibrium screens are exceeded (opx–cpx 80,\n"
            "pyx–gt 60, gt–ol 50 °C); " r"$\sigma_i^2$ = 50$^2$ + excess$^2$;" "\n"
            r"$w_i$ = 1/$\sigma_i^2$, shown as % of $\Sigma w$." "\n"
            "All 13 runs fail at least one screen.",
            fontsize=5.8, va="top", ha="left", transform=tx.transAxes, linespacing=1.3)
    r1.mf._save(fig, "Fig04_TA98_WZ07_measured_Fe3", tight=False)

    out = pd.DataFrame({
        "SampleID": e["SampleID"], "Reference": e["Reference"],
        "T_TA98_C": e["T_TA98"].round(1), "dT_TA98_minus_WZ07_Fe2only_C": e["dT_Fe2"].round(1),
        "dT_opx_cpx_C": e["dT_opx_cpx"].round(1), "dT_pyx_gt_C": e["dT_pyx_gt"].round(1),
        "dT_gt_ol_C": e["dT_gt_ol"].round(1), "screen_excess_C": excess.round(1),
        "sigma_i_C": sigma.round(1), "w_raw_per_C2": w, "weight_pct": pct.round(2),
        "circle_area_pt2": (AREA_PER_PCT * pct).round(1),
        "quality_factor_sensitivity_only": ref["quality_factor"].to_numpy(),
    }).iloc[order]
    out.to_csv(r1.OUT / "Fig04_weights.csv", index=False)
    with open(r1.OUT / "Fig04_weights.csv", "a") as fh:
        fh.write(f"# n_eff=(sum w)^2/sum w^2={n_eff:.3f}; sigma0=50 C; weights = primary item-1 inverse-variance "
                 f"(disequilibrium only); quality_factor is the oxide-total placeholder used only in the sensitivity run\n")
    print("n_eff", n_eff, "pct range", pct.min(), pct.max())


if __name__ == "__main__":
    main()
