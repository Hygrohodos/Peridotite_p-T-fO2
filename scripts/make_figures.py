#!/usr/bin/env python3
"""
Rebuild chapter figures from the calculated xenolith table.

Usage (from the project root):
    python scripts/make_figures.py

Outputs PNG + PDF files to Figures/.
Edit symbol / colour choices below if you want a closer Grapher match.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def find_root() -> Path:
    here = Path(__file__).resolve().parent
    for root in (here, here.parent):
        if (root / "calculations").is_dir():
            return root
    return here


ROOT = find_root()
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from calculations.pt_fo2.params import PARAMS
from calculations.pt_fo2 import calibration as calib

FULL = ROOT / "data" / "derived" / "screened_PT.csv"
OUT = ROOT / "Figures"
_LOC_CANDIDATES = (
    ROOT / "data" / "master" / "kimberlite_locations.csv",
    ROOT / "data" / "kimberlite_locations.csv",
)
LOC = next((p for p in _LOC_CANDIDATES if p.is_file()), _LOC_CANDIDATES[0])


CRATON_STYLE = {
    "Slave Craton": dict(marker="o", color="#e07b00", label="Slave"),
    "Siberian Craton": dict(marker="s", color="#4c72b0", label="Siberian"),
    "Kaapvaal Craton": dict(marker="^", color="#c44e9a", label="Kaapvaal"),
}
OTHER = dict(marker="o", color="#8c8c8c", label="Other", facecolors="none")


GREEN_OUTER = "#b7dfb7"   # ±0.05 Fe3+
GREEN_INNER = "#3d9a4a"   # ±0.02 Fe3+ (darker)
# ΔlogfO2(FMQ) equivalent of a 0.02 / 0.05 change in garnet Fe3+/ΣFe, as
# quoted in the manuscript (0.44 log units per 0.02).  These constants are
# used for every Fe3+ → log-unit conversion (Fig. 10 bars, Table 1 notes).
# main() prints fo2_equivalents(df), the median shift for +0.02 / +0.05 on
# the current screened samples, as a check; it is not substituted.
FO2_EQ_02 = 0.44
FO2_EQ_05 = 0.95


def fo2_equivalents(df):
    """Median ΔFMQ change for +0.02 and +0.05 garnet Fe3+/ΣFe (screened)."""
    use = df[df["screen_pass"] == True]
    x = pd.to_numeric(use["XFe3_gt_corr"], errors="coerce")
    f0 = _fo2_from_xfe3(use, x)
    out = []
    for d in (0.02, 0.05):
        f1 = _fo2_from_xfe3(use, x + d)
        out.append(float(np.nanmedian((f1 - f0).to_numpy(dtype=float))))
    return out[0], out[1]


def _gpa_axis(ax, which="y"):
    """Label a pressure axis held in kbar internally as GPa (1 GPa = 10 kbar)."""
    from matplotlib.ticker import FuncFormatter, MultipleLocator
    axis = ax.yaxis if which == "y" else ax.xaxis
    axis.set_major_locator(MultipleLocator(5.0))
    axis.set_major_formatter(FuncFormatter(lambda v, _: f"{v / 10:.1f}"))
    (ax.set_ylabel if which == "y" else ax.set_xlabel)("P (GPa)")


_GEOTHERM = None


def _load_geotherm():
    """P (kbar) → depth (km) from data/geotherm.xlsx."""
    global _GEOTHERM
    if _GEOTHERM is not None:
        return _GEOTHERM
    path = ROOT / "data" / "geotherm.xlsx"
    g = pd.read_excel(path)
    p = pd.to_numeric(g["P (kbar)"], errors="coerce")
    z = pd.to_numeric(g["z (km)"], errors="coerce")
    m = p.notna() & z.notna()
    order = np.argsort(p[m].to_numpy())
    _GEOTHERM = (p[m].to_numpy()[order], z[m].to_numpy()[order])
    return _GEOTHERM


def depth_km(P_kbar):
    P = np.asarray(P_kbar, dtype=float)
    p, z = _load_geotherm()
    out = np.interp(P, p, z)
    slope_hi = (z[-1] - z[-2]) / (p[-1] - p[-2])
    slope_lo = z[0] / p[0] if p[0] else 0.0
    out = np.where(P > p[-1], z[-1] + slope_hi * (P - p[-1]), out)
    out = np.where(P < p[0], slope_lo * P, out)
    return out


def pressure_from_depth(z_km):
    z_km = np.asarray(z_km, dtype=float)
    p, z = _load_geotherm()
    out = np.interp(z_km, z, p)
    slo = p[0] / z[0] if z[0] else 0.0
    shi = (p[-1] - p[-2]) / (z[-1] - z[-2])
    out = np.where(z_km < z[0], slo * z_km, out)
    out = np.where(z_km > z[-1], p[-1] + shi * (z_km - z[-1]), out)
    return out


def _add_depth_axis(ax, km_interval=None):
    """Right-hand depth (km) axis. km_interval uses independent round ticks."""
    axr = ax.twinx()
    axr.set_ylim(ax.get_ylim())
    if km_interval:
        zlo = float(np.min(depth_km(list(ax.get_ylim()))))
        zhi = float(np.max(depth_km(list(ax.get_ylim()))))
        ticks_km = np.arange(
            np.ceil(zlo / km_interval) * km_interval,
            np.floor(zhi / km_interval) * km_interval + 0.1,
            km_interval,
        )
        axr.set_yticks(pressure_from_depth(ticks_km))
        axr.set_yticklabels([f"{int(z)}" for z in ticks_km])
        minor_km = np.arange(
            np.ceil(zlo / 10.0) * 10.0,
            np.floor(zhi / 10.0) * 10.0 + 0.1,
            10.0,
        )
        axr.set_yticks(pressure_from_depth(minor_km), minor=True)
        axr.tick_params(which="minor", direction="in", length=2.0)
    else:
        ticks = [t for t in ax.get_yticks() if
                 min(ax.get_ylim()) <= t <= max(ax.get_ylim())]
        axr.set_yticks(ticks)
        axr.set_yticklabels([f"{float(depth_km(t)):.0f}" for t in ticks])
    axr.set_ylabel("Depth (km)")
    axr.tick_params(direction="in")
    axr.spines["top"].set_visible(True)
    return axr


def _label_on_curve(ax, x, y, text, color, y_at, dx=0.10):
    """Place a label immediately to the right of a curve, at the nearest y."""
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if not m.any():
        return
    i = int(np.argmin(np.abs(y[m] - y_at)))
    ax.text(float(x[m][i]) + dx, float(y[m][i]), text,
            color=color, fontsize=8, va="center", ha="left", clip_on=True, zorder=8)


def _panel_label(ax, letter):
    ax.text(0.02, 0.97, f"{letter}.", transform=ax.transAxes,
            fontsize=9, fontweight="bold", va="top", ha="left", zorder=10)


def _style_ax(ax, xlabel, ylabel, title=None):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.tick_params(direction="in", top=True, right=True)
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    if title:
        letter = title.strip()[:2]
        if letter[:1].isalpha() and letter.endswith("."):
            _panel_label(ax, letter[0])


def _fe3_match_band(ax, kind="1to1"):
    """Light ±0.05 and darker ±0.02 agreement bands."""
    if kind == "1to1":
        ax.fill_between([0, 0.35], [-0.05, 0.30], [0.05, 0.40],
                        color=GREEN_OUTER, lw=0, zorder=0, label=r"$\pm$0.05")
        ax.fill_between([0, 0.35], [-0.02, 0.33], [0.02, 0.37],
                        color=GREEN_INNER, lw=0, zorder=0, label=r"$\pm$0.02")
    else:
        ax.axhspan(-0.05, 0.05, color=GREEN_OUTER, lw=0, zorder=0, label=r"$\pm$0.05")
        ax.axhspan(-0.02, 0.02, color=GREEN_INNER, lw=0, zorder=0, label=r"$\pm$0.02")


def _fo2_equiv_scale(ax, x0=-0.2, p=24.0):
    """Horizontal scale: fO2 width equivalent to Fe3+ ±0.02 and ±0.05."""
    ax.plot([x0 - FO2_EQ_05, x0 + FO2_EQ_05], [p + 2.2, p + 2.2],
            color=GREEN_OUTER, lw=7, solid_capstyle="butt", zorder=5)
    ax.plot([x0 - FO2_EQ_02, x0 + FO2_EQ_02], [p, p],
            color=GREEN_INNER, lw=7, solid_capstyle="butt", zorder=5)
    ax.text(x0 + FO2_EQ_05 + 0.08, p + 2.2,
            rf"$\Delta$Fe$^{{3+}}$ $\pm$0.05 $\approx$ {FO2_EQ_05:.2f} log",
            fontsize=6, va="center", color="0.2")
    ax.text(x0 + FO2_EQ_02 + 0.08, p,
            rf"$\Delta$Fe$^{{3+}}$ $\pm$0.02 $\approx$ {FO2_EQ_02:.2f} log",
            fontsize=6, va="center", color="0.2")


def _scatter_by_craton(ax, df, x, y, screened_only=False, legend_loc="lower left"):
    data = df.copy()
    if screened_only and "screen_pass" in data.columns:
        data = data[data["screen_pass"] == True]
    used = set()
    for craton, sty in CRATON_STYLE.items():
        sub = data[data["Craton"].fillna("").astype(str) == craton]
        if sub.empty:
            continue
        ax.scatter(sub[x], sub[y], s=28, edgecolors="k", linewidths=0.4,
                   marker=sty["marker"], c=sty["color"], label=sty["label"], zorder=3)
        used.add(craton)
    other = data[~data["Craton"].fillna("").astype(str).isin(CRATON_STYLE)]
    if not other.empty:
        ax.scatter(other[x], other[y], s=18, edgecolors=OTHER["color"],
                   linewidths=0.6, marker="o", facecolors="none", label="Other", zorder=2)
    ax.legend(frameon=True, facecolor="white", edgecolor="0.55", framealpha=0.92, fontsize=7, loc=legend_loc)


def _apply_journal_style():
    """CMP-ready fonts and line weights at print size (no later rescaling)."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 7,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "lines.linewidth": 1.15,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def _save(fig, stem, tight=True, dpi=600):
    OUT.mkdir(parents=True, exist_ok=True)
    if tight:
        fig.tight_layout()
    fig.savefig(OUT / f"{stem}.png", dpi=dpi)
    fig.savefig(OUT / f"{stem}.pdf")
    plt.close(fig)
    print(f"  wrote {stem}.png / .pdf")


def _wz07_sweep(r, P, xs):
    from calculations.pt_fo2.fe3_correction import _wz07corr_fe2only
    kd = (float(r["garnet_Fe"]) * (1.0 - xs) / float(r["garnet_Mg"])) / (
        float(r["olivine_Fe"]) / float(r["olivine_Mg"])
    )
    T = _wz07corr_fe2only(
        np.full_like(xs, float(r["olivine_Mg"])),
        np.full_like(xs, float(r["olivine_Fe"])),
        np.full_like(xs, float(r["garnet_Mg"])),
        float(r["garnet_Fe"]) * (1.0 - xs),
        np.full_like(xs, float(r["garnet_Ca"])),
        np.full_like(xs, float(r["garnet_Mn"])),
        np.full_like(xs, P),
        offset_c=0.0,
    )
    return kd, T


def fig02_kd_fe3(df):
    """Fig. 2 stacked: (a) KD and ΔT_WZ07; (b) TA98 − WZ07 for 21-2 and 22-1.

    CMP 119 mm wide so type and lines print at this size.
    """
    row_hi = df[df["SampleID"].astype(str).str.replace(" ", "") == "21-2"]
    row_lo = df[df["SampleID"].astype(str).str.replace(" ", "") == "22-1"]
    if row_hi.empty or row_lo.empty:
        print("  skip Fig02 (21-2 or 22-1 missing)")
        return
    hi, lo = row_hi.iloc[0], row_lo.iloc[0]
    xs = np.linspace(0.0, 0.20, 81)
    kd_hi, T_hi = _wz07_sweep(hi, 55.0, xs)
    kd_lo, T_lo = _wz07_sweep(lo, float(lo["P_CA91"]), xs)
    terr = T_hi - T_hi[0]
    dT_hi = float(hi["T_TA98"]) - T_hi
    dT_lo = float(lo["T_TA98"]) - T_lo

    # CMP single column: 84 mm wide
    fig, (ax_a, ax_b) = plt.subplots(
        2, 1, figsize=(84 / 25.4, 107 / 25.4), sharex=True,
    )
    fig.subplots_adjust(hspace=0.18, left=0.22, right=0.76, top=0.96, bottom=0.12)

    # --- a ---
    ax_a.plot(xs, kd_hi, color="k", lw=1.15)
    ax_a.set_ylim(1.60, 2.20)
    ax_a.set_ylabel(r"$K_{\mathrm{D}}$", labelpad=6)
    ax_a2 = ax_a.twinx()
    ax_a2.plot(xs, terr, color="#5e3c99", lw=1.15)
    ax_a2.set_ylim(0, 250)
    ax_a2.set_yticks([0, 50, 100, 150, 200, 250])
    ax_a2.set_ylabel(r"$T(x)-T(x=0)$ (°C)", labelpad=6)
    ax_a.tick_params(direction="in", top=True, pad=4)
    ax_a2.tick_params(direction="in", pad=3)
    ax_a.text(0.02, 2.08, r"$K_{\mathrm{D}}$", color="k", fontsize=8)
    ax_a.text(0.105, 1.68, r"$T(x)-T(x=0)$", color="#5e3c99", fontsize=8)
    ax_a.text(0.02, 0.96, "a.", transform=ax_a.transAxes, fontsize=9,
              fontweight="bold", va="top")

    # --- b ---
    ax_b.plot(xs, kd_hi, color="k", lw=1.15, ls="--")
    ax_b.plot(xs, kd_lo, color="#b35806", lw=1.15, ls="--")
    kd_min = min(kd_hi.min(), kd_lo.min())
    kd_max = max(kd_hi.max(), kd_lo.max())
    pad = 0.08 * (kd_max - kd_min)
    ax_b.set_ylim(kd_min - pad, kd_max + pad)
    ax_b.set_ylabel(r"$K_{\mathrm{D}}$")
    ax_b2 = ax_b.twinx()
    ax_b2.axhline(0.0, color="0.55", lw=0.6, ls=":")
    ax_b2.plot(xs, dT_hi, color="k", lw=1.15)
    ax_b2.plot(xs, dT_lo, color="#b35806", lw=1.15)
    ax_b2.set_ylim(-300, 140)
    ax_b2.set_ylabel(r"$T_{\mathrm{TA98}}-T_{\mathrm{WZ07}}$ (°C)")
    ax_b.set_xlabel(r"Fe$^{3+}$/$\Sigma$Fe removed", labelpad=6)
    ax_b.set_xlim(0, 0.20)
    ax_b.set_xticks([0, 0.05, 0.10, 0.15, 0.20])
    ax_b.tick_params(direction="in", top=True, pad=4)
    ax_b2.tick_params(direction="in", pad=3)
    ax_b.xaxis.labelpad = 6
    ax_b2.text(0.045, float(dT_hi[18]) + 8, f"{float(hi['T_TA98']):.0f} °C",
               color="k", fontsize=8, ha="left", va="bottom")
    ax_b2.text(0.055, float(dT_lo[22]) - 28, f"{float(lo['T_TA98']):.0f} °C",
               color="#b35806", fontsize=8, ha="left", va="top")
    ax_b.text(0.015, 0.97, "b.", transform=ax_b.transAxes, fontsize=9,
              fontweight="bold", va="top")

    _save(fig, "Fig02_KD_WZ07_Fe3", tight=False)


def fig02b_ta98_minus_wz07(df):
    """Panel b of Fig. 2; kept as a no-op for older calls."""
    print("  Fig02b is panel b of Fig02")
    return
    from calculations.pt_fo2.fe3_correction import _wz07corr_fe2only

    specs = [
        ("21-2", 55.0, "#2ca02c"),
        ("22-1", None, "#d62728"),
    ]
    xs = np.linspace(0.0, 0.20, 81)
    series = []
    for sid, P_fix, color in specs:
        row = df[df["SampleID"].astype(str).str.replace(" ", "") == sid]
        if row.empty:
            print(f"  skip Fig02b ({sid} not found)")
            return
        r = row.iloc[0]
        P = float(P_fix) if P_fix is not None else float(r["P_CA91"])
        T_ta98 = float(r["T_TA98"])
        kd = (float(r["garnet_Fe"]) * (1.0 - xs) / float(r["garnet_Mg"])) / (
            float(r["olivine_Fe"]) / float(r["olivine_Mg"])
        )
        T_wz = _wz07corr_fe2only(
            np.full_like(xs, float(r["olivine_Mg"])),
            np.full_like(xs, float(r["olivine_Fe"])),
            np.full_like(xs, float(r["garnet_Mg"])),
            float(r["garnet_Fe"]) * (1.0 - xs),
            np.full_like(xs, float(r["garnet_Ca"])),
            np.full_like(xs, float(r["garnet_Mn"])),
            np.full_like(xs, P),
            offset_c=0.0,
        )
        series.append((sid, P, T_ta98, kd, T_ta98 - T_wz, color))

    dT_all = np.concatenate([s[4] for s in series])
    kd_all = np.concatenate([s[3] for s in series])
    dT_lo = np.floor(np.nanmin(dT_all) / 50.0) * 50.0
    dT_hi = np.ceil(np.nanmax(dT_all) / 50.0) * 50.0
    kd_pad = 0.08 * (float(np.nanmax(kd_all)) - float(np.nanmin(kd_all)))

    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    ax2 = ax.twinx()
    ax2.axhline(0.0, color="0.5", lw=0.6, ls=":")
    for sid, P, T_ta98, kd, dT, color in series:
        ax.plot(xs, kd, color=color, lw=1.3, ls="--")
        ax2.plot(xs, dT, color=color, lw=1.4)
        ax2.text(0.02, dT[8], f"{T_ta98:.0f} °C", color=color,
                 fontsize=9, ha="left", va="bottom")

    ax.set_xlim(0, 0.20)
    ax.set_ylim(float(np.nanmin(kd_all)) - kd_pad,
                float(np.nanmax(kd_all)) + kd_pad)
    ax2.set_ylim(dT_lo, dT_hi)
    ax.set_xlabel(r"Fe$^{3+}$/$\Sigma$Fe removed")
    ax.set_ylabel(r"$K_{\mathrm{D}}$")
    ax2.set_ylabel(r"$T_{\mathrm{TA98}}-T_{\mathrm{WZ07}}$ (°C)", color="k")
    ax.set_xticks([0, 0.05, 0.10, 0.15, 0.20])
    ax.tick_params(direction="in", top=True, pad=10)
    ax.xaxis.labelpad = 8
    ax2.tick_params(direction="in", colors="k")
    ax.set_title("TA98 $-$ WZ07 vs Fe$^{3+}$ removed", loc="left", fontsize=10)
    _save(fig, "Fig02b_TA98_WZ07_Fe3")


def fig04_screens(df):
    """Each panel shows only samples that passed the previous screen."""
    fig, axes = plt.subplots(1, 3, figsize=(174 / 25.4, 62 / 25.4))
    sp = PARAMS.screening
    lim = [700, 1500]
    xx = np.linspace(lim[0], lim[1], 200)

    four = df["screen_4phases"] == True if "screen_4phases" in df.columns else pd.Series(True, index=df.index)
    passed_T = four & (pd.to_numeric(df["T_TA98"], errors="coerce") > sp.t_min_c)
    dT_px = (pd.to_numeric(df["T_TA98"], errors="coerce")
             - pd.to_numeric(df["T_NG10"], errors="coerce")).abs()
    passed_opx_cpx = passed_T & (dT_px < sp.dt_opx_cpx_c)

    TA = pd.to_numeric(df["T_TA98"], errors="coerce")
    HA = pd.to_numeric(df["T_HA84"], errors="coerce")
    WZ = pd.to_numeric(df["T_WZ07"], errors="coerce")
    HA_pred = sp.ha84_poly_a + sp.ha84_poly_b * TA + sp.ha84_poly_c * TA ** 2
    passed_pyx_gt = passed_opx_cpx & ((HA - HA_pred).abs() < sp.dt_pyx_gt_c)

    WZ_pred = 88.99741613 + 0.9864930064 * HA
    passed_gt_ol = passed_pyx_gt & ((WZ - WZ_pred).abs() < sp.dt_gt_ol_c)

    print(f"  Fig04 live screens: T>{sp.t_min_c:.0f}={int(passed_T.sum())}  "
          f"opx–cpx={int(passed_opx_cpx.sum())}  "
          f"pyx–gt={int(passed_pyx_gt.sum())}  "
          f"gt–ol={int(passed_gt_ol.sum())}")

    # a. samples that passed T>900; X = fail opx–cpx
    ax = axes[0]
    use = df[passed_T]
    ax.scatter(use["T_TA98"], use["T_NG10"], s=16, c="#4c72b0",
               edgecolors="k", linewidths=0.25)
    ax.plot(lim, lim, "k-", lw=0.8)
    ax.plot(lim, [v + sp.dt_opx_cpx_c for v in lim], "k--", lw=0.6)
    ax.plot(lim, [v - sp.dt_opx_cpx_c for v in lim], "k--", lw=0.6)
    fail = passed_T & ~passed_opx_cpx
    if fail.any():
        ax.scatter(df.loc[fail, "T_TA98"], df.loc[fail, "T_NG10"],
                   marker="x", c="r", s=28)
    _style_ax(ax, "Two-pyroxene (TA98; °C)", "Ca-in-opx (NG10; °C)", "a. Opx–cpx")
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.text(0.97, 0.04, f"pass {int(passed_opx_cpx.sum())}\nfail {int(fail.sum())}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8)

    # b. samples that passed opx–cpx; X = fail pyx–gt
    ax = axes[1]
    use = df[passed_opx_cpx]
    ax.scatter(use["T_TA98"], use["T_HA84"], s=16, c="#55a868",
               edgecolors="k", linewidths=0.25)
    yy = sp.ha84_poly_a + sp.ha84_poly_b * xx + sp.ha84_poly_c * xx ** 2
    ax.plot(xx, yy, "k-", lw=0.8)
    ax.plot(xx, yy + sp.dt_pyx_gt_c, "k--", lw=0.6)
    ax.plot(xx, yy - sp.dt_pyx_gt_c, "k--", lw=0.6)
    fail = passed_opx_cpx & ~passed_pyx_gt
    if fail.any():
        ax.scatter(df.loc[fail, "T_TA98"], df.loc[fail, "T_HA84"],
                   marker="x", c="r", s=28)
    _style_ax(ax, "Two-pyroxene (TA98; °C)", "Opx–garnet (HA84; °C)", "b. Pyx–garnet")
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.text(0.97, 0.04, f"pass {int(passed_pyx_gt.sum())}\nfail {int(fail.sum())}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8)

    # c. samples that passed pyx–gt; X = fail gt–ol
    ax = axes[2]
    use = df[passed_pyx_gt]
    ax.scatter(use["T_HA84"], use["T_WZ07"], s=16, c="#c44e52",
               edgecolors="k", linewidths=0.25)
    yy = 88.99741613 + 0.9864930064 * xx
    ax.plot(xx, yy, "k-", lw=0.8)
    ax.plot(xx, yy + sp.dt_gt_ol_c, "k--", lw=0.6)
    ax.plot(xx, yy - sp.dt_gt_ol_c, "k--", lw=0.6)
    fail = passed_pyx_gt & ~passed_gt_ol
    if fail.any():
        ax.scatter(df.loc[fail, "T_HA84"], df.loc[fail, "T_WZ07"],
                   marker="x", c="r", s=28)
    _style_ax(ax, "Opx–garnet (HA84; °C)", "Olivine–garnet (WZ07; °C)", "c. Garnet–olivine")
    lim_c = [880, 1480]
    ax.set_xlim(lim_c); ax.set_ylim(lim_c)
    ax.text(0.97, 0.04, f"pass {int(passed_gt_ol.sum())}\nfail {int(fail.sum())}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8)
    _save(fig, "Fig03_equilibrium_screens")


def _fe3_at_dt(df, dt):
    from calculations.pt_fo2.fe3_correction import iterate_garnet_fe3
    from calculations.pt_fo2.params import Fe3Params, CalculationParams
    fp = Fe3Params(
        xfe3_start=PARAMS.fe3.xfe3_start,
        dxfe3=PARAMS.fe3.dxfe3,
        dt_match_c=float(dt),
        wz07_offset_c=PARAMS.fe3.wz07_offset_c,
        max_steps=PARAMS.fe3.max_steps,
        apply_empirical_correction=False,
    )
    out = iterate_garnet_fe3(
        df, verbose=False,
        params=CalculationParams(
            screening=PARAMS.screening, fe3=fp, export=PARAMS.export
        ),
    )
    return out["XFe3_gt_iter"]


def _apply_emp(est, T, slope=None, intercept=None):
    """Remove high-T residual slope; mean residual matches T < t_emp_min_c."""
    est = pd.to_numeric(est, errors="coerce")
    T = pd.to_numeric(T, errors="coerce")
    sl = PARAMS.fe3.emp_slope if slope is None else slope
    ic = PARAMS.fe3.emp_intercept if intercept is None else intercept
    out = np.where(T >= PARAMS.fe3.t_emp_min_c, est + sl * T + ic, est)
    return pd.Series(np.clip(out, 0.0, 1.0), index=est.index)


def _error_groups(df):
    """Fig. 5 blue / red groups and the remainder (calibration.error_groups)."""
    base, meas, blue, red, rem, _ = calib.error_groups(df)
    return base, meas, blue, red, rem


REF_STYLE = {
    "Kopylova et al. 1999": dict(marker="s", label="MC04 / K99"),
    "Creighton et al. 2010": dict(marker="D", label="CR10"),
    "Yaxley et al. 2017": dict(marker="^", label="YA17"),
}
RED_IDS = set(PARAMS.fe3.emp_red_ids)


def _slave_fe3_mask(df):
    return calib.calibration_mask(df), calib.measured_fe3(df)


def _norm_sid(s):
    # Sample IDs are stored un-mangled (e.g. 21-6, not Excel's 21-Jun);
    # the replacement keeps older external tables joinable.
    return (
        s.astype(str).str.strip()
        .str.replace(r"^(\d+)-Jun$", r"\1-6", regex=True)
        .str.replace(r"^(\d+)-Mar$", r"\1-3", regex=True)
        .str.replace(r"^(\d+)-Aug$", r"\1-8", regex=True)
    )


def _fe3_sigma(df):
    path = ROOT / "data" / "derived" / "garnet_Fe3_1sigma.csv"
    sig = pd.Series(index=df.index, dtype=float)
    if not path.exists():
        return sig
    tab = pd.read_csv(path)
    tab["sid"] = _norm_sid(tab["SampleID"])
    m = dict(zip(tab["sid"], pd.to_numeric(tab["garnet_Fe3+_1sigma"], errors="coerce")))
    sid = _norm_sid(df["SampleID"])
    sig = sid.map(m)
    return pd.to_numeric(sig, errors="coerce")


def _dash_hull(ax, x, y, pad=1.18):
    """Dashed convex hull around the regression-sample cloud."""
    x = np.asarray(pd.to_numeric(x, errors="coerce"), dtype=float)
    y = np.asarray(pd.to_numeric(y, errors="coerce"), dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 3:
        return
    try:
        from scipy.spatial import ConvexHull
        pts = np.column_stack([x, y])
        hull = ConvexHull(pts)
        verts = pts[hull.vertices]
        c = verts.mean(axis=0)
        verts = c + pad * (verts - c)
        verts = np.vstack([verts, verts[0]])
        ax.plot(verts[:, 0], verts[:, 1], color="0.25", ls="--", lw=0.9, zorder=2)
    except Exception:
        ax.plot(x, y, color="0.25", ls="--", lw=0.6, zorder=2)


def _plot_fe3_row(axL, axR, df, meas, T, est, ok, title,
                  blue=None, red=None, fit=None, legend=False, hull=None):
    _fe3_match_band(axL, "1to1")
    axL.plot([0, 0.35], [0, 0.35], color="0.35", lw=0.9, zorder=1)
    sigma = _fe3_sigma(df)
    has_e = ok & sigma.notna() & (sigma > 0)
    if has_e.any():
        axL.errorbar(meas[has_e], est[has_e], xerr=sigma[has_e],
                     fmt="none", ecolor="0.45", elinewidth=0.7,
                     capsize=2, zorder=2)
    plotted = pd.Series(False, index=df.index)
    for craton, sty in CRATON_STYLE.items():
        sel = ok & (df["Craton"].fillna("").astype(str) == craton)
        if not sel.any():
            continue
        plotted = plotted | sel
        axL.scatter(meas[sel], est[sel], s=36, c=sty["color"],
                    marker=sty["marker"], edgecolors="k", linewidths=0.4,
                    zorder=3, label=sty["label"] if legend else None)
    other = ok & ~plotted
    if other.any():
        axL.scatter(meas[other], est[other], s=28, facecolors="none",
                    edgecolors=OTHER["color"], marker="o", linewidths=0.6,
                    zorder=2, label="Other" if legend else None)
    if blue is not None and blue.any():
        axL.scatter(meas[blue], est[blue], s=90, facecolors="none",
                    edgecolors="#1f9ad6", linewidths=1.6, zorder=4)
    if red is not None and red.any():
        axL.scatter(meas[red], est[red], s=90, facecolors="none",
                    edgecolors="#e31a1c", linewidths=1.6, zorder=4)
    if hull is not None and hull.any():
        _dash_hull(axL, meas[hull], est[hull])
    axL.set_xlim(0, 0.35)
    axL.set_ylim(0, 0.20)
    axL.set_ylabel(r"Garnet Fe$^{3+}$/$\Sigma$Fe (estimated)")
    _style_ax(axL, r"Garnet Fe$^{3+}$/$\Sigma$Fe (measured)",
              axL.get_ylabel(), title)
    if title and "Delta" in title or (title and "T =" in title):
        axL.text(0.10, 0.97, title.split(".", 1)[-1].strip(),
                 transform=axL.transAxes, fontsize=8, va="top", ha="left")
    if legend:
        axL.legend(frameon=True, facecolor="white", edgecolor="0.55", framealpha=0.92, fontsize=6.5, loc="upper right")

    _fe3_match_band(axR, "resid")
    axR.axhline(0, color="0.35", lw=0.8, zorder=1)
    plotted = pd.Series(False, index=df.index)
    resid = meas - est
    for craton, sty in CRATON_STYLE.items():
        sel = ok & (df["Craton"].fillna("").astype(str) == craton)
        if not sel.any():
            continue
        plotted = plotted | sel
        axR.scatter(T[sel], resid[sel], s=36, c=sty["color"],
                    marker=sty["marker"], edgecolors="k", linewidths=0.4,
                    zorder=3)
    other = ok & ~plotted
    if other.any():
        axR.scatter(T[other], resid[other], s=28, facecolors="none",
                    edgecolors=OTHER["color"], marker="o", linewidths=0.6,
                    zorder=2)
    if blue is not None and blue.any():
        axR.scatter(T[blue], resid[blue], s=90, facecolors="none",
                    edgecolors="#1f9ad6", linewidths=1.6, zorder=4)
    if red is not None and red.any():
        axR.scatter(T[red], resid[red], s=90, facecolors="none",
                    edgecolors="#e31a1c", linewidths=1.6, zorder=4)
    if hull is not None and hull.any():
        _dash_hull(axR, T[hull], resid[hull])
    if fit is not None:
        slope, intercept = fit
        tt = np.linspace(1130, 1350, 50)
        axR.plot(tt, slope * tt + intercept, color="k", lw=1.0, zorder=2)
    axR.set_xlim(900, 1400)
    axR.set_ylim(-0.12, 0.20)
    axR.set_ylabel(r"Fe$^{3+}$/$\Sigma$Fe (measured−estimated)")
    _style_ax(axR, "Temperature (°C)", axR.get_ylabel())


def fig05_dt_set(df):
    """Fig. 5 — ΔT = 50, 25, 10 °C, uncorrected iterator (calibration set)."""
    if "garnet_Fe3+" not in df.columns:
        print("  skip Fig05 (no measured garnet_Fe3+)")
        return None
    base, meas = _slave_fe3_mask(df)
    T = pd.to_numeric(df["T_TA98"], errors="coerce")
    dts = (50, 25, 10)
    estimates = {dt: _fe3_at_dt(df, dt) for dt in dts}

    fit = calib.fit_empirical_correction(df)
    assert int(PARAMS.fe3.emp_fit_dt_c) in dts
    blue, red = fit["masks"]["blue"], fit["masks"]["red"]
    rem = fit["masks"]["reg"]          # dashed envelope = Eq. 5 regression set
    slope, intercept = fit["slope"], fit["intercept"]
    print(
        f"  Fig05 Eq. 5 fit (ΔT={fit['fit_dt_c']:.0f}, T≥{fit['t_min_c']:.0f} °C, "
        f"n={fit['n_regression']}): meas−est = {slope:.4e}*T + ({intercept:.5f})"
        f"   blue {fit['n_blue']}  red {fit['n_red']}"
    )
    if (abs(slope - PARAMS.fe3.emp_slope) > 5e-8
            or abs(intercept - PARAMS.fe3.emp_intercept) > 5e-5):
        print("  WARNING: params.py Eq. 5 coefficients differ from the Fig. 5 fit; "
              "run scripts/fit_calibration.py --update-params")

    fig, axes = plt.subplots(3, 2, figsize=(174 / 25.4, 210 / 25.4))
    letters = "abc"
    for i, dt in enumerate(dts):
        est = estimates[dt]
        ok = base & est.notna() & T.notna()
        _plot_fe3_row(
            axes[i, 0], axes[i, 1], df, meas, T, est, ok,
            title=f"{letters[i]}. $\\Delta$T = {dt:.0f} °C",
            blue=blue if dt == PARAMS.fe3.emp_fit_dt_c else None,
            red=red if dt == PARAMS.fe3.emp_fit_dt_c else None,
            fit=(slope, intercept) if dt == PARAMS.fe3.emp_fit_dt_c else None,
            legend=(i == 0),
            hull=rem if dt == PARAMS.fe3.emp_fit_dt_c else None,
        )
        if i < 2:
            axes[i, 0].set_xlabel("")
            axes[i, 1].set_xlabel("")
    fig.tight_layout()
    _save(fig, "Fig05_estimated_vs_measured_Fe3")
    return slope, intercept, blue, red


def _fo2_from_xfe3(df, xfe3):
    """Stagno ΔFMQ using a supplied garnet XFe3 (cation Fe split from garnet_Fe)."""
    from calculations.pt_fo2.fo2 import calculate_fo2
    tmp = df.copy()
    Fe_tot = pd.to_numeric(tmp.get("garnet_Fe"), errors="coerce")
    if "Fe3_gt_corr" in tmp.columns and "Fe2_gt_corr" in tmp.columns:
        Fe_tot = Fe_tot.fillna(
            pd.to_numeric(tmp["Fe3_gt_corr"], errors="coerce")
            + pd.to_numeric(tmp["Fe2_gt_corr"], errors="coerce")
        )
    x = pd.to_numeric(xfe3, errors="coerce")
    tmp["Fe3_gt_corr"] = x * Fe_tot
    tmp["Fe2_gt_corr"] = (1.0 - x) * Fe_tot
    return calculate_fo2(tmp, verbose=False)["delta_logfO2_FMQ"]


def fig10_fo2(df):
    """Estimated vs measured ΔQFM, Grapher-style (Slave / Siberian / Kaapvaal)."""
    if "garnet_Fe3+" not in df.columns or "delta_logfO2_FMQ" not in df.columns:
        print("  skip Fig10")
        return
    use = df.copy()
    if "screen_pass" in use.columns:
        use = use[use["screen_pass"] == True].copy()
    meas_x = pd.to_numeric(use["garnet_Fe3+"], errors="coerce")
    fo2_est = pd.to_numeric(use["delta_logfO2_FMQ"], errors="coerce")
    fo2_meas = _fo2_from_xfe3(use, meas_x)
    ok = meas_x.notna() & (meas_x > 0) & fo2_est.notna() & fo2_meas.notna()
    if not ok.any():
        print("  skip Fig10 (no measured vs estimated fO2 pairs)")
        return

    fig, ax = plt.subplots(figsize=(4.6, 4.4))
    ax.fill_between([-6, 1], [-6 - FO2_EQ_05, 1 - FO2_EQ_05],
                    [-6 + FO2_EQ_05, 1 + FO2_EQ_05],
                    color=GREEN_OUTER, lw=0, zorder=0)
    ax.fill_between([-6, 1], [-6 - FO2_EQ_02, 1 - FO2_EQ_02],
                    [-6 + FO2_EQ_02, 1 + FO2_EQ_02],
                    color=GREEN_INNER, lw=0, zorder=0)
    ax.plot([-6, 1], [-6, 1], color="0.25", lw=1.0, zorder=1)

    cr = use["Craton"].fillna("").astype(str)
    plotted = pd.Series(False, index=use.index)
    for name, sty in CRATON_STYLE.items():
        sel = ok & (cr == name)
        plotted = plotted | sel
        if sel.any():
            ax.scatter(fo2_meas[sel], fo2_est[sel], marker=sty["marker"],
                       c=sty["color"], s=28, edgecolors="k", linewidths=0.4,
                       zorder=3, label=sty["label"])
    other = ok & ~plotted
    if other.any():
        ax.scatter(fo2_meas[other], fo2_est[other], marker=OTHER["marker"],
                   facecolors="none", edgecolors=OTHER["color"], s=22,
                   linewidths=0.6, zorder=2, label=OTHER["label"])

    _, _, blue, red, _ = _error_groups(use)
    if (ok & blue).any():
        ax.scatter(fo2_meas[ok & blue], fo2_est[ok & blue], s=90,
                   facecolors="none", edgecolors="#1f9ad6", linewidths=1.6, zorder=4)
    if (ok & red).any():
        ax.scatter(fo2_meas[ok & red], fo2_est[ok & red], s=90,
                   facecolors="none", edgecolors="#e31a1c", linewidths=1.6, zorder=4)

    ax.set_xlim(-4.5, -0.5)
    ax.set_ylim(-4.2, -0.6)
    ax.legend(frameon=True, facecolor="white", edgecolor="0.55", framealpha=0.92, fontsize=7, loc="lower right")
    _style_ax(ax,
              r"log $f$O$_2$ ($\Delta$QFM; measured)",
              r"log $f$O$_2$ ($\Delta$QFM; estimated)")
    n = int(ok.sum())
    d = (fo2_est[ok] - fo2_meas[ok]).abs()
    print(f"  Fig10 n={n}  ±{FO2_EQ_02:.2f}: {(d<=FO2_EQ_02).mean()*100:.0f}%  "
          f"±{FO2_EQ_05:.2f}: {(d<=FO2_EQ_05).mean()*100:.0f}%")
    _save(fig, "Fig10_fO2_estimated")


def fig11_p_fo2_measured(df):
    """Fig. 12 layout, ΔFMQ from measured garnet Fe³⁺ only."""
    if "garnet_Fe3+" not in df.columns:
        print("  skip Fig11")
        return
    from calculations.pt_fo2.buffers import (
        buffer_table, nimis2015_band, verify_against_literature,
    )
    meas = pd.to_numeric(df["garnet_Fe3+"], errors="coerce")
    use = df.loc[meas.notna() & (meas > 0) & (meas < 0.5)].copy()
    if use.empty:
        print("  skip Fig11 (no measured Fe3+)")
        return
    fo2 = _fo2_from_xfe3(use, use["garnet_Fe3+"])
    use = use.assign(fo2_meas=fo2)
    use = use.loc[use["fo2_meas"].notna()]

    fig, ax = plt.subplots(figsize=(174 / 25.4, 120 / 25.4))
    P = np.linspace(30, 70, 81)
    buf = buffer_table(P, df=use)
    lo, hi = nimis2015_band(P)
    ax.fill_betweenx(P, lo, hi, color="0.85", lw=0, zorder=0)
    ax.plot(buf["dEMOD"], P, color="#d4a017", lw=1.15, ls="--", zorder=1)
    ax.plot(buf["dFeNi"], P, color="#4c72b0", lw=1.15, ls="-.", zorder=1)
    ax.plot(buf["dBSE"], P, color="k", lw=1.0, zorder=1)
    ax.axvline(0, color="0.4", lw=0.6, zorder=1)
    _label_on_curve(ax, buf["dEMOD"], P, "EMOG/D", "#d4a017", 42)
    _label_on_curve(ax, buf["dFeNi"], P, "Fe–Ni", "#4c72b0", 58)
    _label_on_curve(ax, buf["dBSE"], P, "BSE", "k", 48)

    sig = _fe3_sigma(use)
    has_s = sig.notna() & (sig > 0)
    if has_s.any():
        fo2_lo = _fo2_from_xfe3(use.loc[has_s], use.loc[has_s, "garnet_Fe3+"] - sig[has_s])
        fo2_hi = _fo2_from_xfe3(use.loc[has_s], use.loc[has_s, "garnet_Fe3+"] + sig[has_s])
        x = use.loc[has_s, "fo2_meas"]
        xerr = np.vstack([
            (x - fo2_lo).abs().to_numpy(),
            (fo2_hi - x).abs().to_numpy(),
        ])
        ax.errorbar(
            x, use.loc[has_s, "P_CA91"], xerr=xerr, fmt="none",
            ecolor="0.45", elinewidth=0.6, capsize=1.5, zorder=2,
        )

    _scatter_by_craton(ax, use, "fo2_meas", "P_CA91", legend_loc="upper left")

    # X only if a mineral-equilibrium test fails. T < 900 °C alone is not an X.
    sp = PARAMS.screening
    d_px = pd.to_numeric(use.get("dT_opx_cpx"), errors="coerce")
    d_gt = pd.to_numeric(use.get("dT_pyx_gt"), errors="coerce")
    if not isinstance(d_gt, pd.Series) or d_gt.isna().all():
        d_gt = pd.to_numeric(use.get("dT_KR00"), errors="coerce")
    d_ol = pd.to_numeric(use.get("dT_gt_ol"), errors="coerce")
    eq_ok = (
        d_px.notna() & (d_px < sp.dt_opx_cpx_c)
        & d_gt.notna() & (d_gt < sp.dt_pyx_gt_c)
        & d_ol.notna() & (d_ol < sp.dt_gt_ol_c)
    )
    fail = ~eq_ok
    if fail.any():
        ax.scatter(use.loc[fail, "fo2_meas"], use.loc[fail, "P_CA91"],
                   marker="o", s=46, c="0.62", edgecolors="0.25",
                   linewidths=0.5, zorder=4, label="_nolegend_")
        ax.scatter(use.loc[fail, "fo2_meas"], use.loc[fail, "P_CA91"],
                   marker="x", s=28, c="w", linewidths=1.15, zorder=5,
                   label="_nolegend_")

    _, _, blue, red, _ = _error_groups(df)
    blue_u = use.index.isin(df.index[blue])
    red_u = use.index.isin(df.index[red])
    if blue_u.any():
        ax.scatter(use.loc[blue_u, "fo2_meas"], use.loc[blue_u, "P_CA91"],
                   s=140, facecolors="none", edgecolors="#1f9ad6",
                   linewidths=1.6, zorder=6, label="Blue error group")
    if red_u.any():
        ax.scatter(use.loc[red_u, "fo2_meas"], use.loc[red_u, "P_CA91"],
                   s=140, facecolors="none", edgecolors="#e31a1c",
                   linewidths=1.6, zorder=6, label="Red error group")
    p_ma, f_ma = _fo2_depth_mean(use["P_CA91"], use["fo2_meas"])
    if len(p_ma):
        ax.plot(f_ma, p_ma, color="#333333", lw=1.8, zorder=6, label="Depth mean")
        _label_on_curve(ax, f_ma, p_ma, "mean", "#333333", 55)
    handles, labels = ax.get_legend_handles_labels()
    if fail.any():
        from matplotlib.legend_handler import HandlerBase
        from matplotlib.lines import Line2D

        class HandlerXOnCircle(HandlerBase):
            def create_artists(self, legend, orig_handle, xdescent, ydescent,
                               width, height, fontsize, trans):
                cx = 0.5 * width - xdescent
                cy = 0.5 * height - ydescent
                circ = Line2D([cx], [cy], marker="o", ls="none",
                              markerfacecolor="0.62", markeredgecolor="0.25",
                              markeredgewidth=0.6, markersize=8)
                cross = Line2D([cx], [cy], marker="x", ls="none", color="w",
                               markeredgecolor="w", markeredgewidth=1.3, markersize=5.5)
                circ.set_transform(trans)
                cross.set_transform(trans)
                return [circ, cross]

        marker = Line2D([0], [0], ls="none")
        handles.insert(0, marker)
        labels.insert(0, "Equilibrium fail")
        handler_map = {marker: HandlerXOnCircle()}
    else:
        handler_map = {}
    ax.legend(handles, labels, handler_map=handler_map,
              frameon=True, facecolor="white", edgecolor="0.55", framealpha=0.95,
              fontsize=6.5, loc="upper left", labelspacing=0.7,
              handletextpad=0.5, borderpad=0.4)
    ax.text(0.98, 0.03, f"n = {len(use)}", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=8)

    ax.invert_yaxis()
    _style_ax(ax, r"log $f$O$_2$ ($\Delta$FMQ)", "P (kbar)")
    ax.set_xlim(-6, 0)
    ax.set_ylim(70, 30)
    _gpa_axis(ax)
    _add_depth_axis(ax, km_interval=50)
    print(f"  Fig11 measured Fe3+ fO2 n={len(use)}  "
          f"blue={int(blue_u.sum())} red={int(red_u.sum())}  unscreened included")
    _save(fig, "Fig09_P_fO2_measured")


def _fo2_depth_mean(P, fo2, half_km=15, min_n=5):
    """Moving mean of ΔFMQ in a ±half_km window, returned as (P, mean)."""
    P = pd.to_numeric(P, errors="coerce")
    fo2 = pd.to_numeric(fo2, errors="coerce")
    m = P.notna() & fo2.notna()
    P, fo2 = P[m], fo2[m]
    if len(P) < min_n:
        return np.array([]), np.array([])
    # approximate 1 kbar ≈ 3.2 km for the window
    half_p = half_km / 3.2
    centres = np.arange(32, 69, 2.0)
    out_p, out_f = [], []
    for c in centres:
        sel = (P - c).abs() <= half_p
        if int(sel.sum()) >= min_n:
            out_p.append(c)
            out_f.append(float(fo2[sel].mean()))
    return np.array(out_p), np.array(out_f)


def fig12_p_fo2(df, mark_cpx_opx=False, stem="Fig10_P_fO2_estimated"):
    if "delta_logfO2_FMQ" not in df.columns:
        print("  skip Fig12")
        return
    from calculations.pt_fo2.buffers import (
        buffer_table, nimis2015_band, verify_against_literature,
    )
    use = df[df["screen_pass"] == True].copy() if "screen_pass" in df.columns else df.copy()
    fig, ax = plt.subplots(figsize=(174 / 25.4, 190 / 25.4))

    P = np.linspace(20, 80, 121)
    buf = buffer_table(P, df=use)
    lo, hi = nimis2015_band(P)
    ax.fill_betweenx(P, lo, hi, color="0.85", lw=0, zorder=0)
    ax.plot(buf["dEMOD"], P, color="#d4a017", lw=1.3, ls="--", zorder=1)
    ax.plot(buf["dFeNi"], P, color="#4c72b0", lw=1.3, ls="-.", zorder=1)
    ax.plot(buf["dBSE"], P, color="k", lw=1.1, zorder=1)
    ax.axvline(0, color="0.4", lw=0.6, zorder=1)
    _label_on_curve(ax, buf["dEMOD"], P, "EMOG/D", "#d4a017", 44)
    _label_on_curve(ax, buf["dFeNi"], P, "Fe–Ni", "#4c72b0", 52)
    _label_on_curve(ax, buf["dBSE"], P, "BSE", "k", 46)

    _scatter_by_craton(ax, use, "delta_logfO2_FMQ", "P_CA91", legend_loc="upper left")

    meas = pd.to_numeric(use["garnet_Fe3+"], errors="coerce")
    slave = use["Craton"].fillna("").astype(str).str.contains("Slave")
    sib = use["Craton"].fillna("").astype(str).str.contains("Siberian")
    kaa = use["Craton"].fillna("").astype(str).str.contains("Kaapvaal")
    has = meas.notna() & (meas > 0) & (meas < 0.5) & use["delta_logfO2_FMQ"].notna()
    cal = has & slave
    val = has & (sib | kaa)
    if cal.any():
        ax.scatter(use.loc[cal, "delta_logfO2_FMQ"], use.loc[cal, "P_CA91"],
                   marker="*", s=90, c="k", edgecolors="k", linewidths=0.3,
                   zorder=5, label="Calibration")
    if val.any():
        ax.scatter(use.loc[val, "delta_logfO2_FMQ"], use.loc[val, "P_CA91"],
                   marker="*", s=90, facecolors="none", edgecolors="k",
                   linewidths=0.8, zorder=5, label="Validation")
    # Fe3+/ΣFe = 0.01 isopleth: median chemistry of samples that sit on
    # the 0.01 floor, T from the xenolith geotherm. Draw only the part
    # inside the axes and label that segment.
    try:
        from calculations.pt_fo2.buffers import xenolith_geotherm_C
        p_iso = np.linspace(30, 70, 81)
        xcol = "XFe3_gt_corr" if "XFe3_gt_corr" in use.columns else "XGtFe3+_corr"
        floor = use[pd.to_numeric(use.get(xcol), errors="coerce") <= 0.012]
        base = floor if len(floor) >= 8 else use
        chem = base.median(numeric_only=True)
        rows = pd.DataFrame([chem] * len(p_iso))
        rows["P_CA91"] = p_iso
        rows["T_TA98"] = xenolith_geotherm_C(p_iso, use)
        iso = _fo2_from_xfe3(rows, np.full(len(p_iso), 0.01)).to_numpy(dtype=float)
        vis = np.isfinite(iso) & (iso >= -6.85) & (iso <= -0.15)
        if vis.any():
            ax.plot(iso[vis], p_iso[vis], color="#8c5a2b", lw=1.15, ls="--", zorder=2)
            ax.text(-4.15, 31.2, r"Fe$^{3+}$/$\Sigma$Fe = 0.01",
                    color="#8c5a2b", fontsize=8, va="center", ha="left",
                    clip_on=True, zorder=8)
            print(f"  isopleth visible P {p_iso[vis][0]:.0f}–{p_iso[vis][-1]:.0f} "
                  f"ΔFMQ {iso[vis][0]:.2f} to {iso[vis][-1]:.2f}")
    except Exception as exc:
        print("  isopleth skipped", exc)
    p_ma, f_ma = _fo2_depth_mean(use["P_CA91"], use["delta_logfO2_FMQ"])
    if len(p_ma):
        ax.plot(f_ma, p_ma, color="#333333", lw=1.8, zorder=6, label="Depth mean")
        _label_on_curve(ax, f_ma, p_ma, "mean", "#333333", 40)
    # visual lower-left: pressure axis is inverted, so high P is the bottom
    # Combined residual 1σ, with and without the calibration error groups.
    # Same 1σ values as Table 1 / text (calibration.calibration_stats),
    # rounded to 3 dp as reported, converted with FO2_EQ_02 per 0.02.
    st = calib.calibration_stats(df, fo2_per_002=FO2_EQ_02).set_index("dataset")
    s_in, w_in = st.loc["Calibration + validation", ["sigma", "sigma_log"]]
    s_ex, w_ex = st.loc["Combined, error groups excluded", ["sigma", "sigma_log"]]
    print(f"  Fig10 bars: 1σ {s_in:.3f} → {w_in:.2f} log; {s_ex:.3f} → {w_ex:.2f} log")
    x0, p_in, p_ex = -6.85, 68.4, 65.6
    ax.plot([x0, x0 + w_in], [p_in, p_in], color="#6b4c9a", lw=3.0, solid_capstyle="butt", zorder=7)
    ax.plot([x0, x0 + w_ex], [p_ex, p_ex], color="#1b7f3a", lw=3.0, solid_capstyle="butt", zorder=7)
    ax.text(x0 + w_in + 0.08, p_in, f"{w_in:.2f} log",
            ha="left", va="center", fontsize=8, color="#6b4c9a", clip_on=True, zorder=8)
    ax.text(x0 + w_ex + 0.08, p_ex, f"{w_ex:.2f} log",
            ha="left", va="center", fontsize=8, color="#1b7f3a", clip_on=True, zorder=8)
    if mark_cpx_opx:
        na = pd.to_numeric(use["cpx_Na"], errors="coerce")
        ti = pd.to_numeric(use["cpx_TiO2"], errors="coerce")
        al = pd.to_numeric(use["opx_Al2O3"], errors="coerce")
        # Enriched cpx (Na and Ti above the screened median) with depleted
        # opx (Al2O3 at or below the five-sample range, 0.9 wt%).
        pat = (na >= 0.10) & (ti >= 0.20) & (al <= 0.90) & use["delta_logfO2_FMQ"].notna()
        ax.scatter(use.loc[pat, "delta_logfO2_FMQ"], use.loc[pat, "P_CA91"],
                   marker="x", s=46, c="#c0392b", linewidths=1.8,
                   zorder=9, label="Enriched cpx, depleted opx")
        print(f"  pattern marks n={int(pat.sum())} of {len(use)}")
    ax.legend(frameon=True, facecolor="white", edgecolor="0.55", framealpha=0.95,
              fontsize=6.5, loc="upper left", labelspacing=0.45, borderpad=0.35)
    plotted = (
        use["delta_logfO2_FMQ"].between(-7, 0)
        & use["P_CA91"].between(30, 70)
    )
    ax.text(0.98, 0.03, f"n = {int(plotted.sum())}", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=8)
    print(f"  Fig12 stars: calibration n={int(cal.sum())}  validation n={int(val.sum())}")
    ax.invert_yaxis()
    _style_ax(ax, r"log $f$O$_2$ ($\Delta$FMQ)", "P (kbar)")
    ax.set_ylim(70, 30)
    _gpa_axis(ax)
    _add_depth_axis(ax, km_interval=50)
    ax.set_xlim(-7, 0)
    _save(fig, stem)
    verify_against_literature(verbose=True)



def _dt_sweep_corrected(df, force=True):
    """Fig. 6 sweep: locked Eq. 5 held fixed, Eq. 5 regression samples only."""
    path = ROOT / "data" / "derived" / "dt_match_sweep_corr.csv"
    if path.exists() and not force:
        return pd.read_csv(path)
    print("  Recomputing ΔT sweep (Eq. 5 fixed; regression samples only)…")
    reg = calib.error_groups(df)[5]
    s = calib.dt_sweep(df, PARAMS.fe3.emp_slope, PARAMS.fe3.emp_intercept, mask=reg)
    path.parent.mkdir(parents=True, exist_ok=True)
    s.to_csv(path, index=False)
    print(f"  wrote {path}")
    return s


def fig06_dt_sweep(df):
    """Fig. 6 — choose ΔT on the Eq. 5 regression samples, Eq. 5 fixed."""
    s = _dt_sweep_corrected(df, force=True)
    # Highest ±0.02 rate; ties broken by smallest |mean|, then ±0.05
    opt = calib.select_dt(s)
    print(f"  Fig06 selected ΔT (Eq. 5 regression samples, n={int(s['n'].max())}): {opt} °C")
    print(s.loc[s["dt"] == opt].to_string(index=False))
    if opt != PARAMS.fe3.dt_match_c:
        print(f"  WARNING: params.py dt_match_c = {PARAMS.fe3.dt_match_c} differs; "
              "run scripts/fit_calibration.py --update-params")

    fig, ax1 = plt.subplots(figsize=(119 / 25.4, 78 / 25.4))
    ax1.plot(s["dt"], s["mean"], color="k", lw=1.15, label=r"$\mu$ (est − meas)")
    ax1.axhline(0, color="0.45", lw=0.7)
    ax1.axvline(opt, color="0.4", ls="--", lw=0.8)
    ax1.set_xlabel(r"$\Delta$T match threshold (°C)")
    ax1.set_ylabel(r"mean Fe$^{3+}$/$\Sigma$Fe (est − meas)")
    ax1.set_ylim(-0.05, 0.05)
    ax2 = ax1.twinx()
    ax2.plot(s["dt"], 100 * s["p02"], color="#5e3c99", lw=1.15, label="% within ±0.02")
    ax2.plot(s["dt"], 100 * s["p05"], color="#b35806", lw=1.15, ls="--",
             label="% within ±0.05")
    ax2.set_ylabel("Success rate (%)")
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    # Drawn on ax2 (the top axes) so no curve is drawn over it; lower left
    # is clear of all three curves.
    ax2.legend(h1 + h2, l1 + l2, frameon=True, facecolor="white", edgecolor="0.55",
               framealpha=0.92, fontsize=7, loc="lower left")
    ax1.tick_params(direction="in")
    ax2.tick_params(direction="in")
    _save(fig, "Fig06_dt_match_calibration")
    return opt


def fig07_optimal_dt(df, opt_dt):
    """Fig. 7 — same layout as Fig. 5, only the optimal ΔT, after correction."""
    if "garnet_Fe3+" not in df.columns:
        print("  skip Fig07")
        return
    base, meas, blue, red, rem = _error_groups(df)
    T = pd.to_numeric(df["T_TA98"], errors="coerce")
    raw = _fe3_at_dt(df, opt_dt)
    sl, ic = PARAMS.fe3.emp_slope, PARAMS.fe3.emp_intercept   # locked Eq. 5
    est = _apply_emp(raw, T, slope=sl, intercept=ic)
    ok = base & est.notna() & T.notna()
    fig, axes = plt.subplots(1, 2, figsize=(174 / 25.4, 72 / 25.4))
    _plot_fe3_row(
        axes[0], axes[1], df, meas, T, est, ok,
        title="a.",
        blue=blue, red=red, fit=None, legend=False,
    )
    _panel_label(axes[1], "b")
    diff_all = (est[ok] - meas[ok]).abs()
    diff_rem = (est[rem & ok] - meas[rem & ok]).abs()
    axes[0].text(
        0.98, 0.04,
        f"all n = {int(ok.sum())}\n"
        f"±0.02: {(diff_all<=0.02).mean()*100:.0f}%\n"
        f"±0.05: {(diff_all<=0.05).mean()*100:.0f}%\n"
        f"remainder n = {int((rem&ok).sum())}\n"
        f"±0.02: {(diff_rem<=0.02).mean()*100:.0f}%\n"
        f"±0.05: {(diff_rem<=0.05).mean()*100:.0f}%",
        transform=axes[0].transAxes, fontsize=7.5, ha="right", va="bottom",
    )
    print(f"  Fig07 ΔT={opt_dt:g}  Eq. 5 = {sl:.4e}*T + ({ic:.5f})  "
          f"remainder ±0.02={(diff_rem<=0.02).mean()*100:.0f}%")
    fig.tight_layout()
    _save(fig, "Fig07_corrected_Fe3")


def fig08_stoichiometry(df):
    """Stoichiometric Fe3+ vs measured: Droop (1987) and Quinn et al. (2016)."""
    if "garnet_Fe3+" not in df.columns:
        print("  skip Fig08")
        return
    meas = pd.to_numeric(df["garnet_Fe3+"], errors="coerce")
    slave = df["Craton"].fillna("").astype(str).str.contains("Slave")
    refs = df["Reference"].fillna("").astype(str)
    sigma = _fe3_sigma(df)
    est = pd.to_numeric(df["XFe3_quinn"], errors="coerce")
    ok = slave & meas.notna() & (meas > 0) & (meas < 0.5) & est.notna()
    fig, ax = plt.subplots(figsize=(4.4, 4.2))
    ax.fill_between([0, 0.40], [-0.05, 0.35], [0.05, 0.45],
                    color=GREEN_OUTER, lw=0, zorder=0)
    ax.fill_between([0, 0.40], [-0.02, 0.38], [0.02, 0.42],
                    color=GREEN_INNER, lw=0, zorder=0)
    ax.plot([0, 0.40], [0, 0.40], color="0.35", lw=0.9, zorder=1)
    has_e = ok & sigma.notna() & (sigma > 0)
    if has_e.any():
        ax.errorbar(meas[has_e], est[has_e], xerr=sigma[has_e],
                    fmt="none", ecolor="0.45", elinewidth=0.7,
                    capsize=2, zorder=2)
    plotted = pd.Series(False, index=df.index)
    for ref, sty in REF_STYLE.items():
        sel = ok & refs.str.contains(ref.split()[0])
        plotted = plotted | sel
        ax.scatter(meas[sel], est[sel], s=36, c="k",
                   marker=sty["marker"], zorder=3)
    other = ok & ~plotted
    if other.any():
        ax.scatter(meas[other], est[other], s=28, c="k", marker="o", zorder=3)
    _, _, blue, red, _ = _error_groups(df)
    if (ok & blue).any():
        ax.scatter(meas[ok & blue], est[ok & blue], s=90, facecolors="none",
                   edgecolors="#1f9ad6", linewidths=1.6, zorder=4)
    if (ok & red).any():
        ax.scatter(meas[ok & red], est[ok & red], s=90, facecolors="none",
                   edgecolors="#e31a1c", linewidths=1.6, zorder=4)
    ax.set_xlim(0, 0.40)
    ax.set_ylim(0, 0.40)
    _style_ax(ax,
              r"Garnet $X$Fe$^{3+}$ (measured)",
              r"Garnet $X$Fe$^{3+}$ (QU16)",
              r"Stoichiometric Fe$^{3+}$ estimation")
    diff = (est[ok] - meas[ok]).abs()
    ax.text(
        0.97, 0.04,
        f"n = {int(ok.sum())}\n±0.02: {(diff<=0.02).mean()*100:.0f}%\n"
        f"±0.05: {(diff<=0.05).mean()*100:.0f}%",
        transform=ax.transAxes, fontsize=8, ha="right", va="bottom",
    )
    print(f"  Fig08 Quinn: n={int(ok.sum())}  "
          f"±0.02={(diff<=0.02).mean()*100:.0f}%  "
          f"±0.05={(diff<=0.05).mean()*100:.0f}%")
    _save(fig, "Fig08_stoichiometric_Fe3")


def fig09_nonslave_fe3(df):
    """Fig. 9 — measured vs estimated Fe3+ outside the Slave craton.

    Same craton symbols/colours as Fig. 12. Estimated Fe3+ is the
    pipeline value after the high-T empirical correction.
    """
    if "garnet_Fe3+" not in df.columns:
        print("  skip Fig09")
        return
    meas = pd.to_numeric(df["garnet_Fe3+"], errors="coerce")
    if "XFe3_gt_corr" in df.columns:
        est = pd.to_numeric(df["XFe3_gt_corr"], errors="coerce")
    else:
        T = pd.to_numeric(df["T_TA98"], errors="coerce")
        est = _apply_emp(_fe3_at_dt(df, PARAMS.fe3.dt_match_c), T)
    T = pd.to_numeric(df["T_TA98"], errors="coerce")
    slave = df["Craton"].fillna("").astype(str).str.contains("Slave")
    passed = df["screen_pass"] == True if "screen_pass" in df.columns else True
    ok = (~slave) & passed & meas.notna() & (meas > 0) & (meas < 0.5) & est.notna()
    if not ok.any():
        print("  skip Fig09 (no non-Slave measured Fe3+)")
        return

    fig, axes = plt.subplots(1, 2, figsize=(174 / 25.4, 72 / 25.4))
    axL, axR = axes
    resid = meas - est
    x_max = float(max(meas[ok].max(), est[ok].max()) + 0.02)
    t_min = float(T[ok].min() - 40)
    t_max = float(T[ok].max() + 40)
    r_min = float(resid[ok].min() - 0.02)
    r_max = float(resid[ok].max() + 0.02)

    axL.fill_between([0, x_max], [-0.05, x_max - 0.05], [0.05, x_max + 0.05],
                     color=GREEN_OUTER, lw=0, zorder=0)
    axL.fill_between([0, x_max], [-0.02, x_max - 0.02], [0.02, x_max + 0.02],
                     color=GREEN_INNER, lw=0, zorder=0)
    axL.plot([0, x_max], [0, x_max], color="0.35", lw=0.9, zorder=1)
    _fe3_match_band(axR, "resid")
    axR.axhline(0, color="0.35", lw=0.8, zorder=1)

    plotted = pd.Series(False, index=df.index)
    for craton, sty in CRATON_STYLE.items():
        if "Slave" in craton:
            continue
        sel = ok & (df["Craton"].fillna("").astype(str) == craton)
        if not sel.any():
            continue
        plotted = plotted | sel
        axL.scatter(meas[sel], est[sel], s=36, c=sty["color"],
                    marker=sty["marker"], edgecolors="k", linewidths=0.4,
                    label=sty["label"], zorder=3)
        axR.scatter(T[sel], resid[sel], s=36, c=sty["color"],
                    marker=sty["marker"], edgecolors="k", linewidths=0.4,
                    zorder=3)
    other = ok & ~plotted
    if other.any():
        axL.scatter(meas[other], est[other], s=28, facecolors="none",
                    edgecolors=OTHER["color"], marker="o", linewidths=0.7,
                    label=OTHER["label"], zorder=2)
        axR.scatter(T[other], resid[other], s=28, facecolors="none",
                    edgecolors=OTHER["color"], marker="o", linewidths=0.7,
                    zorder=2)

    axL.set_xlim(0, x_max)
    axL.set_ylim(0, x_max)
    axL.legend(frameon=True, facecolor="white", edgecolor="0.55", framealpha=0.92,
               fontsize=7, loc="upper right")
    _style_ax(axL,
              r"Garnet Fe$^{3+}$/$\Sigma$Fe (measured)",
              r"Garnet Fe$^{3+}$/$\Sigma$Fe (estimated)",
              "a.")
    _panel_label(axR, "b")
    diff = (est[ok] - meas[ok]).abs()
    axL.text(
        0.98, 0.04,
        f"n = {int(ok.sum())}\n±0.02: {(diff<=0.02).mean()*100:.0f}%\n"
        f"±0.05: {(diff<=0.05).mean()*100:.0f}%",
        transform=axL.transAxes, fontsize=8, ha="right", va="bottom",
    )

    axR.set_xlim(t_min, t_max)
    axR.set_ylim(r_min, r_max)
    _style_ax(axR, "Temperature (°C)",
              r"Fe$^{3+}$/$\Sigma$Fe (measured−estimated)")
    print(f"  Fig09 non-Slave: n={int(ok.sum())}  "
          f"±0.02={(diff<=0.02).mean()*100:.0f}%  "
          f"±0.05={(diff<=0.05).mean()*100:.0f}%")
    fig.tight_layout()
    _save(fig, "Fig08_nonslave_Fe3")


def _geo_rings(geom):
    """Yield exterior rings (Nx2) from a GeoJSON geometry."""
    if not geom:
        return
    gtype = geom.get("type")
    coords = geom.get("coordinates") or []
    if gtype == "Polygon":
        polys = [coords]
    elif gtype == "MultiPolygon":
        polys = coords
    else:
        return
    for poly in polys:
        if not poly:
            continue
        ring = np.asarray(poly[0], dtype=float)
        if ring.ndim == 2 and len(ring) >= 3:
            yield ring[:, :2]


def _add_geojson(ax, path, facecolor, edgecolor="none", lw=0.0, zorder=0, alpha=1.0):
    import json
    data = json.loads(Path(path).read_text())
    for feat in data.get("features", []):
        for ring in _geo_rings(feat.get("geometry")):
            ax.add_patch(plt.Polygon(
                ring, closed=True, facecolor=facecolor, edgecolor=edgecolor,
                linewidth=lw, zorder=zorder, alpha=alpha, joinstyle="round",
            ))


def _label_xy(feat):
    rings = list(_geo_rings(feat.get("geometry")))
    if not rings:
        return None
    ring = max(rings, key=len)
    return float(np.mean(ring[:, 0])), float(np.mean(ring[:, 1]))


def fig01_map():
    """Wide geologic base map (top third of a page) in the Hygrohodos palette."""
    web = ROOT / "data" / "spatial" / "web"
    broad = web / "geology_broad.geojson"
    cratons = web / "geology_cratons.geojson"
    pipes = web / "kimberlite_localities.geojson"
    if not broad.is_file():
        print("  skip Fig01 (no geology_broad.geojson)")
        return

    colors = {
        "oceanic": "#d5e4ee",
        "other": "#cfc8bc",
        "orogen": "#9aa0c4",
        "sedimentary_basin": "#8fb9b0",
        "volcanic_province": "#c47a62",
        "craton": "#c4a574",
    }
    labels = {
        "oceanic": "Oceanic crust",
        "sedimentary_basin": "Sedimentary basin",
        "volcanic_province": "Volcanic / magmatic",
        "orogen": "Orogen",
        "craton": "Craton / shield",
    }

    fig, ax = plt.subplots(figsize=(174 / 25.4, 82 / 25.4))
    ax.set_facecolor(colors["oceanic"])

    import json
    broad_data = json.loads(broad.read_text())
    draw_order = ["oceanic", "other", "orogen", "sedimentary_basin",
                  "volcanic_province", "craton"]
    by_class = {f["properties"].get("broad"): f for f in broad_data["features"]}
    for key in draw_order:
        feat = by_class.get(key)
        if feat is None:
            continue
        for ring in _geo_rings(feat.get("geometry")):
            ax.add_patch(plt.Polygon(
                ring, closed=True, facecolor=colors[key], edgecolor="none",
                linewidth=0, zorder=1, joinstyle="round",
            ))

    if cratons.is_file():
        craton_data = json.loads(cratons.read_text())
        for feat in craton_data["features"]:
            for ring in _geo_rings(feat.get("geometry")):
                ax.add_patch(plt.Polygon(
                    ring, closed=True, facecolor="none", edgecolor="#6b4f2a",
                    linewidth=0.25, zorder=2, alpha=0.7,
                ))
        # Named craton labels only
        wanted = {
            "Slave Craton": (8, 6),
            "Siberian Craton": (4, -8),
            "Kaapvaal Craton": (10, -8),
        }
        placed = set()
        # Prefer the named polygon; fall back to group centroid
        for feat in craton_data["features"]:
            name = feat["properties"].get("name") or ""
            group = feat["properties"].get("group") or ""
            key = name if name in wanted else (group if group in wanted else None)
            if key is None or key in placed:
                continue
            if key == "Siberian Craton" and name != "Anabar Shield":
                continue
            if key == "Slave Craton" and name != "Slave Craton":
                continue
            if key == "Kaapvaal Craton" and name != "Kaapvaal Craton":
                continue
            xy = _label_xy(feat)
            if xy is None:
                continue
            dx, dy = wanted[key]
            ax.text(xy[0] + dx, xy[1] + dy, key.replace(" Craton", ""),
                    fontsize=8, fontweight="bold", color="#3b2a14",
                    ha="left", va="center", zorder=5,
                    bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                              edgecolor="none", alpha=0.7))
            placed.add(key)

    if pipes.is_file():
        pipe_data = json.loads(pipes.read_text())
        by_c = {k: {"lon": [], "lat": []} for k in list(CRATON_STYLE) + ["Other"]}
        for feat in pipe_data["features"]:
            geom = feat.get("geometry") or {}
            if geom.get("type") != "Point":
                continue
            lon, lat = geom["coordinates"][:2]
            cr = feat.get("properties", {}).get("craton") or "Other"
            bucket = cr if cr in CRATON_STYLE else "Other"
            by_c[bucket]["lon"].append(lon)
            by_c[bucket]["lat"].append(lat)
        for craton, sty in CRATON_STYLE.items():
            pts = by_c[craton]
            if not pts["lon"]:
                continue
            ax.scatter(pts["lon"], pts["lat"], s=14, marker=sty["marker"],
                       c=sty["color"], edgecolors="k", linewidths=0.25,
                       zorder=4, label=sty["label"])
        if by_c["Other"]["lon"]:
            ax.scatter(by_c["Other"]["lon"], by_c["Other"]["lat"], s=10,
                       marker="o", facecolors="none", edgecolors="#1f1a17",
                       linewidths=0.4, zorder=3, label="Other")

    from matplotlib.patches import Patch
    handles = [Patch(facecolor=colors[k], edgecolor="0.4", linewidth=0.3,
                     label=labels[k])
               for k in ("craton", "sedimentary_basin", "orogen",
                         "volcanic_province", "oceanic")]
    geo_leg = ax.legend(handles=handles, loc="lower left", fontsize=6.5,
                        frameon=True, facecolor="white", edgecolor="0.55",
                        framealpha=0.9, title="Geology", title_fontsize=7)
    ax.add_artist(geo_leg)
    ax.legend(loc="lower right", fontsize=6.5, frameon=True, facecolor="white",
              edgecolor="0.55", framealpha=0.9, title="Kimberlites",
              title_fontsize=7, scatterpoints=1)

    ax.set_xlim(-170, 190)
    ax.set_ylim(-58, 78)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.tick_params(direction="in", top=True, right=True)
    _save(fig, "Fig01_locality_map")


def main():
    if not FULL.exists():
        print(f"Missing {FULL} — run: python scripts/run_pipeline.py")
        sys.exit(1)
    df = pd.read_csv(FULL, encoding="utf-8-sig", low_memory=False)
    df = df.rename(columns={df.columns[0]: "SampleID"})
    print(f"Loaded {len(df)} rows from {FULL.name}")
    tab = calib.calibration_stats(df, fo2_per_002=FO2_EQ_02)
    tab.to_csv(ROOT / "data" / "derived" / "Table1_calibration_stats.csv", index=False)
    print(tab.drop(columns="notes").round(4).to_string(index=False))
    e02, e05 = fo2_equivalents(df)
    print(f"  fO2 equivalents used: 0.02 → {FO2_EQ_02:.2f} log, 0.05 → {FO2_EQ_05:.2f} log"
          f"  (median on current screened samples: {e02:.3f}, {e05:.3f})")
    print(f"Writing figures to {OUT}")
    _apply_journal_style()
    fig01_map()
    fig02_kd_fe3(df)
    fig04_screens(df)
    try:
        from plot_exp_ta98_wz07 import main as plot_fig03_experiments
        plot_fig03_experiments()
    except Exception as exc:
        print(f"  Fig03 experiment plot skipped: {exc}")
    fig05_dt_set(df)
    fig06_dt_sweep(df)
    fig07_optimal_dt(df, PARAMS.fe3.dt_match_c)
    fig09_nonslave_fe3(df)
    fig11_p_fo2_measured(df)
    fig12_p_fo2(df)
    print("Done.")


if __name__ == "__main__":
    main()
