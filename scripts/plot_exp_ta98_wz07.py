"""TA98 vs WZ07 on experiments that report measured garnet Fe3+/ΣFe.

Sources
-------
* Matjuschkin et al. 2014 — full mineral analyses from
  data/incoming/melt_experiments.xlsx (the 2026-09-13 compilation is
  missing MgO/CaO for this study).
* Rzehak et al. 2020 — data/incoming/experiments_20260913.csv.

WZ07 is the olivine–garnet Fe–Mg thermometer.  "Including measured
garnet Fe3+" uses Fe2+ = Fe_tot × (1 − XFe3+) on the garnet site.
"Not including" treats all garnet Fe as Fe2+.  TA98 is two-pyroxene
and does not use garnet Fe3+.  Both thermometers are evaluated at the
experimental pressure.  The pipeline's flat +114 °C WZ07 offset is
NOT applied.

MaLm2 (Rzehak) has no olivine and no cpx, so it cannot be plotted.
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

from calculations.pt_fo2.cations import convert_oxides_to_cations
from calculations.pt_fo2.screening import apply_equilibrium_screening
from calculations.pt_fo2.thermobarometry import (
    ta98_temperature, wz07_temperature,
    bk90_temperature, ng10_temperature,
    ha84_temperature, kr00_temperature, na09_temperature,
)

IN_COMP = ROOT / "data" / "incoming" / "experiments_20260913.csv"
OUT = ROOT / "docs" / "figures"
OUT.mkdir(parents=True, exist_ok=True)


def _clean_cols(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [str(c).replace("\xa0", "").strip() for c in df.columns]
    return df


def assemble() -> pd.DataFrame:
    raw = _clean_cols(pd.read_csv(IN_COMP))
    xfe3 = pd.to_numeric(raw.get("garnet_Fe3+_Fetot"), errors="coerce")
    raw = raw.loc[xfe3.notna()].copy()
    raw = raw.rename(columns={
        "Sample ID": "SampleID",
        "Temperature": "T_exp",
        "Pressure_kbar": "P_kbar",
        "garnet_Fe3+_Fetot": "XFe3",
    })
    raw["XFe3"] = pd.to_numeric(raw["XFe3"], errors="coerce")
    return raw.reset_index(drop=True)


def temperatures(df: pd.DataFrame) -> pd.DataFrame:
    cat = convert_oxides_to_cations(df, verbose=False)
    P = pd.to_numeric(cat["P_kbar"], errors="coerce").to_numpy(dtype=float)
    cat["T_TA98"] = ta98_temperature(cat, P)
    cat["T_WZ07_allFe"] = wz07_temperature(cat, P)
    cat["T_WZ07"] = cat["T_WZ07_allFe"]
    cat["T_BK90"] = bk90_temperature(cat, P)
    cat["T_NG10"] = ng10_temperature(cat["T_BK90"])
    cat["T_HA84"] = ha84_temperature(cat, P)
    cat["T_KR00"] = kr00_temperature(cat, P)
    cat["T_NA09"] = na09_temperature(cat, P)

    fe2 = cat.copy()
    xfe3 = pd.to_numeric(fe2["XFe3"], errors="coerce").fillna(0.0)
    fe2["garnet_Fe"] = pd.to_numeric(fe2["garnet_Fe"], errors="coerce") * (1.0 - xfe3)
    cat["T_WZ07_Fe2"] = wz07_temperature(fe2, P)

    has_ol = pd.to_numeric(cat.get("olivine_Mg"), errors="coerce").notna() & (
        pd.to_numeric(cat.get("olivine_Fe"), errors="coerce") > 0
    )
    has_cpx = pd.to_numeric(cat.get("cpx_Ca"), errors="coerce").notna()
    has_opx = pd.to_numeric(cat.get("opx_Ca"), errors="coerce").notna()
    cat["usable"] = has_ol & has_cpx & has_opx & cat["XFe3"].notna()
    cat["dT_Fe2"] = cat["T_TA98"] - cat["T_WZ07_Fe2"]
    cat["dT_allFe"] = cat["T_TA98"] - cat["T_WZ07_allFe"]

    cat, stats = apply_equilibrium_screening(cat, verbose=True)
    reasons = []
    for _, r in cat.iterrows():
        if r.get("screen_pass") is True:
            reasons.append("")
            continue
        bits = []
        if not r.get("screen_4phases"):
            bits.append("4-phase")
        elif not r.get("screen_T900"):
            bits.append("T<900")
        elif not r.get("screen_opx_cpx"):
            bits.append("opx–cpx")
        elif not r.get("screen_pyx_gt"):
            bits.append("pyx–gt")
        elif not r.get("screen_gt_ol"):
            bits.append("gt–ol")
        reasons.append("+".join(bits))
    cat["fail_reason"] = reasons
    return cat


STYLE = {
    "Matjuschkin et al. 2014": dict(marker="o", color="#333333", label="Matjuschkin 2014"),
    "Rzehak et al. 2020": dict(marker="s", color="#4c72b0", label="Rzehak 2020"),
}

# Annotated groups on the Fe³⁺-included figure
RED_GROUP = {
    "F1100r", "F1100o", "F1200o", "F1200r",
    "F1300r", "F14002r", "F14002o",
    "MaLm6", "MaLm11", "MaLm1",
}
GREEN_GROUP = {"F1300o", "MaLm12", "MaLm9"}


def _sid(s) -> str:
    return str(s).replace(" ", "")


def _hull_xy(x, y, pad=0.12):
    """Closed polygon around points, expanded from the centroid."""
    pts = np.column_stack([np.asarray(x, float), np.asarray(y, float)])
    if len(pts) == 0:
        return None
    if len(pts) == 1:
        cx, cy = pts[0]
        return np.array([[cx - 20, cy - 20], [cx + 20, cy - 20],
                         [cx + 20, cy + 20], [cx - 20, cy + 20],
                         [cx - 20, cy - 20]])
    if len(pts) == 2:
        (x0, y0), (x1, y1) = pts
        dx, dy = x1 - x0, y1 - y0
        nrm = np.hypot(dx, dy) + 1e-9
        px, py = -dy / nrm * 25, dx / nrm * 25
        return np.array([[x0 + px, y0 + py], [x1 + px, y1 + py],
                         [x1 - px, y1 - py], [x0 - px, y0 - py],
                         [x0 + px, y0 + py]])
    try:
        from scipy.spatial import ConvexHull
        hull = ConvexHull(pts)
        vert = pts[hull.vertices]
    except Exception:
        order = np.argsort(np.arctan2(pts[:, 1] - pts[:, 1].mean(),
                                      pts[:, 0] - pts[:, 0].mean()))
        vert = pts[order]
    c = vert.mean(axis=0)
    vert = c + (1.0 + pad) * (vert - c)
    return np.vstack([vert, vert[0]])


def _plot(df: pd.DataFrame, ycol: str, ylabel: str, title: str, stem: str,
          href=None, href_label=None, groups=False) -> None:
    from matplotlib.patches import Polygon

    use = df[df["usable"]].copy()
    use["_sid"] = use["SampleID"].map(_sid)
    fig, ax = plt.subplots(figsize=(119 / 25.4, 95 / 25.4))
    ax.axhline(0, color="0.6", lw=0.7)
    if href is not None:
        ax.axhline(href, color="0.45", lw=0.8, ls=":")
        ax.text(use["T_TA98"].min(), href + (12 if href >= 0 else -18),
                href_label or f"{href:.0f} °C",
                fontsize=7, color="0.35")

    if groups:
        group_spec = [
            (RED_GROUP, 1, "#5e3c99", (0, (6, 3)), "1"),
            (GREEN_GROUP, 2, "#b35806", (0, (2, 2)), "2"),
        ]
        xx_all = np.linspace(use["T_TA98"].min() - 30, use["T_TA98"].max() + 30, 80)
        for ids, num, color, dashes, label in group_spec:
            g = use[use["_sid"].isin(ids)]
            if g.empty:
                continue
            poly = _hull_xy(g["T_TA98"], g[ycol], pad=0.28)
            if poly is not None:
                ax.add_patch(Polygon(
                    poly, closed=True, fill=False, edgecolor=color,
                    linewidth=1.3, linestyle=dashes, zorder=2,
                ))
            if len(g) >= 2:
                a, b = np.polyfit(g["T_TA98"].to_numpy(float),
                                  g[ycol].to_numpy(float), 1)
                yy = a * xx_all + b
                ax.plot(xx_all, yy, color=color, lw=1.4, zorder=2)
                # label 1 / 2 near the mid-T of the group
                tx = float(g["T_TA98"].mean())
                ty = a * tx + b
                ax.text(tx, ty + (18 if num == 1 else -22), label,
                        color=color, fontsize=11, fontweight="bold",
                        ha="center", va="center", zorder=5,
                        bbox=dict(boxstyle="circle,pad=0.15",
                                  facecolor="white", edgecolor=color, lw=0.8))
                eq = rf"{label}:  $T_{{\mathrm{{TA98}}}}-T_{{\mathrm{{WZ07}}}}$ = {a:.3f}$\,T$ {b:+.0f}"
                if num == 1:
                    ax.text(0.03, 0.97, eq, transform=ax.transAxes,
                            color=color, fontsize=7.5, va="top", ha="left",
                            bbox=dict(boxstyle="round,pad=0.25",
                                      facecolor="white", edgecolor=color,
                                      lw=0.6, alpha=0.92))
                else:
                    ax.text(0.03, 0.03, eq, transform=ax.transAxes,
                            color=color, fontsize=7.5, va="bottom", ha="left",
                            bbox=dict(boxstyle="round,pad=0.25",
                                      facecolor="white", edgecolor=color,
                                      lw=0.6, alpha=0.92))
                print(f"  group {label}: n={len(g)}  "
                      f"TA98−WZ07 = {a:.3f}·T + {b:.1f}")

    for ref, sty in STYLE.items():
        sub = use[use["Reference"] == ref]
        if sub.empty:
            continue
        ax.scatter(sub["T_TA98"], sub[ycol], s=42, edgecolors="k",
                   linewidths=0.4, zorder=3, **sty)
        fail = sub["screen_pass"] != True
        if fail.any():
            ax.scatter(sub.loc[fail, "T_TA98"], sub.loc[fail, ycol],
                       marker="x", s=70, c="k", linewidths=1.2, zorder=4)
        for _, r in sub.iterrows():
            ax.annotate(str(r["SampleID"]), (r["T_TA98"], r[ycol]),
                        fontsize=6, xytext=(4, 3), textcoords="offset points")

    ax.set_xlabel("TA98 temperature (°C)")
    ax.set_ylabel(ylabel)
    ax.legend(frameon=True, facecolor="white", edgecolor="0.55",
              framealpha=0.92, fontsize=8, loc="upper right")
    ax.tick_params(direction="in", top=True, right=True)
    fig.tight_layout()
    fig.savefig(OUT / f"{stem}.png", dpi=600)
    fig.savefig(OUT / f"{stem}.pdf")
    plt.close(fig)
    print(f"  wrote {OUT / (stem + '.png')}")


def main() -> None:
    df = temperatures(assemble())
    use = df[df["usable"]]
    cols = ["SampleID", "Reference", "T_exp", "P_kbar", "XFe3",
            "T_TA98", "T_WZ07_allFe", "T_WZ07_Fe2", "dT_allFe", "dT_Fe2",
            "screen_pass", "fail_reason", "dT_opx_cpx", "dT_KR00", "dT_gt_ol"]
    show = use[cols].copy()
    for c in ["T_exp", "P_kbar", "XFe3", "T_TA98", "T_WZ07_allFe", "T_WZ07_Fe2",
              "dT_allFe", "dT_Fe2", "dT_opx_cpx", "dT_KR00", "dT_gt_ol"]:
        show[c] = pd.to_numeric(show[c], errors="coerce").round(2)
    print(show.to_string(index=False))
    table = ROOT / "data" / "derived" / "exp_TA98_WZ07_Fe3.csv"
    table.parent.mkdir(parents=True, exist_ok=True)
    use[cols].to_csv(table, index=False)
    print(f"  table {table}")

    f13 = use[use["SampleID"].astype(str).str.replace(" ", "") == "F1300o"]
    href = float(f13["dT_Fe2"].iloc[0]) if len(f13) else -114.0
    print(f"F1300o TA98−WZ07 (Fe2+ only) = {href:.2f} °C  "
          f"(flat offset locked at {118} °C)")
    _plot(df, "dT_Fe2",
          r"$T_{\mathrm{TA98}}-T_{\mathrm{WZ07}}$ (°C), Fe$^{2+}$ only",
          "",
          "Fig04_TA98_WZ07_measured_Fe3",
          href=href,
          href_label=f"F1300o offset = {href:.0f} °C",
          groups=True)


if __name__ == "__main__":
    main()
