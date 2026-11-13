#!/usr/bin/env python3
"""Table 1 and headline numbers: original (main 87024f0) / rev. v2 (0f43eb7, Slave-25 mean + slope) / this revision.
Reads data/reference/{Table1,p-T-fO2}_{original,v2}.csv and the current data/derived outputs; writes Markdown to --out or stdout."""
import argparse
from pathlib import Path
import pandas as pd
R = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser(); ap.add_argument("--out")
ap.add_argument("--derived", default=str(R / "data" / "derived"), help="folder with Table1_calibration_stats.csv and p-T-fO2.csv")
a = ap.parse_args(); D = Path(a.derived)
src = {"Original": (R/"data/reference/Table1_original.csv", R/"data/reference/p-T-fO2_original.csv"),
       "Rev. v2 (Slave-25 mean + slope)": (R/"data/reference/Table1_v2.csv", R/"data/reference/p-T-fO2_v2.csv"),
       "This revision": (D/"Table1_calibration_stats.csv", D/"p-T-fO2.csv")}
rows = ["Slave calibration", "Slave, error groups excluded", "Siberia + Kaapvaal validation", "Calibration + validation", "Combined, error groups excluded"]
L = ["| Dataset | n | " + " | ".join(src) + " |", "|---|---|" + "---|" * len(src)]
T = {k: pd.read_csv(v[0]).set_index("dataset") for k, v in src.items()}
for r in rows:
    cells = [f"{T[k].loc[r,'p02']*100:.0f} / {T[k].loc[r,'p05']*100:.0f} %, 1σ {T[k].loc[r,'sigma']:.3f}, {T[k].loc[r,'sigma_log']:.2f} log" for k in src]
    L.append(f"| {r} | {int(T['Original'].loc[r,'n'])} | " + " | ".join(cells) + " |")
L += ["", "Cells: within ±0.02 / ±0.05 Fe3+/ΣFe; 1σ Fe3+/ΣFe; 1σ ΔFMQ (log units). Charge-balance rows are unchanged (23/60 %, 0.085; 25/53 %, 0.087).", "",
      "| Fully screened xenoliths (n = 450) | " + " | ".join(src) + " |", "|---|" + "---|" * len(src)]
P = {k: pd.read_csv(v[1]) for k, v in src.items()}
def m(k):
    s = P[k][P[k].Screen_pass == True]; fl = s["XGtFe3+_corr"] <= 0.0105; F = s["logfO2(ΔFMQ)"]
    sh = (s.P_GPa < 4.5) & (F <= -4.5)
    return dict(n=len(s), floor=f"{fl.sum()} ({fl.mean()*100:.1f} %)", med=f"{F.median():.2f}", mednf=f"{F[~fl].median():.2f}",
                sh=f"{(sh&fl).sum()} of {sh.sum()}")
M = {k: m(k) for k in src}
for lab, key in [("At the 0.01 starting value (floor)", "floor"), ("Median ΔFMQ, all", "med"), ("Median ΔFMQ, floor excluded", "mednf"),
                 ("Shallow (< 4.5 GPa) and ΔFMQ ≤ −4.5: at floor / total", "sh")]:
    L.append(f"| {lab} | " + " | ".join(str(M[k][key]) for k in src) + " |")
txt = "\n".join(L) + "\n"
(Path(a.out).write_text(txt) if a.out else print(txt))
