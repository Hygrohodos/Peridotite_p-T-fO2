#!/usr/bin/env python3
"""Online Resource 1 material: diagnostics of the WITHDRAWN rev. v2 correction (level anchored to the Slave-25 mean plus
the weighted experimental slope; commit 0f43eb7), quoted in manuscript Section 3.3 as the reason the slope is not applied.
Inputs: data/reference/fO2_residuals_v2.csv and p-T-fO2_v2.csv (git show 0f43eb7:data/derived/...), and TA98 from
data/derived/screened_PT.csv (TA98 does not depend on the Fe3+ correction). Writes JSON to --out or stdout."""
import argparse, json
from pathlib import Path
import pandas as pd
from scipy import stats
R = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser(); ap.add_argument("--out"); a = ap.parse_args()
r = pd.read_csv(R / "data/reference/fO2_residuals_v2.csv")
scr = pd.read_csv(R / "data/derived/screened_PT.csv", low_memory=False); scr = scr.rename(columns={scr.columns[0]: "SampleID"})
T = scr.drop_duplicates("SampleID").set_index("SampleID")["T_TA98"]
r["T"] = r.SampleID.map(T); r["d"] = r.XFe3_est - r.XFe3_meas
out = {"source": "rev. v2 (0f43eb7), withdrawn"}
for lab, m in [("excl_error_groups_52", ~r.error_group), ("all_57", r.d.notna())]:
    q = r[m]; lr = stats.linregress(q["T"], q.d)
    out[f"resid_vs_T_{lab}"] = dict(n=int(len(q)), slope_per_100C=100 * lr.slope, se_per_100C=100 * lr.stderr, p=lr.pvalue)
p = pd.read_csv(R / "data/reference/p-T-fO2_v2.csv"); s = p[p.Screen_pass == True]; fl = s["XGtFe3+_corr"] <= 0.0105
out["floor"] = dict(n=int(fl.sum()), of=int(len(s)), share=float(fl.mean()), below_1162C=int((fl & (s.T_TA98 < 1162)).sum()))
v = r[(r.dataset != "Slave calibration") & (r.XFe3_est <= 0.0105)]
out["validation_floor_samples"] = v[["SampleID", "XFe3_meas", "res_log"]].to_dict("records")
txt = json.dumps(out, indent=1, default=float)
(Path(a.out).write_text(txt) if a.out else print(txt))
