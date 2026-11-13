#!/usr/bin/env python3
"""Numbers quoted in the text that are not in Table 1 (writes JSON to --out or stdout).
Floor = garnet Fe3+/SumFe at the 0.01 starting value (XGtFe3+_corr <= 0.0105)."""
import argparse, json, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "scripts")]
ap = argparse.ArgumentParser(); ap.add_argument("--out")
ap.add_argument("--derived", default=str(ROOT / "data" / "derived"), help="folder with screened_PT.csv, fO2_residuals.csv, fO2_uncertainty.json")
a_ = ap.parse_args()
D = Path(a_.derived)
# full-precision working table (Geology round 2, N9: no binning on the rounded publication columns)
s = pd.read_csv(D / "screened_PT.csv", low_memory=False); s = s.rename(columns={s.columns[0]: "SampleID"})
s = s[s["screen_pass"] == True].copy()
s["P_GPa"] = pd.to_numeric(s["P_CA91"], errors="coerce") / 10.0
F = pd.to_numeric(s["delta_logfO2_FMQ"], errors="coerce"); X = pd.to_numeric(s["XFe3_gt_corr"], errors="coerce"); fl = X <= 0.0105
out = {}
out["n_screened"] = int(len(s)); out["n_floor"] = int(fl.sum()); out["floor_share"] = float(fl.mean())
out["median_dFMQ_all"] = float(F.median()); out["median_dFMQ_nonfloor"] = float(F[~fl].median()); out["median_dFMQ_floor"] = float(F[fl].median())
out["mean_dFMQ_all"] = float(F.mean()); out["mean_dFMQ_nonfloor"] = float(F[~fl].mean())
T = s.T_TA98; P = s.P_GPa
out["floor_by_T"] = {lab: [int((fl & m).sum()), int(m.sum())] for lab, m in
    [("<1000", T < 1000), ("1000-1100", (T >= 1000) & (T < 1100)), ("1100-1200", (T >= 1100) & (T < 1200)), (">=1200", T >= 1200)]}
out["floor_by_P"] = {lab: [int((fl & m).sum()), int(m.sum())] for lab, m in
    [("<3.5", P < 3.5), ("3.5-4.5", (P >= 3.5) & (P < 4.5)), ("4.5-5.5", (P >= 4.5) & (P < 5.5)), (">=5.5", P >= 5.5)]}
out["floor_below_1100"] = int((fl & (T < 1100)).sum()); out["floor_below_1000"] = int((fl & (T < 1000)).sum())
sh = (P < 4.5) & (F <= -4.5); out["shallow_reduced"] = [int((sh & fl).sum()), int(sh.sum())]
out["shallow_reduced_kaapvaal"] = int((sh & s.Craton.str.contains("Kaap", na=False)).sum())
out["P_below_4.5_floor"] = [int((fl & (P < 4.5)).sum()), int((P < 4.5).sum())]
out["floor_P_ge_4.5"] = int((fl & (P >= 4.5)).sum())
# deep reduced (BSE-relative check is done in the figure; here dFMQ < -3.5 at P >= 5 GPa as a simple count)
deep = P >= 5.0
out["deep_P_ge_5"] = int(deep.sum()); out["deep_dFMQ_lt_m3.5"] = [int((deep & (F < -3.5)).sum()), int((deep & (F < -3.5) & fl).sum())]
# deep (>= 5 GPa) samples vs BSE, and a simple scatter test (Geology round 2, N2):
# if every off-floor deep sample sat at the deep median and differed only by the reduced-half fO2 1sigma,
# how many would fall below BSE by chance?
from calculations.pt_fo2.buffers import buffer_table
from scipy.stats import norm
bse = np.asarray(buffer_table(s["P_CA91"].to_numpy(float), df=s)["dBSE"]); Fv = F.to_numpy(float)
deep5 = (s["P_GPa"] >= 5.0).to_numpy(); flv = fl.to_numpy(); sig = json.load(open(D / "fO2_uncertainty.json"))["reduced"]["sigma_log"]
m = deep5 & ~flv; med = float(np.median(Fv[m])); marg = bse - Fv
out["deep_vs_BSE"] = dict(n_deep=int(deep5.sum()), n_below=int((deep5 & (marg > 0)).sum()), n_below_floor=int((deep5 & (marg > 0) & flv).sum()),
    n_offfloor=int(m.sum()), offfloor_below=int((m & (marg > 0)).sum()), deep_median_offfloor=med, sigma_log=float(sig),
    expected_below_from_scatter=float(norm.cdf((bse[m] - med) / sig).sum()),
    below_within_1sigma=int((deep5 & (marg > 0) & (marg <= sig)).sum()), below_more_than_2sigma=int((deep5 & (marg > 2 * sig)).sum()),
    below_more_than_2sigma_offfloor=int((m & (marg > 2 * sig)).sum()),
    ids_below=s.loc[deep5 & (marg > 0), "SampleID"].astype(str).tolist())
plotted = (F.between(-7, 0) & s["P_CA91"].between(30, 70)).to_numpy()
out["fig11_floor_plotted"] = [int((plotted & flv).sum()), int(flv.sum())]
# residual trend with T on the 52 / 57
r = pd.read_csv(D / "fO2_residuals.csv"); scr = pd.read_csv(D / "screened_PT.csv", low_memory=False)
scr = scr.rename(columns={scr.columns[0]: "SampleID"})
tt = scr.drop_duplicates("SampleID").set_index("SampleID")["T_TA98"]
r["T"] = r.SampleID.map(tt); r["d"] = r.XFe3_est - r.XFe3_meas
slave25 = (r.dataset == "Slave calibration") & ~r.error_group
for lab, m in [("all57", r.d.notna()), ("excl52", ~r.error_group), ("slave25", slave25), ("slave25_T_ge_1130", slave25 & (r["T"] >= 1130))]:
    q = r[m]; lr = stats.linregress(q["T"], q.d)
    out[f"resid_vs_T_{lab}"] = dict(n=int(len(q)), slope_per_100C=float(100 * lr.slope), se_per_100C=float(100 * lr.stderr), p=float(lr.pvalue))
out["resid_by_T_excl52"] = {lab: [float(r[~r.error_group & m].d.mean()), int((~r.error_group & m).sum())] for lab, m in
    [("<1100", r["T"] < 1100), ("1100-1200", (r["T"] >= 1100) & (r["T"] < 1200)), ("1200-1300", (r["T"] >= 1200) & (r["T"] < 1300)), (">=1300", r["T"] >= 1300)]}
out["cal_val_T_coverage"] = dict(below_1000=int((r["T"] < 1000).sum()), below_1100=int((r["T"] < 1100).sum()), min=float(r["T"].min()), max=float(r["T"].max()))
fl57 = r.XFe3_est <= 0.0105
out["cal_val_floor"] = r.loc[fl57, ["SampleID", "dataset", "XFe3_meas", "res_log"]].to_dict("records")
out["fo2_sigma_nonfloor57"] = [float(r.loc[~fl57, "res_log"].std(ddof=1)), int((~fl57).sum())]
ou = r[r.dataset != "Slave calibration"]; out["outside_005_validation"] = ou.loc[(ou.d.abs() > 0.05 + 1e-9), ["SampleID", "XFe3_meas", "XFe3_est", "d"]].to_dict("records")
c = r[r.dataset == "Slave calibration"]; out["outside_005_calibration"] = c.loc[(c.d.abs() > 0.05 + 1e-9), ["SampleID", "error_group", "d"]].to_dict("records")
eg = r[r.error_group]; out["error_groups"] = eg[["SampleID", "XFe3_meas", "XFe3_est", "d", "dFMQ_est"]].to_dict("records")
out["sample_21_6"] = float(r.loc[r.SampleID == "21-6", "d"].iloc[0]); out["Uv-50-04"] = float(r.loc[r.SampleID == "Uv-50-04", "d"].iloc[0])
txt = json.dumps(out, indent=1, default=float)
(Path(a_.out).write_text(txt) if a_.out else print(txt))
