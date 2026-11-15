#!/usr/bin/env python3
"""Reproduce the coefficients of the TA98 + garnet-Cr2O3 correction (Eq. 4, revision v7_cr) from the pipeline.

    python3 scripts/fit_cr_correction.py [--out cr_fit.json] [--check]

Method (full description in calculations/pt_fo2/cr_correction.py):
  natural: screened Slave calibration + Siberia/Kaapvaal validation xenoliths with measured garnet Fe3+,
           D = T_TA98 - T_WZ07(measured Fe3+), T_TA98 >= 1120 C, without the error groups and Uv-50-04 (n = 37);
  experiments: the 8 Matjuschkin et al. (2014) runs (garnet Cr2O3 known), sigma_i from
           scripts/weighted_thermometer_regression.py;
  joint WLS  y = a*nat + e*exp + b*t + c*(Cr2O3 - Cr_ref)*t,  t = T_TA98 - pivot, weights 1/sigma^2, fixed scale,
           natural sigma = residual SD of the natural-only OLS fit of y ~ t + (CaO - CaO_ref)*t (17.43 C; the
           same natural weighting was used for every composition model compared).
--check exits non-zero unless the fitted a, b, c and Cr_ref equal calculations/pt_fo2/params.py to 1e-6.
"""
import argparse, json, sys
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "scripts")]
from calculations import run_full_pipeline                       # noqa: E402
from calculations.pt_fo2 import calibration as calib             # noqa: E402
from calculations.pt_fo2 import cr_correction as crc             # noqa: E402
from calculations.pt_fo2.level_anchor import discrepancy_measured  # noqa: E402
from calculations.pt_fo2.params import PARAMS                    # noqa: E402
import weighted_thermometer_regression as wtr                    # noqa: E402

INPUT = ROOT / "data" / "master" / "xenolith_geochem.csv"
EXPERIMENTS = ROOT / "data" / "incoming" / "experiments_20260913.csv"


def natural_samples(table):
    """The 37 natural fit samples: id, craton, T, t, y = D, garnet Cr2O3 and CaO."""
    m = (calib.calibration_mask(table) | calib.validation_mask(table)) & calib.measured_fe3(table).notna()
    t = table[m]
    g = pd.DataFrame({"id": t["SampleID"].astype(str).str.strip().values,
                      "craton": t["Craton"].astype(str).values,
                      "T": pd.to_numeric(t["T_TA98"], errors="coerce").values,
                      "y": discrepancy_measured(table)[m].values,
                      "Cr2O3": pd.to_numeric(t["garnet_Cr2O3"], errors="coerce").values,
                      "CaO": pd.to_numeric(t["garnet_CaO"], errors="coerce").values})
    g = g[(g["T"] >= crc.NATURAL_T_MIN_FIT) & ~g["id"].isin(crc.FIT_EXCLUDED)].reset_index(drop=True)
    g["t"] = g["T"] - crc.CR_PIVOT
    return g


def experiment_runs():
    """The 8 Matjuschkin runs with garnet Cr2O3, their y, t and sigma_i (paper weights)."""
    runs = wtr.load_runs()
    wt, _ = wtr.fit(runs, t_ref=crc.CR_PIVOT)
    runs = runs.merge(wt[["SampleID", "sigma_i_C", "weight_pct"]], on="SampleID")
    raw = pd.read_csv(EXPERIMENTS).drop_duplicates("Sample ID").set_index("Sample ID")
    runs["Cr2O3"] = pd.to_numeric(raw.loc[runs.SampleID, "garnet_Cr2O3"], errors="coerce").values
    runs = runs[runs.SampleID.str.startswith("F") & runs.Cr2O3.notna()].reset_index(drop=True)
    return pd.DataFrame({"id": runs.SampleID, "T": runs.T_TA98, "t": runs.T_TA98 - crc.CR_PIVOT,
                         "y": runs.y_TA98_minus_WZ07_C, "Cr2O3": runs.Cr2O3, "sig": runs.sigma_i_C})


def fit(table):
    g = natural_samples(table)
    e = experiment_runs()
    cr_ref = float(g.Cr2O3.mean()); ca_ref = float(g.CaO.mean())
    nat_sd = float(np.sqrt(smf.ols("y ~ t + t:dCaO", g.assign(dCaO=g.CaO - ca_ref)).fit().scale))
    D = pd.concat([g.assign(sig=nat_sd, nat=1, exp=0), e.assign(nat=0, exp=1)], ignore_index=True)
    D["dCr"] = D.Cr2O3 - cr_ref
    m = smf.wls("y ~ 0 + nat + exp + t + t:dCr", D, weights=1 / D.sig ** 2).fit(cov_type="fixed scale")
    nat_only = smf.ols("y ~ t + t:dCr", g.assign(dCr=g.Cr2O3 - cr_ref)).fit()
    p = lambda k: dict(value=float(m.params[k]), se=float(m.bse[k]), p=float(m.pvalues[k]))
    cr_exp = float(e.Cr2O3.mean())
    out = dict(
        n_natural=int(len(g)), n_by_craton=g.craton.value_counts().to_dict(), n_experiments=int(len(e)),
        natural_T_range_C=[float(g["T"].min()), float(g["T"].max())],
        natural_Cr2O3_range_wt=[float(g.Cr2O3.min()), float(g.Cr2O3.max())],
        experiment_T_range_C=[float(e["T"].min()), float(e["T"].max())],
        experiment_Cr2O3_range_wt=[float(e.Cr2O3.min()), float(e.Cr2O3.max())],
        pivot_C=crc.CR_PIVOT, cr_ref_wt=cr_ref, natural_sigma_C=nat_sd,
        a=p("nat"), experiment_level=p("exp"), b=p("t"), c=p("t:dCr"),
        chi2=float(np.sum(m.resid ** 2 / D.sig ** 2)), n_points=int(len(D)),
        natural_only=dict(a=float(nat_only.params.Intercept), b=float(nat_only.params.t),
                          c=float(nat_only.params["t:dCr"]), c_se=float(nat_only.bse["t:dCr"])),
        slope_check=dict(mean_experiment_Cr2O3_wt=cr_exp,
                         predicted_slope=float(m.params.t + m.params["t:dCr"] * (cr_exp - cr_ref)),
                         measured_weighted_slope_13_runs=float(wtr.fit(wtr.load_runs(), t_ref=crc.CR_PIVOT)[1].b1)),
    )
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path)
    ap.add_argument("--check", action="store_true")
    a_ = ap.parse_args()
    table, _ = run_full_pipeline(INPUT, verbose=False)
    res = fit(table)
    print(json.dumps(res, indent=1, default=float))
    if a_.out:
        a_.out.write_text(json.dumps(res, indent=1, default=float))
    if a_.check:
        fp = PARAMS.fe3
        ok = all(abs(x - y) < 1e-6 for x, y in ((res["a"]["value"], fp.cr_a), (res["b"]["value"], fp.cr_b),
                                                (res["c"]["value"], fp.cr_c), (res["cr_ref_wt"], fp.cr_ref)))
        sys.exit(0 if ok else "fitted coefficients differ from calculations/pt_fo2/params.py")
