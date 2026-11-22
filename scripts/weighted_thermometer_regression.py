#!/usr/bin/env python3
"""Weighted TA98 - WZ07 regression on experiments (manuscript Fig. 5, Sect. 3.3).

Stand-alone module and CLI. It reads a CSV of experimental run products,
weights each run by how far it is from equilibrium, and fits

    y_i = T_TA98,i - T_WZ07,i  =  b0 + b1 (T_TA98,i - T_ref)        [deg C]

by weighted least squares (WLS). Adding a run = adding a row to the CSV; the
same weights and fit are then applied to every run consistently.

Dependencies: Python >= 3.9, numpy >= 1.24, pandas >= 2.0 (matplotlib only
for the figure script; pytest optional for the test). No other repo imports,
so it runs from a clean checkout:

    python3 scripts/weighted_thermometer_regression.py            # defaults
    python3 scripts/weighted_thermometer_regression.py --input my_runs.csv \
        --outdir results/ --t-ref 1350 --sigma0 50

Input CSV schema (one row per run; extra columns are ignored; all T in deg C)
-----------------------------------------------------------------------------
    SampleID            run identifier, unique               (required)
    T_TA98              two-pyroxene temperature, deg C       (required)
    T_WZ07_Fe2          olivine-garnet T with measured Fe3+ removed (Fe2+
                        only), deg C                          (required)
    dT_opx_cpx          TA98 - NG10, deg C                    (required)
    dT_pyx_gt           HA84 - quadratic(TA98), deg C         (required)
    dT_gt_ol            WZ07pred - WZ07, deg C                (required)
    Reference           source publication                    (optional)
    sigma_analytical_C  per-run analytical 1-sigma on y, deg C (optional;
                        blank / missing -> SIGMA0 = 50 deg C. This is a
                        PLACEHOLDER until analytical quality is assessed.
                        Non-numeric, non-finite or <= 0 values are errors.)
    include             use the run in the fit (optional). Blank = included;
                        1/true/yes/y = included; 0/false/no/n = excluded;
                        anything else is an error.
Required values must be finite numbers: a blank or NaN T or dT is an error
(the offending rows are listed); nothing is silently dropped.
Column names can be remapped with --col NAME=CSVCOLUMN (NAME must be one of
the canonical names above; clashes with an existing column are errors).
Sanity warnings (no effect on the fit): median temperature > 1600 (Kelvin?),
any T > 2000 or < 0, and near-duplicate runs (all five temperatures within
DUP_TOL_C of another run). At least MIN_RUNS = 3 included runs are needed.
The default input, data/derived/exp_TA98_WZ07_Fe3.csv, is written by
scripts/plot_exp_ta98_wz07.py from the raw experimental analyses.

Weights
-------
    excess_i  = sqrt( sum_j max(0, |dT_ij| - thr_j)^2 )   screens j (a run
                fails screen j when |dT_ij| > thr_j, i.e. exactly when that
                term is non-zero):
                opx-cpx 80, pyx-gt 60, gt-ol 50 deg C (SCREEN_THRESHOLDS_C)
    sigma_i^2 = sigma_analytical_i^2 + excess_i^2         (sigma_analytical
                defaults to SIGMA0 = 50 deg C, the gt-ol acceptance width)
    w_i       = 1 / sigma_i^2;  shown as w_i / sum(w) (per cent)
    n_eff     = (sum w)^2 / sum w^2   (Kish effective sample size)

Fit and uncertainty convention
------------------------------
WLS with weights w_i. Parameter covariance = (X^T W X)^-1 * s, where
s = max(1, chi2_red) and chi2_red = sum w_i r_i^2 / (n - p). That is,
sigma_i are taken as absolute 1-sigma uncertainties, and the covariance is
inflated only if the scatter exceeds them (never shrunk). For the current
13 runs chi2_red = 0.84, so s = 1 and the SEs are the unscaled ones
(-14 +/- 27, -0.555 +/- 0.216). This is the convention of
calculations/pt_fo2/exp_calibration.py (band_cov), which the pipeline uses.
The 95 % band is fit +/- 1.96 * SE(prediction).

Methods caveats: sigma_i is treated as an absolute 1-sigma but sigma0 is a
placeholder and the excess term is heuristic, so the absolute scale is an
assumption. The weights are built from the same thermometers as y (dT_opx_cpx
contains TA98; for the 13 current runs r(dT_opx_cpx, y) = -0.48), so they are
not independent of the response; the case bootstrap carries weights with runs.
Exit status of the CLI: 0 OK, 2 input/validation error.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import warnings
from dataclasses import dataclass, field, asdict
from pathlib import Path

import numpy as np
import pandas as pd

__version__ = "1.0"

# ---- named constants (all overridable from the CLI / function arguments) ----
SCREEN_THRESHOLDS_C = {"dT_opx_cpx": 80.0, "dT_pyx_gt": 60.0, "dT_gt_ol": 50.0}
SIGMA0_C = 50.0          # analytical / equilibrium-floor 1-sigma (placeholder)
T_REF_C = 1350.0         # centring temperature for the slope
Z95 = 1.959964

MIN_RUNS = 3             # minimum included runs for a fit with a usable chi2
KELVIN_WARN_MEDIAN_C = 1600.0   # median T above this looks like Kelvin
T_PLAUSIBLE_C = (0.0, 2000.0)
DUP_TOL_C = 0.5          # near-duplicate: all five temperatures within this
TRUE_VALUES = {"1", "1.0", "true", "yes", "y", "t"}
FALSE_VALUES = {"0", "0.0", "false", "no", "n", "f"}

REQUIRED = ("SampleID", "T_TA98", "T_WZ07_Fe2", "dT_opx_cpx", "dT_pyx_gt", "dT_gt_ol")
OPTIONAL = ("Reference", "sigma_analytical_C", "include")
REPO = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = REPO / "data" / "derived" / "exp_TA98_WZ07_Fe3.csv"


class InputError(ValueError):
    """Invalid input data or options (CLI exits with status 2)."""


class DataSanityWarning(UserWarning):
    """Input looks suspicious (units, duplicates); the fit is still run."""


def _rows(df, mask):
    ids = df.loc[mask, "SampleID"].astype(str).tolist() if "SampleID" in df else []
    idx = [int(i) + 2 for i in np.flatnonzero(np.asarray(mask))]   # CSV line numbers
    return ", ".join(f"{a} (line {b})" for a, b in zip(ids, idx)) if ids else str(idx)


def parse_colmap(items):
    """['NAME=CSVCOL', ...] -> {NAME: CSVCOL}, with validation."""
    out = {}
    for it in items or []:
        if it.count("=") != 1:
            raise InputError(f"--col {it!r}: expected NAME=CSVCOLUMN")
        k, v = (x.strip() for x in it.split("="))
        if not k or not v:
            raise InputError(f"--col {it!r}: NAME and CSVCOLUMN must be non-empty")
        if k not in REQUIRED + OPTIONAL:
            raise InputError(f"--col {it!r}: unknown NAME {k!r}; allowed: {', '.join(REQUIRED + OPTIONAL)}")
        if k in out:
            raise InputError(f"--col: {k!r} mapped twice")
        out[k] = v
    if len(set(out.values())) != len(out):
        raise InputError(f"--col: two names mapped to the same CSV column ({out})")
    return out


def _apply_colmap(df, colmap):
    for k, v in colmap.items():
        if v not in df.columns:
            raise InputError(f"--col {k}={v}: CSV has no column {v!r}")
        if k != v and k in df.columns:
            raise InputError(f"--col {k}={v}: CSV already has a column {k!r}; "
                             f"drop or rename it, or omit the mapping")
    return df.rename(columns={v: k for k, v in colmap.items()})


def _parse_include(col):
    v = col.astype("string").str.strip().str.lower().fillna("")
    bad = ~(v.isin(TRUE_VALUES | FALSE_VALUES) | (v == ""))
    return v, bad


def sanity_warnings(df):
    """Warn on Kelvin-looking temperatures and near-duplicate runs."""
    msgs = []
    T = df[["T_TA98", "T_WZ07_Fe2"]].to_numpy(float)
    if len(df) and np.median(T) > KELVIN_WARN_MEDIAN_C:
        msgs.append(f"median temperature {np.median(T):.0f} > {KELVIN_WARN_MEDIAN_C:.0f}: "
                    "inputs must be in deg C (Kelvin?)")
    lo, hi = T_PLAUSIBLE_C
    out = ((T < lo) | (T > hi)).any(axis=1)
    if out.any():
        msgs.append(f"temperature outside {lo:.0f}-{hi:.0f} deg C: {_rows(df, out)}")
    A = df[list(REQUIRED[1:])].to_numpy(float)
    for i in range(len(A)):
        for j in range(i + 1, len(A)):
            if np.all(np.abs(A[i] - A[j]) <= DUP_TOL_C):
                msgs.append(f"near-duplicate runs {df.SampleID.iloc[i]} and {df.SampleID.iloc[j]} "
                            f"(all temperatures within {DUP_TOL_C} deg C)")
    for m in msgs:
        warnings.warn(m, DataSanityWarning, stacklevel=3)
    return msgs


@dataclass
class FitResult:
    b0: float
    b1: float
    se_b0: float
    se_b1: float
    cov: list
    cov_scale: float
    chi2: float
    chi2_red: float
    dof: int
    n: int
    n_eff: float
    t_ref: float
    sigma0: float
    thresholds: dict = field(default_factory=dict)

    def predict(self, T):
        return self.b0 + self.b1 * (np.asarray(T, float) - self.t_ref)

    def band(self, T, z=Z95):
        """Pointwise confidence band of the fitted line (z * SE)."""
        dx = np.asarray(T, float) - self.t_ref
        C = np.asarray(self.cov)
        se = np.sqrt(C[0, 0] + 2 * dx * C[0, 1] + dx ** 2 * C[1, 1])
        y = self.predict(T)
        return y - z * se, y + z * se

    def equation(self, digits_b0=0, digits_b1=3):
        return (f"T_TA98 - T_WZ07 = {self.b0:.{digits_b0}f} (+/-{self.se_b0:.{digits_b0}f}) "
                f"{self.b1:+.{digits_b1}f} (+/-{self.se_b1:.{digits_b1}f}) (T_TA98 - {self.t_ref:.0f})")


# ------------------------------------------------------------------ input ----
def load_runs(path=DEFAULT_INPUT, colmap=None):
    """Read the run CSV, rename columns per ``colmap`` {canonical: csv_name},
    validate the schema and values, and return the included runs with
    y = T_TA98 - T_WZ07.  Raises InputError with the offending rows listed."""
    df = pd.read_csv(path, comment="#", dtype={"SampleID": str})
    if isinstance(colmap, (list, tuple)):
        colmap = parse_colmap(colmap)
    if colmap:
        df = _apply_colmap(df, colmap)
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        raise InputError(f"{path}: missing required column(s) {missing}")
    df = df.copy()
    blank_id = df["SampleID"].isna() | (df["SampleID"].astype(str).str.strip() == "")
    if blank_id.any():
        raise InputError(f"{path}: blank SampleID on line(s) {[int(i) + 2 for i in np.flatnonzero(blank_id)]}")
    df["SampleID"] = df["SampleID"].astype(str).str.strip()
    if df["SampleID"].duplicated().any():
        raise InputError(f"duplicate SampleID: {df.loc[df.SampleID.duplicated(), 'SampleID'].tolist()}")
    if "include" in df:
        v, bad = _parse_include(df["include"])
        if bad.any():
            raise InputError(f"unrecognised 'include' value(s) {sorted(set(df.loc[bad, 'include'].astype(str)))} "
                             f"in {_rows(df, bad)}; use blank/1/true/yes or 0/false/no")
        df = df[~v.isin(FALSE_VALUES)].copy()
    problems = []
    for c in REQUIRED[1:]:
        num = pd.to_numeric(df[c], errors="coerce")
        bad = ~np.isfinite(num.to_numpy(float))
        if bad.any():
            problems.append(f"{c}: blank/NaN/non-numeric in {_rows(df, bad)}")
        df[c] = num.astype(float)
    if "sigma_analytical_C" in df:
        raw = df["sigma_analytical_C"]
        blank = raw.isna() | (raw.astype("string").str.strip().fillna("") == "")
        num = pd.to_numeric(raw.where(~blank), errors="coerce")
        nonnum = (~blank) & num.isna()
        if nonnum.any():
            problems.append(f"sigma_analytical_C: non-numeric in {_rows(df, nonnum)}")
        badval = (~blank) & num.notna() & ~((num > 0) & np.isfinite(num))
        if badval.any():
            problems.append(f"sigma_analytical_C: must be finite and > 0 in {_rows(df, badval)}")
        df["sigma_analytical_C"] = num
    if problems:
        raise InputError(f"{path}: invalid values (nothing was fitted):\n  " + "\n  ".join(problems))
    if "Reference" not in df:
        df["Reference"] = ""
    df["Reference"] = df["Reference"].fillna("").astype(str)
    df["y_TA98_minus_WZ07_C"] = df["T_TA98"] - df["T_WZ07_Fe2"]
    df = df.reset_index(drop=True)
    sanity_warnings(df)
    return df


# ---------------------------------------------------------------- weights ----
def screen_excess(df, thresholds=SCREEN_THRESHOLDS_C):
    ex2 = np.zeros(len(df))
    for col, thr in thresholds.items():
        ex2 += np.maximum(0.0, np.abs(df[col].to_numpy(float)) - thr) ** 2
    return np.sqrt(ex2)


def analytical_sigma(df, sigma0=SIGMA0_C):
    """Per-run analytical sigma; blank/missing -> sigma0 (flagged as default)."""
    if "sigma_analytical_C" in df:
        s = pd.to_numeric(df["sigma_analytical_C"], errors="coerce").to_numpy(float)
    else:
        s = np.full(len(df), np.nan)
    is_default = ~np.isfinite(s)
    return np.where(is_default, sigma0, s), is_default


def compute_weights(df, thresholds=SCREEN_THRESHOLDS_C, sigma0=SIGMA0_C):
    """Return a DataFrame of per-run weight terms (same row order as df)."""
    for c in thresholds:
        bad = ~np.isfinite(df[c].to_numpy(float))
        if bad.any():
            raise InputError(f"{c}: non-finite screen value in {_rows(df, bad)}")
    if not (np.isfinite(sigma0) and sigma0 > 0):
        raise InputError(f"sigma0 must be finite and > 0 (got {sigma0})")
    ex = screen_excess(df, thresholds)
    sa, dflt = analytical_sigma(df, sigma0)
    sigma = np.sqrt(sa ** 2 + ex ** 2)
    w = 1.0 / sigma ** 2
    out = pd.DataFrame({
        "SampleID": df["SampleID"], "Reference": df["Reference"],
        "T_TA98_C": df["T_TA98"], "y_TA98_minus_WZ07_C": df["y_TA98_minus_WZ07_C"],
        **{c: df[c] for c in thresholds},
        "screens_failed": [", ".join(c[3:].replace("_", "-") for c, t in thresholds.items()
                                     if abs(r[c]) > t) for _, r in df.iterrows()],
        "screen_excess_C": ex, "sigma_analytical_C": sa,
        "sigma_analytical_is_default": dflt, "sigma_i_C": sigma,
        "w_raw_per_C2": w, "weight_frac": w / w.sum(), "weight_pct": 100 * w / w.sum(),
    })
    return out


def n_eff(w):
    w = np.asarray(w, float)
    return float(w.sum() ** 2 / (w ** 2).sum())


# -------------------------------------------------------------------- fit ----
def wls_line(T, y, w, t_ref=T_REF_C):
    T, y, w = (np.asarray(a, float) for a in (T, y, w))
    X = np.column_stack([np.ones_like(T), T - t_ref])
    A = X.T @ (X * w[:, None])
    cov0 = np.linalg.inv(A)
    b = cov0 @ (X.T @ (w * y))
    r = y - X @ b
    dof = len(y) - 2
    chi2 = float(np.sum(w * r ** 2))
    chi2_red = chi2 / dof if dof > 0 else float("nan")
    s = max(1.0, chi2_red) if math.isfinite(chi2_red) else 1.0
    return b, cov0 * s, s, chi2, chi2_red, dof


def fit(df, thresholds=SCREEN_THRESHOLDS_C, sigma0=SIGMA0_C, t_ref=T_REF_C):
    """Weights + WLS fit. Returns (weights DataFrame, FitResult)."""
    n = len(df)
    if n == 0:
        raise InputError("no runs to fit (empty input, or every run has include = 0)")
    if n < MIN_RUNS:
        raise InputError(f"only {n} run(s) included; at least {MIN_RUNS} are needed "
                         f"(with n <= 2 the line has no degrees of freedom and its SEs are meaningless)")
    if np.unique(df["T_TA98"].to_numpy(float)).size < 2:
        raise InputError("all included runs have the same T_TA98; the slope is undefined")
    wt = compute_weights(df, thresholds, sigma0)
    b, cov, s, chi2, chi2_red, dof = wls_line(df["T_TA98"], df["y_TA98_minus_WZ07_C"],
                                             wt["w_raw_per_C2"], t_ref)
    res = FitResult(b0=float(b[0]), b1=float(b[1]), se_b0=float(np.sqrt(cov[0, 0])),
                    se_b1=float(np.sqrt(cov[1, 1])), cov=cov.tolist(), cov_scale=float(s),
                    chi2=chi2, chi2_red=chi2_red, dof=dof, n=len(df),
                    n_eff=n_eff(wt["w_raw_per_C2"]), t_ref=t_ref, sigma0=sigma0,
                    thresholds=dict(thresholds))
    return wt, res


# ------------------------------------------------------- robustness checks ----
def ols(df, t_ref=T_REF_C):
    b, cov, *_ = wls_line(df["T_TA98"], df["y_TA98_minus_WZ07_C"], np.ones(len(df)), t_ref)
    # OLS: unit weights carry no absolute scale -> use the residual-variance
    # covariance (scale = chi2_red, not max(1, .)).
    T = df["T_TA98"].to_numpy(float); y = df["y_TA98_minus_WZ07_C"].to_numpy(float)
    X = np.column_stack([np.ones_like(T), T - t_ref])
    r = y - X @ b
    s2 = float(r @ r / (len(y) - 2))
    C = np.linalg.inv(X.T @ X) * s2
    return dict(b0=float(b[0]), b1=float(b[1]), se_b0=float(np.sqrt(C[0, 0])),
                se_b1=float(np.sqrt(C[1, 1])), resid_sd=math.sqrt(s2))


def leave_one_out(df, **kw):
    rows = []
    for i in range(len(df)):
        _, r = fit(df.drop(index=i).reset_index(drop=True), **kw)
        rows.append(dict(dropped=df.loc[i, "SampleID"], b0=r.b0, b1=r.b1, se_b1=r.se_b1, n_eff=r.n_eff))
    return pd.DataFrame(rows)


def bootstrap_slope(df, n=10000, seed=1, **kw):
    """Case-resampling bootstrap of the weighted fit (weights travel with runs)."""
    rng = np.random.default_rng(seed)
    wt = compute_weights(df, kw.get("thresholds", SCREEN_THRESHOLDS_C), kw.get("sigma0", SIGMA0_C))
    T = df["T_TA98"].to_numpy(float); y = df["y_TA98_minus_WZ07_C"].to_numpy(float)
    w = wt["w_raw_per_C2"].to_numpy(float)
    t_ref = kw.get("t_ref", T_REF_C)
    b1 = []
    for _ in range(n):
        i = rng.integers(0, len(y), len(y))
        if len(np.unique(T[i])) < 2:
            continue
        b1.append(wls_line(T[i], y[i], w[i], t_ref)[0][1])
    b1 = np.array(b1)
    return dict(n=len(b1), median=float(np.median(b1)), sd=float(b1.std(ddof=1)),
                ci95=[float(v) for v in np.percentile(b1, [2.5, 97.5])],
                frac_negative=float(np.mean(b1 < 0)))


# ----------------------------------------------------------------- output ----
def write_outputs(wt, res, outdir, stem="weighted_regression", extras=None):
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    wt.sort_values("w_raw_per_C2", ascending=False).to_csv(outdir / f"{stem}_weights.csv",
                                                           index=False, float_format="%.6g")
    summ = asdict(res); summ["equation"] = res.equation(); summ["version"] = __version__
    summ["se_convention"] = "cov = inv(X^T W X) * max(1, chi2_red); w = 1/sigma_i^2 (absolute)"
    summ["n_runs_default_sigma"] = int(wt["sigma_analytical_is_default"].sum())
    if extras:
        summ.update(extras)
    (outdir / f"{stem}_fit.json").write_text(json.dumps(summ, indent=2, default=float))
    lines = [res.equation(),
             f"n = {res.n}, n_eff = (sum w)^2/sum w^2 = {res.n_eff:.2f}",
             f"chi2 = {res.chi2:.2f}, dof = {res.dof}, chi2_red = {res.chi2_red:.3f}, cov scale = {res.cov_scale:.3f}",
             f"sigma0 = {res.sigma0} C ({summ['n_runs_default_sigma']} of {res.n} runs use the default/placeholder)",
             f"screens (C): {res.thresholds}"]
    (outdir / f"{stem}_fit.txt").write_text("\n".join(lines) + "\n")
    return lines


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0] +
                                " All temperatures in deg C.")
    p.add_argument("--input", default=str(DEFAULT_INPUT), help="run CSV (deg C)")
    p.add_argument("--outdir", default=str(REPO / "data" / "derived" / "weighted_regression"))
    p.add_argument("--t-ref", type=float, default=T_REF_C)
    p.add_argument("--sigma0", type=float, default=SIGMA0_C)
    p.add_argument("--thr-opx-cpx", type=float, default=SCREEN_THRESHOLDS_C["dT_opx_cpx"])
    p.add_argument("--thr-pyx-gt", type=float, default=SCREEN_THRESHOLDS_C["dT_pyx_gt"])
    p.add_argument("--thr-gt-ol", type=float, default=SCREEN_THRESHOLDS_C["dT_gt_ol"])
    p.add_argument("--col", action="append", default=[], metavar="NAME=CSVCOL",
                   help="map a canonical column name to a CSV column")
    p.add_argument("--robustness", action="store_true",
                   help="also report OLS, leave-one-out and bootstrap slope")
    p.add_argument("--n-boot", type=int, default=10000)
    a = p.parse_args(argv)
    thr = {"dT_opx_cpx": a.thr_opx_cpx, "dT_pyx_gt": a.thr_pyx_gt, "dT_gt_ol": a.thr_gt_ol}
    try:
        df = load_runs(a.input, parse_colmap(a.col))
        kw = dict(thresholds=thr, sigma0=a.sigma0, t_ref=a.t_ref)
        wt, res = fit(df, **kw)
    except InputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2)
    extras = None
    if a.robustness:
        loo = leave_one_out(df, **kw)
        extras = dict(ols=ols(df, a.t_ref),
                      loo_b1_min=float(loo.b1.min()), loo_b1_max=float(loo.b1.max()),
                      loo_b1_min_dropped=loo.loc[loo.b1.idxmin(), "dropped"],
                      loo_b1_max_dropped=loo.loc[loo.b1.idxmax(), "dropped"],
                      bootstrap_b1=bootstrap_slope(df, n=a.n_boot, **kw))
    for line in write_outputs(wt, res, a.outdir, extras=extras):
        print(line)
    if extras:
        loo.to_csv(Path(a.outdir) / "weighted_regression_loo.csv",
                                       index=False, float_format="%.6g")
        print(json.dumps(extras, indent=2, default=float))
    return res


if __name__ == "__main__":
    main()
