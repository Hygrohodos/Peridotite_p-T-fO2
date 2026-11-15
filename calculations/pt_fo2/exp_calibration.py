"""
exp_calibration.py
------------------
Weighted calibration of the TA98 − WZ07 discrepancy on ALL experiments with
measured garnet Fe3+/ΣFe (manuscript Section 3.3; reviewer item 1).

Replaces the single-experiment (F1300o) flat −118 °C offset by a weighted
best-fit curve through every usable run product in
``data/derived/exp_TA98_WZ07_Fe3.csv`` (Matjuschkin et al. 2014, n = 8;
Rzehak et al. 2020, n = 5), written by ``scripts/plot_exp_ta98_wz07.py``.

y = TA98 − WZ07(Fe2+ only, measured Fe3+ removed)  [°C]

Disequilibrium weight (primary, "inverse variance")
    The experiments are scored with the same three equilibrium screens as
    the xenoliths (Section 3.2; thresholds in params.ScreeningParams):
        opx–cpx  |TA98 − NG10|            < 80 °C
        pyx–gt   |HA84 − quadratic(TA98)| < 60 °C
        gt–ol    |WZ07pred − WZ07|        < 50 °C
    The amount by which a run exceeds each threshold is a temperature
    misfit between thermometers that should agree at equilibrium, so it is
    treated as an extra 1σ temperature uncertainty on y:
        σ_i² = σ0² + Σ_j max(0, |dT_ij| − thr_j)²
        w_i  = 1 / σ_i²
    σ0 (default 50 °C, the gt–ol acceptance width) is the scatter expected
    for an equilibrated run; the reduced χ² of the fit is reported so the
    choice can be checked (≈ 1 means σ_i are consistent with the scatter).

Alternatives for sensitivity checks: ``ratio`` (w = 1/D², D = largest
|dT_j|/thr_j), ``unweighted``, and ``quality`` (inverse variance ×
oxide-total quality factor).

Model: weighted least squares, y = b0 + b1 (T − T_ref) (``deg`` 0/1/2);
uncertainty from the parameter covariance (scaled by max(1, χ²_red)) and
from a case-resampling bootstrap.

Pipeline mapping (fe3_correction.wz07_offset):
    offset(T) = −y(T) = wz07_offset_c + wz07_offset_slope · T
"""

from __future__ import annotations

import numpy as np
import pandas as pd

SCREEN_THR = {"dT_opx_cpx": 80.0, "dT_pyx_gt": 60.0, "dT_gt_ol": 50.0}
PHASES = ("cpx", "opx", "garnet", "olivine")
T_REF = 1350.0


def screen_excess(e):
    """Per-run quadrature sum of the °C by which each screen is exceeded."""
    ex2 = sum(np.maximum(0.0, e[c].abs() - t) ** 2 for c, t in SCREEN_THR.items())
    return np.sqrt(ex2.to_numpy(dtype=float))


def screen_ratio(e):
    return np.max([e[c].abs().to_numpy(float) / t for c, t in SCREEN_THR.items()], axis=0)


def quality_factor(e, raw, window=1.5):
    """Analytical-quality factor from oxide totals (1 = all totals within
    100 ± ``window`` wt%; decreases as 1 / (1 + Σ excess wt%))."""
    tot = raw.set_index(raw["Sample ID"].astype(str).str.replace(" ", ""))
    q = []
    for sid in e["SampleID"].astype(str).str.replace(" ", ""):
        ex = 0.0
        for ph in PHASES:
            v = pd.to_numeric(tot.loc[sid, f"{ph}_Total"], errors="coerce")
            v = float(np.atleast_1d(v)[0])
            if np.isfinite(v):
                ex += max(0.0, abs(v - 100.0) - window)
        q.append(1.0 / (1.0 + ex))
    return np.asarray(q)


def weights(e, scheme="invvar", sigma0=50.0, raw=None):
    if scheme == "invvar":
        return 1.0 / (sigma0 ** 2 + screen_excess(e) ** 2)
    if scheme == "ratio":
        return 1.0 / screen_ratio(e) ** 2
    if scheme == "unweighted":
        return np.ones(len(e))
    if scheme == "quality":
        return weights(e, "invvar", sigma0) * quality_factor(e, raw)
    raise ValueError(scheme)


def _design(T, deg):
    return np.vander(np.asarray(T, float) - T_REF, deg + 1, increasing=True)


def wls(T, y, w, deg=1):
    X = _design(T, deg)
    W = np.asarray(w, float)
    A = X.T @ (X * W[:, None])
    cov = np.linalg.inv(A)
    b = cov @ (X.T @ (W * y))
    r = y - X @ b
    dof = len(y) - (deg + 1)
    chi2 = float(np.sum(W * r ** 2))
    k = deg + 1
    return dict(beta=b, cov=cov, chi2=chi2, dof=dof,
                chi2_red=chi2 / dof if dof > 0 else np.nan,
                aicc=chi2 + 2 * k + 2 * k * (k + 1) / max(len(y) - k - 1, 1),
                deg=deg)


def predict(fit, T):
    return _design(T, fit["deg"]) @ fit["beta"]


def band_cov(fit, T, z=1.96):
    X = _design(T, fit["deg"])
    s = max(1.0, fit["chi2_red"]) if np.isfinite(fit["chi2_red"]) else 1.0
    se = np.sqrt(np.einsum("ij,jk,ik->i", X, fit["cov"] * s, X))
    return predict(fit, T) - z * se, predict(fit, T) + z * se


def bootstrap(T, y, w, deg=1, n=4000, seed=1):
    rng = np.random.default_rng(seed)
    T, y, w = map(lambda a: np.asarray(a, float), (T, y, w))
    out = []
    for _ in range(n):
        i = rng.integers(0, len(y), len(y))
        if len(np.unique(T[i])) <= deg:
            continue
        out.append(wls(T[i], y[i], w[i], deg)["beta"])
    return np.array(out)


def band_boot(betas, T, deg, q=(2.5, 97.5)):
    Y = _design(T, deg) @ betas.T
    return np.percentile(Y, q[0], axis=1), np.percentile(Y, q[1], axis=1)


def to_offset_params(fit):
    """Linear fit → (wz07_offset_c, wz07_offset_slope) for params.Fe3Params."""
    b = fit["beta"]
    if fit["deg"] == 0:
        return float(-b[0]), 0.0
    if fit["deg"] == 1:
        return float(-(b[0] - b[1] * T_REF)), float(-b[1])
    raise ValueError("pipeline offset supports deg ≤ 1")
