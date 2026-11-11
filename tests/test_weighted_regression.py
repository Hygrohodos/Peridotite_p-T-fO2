"""Reproduces the revision-v1 Fig. 4 fit from the 13 current experiments.
Run:  python -m pytest tests/test_weighted_regression.py   (or python tests/test_weighted_regression.py)
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import weighted_thermometer_regression as wtr


def _current():
    df = wtr.load_runs()  # data/derived/exp_TA98_WZ07_Fe3.csv
    return df, *wtr.fit(df)


def test_reproduces_current_fit():
    df, wt, r = _current()
    assert r.n == 13
    assert round(r.b0) == -14 and round(r.se_b0) == 27
    assert round(r.b1, 3) == -0.555 and round(r.se_b1, 3) == 0.216
    assert round(r.n_eff, 1) == 8.4
    assert r.cov_scale == 1.0 and abs(r.chi2_red - 0.839) < 1e-3


def test_weights_match_figure_table():
    _, wt, _ = _current()
    pct = wt.set_index("SampleID")["weight_pct"]
    assert round(pct["F1300o"], 1) == 24.0 and round(pct["MaLm1"], 1) == 0.4
    assert abs(wt["weight_frac"].sum() - 1) < 1e-12
    # sigma_i^2 = 50^2 + excess^2 for every run
    assert np.allclose(wt.sigma_i_C ** 2, 50 ** 2 + wt.screen_excess_C ** 2)


def test_excess_formula_by_hand():
    d = pd.DataFrame({"dT_opx_cpx": [100.0, -10], "dT_pyx_gt": [-90.0, 0], "dT_gt_ol": [50.0, 0]})
    assert np.allclose(wtr.screen_excess(d), [np.hypot(20, 30), 0.0])


def test_explicit_sigma_and_include_and_centring(tmp_path):
    df = pd.read_csv(wtr.DEFAULT_INPUT, comment="#")
    df["sigma_analytical_C"] = 50.0
    df["include"] = 1
    f = tmp_path / "runs.csv"; df.to_csv(f, index=False)
    _, r = wtr.fit(wtr.load_runs(f))
    _, _, r0 = _current()
    assert np.isclose(r.b1, r0.b1) and np.isclose(r.b0, r0.b0)
    _, r2 = wtr.fit(wtr.load_runs(f), t_ref=0.0)   # centring changes b0, not slope
    assert np.isclose(r2.b1, r0.b1) and np.isclose(r2.predict(1350), r0.b0)
    df.loc[0, "include"] = 0; df.to_csv(f, index=False)
    assert wtr.fit(wtr.load_runs(f))[1].n == 12


if __name__ == "__main__":
    import tempfile
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn(Path(tempfile.mkdtemp())) if "tmp_path" in fn.__code__.co_varnames else fn()
            print("ok", name)
