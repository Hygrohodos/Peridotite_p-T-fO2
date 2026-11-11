"""Reproduces the revision-v1 Fig. 4 fit from the 13 current experiments.
Run:  python3 -m pytest tests/test_weighted_regression.py   (or python3 tests/test_weighted_regression.py)
Needs pytest.
"""
import subprocess
import sys
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

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


# ------------------------------------------------ validation (review of 34c7e65) ----
SCRIPT = Path(wtr.__file__)


def _csv(tmp_path, edit=None, name="runs.csv"):
    df = pd.read_csv(wtr.DEFAULT_INPUT, comment="#")
    if edit:
        r = edit(df)
        df = df if r is None else r
    f = tmp_path / name
    df.to_csv(f, index=False)
    return f


def _cli(*args):
    return subprocess.run([sys.executable, str(SCRIPT), *map(str, args)], capture_output=True, text=True)


@pytest.mark.parametrize("col", ["T_TA98", "T_WZ07_Fe2", "dT_opx_cpx", "dT_pyx_gt", "dT_gt_ol"])
def test_defect1_blank_required_value_raises_and_names_row(tmp_path, col):
    def ed(df):
        df[col] = df[col].astype(object); df.loc[3, col] = ""
    f = _csv(tmp_path, ed)
    with pytest.raises(wtr.InputError, match=rf"{col}.*F1200r \(line 5\)"):
        wtr.load_runs(f)
    r = _cli("--input", f, "--outdir", tmp_path / "o")
    assert r.returncode == 2 and "F1200r" in r.stderr


def test_defect1_nan_and_text_values_raise(tmp_path):
    def ed(df):
        df["dT_gt_ol"] = df["dT_gt_ol"].astype(object)
        df.loc[0, "dT_gt_ol"] = "nan"; df.loc[1, "dT_gt_ol"] = "abc"
    with pytest.raises(wtr.InputError, match=r"F1100o.*F1100r"):
        wtr.load_runs(_csv(tmp_path, ed))


@pytest.mark.parametrize("val", ["", "1", "true", "Yes", "Y"])
def test_defect2_blank_or_true_include_is_included(tmp_path, val):
    def ed(df):
        df["include"] = "1"; df.loc[4, "include"] = val
    assert len(wtr.load_runs(_csv(tmp_path, ed))) == 13


@pytest.mark.parametrize("val", ["0", "false", "no", "N"])
def test_defect2_false_include_excludes(tmp_path, val):
    def ed(df):
        df["include"] = ""; df.loc[4, "include"] = val
    d = wtr.load_runs(_csv(tmp_path, ed))
    assert len(d) == 12 and "F1300o" not in set(d.SampleID)


@pytest.mark.parametrize("val", ["x", "included", "2", "nan?"])
def test_defect2_unrecognised_include_raises(tmp_path, val):
    def ed(df):
        df["include"] = "1"; df.loc[4, "include"] = val
    with pytest.raises(wtr.InputError, match=r"include.*F1300o"):
        wtr.load_runs(_csv(tmp_path, ed))


@pytest.mark.parametrize("val,msg", [("abc", "non-numeric"), ("0", "> 0"), ("-5", "> 0"), ("inf", "> 0")])
def test_defect5_bad_sigma_raises(tmp_path, val, msg):
    def ed(df):
        df["sigma_analytical_C"] = ""; df.loc[2, "sigma_analytical_C"] = val
    with pytest.raises(wtr.InputError, match=rf"sigma_analytical_C.*{msg}.*F1200o"):
        wtr.load_runs(_csv(tmp_path, ed))


def test_defect5_blank_sigma_uses_default(tmp_path):
    def ed(df):
        df["sigma_analytical_C"] = ""; df.loc[2, "sigma_analytical_C"] = "70"
    wt, _ = wtr.fit(wtr.load_runs(_csv(tmp_path, ed)))
    assert wt.sigma_analytical_is_default.sum() == 12
    assert wt.set_index("SampleID").loc["F1200o", "sigma_analytical_C"] == 70


@pytest.mark.parametrize("n", [0, 1, 2])
def test_defect6_too_few_runs_raises(tmp_path, n):
    def ed(df):
        df["include"] = 0; df.loc[: n - 1, "include"] = 1
    d = wtr.load_runs(_csv(tmp_path, ed))
    assert len(d) == n
    with pytest.raises(wtr.InputError, match="no runs" if n == 0 else "at least 3"):
        wtr.fit(d)
    assert _cli("--input", _csv(tmp_path, ed, "b.csv"), "--outdir", tmp_path / "o").returncode == 2


def test_defect6_three_runs_ok(tmp_path):
    d = wtr.load_runs(_csv(tmp_path, lambda df: df.head(3)))
    assert wtr.fit(d)[1].dof == 1


def test_defect7_screens_failed_matches_excess():
    d = pd.DataFrame({"SampleID": ["a", "b", "c"], "Reference": "", "T_TA98": [1200.0, 1300, 1400],
                      "y_TA98_minus_WZ07_C": [0.0, 1, 2],
                      "dT_opx_cpx": [80.0, 80.5, 0], "dT_pyx_gt": [-60.0, 0, 0], "dT_gt_ol": [50.0, 0, 0]})
    wt = wtr.compute_weights(d)
    assert wt.screens_failed.tolist() == ["", "opx-cpx", ""]      # exactly at threshold = pass
    assert ((wt.screen_excess_C > 0) == (wt.screens_failed != "")).all()
    _, wt0, _ = _current()
    assert ((wt0.screen_excess_C > 0) == (wt0.screens_failed != "")).all()
    d.loc[2, "dT_gt_ol"] = np.nan
    with pytest.raises(wtr.InputError, match="dT_gt_ol.*c"):     # NaN cannot pass
        wtr.compute_weights(d)


def test_defect8_colmap_ok_and_errors(tmp_path):
    f = _csv(tmp_path, lambda df: df.rename(columns={"SampleID": "Run"}))
    assert len(wtr.load_runs(f, wtr.parse_colmap(["SampleID=Run"]))) == 13
    for bad in (["BAD"], ["a=b=c"], ["=Run"], ["SampleID="], ["Nonsense=Run"],
                ["SampleID=Run", "SampleID=Run2"], ["SampleID=Run", "Reference=Run"]):
        with pytest.raises(wtr.InputError):
            wtr.parse_colmap(bad)
    with pytest.raises(wtr.InputError, match="no column"):
        wtr.load_runs(f, {"SampleID": "Missing"})
    g = _csv(tmp_path, lambda df: df.assign(Run=df.SampleID), "g.csv")   # SampleID AND Run present
    with pytest.raises(wtr.InputError, match="already has a column 'SampleID'"):
        wtr.load_runs(g, {"SampleID": "Run"})
    for args in (["--col", "BAD"], ["--col", "SampleID=Run", "--input", g]):
        r = _cli(*args, "--outdir", tmp_path / "o")
        assert r.returncode == 2 and r.stderr.startswith("error:")


def test_defect9_sanity_warnings(tmp_path):
    with warnings.catch_warnings():
        warnings.simplefilter("error", wtr.DataSanityWarning)
        wtr.load_runs()                                            # current data: no warning
    def kelvin(df):
        df["T_TA98"] += 273.15; df["T_WZ07_Fe2"] += 273.15
    with pytest.warns(wtr.DataSanityWarning) as rec:
        wtr.load_runs(_csv(tmp_path, kelvin))
    msgs = [str(r.message) for r in rec]
    assert any("Kelvin" in m for m in msgs) and any("outside 0-2000" in m for m in msgs)
    def dup(df):
        r = df.iloc[[4]].copy(); r["SampleID"] = "F1300o_copy"; r["T_TA98"] += 0.1
        return pd.concat([df, r])
    with pytest.warns(wtr.DataSanityWarning, match="near-duplicate runs F1300o and F1300o_copy"):
        wtr.load_runs(_csv(tmp_path, dup, "d.csv"))


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))

