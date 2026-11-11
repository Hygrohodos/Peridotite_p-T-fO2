"""Revision v2 (option B): the WZ07 match level is anchored to the natural Slave calibration mean."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "scripts")]
import make_figures as mf  # noqa: E402
import weighted_thermometer_regression as wtr  # noqa: E402
from calculations.pt_fo2.params import PARAMS  # noqa: E402
from calculations.pt_fo2 import level_anchor as la  # noqa: E402
from calculations.pt_fo2 import calibration as calib  # noqa: E402


@pytest.fixture(scope="module")
def df():
    d = pd.read_csv(mf.FULL, encoding="utf-8-sig", low_memory=False)
    return d.rename(columns={d.columns[0]: "SampleID"})


def test_params_level_is_the_anchor(df):
    a = la.anchor_level(df)
    assert a["n"] == 25
    assert PARAMS.fe3.dt_match_c == pytest.approx(a["dt_level_C"], abs=0.05)
    assert a["mean_D_C"] == pytest.approx(-91.89, abs=0.05)


def test_anchor_uses_no_validation_samples(df):
    m = la.anchor_set(df)
    assert not (m & calib.validation_mask(df)).any()
    assert not df.loc[m, "SampleID"].astype(str).isin(la.BLUE_IDS + la.RED_IDS).any()


def test_error_groups_unchanged_at_adopted_level(df):
    _, _, blue, red, rem, _ = calib.error_groups(df)
    assert sorted(df.SampleID[blue].astype(str)) == sorted(la.BLUE_IDS)
    assert sorted(df.SampleID[red].astype(str)) == sorted(la.RED_IDS)
    assert (rem == la.anchor_set(df)).all()


def test_pipeline_curve_is_the_standalone_fit():
    _, f = wtr.fit(wtr.load_runs())
    T = np.array([1100.0, 1162.03, 1250.0, 1350.0, 1500.0, 1619.8, 1700.0])
    exp = f.predict(np.clip(T, PARAMS.fe3.wz07_offset_tmin_c, PARAMS.fe3.wz07_offset_tmax_c))
    assert np.allclose(la.curve(T), exp, atol=0.01)
    assert la.curve([900.0])[0] == pytest.approx(la.curve([1162.03])[0])  # held flat below the experiments


def test_anchor_reproduces_mean_measured_fe3(df):
    """At the anchored level the anchor set's mean est - meas is ~0 (definition of the anchor)."""
    m = la.anchor_set(df)
    est = pd.to_numeric(df["XFe3_gt_iter"], errors="coerce")
    meas = calib.measured_fe3(df)
    assert abs(float((est - meas)[m].mean())) < 0.003
    assert not PARAMS.fe3.apply_empirical_correction
