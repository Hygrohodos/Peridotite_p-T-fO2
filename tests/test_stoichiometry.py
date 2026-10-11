"""Charge-balance garnet Fe3+ against Droop (1987) Table 3.

Run with ``python -m pytest tests`` or ``python tests/test_stoichiometry.py``.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from calculations.pt_fo2.stoichiometry import calculate_stoichiometric_fe3  # noqa: E402

# Droop (1987) Mineral Mag 51:434, Table 3, garnet analysis A
# (ferrian uvarovite, Edessa).  ZnO 0.19 wt% is not one of the pipeline
# oxides, so it is passed as NiO (same charge and cation/oxygen numbers;
# the molar mass difference changes S by < 0.001).
DROOP_GARNET = dict(
    garnet_SiO2=35.94, garnet_TiO2=0.12, garnet_Al2O3=0.28,
    garnet_Cr2O3=18.11, garnet_FeO=11.37, garnet_MnO=0.03,
    garnet_MgO=0.00, garnet_CaO=33.33, garnet_NiO=0.19 * 74.692 / 81.379,
)


def _calc(rows):
    return calculate_stoichiometric_fe3(pd.DataFrame(rows), verbose=False)


def test_droop_table3_garnet():
    r = _calc([DROOP_GARNET]).iloc[0]
    assert abs(r["cation_sum_12O"] - 8.268) < 0.001      # S, Table 3 A
    assert abs(r["Fe3_droop"] - 0.779) < 0.001           # Fe3+, Table 3 B
    assert abs(r["Fe2_droop"] - 0.013) < 0.001           # Fe2+, Table 3 B
    assert abs(r["XFe3_droop"] - 0.779 / 0.792) < 0.002


def test_droop_equals_quinn():
    r = _calc([DROOP_GARNET]).iloc[0]
    assert abs(r["Fe3_droop"] - r["Fe3_quinn"]) < 1e-9
    assert abs(r["Fe2_droop"] - r["Fe2_quinn"]) < 1e-9
    assert abs(r["XFe3_droop"] - r["XFe3_quinn"]) < 1e-12


def test_no_fe3_when_cations_le_8_and_cap_at_total_fe():
    # Ideal pyrope-almandine: S = 8 exactly -> no Fe3+.
    ideal = dict(garnet_SiO2=3 * 60.083, garnet_Al2O3=101.961,
                 garnet_MgO=1.5 * 40.304, garnet_FeO=1.5 * 71.844)
    # Fe-poor analysis with a large cation excess -> Fe3+ capped at ΣFe.
    capped = dict(garnet_SiO2=38.0, garnet_Al2O3=21.0, garnet_MgO=36.0,
                  garnet_FeO=0.5, garnet_CaO=4.0)
    r = _calc([ideal, capped])
    assert abs(r.loc[0, "XFe3_droop"]) < 1e-9 and abs(r.loc[0, "XFe3_quinn"]) < 1e-9
    assert r.loc[1, "XFe3_droop"] == 1.0 and r.loc[1, "XFe3_quinn"] == 1.0
    assert r.loc[1, "Fe2_droop"] == 0.0 and r.loc[1, "Fe2_quinn"] == 0.0


def test_no_garnet_analysis_is_nan():
    r = _calc([dict(garnet_SiO2=np.nan, garnet_FeO=np.nan)]).iloc[0]
    for c in ("Fe3_droop", "Fe2_droop", "XFe3_droop",
              "Fe3_quinn", "Fe2_quinn", "XFe3_quinn"):
        assert np.isnan(r[c]), c


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
