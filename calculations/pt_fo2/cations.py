"""
cations.py
----------
Convert mineral oxide weight-percent analyses to cations (atoms per formula unit).

This is Stage 1 of the oxythermobarometric pipeline.
"""

import pandas as pd
import numpy as np

# Molecular weights and stoichiometry used in the original SD1 workbook
OXIDE_INFO = {
    "SiO2":  {"mw": 60.083,  "n_cat": 1, "n_O": 2, "cation": "Si"},
    "TiO2":  {"mw": 79.865,  "n_cat": 1, "n_O": 2, "cation": "Ti"},
    "Al2O3": {"mw": 101.961, "n_cat": 2, "n_O": 3, "cation": "Al"},
    "Cr2O3": {"mw": 151.989, "n_cat": 2, "n_O": 3, "cation": "Cr"},
    "FeO":   {"mw": 71.844,  "n_cat": 1, "n_O": 1, "cation": "Fe"},
    "MnO":   {"mw": 70.937,  "n_cat": 1, "n_O": 1, "cation": "Mn"},
    "MgO":   {"mw": 40.304,  "n_cat": 1, "n_O": 1, "cation": "Mg"},
    "CaO":   {"mw": 56.077,  "n_cat": 1, "n_O": 1, "cation": "Ca"},
    "Na2O":  {"mw": 61.97,   "n_cat": 2, "n_O": 1, "cation": "Na"},
    "NiO":   {"mw": 74.692,  "n_cat": 1, "n_O": 1, "cation": "Ni"},
}

# Which oxide columns belong to each mineral and the ideal oxygen basis
MINERAL_SETTINGS = {
    "cpx": {
        "oxide_cols": [
            "cpx_SiO2", "cpx_TiO2", "cpx_Al2O3", "cpx_Cr2O3", "cpx_FeO",
            "cpx_MnO", "cpx_MgO", "cpx_CaO", "cpx_Na2O", "cpx_NiO",
        ],
        "ideal_oxygens": 6,
        "prefix": "cpx_",
    },
    "opx": {
        "oxide_cols": [
            "opx_SiO2", "opx_TiO2", "opx_Al2O3", "opx_Cr2O3", "opx_FeO",
            "opx_MnO", "opx_MgO", "opx_CaO", "opx_Na2O", "opx_NiO",
        ],
        "ideal_oxygens": 6,
        "prefix": "opx_",
    },
    "garnet": {
        "oxide_cols": [
            "garnet_SiO2", "garnet_TiO2", "garnet_Al2O3", "garnet_Cr2O3",
            "garnet_FeO", "garnet_MnO", "garnet_MgO", "garnet_CaO",
            "garnet_Na2O", "garnet_NiO",
        ],
        "ideal_oxygens": 12,
        "prefix": "garnet_",
    },
    "olivine": {
        "oxide_cols": [
            "olivine_SiO2", "olivine_TiO2", "olivine_Al2O3", "olivine_Cr2O3",
            "olivine_FeO", "olivine_MnO", "olivine_MgO", "olivine_CaO",
            "olivine_Na2O", "olivine_NiO",
        ],
        "ideal_oxygens": 4,
        "prefix": "olivine_",
    },
}


def oxides_to_apfu(df, oxide_cols, ideal_oxygens=None, round_digits=4):
    """
    Convert a set of oxide columns (wt%) into cations on an ideal-oxygen basis.

    Parameters
    ----------
    df : pandas.DataFrame
        Input data containing the oxide columns.
    oxide_cols : list of str
        Column names of the oxides for one mineral.
    ideal_oxygens : float or None
        Number of oxygens to normalise to (6 for pyroxenes, 12 for garnet, 4 for olivine).
    round_digits : int
        Decimal places for the final cation values.

    Returns
    -------
    pandas.DataFrame
        Columns for each cation + Sum_Cations (+ Norm_O if requested).
    """
    if not oxide_cols:
        return pd.DataFrame(index=df.index)

    cat_df = pd.DataFrame(index=df.index)
    total_O = pd.Series(0.0, index=df.index, dtype=float)
    cat_moles = {}

    for full_col in oxide_cols:
        base_ox = full_col.split("_")[-1]
        if base_ox not in OXIDE_INFO or full_col not in df.columns:
            continue
        info = OXIDE_INFO[base_ox]
        wt = pd.to_numeric(df[full_col], errors="coerce").fillna(0.0)
        moles_ox = wt / info["mw"]
        cat_name = info["cation"]
        if cat_name not in cat_moles:
            cat_moles[cat_name] = pd.Series(0.0, index=df.index, dtype=float)
        cat_moles[cat_name] += moles_ox * info["n_cat"]
        total_O += moles_ox * info["n_O"]

    if not cat_moles:
        return pd.DataFrame(index=df.index)

    if ideal_oxygens is not None and total_O.sum() > 0:
        scale = ideal_oxygens / total_O.replace(0, np.nan)
    else:
        scale = pd.Series(1.0, index=df.index)

    for cat_name, moles in cat_moles.items():
        cat_df[cat_name] = (moles * scale).round(round_digits)
    cat_df["Sum_Cations"] = cat_df.sum(axis=1).round(round_digits)
    if ideal_oxygens is not None:
        cat_df["Norm_O"] = ideal_oxygens
    return cat_df


def convert_oxides_to_cations(df, verbose=True):
    """
    Run the full oxide-to-cation conversion for cpx, opx, garnet and olivine.

    Parameters
    ----------
    df : pandas.DataFrame
        Raw geochemical data (oxide wt% columns).
    verbose : bool
        Print progress messages.

    Returns
    -------
    pandas.DataFrame
        Original columns plus all new cation columns (prefixed).
    """
    result = df.copy()

    if verbose:
        print("Converting oxides to cations …")

    for name, config in MINERAL_SETTINGS.items():
        apfu = oxides_to_apfu(
            df,
            config["oxide_cols"],
            ideal_oxygens=config["ideal_oxygens"],
        )
        if apfu.empty:
            if verbose:
                print(f"  Warning: no columns found for {name}")
            continue
        apfu = apfu.add_prefix(config["prefix"])
        result = pd.concat([result, apfu], axis=1)

    if verbose:
        print(f"  Total samples processed: {len(result)}")

    return result
