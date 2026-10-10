"""
pipeline.py
-----------
High-level orchestration of the full oxythermobarometric calculation.

This is the single entry point that the website (or any other script)
should call.  It chains:

  1. Load & clean
  2. Oxide → cation conversion
  3. Iterative CA91 P + TA98 T  (+ supporting thermometers)
  4. Equilibrium screening
  5. Garnet Fe³⁺ iteration + empirical correction
  6. (Future) Oxygen fugacity

Usage example
-------------
from calculations.pipeline import run_full_pipeline
result, stats = run_full_pipeline("data/xenolith_geochem.csv")
"""

from pathlib import Path
import pandas as pd
import numpy as np

from .cations import convert_oxides_to_cations
from .thermobarometry import iterate_PT, calculate_all_temperatures
from .screening import apply_equilibrium_screening
from .fe3_correction import iterate_garnet_fe3, apply_empirical_fe3_correction
from .fo2 import calculate_fo2
from .params import PARAMS
from .stoichiometry import calculate_stoichiometric_fe3


def load_and_clean(csv_path, verbose=True):
    """
    Read the geochemical CSV and perform the same column cleaning that
    was used in the original monolithic script.
    """
    csv_path = Path(csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(f"Cannot find input file: {csv_path}")

    df = pd.read_csv(
        csv_path,
        encoding="cp1252",
        encoding_errors="ignore",
        on_bad_lines="skip",
    )

    # Clean column names (remove non-breaking spaces, extra whitespace)
    df.columns = (
        df.columns
        .str.strip()
        .str.replace("\xa0", "", regex=False)
        .str.replace(r"\s+", "", regex=True)
    )

    # Ensure a consistent SampleID column
    if "SampleID" not in df.columns:
        for c in df.columns:
            if "sample" in c.lower() and "id" in c.lower():
                df = df.rename(columns={c: "SampleID"})
                break

    # Force a few known columns to numeric
    cols_to_numeric = [
        "cpx_TiO2", "opx_TiO2", "opx_MnO", "opx_Na2O", "opx_NiO",
        "garnet_TiO2", "garnet_Fe3+", "olivine_MnO",
    ]
    for c in cols_to_numeric:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    df = fill_craton_from_locations(df, csv_path, verbose=verbose)

    if verbose:
        print(f"Loaded {len(df)} samples from {csv_path.name}")

    return df


def _loc_key(val):
    return (
        str(val).strip().lower()
        .replace("_", " ")
        .replace("-", " ")
        .replace("/", " ")
        .replace("  ", " ")
    )


def _location_tables(csv_path):
    """Resolve data/locations.csv and data/master/kimberlite_locations.csv."""
    root = Path(csv_path).resolve()
    for _ in range(4):
        if (root / "data").is_dir() or (root / "calculations").is_dir():
            break
        root = root.parent
    if not (root / "data").is_dir():
        root = Path(csv_path).resolve().parent
        while root != root.parent and not (root / "data").is_dir():
            root = root.parent
    cands = [
        root / "data" / "locations.csv",
        root / "data" / "master" / "kimberlite_locations.csv",
        root / "data" / "master" / "locations.csv",
    ]
    tables = []
    for p in cands:
        if p.is_file():
            tables.append(pd.read_csv(p, encoding="utf-8-sig"))
    return tables


def fill_craton_from_locations(df, csv_path=None, verbose=True):
    """Fill blank Craton using locations.csv / kimberlite_locations.csv."""
    if "Location" not in df.columns:
        return df
    result = df.copy()
    if "Craton" not in result.columns:
        result["Craton"] = ""
    blank = result["Craton"].isna() | (result["Craton"].astype(str).str.strip() == "")

    lookup = {}
    tables = _location_tables(csv_path or ".")
    for tab in tables:
        loc_col = next((c for c in tab.columns if c.lower() == "location"), None)
        cr_col = next((c for c in tab.columns if c.lower() == "craton"), None)
        if loc_col is None or cr_col is None:
            continue
        for _, row in tab.iterrows():
            cr = str(row[cr_col]).strip() if pd.notna(row[cr_col]) else ""
            if cr and cr.lower() not in {"nan", "none"}:
                lookup[_loc_key(row[loc_col])] = cr

    # Notes / country hints from data/locations.csv (first table if present)
    notes_map = {}
    country_map = {}
    if tables:
        tab = tables[0]
        loc_col = next((c for c in tab.columns if c.lower() == "location"), None)
        notes_col = next((c for c in tab.columns if c.lower() == "notes"), None)
        ctry_col = next((c for c in tab.columns if c.lower() == "country"), None)
        if loc_col is not None:
            for _, row in tab.iterrows():
                k = _loc_key(row[loc_col])
                if notes_col is not None and pd.notna(row.get(notes_col)):
                    notes_map[k] = str(row[notes_col]).lower()
                if ctry_col is not None and pd.notna(row.get(ctry_col)):
                    country_map[k] = str(row[ctry_col]).strip()

    locality_rules = {
        "gahcho kué": "Slave Craton",
        "gahcho kue": "Slave Craton",
        "obnazhennaya": "Siberian Craton",
        "lashaine tanzania": "Tanzania Craton",
        "labait tanzania": "Tanzania Craton",
        "nikos somerset islnd": "Rae Craton",
        "batty bay somerset island": "Rae Craton",
        "nanorluk somerset island": "Rae Craton",
        "ham somerset islnd": "Rae Craton",
        "amayersuk somerset islnd": "Rae Craton",
        "elwin bay somerset islnd": "Rae Craton",
        "jp somerset islnd": "Rae Craton",
        "kuopio kaavi pipe 14": "Karelian Craton",
        "kuopio kaavi pipe 9": "Karelian Craton",
        "lahtojoki pipe 7": "Karelian Craton",
    }
    note_rules = (
        ("slave", "Slave Craton"),
        ("lac de gras", "Slave Craton"),
        ("ekati", "Slave Craton"),
        ("kaapvaal", "Kaapvaal Craton"),
        ("sakha", "Siberian Craton"),
        ("ontario", "Superior Craton"),
        ("free state", "Kaapvaal Craton"),
        ("northern cape", "Kaapvaal Craton"),
        ("limpopo", "Kaapvaal Craton"),
        ("gauteng", "Kaapvaal Craton"),
        ("maloti", "Kaapvaal Craton"),
        ("lesotho", "Kaapvaal Craton"),
        ("botswana", "Kaapvaal Craton"),
    )
    country_rules = {
        "south africa": "Kaapvaal Craton",
        "lesotho": "Kaapvaal Craton",
        "botswana": "Kaapvaal Craton",
        "eswatini": "Kaapvaal Craton",
        "swaziland": "Kaapvaal Craton",
    }

    filled = 0
    for i in result.index[blank]:
        k = _loc_key(result.at[i, "Location"])
        cr = lookup.get(k, "") or locality_rules.get(k, "")
        if not cr:
            notes = notes_map.get(k, "")
            for needle, name in note_rules:
                if needle in notes or needle in k:
                    cr = name
                    break
        if not cr:
            cr = country_rules.get(country_map.get(k, "").lower(), "")
        if cr:
            result.at[i, "Craton"] = cr
            filled += 1

    if verbose:
        still = int((result["Craton"].isna() | (result["Craton"].astype(str).str.strip() == "")).sum())
        print(f"[craton] filled {filled} blank rows from locations tables; {still} still blank")
        print(result["Craton"].fillna("(blank)").value_counts().to_string())
    return result


def run_full_pipeline(csv_path_or_df, verbose=True, write_csv=None, write_full_csv=None, params=None):
    """
    Execute the complete calculation sequence that is currently implemented.

    Parameters
    ----------
    csv_path_or_df : str, Path, or pandas.DataFrame
        Either a path to the input CSV or an already-loaded DataFrame
        (useful for the website upload page).
    verbose : bool
        Print progress and screening statistics.
    write_csv : str or Path or None
        If given, write the final result table to this path.

    Returns
    -------
    result : pandas.DataFrame
        Full table with cations, P-T, all thermometers and screen flags.
    stats : dict
        Counts of how many samples passed each screening step.
    """
    if params is None:
        params = PARAMS
    ep = params.export

    # 0. Load
    if isinstance(csv_path_or_df, (str, Path)):
        df = load_and_clean(csv_path_or_df, verbose=verbose)
    else:
        df = csv_path_or_df.copy()

    # 1. Cations
    result = convert_oxides_to_cations(df, verbose=verbose)

    # 2. Iterative P–T + supporting thermometers
    if verbose:
        print("\nCalculating iterated TA98 temperatures and CA91 pressures …")
    T_TA98, P_CA91 = iterate_PT(result, verbose=verbose)
    result["T_TA98"] = T_TA98
    result["P_CA91"] = P_CA91
    result = calculate_all_temperatures(result, P_CA91, verbose=verbose)

    # 3. Equilibrium screening
    result, stats = apply_equilibrium_screening(result, verbose=verbose, params=params)

    # 4. Garnet Fe³⁺ iteration & empirical correction
    result = iterate_garnet_fe3(result, verbose=verbose, params=params)
    result = apply_empirical_fe3_correction(result, verbose=verbose, params=params)

    # 5. Stoichiometric Fe³⁺ (Droop 1987 / Quinn 2016)
    result = calculate_stoichiometric_fe3(result, verbose=verbose)

    # 6. Oxygen fugacity (Stagno 2013)
    result = calculate_fo2(result, verbose=verbose)

    # Optional write-out (slim public columns only)
    if write_csv is not None:
        out = Path(write_csv)
        export_map = [
            ("SampleID", "SampleID"),
            ("Location", "Location"),
            ("Craton", "Craton"),
            ("Lithology", "Lithology"),
            ("Reference", "Reference"),
            ("P_CA91", "P_CA91"),
            ("T_TA98", "T_TA98"),
            ("screen_pass", "Screen_pass"),
            ("XFe3_gt_iter", "XGtFe3+_iter"),
            ("XFe3_gt_corr", "XGtFe3+_corr"),
            ("delta_logfO2_FMQ", "logfO2(\u0394FMQ)"),
        ]
        export_cols = []
        rename = {}
        for src, dst in export_map:
            if src in result.columns:
                export_cols.append(src)
                if src != dst:
                    rename[src] = dst
        slim = result[export_cols].rename(columns=rename)

        # Rounding for publication-ready output
        if "P_CA91" in slim.columns:
            slim["P_CA91"] = slim["P_CA91"].round(ep.p_decimals)
        if "T_TA98" in slim.columns:
            slim["T_TA98"] = slim["T_TA98"].round(ep.t_decimals)
        if "XGtFe3+_iter" in slim.columns:
            slim["XGtFe3+_iter"] = slim["XGtFe3+_iter"].round(ep.xfe3_decimals)
        if "XGtFe3+_corr" in slim.columns:
            slim["XGtFe3+_corr"] = slim["XGtFe3+_corr"].round(ep.xfe3_decimals)
        if "logfO2(\u0394FMQ)" in slim.columns:
            slim["logfO2(\u0394FMQ)"] = slim["logfO2(\u0394FMQ)"].round(ep.fo2_decimals)

        # utf-8-sig (BOM) so Excel on Windows shows Unicode correctly (e.g. Δ)
        if write_full_csv is not None:
            Path(write_full_csv).parent.mkdir(parents=True, exist_ok=True)
            result.to_csv(write_full_csv, index=False, encoding="utf-8-sig")
            if verbose:
                print(f"Full working table written to {Path(write_full_csv).resolve()}")
        slim.to_csv(out, index=False, encoding="utf-8-sig")
        if verbose:
            print(f"\nResults written to {out.resolve()}")
            print(f"  Columns: {list(slim.columns)}")

    return result, stats
