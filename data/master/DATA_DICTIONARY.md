# Data Dictionary – Master Xenolith Table

File: `data/master/xenolith_geochem.csv`

This is the authoritative source table. All public results on Hygrohodos are derived from it.

## Required columns (minimum for the pipeline)

| Column | Description | Units / notes |
|--------|-------------|----------------|
| Sample ID | Unique sample identifier | text |
| Location | Kimberlite pipe / locality name | text (must match a name in `kimberlite_locations.csv` where possible) |
| Craton | Craton or tectonic setting | text |
| Lithology | Rock type | text |
| Reference | Literature source | text (include year) |
| cpx_SiO2 … cpx_Total | Clinopyroxene major elements | wt% |
| opx_SiO2 … opx_Total | Orthopyroxene major elements | wt% |
| garnet_SiO2 … garnet_Total | Garnet major elements | wt% |
| olivine_SiO2 … olivine_Total | Olivine major elements | wt% |

Optional but useful:
- `garnet_Fe3+`, `garnet_Fe2+` (if measured)
- Any additional oxides the original authors reported

## Notes for adding new data

- Keep the exact column names used in the existing master file.
- Prefer one row per sample.
- Always record the original paper / DOI in the Reference column.
- Coordinates belong in `kimberlite_locations.csv`, not in every sample row.
