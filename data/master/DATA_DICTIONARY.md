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
- One row per analysis. Do not add a second copy of a published analysis (e.g. a compilation's re-listing of the original paper): keep the row with the real citation, valid Location/Craton and measured garnet Fe3+.
- Text fields carry no leading/trailing spaces; the header has no non-breaking spaces.
- Always record the original paper / DOI in the Reference column.
- Coordinates belong in `kimberlite_locations.csv`, not in every sample row.

## Encoding

The file is Windows-1252 (cp1252), LF line endings, no quoted fields. The only non-ASCII character is `é` (0xE9) in `Gahcho_Kué`. The pipeline reads it with `encoding="cp1252"`; do not re-save it as UTF-8 without also changing the loader and the location tables.

## Quality control log (`qc/`)

2026-10-10 de-duplication (1171 → 1100 rows):

- `qc/duplicates_removed_2026-10-10.csv`: the 71 rows removed (including 38-Mx1, a probable transcription duplicate of 38-1 whose garnet CaO and olivine FeO differ by 0.11 wt%, removed on author approval). A row was treated as a duplicate when every oxide reported in both rows agreed within 0.05 wt% (≥ 26 common oxides), whether or not the IDs matched. The kept row was chosen by, in order: measured garnet Fe3+, a real citation over a `RefNNN` placeholder, a non-blank craton, a valid location, more oxides reported, more decimal places. No measured Fe3+ value was lost (all were already on the kept row). Kept rows were not otherwise edited.
- `qc/duplicates_flagged_2026-10-10.csv`: same-ID (or near-identical) pairs that were **kept** because the analyses differ by more than 0.05 wt% (different grains, re-analyses or different studies). These need an author decision.
- `qc/whitespace_fixes_2026-10-10.csv`: leading spaces stripped from 69 Sample IDs; non-breaking spaces in `vr19674 pg2` / `vr19674 pg3` replaced by ordinary spaces.
