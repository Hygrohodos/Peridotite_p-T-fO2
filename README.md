# Peridotite P–T–fO2 by comparative geothermometry

Oxygen fugacity of peridotite xenoliths in kimberlite, estimated by matching
two-pyroxene (TA98) and olivine–garnet (WZ07) temperatures while iterating
garnet Fe3+/ΣFe, then applying the Stagno et al. (2013) oxybarometer.

Locked parameters (2026-10-10):

- starting garnet Fe3+/ΣFe = 0.01 (`calculations/pt_fo2/params.py`, `xfe3_start`)
- WZ07 offset = flat −118 °C (`wz07_offset_c = 118`, slope = 0)
- thermometer agreement ΔT = 26 °C
- pyroxene–garnet screen = |HA84 − quadratic in TA98| < 60 °C
- high-T empirical correction applied at T ≥ 1130 °C

Calibration (Slave, measured Fe3+, screened): 58 % within ±0.02, 81 % within ±0.05 (n = 31).
Validation (Siberian + Kaapvaal, measured Fe3+, screened): 63 % / 93 % (n = 27).
Residual 1σ calibration + validation, error groups included: 0.039 (≈ 0.86 log ΔFMQ).
Residual 1σ with the 3 red and 3 blue error groups excluded: 0.027 (≈ 0.59 log ΔFMQ).

## Run

```bash
pip install -r requirements.txt
python run_pipeline.py
python scripts/make_figures.py
```

Input: `data/master/xenolith_geochem.csv` (craton filled from `data/locations.csv` and `data/master/kimberlite_locations.csv`).
Outputs written by the pipeline: `data/derived/p-T-fO2.csv`, `data/derived/screened_PT.csv`.
Figures: `Figures/Fig01`–`Fig10`.

`supplement/thermometer_choice/` is a standalone comparison calculator for the calibrations not used in the paper. It is not imported by the pipeline.

This bundle is the P–T–fO2 calculation only. Olivine FC/DFC, garnet diffusion, and geophysical inversion are not included.

## Repository layout

| Path | Contents |
|---|---|
| `run_pipeline.py` | Runs the full P–T–fO2 calculation |
| `calculations/pt_fo2/` | Calculation modules (stoichiometry, Fe3+ correction, thermobarometry, screening, fO2, buffers, locked parameters in `params.py`) |
| `scripts/` | Figure scripts |
| `data/master/` | Input xenolith geochemistry and kimberlite locations (see `DATA_DICTIONARY.md`) |
| `data/derived/` | Calculated P–T–fO2 results |
| `data/spatial/web/` | GeoJSON layers for the locality map |
| `Figures/` | Paper figures (PNG and PDF) |
| `supplement/thermometer_choice/` | Standalone thermometer comparison calculator |

To get everything at once, use **Code → Download ZIP** at the top of this page, or clone:

```bash
git clone https://github.com/Hygrohodos/Peridotite_p-T-fO2.git
```

## License

MIT, see [`LICENSE`](LICENSE).

