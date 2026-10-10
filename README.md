# Peridotite P–T–fO2 by comparative geothermometry

Oxygen fugacity of peridotite xenoliths in kimberlite, estimated by matching
two-pyroxene (TA98) and olivine–garnet (WZ07) temperatures while iterating
garnet Fe3+/ΣFe, then applying the Stagno et al. (2013) oxybarometer.

Locked parameters (2026-10-10; `calculations/pt_fo2/params.py` is authoritative):

- starting garnet Fe3+/ΣFe = 0.01 (`xfe3_start`)
- WZ07 offset = flat −118 °C (`wz07_offset_c = 118`, slope = 0), from experiment F1300o
- thermometer agreement ΔT = 29 °C (`dt_match_c`), chosen in Fig. 6
- pyroxene–garnet screen = |HA84 − quadratic in TA98| < 60 °C
- high-T empirical correction (Eq. 5), applied at T ≥ 1130 °C:
  XFe3,corr = XFe3,iter + 3.8319×10⁻⁴·T − 0.46673 (T in °C, clipped to 0–1)

How Eq. 5 and ΔT are derived (`calculations/pt_fo2/calibration.py`, run by `scripts/fit_calibration.py`):

1. Eq. 5 is the least-squares line of (measured − estimated) Fe3+/ΣFe on TA98
   temperature at ΔT = 25 °C, through the Slave samples inside the dashed
   envelope of Fig. 5b: the 31 screened Slave samples less the blue (3) and
   red (3) error groups, T ≥ 1130 °C (n = 21). No level shift is applied.
2. Fig. 6 holds Eq. 5 fixed, sweeps ΔT from 0 to 75 °C and scores only those
   21 samples. ΔT = 29 °C has the highest ±0.02 rate (76 %; ties broken by the
   smallest mean bias).

Results (Table 1; `data/derived/Table1_calibration_stats.csv`). 1σ is the residual
standard deviation of estimated − measured Fe3+/ΣFe. Log units use 0.44 log units
ΔFMQ per 0.02 Fe3+/ΣFe.

| Dataset | n | ±0.02 | ±0.05 | 1σ | ≈ log units |
|---|---|---|---|---|---|
| Slave calibration | 31 | 65 % | 81 % | 0.045 | 0.99 |
| Slave, error groups excluded | 25 | 80 % | 100 % | 0.017 | 0.37 |
| Siberian + Kaapvaal validation | 27 | 63 % | 85 % | 0.035 | 0.77 |
| Calibration + validation | 58 | 64 % | 83 % | 0.040 | 0.88 |
| Combined, error groups excluded | 52 | 71 % | 92 % | 0.028 | 0.62 |

The Fig. 10 bars are computed from the last two rows.

## Run

```bash
pip install -r requirements.txt
python run_pipeline.py
python scripts/fit_calibration.py     # checks params.py; add --update-params to rewrite it
python scripts/make_figures.py
```

If `fit_calibration.py --update-params` changes `params.py`, run the pipeline and figures again.

Input: `data/master/xenolith_geochem.csv` (craton filled from `data/locations.csv` and `data/master/kimberlite_locations.csv`).
Outputs: `data/derived/p-T-fO2.csv` (publication table; `P_CA91` in kbar with an SI copy in `P_GPa`), `data/derived/screened_PT.csv` (full working table, pressures in kbar), `data/derived/emp_correction_fit.json`, `data/derived/dt_match_sweep_corr.csv`, `data/derived/Table1_calibration_stats.csv`, `data/derived/exp_TA98_WZ07_Fe3.csv`.
Figures: `Figures/Fig01`–`Fig10` (pressure axes in GPa). Fig. 4 is written by `scripts/plot_exp_ta98_wz07.py`, which `make_figures.py` calls.

`supplement/thermometer_choice/` is a standalone comparison calculator for the calibrations not used in the paper. It is not imported by the pipeline.

This bundle is the P–T–fO2 calculation only. Olivine FC/DFC, garnet diffusion, and geophysical inversion are not included.

## Repository layout

| Path | Contents |
|---|---|
| `run_pipeline.py` | Runs the full P–T–fO2 calculation |
| `calculations/pt_fo2/` | Calculation modules (stoichiometry, Fe3+ correction, thermobarometry, screening, fO2, buffers, locked parameters in `params.py`) |
| `scripts/` | Figure scripts and `fit_calibration.py` (derives Eq. 5 and ΔT) |
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

