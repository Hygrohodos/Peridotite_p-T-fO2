# Peridotite P–T–fO2 by comparative geothermometry

Oxygen fugacity of peridotite xenoliths in kimberlite, estimated by matching
two-pyroxene (TA98) and olivine–garnet (WZ07) temperatures while iterating
garnet Fe3+/ΣFe, then applying the Stagno et al. (2013) oxybarometer.

> Revision v7_cr (15 Nov 2026) is the current method and the default of `calculations/pt_fo2/params.py`. The
> sections after "Revision v7_cr" (locked parameters, level, results, fO2 uncertainty) describe revision v3 and are
> kept for history. Figure numbers follow manuscript v7 (Figs 1–10).

## Revision v7_cr: TA98 + garnet-Cr2O3 calibration correction (Eq. 4)

Full method, coefficients and validity range: `calculations/pt_fo2/cr_correction.py` (docstring).
Coefficients are reproduced from the pipeline by `python3 scripts/fit_cr_correction.py --check`.

Garnet Fe3+/ΣFe is raised from 0.01 in 0.001 steps until

    T_TA98 − T_WZ07(Fe2+)  ≤  a + b·t + c·(Cr2O3 − 5.364)·t,      t = T_TA98 − 1207.02 °C

i.e. T_corr = T_WZ07 − 102 − 0.285·t + 0.067·(Cr2O3 − 5.36)·t, with garnet Cr2O3 in wt%.

| Coefficient | Value | 1σ |
|---|---|---|
| a (adopted level at the pivot, ΔTmatch; Fig. 6 ±0.02 optimum) | −102 °C | — |
| a of the joint WLS fit (not used as ΔTmatch) | −92.77 °C | 2.92 |
| b (slope at Cr2O3 = 5.364 wt%) | −0.2850 °C/°C | 0.0472 |
| c (Cr2O3 × T term) | +0.0673 °C/°C per wt% | 0.0199 (p = 0.0007) |

- **Joint fit.** 37 natural xenoliths with measured garnet Fe3+ (15 Slave, 13 Kaapvaal, 9 Siberia; screened
  calibration + validation sets, T_TA98 ≥ 1120 °C, without the four error-group samples and Uv-50-04; y = TA98 −
  WZ07 at the measured Fe3+) plus the 8 Matjuschkin et al. (2014) runs (paper σ_i, own level, so they constrain only
  b and c). Weighted least squares, natural σ = 17.43 °C. The natural samples alone give the same coefficients
  within error (c = 0.065 ± 0.020). The adopted level a is the Fig. 6 ±0.02 optimum (−102 °C), not the WLS intercept (−92.77 °C); b and c are held fixed.
- **Slope check.** At the mean Matjuschkin garnet Cr2O3 (1.33 wt%) the correction predicts a TA98 − WZ07 slope of
  −0.557 °C/°C; the weighted fit of the 13 experiments (Fig. 5) gives −0.555 ± 0.216.
- **Validity range and clamping.** Cr2O3 1.3–9.8 wt% and T_TA98 1124–1394 °C. Outside, the input is held at the
  limit (no extrapolation); samples without garnet Cr2O3 use 5.364 wt%. Flagged samples are listed in
  `cr_correction_flags.csv` (with `screen_pass`). Of the 424 screened xenoliths, 231 are below 1124 °C, 6 above
  1394 °C, 2 below 1.3 wt% Cr2O3 and 2 above 9.8 wt%.
- **Where it is applied.** After the analytical-quality, four-phase, T > 900 °C and equilibrium screens (which use
  uncorrected temperatures), in the garnet Fe3+ iteration (`fe3_correction.iterate_garnet_fe3`). P and TA98 are
  unchanged.
- **Build.** `PAPER_FIG_OUT=/some/dir python3 scripts/paper_figures/make_paper_figures.py` writes Figs 1–10 (Fig. 5 =
  natural samples split at 4.9 wt% Cr2O3 plus the experiments, `fig5_experiments_cr.py`), Table 1, the text numbers,
  `cr_correction_flags.csv` and `cr_correction.json`. `--v7-eq4` rebuilds manuscript v7 (experimental slope 0.555,
  Slave-median level −103.8 °C).

Results (Table 1, v7_cr):

| Dataset | n | ±0.02 | ±0.05 | 1σ Fe3+/ΣFe | 1σ log units |
|---|---|---|---|---|---|
| Slave calibration | 23 | 61 % | 83 % | 0.046 | 0.75 |
| Slave, error groups excluded | 19 | 74 % | 100 % | 0.020 | 0.43 |
| Siberian + Kaapvaal validation | 26 | 69 % | 85 % | 0.035 | 1.25 |
| Calibration + validation | 49 | 65 % | 84 % | 0.041 | 1.04 |
| Combined, error groups excluded | 45 | 71 % | 91 % | 0.030 | 0.99 |

The Siberian + Kaapvaal samples above 1120 °C are part of the correction fit, so the validation row is not an
independent test in v7_cr. Of the 424 screened xenoliths, 136 (32 %) stay at the 0.01 starting value.

## Revision v3 (history)

Locked parameters (revision v3, 2026-11-11; `calculations/pt_fo2/params.py` is authoritative):

- starting garnet Fe3+/ΣFe = 0.01 (`xfe3_start`), step 0.001
- no WZ07 correction (`wz07_offset_c = 0`, `wz07_offset_slope = 0`); Eq. 5 off (`apply_empirical_correction = False`)
- stopping rule: iterate until TA98 − WZ07(Fe2+) ≤ −96 °C (`dt_match_c`)
- pyroxene–garnet screen = |HA84 − quadratic in TA98| < 60 °C

How the −96 °C level is set (`calculations/pt_fo2/level_anchor.py`, `scripts/anchor_level.py`):
it is the median of TA98 − WZ07 at the measured garnet Fe3+/ΣFe over all 30 screened Slave
calibration samples (bootstrap 95 % CI −108 to −81 °C). It uses no validation sample and no
outcome-based exclusion, and was fixed before any scoring. For reference, the 52 calibration +
validation xenoliths outside the error groups give −96 ± 32 °C (revision v3; no v7 figure). The 13 experiments
show a different level and a slope (−0.55 ± 0.22 per °C, Fig. 5); the natural samples do not
(−0.06 ± 0.05 over the 52), so the experimental fit is shown but not applied.
`python3 scripts/anchor_level.py --check` fails if `params.py` and the anchor disagree.

History: the original paper used a flat −118 °C offset (experiment F1300o), ΔT = 29 °C and the
high-T Eq. 5 (kept in `params.py` for reproducing it, and in `data/derived/legacy_original/`).
Revision v2 used the weighted experimental slope; it was withdrawn because it biased the
residuals with temperature (+0.023 Fe3+/ΣFe per 100 °C) and put 35 % of samples at the floor.

Results (Table 1; `data/derived/Table1_calibration_stats.csv`). 1σ is the residual
standard deviation of estimated − measured Fe3+/ΣFe. The log-unit 1σ is the
standard deviation of ΔFMQ(estimated) − ΔFMQ(measured), with ΔFMQ calculated
(Stagno et al. 2013) from the estimated and from the measured garnet Fe3+/ΣFe at
the same P, T and mineral chemistry (`calibration.fo2_residuals`).

| Dataset | n | ±0.02 | ±0.05 | 1σ Fe3+/ΣFe | 1σ log units |
|---|---|---|---|---|---|
| Slave calibration | 30 | 50 % | 80 % | 0.044 | 0.63 |
| Slave, error groups excluded | 25 | 60 % | 96 % | 0.024 | 0.42 |
| Siberian + Kaapvaal validation | 27 | 56 % | 93 % | 0.032 | 0.98 |
| Calibration + validation | 57 | 53 % | 86 % | 0.039 | 0.81 |
| Combined, error groups excluded | 52 | 58 % | 94 % | 0.029 | 0.76 |
| Charge balance (Droop 1987; Quinn et al. 2016), screened Slave | 30 | 23 % | 60 % | 0.085 | – |
| Charge balance (Droop 1987; Quinn et al. 2016), all Slave before screening | 57 | 25 % | 53 % | 0.087 | – |

The residuals show no trend with temperature (−0.004 ± 0.005 per 100 °C, p = 0.46, n = 52).
Of the 450 fully screened xenoliths, 89 (19.8 %) stay at the 0.01 starting value (the "floor");
86 of them are below 1100 °C, where only 8 of the 57 calibration/validation samples lie, and
none of the 57 is at the floor. `scripts/v3_text_numbers.py` and `scripts/v3_level_sensitivity.py`
write these and the other numbers quoted in the text; `scripts/table1_three_way.py` compares
original / v2 / v3 using the reference outputs in `data/reference/`.

Charge-balance stoichiometry (`calculations/pt_fo2/stoichiometry.py`) implements
Droop (1987) and Quinn et al. (2016) separately (`XFe3_droop`, `XFe3_quinn` in
`screened_PT.csv`). For anhydrous garnet the two procedures are algebraically
identical (Fe3+ = 24(1 − 8/S) per 8-cation formula), so the columns agree to
rounding and Table 1 reports one charge-balance row per sample set. Mean bias
(estimated − measured) is +0.017 for the 30 screened Slave samples, against
−0.008 for comparative geothermometry. `tests/test_stoichiometry.py` reproduces
the Droop (1987) Table 3 garnet (Fe3+ 0.779, Fe2+ 0.013 per 12 O).

fO2 uncertainty (`python scripts/fo2_uncertainty.py`; `data/derived/fO2_uncertainty.json`,
`data/derived/fO2_residuals.csv`, `Figures/diagnostics/fO2_residual_vs_Fe3.png`):
the 57 calibration + validation samples (error groups retained) are split at their
median estimated ΔFMQ (−2.45). The 1σ is 0.83 log units for the reduced half
(n = 28) and 0.71 for the oxidised half (n = 29), 0.81 overall. These are the two
Fig. 10 bars. They do not apply to floor samples, whose estimates are censored at 0.01.

Screening (1100 samples): four phases 747 → analytical quality (Nimis and Grütter 2010: oxide totals 98.5–101.5 wt%, pyroxene cation sums ≥ 3.990 on 6 O, garnet Cr2O3 > 1 wt% or Mg# > 0.6) 695 → TA98 > 900 °C 565 → opx–cpx 545 → pyx–gt 467 → gt–ol 424 fully screened. Calibration (Slave, measured Fe3+) 57 → 23; validation (Siberian + Kaapvaal) 52 → 26; neither set was affected by the de-duplication.

## Run

```bash
pip install -r requirements.txt
python3 run_pipeline.py                 # data/derived/p-T-fO2.csv, screened_PT.csv
python3 scripts/anchor_level.py --check # level check (writes data/derived/level_anchor.json)
python3 scripts/table1_stats.py         # Table 1
python3 scripts/fo2_uncertainty.py      # log-unit column of Table 1, Fig. 10 bars
PAPER_FIG_OUT=/some/dir bash scripts/paper_figures/make_all.sh   # Figs 1-10 + text numbers
python3 -m pytest tests
```

`scripts/make_figures.py` now writes Table 1 only (its legacy figure functions draw the original
Eq. 5 calibration and run only if `apply_empirical_correction` is switched back on).
`scripts/fit_calibration.py` belongs to the original calibration and is not used.

Input: `data/master/xenolith_geochem.csv` (1100 rows after the 2026-10-10 de-duplication, see `data/master/DATA_DICTIONARY.md` and `data/master/qc/`; craton filled from `data/locations.csv` and `data/master/kimberlite_locations.csv`).
Outputs: `data/derived/p-T-fO2.csv` (publication table; `P_CA91` in kbar with an SI copy in `P_GPa`), `data/derived/screened_PT.csv` (full working table, pressures in kbar), `data/derived/level_anchor.json`, `data/derived/Table1_calibration_stats.csv`, `data/derived/exp_TA98_WZ07_Fe3.csv`. `data/derived/legacy_original/` holds the stale Eq. 5 outputs for reference.
Figures: `make_all.sh` writes Fig1–Fig10 (`FigN_final.png`, numbered as in manuscript v7), `Fig4_per_sample.csv`, `Fig5_weights.csv`, `Fig5_equation.json` (fit coefficients at 1350 °C and at the 1207.0 °C pivot), `Online_Resource_1/` (the run/weight table moved out of Fig. 5), `Fig6_level_sweep.csv`, `level_anchor.json`, `figure_stats.json`, `text_numbers.json`, `level_sensitivity.csv` and `Table1_three_way.md` to `$PAPER_FIG_OUT` (default `scripts/paper_figures/out`, git-ignored); it writes nothing into `data/`.

Figure style: every figure takes its style from `scripts/paper_figures/unified.py`, the single house-style module (approved 14 Nov 2026). It sets CMP widths of 84 or 174 mm, Arial-metric 8 pt lettering, bold lower-case panel letters outside the axes, frameless keys, a colour-blind-checked palette and 600 dpi output. Online Resource 2 (`figS_slope_extension.py`) and the legacy/diagnostic plots use it too. A stand-alone copy (`house_style.py`, `house_style.mplstyle`, README of the specifications) is kept in OneDrive `Comparative Geothermometry/figure_style/` for other papers.

`supplement/thermometer_choice/` is a standalone comparison calculator for the calibrations not used in the paper. It is not imported by the pipeline.

This bundle is the P–T–fO2 calculation only. Olivine FC/DFC, garnet diffusion, and geophysical inversion are not included.

## Weighted TA98 − WZ07 regression on experiments (Fig. 5; its slope was Eq. 4 in manuscript v7; v7_cr uses it as the slope check)

`scripts/weighted_thermometer_regression.py` is a stand-alone module and CLI
(Python ≥ 3.9, numpy, pandas; no other repository imports; all temperatures in °C).
It weights each experiment by its apparent disequilibrium and fits
TA98 − WZ07(Fe2+ only) = b0 + b1 (T_TA98 − 1350) by weighted least squares:

- excess_i = quadrature sum of max(0, |ΔT_j| − t_j) over the three screens
  (opx–cpx 80, pyx–gt 60, gt–ol 50 °C; `SCREEN_THRESHOLDS_C`). A run fails
  screen j when |ΔT_j| > t_j, i.e. exactly when that term is non-zero.
- σ_i² = σ_analytical,i² + excess_i², with σ_analytical defaulting to
  `SIGMA0_C` = 50 °C
- w_i = 1/σ_i²; n_eff = (Σw)²/Σw²
- parameter covariance = (XᵀWX)⁻¹ × max(1, χ²_red)

**Methods notes.**
- σ_i is treated as an absolute 1σ uncertainty on y. But σ0 = 50 °C is a
  placeholder (the gt–ol acceptance width, pending the analytical-quality
  assessment) and the excess term is heuristic, so the absolute scale is an
  assumption, not a measurement.
- The covariance is scaled by max(1, χ²_red): inflated if the scatter
  exceeds the σ_i, never shrunk. For the 13 current runs χ²_red = 0.84
  (p = 0.60), so the SEs are unscaled; always scaling would give ±24.5 and
  ±0.198 instead of ±26.8 and ±0.216.
- The weights come from the same thermometers as y: ΔT_opx–cpx contains
  TA98 and correlates with y (r = −0.48 for the 13 runs); ΔT_gt–ol contains
  WZ07 (r = 0.05). The weights are therefore not independent of the
  response. The case-resampling bootstrap carries each run's weight with it.

```bash
python3 scripts/weighted_thermometer_regression.py               # current 13 runs
python3 scripts/weighted_thermometer_regression.py --robustness  # + OLS, leave-one-out, bootstrap
python3 -m pytest tests/test_weighted_regression.py               # or: python3 tests/test_weighted_regression.py
python3 scripts/make_fig05_weightcircles.py --pivot 1207.0        # Fig. 5 + Online Resource 1 weight table (uses the module; without --pivot the pivot comes from a pipeline run)
python3 scripts/make_fig05_weightcircles.py --input my_runs.csv --outdir out/
```

Outputs (default `data/derived/weighted_regression/`): `weighted_regression_weights.csv`,
`weighted_regression_fit.json`, `weighted_regression_fit.txt` (and `weighted_regression_loo.csv`
with `--robustness`). Options: `--input`, `--outdir`, `--t-ref`, `--sigma0`,
`--thr-opx-cpx`, `--thr-pyx-gt`, `--thr-gt-ol`, `--col NAME=CSVCOLUMN`.
Invalid input stops with a message naming the rows and exit status 2. The
checks cover blank, NaN or non-numeric T or ΔT; bad `include` or
`sigma_analytical_C` values; fewer than 3 included runs; and bad or clashing
`--col`. Suspicious input gives a warning but is still fitted: a median
temperature above 1600 (Kelvin?), T outside 0–2000 °C, or near-duplicate runs.

The Fig. 5 script runs from a clean checkout (needs numpy, pandas,
matplotlib and this repository). Its two cross-checks are optional and are
skipped with a notice when they cannot apply (`--no-crosscheck` skips both):
the raw-data pipeline (`calculations/pt_fo2/exp_calibration`, only when the
input runs equal the pipeline's 13), and the item-1 results file (in-repo copy
`data/reference/item1_experiment_weights.csv`, which also supplies the sensitivity-only
quality-factor column). New
references are styled from a fallback cycle of the paper palette. New runs
get a default label offset; hand-tune `LABEL_OFF` if labels collide.

### Adding experiments

The input CSV has one row per run with columns `SampleID`, `T_TA98`, `T_WZ07_Fe2`
(WZ07 with the measured garnet Fe3+ removed), `dT_opx_cpx`, `dT_pyx_gt`, `dT_gt_ol`
(all °C, computed exactly as for the xenoliths in Section 3.2), and optionally
`Reference`, `sigma_analytical_C` (blank = 50 °C default; otherwise finite and > 0)
and `include` (blank or 1/true/yes = used; 0/false/no = excluded; anything else is an error).

1. Preferred: add the new run's mineral analyses to the experimental input read by
   `scripts/plot_exp_ta98_wz07.py` and rerun it; it rewrites
   `data/derived/exp_TA98_WZ07_Fe3.csv` with the thermometers and screens computed
   the same way as for every other run.
2. Or copy that CSV, append rows with the columns above, and pass `--input your.csv`.
3. Rerun the CLI and `make_fig05_weightcircles.py` (with the same `--input`). The test pins the published
   13-run result (−14 ± 27, −0.555 ± 0.216, n_eff 8.4); when the data set is
   deliberately changed, update those numbers in the test together with the paper.

## Repository layout

| Path | Contents |
|---|---|
| `run_pipeline.py` | Runs the full P–T–fO2 calculation |
| `calculations/pt_fo2/` | Calculation modules (stoichiometry, Fe3+ correction, thermobarometry, screening, fO2, buffers, locked parameters in `params.py`) |
| `scripts/` | Level anchor, Table 1, text-number and figure scripts (`paper_figures/make_all.sh`); legacy `fit_calibration.py` |
| `data/master/` | Input xenolith geochemistry and kimberlite locations (see `DATA_DICTIONARY.md`) |
| `data/derived/` | Calculated P–T–fO2 results |
| `data/spatial/web/` | GeoJSON layers for the locality map |
| `Figures/` | `original_submission/`: original-submission figures under their old numbers (not v7); `diagnostics/`: check plots. The v7 paper figures are built into `$PAPER_FIG_OUT` |
| `supplement/thermometer_choice/` | Standalone thermometer comparison calculator |
| `tests/` | Unit tests (charge-balance garnet Fe3+ against Droop 1987 Table 3; weighted regression reproduces the Fig. 5 fit; level anchor = params) |

To get everything at once, use **Code → Download ZIP** at the top of this page, or clone:

```bash
git clone https://github.com/Hygrohodos/Peridotite_p-T-fO2.git
```

## License

MIT, see [`LICENSE`](LICENSE).

