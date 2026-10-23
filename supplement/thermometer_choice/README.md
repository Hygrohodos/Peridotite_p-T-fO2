# Thermometer comparison (not the manuscript workflow)

Standalone calculator for the calibrations considered and not used in Greene et al. (in prep.). It is not imported by `run_pipeline.py` and it does not calculate garnet Fe³⁺ or fO₂.

Use it to support the choice of TA98, WZ07, HA84 and CA91. Preferred peridotite pair is `-t TA98 -b CA91`. KR00, NA09, LG88, NG09, OW79 and the other codes are comparisons. SS81 and NB84 are not verified. Sudholz et al. (2022) is not in this file.

`NOTES_FOR_SUPPLEMENT.md` describes the screen (TA98–HA84 quadratic) and the Fe³⁺ iteration (ΔT = 29 °C, Eq. 5 fitted at ΔT = 25 °C) used in the manuscript. `calculations/pt_fo2/params.py` in the main repository is authoritative if they ever differ.

# P–T calculator

Standalone Python 3.9+ script (standard library only) for the thermometers and barometers used in the Hygrohodos P–T–fO₂ workflow. It reads oxide wt% and solves one thermometer with one barometer.

This is a calculation tool, not the equilibrium screen and not the oxybarometer. Read `NOTES_FOR_SUPPLEMENT.md` before describing the methods.

## Run

```text
python pt_calculator.py --list
python pt_calculator.py --self-check
python pt_calculator.py example_jericho_22-4.csv -t TA98 -b CA91
python pt_calculator.py chemistry.csv -t NA09 -b BE15 -o eclogite_pt.csv
python pt_calculator.py chemistry.csv -t BU17 -b MGP -o xenocryst_pt.csv
python pt_calculator.py chemistry.csv --interactive
python pt_calculator.py chemistry.csv -t TA98 -b CA91 --also-all -o all_at_that_pt.csv
```

`--self-check` must print Jericho 22-4 at about 794.2 °C and 37.17 kbar. That is the numerical anchor for the whole file.

## CSV

One row per sample. Oxide weights in wt%, blank if not analysed. Column names are `{phase}_{oxide}` with phases `cpx`, `opx`, `garnet` (or `gt`), `olivine` (or `ol`) and oxides `SiO2 TiO2 Al2O3 Cr2O3 FeO MnO MgO CaO Na2O NiO`. `SampleID`, `Location`, `Craton`, `Lithology` and `Reference` are copied through. Olivine aluminium for Bussweiler may be `olivine_Al2O3` or `olivine_Al_ppm`.

Cations are calculated in the script (6 oxygens for pyroxene, 12 for garnet, 4 for olivine, rounded to 4 decimal places). Do not pre-normalise. All iron is FeO, taken as Fe²⁺.

## Pairs to use

| Purpose | Command |
|---|---|
| Peridotite, the preferred solution | `-t TA98 -b CA91` |
| Eclogite, the preferred solution | `-t NA09 -b BE15` |
| Xenocryst Al-in-olivine at the Mg# pressure | `-t BU17 -b MGP` |

Any other pair is a comparison. `SS81` and `NB84` are not verified against the papers. `WO74` can return a pressure far outside the mantle range. Do not describe those three as results.

`--also-all` evaluates every other calibration at the solved P and T. The solved `T_C` is the last thermometer call, made at the previous pressure. A thermometer re-evaluated at the final pressure can differ by less than the convergence tolerance (0.1 °C / 0.1 kbar, or 0.001 when the barometer is `BE15`). Quote the solved pair, not the re-evaluation.
