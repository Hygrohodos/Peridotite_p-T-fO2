# Notes for the P–T–fO₂ supplementary information

Give this file to the chat that is writing the methods. The script `pt_calculator.py` is the P–T engine. It does not calculate fO₂, Fe³⁺ or the equilibrium screen. Those steps are specified below so they are not invented, and so they are not described as part of the thermometer equations.

Units: the script works in kbar internally; the manuscript reports pressure in GPa (1 GPa = 10 kbar). Quote GPa in the supplement text.

The script matches the website. `python pt_calculator.py --self-check` checks that. Anchor analysis: Jericho 22-4 (Kopylova et al. 1999), in `example_jericho_22-4.csv`.

- Taylor (1998) × Carswell (1991): **794 °C, 3.72 GPa** (37.2 kbar; unrounded 794.14 °C, 37.16 kbar).
- Nakamura (2009) at that pressure: **832 °C** (unrounded 831.58 °C). This is the published thermometer. An earlier workbook had `2·Xprp + Xalm` instead of `2·Xprp·Xalm` and gave 596 °C. Do not use 596 °C, and do not say Nakamura is 300 °C cold.

## What to call the preferred P–T

**Peridotite.** Temperature is Taylor (1998) two-pyroxene. Pressure is Carswell (1991) Al-in-orthopyroxene, the Nickel & Green (1985) formulation. They are solved together, starting at 40 kbar and 1000 °C, stopping when both have changed by less than 0.1, or after 25 cycles. Carswell pressure is limited to 0–80 kbar.

**Eclogite.** Temperature is the published Nakamura (2009) garnet–clinopyroxene thermometer. Pressure is the published Beyer et al. (2015) garnet–clinopyroxene barometer. Solved together from 50 kbar, tolerance 0.001, up to 1000 cycles. This is not the Beyer workbook. That sheet substituted raw XFe, XMg and XCa for Nakamura’s A, B and C terms, omitted −273, and put Ca in the clinopyroxene denominator. It also did not use the published Beyer activities. Neither error is in this script.

A result is kept for an eclogite only when the iteration converges, clinopyroxene Si is below 1.985 cations per 6 oxygens, pressure is between 5 and 80 kbar, and temperature is between 400 and 1600 °C. That P–T then replaces the peridotite solution for that sample.

**Olivine xenocrysts**, where Al-in-olivine is used. Temperature is Bussweiler et al. (2017) only:

T (°C) = (11245 + 46·P) / (13.68 − ln Al_ppm) − 273

with P in kbar and Al in ppm. Natural log. There is no Cr# term. De Hoog et al. (2010) is a different calibration and is not used. Pressure is not a pyroxene barometer. It is an Mg# projection, split at 91.5:

- Mg# ≥ 91.5: P = −16.78·Mg# + 1575.6 kbar
- Mg# < 91.5: P = −14.49·Mg# + 1356.9 kbar

Mg# = 100·Mg/(Mg+Fe) from olivine MgO and FeO.

## Nakamura (2009), the equation to print

T in °C, P in kbar. Garnet fractions are over Fe²⁺+Mn+Mg+Ca. All probe Fe is Fe²⁺.

- A = 0.5·Xgrs·(Xprp − Xalm − Xsps)
- B = 0.5·Xgrs·(Xprp − Xalm + Xsps)
- C = 0.5·(Xgrs + Xsps)·(Xprp − Xalm)
- KD = (Fe/Mg)garnet / (Fe/Mg)cpx
- XMg,cpx = Mg / (Al + Fe + Mg)
- XFe,cpx = Fe / (Al + Fe + Mg)

No calcium in the clinopyroxene denominator. Mn is included when it was analysed. The −273 is outside the fraction, not in the denominator.

T = {2784 + 14.52 P
  + (2601 + 1.44 P)·(2·Xgrs·Xprp − A)
  + (1183 + 6.98 P)·(Xgrs² − A)
  − 105·(2·Xgrs·Xalm + B)
  + (814.6 + 3.61 P)·(Xgrs² + B)
  − (254.6 + 8.42 P)·(2·Xprp·Xalm − Xalm² + C)
  − 83.6·(Xprp² − 2·Xprp·Xalm + C)
  + 1388·Xsps
  − 462·(XMg,cpx − XFe,cpx)}
  / {ln KD + 1.431
  + 0.695·(2·Xgrs·Xprp + Xgrs² − 2A)
  + 0.203·(Xgrs² − 2·Xgrs·Xalm)
  + 0.922·Xsps}
  − 273

The product `2·Xprp·Xalm` is the term that was previously wrong.

## Shared calculation rules

- Pyroxene cations on 6 oxygens, garnet on 12, olivine on 4. Oxide masses are the ones in the script. Cations are rounded to 4 decimal places, half to even, before the thermometers.
- Exchange thermometers use total Fe as Fe²⁺. Droop Fe³⁺ is used only inside the Beyer activity model, and is not fed back into Nakamura or Taylor.
- Kelvin offset is 273.0, except Brey & Köhler (1990) Al-in-opx pressure, which uses 273.15.
- Harley / Lee–Ganguly use R = 1.9872 cal. Nickel–Green / Carswell use R = 1.98726 cal. Brey & Köhler pressure uses 8.3143 J.

## Comparison calibrations in the script

These can be reported as comparisons evaluated at the preferred pressure. They are not the preferred solution.

Thermometers: Nimis & Grütter (2010) and Brey & Köhler (1990) Ca-in-opx; Brey & Köhler (1990) and Bertrand & Mercier (1985) two-pyroxene; Harley (1984), Lee & Ganguly (1988) and Nimis & Grütter (2009) opx–garnet; Krogh Ravna (2000), Krogh (1988), Ai (1994), Ellis & Green (1979), Powell (1985) and Råheim & Green (1974) garnet–cpx; Wu & Zhao (2007) and O’Neill & Wood (1979) olivine–garnet; Nimis & Taylor (2000) and the Ashchepkov correction of it; Boyd (1973), Wood & Banno (1973), Wells (1977) and Nehru & Wyllie (1974).

Barometers: Nickel & Green (1985), Taylor (1998) Al-in-opx, Brey & Köhler (1990), Harley & Green (1982), Harley (1984), Wood (1974).

Nimis & Grütter (2010) is not an independent thermometer. It is

T = −628.7 + 2.069·T_BK90 − 4.53×10⁻⁴·T_BK90²

with T_BK90 the Brey & Köhler Ca-in-opx temperature. Ai (1994) takes garnet Mg# on a 0–100 scale, not a fraction. Krogh Ravna (2000) takes pressure in GPa.

Do not recommend, and do not quietly “correct”:

- Sachtleben & Seck (1981), code `SS81`. The workbook expression on orthopyroxene Ca (cations per 6 oxygens) only spans about 812–850 °C. The published form was not established.
- Nickel & Brey (1984), code `NB84`. The workbook transcription is not a usable thermometer.
- Wood (1974), code `WO74`. Included as implemented. It is a comparison barometer and can return a meaningless pressure.

## What this script does not do

**Garnet Fe³⁺.** For peridotite fO₂, Fe³⁺/ΣFe in garnet is increased from 0.01 in steps of 0.001 until Taylor (1998) minus the corrected Wu & Zhao temperature is 29 °C or less. The correction is a flat −118 °C: T_WZ07,corr = T_WZ07,raw − 118. At T ≥ 1130 °C the empirical correction (manuscript Eq. 5) is added: ΔXFe³⁺ = 3.8319×10⁻⁴·T − 0.46673 (T in °C), clipped so Fe³⁺/ΣFe stays within 0–1. Eq. 5 is the regression at ΔT = 25 °C through the Slave samples inside the dashed envelope of Fig. 5b (n = 21); ΔT = 29 °C is the Fig. 6 optimum on those samples. Both are derived by `scripts/fit_calibration.py` in the main repository; `calculations/pt_fo2/params.py` is authoritative. Measured garnet Fe³⁺ is not required. This iteration is not in the script. `WZ07` here is the raw thermometer.

**Peridotite fO₂.** Stagno et al. (2013) garnet–olivine–orthopyroxene oxybarometer at the Taylor–Carswell P–T, using that iterated Fe³⁺. FMQ is the O’Neill (1987) expression used by Stagno, T in kelvin and P in kbar:

log fO₂ (FMQ) = 82.75 + 0.00484 T − 30681/T − 24.45 log₁₀ T + 94 P/T − 0.002 P

**Eclogite fO₂.** Stagno et al. (2015) andradite oxybarometer at the Nakamura–Beyer P–T. Garnet Fe³⁺ for that barometer is Droop charge balance on an 8-cation basis. It replaces the peridotite ΔFMQ when the eclogite P–T is accepted.

**Equilibrium screen**, applied to peridotite, not to the eclogite solution. A sample passes all of:

1. Olivine, orthopyroxene, clinopyroxene and garnet all present.
2. Taylor temperature above 900 °C.
3. |Taylor − Nimis & Grütter 2010| < 80 °C.
4. |Harley − (1203.7 − 0.9252·T_Taylor + 7.195802×10⁻⁴·T_Taylor²)| < 60 °C (median quadratic of HA84 on TA98 in 50 °C TA98 bins, fitted to samples that passed test 3). This replaced the earlier Krogh Ravna (2000) versus Harley test, which must not be used.
5. |Wu & Zhao − (88.99741613 + 0.9864930064·T_Harley)| < 50 °C.

Harley is Harley (1984). Nakamura is not a pass/fail gate. An older quadratic of Nakamura against Harley was fit to the cold workbook thermometer and must not be revived.

**Depth**, if a figure needs it: z (km) = P (kbar) × 10⁵ / (3300 × 9.80665), about 3.09 km per kbar.

## Wording to avoid

- Nakamura (2009) running ~300 °C below Taylor. That was the product/sum error.
- Al-in-olivine as “Bussweiler + De Hoog Cr#” or with an `1860·Cr#` term.
- The Beyer workbook cell as the barometer.
- Sachtleben & Seck or Nickel & Brey as working thermometers.
- A Nakamura-versus-Harley residual as an equilibrium test.
