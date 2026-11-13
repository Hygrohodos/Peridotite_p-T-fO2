#!/usr/bin/env python3
"""Figures 1-11, Table 1 and the numbers quoted in the text of Greene et al. (revision).

Usage:  python3 scripts/paper_figures/make_paper_figures.py [--extend-slope]
        (output folder: $PAPER_FIG_OUT, default scripts/paper_figures/out; nothing in data/ is changed)
        --extend-slope: Supplementary Data option. The slope of c(T) is extended linearly below the coldest
        experiment (1162 deg C) instead of being held flat there, and the level is re-derived by the same
        Slave-30 rule (MATCH_LEVEL_EXTENDED_C). Use a separate $PAPER_FIG_OUT; the paper uses the default.

Data flow
    experimental runs (scripts/weighted_thermometer_regression.py)
      -> weighted fit c(T) = TA98 - WZ07 (Fig. 4), held flat outside the experimental TA98 range
    data/master/xenolith_geochem.csv
      -> calculation pipeline stopped at MATCH_LEVEL_C, with the slope of c(T) about PIVOT_T_C
                                                            -> screened_PT.csv, p-T-fO2.csv
      -> fO2 residuals and 1 sigma; Table 1                 -> fO2_residuals.csv, fO2_uncertainty.json,
                                                               Table1_calibration_stats.csv
      -> Figs 5-11, and the numbers quoted in the text, all from these files.
Figs 1-4 (localities, Fe3+ effect on WZ07, equilibrium screens, experiments) do not depend on the level.
"""

# ---------------------------------------------------------------------------------------------------
# Match level of Eq. 4. Garnet Fe3+/SumFe is raised from 0.01 until
#     T_TA98 - T_WZ07(Fe2+) - [c(T_TA98) - c(PIVOT_T_C)] <= MATCH_LEVEL_C,
# where c(T) is the weighted experimental fit of Fig. 4 (its slope applied about the pivot, held flat outside
# the experimental TA98 range). Over the 30 Slave calibration samples, PIVOT_T_C is the median T_TA98 and
# MATCH_LEVEL_C the median of the slope-removed T_TA98 - T_WZ07(measured Fe3+) (both checked in
# write_statistics). This is the bold line in Fig. 5; Figs 6-11 and every number in the text use it.
MATCH_LEVEL_C = -99.1
PIVOT_T_C = 1208.5
# The same rule with the slope extended below the experiments (--extend-slope; not used in the paper).
MATCH_LEVEL_EXTENDED_C = -101.3
# ---------------------------------------------------------------------------------------------------

import argparse, json, shutil, subprocess, sys
from pathlib import Path
import pandas as pd

from load import REPO, PARAMS, calib
from calculations import run_full_pipeline
from calculations.pt_fo2.level_anchor import anchor_level
from common import set_match_level, to_pipeline

HERE = Path(__file__).resolve().parent
SCRIPTS = Path(REPO) / "scripts"
sys.path.insert(0, str(SCRIPTS))


def set_wz07_correction(extend_below=False):
    """Write the weighted experimental fit c(T) = b0 + b1 (T - t_ref) into PARAMS.fe3 as the WZ07 correction.

    The pipeline subtracts offset(T) = wz07_offset_c + wz07_offset_slope * T from WZ07, so offset = -c(T);
    T is held within the TA98 range of the experiments, or only below its maximum if extend_below. Only the
    slope matters: the level is set about the pivot (set_match_level), which cancels the intercept.
    """
    import weighted_thermometer_regression as wtr
    runs = wtr.load_runs()
    _, fit = wtr.fit(runs)
    fp = PARAMS.fe3
    fp.wz07_offset_slope = -fit.b1
    fp.wz07_offset_c = -(fit.b0 - fit.b1 * fit.t_ref)
    fp.wz07_offset_tmin_c = None if extend_below else float(runs.T_TA98.min())
    fp.wz07_offset_tmax_c = float(runs.T_TA98.max())
    return fit


def calculate_at_match_level(level_c, pivot_c, out):
    """Run the pipeline at match level level_c about pivot pivot_c; write and return the full table."""
    assert not PARAMS.fe3.apply_empirical_correction, "the method applies no Eq. 5 (calculations/pt_fo2/params.py)"
    set_match_level(level_c, pivot_c)
    run_full_pipeline(Path(REPO) / "data" / "master" / "xenolith_geochem.csv", verbose=False,
                      write_csv=out / "p-T-fO2.csv", write_full_csv=out / "screened_PT.csv")
    table = pd.read_csv(out / "screened_PT.csv", encoding="utf-8-sig", low_memory=False)
    return table.rename(columns={table.columns[0]: "SampleID"})


def write_statistics(table, level_c, pivot_c, out):
    """Write the fO2 uncertainty, Table 1, the Slave-30 anchor and the level sensitivity for the table."""
    from fo2_uncertainty import write_uncertainty
    from level_sensitivity import level_sensitivity
    anchor = anchor_level(table)
    assert abs(anchor["median_T_TA98_C"] - pivot_c) < 0.5, f"PIVOT_T_C is not the Slave-30 median TA98 ({anchor['median_T_TA98_C']:.1f})"
    median_at_pivot = anchor["median_r_C"] - to_pipeline(0.0)
    assert abs(median_at_pivot - level_c) < 0.1, f"MATCH_LEVEL_C is not the Slave-30 median at the pivot ({median_at_pivot:.2f})"
    (out / "level_anchor.json").write_text(json.dumps(anchor, indent=1, default=float))
    write_uncertainty(table, out)
    calib.calibration_stats(table).to_csv(out / "Table1_calibration_stats.csv", index=False)
    level_sensitivity(table, level_c, to_pipeline).to_csv(out / "level_sensitivity.csv", index=False)


def draw_figures(table):
    """Draw Figs 1-3 and 5-11 into the output folder; Figs 5-11 from the table at the match level."""
    from unified import save
    import fig1, fig2, fig3
    from fig5_slave_level import draw_fig5
    from fig6_fe3_match_level import draw_fig6
    from fig7_level_sensitivity import draw_fig7
    from fig8_fe3_match_level import draw_fig8
    from fig9_validation import draw_fig9
    from fig10_fo2_measured import draw_fig10
    from fig11_fo2_estimated import draw_fig11
    for n, draw in ((1, fig1.fig1), (2, fig2.fig2), (3, fig3.fig3)):
        fig, px = draw()
        save(fig, f"Fig{n}_final.png", px)
    for n, draw in ((5, draw_fig5), (6, draw_fig6), (7, draw_fig7), (8, draw_fig8), (9, draw_fig9), (10, draw_fig10), (11, draw_fig11)):
        fig, px = draw(table)
        save(fig, f"Fig{n}_final.png", px)


def draw_fig4(out):
    """Fig. 4 (the 13 experiments and their weighted fit), drawn by scripts/make_fig04_weightcircles.py."""
    tmp = out / "_fig4"
    run("make_fig04_weightcircles.py", "--outdir", tmp)
    (tmp / "Fig04_TA98_WZ07_measured_Fe3.png").replace(out / "Fig4_final.png")
    (tmp / "Fig04_weights.csv").replace(out / "Fig4_weights.csv")
    shutil.rmtree(tmp)


def write_quoted_numbers(out):
    """Numbers quoted in the text and the three-way Table 1, from the output folder."""
    run("text_numbers.py", "--derived", out, "--out", out / "text_numbers.json")
    run("table1_three_way.py", "--derived", out, "--out", out / "Table1_three_way.md")


def run(script, *args):
    """Run one of the repository scripts in scripts/."""
    subprocess.run([sys.executable, str(SCRIPTS / script), *map(str, args)], check=True, stdout=subprocess.DEVNULL)


if __name__ == "__main__":
    from unified import OUTD
    cli = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    cli.add_argument("--extend-slope", action="store_true", help="extend the slope below 1162 deg C (Supplementary Data)")
    extend = cli.parse_args().extend_slope
    level = MATCH_LEVEL_EXTENDED_C if extend else MATCH_LEVEL_C
    out = Path(OUTD)
    (out / "figure_stats.json").unlink(missing_ok=True)
    fit = set_wz07_correction(extend)
    table = calculate_at_match_level(level, PIVOT_T_C, out)
    write_statistics(table, level, PIVOT_T_C, out)
    draw_figures(table)
    draw_fig4(out)
    write_quoted_numbers(out)
    shape = "extended" if extend else "held flat"
    print(f"WZ07 correction c(T) = {fit.b0:.1f} {fit.b1:+.3f} (T - {fit.t_ref:.0f}), {shape} below the experiments; "
          f"match level {level:g} deg C at {PIVOT_T_C:g} deg C; outputs in {out}")
