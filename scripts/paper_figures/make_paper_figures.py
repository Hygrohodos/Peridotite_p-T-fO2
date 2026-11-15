#!/usr/bin/env python3
"""Figures 1-10, Table 1 and the numbers quoted in the text of Greene et al. (revision v7_cr).

Usage:  python3 scripts/paper_figures/make_paper_figures.py            (revision v7_cr, the paper)
        python3 scripts/paper_figures/make_paper_figures.py --v7-eq4   (revision v7: experimental slope, Slave median)
        (output folder: $PAPER_FIG_OUT, default scripts/paper_figures/out; nothing in data/ is changed)
        --extend-slope (with --v7-eq4 only): the v7 slope extended linearly below 1162 deg C.

Data flow (v7_cr)
    data/master/xenolith_geochem.csv
      -> calculation pipeline: cations, P-T, screening (analytical quality of Nimis and Grutter 2010, four-phase,
         T_TA98 > 900 deg C, then the equilibrium screens; calculations/pt_fo2/screening.py)
      -> garnet Fe3+/SumFe iterated with the TA98 + garnet-Cr2O3 correction (Eq. 4;
         calculations/pt_fo2/cr_correction.py, coefficients in calculations/pt_fo2/params.py and reproduced by
         scripts/fit_cr_correction.py)                    -> screened_PT.csv, p-T-fO2.csv, cr_correction_flags.csv
      -> fO2 residuals and 1 sigma; Table 1                -> fO2_residuals.csv, fO2_uncertainty.json,
                                                              Table1_calibration_stats.csv, level_anchor.json
      -> Figs 3-10 and the numbers quoted in the text, all from these files (Fig. 5: fig5_experiments_cr.py).

Eq. 4 (v7_cr). Garnet Fe3+/SumFe is raised from 0.01 until
    T_TA98 - T_WZ07(Fe2+) <= a + b t + c (Cr2O3 - 5.364) t,   t = T_TA98 - 1207.02,
a = -102 (the Fig. 6 +-0.02 optimum; the joint WLS intercept is -92.77), b = -0.2850, c = 0.0673,
with T_TA98 held within 1124-1394 deg C and Cr2O3 within 1.3-9.8 wt%.
a is dT_match (the bold line of Fig. 4 and the adopted level of Figs 6 and 7).
"""

import argparse, json, shutil, subprocess, sys
from pathlib import Path
import pandas as pd

from load import REPO, PARAMS, calib
from calculations import run_full_pipeline
from calculations.pt_fo2.level_anchor import anchor_level, match_level
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


INPUT = Path(REPO) / "data" / "master" / "xenolith_geochem.csv"


def slave_match_level():
    """(dT_match at the pivot, pivot) in deg C from the screened Slave calibration samples. The screens and
    temperatures do not depend on the level, so they are taken from a pipeline run at any level."""
    table, _ = run_full_pipeline(INPUT, verbose=False)
    return match_level(table)


def calculate_at_match_level(level_c, pivot_c, out):
    """Run the pipeline at match level level_c about pivot pivot_c; write and return the full table."""
    assert not PARAMS.fe3.apply_empirical_correction, "the method applies no Eq. 5 (calculations/pt_fo2/params.py)"
    set_match_level(level_c, pivot_c)
    run_full_pipeline(INPUT, verbose=False,
                      write_csv=out / "p-T-fO2.csv", write_full_csv=out / "screened_PT.csv")
    table = pd.read_csv(out / "screened_PT.csv", encoding="utf-8-sig", low_memory=False)
    return table.rename(columns={table.columns[0]: "SampleID"})


def write_statistics(table, level_c, out):
    """Write the fO2 uncertainty, Table 1, the Slave anchor and the level sensitivity for the table."""
    from fo2_uncertainty import write_uncertainty
    from level_sensitivity import level_sensitivity
    (out / "level_anchor.json").write_text(json.dumps(anchor_level(table), indent=1, default=float))
    write_uncertainty(table, out)
    calib.calibration_stats(table).to_csv(out / "Table1_calibration_stats.csv", index=False)
    level_sensitivity(table, level_c, to_pipeline).to_csv(out / "level_sensitivity.csv", index=False)


def draw_figures(table):
    """Draw Figs 1-4 and 6-10 into the output folder; Figs 3, 4 and 6-10 from the table at the match level."""
    from unified import save
    import fig1, fig2, fig3
    from fig4_slave_level import draw_fig4
    from fig6_level_sensitivity import draw_fig6
    from fig7_fe3_match_level import draw_fig7
    from fig8_validation import draw_fig8
    from fig9_fo2_measured import draw_fig9
    from fig10_fo2_estimated import draw_fig10
    for n, draw in ((1, fig1.fig1), (2, fig2.fig2)):
        fig, px = draw()
        save(fig, f"Fig{n}_final.png", px)
    for n, draw in ((3, fig3.fig3), (4, draw_fig4), (6, draw_fig6), (7, draw_fig7), (8, draw_fig8), (9, draw_fig9), (10, draw_fig10)):
        # file and function numbers follow manuscript v7; the former Fig. 6 is retired_former_fig6_fe3_three_levels.py (not drawn)
        fig, px = draw(table)
        save(fig, f"Fig{n}_final.png", px)


def draw_fig5(out, pivot):
    """Fig. 5 (the 13 experiments and their weighted fit, equation written about the pivot) and the Online
    Resource 1 run/weight table, both from scripts/make_fig05_weightcircles.py."""
    tmp = out / "_fig5"
    run("make_fig05_weightcircles.py", "--outdir", tmp, "--pivot", repr(float(pivot)))
    (tmp / "Fig05_TA98_WZ07_measured_Fe3.png").replace(out / "Fig5_final.png")
    (tmp / "Fig05_weights.csv").replace(out / "Fig5_weights.csv")
    (tmp / "Fig05_TA98_WZ07_measured_Fe3_equation.json").replace(out / "Fig5_equation.json")
    shutil.rmtree(out / "Online_Resource_1", ignore_errors=True)
    (tmp / "Online_Resource_1").replace(out / "Online_Resource_1")
    shutil.rmtree(tmp)


def write_quoted_numbers(out):
    """Numbers quoted in the text and the three-way Table 1, from the output folder."""
    run("text_numbers.py", "--derived", out, "--out", out / "text_numbers.json")
    run("table1_three_way.py", "--derived", out, "--out", out / "Table1_three_way.md")


def run(script, *args):
    """Run one of the repository scripts in scripts/."""
    subprocess.run([sys.executable, str(SCRIPTS / script), *map(str, args)], check=True, stdout=subprocess.DEVNULL)


def write_cr_outputs(table, out):
    """Fig. 5 (natural samples split by garnet Cr2O3 + experiments), the range flags of the Cr2O3 correction
    (cr_correction_flags.csv) and its coefficients and counts (cr_correction.json)."""
    from unified import save
    from fig5_experiments_cr import draw_fig5_cr
    from calculations.pt_fo2.cr_correction import write_flags
    fig, px, st = draw_fig5_cr(table)
    save(fig, "Fig5_final.png", px)
    (out / "Fig5_cr_split.json").write_text(json.dumps(st, indent=1, default=float))
    fp = PARAMS.fe3
    counts = write_flags(table, out / "cr_correction_flags.csv", fp)
    scr = table[table["screen_pass"] == True]
    scr_counts = write_flags(scr, out / "_tmp_flags.csv", fp); (out / "_tmp_flags.csv").unlink()
    v7_rule, _ = match_level(table)
    info = dict(a=fp.cr_a, b=fp.cr_b, c=fp.cr_c, cr_ref=fp.cr_ref, pivot=fp.cr_pivot_c,
                cr_range_wt=[fp.cr_min_wt, fp.cr_max_wt], t_range_c=[fp.cr_tmin_c, fp.cr_tmax_c],
                v7_rule_level_with_cr_slope=v7_rule, flags_iterated=counts, flags_screened=scr_counts)
    (out / "cr_correction.json").write_text(json.dumps(info, indent=1, default=float))
    return info


if __name__ == "__main__":
    from unified import OUTD
    cli = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    cli.add_argument("--v7-eq4", action="store_true", help="revision v7 Eq. 4 (experimental slope, Slave median level)")
    cli.add_argument("--extend-slope", action="store_true", help="with --v7-eq4: extend the slope below 1162 deg C")
    a_ = cli.parse_args()
    out = Path(OUTD)
    (out / "figure_stats.json").unlink(missing_ok=True)
    if a_.v7_eq4:
        PARAMS.fe3.wz07_correction = "experimental"
        fit = set_wz07_correction(a_.extend_slope)
        level, pivot = slave_match_level()
    else:
        assert PARAMS.fe3.wz07_correction == "cr2o3" and PARAMS.fe3.wz07_offset_slope == 0 and PARAMS.fe3.wz07_offset_c == 0
        level, pivot = PARAMS.fe3.cr_a, PARAMS.fe3.cr_pivot_c
    table = calculate_at_match_level(level, pivot, out)
    write_statistics(table, level, out)
    draw_figures(table)
    draw_fig5(out, pivot)          # Online Resource 1 run table and the weighted experimental fit (Fig5_equation.json)
    if not a_.v7_eq4:
        info = write_cr_outputs(table, out)   # replaces Fig5_final.png with the Cr2O3 figure
    write_quoted_numbers(out)
    if a_.v7_eq4:
        shape = "extended" if a_.extend_slope else "held flat"
        print(f"v7 Eq. 4: c(T) = {fit.b0:.1f} {fit.b1:+.3f} (T - {fit.t_ref:.0f}), {shape} below the experiments; "
              f"match level {level:.1f} deg C at {pivot:.1f} deg C; outputs in {out}")
    else:
        print(f"Eq. 4 (Cr2O3): {level:.2f} {PARAMS.fe3.cr_b:+.4f} t {PARAMS.fe3.cr_c:+.4f} (Cr2O3 - {PARAMS.fe3.cr_ref:.3f}) t, "
              f"t = T_TA98 - {pivot:.2f}; flagged {info['flags_screened']}; outputs in {out}")
