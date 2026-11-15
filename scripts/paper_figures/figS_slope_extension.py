"""Online Resource 2 figure and table: the experimental slope held flat below 1162 deg C (paper) or extended.

Usage:  PAPER_FIG_OUT=<dir> python3 figS_slope_extension.py <held-flat output> <extended output>
where the two inputs are output folders of make_paper_figures.py without and with --extend-slope.
Draws the screened Slave calibration samples (slope-removed TA98 - WZ07 against TA98, Fig4_per_sample.csv) with
both match levels, and tabulates the statistics of both runs.
"""
import sys
from common import *

CASES = (('held flat', SLAVE), ('extended', '0.2'))   # extended: dark grey open squares (blue squares mean Siberian elsewhere)
SETS = ('Slave calibration', 'Slave, error groups excluded', 'Siberia + Kaapvaal validation', 'Calibration + validation',
        'Combined, error groups excluded')


def summary(d, case):
    """One row of statistics for one output folder."""
    tab = pd.read_csv(os.path.join(d, 'Table1_calibration_stats.csv')).set_index('dataset')
    tn = json.load(open(os.path.join(d, 'text_numbers.json')))
    fs = json.load(open(os.path.join(d, 'figure_stats.json')))['fig4']
    unc = json.load(open(os.path.join(d, 'fO2_uncertainty.json')))
    row = {'slope below 1162 C': case, 'pivot (C)': round(fs['pivot'], 1), 'match level at the pivot (C)': round(fs['level'], 1)}
    for ds in SETS:
        x = tab.loc[ds]
        row[f'{ds}: n'] = int(x.n)
        row[f'{ds}: within 0.02 (%)'] = round(100 * x.p02); row[f'{ds}: within 0.05 (%)'] = round(100 * x.p05)
        row[f'{ds}: 1 sigma Fe3+/SumFe'] = round(x.sigma, 3); row[f'{ds}: 1 sigma log fO2'] = round(x.sigma_log, 2)
    row.update({'screened samples': tn['n_screened'], 'at 0.01 floor': tn['n_floor'], 'at floor, TA98 < 1100 C': tn['floor_below_1100'],
                'Slave, error groups excluded: residual-T trend p': round(tn['resid_vs_T_slave_excl_eg']['p'], 3),
                'Slave, error groups excluded: residual-T trend p, TA98 >= 1130 C': round(tn['resid_vs_T_slave_excl_eg_T_ge_1130']['p'], 3),
                'Combined, error groups excluded: residual-T trend p': float(f"{tn['resid_vs_T_excl_eg']['p']:.2g}"),
                'median dFMQ, all': round(tn['median_dFMQ_all'], 2),
                'fO2 1 sigma reduced / oxidised': f"{unc['reduced']['sigma_log']:.2f} / {unc['oxidised']['sigma_log']:.2f}"})
    return row, fs['level']


def draw(dirs):
    """Supplementary figure from the two output folders."""
    W, H = COL1, 88; fig = newfig(W, H)          # house style: one column, axes as Fig. 4
    a = axmm(fig, W, H, 15.5, 11.0, 66, 66)
    a.axvline(1162, color='0.5', lw=LW_THIN, ls=(0, (1, 1.5)), zorder=1)
    a.text(1166, -335, '1162 \u00b0C', fontsize=FS_ANNOT, color='0.35', va='bottom')
    rows, keys = [], []
    for d, (case, col) in zip(dirs, CASES):
        row, lev = summary(d, case); rows.append(row)
        u = pd.read_csv(os.path.join(d, 'Fig4_per_sample.csv'))
        flat = case == 'held flat'
        a.axhline(lev, color=col, lw=LW_MAIN, ls='-' if flat else '--', zorder=2)
        a.scatter(u.T_TA98, u.dT_corr, s=S_STD, marker='o' if flat else 's', facecolors=col if flat else 'none',
                  edgecolors='k' if flat else col, linewidths=EW_STD if flat else 0.8, zorder=4 if flat else 3)
        keys.append(Line2D([], [], color=col, ls='-' if flat else '--', lw=LW_REF, marker='o' if flat else 's', ms=MS_KEY,
                           mfc=col if flat else 'none', mec='k' if flat else col, label=M(f'slope {case} below 1162 \u00b0C (level {lev:.0f} \u00b0C)')))
    a.set_xlim(950, 1350); a.set_ylim(-350, 50)
    a.xaxis.set_major_locator(MultipleLocator(100)); a.xaxis.set_minor_locator(MultipleLocator(25))
    a.yaxis.set_major_locator(MultipleLocator(100)); a.yaxis.set_minor_locator(MultipleLocator(25))
    a.yaxis.set_major_formatter(FuncFormatter(lambda v, _: M(f'{v:.0f}')))
    a.set_xlabel(XT_TA98); a.set_ylabel('$T_\\mathrm{TA98}-T_\\mathrm{WZ07}$, slope removed (\u00b0C)')
    a.legend(handles=keys, loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=1, fontsize=FS_KEY, borderaxespad=0.5, handlelength=2.2)
    pd.DataFrame(rows).set_index('slope below 1162 C').T.to_csv(os.path.join(OUTD, 'Online_Resource_2_slope_extension_table.csv'))
    return fig, pxsize(W, H)


if __name__ == '__main__':
    f, px = draw(sys.argv[1:3]); save(f, 'Online_Resource_2_slope_extension.png', px)
