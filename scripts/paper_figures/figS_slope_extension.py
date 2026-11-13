"""Supplementary figure and table: the experimental slope held flat below 1162 deg C (paper) or extended.

Usage:  PAPER_FIG_OUT=<dir> python3 figS_slope_extension.py <held-flat output> <extended output>
where the two inputs are output folders of make_paper_figures.py without and with --extend-slope.
Draws the 30 Slave calibration samples (slope-removed TA98 - WZ07 against TA98, Fig5_per_sample.csv) with
both match levels, and tabulates the statistics of both runs.
"""
import sys
from common import *

CASES = (('held flat', SLAVE), ('extended', '#2563eb'))
SETS = (('Slave 30', 'Slave calibration'), ('Slave 25', 'Slave, error groups excluded'),
        ('Validation 27', 'Siberia + Kaapvaal validation'), ('All 57', 'Calibration + validation'),
        ('52, error groups excluded', 'Combined, error groups excluded'))


def summary(d, case):
    """One row of statistics for one output folder."""
    tab = pd.read_csv(os.path.join(d, 'Table1_calibration_stats.csv')).set_index('dataset')
    tn = json.load(open(os.path.join(d, 'text_numbers.json')))
    fs = json.load(open(os.path.join(d, 'figure_stats.json')))['fig5']
    unc = json.load(open(os.path.join(d, 'fO2_uncertainty.json')))
    row = {'slope below 1162 C': case, 'match level at 1208.5 C (C)': round(fs['level'], 1)}
    for lab, ds in SETS:
        x = tab.loc[ds]
        row[f'{lab}: within 0.02 (%)'] = round(100 * x.p02); row[f'{lab}: within 0.05 (%)'] = round(100 * x.p05)
        row[f'{lab}: 1 sigma Fe3+/SumFe'] = round(x.sigma, 3); row[f'{lab}: 1 sigma log fO2'] = round(x.sigma_log, 2)
    row.update({'at 0.01 floor (of 450)': tn['n_floor'], 'at floor, TA98 < 1100 C': tn['floor_below_1100'],
                'Slave 25 residual-T trend p': round(tn['resid_vs_T_slave25']['p'], 3),
                'Slave 25 residual-T trend p, TA98 >= 1130 C': round(tn['resid_vs_T_slave25_T_ge_1130']['p'], 3),
                '52 residual-T trend p': float(f"{tn['resid_vs_T_excl52']['p']:.2g}"),
                'median dFMQ, all': round(tn['median_dFMQ_all'], 2),
                'fO2 1 sigma reduced / oxidised': f"{unc['reduced']['sigma_log']:.2f} / {unc['oxidised']['sigma_log']:.2f}"})
    return row, fs['level']


def draw(dirs):
    """Supplementary figure from the two output folders."""
    W, H = 119, 84; fig = newfig(W, H)
    a = axmm(fig, W, H, 19, 12.5, 96, 68)
    a.axvline(1162, color='0.5', lw=LW_THIN, ls=':', zorder=1)
    a.text(1166, -335, '1162 \u00b0C', fontsize=FS_SMALL, color='0.4', va='bottom')
    rows, keys = [], []
    for d, (case, col) in zip(dirs, CASES):
        row, lev = summary(d, case); rows.append(row)
        u = pd.read_csv(os.path.join(d, 'Fig5_per_sample.csv'))
        flat = case == 'held flat'
        a.axhline(lev, color=col, lw=LW_MAIN, ls='-' if flat else '--', zorder=2)
        a.scatter(u.T_TA98, u.dT_corr, s=S_STD, marker='o' if flat else 's', facecolors=col if flat else 'none',
                  edgecolors='k' if flat else col, linewidths=EW_STD if flat else 1.0, zorder=4 if flat else 3)
        keys.append(Line2D([], [], color=col, ls='-' if flat else '--', lw=LW_REF, marker='o' if flat else 's', ms=MS_KEY,
                           mfc=col if flat else 'none', mec='k' if flat else col, label=M(f'{case} below 1162 \u00b0C (level {lev:.0f} \u00b0C)')))
    a.set_xlim(950, 1350); a.set_ylim(-350, 50)
    a.xaxis.set_major_locator(MultipleLocator(100)); a.xaxis.set_minor_locator(MultipleLocator(25))
    a.yaxis.set_major_locator(MultipleLocator(100)); a.yaxis.set_minor_locator(MultipleLocator(25))
    a.yaxis.set_major_formatter(FuncFormatter(lambda v, _: M(f'{v:.0f}')))
    a.set_xlabel('TA98 temperature (\u00b0C)'); a.set_ylabel('$T_\\mathrm{TA98}-T_\\mathrm{WZ07}$, slope removed (\u00b0C)')
    a.legend(handles=keys, loc='upper left', frameon=False, fontsize=FS_KEY, handletextpad=0.4, labelspacing=0.35, borderaxespad=0.4)
    pd.DataFrame(rows).set_index('slope below 1162 C').T.to_csv(os.path.join(OUTD, 'Supplementary_slope_extension_table.csv'))
    return fig, (2811, 1984)


if __name__ == '__main__':
    f, px = draw(sys.argv[1:3]); save(f, 'Supplementary_slope_extension.png', px)
