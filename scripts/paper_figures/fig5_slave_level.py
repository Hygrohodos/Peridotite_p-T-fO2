"""Fig. 5: TA98 - WZ07 against TA98 for the screened Slave calibration samples, with WZ07 calculated with
the measured garnet Fe3+/SumFe. Open: as calculated; filled: after removing the experimental slope of Fig. 4
about the pivot, c(T) - c(pivot). The bold line is the match level dT_match (Eq. 4) at the pivot
(calculations/pt_fo2/level_anchor.py). Per-sample values are recomputed from the pipeline table.
"""
from common import *
from calculations.pt_fo2.fe3_correction import _wz07corr_fe2only
from calculations.pt_fo2.level_anchor import anchor_level


def per_sample(df):
    """TA98, WZ07 at measured Fe3+, and TA98 - WZ07 before and after removing the slope, for the screened Slave samples."""
    base, meas = mf._error_groups(df)[:2]
    use = base & meas.notna()
    u = df[use].copy(); m = meas[use].to_numpy(float)
    num = lambda c: pd.to_numeric(u[c], errors='coerce').to_numpy(float)
    u['T_TA98'] = num('T_TA98')
    u['WZm'] = _wz07corr_fe2only(num('olivine_Mg'), num('olivine_Fe'), num('garnet_Mg'), num('garnet_Fe') * (1.0 - m),
                                 num('garnet_Ca'), num('garnet_Mn'), num('P_CA91'), offset_c=0.0)
    u['dT_raw'] = u.T_TA98 - u.WZm
    u['dT_corr'] = u.dT_raw - slope_term(u.T_TA98)
    return u


def trend(x, y):
    """Least-squares slope (deg C per deg C), its standard error and p."""
    from scipy import stats
    r = stats.linregress(x, y)
    return dict(slope=float(r.slope), se=float(r.stderr), p=float(r.pvalue), n=int(len(x)))


def draw_fig5(df):
    """Fig. 5 from the pipeline table df; the bold line is the match level used by Figs 6-11."""
    u = per_sample(df); lev = match_level(); A = anchor_level(df)
    W, H = 119, 84; fig = newfig(W, H)
    a = axmm(fig, W, H, 19, 12.5, 96, 62)
    a.axhline(0, color='0.6', lw=LW_THIN, zorder=1)
    a.axhline(lev, color='k', lw=2.6, zorder=2)
    a.plot([pivot_t()], [lev], marker='D', ms=4.5, mfc='white', mec='k', mew=1.2, zorder=7, clip_on=False)
    a.annotate(M(f'\u0394T$_\\mathrm{{match}}$ = {lev:.0f} \u00b0C at {pivot_t():.0f} \u00b0C'), xy=(pivot_t(), lev), xytext=(1020, -280),
               fontsize=FS_ANNOT + 1, fontweight='bold', ha='left', va='center', zorder=6,
               arrowprops=dict(arrowstyle='-', color='k', lw=LW_THIN, shrinkA=2, shrinkB=4))
    a.scatter(u.T_TA98, u.dT_raw, s=S_STD * 0.8, facecolors='none', edgecolors=SLAVE, linewidths=0.9, zorder=3)
    a.scatter(u.T_TA98, u.dT_corr, s=S_STD, c=SLAVE, edgecolors='k', linewidths=EW_STD, zorder=4)
    a.set_xlim(950, 1350); a.set_ylim(-350, 50)
    a.xaxis.set_major_locator(MultipleLocator(100)); a.xaxis.set_minor_locator(MultipleLocator(25))
    a.yaxis.set_major_locator(MultipleLocator(100)); a.yaxis.set_minor_locator(MultipleLocator(25))
    a.yaxis.set_major_formatter(FuncFormatter(lambda v, _: M(f'{v:.0f}')))
    a.set_xlabel('TA98 temperature (\u00b0C)'); a.set_ylabel('$T_\\mathrm{TA98}-T_\\mathrm{WZ07}$ (\u00b0C)')
    key = [Line2D([], [], marker='o', ls='', mfc=SLAVE, mec='k', mew=EW_STD, ms=MS_KEY, label='slope removed'),
           Line2D([], [], marker='o', ls='', mfc='none', mec=SLAVE, mew=0.9, ms=MS_KEY, label='as calculated')]
    a.legend(handles=key, loc='upper left', frameon=False, fontsize=FS_KEY, handletextpad=0.3, labelspacing=0.35, borderaxespad=0.4)
    hi = (u.T_TA98 >= 1130).to_numpy()
    st = dict(n=int(len(u)), level=float(lev), pivot=pivot_t(), median_raw=float(u.dT_raw.median()), median_slope_removed=float(u.dT_corr.median()),
              median_T=float(u.T_TA98.median()), ci95_bootstrap_D=A['ci95_bootstrap_D_C'], n_raw_lt_0=int((u.dT_raw < 0).sum()),
              trend_raw_all=trend(u.T_TA98, u.dT_raw), trend_corr_all=trend(u.T_TA98, u.dT_corr),
              trend_raw_T_ge_1130=trend(u.T_TA98[hi], u.dT_raw[hi]), trend_corr_T_ge_1130=trend(u.T_TA98[hi], u.dT_corr[hi]),
              n_outside_axes=int(((u.dT_raw < -350) | (u.dT_corr < -350) | (u.dT_raw > 50) | (u.T_TA98 < 950) | (u.T_TA98 > 1350)).sum()))
    log('fig5', st); print('Fig5', st)
    u[['SampleID', 'T_TA98', 'WZm', 'dT_raw', 'dT_corr']].to_csv(os.path.join(OUTD, 'Fig5_per_sample.csv'), index=False)
    return fig, (2811, 1984)


if __name__ == '__main__':
    f, px = draw_fig5(load()); save(f, 'Fig5_final.png', px)
