"""Fig. 4: TA98 - WZ07 against TA98 for the screened Slave calibration samples, with WZ07 calculated with
the measured garnet Fe3+/SumFe. Open: as calculated; filled: after removing the experimental slope of Fig. 5
about the pivot, c(T) - c(pivot). The bold line is the match level dT_match (Eq. 4) at the pivot
(calculations/pt_fo2/level_anchor.py). Blue and red rings mark the two error groups of Fig. 7 on the filled
(slope-removed) points, in the ring style of Fig. 7 (fe3rows.rings). Per-sample values are recomputed from the
pipeline table.
"""
from common import *
from calculations.pt_fo2.fe3_correction import _wz07corr_fe2only
from calculations.pt_fo2.level_anchor import anchor_level
from fe3rows import rings


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
    u['dT_corr'] = u.dT_raw - slope_term(u.T_TA98, pd.to_numeric(u['garnet_Cr2O3'], errors='coerce').to_numpy(float))
    return u


def error_group_masks(df, u):
    """Boolean masks over the rows of u (per_sample) for the blue and red error groups of Fig. 7."""
    blue, red = mf._error_groups(df)[2:4]
    sid = u.SampleID.astype(str)
    return (sid.isin(set(df.SampleID[blue].astype(str))).to_numpy(),
            sid.isin(set(df.SampleID[red].astype(str))).to_numpy())


def trend(x, y):
    """Least-squares slope (deg C per deg C), its standard error and p."""
    from scipy import stats
    r = stats.linregress(x, y)
    return dict(slope=float(r.slope), se=float(r.stderr), p=float(r.pvalue), n=int(len(x)))


def draw_fig4(df):
    """Fig. 4 from the pipeline table df; the bold line is the match level used by Figs 6-10."""
    u = per_sample(df); lev = match_level(); A = anchor_level(df)
    W, H = COL1, 80; fig = newfig(W, H)
    a = axmm(fig, W, H, 15.5, 11.0, 66, 66)
    a.axhline(0, color='0.6', lw=LW_THIN, zorder=1)
    a.axhline(lev, color='k', lw=2.2, zorder=2)
    a.plot([pivot_t()], [lev], marker='D', ms=4.5, mfc='white', mec='k', mew=1.2, zorder=7, clip_on=False)
    a.annotate(M(f'\u0394T$_\\mathrm{{match}}$ = {lev:.0f} \u00b0C at {pivot_t():.0f} \u00b0C'), xy=(pivot_t(), lev), xytext=(958, 28),
               fontsize=FS_ANNOT, fontweight='bold', ha='left', va='center', zorder=6,
               arrowprops=dict(arrowstyle='-|>', mutation_scale=7, color='k', lw=LW_THIN, shrinkA=2, shrinkB=4, relpos=(0.985, 0.0),
                               connectionstyle='angle,angleA=90,angleB=135,rad=0'))
    a.scatter(u.T_TA98, u.dT_raw, s=S_STD * 0.8, facecolors='none', edgecolors=SLAVE_T, linewidths=0.8, zorder=3)
    a.scatter(u.T_TA98, u.dT_corr, s=S_STD, c=SLAVE, edgecolors='k', linewidths=EW_STD, zorder=4)
    blue, red = error_group_masks(df, u)
    rings(a, u.T_TA98.to_numpy(), u.dT_corr.to_numpy(), blue, red)
    a.set_xlim(950, 1350); a.set_ylim(-350, 50)
    a.xaxis.set_major_locator(MultipleLocator(100)); a.xaxis.set_minor_locator(MultipleLocator(25))
    a.yaxis.set_major_locator(MultipleLocator(100)); a.yaxis.set_minor_locator(MultipleLocator(25))
    a.yaxis.set_major_formatter(FuncFormatter(lambda v, _: M(f'{v:.0f}')))
    a.set_xlabel(XT_TA98); a.set_ylabel(DT_LABEL)
    key = [Line2D([], [], marker='o', ls='', mfc=SLAVE, mec='k', mew=EW_STD, ms=MS_KEY, label='slope removed'),
           Line2D([], [], marker='o', ls='', mfc='none', mec=SLAVE_T, mew=0.8, ms=MS_KEY, label='as calculated'),
           Line2D([], [], marker='o', ls='', mfc='none', mec=BLUE_RING, mew=LW_RING, ms=MS_KEY + 1.5, label='blue error group'),
           Line2D([], [], marker='o', ls='', mfc='none', mec=RED_RING, mew=LW_RING, ms=MS_KEY + 1.5, label='red error group')]
    a.legend(handles=key, loc='lower center', bbox_to_anchor=(0.43, 0.0), fontsize=FS_KEY, borderaxespad=0.3,
             ncol=2, columnspacing=1.0)
    hi = (u.T_TA98 >= 1130).to_numpy()
    st = dict(n=int(len(u)), level=float(lev), pivot=pivot_t(), median_raw=float(u.dT_raw.median()), median_slope_removed=float(u.dT_corr.median()),
              median_T=float(u.T_TA98.median()), ci95_bootstrap_D=A['ci95_bootstrap_D_C'], n_raw_lt_0=int((u.dT_raw < 0).sum()),
              trend_raw_all=trend(u.T_TA98, u.dT_raw), trend_corr_all=trend(u.T_TA98, u.dT_corr),
              trend_raw_T_ge_1130=trend(u.T_TA98[hi], u.dT_raw[hi]), trend_corr_T_ge_1130=trend(u.T_TA98[hi], u.dT_corr[hi]),
              error_group_blue=sorted(u.SampleID[blue].astype(str)), error_group_red=sorted(u.SampleID[red].astype(str)),
              n_outside_axes=int(((u.dT_raw < -350) | (u.dT_corr < -350) | (u.dT_raw > 50) | (u.T_TA98 < 950) | (u.T_TA98 > 1350)).sum()))
    log('fig4', st); print('Fig4', st)
    u[['SampleID', 'T_TA98', 'WZm', 'dT_raw', 'dT_corr']].to_csv(os.path.join(OUTD, 'Fig4_per_sample.csv'), index=False)
    return fig, pxsize(W, H)


if __name__ == '__main__':
    f, px = draw_fig4(load()); save(f, 'Fig4_final.png', px)
