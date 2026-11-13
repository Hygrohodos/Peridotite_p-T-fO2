"""Fig. 2: effect of garnet Fe3+ on the olivine-garnet (WZ07) temperature and on KD (Eq. 3) for
Jericho peridotite 21-2, and the idea of the calculation: Fe3+/SumFe at which WZ07 meets TA98.

Input: data/fig2_curves.csv (x = garnet Fe3+/SumFe removed from the Fe2+ budget; KD and T_WZ07 of 21-2).
"""
from unified import *

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'fig2_curves.csv')
TA98_21_2 = 1288.893081993608        # TA98 temperature of 21-2 (deg C); independent of garnet Fe3+
COLOUR = TA98C                       # sample colour, as in the other figures
COLOUR_TEXT = '#C25600'


def fig2():
    """Olivine-garnet temperature (top) and KD (bottom) against garnet Fe3+/SumFe for sample 21-2."""
    d = pd.read_csv(DATA)
    x, T_wz, KD = d.x.values, d.T_WZ07_21_2.values, d.KD_21_2.values
    x_match = float(np.interp(0, -d.dT_21_2.values, x))          # where TA98 - WZ07 reaches 0
    T0, T15 = float(T_wz[0]), float(np.interp(0.15, x, T_wz))

    W, H = 84, 100
    fig, (t, k) = plt.subplots(2, 1, figsize=(W * MM, H * MM), sharex=True,
                               gridspec_kw=dict(hspace=0.10, height_ratios=[1.6, 1]))
    fig.subplots_adjust(left=0.20, right=0.97, top=0.97, bottom=0.115)

    t.fill_between(x, T_wz, TA98_21_2, color=COLOUR, alpha=0.16, lw=0)
    t.axhline(T0, color=COLOUR, lw=LW_THIN, ls=':', zorder=1)
    t.axhline(TA98_21_2, color='k', lw=LW_MAIN, ls=(0, (5, 2)), zorder=2)
    t.plot(x, T_wz, color=COLOUR, lw=LW_MAIN, zorder=4)
    t.set_ylim(1130, 1480); t.set_xlim(0, 0.2)
    t.set_ylabel('Olivine\u2013garnet\ntemperature (\u00b0C)', linespacing=1.1)
    t.annotate('', xy=(0.15, T15), xytext=(0.15, T0),
               arrowprops=dict(arrowstyle='<->', color=COLOUR, lw=LW_ONE, shrinkA=0, shrinkB=0))
    t.text(0.155, T0 + 0.28 * (T15 - T0), f'+{T15 - T0:.0f} \u00b0C', color=COLOUR_TEXT,
           fontsize=FS_ANNOT, fontweight='bold', ha='left', va='center')
    t.text(0.003, TA98_21_2 + 7, f'TA98 = {TA98_21_2:.0f} \u00b0C', color='k', ha='left', va='bottom',
           fontsize=FS_ANNOT, fontweight='bold')
    t.plot([x_match], [TA98_21_2], 'o', ms=8, mfc='white', mec='k', mew=LW_RING, zorder=6)
    t.annotate(f'match at {x_match:.3f}', xy=(x_match, TA98_21_2), xytext=(0.10, 1440),
               fontsize=FS_ANNOT, fontweight='bold', ha='center', va='center',
               arrowprops=dict(arrowstyle='-|>', color='k', lw=LW_ONE, shrinkA=6, shrinkB=5))

    k.plot(x, KD, color=COLOUR, lw=LW_MAIN)
    k.set_xticks([0, .05, .1, .15, .2]); k.set_xlabel('Garnet Fe$^{3+}$/\u03a3Fe')
    k.set_ylabel(r'$K_\mathrm{D}$')
    rng = KD.max() - KD.min(); k.set_ylim(KD.min() - 0.15 * rng, KD.max() + 0.15 * rng)

    for ax in (t, k):
        ax.grid(axis='y', color='0.9', lw=0.5, zorder=0); ax.set_axisbelow(True); minor(ax)
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    fig.align_ylabels((t, k))
    return fig, (int(W * MM * 600), int(H * MM * 600))


if __name__ == '__main__':
    f, px = fig2(); save(f, 'Fig2_final.png', px)
