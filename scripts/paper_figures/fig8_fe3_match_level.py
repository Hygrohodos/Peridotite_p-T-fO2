"""Fig. 8: estimated against measured garnet Fe3+/SumFe for the screened Slave calibration samples at the match level,
and the residuals against temperature, with the statistics for all samples and for those outside the error groups."""
from fe3rows import *
from common import *


def draw_fig8(df):
    """Fig. 8 from the pipeline table df, at the match level."""
    base, meas, blue, red, rem = mf._error_groups(df); T = pd.to_numeric(df['T_TA98'], errors='coerce')
    lev = match_level(); est = fe3_at(df, lev); ok = base & est.notna() & T.notna()
    W, H = 174, 72; fig = newfig(W, H); b0 = 12.4; ah = H - b0 - 8.2
    aL = axmm(fig, W, H, 17.4, b0, 64.5, ah); aR = axmm(fig, W, H, 102.4, b0, 64.5, ah)
    row(aL, aR, df, meas, T, est, ok, blue=blue, red=red)
    aL.set_ylim(0, 0.225); fmt(aL, np.arange(0, 0.351, 0.05), np.arange(0, 0.2251, 0.025), '%.2f', '%.3f')
    letter(aL, 'a.'); letter(aR, 'b.')
    da = est[ok] - meas[ok]; dr = est[rem & ok] - meas[rem & ok]
    aL.text(0.98, 0.04, f"all n = {int(ok.sum())}\n\u00b10.02: {pct(da, 0.02) * 100:.0f}%\n\u00b10.05: {pct(da, 0.05) * 100:.0f}%\n"
            f"remainder n = {int((rem & ok).sum())}\n\u00b10.02: {pct(dr, 0.02) * 100:.0f}%\n\u00b10.05: {pct(dr, 0.05) * 100:.0f}%",
            transform=aL.transAxes, fontsize=FS_ANNOT, ha='right', va='bottom', zorder=20, path_effects=halo(2.0))
    log('fig8', dict(dt=lev, n=int(ok.sum()), p02=pct(da, 0.02), p05=pct(da, 0.05), mean=float(da.mean()),
                     n_rem=int((rem & ok).sum()), p02_rem=pct(dr, 0.02), p05_rem=pct(dr, 0.05),
                     blue=sorted(df.SampleID[blue].astype(str)), red=sorted(df.SampleID[red].astype(str))))
    fig.legend(handles=key_handles(), loc='upper center', bbox_to_anchor=(0.5, 1 - 0.8 / H), ncol=5, frameon=False,
               fontsize=FS_KEY, handletextpad=0.4, columnspacing=2.2, handlelength=1.6, borderaxespad=0)
    return fig, (4110, 1700)


if __name__ == '__main__':
    f, px = draw_fig8(load()); save(f, 'Fig8_final.png', px)
