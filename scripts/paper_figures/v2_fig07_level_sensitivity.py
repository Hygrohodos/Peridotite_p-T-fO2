"""Fig. 7 (revision v3): sensitivity of the Fe3+ agreement to the flat match level on the 25 Slave samples
outside the error groups (no WZ07 correction, no Eq. 5). The adopted level is ANCHORED to the median
natural discrepancy of all 30 Slave calibration samples (level_anchor.py), not chosen from this sweep;
the +-0.02 optimum of the sweep is marked for comparison only."""
from v2common import *
import json
def fig6():
    df=load(); base,meas,blue,red,rem=mf._error_groups(df)
    s=calib.dt_sweep(df,0.0,0.0,dts=np.arange(-250,0.5,1.0),mask=rem)
    a=anchor()
    lev=FP.dt_match_c; assert abs(lev-a['dt_level_C'])<0.05
    sa=calib.dt_sweep(df,0.0,0.0,dts=[lev],mask=rem).iloc[0]
    opt=calib.select_dt(s); ro=s[s.dt==opt].iloc[0]; ties=s.loc[s.p02==s.p02.max(),'dt']
    s.to_csv(os.path.join(OUTD,'Fig7_level_sweep.csv'),index=False)
    W,H=119,84; fig=newfig(W,H)
    ax=axmm(fig,W,H,17,12.5,87,61.5); ax2=ax.twinx()
    ax.set_zorder(ax2.get_zorder()+1); ax.patch.set_visible(False)
    ax.axhline(0,color='0.55',lw=LW_THIN,zorder=1)
    ax.axvline(lev,color='k',lw=LW_ONE,zorder=1)
    ax.axvline(opt,color='0.45',lw=LW_THIN,ls=(0,(1.5,1.5)),zorder=1)
    ax.plot(s.dt,s['mean'],color='k',lw=LW_MAIN,zorder=5)
    ax.set_xlim(-250,0); ax.set_ylim(-0.10,0.10)
    ax.set_xlabel('\u0394T match level, T$_\\mathrm{TA98}$ \u2212 T$_\\mathrm{WZ07}$ (\u00b0C)'); ax.set_ylabel('mean Fe$^{3+}$/\u03a3Fe (est \u2212 meas)',labelpad=3)
    yt=[-0.10,-0.05,0,0.05,0.10]; ax.set_yticks(yt); ax.set_yticklabels([M(f'{v:.2f}') for v in yt])
    ax.xaxis.set_major_locator(MultipleLocator(50)); ax.xaxis.set_minor_locator(MultipleLocator(25)); ax.yaxis.set_minor_locator(MultipleLocator(0.025))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v,_:M(f'{v:.0f}')))
    ax2.plot(s.dt,100*s.p05,color=G05,lw=LW_MAIN,ls=(0,(4,1.6)),zorder=3)
    ax2.plot(s.dt,100*s.p02,color=G02,lw=LW_MAIN,zorder=4)
    ax2.set_ylim(0,105); ax2.set_yticks(range(0,101,20)); ax2.yaxis.set_minor_locator(MultipleLocator(10))
    ax2.set_ylabel('Success rate (%)',labelpad=4,color=G02_T)
    ax2.tick_params(axis='y',colors=G02_T,which='both'); ax2.spines['right'].set_color(G02_T)
    ax.tick_params(axis='y',right=False,which='both'); ax.spines['left'].set_color('k')
    ax2.plot([lev],[100*sa.p02],marker='o',ms=5.5,mfc='white',mec=G02_T,mew=LW_RING,zorder=8,clip_on=False)
    ax2.annotate(M(f'adopted level\n{lev:.0f} \u00b0C:\n{100*sa.p02:.0f} % \u00b10.02\n{100*sa.p05:.0f} % \u00b10.05'),xy=(lev,100*sa.p02),xytext=(-233,47),
                 fontsize=FS_ANNOT,color=G02_T,fontweight='bold',ha='left',va='center',zorder=9,linespacing=1.15,
                 path_effects=halo(2.2),arrowprops=dict(arrowstyle='-',color=G02_T,lw=LW_THIN,shrinkA=2,shrinkB=4))
    hs=[Line2D([],[],color='k',lw=LW_MAIN,label='mean (est \u2212 meas)'),Line2D([],[],color=G02,lw=LW_MAIN,label='% within \u00b10.02'),
        Line2D([],[],color=G05,lw=LW_MAIN,ls=(0,(4,1.6)),label='% within \u00b10.05'),
        Line2D([],[],color='k',lw=LW_ONE,label=M(f'adopted level ({lev:.0f} \u00b0C)')),
        Line2D([],[],color='0.45',lw=LW_THIN,ls=(0,(1.5,1.5)),label=M(f'\u00b10.02 optimum ({opt:.0f} \u00b0C)'))]
    ax.legend(handles=hs,loc='lower center',bbox_to_anchor=(0.5,1.01),ncol=3,frameon=False,fontsize=FS_SMALL,handletextpad=0.4,labelspacing=0.3,columnspacing=0.8,handlelength=1.8)
    st=dict(level=float(lev),opt=float(opt),n=int(sa.n),p02=float(sa.p02),p05=float(sa.p05),mean=float(sa['mean']),
            tuned_opt=int(opt),tuned_p02=float(ro.p02),tuned_p05=float(ro.p05),p02_max_range=[float(ties.min()),float(ties.max())])
    log('fig07_sweep',st); print('Fig7(sweep)',st)
    return fig,(2811,1984)
if __name__=='__main__':
    f,px=fig6(); save(f,'Fig7_final.png',px)
