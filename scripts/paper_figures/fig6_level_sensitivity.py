"""Fig. 6: sensitivity of the Fe3+ agreement to the match level (at the pivot, with the experimental slope)
on the screened Slave samples outside the error groups. The adopted level (make_paper_figures.py) is the
+-0.02 optimum of this sweep (ties broken by smallest |mean|)."""
from common import *
def draw_fig6(df):
    """Fig. 6 from the pipeline table df: sweep of the match level, adopted level marked."""
    base,meas,blue,red,rem=mf._error_groups(df)
    lv=np.arange(-250,0.5,1.0)
    s=calib.dt_sweep(df,0.0,0.0,dts=[to_pipeline(v) for v in lv],mask=rem)
    s['dt']=lv.round(0)   # level at the pivot (dt_sweep truncates the pipeline value to int)
    lev=match_level()
    sa=calib.dt_sweep(df,0.0,0.0,dts=[to_pipeline(lev)],mask=rem).iloc[0]
    opt=calib.select_dt(s); ro=s[s.dt==opt].iloc[0]; ties=s.loc[s.p02==s.p02.max(),'dt']
    s.to_csv(os.path.join(OUTD,'Fig6_level_sweep.csv'),index=False)
    W,H=COL1,88; fig=newfig(W,H)
    ax=axmm(fig,W,H,15.5,11.0,57.0,60.0); ax2=ax.twinx()
    ax.set_zorder(ax2.get_zorder()+1); ax.patch.set_visible(False)
    ax.axhline(0,color='0.55',lw=LW_THIN,zorder=1)
    ax.axvline(lev,ymax=100/125,color='k',lw=LW_ONE,zorder=1)   # stop at 100 %, below the read-out
    ax.axvline(opt,ymax=100/125,color='0.45',lw=LW_THIN,ls=(0,(1.5,1.5)),zorder=1)
    ax.plot(s.dt,s['mean'],color='k',lw=LW_MAIN,zorder=5)
    ax.set_xlim(-250,0); ax.set_ylim(-0.10,0.10)
    ax.set_xlabel(M(f'\u0394T$_\\mathrm{{match}}$ at {pivot_t():.0f} \u00b0C (\u00b0C)')); ax.set_ylabel('Mean Fe$^{3{+}}$/\u03a3Fe, estimated \u2212 measured',labelpad=3)
    yt=[-0.10,-0.05,0,0.05,0.10]; ax.set_yticks(yt); ax.set_yticklabels([M(f'{v:.2f}') for v in yt])
    ax.xaxis.set_major_locator(MultipleLocator(50)); ax.xaxis.set_minor_locator(MultipleLocator(10)); ax.yaxis.set_minor_locator(MultipleLocator(0.025))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v,_:M(f'{v:.0f}')))
    ax2.plot(s.dt,100*s.p05,color=G05,lw=LW_MAIN,ls=(0,(4,1.6)),zorder=3)
    ax2.plot(s.dt,100*s.p02,color=G02,lw=LW_MAIN,zorder=4)
    ax2.set_ylim(0,125); ax2.set_yticks(range(0,101,20)); ax2.yaxis.set_minor_locator(mpl.ticker.FixedLocator(range(10,100,20)))
    ax2.set_ylabel('Success rate (%)',labelpad=3,color=G02_T)
    ax2.tick_params(axis='y',colors=G02_T,which='both'); ax2.spines['right'].set_color(G02_T)
    ax.tick_params(axis='y',right=False,which='both'); ax.spines['left'].set_color('k')
    ax2.plot([lev],[100*sa.p02],marker='o',ms=5.0,mfc='white',mec=G02_T,mew=LW_RING,zorder=8,clip_on=False)
    # read-out at the adopted level, in the head-room above 100 %
    ax2.text(-3,121,M(f'adopted level:\n{100*sa.p02:.0f} % \u00b10.02, {100*sa.p05:.0f} % \u00b10.05'),fontsize=FS_ANNOT,color=G02_T,fontweight='bold',
             ha='right',va='top',zorder=9,linespacing=1.15)
    hs=[Line2D([],[],color='k',lw=LW_MAIN,label='mean (estimated \u2212 measured)'),Line2D([],[],color=G02,lw=LW_MAIN,label='% within \u00b10.02'),
        Line2D([],[],color=G05,lw=LW_MAIN,ls=(0,(4,1.6)),label='% within \u00b10.05'),
        Line2D([],[],color='k',lw=LW_ONE,label=M(f'adopted level ({lev:.0f} \u00b0C)')),
        Line2D([],[],color='0.45',lw=LW_THIN,ls=(0,(1.5,1.5)),label=M(f'\u00b10.02 optimum ({opt:.0f} \u00b0C)'))]
    hs=[hs[0],hs[3],hs[4],hs[1],hs[2]]   # column 1: mean and the two levels; column 2: success rates
    ax.legend(handles=hs,loc='lower center',bbox_to_anchor=((W/2-15.5)/57.0,1.01),ncol=2,fontsize=FS_KEY,columnspacing=1.2,handlelength=1.8)
    st=dict(level=float(lev),opt=float(opt),n=int(sa.n),p02=float(sa.p02),p05=float(sa.p05),mean=float(sa['mean']),
            tuned_opt=int(opt),tuned_p02=float(ro.p02),tuned_p05=float(ro.p05),p02_max_range=[float(ties.min()),float(ties.max())])
    log('fig6',st); print('Fig6(sweep)',st)
    return fig,pxsize(W,H)
if __name__=='__main__':
    f,px=draw_fig6(load()); save(f,'Fig6_final.png',px)
