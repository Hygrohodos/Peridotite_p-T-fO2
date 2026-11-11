from unified import *
from matplotlib.ticker import MultipleLocator
FAIL_ALPHA=0.5
def screens(df):
    sp=PARAMS.screening
    four=df["screen_4phases"]==True if "screen_4phases" in df.columns else pd.Series(True,index=df.index)
    TA=pd.to_numeric(df["T_TA98"],errors="coerce"); HA=pd.to_numeric(df["T_HA84"],errors="coerce"); WZ=pd.to_numeric(df["T_WZ07"],errors="coerce"); NG=pd.to_numeric(df["T_NG10"],errors="coerce")
    pT=four&(TA>sp.t_min_c); p1=pT&((TA-NG).abs()<sp.dt_opx_cpx_c)
    HAp=sp.ha84_poly_a+sp.ha84_poly_b*TA+sp.ha84_poly_c*TA**2
    p2=p1&((HA-HAp).abs()<sp.dt_pyx_gt_c); WZp=88.99741613+0.9864930064*HA
    p3=p2&((WZ-WZp).abs()<sp.dt_gt_ol_c)
    return sp,TA,HA,WZ,NG,(pT,p1,p2,p3)
def fig3():
    df=load(); sp,TA,HA,WZ,NG,(pT,p1,p2,p3)=screens(df)
    W,H=174,62; fig=newfig(W,H)
    side=44.0; lefts=[13.0,13+side+13.5,13+2*(side+13.5)]; bot=10.5
    panels=[
     dict(x=TA,y=NG,prev=pT,ok=p1,xl='Two-pyroxene (TA98; \u00b0C)',yl='Ca-in-opx (NG10; \u00b0C)',lim=(700,1500),ttl='Opx\u2013cpx equilibrium',thr=sp.dt_opx_cpx_c,line=lambda t:t),
     dict(x=TA,y=HA,prev=p1,ok=p2,xl='Two-pyroxene (TA98; \u00b0C)',yl='Opx\u2013garnet (HA84; \u00b0C)',lim=(700,1500),ttl='Pyroxene\u2013garnet equilibrium',thr=sp.dt_pyx_gt_c,line=lambda t:sp.ha84_poly_a+sp.ha84_poly_b*t+sp.ha84_poly_c*t**2),
     dict(x=HA,y=WZ,prev=p2,ok=p3,xl='Opx\u2013garnet (HA84; \u00b0C)',yl='Olivine\u2013garnet (WZ07; \u00b0C)',lim=(880,1480),ttl='Garnet\u2013olivine equilibrium',thr=sp.dt_gt_ol_c,line=lambda t:88.99741613+0.9864930064*t)]
    counts=[]
    for i,(p,L) in enumerate(zip(panels,'abc')):
        ax=axmm(fig,W,H,lefts[i],bot,side,side)
        lo,hi=p['lim']; t=np.linspace(lo,hi,200); yy=p['line'](t); ok=p['ok']; fail=p['prev']&~p['ok']
        ax.plot(t,yy,color='k',lw=LW_ONE,zorder=4)
        for sgn in (1,-1): ax.plot(t,yy+sgn*p['thr'],color='k',lw=LW_THIN,ls=(0,(4,2)),zorder=4)
        ax.scatter(p['x'][ok],p['y'][ok],s=S_DENSE,c=PASS,edgecolors=PASS_E,linewidths=EW_DENSE,zorder=3)
        ax.scatter(p['x'][fail],p['y'][fail],marker='x',s=17,c=FAIL,linewidths=LW_THIN+0.1,alpha=FAIL_ALPHA,zorder=1.5)
        ax.set_xlim(lo,hi); ax.set_ylim(lo,hi); ax.set_box_aspect(1)
        tk=list(range(800,1501,200)) if lo==700 else list(range(900,1401,200))
        ax.set_xticks(tk); ax.set_yticks(tk); ax.xaxis.set_minor_locator(MultipleLocator(100)); ax.yaxis.set_minor_locator(MultipleLocator(100))
        ax.set_xlabel(p['xl']); ax.set_ylabel(p['yl'],labelpad=2)
        ax.set_title(p['ttl'],fontsize=FS_ANNOT,pad=3.5,loc='center')
        letter(ax,L+'.',x=0.035,y=0.965)
        n_ok,n_fail=int(ok.sum()),int(fail.sum()); counts.append((n_ok,n_fail))
        hs=[Line2D([],[],marker='o',ls='',mfc=PASS,mec=PASS_E,mew=EW_STD,ms=3.8,label=f'pass {n_ok}'),
            Line2D([],[],marker='x',ls='',mec=mpl.colors.to_rgba(FAIL,FAIL_ALPHA),mew=LW_THIN+0.1,ms=4.2,label=f'fail {n_fail}'),
            Line2D([],[],color='k',lw=LW_THIN,ls=(0,(4,2)),label=f'\u00b1{p["thr"]:.0f} \u00b0C')]
        ax.legend(handles=hs,loc='lower right',frameon=False,fontsize=FS_KEY,handlelength=1.6,handletextpad=0.5,labelspacing=0.25,borderaxespad=0.4)
    print('counts',counts); assert counts==[(583,24),(495,88),(450,45)],counts
    return fig,(4110,1464)
if __name__=='__main__':
    f,px=fig3(); save(f,'Fig3_final.png',px)
