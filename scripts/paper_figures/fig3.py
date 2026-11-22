from unified import *
from common import log
from matplotlib.ticker import MultipleLocator
FAIL_ALPHA=0.5
def screens(df):
    """Thermometers and the successive screen masks of the pipeline table (calculations/pt_fo2/screening.py)."""
    sp=PARAMS.screening; num=lambda c: pd.to_numeric(df[c],errors="coerce")
    TA,HA,WZ,NG=num("T_TA98"),num("T_HA84"),num("T_WZ07"),num("T_NG10")
    masks=tuple(df[c]==True for c in ("screen_T900","screen_opx_cpx","screen_pyx_gt","screen_gt_ol"))
    return sp,TA,HA,WZ,NG,masks
def fig3(df):
    """Fig. 3 from the pipeline table df: samples entering each equilibrium screen (quality and T > 900 C passed)."""
    sp,TA,HA,WZ,NG,(pT,p1,p2,p3)=screens(df)
    W,H=174,64; fig=newfig(W,H)
    side=44.0; lefts=[13.5,13.5+side+13.5,13.5+2*(side+13.5)]; bot=10.5
    panels=[
     dict(x=TA,y=NG,prev=pT,ok=p1,xl=XT_TA98,yl='Ca-in-opx $T_\\mathrm{NG10}$ (\u00b0C)',lim=(700,1500),ttl='Opx\u2013cpx equilibrium',thr=sp.dt_opx_cpx_c,line=lambda t:t),
     dict(x=TA,y=HA,prev=p1,ok=p2,xl=XT_TA98,yl='Opx\u2013garnet $T_\\mathrm{HA84}$ (\u00b0C)',lim=(700,1500),ttl='Pyroxene\u2013garnet equilibrium',thr=sp.dt_pyx_gt_c,line=lambda t:sp.ha84_poly_a+sp.ha84_poly_b*t+sp.ha84_poly_c*t**2),
     dict(x=HA,y=WZ,prev=p2,ok=p3,xl='Opx\u2013garnet $T_\\mathrm{HA84}$ (\u00b0C)',yl='Olivine\u2013garnet $T_\\mathrm{WZ07}$ (\u00b0C)',lim=(880,1480),ttl='Garnet\u2013olivine equilibrium',thr=sp.dt_gt_ol_c,line=lambda t:88.99741613+0.9864930064*t)]
    counts=[]
    for i,(p,L) in enumerate(zip(panels,'abc')):
        ax=axmm(fig,W,H,lefts[i],bot,side,side)
        lo,hi=p['lim']; t=np.linspace(lo,hi,200); yy=p['line'](t); ok=p['ok']; fail=p['prev']&~p['ok']
        ax.plot(t,yy,color='k',lw=LW_ONE,zorder=4)
        for sgn in (1,-1): ax.plot(t,yy+sgn*p['thr'],color='k',lw=LW_THIN,ls=(0,(4,2)),zorder=4)
        ax.scatter(p['x'][ok],p['y'][ok],s=S_DENSE,c=PASS,edgecolors=PASS_E,linewidths=EW_DENSE,zorder=3)
        ax.scatter(p['x'][fail],p['y'][fail],marker='x',s=15,c=FAIL,linewidths=0.8,alpha=FAIL_ALPHA,zorder=1.5)
        ax.set_xlim(lo,hi); ax.set_ylim(lo,hi); ax.set_box_aspect(1)
        tk=list(range(800,1501,200)) if lo==700 else list(range(900,1401,200))
        ax.set_xticks(tk); ax.set_yticks(tk); ax.xaxis.set_minor_locator(MultipleLocator(100)); ax.yaxis.set_minor_locator(MultipleLocator(100))
        ax.set_xlabel(p['xl']); ax.set_ylabel(p['yl'],labelpad=2)
        ax.set_title(p['ttl'],fontsize=FS_LABEL,pad=4,loc='center')
        letter(ax,L)
        n_ok,n_fail=int(ok.sum()),int(fail.sum()); counts.append((n_ok,n_fail))
        hs=[Line2D([],[],marker='o',ls='',mfc=PASS,mec=PASS_E,mew=EW_STD,ms=3.8,label=f'pass {n_ok}'),
            Line2D([],[],marker='x',ls='',mec=mpl.colors.to_rgba(FAIL,FAIL_ALPHA),mew=0.8,ms=4.0,label=f'fail {n_fail}'),
            Line2D([],[],color='k',lw=LW_THIN,ls=(0,(4,2)),label=f'\u00b1{p["thr"]:.0f} \u00b0C')]
        ax.legend(handles=hs,loc='lower right',fontsize=FS_KEY,labelspacing=0.3,borderaxespad=0.3)
    print('counts',counts); log('fig3',dict(entering=int(pT.sum()),pass_fail=counts))
    return fig,pxsize(W,H)
if __name__=='__main__':
    f,px=fig3(load()); save(f,'Fig3_final.png',px)
