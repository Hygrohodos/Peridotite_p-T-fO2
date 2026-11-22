"""Fig. 5 (revision v7_cr): TA98 - WZ07 against T_TA98, calculated with the measured garnet Fe3+, for the natural
calibration + validation xenoliths split at garnet Cr2O3 4.9 wt% (gap 4.49-5.32 wt%; a: low Cr, b: high Cr), each with
its own OLS fit and 95 % band (the other group's fit dotted), colour = garnet Cr2O3, and the 13 experiments with
their weighted fit and 95 % band (circle area = weight; scripts/make_fig05_weightcircles.py). Natural samples with
T_TA98 >= 1110 deg C, excluding the error groups and Uv-50-04 (n = 15 + 23). House style: unified.py.
The Cr2O3 dependence shown here is the basis of the correction in calculations/pt_fo2/cr_correction.py."""
import json
from common import *
import unified as uf
import matplotlib as mpl
import statsmodels.formula.api as smf
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator, FuncFormatter
from scipy import stats
import weighted_thermometer_regression as wtr
from make_fig05_weightcircles import EXP_STYLE, SIZE_PER_PCT
from calculations.pt_fo2.level_anchor import discrepancy_measured
PIV=1207.0242233332644; SPLIT=4.905            # garnet Cr2O3 gap 4.49-5.32 wt%
WS=SIZE_PER_PCT*0.56                            # circle scale as in figures_v7 Fig. 5 (area proportional to weight)


def natural_set(table):
    """Screened Slave calibration + Siberia/Kaapvaal validation samples with measured garnet Fe3+, with
    D = T_TA98 - T_WZ07(measured Fe3+), T, garnet Cr2O3 and craton."""
    m = (calib.calibration_mask(table) | calib.validation_mask(table)) & calib.measured_fe3(table).notna()
    t = table[m]
    cr = t['Craton'].fillna('').astype(str)
    return pd.DataFrame({'id': t.SampleID.astype(str).str.strip().values,
                         'T': pd.to_numeric(t.T_TA98, errors='coerce').values,
                         'D': discrepancy_measured(table)[m].values,
                         'garnet_Cr2O3': pd.to_numeric(t.garnet_Cr2O3, errors='coerce').values,
                         'set': np.where(cr.str.contains('Slave'), 'Slave', np.where(cr.str.contains('Siber'), 'Siberia', 'Kaapvaal'))})


def draw_fig5_cr(table):
    """Return (fig, px, stats) for Fig. 5."""
    CUT=1110.0   # t1876u: natural fit cut lowered from 1120 to 1110 C, which adds only 22-7 (1113.6 C, Cr2O3 8.05 wt%) to the high-Cr fit
    GA=natural_set(table); GA['crat']=GA['set']
    G=GA[(GA['T']>=CUT)&~GA.id.isin(['14-107','MX088','MX044','MX104','Uv-50-04'])].copy(); G['y']=G.D; G['t']=G['T']-PIV   # n = 38
    G['grp']=np.where(G.garnet_Cr2O3<SPLIT,'low','high')
    fits={g:smf.ols('y~t',G[G.grp==g]).fit() for g in ('low','high')}
    MI=smf.ols('y~t*C(grp)',G).fit(); KI=[c for c in MI.params.index if c.startswith('t:')][0]
    # SS scan over every gap between adjacent Cr2O3 values (>= 8 per group): total residual SS of the two-line fit
    u=np.sort(G.garnet_Cr2O3.unique()); SCAN=[]
    for th in (u[:-1]+u[1:])/2:
        lo_=G.garnet_Cr2O3<th
        if lo_.sum()>=8 and (~lo_).sum()>=8: SCAN.append((float(th),int(lo_.sum()),float(sum(smf.ols('y~t',G[m]).fit().ssr for m in (lo_,~lo_)))))
    SS1=float(smf.ols('y~t',G).fit().ssr)
    runs=wtr.load_runs(); wt,fp=wtr.fit(runs,t_ref=PIV); runs=runs.merge(wt[['SampleID','weight_pct']],on='SampleID')
    YLIM=(-180,200); XLIM=(1100,1650)
    clip=runs.y_TA98_minus_WZ07_C<YLIM[0]
    # ---------------- figure ----------------
    W,H=uf.COL2,84.0; fig=uf.newfig(W,H)
    AXW,AXH,B0=62.0,58.0,11.0; A=uf.axmm(fig,W,H,15.0,B0,AXW,AXH); B=uf.axmm(fig,W,H,15.0+AXW+8.0,B0,AXW,AXH); CB=uf.axmm(fig,W,H,15.0+2*AXW+15.0,B0,2.6,AXH)
    cmap=mpl.colormaps['YlOrRd']; cmap=mpl.colors.LinearSegmentedColormap.from_list('ylord',cmap(np.linspace(0.08,0.92,256)))
    norm=mpl.colors.Normalize(vmin=1.5,vmax=10.0)
    MK={'Slave':'o','Siberia':'s','Kaapvaal':'^'}; SZ={'Slave':uf.S_STD,'Siberia':uf.S_STD*0.9,'Kaapvaal':uf.S_STD*1.15}
    tt_e=np.linspace(runs.T_TA98.min(),runs.T_TA98.max(),200); elo,ehi=fp.band(tt_e)
    def band(m,x):
        Xd=np.column_stack([np.ones_like(x),x-PIV]); yh=Xd@m.params.values; se=np.sqrt(np.einsum('ij,jk,ik->i',Xd,m.cov_params().values,Xd)); q=stats.t.ppf(0.975,m.df_resid); return yh,yh-q*se,yh+q*se
    for P,g,s in ((A,'low','a'),(B,'high','b')):
        o='high' if g=='low' else 'low'; d=G[G.grp==g]; do=G[G.grp==o]; m=fits[g]
        P.axhline(0,color='0.6',lw=uf.LW_THIN,zorder=0)
        P.fill_between(tt_e,elo,ehi,color='0.92',lw=0,zorder=0.3); P.plot(tt_e,fp.b0+fp.b1*(tt_e-PIV),color='0.25',lw=uf.LW_ONE,ls=(0,(5,2.5)),zorder=1.5)
        xs=np.linspace(d['T'].min(),d['T'].max(),100); yh,lo,hi=band(m,xs); P.fill_between(xs,lo,hi,color='0.80',lw=0,zorder=0.6); P.plot(xs,yh,color='k',lw=uf.LW_MAIN,zorder=4)
        xo=np.linspace(do['T'].min(),do['T'].max(),2); P.plot(xo,fits[o].params.Intercept+fits[o].params.t*(xo-PIV),color='0.45',lw=uf.LW_ONE,ls=(0,(1,1.6)),zorder=3.5)
        for ref_,sty in EXP_STYLE.items():
            e=runs[(runs.Reference==ref_)&~clip]; col=sty['color']
            P.scatter(e.T_TA98,e.y_TA98_minus_WZ07_C,s=WS*e.weight_pct,marker='o',facecolors=col,alpha=0.18,edgecolors='none',zorder=2)
            P.scatter(e.T_TA98,e.y_TA98_minus_WZ07_C,s=WS*e.weight_pct,marker='o',facecolors='none',edgecolors=col,linewidths=0.6,alpha=0.85,zorder=2)
            P.scatter(e.T_TA98,e.y_TA98_minus_WZ07_C,s=uf.S_STD*0.75,marker=sty['marker'],color=col,edgecolors='k',linewidths=uf.EW_STD,zorder=5)
            ec=runs[(runs.Reference==ref_)&clip]   # off-scale runs: arrow at the lower edge
            for _,r in ec.iterrows():
                P.annotate('',xy=(r.T_TA98,YLIM[0]+3),xytext=(r.T_TA98,YLIM[0]+22),arrowprops=dict(arrowstyle='-|>,head_length=0.35,head_width=0.18',color=col,lw=0.9,shrinkA=0,shrinkB=0),zorder=6)
                right=r.T_TA98>XLIM[1]-60
                P.text(r.T_TA98+(-7 if right else 7),YLIM[0]+13,uf.M(f'{r.y_TA98_minus_WZ07_C:.0f}'),fontsize=uf.FS_SMALL,color=col,va='center',ha='right' if right else 'left',zorder=6)
        for c,mk in MK.items():
            e=d[d.crat==c]; P.scatter(e['T'],e.y,s=SZ[c],marker=mk,c=e.garnet_Cr2O3,cmap=cmap,norm=norm,edgecolors='k',linewidths=uf.EW_STD,zorder=6)
        P.set_xlim(*XLIM); P.set_ylim(*YLIM)
        P.xaxis.set_major_locator(MultipleLocator(100)); P.xaxis.set_minor_locator(MultipleLocator(25))
        P.annotate('weighted fit',xy=(1450,fp.b0+fp.b1*(1450-PIV)),xytext=(6,7),textcoords='offset points',fontsize=uf.FS_ANNOT,color='0.25',ha='left',va='bottom',zorder=7)
        P.yaxis.set_major_locator(MultipleLocator(100)); P.yaxis.set_minor_locator(MultipleLocator(25))
        P.yaxis.set_major_formatter(FuncFormatter(lambda v,_: uf.M(f'{v:.0f}')))
        P.set_xlabel(uf.XT_TA98); uf.letter(P,s)
        head=('Cr$_2$O$_3$ < ' if g=='low' else 'Cr$_2$O$_3$ \u2265 ')+f'{SPLIT:.1f} wt% (n = {int(m.nobs)})'
        P.text(0.035,0.965,head+'\n'+uf.M(f'slope {m.params.t:.2f} \u00b1 {m.bse.t:.2f}'),transform=P.transAxes,va='top',ha='left',fontsize=uf.FS_ANNOT,linespacing=1.3,zorder=8)
    A.set_ylabel(uf.DT_LABEL+', Fe$^{2{+}}$ only'); B.tick_params(labelleft=False)
    cb=fig.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=cmap),cax=CB); cb.set_label('Garnet Cr$_2$O$_3$ (wt%)'); cb.outline.set_linewidth(uf.LW_FRAME)
    CB.tick_params(direction='in',width=uf.LW_TICK,which='both',right=True,left=False); CB.minorticks_off(); CB.yaxis.set_major_locator(MultipleLocator(2))
    CB.axhline(SPLIT,color='k',lw=1.0); CB.text(-0.6,SPLIT,f'{SPLIT:.1f}',transform=CB.get_yaxis_transform(),ha='right',va='center',fontsize=uf.FS_SMALL)
    # ---------------- key: two rows above the panels ----------------
    nk=lambda mk,ms: Line2D([],[],marker=mk,ls='',mfc='0.75',mec='k',mew=uf.EW_STD,ms=ms)
    ek=lambda sty: Line2D([],[],marker=sty['marker'],ls='',mfc=sty['color'],mec='k',mew=uf.EW_STD,ms=uf.MS_KEY*0.9)
    ok=lambda mk,ms: Line2D([],[],marker=mk,ls='',mfc='w',mec='0.35',mew=1.1,ms=ms)
    r1=[nk('o',uf.MS_KEY),nk('s',uf.MS_KEY),nk('^',5.3),Line2D([],[],color='0.45',lw=uf.LW_ONE,ls=(0,(1,1.6)))]
    l1=['Slave','Siberia','Kaapvaal','other Cr$_2$O$_3$ group']
    mat,rz=EXP_STYLE['Matjuschkin et al. 2014'],EXP_STYLE['Rzehak et al. 2020']
    r2=[ek(mat),ek(rz)]
    l2=[f"{mat['label']} (n = {int((runs.Reference=='Matjuschkin et al. 2014').sum())})",f"{rz['label']} (n = {int((runs.Reference=='Rzehak et al. 2020').sum())})"]
    L1=fig.legend(handles=r1,labels=l1,loc='upper left',bbox_to_anchor=(15.0/W,1-1.0/H),ncol=4,frameon=False,handlelength=1.8,columnspacing=1.4)
    L2=fig.legend(handles=r2,labels=l2,loc='upper left',bbox_to_anchor=(15.0/W,1-6.0/H),ncol=2,frameon=False,handlelength=1.8,columnspacing=1.4)
    # weight key (circle area = % of total weight), right of the second key row
    x0=120.0; yk=1-3.9/H
    fig.text(x0/W,yk,'weight',ha='right',va='center',fontsize=uf.FS_KEY)
    for k,dx in zip((5,10,20),(3.0,13.5,25.5)):
        for fc,ec,al,mw in (('0.5','none',0.25,0),('none','0.3',1.0,0.6)):
            fig.add_artist(Line2D([(x0+dx)/W],[yk],marker='o',ls='',ms=np.sqrt(WS*k),mfc=fc,mec=ec,alpha=al,mew=mw,transform=fig.transFigure))
        fig.text((x0+dx+np.sqrt(WS*k)/2/72*25.4+0.8)/W,yk,f'{k}%',va='center',ha='left',fontsize=uf.FS_KEY)
    px=uf.pxsize(W,H)
    st=dict(split='garnet Cr2O3',split_wt=SPLIT,split_gap=[4.49,5.32],groups={g:dict(n=int(m.nobs),slope=m.params.t,se=m.bse.t,p=m.pvalues.t,level=m.params.Intercept) for g,m in fits.items()},
        slope_difference=dict(coef=MI.params[KI],se=MI.bse[KI],p=MI.pvalues[KI]),natural_T_cut=CUT,
        experiments=dict(n=len(runs),slope=fp.b1,se=fp.se_b1,level_at_pivot=fp.b0,se_level=fp.se_b0,n_eff=fp.n_eff,off_scale=list(runs.SampleID[clip])))
    return fig, px, st


if __name__ == '__main__':
    fig, px, st = draw_fig5_cr(load())
    uf.save(fig, 'Fig5_final.png', px); print(json.dumps(st, indent=1, default=float))
