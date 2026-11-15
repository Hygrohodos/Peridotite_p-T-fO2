"""Fig. 8: estimated vs measured garnet Fe3+/SumFe for the screened Siberian and Kaapvaal validation
samples. The estimates are the pipeline column XFe3_gt_corr, i.e. at the adopted match level (make_paper_figures.py)."""
from fe3rows import *
from common import log
def draw_fig8(df):
    """Fig. 8 from the pipeline table df."""
    meas=pd.to_numeric(df["garnet_Fe3+"],errors="coerce")
    est=pd.to_numeric(df["XFe3_gt_corr"],errors="coerce")
    T=pd.to_numeric(df["T_TA98"],errors="coerce"); cr=df["Craton"].fillna("").astype(str); slave=cr.str.contains("Slave")
    passed=df["screen_pass"]==True
    ok=(~slave)&passed&meas.notna()&(meas>0)&(meas<0.5)&est.notna(); resid=meas-est
    x_max=float(max(meas[ok].max(),est[ok].max())+0.02); t_min=float(T[ok].min()-40); t_max=float(T[ok].max()+40)
    r_min=float(resid[ok].min()-0.02); r_max=float(resid[ok].max()+0.02)
    W,H=174,72; fig=newfig(W,H); b0=11.0; ah=H-b0-8.2
    aL=axmm(fig,W,H,17.4,b0,64.5,ah); aR=axmm(fig,W,H,102.4,b0,64.5,ah)
    aL.fill_between([0,x_max],[-0.05,x_max-0.05],[0.05,x_max+0.05],color=GREEN_OUTER,lw=0,zorder=0)
    aL.fill_between([0,x_max],[-0.02,x_max-0.02],[0.02,x_max+0.02],color=GREEN_INNER,lw=0,zorder=0)
    aL.plot([0,x_max],[0,x_max],color='0.35',lw=LW_ONE,zorder=1); band_fill(aR,'resid'); aR.axhline(0,color='0.35',lw=LW_THIN,zorder=1)
    # non-Slave craton symbols; 'Other' open circles beneath (same rule as the repo)
    for k,mk,c,lab,_ in CRATON:
        if 'Slave' in k: continue
        sel=ok&(cr==k)
        if sel.any():
            aL.scatter(meas[sel],est[sel],s=S_STD,c=c,marker=mk,edgecolors='k',linewidths=EW_STD,zorder=3)
            aR.scatter(T[sel],resid[sel],s=S_STD,c=c,marker=mk,edgecolors='k',linewidths=EW_STD,zorder=3)
    other=ok&~cr.isin([c[0] for c in CRATON])
    if other.any():
        aL.scatter(meas[other],est[other],s=S_OTHER,facecolors='none',edgecolors=OTHER_EC,linewidths=EW_OTHER,zorder=2); aR.scatter(T[other],resid[other],s=S_OTHER,facecolors='none',edgecolors=OTHER_EC,linewidths=EW_OTHER,zorder=2)
    aL.set_xlim(0,x_max); aL.set_ylim(0,x_max); aR.set_xlim(900,1400)   # as Fig. 7b; data span t_min+40 to t_max-40; aR.set_ylim(r_min,r_max)
    aL.set_xlabel(XM); aL.set_ylabel(YE); aR.set_xlabel(XT_TA98); aR.set_ylabel(YR)
    fmt(aL,np.arange(0,x_max,0.05),np.arange(0,x_max,0.05),'%.2f','%.2f'); fmt(aR,np.arange(900,1401,100),np.arange(-0.15,0.1001,0.05),'%d','%.2f')
    letter(aL,'a'); letter(aR,'b')
    diff=(est[ok]-meas[ok]).abs()
    aL.text(0.98,0.04,f"n = {int(ok.sum())}\n\u00b10.02: {(diff<=0.02).mean()*100:.0f}%\n\u00b10.05: {(diff<=0.05).mean()*100:.0f}%",transform=aL.transAxes,fontsize=FS_ANNOT,ha='right',va='bottom',zorder=20,path_effects=halo(2.0))
    log('fig8_validation',dict(n=int(ok.sum()),p02=float((diff<=0.02).mean()),p05=float((diff<=0.05).mean())))
    print('Fig8 n',int(ok.sum()),(diff<=0.02).mean(),(diff<=0.05).mean())
    h=[kh('sib','Siberian'),kh('kaa','Kaapvaal')]+([kh('oth','Other')] if other.any() else [])+[Patch(fc=GREEN_INNER,ec='none',label=r'within $\pm$0.02'),Patch(fc=GREEN_OUTER,ec='none',label=r'within $\pm$0.05')]
    fig.legend(handles=h,loc='upper center',bbox_to_anchor=(0.5,1-0.8/H),ncol=len(h),fontsize=FS_KEY,columnspacing=2.0,borderaxespad=0)
    return fig,pxsize(W,H)
if __name__=='__main__':
    f,px=draw_fig8(load()); save(f,'Fig8_final.png',px)
