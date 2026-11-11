"""Fig. 5 (revision v2): estimated vs measured garnet Fe3+/SumFe at three match levels
(adopted level +-25 deg C), weighted experimental correction, no Eq. 5, no regression envelope."""
from fe3rows import *
from v2common import *
def fig5():
    df=load(); base,meas,blue,red,rem=mf._error_groups(df); T=pd.to_numeric(df["T_TA98"],errors="coerce")
    opt=FP.dt_match_c; dts=(opt+25,opt,opt-25)
    W,H=174,210; fig=newfig(W,H)
    top=8.5; bot=12.0; gap=8.4; h=(H-top-bot-2*gap)/3; wA=64.5; lL=17.4; lR=102.4
    out={}
    for i,dt in enumerate(dts):
        y0=H-top-(i+1)*h-i*gap
        aL=axmm(fig,W,H,lL,y0,wA,h); aR=axmm(fig,W,H,lR,y0,wA,h)
        e=fe3_at(df,dt); ok=base&e.notna()&T.notna(); mid=(dt==opt)
        row(aL,aR,df,meas,T,e,ok,blue=blue if mid else None,red=red if mid else None)
        aL.set_ylim(0,0.225); fmt(aL,np.arange(0,0.351,0.05),np.arange(0,0.2251,0.025),'%.2f','%.3f')
        letter(aL,'abc'[i]+'.',x=0.03,y=0.97); aL.text(0.12,0.965,M(f'\u0394T = {dt+MEANC:.0f} \u00b0C'),transform=aL.transAxes,fontsize=FS_ANNOT,va='top',ha='left',zorder=20)
        d=(e[ok]-meas[ok]); r=ok&~blue&~red; dr=(e[r]-meas[r])
        aL.text(0.98,0.04,f"all n = {int(ok.sum())}\n\u00b10.02: {pct(d,0.02)*100:.0f}%\n\u00b10.05: {pct(d,0.05)*100:.0f}%",transform=aL.transAxes,fontsize=FS_ANNOT,ha='right',va='bottom',zorder=20,path_effects=halo(2.0))
        res=(meas-e)[ok]; Tl=T[ok]
        clip=int(((res<-0.12)|(res>0.20)|(Tl<900)|(Tl>1400)|(e[ok]>0.20)).sum())
        out[round(float(dt),1)]=dict(n=int(ok.sum()),p02=pct(d,0.02),p05=pct(d,0.05),n_rem=int(r.sum()),p02_rem=pct(dr,0.02),p05_rem=pct(dr,0.05),mean=float(d.mean()),points_outside_axes=clip)
        if i<2: aL.set_xlabel(''); aR.set_xlabel('')
    fig.legend(handles=key_handles(),loc='upper center',bbox_to_anchor=(0.5,1-0.8/H),ncol=5,frameon=False,fontsize=FS_KEY,handletextpad=0.4,columnspacing=2.2,handlelength=1.6,borderaxespad=0)
    log('fig06_rows',out); print('Fig5',out)
    return fig,(4110,4960)
if __name__=='__main__':
    f,px=fig5(); save(f,'Fig6_final.png',px)
