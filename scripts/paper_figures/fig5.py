from fe3rows import *
def fig5():
    df=load(); base,meas=mf._slave_fe3_mask(df); T=pd.to_numeric(df["T_TA98"],errors="coerce")
    dts=(50,25,10); est={dt:mf._fe3_at_dt(df,dt) for dt in dts}
    fit=calib.fit_empirical_correction(df); assert int(PARAMS.fe3.emp_fit_dt_c) in dts
    blue,red,rem=fit["masks"]["blue"],fit["masks"]["red"],fit["masks"]["reg"]; sl,ic=fit["slope"],fit["intercept"]
    print('Fig5 Eq5 slope',sl,'icpt',ic,'n_reg',fit['n_regression'],'blue',fit['n_blue'],'red',fit['n_red'])
    W,H=174,210; fig=newfig(W,H)
    top=8.5; bot=12.0; gap=8.4; h=(H-top-bot-2*gap)/3; wA=64.5; lL=17.4; lR=102.4
    for i,dt in enumerate(dts):
        y0=H-top-(i+1)*h-i*gap
        aL=axmm(fig,W,H,lL,y0,wA,h); aR=axmm(fig,W,H,lR,y0,wA,h)
        e=est[dt]; ok=base&e.notna()&T.notna(); fitrow=(dt==PARAMS.fe3.emp_fit_dt_c)
        row(aL,aR,df,meas,T,e,ok,blue=blue if fitrow else None,red=red if fitrow else None,fit=(sl,ic) if fitrow else None,hull=rem if fitrow else None)
        letter(aL,'abc'[i]+'.',x=0.03,y=0.97); aL.text(0.12,0.965,f'\u0394T = {dt:.0f} \u00b0C',transform=aL.transAxes,fontsize=FS_ANNOT,va='top',ha='left',zorder=20)
        if i<2: aL.set_xlabel(''); aR.set_xlabel('')
        print('  dT',dt,'n',int(ok.sum()))
    fig.legend(handles=key_handles(),loc='upper center',bbox_to_anchor=(0.5,1-0.8/H),ncol=5,frameon=False,fontsize=FS_KEY,handletextpad=0.4,columnspacing=2.2,handlelength=1.6,borderaxespad=0)
    return fig,(4110,4960)
if __name__=='__main__':
    f,px=fig5(); save(f,'Fig5_final.png',px)
