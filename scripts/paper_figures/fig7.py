from fe3rows import *
def fig7():
    df=load(); base,meas,blue,red,rem=mf._error_groups(df); T=pd.to_numeric(df["T_TA98"],errors="coerce")
    opt=PARAMS.fe3.dt_match_c; raw=mf._fe3_at_dt(df,opt); sl,ic=PARAMS.fe3.emp_slope,PARAMS.fe3.emp_intercept
    est=mf._apply_emp(raw,T,slope=sl,intercept=ic); ok=base&est.notna()&T.notna()
    W,H=174,72; fig=newfig(W,H); b0=12.4; ah=H-b0-8.2
    aL=axmm(fig,W,H,17.4,b0,64.5,ah); aR=axmm(fig,W,H,102.4,b0,64.5,ah)
    row(aL,aR,df,meas,T,est,ok,blue=blue,red=red)
    letter(aL,'a.'); letter(aR,'b.')
    da=(est[ok]-meas[ok]).abs(); dr=(est[rem&ok]-meas[rem&ok]).abs()
    aL.text(0.98,0.04,f"all n = {int(ok.sum())}\n\u00b10.02: {(da<=0.02).mean()*100:.0f}%\n\u00b10.05: {(da<=0.05).mean()*100:.0f}%\nremainder n = {int((rem&ok).sum())}\n\u00b10.02: {(dr<=0.02).mean()*100:.0f}%\n\u00b10.05: {(dr<=0.05).mean()*100:.0f}%",
            transform=aL.transAxes,fontsize=FS_ANNOT,ha='right',va='bottom',zorder=20,path_effects=halo(2.0))
    print('Fig7 dT',opt,'all n',int(ok.sum()),(da<=0.02).mean(),(da<=0.05).mean(),'rem n',int((rem&ok).sum()),(dr<=0.02).mean(),(dr<=0.05).mean())
    fig.legend(handles=key_handles(),loc='upper center',bbox_to_anchor=(0.5,1-0.8/H),ncol=5,frameon=False,fontsize=FS_KEY,handletextpad=0.4,columnspacing=2.2,handlelength=1.6,borderaxespad=0)
    return fig,(4110,1700)
if __name__=='__main__':
    f,px=fig7(); save(f,'Fig7_final.png',px)
