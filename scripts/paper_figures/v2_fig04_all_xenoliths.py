"""Fig. X (new; Paolo Nimis, line 255): TA98 vs WZ07 for all 57 xenoliths with measured garnet Fe3+.
(a) WZ07 vs TA98 with Fe3+ = 0 (open) and measured Fe3+ (filled); (b) TA98 - WZ07 (measured Fe3+) vs TA98,
with the 13 experiments, the weighted regression (stand-alone module) and its 95 % band, and the
xenolith mean +- 1 sigma (error groups excluded). Per-sample WZ07 values are recomputed here from the
pipeline output, not read from a side file."""
from v2common import *
from calculations.pt_fo2.fe3_correction import _wz07corr_fe2only
EXP_M='#5e3c99'; EXP_R='#4c72b0'; XMEAN='#b35806'; XBAND='#fde0c5'
def per_sample(df):
    base,meas,blue,red,rem=mf._error_groups(df)
    val=calib.validation_mask(df) & meas.notna()
    use=(base|val)&meas.notna()
    u=df[use].copy(); m=meas[use]
    num=lambda c: pd.to_numeric(u[c],errors='coerce').to_numpy(float)
    P=num('P_CA91'); gFe=num('garnet_Fe')
    def wz(x): return _wz07corr_fe2only(num('olivine_Mg'),num('olivine_Fe'),num('garnet_Mg'),gFe*(1.0-x),num('garnet_Ca'),num('garnet_Mn'),P,offset_c=0.0)
    u['T_TA98']=num('T_TA98'); u['WZ0']=wz(np.zeros(len(u))); u['WZm']=wz(m.to_numpy(float))
    u['dTm']=u.T_TA98-u.WZm; u['grp']=np.where(blue[use],'blue',np.where(red[use],'red',''))
    return u
def figX():
    df=load(); u=per_sample(df); cr=u['Craton'].fillna('').astype(str)
    W,H=174,84; fig=newfig(W,H)
    a=axmm(fig,W,H,14.5,11.5,62,62); b=axmm(fig,W,H,98.5,11.5,72,62)
    lim=(900,1450); a.plot(lim,lim,color='0.35',lw=LW_ONE,zorder=1)
    for k,mk,c,lab,_ in CRATON:
        s=(cr==k).to_numpy()
        a.scatter(u.T_TA98[s],u.WZ0[s],s=S_STD*0.8,facecolors='none',edgecolors=c,marker=mk,linewidths=0.8,zorder=2)
        a.scatter(u.T_TA98[s],u.WZm[s],s=S_STD,c=c,marker=mk,edgecolors='k',linewidths=EW_STD,zorder=3)
        b.scatter(u.T_TA98[s],u.dTm[s],s=S_STD,c=c,marker=mk,edgecolors='k',linewidths=EW_STD,zorder=4)
    for g,col in (('blue',BLUE_RING),('red',RED_RING)):
        s=(u.grp==g).to_numpy(); b.scatter(u.T_TA98[s],u.dTm[s],s=S_RING,facecolors='none',edgecolors=col,linewidths=LW_RING,zorder=5)
    a.set_xlim(lim); a.set_ylim(lim); a.set_xticks(range(900,1451,100)); a.set_yticks(range(900,1451,100)); minor(a)
    a.set_xlabel('TA98 temperature (\u00b0C)'); a.set_ylabel('WZ07 temperature (\u00b0C)'); letter(a,'a.')
    a.text(0.05,0.80,'1:1',transform=a.transAxes,fontsize=FS_ANNOT,color='0.35')
    ka=[kh('slave','Slave'),kh('sib','Siberian'),kh('kaa','Kaapvaal'),
        Line2D([],[],marker='o',ls='',mfc='k',mec='k',ms=MS_KEY,label='measured Fe$^{3+}$'),
        Line2D([],[],marker='o',ls='',mfc='none',mec='k',mew=0.8,ms=MS_KEY,label='Fe$^{3+}$ = 0')]
    a.legend(handles=ka,loc='lower right',frameon=False,fontsize=FS_KEY,handletextpad=0.3,labelspacing=0.3,borderaxespad=0.3)
    # panel b
    b.axhline(0,color='0.55',lw=LW_THIN,zorder=1)
    tt=np.linspace(TLO,THI,200); lo,hi=FIT.band(tt)
    b.fill_between(tt,lo,hi,color='0.85',lw=0,zorder=0.5); b.plot(tt,FIT.predict(tt),color='k',lw=LW_MAIN,zorder=3)
    for refn,mk,c in (('Matjuschkin','o',EXP_M),('Rzehak','s',EXP_R)):
        s=RUNS.Reference.str.startswith(refn)
        b.scatter(RUNS.T_TA98[s],RUNS.y_TA98_minus_WZ07_C[s],marker=mk,s=S_STD,facecolors='white',edgecolors=c,linewidths=1.1,zorder=4)
    ex=(u.grp=='').to_numpy(); mu=float(u.dTm[ex].mean()); sd=float(u.dTm[ex].std(ddof=1)); n_ex=int(ex.sum())
    b.axhspan(mu-sd,mu+sd,color=XBAND,lw=0,zorder=0.3); b.axhline(mu,color=XMEAN,lw=LW_REF,ls=(0,(5,2)),zorder=2)
    b.set_xlim(900,1650); b.set_ylim(-330,230); b.set_xticks(range(900,1651,100)); minor(b)
    b.set_xlabel('TA98 temperature (\u00b0C)'); b.set_ylabel('$T_\\mathrm{TA98}-T_\\mathrm{WZ07}$ (\u00b0C), measured Fe$^{3+}$'); letter(b,'b.')
    kb=[Line2D([],[],color='k',lw=LW_MAIN,label='Experiments: fit, 95 % band'),
        Line2D([],[],marker='o',ls='',mfc='white',mec=EXP_M,mew=1.1,ms=MS_KEY,label='Matjuschkin et al. (2014)'),
        Line2D([],[],marker='s',ls='',mfc='white',mec=EXP_R,mew=1.1,ms=MS_KEY,label='Rzehak et al. (2020)'),
        Line2D([],[],color=XMEAN,lw=LW_REF,ls=(0,(5,2)),label=M(f'Xenoliths: {mu:.0f} \u00b1 {sd:.0f} \u00b0C (n = {n_ex})')),
        Line2D([],[],marker='o',ls='',mfc='none',mec=BLUE_RING,mew=LW_RING,ms=6,label='Blue / red error groups')]
    b.legend(handles=kb,loc='lower center',bbox_to_anchor=(0.5,1.0),ncol=2,columnspacing=0.8,handlelength=1.6,frameon=False,fontsize=FS_SMALL,handletextpad=0.4,labelspacing=0.25,borderaxespad=0.2)
    st=dict(n=int(len(u)),n_excl=n_ex,mean_excl=mu,sd_excl=sd,n_WZ07_gt_TA98=int((u.dTm<0).sum()),
            slave25_mean=float(u.dTm[ex&(u['Craton']=='Slave Craton').to_numpy()].mean()),
            fit=dict(b0=FIT.b0,b1=FIT.b1,se_b0=FIT.se_b0,se_b1=FIT.se_b1,n_eff=FIT.n_eff,chi2_red=FIT.chi2_red))
    log('fig04_all57',st); print('FigX',st)
    u[['SampleID','Craton','T_TA98','WZ0','WZm','dTm','grp']].to_csv(os.path.join(OUTD,'Fig4_per_sample.csv'),index=False)
    return fig,(4110,1984)
if __name__=='__main__':
    f,px=figX(); save(f,'Fig4_final.png',px)
