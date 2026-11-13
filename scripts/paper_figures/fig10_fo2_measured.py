"""Fig. 10: pressure vs fO2 calculated from MEASURED garnet Fe3+/SumFe (n = 103). Independent of the match level."""
from unified import *
from common import log
from fig11_fo2_estimated import depth_axis, gpa
from matplotlib.legend_handler import HandlerBase
class HandlerXOnCircle(HandlerBase):
    def create_artists(self,legend,orig,xd,yd,width,height,fontsize,trans):
        cx=0.5*width-xd; cy=0.5*height-yd
        a=Line2D([cx],[cy],marker='o',ls='none',mfc='0.62',mec='0.25',mew=EW_STD,ms=6.5)
        b=Line2D([cx],[cy],marker='x',ls='none',mec='w',mew=LW_THIN+0.1,ms=4.2)
        a.set_transform(trans); b.set_transform(trans); return [a,b]
FAILH=Line2D([0],[0],ls='none')
def draw_fig10(df):
    """Fig. 10 from the pipeline table df (measured Fe3+ only)."""
    from calculations.pt_fo2.buffers import buffer_table, nimis2015_band
    meas=pd.to_numeric(df["garnet_Fe3+"],errors="coerce")
    use=df.loc[meas.notna()&(meas>0)&(meas<0.5)].copy()
    fo2=mf._fo2_from_xfe3(use,use["garnet_Fe3+"]); use=use.assign(fo2_meas=fo2); use=use.loc[use["fo2_meas"].notna()]
    W,H=174,120; fig=newfig(W,H); L,B,AW,AH=14.5,12,144.5,96; ax=axmm(fig,W,H,L,B,AW,AH)
    P=np.linspace(30,70,81); buf=buffer_table(P,df=use); lo,hi=nimis2015_band(P)
    ax.fill_betweenx(P,lo,hi,color='0.87',lw=0,zorder=0)
    ax.plot(buf['dEMOD'],P,color=EMOD,lw=LW_REF,ls='--',zorder=1)
    ax.plot(buf['dFeNi'],P,color=FENI,lw=LW_REF,ls='-.',zorder=1)
    ax.plot(buf['dBSE'],P,color='k',lw=LW_REF,zorder=1)
    ax.axvline(0,color='0.4',lw=LW_THIN,zorder=1)
    cr=use['Craton'].fillna('').astype(str)
    scatter_cratons(ax,use.fo2_meas,use.P_CA91,cr)
    has_other=bool((~cr.isin([c[0] for c in CRATON])).any())
    sp=PARAMS.screening
    d_px=pd.to_numeric(use.get("dT_opx_cpx"),errors="coerce"); d_gt=pd.to_numeric(use.get("dT_pyx_gt"),errors="coerce")
    if not isinstance(d_gt,pd.Series) or d_gt.isna().all(): d_gt=pd.to_numeric(use.get("dT_KR00"),errors="coerce")
    d_ol=pd.to_numeric(use.get("dT_gt_ol"),errors="coerce")
    eq_ok=(d_px.notna()&(d_px<sp.dt_opx_cpx_c)&d_gt.notna()&(d_gt<sp.dt_pyx_gt_c)&d_ol.notna()&(d_ol<sp.dt_gt_ol_c))
    q_fail=(use['screen_4phases']==True)&(use['screen_quality']!=True)   # analytical quality (screening.py)
    fail=~eq_ok|q_fail
    ax.scatter(use.loc[fail,"fo2_meas"],use.loc[fail,"P_CA91"],marker="o",s=40,c="0.62",edgecolors="0.25",linewidths=EW_STD,zorder=4)
    ax.scatter(use.loc[fail,"fo2_meas"],use.loc[fail,"P_CA91"],marker="x",s=22,c="w",linewidths=LW_THIN+0.1,zorder=5)
    _,_,blue,red,_=mf._error_groups(df); blue_u=use.index.isin(df.index[blue]); red_u=use.index.isin(df.index[red])
    ax.scatter(use.loc[blue_u,"fo2_meas"],use.loc[blue_u,"P_CA91"],s=S_RING,facecolors="none",edgecolors=BLUE_RING,linewidths=LW_RING,zorder=6)
    ax.scatter(use.loc[red_u,"fo2_meas"],use.loc[red_u,"P_CA91"],s=S_RING,facecolors="none",edgecolors=RED_RING,linewidths=LW_RING,zorder=6)
    p_ma,f_ma=mf._fo2_depth_mean(use["P_CA91"],use["fo2_meas"])
    ax.plot(f_ma,p_ma,color=MEANC,lw=LW_MAIN,zorder=7)
    ax.invert_yaxis(); ax.set_xlabel('log $f$O$_2$ (\u0394FMQ)'); ax.set_xlim(-6,0); ax.set_ylim(70,30)
    gpa(ax); ax.set_xticks(range(-6,1)); ax.set_xticklabels([M(str(v)) for v in range(-6,1)]); ax.xaxis.set_minor_locator(MultipleLocator(0.5))
    depth_axis(ax)
    ix=lambda arr,p:float(np.interp(p,P,np.asarray(arr,float)))
    ax.text(ix(buf['dFeNi'],58)-0.12,58,'Fe\u2013Ni',color=FENI,fontsize=FS_ANNOT,ha='right',va='center',path_effects=halo(),zorder=9)
    ax.text(ix(buf['dEMOD'],39.5)+0.12,39.5,'EMOG/D',color=EMOD_T,fontsize=FS_ANNOT,ha='left',va='center',path_effects=halo(),zorder=9)
    ax.text(ix(buf['dBSE'],68.2)-0.1,68.2,'BSE',color='k',fontsize=FS_ANNOT,ha='right',va='center',path_effects=halo(2.4,'0.87'),zorder=9)
    ax.text(0.985,0.015,f'n = {len(use)}',transform=ax.transAxes,ha='right',va='bottom',fontsize=FS_ANNOT,path_effects=halo(),zorder=9)
    failh=FAILH
    r1=[kh('slave','Slave'),kh('sib','Siberian'),kh('kaa','Kaapvaal')]+([kh('oth','Other')] if has_other else [])
    bl=Line2D([],[],marker='o',ls='',mfc='none',mec=BLUE_RING,mew=LW_RING,ms=7,label='Blue error group')
    rd=Line2D([],[],marker='o',ls='',mfc='none',mec=RED_RING,mew=LW_RING,ms=7,label='Red error group')
    handles=[];labels=[]
    cells=[[kh('slave','Slave'),bl],[kh('sib','Siberian'),rd],[kh('kaa','Kaapvaal'),kh('mean','Depth Mean')],[failh,kh('blank')]]
    names=[['Slave','Blue error group'],['Siberian','Red error group'],['Kaapvaal','Depth Mean'],['Quality or equilibrium fail',' ']]
    for c in range(4):
        for r in range(2): handles.append(cells[c][r]); labels.append(names[c][r])
    fig.legend(handles=handles,labels=labels,handler_map={FAILH:HandlerXOnCircle()},loc='upper center',bbox_to_anchor=((L+AW/2)/W,1-0.6/H),ncol=4,frameon=False,fontsize=FS_KEY,handletextpad=0.4,columnspacing=2.2,labelspacing=0.5,handlelength=1.6,borderaxespad=0)
    print('Fig10 n',len(use),'blue',int(blue_u.sum()),'red',int(red_u.sum()),'fail',int(fail.sum()),'other',has_other)
    log('fig10',dict(n=int(len(use)),fail=int(fail.sum()),quality_fail=int(q_fail.sum()),equilibrium_fail=int((~eq_ok).sum())))
    assert len(use)==103
    return fig,(4110,2834)
if __name__=='__main__':
    f,px=draw_fig10(load()); save(f,'Fig10_final.png',px)
