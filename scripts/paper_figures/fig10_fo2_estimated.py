"""Fig. 10: estimated fO2 vs P for the fully screened xenoliths. Floor samples (garnet Fe3+/SumFe at the 0.01 starting
value) are plotted with the same symbols as the others. Depth means with (solid) and without (dashed) the floor samples.
The estimates are the pipeline columns XFe3_gt_corr and delta_logfO2_FMQ, i.e. at the adopted match level
(make_paper_figures.py)."""
from unified import *
from common import log
def depth_axis(ax):
    """Depth (km) axis on the right of a pressure axis."""
    axr=mf._add_depth_axis(ax,km_interval=50); axr.set_ylabel('Depth (km)',fontsize=FS_LABEL)
    axr.tick_params(which='major',length=3,width=LW_TICK,direction='in'); axr.tick_params(which='minor',length=1.6,width=0.5,direction='in')
    return axr
def gpa(ax):
    """Pressure axis in GPa (data in kbar)."""
    ax.yaxis.set_major_locator(MultipleLocator(5)); ax.yaxis.set_major_formatter(FuncFormatter(lambda v,_:f'{v/10:.1f}')); ax.set_ylabel('Pressure (GPa)')
    ax.yaxis.set_minor_locator(MultipleLocator(1))
def draw_fig10(df):
    """Fig. 10 from the pipeline table df."""
    from calculations.pt_fo2.buffers import buffer_table, nimis2015_band, xenolith_geotherm_C
    use=df[df["screen_pass"]==True].copy()
    W,H=174,190; fig=newfig(W,H); ax=axmm(fig,W,H,14.5,12,144.5,163)
    P=np.linspace(20,80,121); buf=buffer_table(P,df=use); lo,hi=nimis2015_band(P)
    ax.fill_betweenx(P,lo,hi,color='0.87',lw=0,zorder=0)
    ax.plot(buf['dEMOD'],P,color=EMOD,lw=LW_REF,ls='--',zorder=1)
    ax.plot(buf['dFeNi'],P,color=FENI,lw=LW_REF,ls='-.',zorder=1)
    ax.plot(buf['dBSE'],P,color='k',lw=LW_REF,zorder=1)
    ax.axvline(0,color='0.4',lw=LW_THIN,zorder=1)
    cr=use['Craton'].fillna('').astype(str)
    scatter_cratons(ax,use.delta_logfO2_FMQ,use.P_CA91,cr)
    meas=pd.to_numeric(use['garnet_Fe3+'],errors='coerce')
    slave=cr.str.contains('Slave'); sib=cr.str.contains('Siberian'); kaa=cr.str.contains('Kaapvaal')
    has=meas.notna()&(meas>0)&(meas<0.5)&use['delta_logfO2_FMQ'].notna(); cal=has&slave; val=has&(sib|kaa)
    ax.scatter(use.loc[cal,'delta_logfO2_FMQ'],use.loc[cal,'P_CA91'],marker='*',s=S_STAR,c='k',edgecolors='white',linewidths=0.35,zorder=5)
    ax.scatter(use.loc[val,'delta_logfO2_FMQ'],use.loc[val,'P_CA91'],marker='*',s=S_STAR,facecolors='none',edgecolors='k',linewidths=LW_THIN,zorder=5)
    p_iso=np.linspace(30,70,81); xcol='XFe3_gt_corr'
    floor=use[pd.to_numeric(use.get(xcol),errors='coerce')<=0.012]; base=floor if len(floor)>=8 else use
    chem=base.median(numeric_only=True); rows=pd.DataFrame([chem]*len(p_iso)); rows['P_CA91']=p_iso; rows['T_TA98']=xenolith_geotherm_C(p_iso,use)
    iso=mf._fo2_from_xfe3(rows,np.full(len(p_iso),0.01)).to_numpy(dtype=float)
    vis=np.isfinite(iso)&(iso>=-6.85)&(iso<=-0.15)
    ax.plot(iso[vis],p_iso[vis],color=ISO,lw=LW_REF,ls='--',zorder=2)
    xf=pd.to_numeric(use[xcol],errors='coerce'); fl=(xf<=0.0105)
    p_ma,f_ma=mf._fo2_depth_mean(use['P_CA91'],use['delta_logfO2_FMQ'])
    ax.plot(f_ma,p_ma,color=MEANC,lw=LW_MAIN,zorder=6)
    p_mb,f_mb=mf._fo2_depth_mean(use.loc[~fl,'P_CA91'],use.loc[~fl,'delta_logfO2_FMQ'])
    ax.plot(f_mb,p_mb,color=MEANC,lw=LW_MAIN,ls=(0,(3,1.5)),zorder=6)
    log('fig10',dict(n_floor=int(fl.sum()),n=int(len(use)),mean_all=dict(zip(map(float,p_ma),map(float,f_ma))),mean_nonfloor=dict(zip(map(float,p_mb),map(float,f_mb)))))
    unc=calib.fo2_uncertainty(calib.fo2_residuals(df)); cut,w_lo,w_hi=unc['split'],unc['reduced']['sigma_log'],unc['oxidised']['sigma_log']
    p_bar,p_lab=67.6,65.55
    for x0,w,side,lab in ((-6.85,w_lo,'left',f'1\u03c3, \u0394FMQ < {M(f"{cut:.2f}")}'),(-0.15,w_hi,'right',f'1\u03c3, \u0394FMQ \u2265 {M(f"{cut:.2f}")}')):
        xs=(x0,x0+w) if side=='left' else (x0-w,x0)
        ax.plot(xs,[p_bar,p_bar],color=BAR,lw=3.0,solid_capstyle='butt',zorder=7)
        for xe in xs: ax.plot([xe,xe],[p_bar-0.45,p_bar+0.45],color=BAR,lw=LW_ONE+0.2,zorder=7)
        ax.text(xs[1]+0.1 if side=='left' else xs[0]-0.1,p_bar,f'{w:.2f} log',ha='left' if side=='left' else 'right',va='center',fontsize=FS_ANNOT,fontweight='bold',color=BAR,zorder=8,path_effects=halo())
        ax.text(x0,p_lab,lab,ha=side,va='center',fontsize=FS_ANNOT,color='0.2',zorder=8)
    nplot=int(use['delta_logfO2_FMQ'].between(-7,0).__and__(use['P_CA91'].between(30,70)).sum())
    ax.text(0.985,0.992,f"n = {nplot}",transform=ax.transAxes,ha='right',va='top',fontsize=FS_ANNOT,path_effects=halo())
    ax.invert_yaxis(); ax.set_xlabel('log $f$O$_2$ (\u0394FMQ)'); ax.set_ylim(70,30); ax.set_xlim(-7,0)
    gpa(ax); ax.set_xticks(range(-7,1)); ax.set_xticklabels([M(str(v)) for v in range(-7,1)]); ax.xaxis.set_minor_locator(MultipleLocator(0.5))
    depth_axis(ax)
    ax.text(-0.82,43.5,'EMOG/D',color=EMOD_T,fontsize=FS_ANNOT,ha='left',va='center',path_effects=halo(),zorder=9)
    ax.text(-5.2,61.0,'Fe\u2013Ni',color=FENI,fontsize=FS_ANNOT,ha='right',va='center',path_effects=halo(),zorder=9)
    ax.text(-4.35,68.3,'BSE',color='k',fontsize=FS_ANNOT,ha='right',va='center',path_effects=halo(2.4,'0.87'),zorder=9)
    ax.text(-4.33,31.0,'Fe$^{3{+}}$/\u03a3Fe = 0.01',color=ISO,fontsize=FS_ANNOT,ha='left',va='center',path_effects=halo(),zorder=9)
    key_rows(fig,[[kh('slave','Slave'),kh('sib','Siberian'),kh('kaa','Kaapvaal'),kh('oth','Other')],
                  [kh('cal','calibration'),kh('val','validation'),kh('mean','depth mean'),kh('mean','depth mean, floor excluded',ls=(0,(3,1.5)))]],W,H,cx_mm=14.5+144.5/2,top_mm=0.6)
    print('Fig10 n',nplot,'cal',int(cal.sum()),'val',int(val.sum()),'sig',round(w_lo,2),round(w_hi,2),'cut',cut)
    return fig,pxsize(W,H)
if __name__=='__main__':
    f,px=draw_fig10(load()); save(f,'Fig10_final.png',px)
