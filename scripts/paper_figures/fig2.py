from unified import *
import os; D=os.path.join(os.path.dirname(os.path.abspath(__file__)),'data')+'/'
d=pd.read_csv(D+'fig2_curves.csv'); x=d.x.values
TA_HI,TA_LO=1288.893081993608,661.4758937660928
xc=float(np.interp(0,-d.dT_21_2.values,x))
C_HI,C_LO=TA98C,WZ07; C_HI_T,C_LO_T='#C25600',WZ07      # colour = sample (21-2 vermilion, 22-1 blue); TA98 black dashed
def fig2():
    W,H=174,100
    fig,axs=plt.subplots(2,2,figsize=(W*MM,H*MM),sharex=True,gridspec_kw=dict(hspace=0.10,wspace=0.20,height_ratios=[1.6,1]))
    fig.subplots_adjust(left=0.10,right=0.985,top=0.915,bottom=0.115)
    cw=[C_HI,C_LO]; cwt=[C_HI_T,C_LO_T]
    cfg=[('21-2',d.T_WZ07_21_2,d.KD_21_2,TA_HI,(1130,1480)),('22-1',d.T_WZ07_22_1,d.KD_22_1,TA_LO,(600,950))]
    for j,(sid,Tw,KD,ta,yl) in enumerate(cfg):
        t,k=axs[0,j],axs[1,j]; col=cw[j]; colt=cwt[j]
        T0=float(Tw.iloc[0]); T15=float(np.interp(0.15,x,Tw))
        t.fill_between(x,Tw,ta,color=col,alpha=0.16,lw=0)
        t.axhline(T0,color=col,lw=LW_THIN,ls=':',zorder=1)
        t.axhline(ta,color='k',lw=LW_MAIN,ls=(0,(5,2)),zorder=2)
        t.plot(x,Tw,color=col,lw=LW_MAIN,zorder=4)
        t.set_ylim(*yl); t.set_xlim(0,0.2)
        if j==0: t.set_ylabel('Olivine\u2013garnet\ntemperature (\u00b0C)',linespacing=1.1)
        t.annotate('',xy=(0.15,T15),xytext=(0.15,T0),arrowprops=dict(arrowstyle='<->',color=col,lw=LW_ONE,shrinkA=0,shrinkB=0))
        t.text(0.155,T0+0.28*(T15-T0),f'+{T15-T0:.0f} \u00b0C\nat x = 0.15',color=colt,fontsize=FS_ANNOT,fontweight='bold',ha='left',va='center')
        if j==0:
            t.text(0.003,ta+7,f'TA98 = {ta:.0f} \u00b0C\n(not affected by Fe$^{{3+}}$)',color='k',ha='left',va='bottom',fontsize=FS_ANNOT,fontweight='bold',linespacing=1.1)
            t.plot([xc],[ta],'o',ms=8,mfc='white',mec='k',mew=LW_RING,zorder=6)
            t.annotate(f'match at Fe$^{{3+}}$/\u03a3Fe = {xc:.3f}',xy=(xc,ta),xytext=(0.094,1425),fontsize=FS_ANNOT,fontweight='bold',ha='center',va='center',
                       arrowprops=dict(arrowstyle='-|>',color='k',lw=LW_ONE,shrinkA=6,shrinkB=5))
            t.set_title('21-2 (TA98 = %.0f \u00b0C): the thermometers can be matched'%ta,fontsize=FS_ANNOT,fontweight='bold',color=colt,loc='left',pad=5)
        else:
            t.text(0.195,ta-7,f'TA98 = {ta:.0f} \u00b0C (not affected by Fe$^{{3+}}$)',color='k',ha='right',va='top',fontsize=FS_ANNOT,fontweight='bold')
            t.text(0.195,742,'the gap only widens',color=colt,fontsize=FS_ANNOT,fontweight='bold',ha='right',va='center')
            t.set_title('22-1 (TA98 = %.0f \u00b0C): no match is possible'%ta,fontsize=FS_ANNOT,fontweight='bold',color=colt,loc='left',pad=5)
        k.plot(x,KD,color=col,lw=LW_MAIN); k.set_xticks([0,.05,.1,.15,.2]); k.set_xlabel('Fe$^{3+}$/\u03a3Fe removed from garnet (x)')
        if j==0: k.set_ylabel(r'$K_\mathrm{D}$')
        rng=KD.max()-KD.min(); k.set_ylim(KD.min()-0.55*rng,KD.max()+0.12*rng)
        k.text(0.003,KD.min()-0.5*rng,'$K_\\mathrm{D}$ falls as Fe$^{3+}$ is removed from the Fe budget',color=colt,fontsize=FS_SMALL,fontweight='bold',ha='left',va='bottom')
        for ax in (t,k):
            ax.grid(axis='y',color='0.9',lw=0.5,zorder=0); ax.set_axisbelow(True); minor(ax)
        t.xaxis.set_minor_locator(AutoMinorLocator(5)); k.xaxis.set_minor_locator(AutoMinorLocator(5))
        letter(t,'ab'[j]+'.',y=1.09,x=(-0.19 if j==0 else -0.12))
    fig.align_ylabels(axs[:,0])
    return fig,(4110,2362)
if __name__=='__main__':
    f,px=fig2(); save(f,'Fig2_final.png',px)
