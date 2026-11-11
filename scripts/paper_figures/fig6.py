from unified import *
import os
def fig6():
    s=pd.read_csv(os.path.join(REPO,'data','derived','dt_match_sweep_corr.csv'))   # read-only
    opt=29; r=s[s.dt==opt].iloc[0]
    W,H=119,78; fig=newfig(W,H)
    ax=axmm(fig,W,H,17,12.5,87,61.5); ax2=ax.twinx()
    ax.set_zorder(ax2.get_zorder()+1); ax.patch.set_visible(False)
    ax.axhline(0,color='0.55',lw=LW_THIN,zorder=1)
    ax.axvline(opt,color='0.35',lw=LW_THIN,ls=(0,(4,2)),zorder=1)
    ax.plot(s.dt,s['mean'],color='k',lw=LW_MAIN,zorder=5)
    ax.set_xlim(-3,78); ax.set_ylim(-0.05,0.05)
    ax.set_xlabel('\u0394T match threshold (\u00b0C)'); ax.set_ylabel('mean Fe$^{3+}$/\u03a3Fe (est \u2212 meas)',labelpad=3)
    ax.set_yticks([-0.04,-0.02,0,0.02,0.04]); ax.set_yticklabels([M(f'{v:.2f}') for v in [-0.04,-0.02,0,0.02,0.04]])
    ax.xaxis.set_major_locator(MultipleLocator(10)); ax.xaxis.set_minor_locator(MultipleLocator(5)); ax.yaxis.set_minor_locator(MultipleLocator(0.01))
    ax2.plot(s.dt,100*s.p05,color=G05,lw=LW_MAIN,ls=(0,(4,1.6)),zorder=3)
    ax2.plot(s.dt,100*s.p02,color=G02,lw=LW_MAIN,zorder=4)
    ax2.set_ylim(0,105); ax2.set_yticks(range(0,101,20)); ax2.yaxis.set_minor_locator(MultipleLocator(10))
    ax2.set_ylabel('Success rate (%)',labelpad=4,color=G02_T)
    ax2.tick_params(axis='y',colors=G02_T,which='both'); ax2.spines['right'].set_color(G02_T)
    ax.tick_params(axis='y',right=False,which='both'); ax.spines['left'].set_color('k')
    ax2.plot([opt],[100*r.p02],marker='o',ms=5.5,mfc='white',mec=G02_T,mew=LW_RING,zorder=8,clip_on=False)
    ax2.annotate(f'{opt} \u00b0C: {100*r.p02:.0f} %',xy=(opt,100*r.p02),xytext=(opt+3.5,79.5),fontsize=FS_ANNOT,color=G02_T,fontweight='bold',ha='left',va='center',zorder=9)
    ax2.text(33,90,'% within \u00b10.05',color=G05_T,fontsize=FS_ANNOT,fontweight='bold',ha='left',va='bottom',path_effects=halo(2.2),zorder=9)
    ax2.text(1,80,'% within \u00b10.02',color=G02_T,fontsize=FS_ANNOT,fontweight='bold',ha='left',va='center',path_effects=halo(2.2),zorder=9)
    x0=38; y0=s['mean'][s.dt==x0].iloc[0]
    ax.text(x0,y0-0.0045,'mean (est \u2212 meas)',color='k',fontsize=FS_ANNOT,fontweight='bold',ha='left',va='top',rotation=-26.7,rotation_mode='anchor',path_effects=halo(2.2),zorder=9)
    ax.text(opt-1.2,-0.0472,'\u0394T = 29 \u00b0C',fontsize=FS_ANNOT,color='0.25',ha='right',va='bottom',zorder=9)
    print('Fig6 checks: mean@0',round(float(s['mean'][s.dt==0].iloc[0]),4),'p02@0',float(s.p02[s.dt==0].iloc[0]),'p05@0',float(s.p05[s.dt==0].iloc[0]),'p02@29',float(r.p02),'mean@29',float(r['mean']))
    return fig,(2811,1842)
if __name__=='__main__':
    f,px=fig6(); save(f,'Fig6_final.png',px)
