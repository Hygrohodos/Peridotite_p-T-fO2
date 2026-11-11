from unified import *
import os; D=os.path.join(os.path.dirname(os.path.abspath(__file__)),'data')+'/'
T1,T1_T,T2,T2_T='#009E73','#00785A','#B04A88','#9E3A78'
d=pd.read_csv(D+'fig4_exp_TA98_WZ07_Fe3.csv')
R1={"F1100r","F1100o","F1200o","F1200r","F1300r","F14002r","F14002o","MaLm6","MaLm11","MaLm1"}
R2={"F1300o","MaLm12","MaLm9"}
d['trend']=np.where(d.SampleID.isin(R1),1,2); d['study']=np.where(d.Reference.str.startswith('Matj'),'M','R')
fit={}
for n in (1,2):
    g=d[d.trend==n]; a,b=np.polyfit(g.T_TA98,g.dT_Fe2,1); r=g.dT_Fe2-(a*g.T_TA98+b); fit[n]=(a,b,r.std(ddof=2),len(g))
OFFS=-117.63138111725789
LAB={'F1100o':(6,-9,'left'),'F1100r':(6,3,'left'),'F1200o':(6,3,'left'),'F1200r':(6,-9,'left'),'F1300r':(-7,6,'right'),
 'F14002o':(7,-9,'left'),'F14002r':(7,3,'left'),'MaLm6':(7,3,'left'),'MaLm11':(-8,3,'right'),'MaLm1':(7,-7,'left'),
 'F1300o':(8,-14,'left'),'MaLm12':(8,-3,'left'),'MaLm9':(-8,7,'right')}
CT={1:T1,2:T2}
def fig4():
    W,H=119,95
    fig=newfig(W,H); b=fig.add_axes([0.125,0.115,0.855,0.775])
    b.set_xlim(1100,1690); b.set_ylim(-325,185)
    b.axhline(0,color='k',lw=LW_ONE,zorder=1); b.text(1106,8,'TA98 = WZ07',fontsize=FS_ANNOT,fontweight='bold',ha='left',va='bottom')
    b.axhline(OFFS,color='0.25',lw=LW_REF,ls=(0,(1,2)),zorder=1); b.text(1106,OFFS+5,M(f'F1300o: {OFFS:.0f} \u00b0C (dotted)'),fontsize=FS_SMALL,color='0.2',ha='left',va='bottom')
    for n,ls in ((1,'-'),(2,(0,(5,2)))):
        g=d[d.trend==n]; a_,b_,sd,nn=fit[n]; xx=np.linspace(g.T_TA98.min()-40,g.T_TA98.max()+40,50)
        b.fill_between(xx,a_*xx+b_-sd,a_*xx+b_+sd,color=CT[n],alpha=0.20 if n==1 else 0.13,lw=0,zorder=0.6)
        b.plot(xx,a_*xx+b_,color=CT[n],lw=LW_MAIN,ls=ls,zorder=3)
    for st,mk in (('M','o'),('R','s')):
        s=d[d.study==st]
        b.scatter(s.T_TA98,s.dT_Fe2,marker=mk,s=S_LARGE if mk=='o' else S_LARGE*0.9,c=[CT[t] for t in s.trend],edgecolors='k',linewidths=EW_LARGE,zorder=5)
    for _,r in d.iterrows():
        dx,dy,ha=LAB[r.SampleID]
        b.annotate(r.SampleID,(r.T_TA98,r.dT_Fe2),xytext=(dx,dy),textcoords='offset points',fontsize=FS_SMALL,ha=ha,va='center',color='k')
    a1,b1,sd1,n1=fit[1]; a2,b2,sd2,n2=fit[2]
    b.text(1685,172,M(f'Trend 1 (n = {n1}): slope {a1:.3f}; \u00b11\u03c3 band = {sd1:.0f} \u00b0C'),color=T1_T,fontsize=FS_ANNOT,fontweight='bold',ha='right',va='top')
    b.text(1108,-236,M(f'Trend 2 (n = {n2}): slope {a2:.3f};\n\u00b11\u03c3 band = {sd2:.0f} \u00b0C (1 d.f.)'),color=T2_T,fontsize=FS_ANNOT,fontweight='bold',ha='left',va='top',linespacing=1.2)
    b.set_xlabel('TA98 temperature (\u00b0C)'); b.set_ylabel('$T_\\mathrm{TA98}-T_\\mathrm{WZ07}$ (\u00b0C)\nWZ07 with measured Fe$^{3+}$')
    b.set_xticks(range(1100,1700,100)); minor(b)
    Mj=Line2D([],[],marker='o',ls='',mfc='0.75',mec='k',mew=EW_LARGE,ms=6,label='Matjuschkin et al. (2014), 5 GPa')
    Rz=Line2D([],[],marker='s',ls='',mfc='0.75',mec='k',mew=EW_LARGE,ms=6,label='Rzehak et al. (2020), 10 GPa')
    t1=Line2D([],[],color=T1,lw=LW_MAIN,marker='o',mfc=T1,mec='k',mew=EW_LARGE,ms=6,label='Trend 1')
    t2=Line2D([],[],color=T2,lw=LW_MAIN,ls=(0,(3,1.5)),marker='o',mfc=T2,mec='k',mew=EW_LARGE,ms=6,label='Trend 2')
    key_rows(fig,[[Mj,t1],[Rz,t2]],W,H,cx_mm=0.55*W,top_mm=0.5,gap=1.6)
    return fig,(2811,2244)
if __name__=='__main__':
    f,px=fig4(); save(f,'Fig4_final.png',px); print(fit)
