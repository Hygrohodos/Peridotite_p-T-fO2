from unified import *
from matplotlib.ticker import FormatStrFormatter
from matplotlib.patches import Patch
XM=r'Garnet Fe$^{3{+}}$/$\Sigma$Fe, measured'; YE=r'Garnet Fe$^{3{+}}$/$\Sigma$Fe, estimated'; YR='Fe$^{3{+}}$/$\\Sigma$Fe, measured \u2212 estimated'
def points(ax,xv,yv,df,ok,meas_unused=None):
    cr=df['Craton'].fillna('').astype(str)
    known=cr.isin([c[0] for c in CRATON]); oth=ok&~known
    if oth.any(): ax.scatter(xv[oth],yv[oth],s=S_OTHER,facecolors='none',edgecolors=OTHER_EC,linewidths=EW_OTHER,marker='o',zorder=2)
    for k,mk,c,lab,_ in CRATON:
        sel=ok&(cr==k)
        if sel.any(): ax.scatter(xv[sel],yv[sel],s=S_STD,c=c,marker=mk,edgecolors='k',linewidths=EW_STD,zorder=3)
def rings(ax,xv,yv,blue,red):
    for m,c in ((blue,BLUE_RING),(red,RED_RING)):
        if m is not None and m.any(): ax.scatter(xv[m],yv[m],s=S_RING,facecolors='none',edgecolors=c,linewidths=LW_RING,zorder=4)
def fmt(ax,xt,yt,xf,yf):
    xt=np.round(xt,4)+0.0; yt=np.round(yt,4)+0.0; ax.set_xticks(xt); ax.set_yticks(yt); ax.xaxis.set_major_formatter(FormatStrFormatter(xf)); ax.yaxis.set_major_formatter(FormatStrFormatter(yf))
    minor(ax)
def row(axL,axR,df,meas,T,est,ok,blue=None,red=None):
    band_fill(axL,'1to1'); axL.plot([0,0.35],[0,0.35],color='0.35',lw=LW_ONE,zorder=1)
    sigma=mf._fe3_sigma(df); has_e=ok&sigma.notna()&(sigma>0)
    if has_e.any(): axL.errorbar(meas[has_e],est[has_e],xerr=sigma[has_e],fmt='none',ecolor='0.45',elinewidth=LW_THIN,capsize=1.5,zorder=2)
    points(axL,meas,est,df,ok); rings(axL,meas,est,blue,red)
    axL.set_xlim(0,0.35); axL.set_ylim(0,0.20); axL.set_xlabel(XM); axL.set_ylabel(YE)
    fmt(axL,np.arange(0,0.351,0.05),np.arange(0,0.2001,0.05),'%.2f','%.2f')
    band_fill(axR,'resid'); axR.axhline(0,color='0.35',lw=LW_THIN,zorder=1)
    resid=meas-est; points(axR,T,resid,df,ok); rings(axR,T,resid,blue,red)
    axR.set_xlim(900,1400); axR.set_ylim(-0.12,0.20); axR.set_xlabel(XT_TA98); axR.set_ylabel(YR)
    fmt(axR,np.arange(900,1401,100),np.arange(-0.10,0.2001,0.05),'%d','%.2f')
def key_handles(rings_too=True):
    h=[kh('slave','Slave'),Patch(fc=GREEN_INNER,ec='none',label=r'within $\pm$0.02'),Patch(fc=GREEN_OUTER,ec='none',label=r'within $\pm$0.05')]
    if rings_too:
        h+=[kh('blue','blue error group'),kh('red','red error group')]
    return h
