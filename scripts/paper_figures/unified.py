"""Shared style for every figure of the manuscript (one set of numbers, imported by all FigN scripts).
Read-only on the repo; output goes to PAPER_FIG_OUT (default scripts/paper_figures/out)."""
import sys; sys.dont_write_bytecode=True
import os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from load import *            # REPO paths, mf (read-only helpers), PARAMS, calib, load()
import matplotlib as mpl; mpl.use('Agg')
import matplotlib.pyplot as plt, numpy as np, pandas as pd
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, FuncFormatter
MM=1/25.4
OUTD=os.environ.get('PAPER_FIG_OUT',os.path.join(os.path.dirname(os.path.abspath(__file__)),'out')); os.makedirs(OUTD,exist_ok=True)

# ---------------- type (Liberation Sans everywhere) ----------------
FS_TICK=7.0; FS_LABEL=7.5; FS_KEY=7.5; FS_ANNOT=7.5; FS_SMALL=7.0; FS_LETTER=9.0
# ---------------- line widths (pt) ----------------
LW_FRAME=0.7      # axes frame
LW_TICK=0.6       # tick marks (major); minor 0.5
LW_MAIN=2.0       # primary data curves: fits, means, matched/unmatched temperature curves, sweep curves
LW_REF=1.4        # labelled reference curves (BSE, Fe-Ni, EMOG/D, isopleth, band edges, TA98 line is MAIN)
LW_ONE=1.0        # 1:1 lines and regression/fit lines inside scatter panels
LW_THIN=0.8       # zero lines, thresholds, dashed envelopes, error bars, dotted guides
LW_RING=1.6       # highlight rings
# ---------------- markers ----------------
S_STD=28; EW_STD=0.4          # standard scatter symbols (black edge)
S_DENSE=12; EW_DENSE=0.3      # dense clouds (Fig 3, map)
S_LARGE=55; EW_LARGE=0.8      # few, individually labelled points (Fig 4)
S_OTHER=18; EW_OTHER=0.6      # 'Other' = open circle
S_RING=100                    # highlight ring size
S_STAR=100                    # calibration / validation stars
MS_KEY=5.0                    # marker size in keys (triangle 5.5, star 8.5)
# ---------------- colour semantics ----------------
SLAVE='#e07b00'; SIB='#4c72b0'; KAA='#c44e9a'; OTHER_EC='#7a7a7a'
GREEN_OUTER='#b7dfb7'; GREEN_INNER='#3d9a4a'          # +-0.05 / +-0.02 bands
G02='#2E8B3E'; G05='#79B97F'; G02_T='#25702F'; G05_T='#3F8A4A'   # same greens as lines/labels
FENI='#1F8A99'; EMOD='#C99A00'; EMOD_T='#9A7400'; ISO='#8c5a2b'; BAR='#6b4c9a'; MEANC='#222222'
BLUE_RING='#1f9ad6'; RED_RING='#e31a1c'
PASS='#4A6FA5'; PASS_E='#1f2f4a'; FAIL='#D7191C'      # Fig 3 screens
WZ07='#0072B2'; TA98C='#D55E00'
CRATON=[('Slave Craton','o',SLAVE,'Slave',MS_KEY),('Siberian Craton','s',SIB,'Siberian',MS_KEY),('Kaapvaal Craton','^',KAA,'Kaapvaal',5.5)]

def style():
    mpl.rcParams.update({
     'font.family':'sans-serif','font.sans-serif':['Liberation Sans','Arial','DejaVu Sans'],
     'mathtext.fontset':'custom','mathtext.rm':'Liberation Sans','mathtext.it':'Liberation Sans:italic','mathtext.bf':'Liberation Sans:bold','mathtext.default':'regular',
     'font.size':FS_TICK,'axes.labelsize':FS_LABEL,'xtick.labelsize':FS_TICK,'ytick.labelsize':FS_TICK,'legend.fontsize':FS_KEY,
     'axes.linewidth':LW_FRAME,'xtick.major.width':LW_TICK,'ytick.major.width':LW_TICK,'xtick.major.size':3,'ytick.major.size':3,
     'xtick.minor.width':0.5,'ytick.minor.width':0.5,'xtick.minor.size':1.6,'ytick.minor.size':1.6,
     'xtick.direction':'in','ytick.direction':'in','xtick.top':True,'ytick.right':True,'xtick.minor.visible':True,'ytick.minor.visible':True,
     'axes.labelpad':3.0,'pdf.fonttype':42,'ps.fonttype':42,
     'lines.solid_capstyle':'butt','figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white'})
style()
def axmm(fig,W,H,l,b,w,h): return fig.add_axes([l/W,b/H,w/W,h/H])
def newfig(W,H): return plt.figure(figsize=(W*MM,H*MM))
def letter(ax,s,x=0.03,y=0.97,**k):
    return ax.text(x,y,s,transform=ax.transAxes,fontsize=FS_LETTER,fontweight='bold',va='top',ha='left',zorder=20,**k)
def save(fig,name,px):
    path=f'{OUTD}/{name}'; fig.savefig(path,dpi=600)
    from PIL import Image
    im=Image.open(path); assert im.size==px,(im.size,px); im.save(path,dpi=(600,600)); print('wrote',path,im.size); plt.close(fig)
def M(t): return t.replace('-','\u2212')
halo=lambda lw=2.4,fc='white':[pe.withStroke(linewidth=lw,foreground=fc)]
def minor(ax,nx=2,ny=2,which='both'):
    if which in('both','x'): ax.xaxis.set_minor_locator(AutoMinorLocator(nx))
    if which in('both','y'): ax.yaxis.set_minor_locator(AutoMinorLocator(ny))
def band_fill(ax,kind):
    """+-0.05 (light) and +-0.02 (dark) agreement bands, drawn without outlines."""
    if kind=='1to1':
        ax.fill_between([0,0.35],[-0.05,0.30],[0.05,0.40],color=GREEN_OUTER,lw=0,zorder=0)
        ax.fill_between([0,0.35],[-0.02,0.33],[0.02,0.37],color=GREEN_INNER,lw=0,zorder=0)
    else:
        ax.axhspan(-0.05,0.05,color=GREEN_OUTER,lw=0,zorder=0); ax.axhspan(-0.02,0.02,color=GREEN_INNER,lw=0,zorder=0)
# ---- key handles (identical look in every figure) ----
def kh(kind,label=None,**k):
    d={'slave':dict(marker='o',mfc=SLAVE,mec='k',mew=EW_STD,ms=MS_KEY),'sib':dict(marker='s',mfc=SIB,mec='k',mew=EW_STD,ms=MS_KEY),
       'kaa':dict(marker='^',mfc=KAA,mec='k',mew=EW_STD,ms=5.5),'oth':dict(marker='o',mfc='none',mec=OTHER_EC,mew=EW_OTHER,ms=MS_KEY),
       'cal':dict(marker='*',mfc='k',mec='k',mew=0.3,ms=8.5),'val':dict(marker='*',mfc='none',mec='k',mew=0.8,ms=8.5),
       'mean':dict(color=MEANC,lw=LW_MAIN),'blank':dict()}[kind]
    d.update(k); d.setdefault('ls','' if 'marker' in d or kind=='blank' else '-')
    if kind=='blank': label=' '
    return Line2D([],[],label=label,**d)
def key_rows(fig,rows,W,H,cx_mm,top_mm,gap=2.2):
    """Two-row (or one-row) key above the axes. rows = list of lists of handles; columns align across rows."""
    nc=max(len(r) for r in rows); order=[]
    for c in range(nc):
        for r in rows: order.append(r[c] if c<len(r) else kh('blank'))
    return fig.legend(handles=order,loc='upper center',bbox_to_anchor=(cx_mm/W,1-top_mm/H),ncol=nc,frameon=False,fontsize=FS_KEY,
        handletextpad=0.4,columnspacing=gap,labelspacing=0.5,handlelength=1.6,borderaxespad=0)
def scatter_cratons(ax,x,y,cr,use=None,zbase=3):
    """Craton symbols exactly as the repo (colour, marker); 'Other' open circles beneath."""
    cr=cr.fillna('').astype(str)
    known=cr.isin([c[0] for c in CRATON])
    ax.scatter(x[~known],y[~known],s=S_OTHER,facecolors='none',edgecolors=OTHER_EC,linewidths=EW_OTHER,marker='o',zorder=zbase-1)
    for k,mk,c,lab,_ in CRATON:
        m=cr==k; ax.scatter(x[m],y[m],s=S_STD,c=c,marker=mk,edgecolors='k',linewidths=EW_STD,zorder=zbase)
