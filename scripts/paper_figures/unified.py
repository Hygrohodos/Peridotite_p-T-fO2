"""House style: the single default style module for every figure in this repo (approved 14 Nov 2026, t1806u).
Imported by all FigN scripts, Fig. 5 (make_fig05_weightcircles.py), Online Resource 2 (figS_slope_extension.py), and,
through style(), by the legacy and diagnostic plots (make_figures.py, make_figures_rev1.py, fo2_uncertainty.py,
make_buffers_check.py, plot_exp_ta98_wz07.py). A stand-alone copy for other papers (house_style.py, .mplstyle and a
README of the specifications) is in OneDrive Comparative Geothermometry/figure_style/.
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

# =====================================================================================================
# HOUSE STYLE (14 Nov 2026, t1804u). One set of numbers for Figs 1-10 and Online Resource 2.
#   Journal: CMP / Springer artwork rules - widths 84 mm (single column) or 174 mm (double column);
#   sans-serif lettering 8-12 pt with minimal size variance; lines >= 0.3 pt; combination art >= 600 dpi;
#   figure parts labelled with lower-case letters (bold a, b, c, outside the axes at the top left).
#   Type: Liberation Sans (metric clone of Arial). 8 pt for everything (ticks, axis labels, keys,
#   annotations); 9 pt bold for panel letters; 7 pt only for sample IDs in Fig. 5.
#   Colour: Okabe-Ito based and checked for deuteranopia, protanopia and tritanopia (min CIELAB dE >= 17
#   between any two categories drawn together; marker shapes differ as well). Meanings are unchanged:
#   orange = Slave, blue squares = Siberian, purple triangles = Kaapvaal, sky-blue / red rings = error
#   groups, greens = +-0.02 / +-0.05 bands, grey field = Nimis et al. (2015) range.
#   Text in a colour uses a darker shade of it (contrast >= 4.5:1 on white).
# =====================================================================================================
FS_TICK=8.0; FS_LABEL=8.0; FS_KEY=8.0; FS_ANNOT=8.0; FS_SMALL=7.0; FS_LETTER=9.0
# ---------------- line widths (pt) ----------------
LW_FRAME=0.6      # axes frame
LW_TICK=0.6       # tick marks (major); minor 0.5
LW_MAIN=1.5       # primary data curves: fits, means, temperature curves, sweep curves
LW_REF=1.1        # labelled reference curves (BSE, Fe-Ni, EMOG/D, isopleth)
LW_ONE=0.9        # 1:1 lines and regression/fit lines inside scatter panels
LW_THIN=0.6       # zero lines, thresholds, dashed envelopes, error bars, leaders
LW_RING=1.3       # highlight rings
# ---------------- markers ----------------
S_STD=24; EW_STD=0.4          # standard scatter symbols (black edge)
S_DENSE=10; EW_DENSE=0.3      # dense clouds (Fig 3, map)
S_OTHER=16; EW_OTHER=0.6      # 'Other' = open circle
S_RING=80                     # highlight ring size
S_STAR=90                     # calibration / validation stars
MS_KEY=4.8                    # marker size in keys (triangle 5.3, star 8)
# ---------------- colour semantics (colour-blind checked) ----------------
SLAVE='#E69F00'; SLAVE_T='#A36A00'                    # Slave (orange); text shade
SIB='#2C5FA8'; KAA='#CC79A7'; OTHER_EC='#7a7a7a'
GREEN_OUTER='#b7dfb7'; GREEN_INNER='#3d9a4a'          # +-0.05 / +-0.02 bands
G02='#2E8B3E'; G05='#79B97F'; G02_T='#25702F'; G05_T='#3F8A4A'   # same greens as lines/labels
FENI='#1F8A99'; FENI_T='#17707C'; EMOD='#9A7400'; EMOD_T='#876500'; ISO='#8c5a2b'; BAR='#6b4c9a'; MEANC='#222222'
BLUE_RING='#56B4E9'; RED_RING='#E31A1C'
PASS='#4A6FA5'; PASS_E='#1f2f4a'; FAIL='#D7191C'      # Fig 3 screens
WZ07='#0072B2'; TA98C=SLAVE
CRATON=[('Slave Craton','o',SLAVE,'Slave',MS_KEY),('Siberian Craton','s',SIB,'Siberian',MS_KEY),('Kaapvaal Craton','^',KAA,'Kaapvaal',5.3)]
COL1, COL2 = 84.0, 174.0      # CMP single / double column widths (mm)
# axis-label wording used everywhere (temperatures as T with the thermometer as subscript)
XT_TA98='Two-pyroxene $T_\\mathrm{TA98}$ (\u00b0C)'
DT_LABEL='$T_\\mathrm{TA98}-T_\\mathrm{WZ07}$ (\u00b0C)'

def style():
    mpl.rcParams.update({
     'font.family':'sans-serif','font.sans-serif':['Liberation Sans','Arial','DejaVu Sans'],
     'mathtext.fontset':'custom','mathtext.rm':'Liberation Sans','mathtext.it':'Liberation Sans:italic','mathtext.bf':'Liberation Sans:bold','mathtext.default':'regular',
     'font.size':FS_TICK,'axes.labelsize':FS_LABEL,'xtick.labelsize':FS_TICK,'ytick.labelsize':FS_TICK,'legend.fontsize':FS_KEY,
     'axes.titlesize':FS_LABEL,'axes.linewidth':LW_FRAME,'xtick.major.width':LW_TICK,'ytick.major.width':LW_TICK,'xtick.major.size':3,'ytick.major.size':3,
     'xtick.minor.width':0.5,'ytick.minor.width':0.5,'xtick.minor.size':1.6,'ytick.minor.size':1.6,'xtick.major.pad':2.5,'ytick.major.pad':2.5,
     'xtick.direction':'in','ytick.direction':'in','xtick.top':True,'ytick.right':True,'xtick.minor.visible':True,'ytick.minor.visible':True,
     'axes.labelpad':3.0,'pdf.fonttype':42,'ps.fonttype':42,'axes.unicode_minus':True,
     'legend.frameon':False,'legend.handletextpad':0.4,'legend.columnspacing':1.4,'legend.labelspacing':0.4,'legend.handlelength':1.5,
     'legend.borderaxespad':0.4,
     'lines.solid_capstyle':'butt','figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white'})
style()
def axmm(fig,W,H,l,b,w,h): return fig.add_axes([l/W,b/H,w/W,h/H])
def newfig(W,H): return plt.figure(figsize=(W*MM,H*MM))
def pxsize(W,H): return (int(W*MM*600),int(H*MM*600))   # pixel size of a W x H mm figure at 600 dpi
def letter(ax,s,x=None,y=None,dx_mm=None,**k):
    """Panel letter (bold, lower case, no period) just outside the axes at the top left, level with the top
    of the y tick labels. dx_mm: distance of the letter's left edge from the axes' left edge (default -11.5 mm)."""
    s=s.rstrip('.')
    from matplotlib.transforms import offset_copy
    tr=offset_copy(ax.transAxes,fig=ax.figure,x=(dx_mm if dx_mm is not None else -11.5)/25.4*72,y=1.2/25.4*72,units='points')
    return ax.text(0,1,s,transform=tr,fontsize=FS_LETTER,fontweight='bold',va='bottom',ha='left',zorder=20,clip_on=False,**k)
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
       'kaa':dict(marker='^',mfc=KAA,mec='k',mew=EW_STD,ms=5.3),'oth':dict(marker='o',mfc='none',mec=OTHER_EC,mew=EW_OTHER,ms=MS_KEY),
       'cal':dict(marker='*',mfc='k',mec='k',mew=0.3,ms=8.0),'val':dict(marker='*',mfc='none',mec='k',mew=0.8,ms=8.0),
       'blue':dict(marker='o',mfc='none',mec=BLUE_RING,mew=LW_RING,ms=6.6),'red':dict(marker='o',mfc='none',mec=RED_RING,mew=LW_RING,ms=6.6),
       'mean':dict(color=MEANC,lw=LW_MAIN),'blank':dict()}[kind]
    d.update(k); d.setdefault('ls','' if 'marker' in d or kind=='blank' else '-')
    if kind=='blank': label=' '
    return Line2D([],[],label=label,**d)
def key_rows(fig,rows,W,H,cx_mm,top_mm,gap=1.8):
    """Two-row (or one-row) key above the axes. rows = list of lists of handles; columns align across rows."""
    nc=max(len(r) for r in rows); order=[]
    for c in range(nc):
        for r in rows: order.append(r[c] if c<len(r) else kh('blank'))
    return fig.legend(handles=order,loc='upper center',bbox_to_anchor=(cx_mm/W,1-top_mm/H),ncol=nc,frameon=False,fontsize=FS_KEY,
        handletextpad=0.4,columnspacing=gap,labelspacing=0.45,handlelength=1.5,borderaxespad=0)
def scatter_cratons(ax,x,y,cr,use=None,zbase=3):
    """Craton symbols exactly as the repo (colour, marker); 'Other' open circles beneath."""
    cr=cr.fillna('').astype(str)
    known=cr.isin([c[0] for c in CRATON])
    ax.scatter(x[~known],y[~known],s=S_OTHER,facecolors='none',edgecolors=OTHER_EC,linewidths=EW_OTHER,marker='o',zorder=zbase-1)
    for k,mk,c,lab,_ in CRATON:
        m=cr==k; ax.scatter(x[m],y[m],s=S_STD,c=c,marker=mk,edgecolors='k',linewidths=EW_STD,zorder=zbase)
