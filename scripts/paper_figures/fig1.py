from unified import *
import json
from pathlib import Path
from matplotlib.patches import Patch
ROOT=Path(mf.ROOT); web=ROOT/'data'/'spatial'/'web'
COL={"oceanic":"#d5e4ee","other":"#cfc8bc","orogen":"#9aa0c4","sedimentary_basin":"#8fb9b0","volcanic_province":"#c47a62","craton":"#c4a574"}
LAB={"craton":"Craton / shield","sedimentary_basin":"Sedimentary basin","orogen":"Orogen","volcanic_province":"Volcanic / magmatic","oceanic":"Oceanic crust"}
def fig1():
    W,H=174,82; fig=newfig(W,H)
    l,b,w=14.0,11.5,156.0; h=w*136/360; ax=axmm(fig,W,H,l,b,w,h)
    ax.set_facecolor(COL['oceanic'])
    broad=json.loads((web/'geology_broad.geojson').read_text()); by={f['properties'].get('broad'):f for f in broad['features']}
    for key in ["oceanic","other","orogen","sedimentary_basin","volcanic_province","craton"]:
        f=by.get(key)
        if f is None: continue
        for ring in mf._geo_rings(f.get('geometry')): ax.add_patch(plt.Polygon(ring,closed=True,facecolor=COL[key],edgecolor='none',linewidth=0,zorder=1,joinstyle='round'))
    cr=json.loads((web/'geology_cratons.geojson').read_text())
    for f in cr['features']:
        for ring in mf._geo_rings(f.get('geometry')): ax.add_patch(plt.Polygon(ring,closed=True,facecolor='none',edgecolor='#6b4f2a',linewidth=0.3,zorder=2,alpha=0.7))
    wanted={"Slave Craton":(8,6),"Siberian Craton":(4,-8),"Kaapvaal Craton":(10,-8)}; placed=set()
    for f in cr['features']:
        name=f['properties'].get('name') or ''; group=f['properties'].get('group') or ''
        key=name if name in wanted else (group if group in wanted else None)
        if key is None or key in placed: continue
        if key=="Siberian Craton" and name!="Anabar Shield": continue
        if key=="Slave Craton" and name!="Slave Craton": continue
        if key=="Kaapvaal Craton" and name!="Kaapvaal Craton": continue
        xy=mf._label_xy(f)
        if xy is None: continue
        dx,dy=wanted[key]
        ax.text(xy[0]+dx,xy[1]+dy,key.replace(' Craton',''),fontsize=FS_ANNOT,fontweight='bold',color='#3b2a14',ha='left',va='center',zorder=5,path_effects=halo(2.6))
        placed.add(key)
    pipes=json.loads((web/'kimberlite_localities.geojson').read_text()); by_c={k:{'lon':[],'lat':[]} for k in [c[0] for c in CRATON]+['Other']}; npts=0
    for f in pipes['features']:
        g=f.get('geometry') or {}
        if g.get('type')!='Point': continue
        lon,lat=g['coordinates'][:2]; c=f.get('properties',{}).get('craton') or 'Other'; bk=c if c in by_c else 'Other'
        by_c[bk]['lon'].append(lon); by_c[bk]['lat'].append(lat); npts+=1
    for k,mk,c,lab,_ in CRATON:
        p=by_c[k]
        if p['lon']: ax.scatter(p['lon'],p['lat'],s=S_DENSE,marker=mk,c=c,edgecolors='k',linewidths=EW_DENSE,zorder=4)
    if by_c['Other']['lon']: ax.scatter(by_c['Other']['lon'],by_c['Other']['lat'],s=S_DENSE,marker='o',facecolors='none',edgecolors='#2b2623',linewidths=EW_OTHER,zorder=3)
    ax.set_xlim(-170,190); ax.set_ylim(-58,78); ax.set_aspect('equal',adjustable='box')
    ax.set_xlabel('Longitude'); ax.set_ylabel('Latitude'); minor(ax)
    # key above the map: row 1 geology, row 2 kimberlites (same key style as Figs 6-11)
    kw=dict(loc='center left',frameon=False,fontsize=FS_KEY,handletextpad=0.4,columnspacing=1.2,handlelength=1.6,borderaxespad=0,borderpad=0)
    geo=[Patch(facecolor=COL[k],edgecolor='0.4',linewidth=0.3,label=LAB[k]) for k in ("craton","sedimentary_basin","orogen","volcanic_province","oceanic")]
    kim=[kh('slave','Slave'),kh('sib','Siberian'),kh('kaa','Kaapvaal'),kh('oth','Other')]
    xh,xl_=30.0,32.0
    for yc,hd,hs in ((2.8,'Geology',geo),(7.6,'Kimberlites',kim)):
        fig.text(xh/W,1-yc/H,hd,ha='right',va='center',fontsize=FS_KEY,fontweight='bold')
        fig.legend(handles=hs,ncol=len(hs),bbox_to_anchor=(xl_/W,1-yc/H),**kw)
    print('Fig1 pipes',npts,{k:len(v['lon']) for k,v in by_c.items()},'axes mm',w,round(h,1))
    return fig,(4110,1937)
if __name__=='__main__':
    f,px=fig1(); save(f,'Fig1_final.png',px)
