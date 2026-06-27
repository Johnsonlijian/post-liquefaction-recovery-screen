"""
R39 Module 6 -within-location, demand-referenced recurrence (manifestation-only,
large-N, no per-site CPT needed).

Logic: at a fixed location geology is constant, so comparing the two events
controls static susceptibility. If a cell liquefied at the LOWER demand of one
event but NOT at the HIGHER demand of the other, susceptibility must have CHANGED
between events. Tally the direction:
  WEAKENED  : manifest at lower-demand event, not at higher-demand event
  STRENGTHENED/RECOVERED: manifest at higher-demand event only is demand-consistent;
              the diagnostic of strengthening is liquefied-then-not at HIGHER later demand.
Honest confound (reported): Feb 2011 was mapped more thoroughly than Sept 2010,
so Sept "absence" is partly mapping incompleteness. We bound this.
"""
import numpy as np, pandas as pd, json
import xml.etree.ElementTree as ET
from pathlib import Path
from shapely.geometry import shape, Point
from shapely.strtree import STRtree
from shapely import prepared
from pyproj import Transformer
BASE = Path(__file__).resolve().parents[1]
R39=BASE; OUT=R39/"outputs"; CACHE=BASE/"data/raw/public_demand_manifestation_cache_2026-06-08"
rep=[]; w=lambda s="":(rep.append(s),print(s))
to_nztm=Transformer.from_crs("EPSG:4326","EPSG:2193",always_xy=True)

def load_shakemap(fp):
    root=ET.fromstring(Path(fp).read_text(encoding="utf-8",errors="replace")); ln=lambda t:t.split('}',1)[-1]
    fields=[]; grid=None
    for el in root.iter():
        n=ln(el.tag)
        if n=="grid_field": fields.append((int(el.get("index")),el.get("name")))
        elif n=="grid_data": grid=el.text
    cols=[nm for _,nm in sorted(fields)]; arr=np.array([float(x) for x in grid.split()]).reshape(-1,len(cols))
    df=pd.DataFrame(arr,columns=cols); x,y=to_nztm.transform(df["LON"].values,df["LAT"].values)
    df["E"]=x; df["N"]=y; pc=[c for c in df.columns if c.upper()=="PGA"][0]
    from scipy.spatial import cKDTree
    tree=cKDTree(df[["E","N"]].values); v=df[pc].values
    return lambda E,N: v[tree.query(np.column_stack([E,N]),k=1)[1]]

pgaS=load_shakemap(CACHE/"darfield_2010_grid.xml"); pgaF=load_shakemap(CACHE/"christchurch_feb2011_grid.xml")
def load_prep(name):
    gj=json.loads((R39/f"data/ecan_mapped_liquefaction_{name}.geojson").read_text(encoding="utf-8"))
    polys=[shape(f["geometry"]) for f in gj["features"] if f.get("geometry")]
    return STRtree(polys), polys
tS,pS=load_prep("sept2010"); tF,pF=load_prep("feb2011")
def inpoly(tree,polys,E,N):
    pt=Point(E,N)
    for idx in tree.query(pt):
        if polys[idx].contains(pt): return 1
    return 0

# grid over union bounds, 100 m
allb=[p.bounds for p in pS+pF]
xmin=min(b[0] for b in allb); xmax=max(b[2] for b in allb); ymin=min(b[1] for b in allb); ymax=max(b[3] for b in allb)
xs=np.arange(xmin,xmax,100.0); ys=np.arange(ymin,ymax,100.0)
w("# R39 Module 6 -demand-referenced recurrence (manifestation-only)\n")
w(f"grid: {len(xs)}x{len(ys)} = {len(xs)*len(ys)} cells @100 m over mapped domain")
rows=[]
for X in xs:
    Es=np.full_like(ys,X)
    pgs=pgaS(Es,ys); pgf=pgaF(Es,ys)
    for k,Y in enumerate(ys):
        mS=inpoly(tS,pS,X,Y); mF=inpoly(tF,pF,X,Y)
        rows.append((X,Y,mS,mF,float(pgs[k]),float(pgf[k])))
g=pd.DataFrame(rows,columns=["E","N","mS","mF","pgaS","pgaF"])
# restrict to cells with non-trivial demand in at least one event (in the affected zone)
g=g[(g["pgaS"]>5)|(g["pgaF"]>5)].copy()   # PGA in %g
w(f"analysis cells (PGA>0.05g in either event): {len(g)}")
w(f"  manifest counts: Sept {int(g['mS'].sum())}, Feb {int(g['mF'].sum())}")
w(f"  PGA (%g): Sept median {g['pgaS'].median():.0f}, Feb median {g['pgaF'].median():.0f}")

# demand-referenced diagnostic
g["feb_higher"]=g["pgaF"]>g["pgaS"]
w("\n## Demand-referenced recurrence diagnostic")
# Among cells where Sept demand was HIGHER than Feb (feb_higher=False):
sub=g[~g["feb_higher"]]
w(f"\ncells with Sept PGA > Feb PGA (n={len(sub)}): a fair test of strengthening/recovery")
w(f"   of those that liquefied in Sept ({int(sub['mS'].sum())}), fraction that did NOT re-liquefy in Feb (lower demand): "
  f"{1-sub[sub.mS==1]['mF'].mean():.2f}  <- demand-consistent non-recurrence (not necessarily recovery)")
# WEAKENING diagnostic: liquefied at LOWER demand event but not at HIGHER demand event
weak=g[(g.mS==0)&(g.mF==1)&(g.pgaF<g.pgaS)]      # Feb liq at lower Feb demand, no Sept liq at higher Sept demand
strong=g[(g.mS==1)&(g.mF==0)&(g.pgaF>g.pgaS)]    # Sept liq, no Feb liq at higher Feb demand -> strengthened/recovered
w(f"\nWEAKENING cells (Feb-only, Feb PGA<Sept PGA): {len(weak)}")
w(f"STRENGTHENED cells (Sept-only, Feb PGA>Sept PGA): {len(strong)}")
w(f"   ratio weak:strong = {len(weak)}:{len(strong)}")
w("\nCAVEAT (honest): Feb 2011 was mapped more thoroughly than Sept 2010, so")
w(f"Sept absence is partly mapping incompleteness -> the WEAKENING count is an")
w("UPPER bound and the STRENGTHENED count a lower bound; the asymmetry cannot be")
w("cleanly attributed to physical change without per-event mapping-effort control.")
g.to_csv(OUT/"mod06_recurrence_grid.csv",index=False)
(OUT/"mod06_SUMMARY.md").write_text("\n".join(rep),encoding="utf-8")
w(f"\n[saved] {OUT/'mod06_SUMMARY.md'}")
