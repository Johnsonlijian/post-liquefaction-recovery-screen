"""
R39 Module 5 -Field-scale re-liquefaction susceptibility test (clean observable).

Sidesteps the recovery-rate identifiability problem by using OBSERVED manifestation
(clean, per-event, known demand) instead of a contaminated CPT slope.

Question: does liquefying in Sept 2010 (Darfield) make a site MORE or LESS likely
to re-liquefy in Feb 2011 (Christchurch), controlling for (i) the Feb demand and
(ii) baseline susceptibility (independent CPT qc1Ncs where available)- This is the
field-scale version of the density-paradox controversy (Darby: strengthen;
Ye/Ha: weaken).

Data (all public): ECan mapped-liquefaction polygons (Sept2010, Feb2011);
USGS ShakeMaps (cached grid.xml); NZGD CPT critical-layer qc1Ncs (104 parsed).
"""
import numpy as np, pandas as pd, json
import xml.etree.ElementTree as ET
from pathlib import Path
from shapely.geometry import shape, Point
from shapely.strtree import STRtree
from pyproj import Transformer
import statsmodels.formula.api as smf

BASE = Path(__file__).resolve().parents[1]
PROC=BASE/"data/processed/nzgd_profile_parse_2026-06-08"; R39=BASE
OUT=R39/"outputs"; CACHE=BASE/"data/raw/public_demand_manifestation_cache_2026-06-08"
rep=[]; w=lambda s="":(rep.append(s),print(s))
to_nztm=Transformer.from_crs("EPSG:4326","EPSG:2193",always_xy=True)

def load_shakemap(fp):
    root=ET.fromstring(Path(fp).read_text(encoding="utf-8",errors="replace"))
    ln=lambda t:t.split('}',1)[-1]
    fields=[]; grid=None
    for el in root.iter():
        n=ln(el.tag)
        if n=="grid_field": fields.append((int(el.get("index")),el.get("name")))
        elif n=="grid_data": grid=el.text
    cols=[nm for _,nm in sorted(fields)]
    arr=np.array([float(x) for x in grid.split()]).reshape(-1,len(cols))
    df=pd.DataFrame(arr,columns=cols)
    # columns include LON LAT PGA (%g) ...
    x,y=to_nztm.transform(df["LON"].values,df["LAT"].values)
    df["E"]=x; df["N"]=y
    pgacol=[c for c in df.columns if c.upper()=="PGA"][0]
    return df[["E","N",pgacol]].rename(columns={pgacol:"PGA"})

def sampler(grid_df):
    from scipy.spatial import cKDTree
    tree=cKDTree(grid_df[["E","N"]].values); vals=np.asarray(grid_df["PGA"].values)
    def f(E,N):
        d,i=tree.query([float(E),float(N)],k=1); return float(vals[int(i)])
    return f

w("# R39 Module 5 -field-scale re-liquefaction test\n")
pga_S=sampler(load_shakemap(CACHE/"darfield_2010_grid.xml"))
pga_F=sampler(load_shakemap(CACHE/"christchurch_feb2011_grid.xml"))
w("ShakeMaps loaded (Darfield 2010, Christchurch Feb 2011).")

# manifestation polygons -> STRtree for fast point-in-polygon
def load_polys(name):
    gj=json.loads((R39/f"data/ecan_mapped_liquefaction_{name}.geojson").read_text(encoding="utf-8"))
    polys=[shape(f["geometry"]) for f in gj["features"] if f.get("geometry")]
    return polys, STRtree(polys)
pS,tS=load_polys("sept2010"); pF,tF=load_polys("feb2011")
w(f"polygons: Sept2010 {len(pS)}, Feb2011 {len(pF)}")
# bounding region of mapped liquefaction (analysis domain) from Feb polygons
fxmin=min(p.bounds[0] for p in pF); fxmax=max(p.bounds[2] for p in pF)
fymin=min(p.bounds[1] for p in pF); fymax=max(p.bounds[3] for p in pF)

def in_any(tree, polys, E, N):
    pt=Point(E,N)
    for idx in tree.query(pt):
        g=polys[idx] if not hasattr(idx,'area') else idx
        if g.contains(pt): return 1
    return 0

# -- (A) CPT-controlled sites (independent susceptibility = qc1Ncs critical layer)
inv=pd.read_csv(PROC/"profile_parse_inventory.csv")
inv=inv.dropna(subset=["easting_ags","northing_ags","qc1Ncs_median_1_10m"])
rows=[]
for _,r in inv.iterrows():
    E,N=r["easting_ags"],r["northing_ags"]
    if not (fxmin-2000<E<fxmax+2000 and fymin-2000<N<fymax+2000): continue
    rows.append(dict(nzgd=r["NZGD_ID"],E=E,N=N,qc1Ncs=r["qc1Ncs_median_1_10m"],
        manifest_S=in_any(tS,pS,E,N),manifest_F=in_any(tF,pF,E,N),
        pga_S=float(pga_S(E,N)),pga_F=float(pga_F(E,N))))
cpt=pd.DataFrame(rows)
w(f"\n## (A) CPT-controlled sites in mapped domain: {len(cpt)}")
if len(cpt)>20:
    w(cpt[["manifest_S","manifest_F"]].value_counts().to_string())
    w(f"  P(Feb manifest | Sept manifest)   = {cpt[cpt.manifest_S==1]['manifest_F'].mean():.2f}")
    w(f"  P(Feb manifest | NO Sept manifest)= {cpt[cpt.manifest_S==0]['manifest_F'].mean():.2f}")
    cpt["lnpgaF"]=np.log(cpt["pga_F"].clip(0.05)); cpt["q"]=cpt["qc1Ncs"]
    try:
        m=smf.logit("manifest_F ~ q + lnpgaF + manifest_S",data=cpt).fit(disp=0)
        w("\n  logit P(Feb manifest) ~ qc1Ncs + ln(Feb PGA) + Sept_manifest:")
        for t in ["q","lnpgaF","manifest_S"]:
            if t in m.params.index:
                w(f"    {t:12s} coef={m.params[t]:+.3f}  OR={np.exp(m.params[t]):.2f}  p={m.pvalues[t]:.3g}")
        w("    -> Sept_manifest coef>0 (controlling demand+resistance) = prior liquefaction")
        w("       ADDS re-liquefaction risk beyond static susceptibility (net weakening);")
        w("       <0 = net strengthening; ~0 = susceptibility persists, no extra event effect.")
    except Exception as e: w(f"  logit failed: {e}")
    cpt.to_csv(OUT/"mod05_cpt_controlled_sites.csv",index=False)

(OUT/"mod05_SUMMARY.md").write_text("\n".join(rep),encoding="utf-8")
w(f"\n[saved] {OUT/'mod05_SUMMARY.md'}")
