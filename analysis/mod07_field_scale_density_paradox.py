"""
R39 Module 7 -FIELD-SCALE density-paradox test from FULLY OPEN data.

Pipeline (no auth anywhere):
  open CPTu profiles (Zenodo 10.5281/zenodo.20839217) -(Year,FinalDepth bridge,
  validated 0.0 m vs 104 ground-truth)-> coordinates (NZGD urban_all) -> join to
  ECan manifestation polygons (Sept2010, Feb2011) + USGS ShakeMaps (PGA).

Question (resolves the density paradox at field scale): controlling for the
critical-layer resistance (qc1Ncs) and the Feb-2011 demand, does liquefying in
Sept 2010 change the probability of re-liquefying in Feb 2011-
  Sept_manifest coef > 0  -> prior liquefaction ADDS risk (net weakening / fabric)
                  < 0  -> net strengthening;  ~0 -> susceptibility persists only.
"""
import numpy as np, pandas as pd, json
import pyarrow.parquet as pq
import xml.etree.ElementTree as ET
from pathlib import Path
from shapely.geometry import shape, Point
from shapely.strtree import STRtree
from pyproj import Transformer
import statsmodels.formula.api as smf
from scipy.spatial import cKDTree

P = Path(__file__).resolve().parents[1]
R39=P; OUT=R39/"outputs"; CACHE=P/"data/raw/public_demand_manifestation_cache_2026-06-08"
rep=[]; w=lambda s="":(rep.append(s),print(s)); Pa=101.325

# ----- 1. CPT processing: Ic + qc1Ncs (Robertson 2009 + Boulanger-Idriss 2014)
def process_sounding(df):
    qc=df["qc (MPa)"].values*1000.0          # kPa
    fs=df["fs (kPa)"].values
    u2=df["u2 (kPa)"].values; a=df["Area ratio (-)"].values
    sv=df["σv_tot (kPa)"].values; svp=df["σv_eff (kPa)"].values
    z=df["Depth (m)"].values; u0=df["u0 (kPa)"].values
    qt=qc+np.where(u2>-998,u2,0)*(1-a)        # cone area correction (u2=-1 sentinel -> skip)
    svp=np.clip(svp,5,None)
    Fr=np.clip(fs/np.clip(qt-sv,1e-3,None)*100,0.1,10)
    # iterative Ic (Robertson)
    n=np.full_like(z,0.5)
    for _ in range(6):
        Qtn=np.clip((qt-sv)/Pa*(Pa/svp)**n,1,1000)
        Ic=np.sqrt((3.47-np.log10(Qtn))**2+(1.22+np.log10(Fr))**2)
        n=np.clip(0.381*Ic+0.05*(svp/Pa)-0.15,0.5,1.0)
    # qc1N, qc1Ncs (B-I 2014, iterate m)
    qc1N=np.clip(qc,1,None)/Pa
    for _ in range(4):
        m=1.338-0.249*np.clip(qc1N,21,254)**0.264
        Cn=np.clip((Pa/svp)**m,None,1.7)
        qc1N=Cn*(np.clip(qc,1,None)/Pa)
    FC=np.clip(80*(Ic+0.0)-137,0,100)
    dq=(11.9+qc1N/14.6)*np.exp(1.63-9.7/(FC+2)-(15.7/(FC+2))**2)
    qc1Ncs=qc1N+dq
    return pd.DataFrame(dict(z=z,Ic=Ic,qc1Ncs=qc1Ncs,u0=u0,svp=svp))

def critical_layer(prof):
    # liquefiable: Ic<=2.6, below water table (u0>0), 1<=z<=15 m
    m=(prof["Ic"]<=2.6)&(prof["u0"]>0)&(prof["z"]>=1)&(prof["z"]<=15)
    s=prof[m]
    if len(s)<3: return np.nan
    return s["qc1Ncs"].quantile(0.20)   # representative loose critical layer (20th pct)

# ----- 2. demand (ShakeMap) + manifestation (polygons)
to_nztm=Transformer.from_crs("EPSG:4326","EPSG:2193",always_xy=True)
def shakemap(fp):
    root=ET.fromstring(Path(fp).read_text(encoding="utf-8",errors="replace")); ln=lambda t:t.split('}',1)[-1]
    fld=[];g=None
    for el in root.iter():
        nm=ln(el.tag)
        if nm=="grid_field": fld.append((int(el.get("index")),el.get("name")))
        elif nm=="grid_data": g=el.text
    cols=[n for _,n in sorted(fld)]; arr=np.array([float(x) for x in g.split()]).reshape(-1,len(cols))
    d=pd.DataFrame(arr,columns=cols); x,y=to_nztm.transform(d["LON"].values,d["LAT"].values)
    pc=[c for c in d.columns if c.upper()=="PGA"][0]
    tree=cKDTree(np.c_[x,y]); v=np.asarray(d[pc].values)
    return lambda E,N: float(v[int(tree.query([float(E),float(N)])[1])])
pgaS=shakemap(CACHE/"darfield_2010_grid.xml"); pgaF=shakemap(CACHE/"christchurch_feb2011_grid.xml")
def polys(name):
    gj=json.loads((R39/f"data/ecan_mapped_liquefaction_{name}.geojson").read_text(encoding="utf-8"))
    pl=[shape(f["geometry"]) for f in gj["features"] if f.get("geometry")]; return STRtree(pl),pl
tS,pS=polys("sept2010"); tF,pF=polys("feb2011")
def inpoly(tree,pl,E,N):
    pt=Point(E,N)
    for i in tree.query(pt):
        if pl[i].contains(pt): return 1
    return 0

# ----- 3. build site table for the 857 bridged soundings
bridge=pd.read_csv(R39/"data/parquet_coord_bridge.csv")  # oid, NZGD ID, E, N, year
w("# R39 Module 7 -field-scale density-paradox test (open data)\n")
w(f"bridged CPTs with coordinates: {len(bridge)}")
need=set(bridge["oid"])
cols=["Sorted ID","Original ID","Year","Area ratio (-)","Depth (m)","qc (MPa)","fs (kPa)","u2 (kPa)","u0 (kPa)","σv_tot (kPa)","σv_eff (kPa)"]
tab=pq.read_table(R39/"data/canterbury_cptu.parquet",columns=cols).to_pandas()
tab=tab[tab["Original ID"].isin(need)]
rows=[]
for oid,df in tab.groupby("Original ID"):
    prof=process_sounding(df.sort_values("Depth (m)"))
    q=critical_layer(prof)
    if not np.isfinite(q): continue
    br=bridge[bridge["oid"]==oid].iloc[0]; E,N=br["NZTM Easting"],br["NZTM Northing"]
    rows.append(dict(oid=oid,year=int(br["year"]),E=E,N=N,qc1Ncs=q,
        mS=inpoly(tS,pS,E,N),mF=inpoly(tF,pF,E,N),
        pgaS=float(pgaS(E,N)),pgaF=float(pgaF(E,N))))
d=pd.DataFrame(rows)
d.to_csv(OUT/"mod07_site_table.csv",index=False)
w(f"sites with computable critical-layer qc1Ncs + spatial join: {len(d)}")
w(f"  qc1Ncs: median {d['qc1Ncs'].median():.0f} (IQR {d['qc1Ncs'].quantile(.25):.0f}-{d['qc1Ncs'].quantile(.75):.0f})")
w(f"  PGA(%g): Sept median {d['pgaS'].median():.0f}, Feb median {d['pgaF'].median():.0f}")
w(f"  manifest: Sept {int(d['mS'].sum())}, Feb {int(d['mF'].sum())}")
w("\ncontingency (Sept x Feb):")
w(d.groupby(['mS','mF']).size().to_string())

# ----- 4. the density-paradox logit
def run(sub,label):
    sub=sub.copy(); sub["lnF"]=np.log(np.clip(sub["pgaF"],3,None)); sub["q"]=sub["qc1Ncs"]
    if sub["mF"].nunique()<2 or sub["mS"].sum()<5: w(f"\n[{label}] insufficient variation (n={len(sub)}, Sept+={int(sub['mS'].sum())})"); return
    try:
        m=smf.logit("mF ~ q + lnF + mS",data=sub).fit(disp=0)
        w(f"\n[{label}] n={len(sub)}  P(Feb|Sept)={sub[sub.mS==1]['mF'].mean():.2f}  P(Feb|noSept)={sub[sub.mS==0]['mF'].mean():.2f}")
        for t in ["q","lnF","mS"]:
            if t in m.params.index:
                ci=m.conf_int().loc[t]
                w(f"    {t:4s} coef={m.params[t]:+.3f} OR={np.exp(m.params[t]):.2f} (95%CI {np.exp(ci[0]):.2f}-{np.exp(ci[1]):.2f}) p={m.pvalues[t]:.3g}")
    except Exception as e: w(f"[{label}] logit err {e}")
run(d,"ALL bridged CPTs (mixed measurement year)")
run(d[d["year"]<=2010],"pre-Feb subset (2010 soundings)")
run(d[d["year"]<=2011],"<=2011 soundings")
w("\n-> mS (Sept_manifest) OR>1 controlling resistance+demand = field evidence for")
w("   net weakening (prior liquefaction adds re-liquefaction risk); OR<1 = strengthening.")
(OUT/"mod07_SUMMARY.md").write_text("\n".join(rep),encoding="utf-8")
w(f"\n[saved] {OUT/'mod07_SUMMARY.md'} + mod07_site_table.csv")
