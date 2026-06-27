"""
R39 Module 8 -density-paradox test with a PROPER susceptibility control (LSN).

mod07 confound: critical-layer qc1Ncs barely predicted manifestation and Feb PGA
had the wrong sign (Port-Hills rock at high PGA didn't liquefy; floodplain sand at
moderate PGA did). The fix is the Christchurch-standard CPT severity index LSN
(van Ballegooy et al. 2014), which integrates the full-profile triggering response
at each event's demand and groundwater. We:
  1. compute LSN per site for the Feb-2011 demand (and Sept-2010);
  2. validate LSN discriminates observed manifestation (it must, or the proxy is bad);
  3. test whether prior (Sept) liquefaction adds re-liquefaction risk BEYOND LSN:
       P(Feb manifest) ~ LSN_Feb + Sept_manifest.
"""
import numpy as np, pandas as pd, json
import pyarrow.parquet as pq
import xml.etree.ElementTree as ET
from pathlib import Path
from shapely.geometry import shape, Point
from shapely.strtree import STRtree
from pyproj import Transformer
from scipy.spatial import cKDTree
import statsmodels.formula.api as smf
from sklearn.metrics import roc_auc_score

P = Path(__file__).resolve().parents[1]
R39=P; OUT=R39/"outputs"; CACHE=P/"data/raw/public_demand_manifestation_cache_2026-06-08"
rep=[]; w=lambda s="":(rep.append(s),print(s)); Pa=101.325

def crr75(q): q=np.clip(q,1,250); return np.exp(q/113+(q/1000)**2-(q/140)**3+(q/137)**4-2.80)
def rd(z,M):
    a=-1.012-1.126*np.sin(z/11.73+5.133); b=0.106+0.118*np.sin(z/11.28+5.142); return np.exp(a+b*M)
def MSF(q,M):
    mm=np.minimum(2.2,1.09+(np.clip(q,1,250)/180)**3); return 1+(mm-1)*(8.64*np.exp(-M/4)-1.325)
def Ksig(q,sp):
    Cs=np.minimum(0.3,1/(37.3-8.27*np.clip(q,1,211)**0.264)); return np.minimum(1.1,1-Cs*np.log(np.clip(sp,30,None)/Pa))

def site_profile(df):
    qc=df["qc (MPa)"].values*1000.0; fs=df["fs (kPa)"].values; u2=df["u2 (kPa)"].values
    a=df["Area ratio (-)"].values; sv=df["σv_tot (kPa)"].values; svp=np.clip(df["σv_eff (kPa)"].values,5,None)
    z=df["Depth (m)"].values; u0=df["u0 (kPa)"].values
    qt=qc+np.where(u2>-998,u2,0)*(1-a); Fr=np.clip(fs/np.clip(qt-sv,1e-3,None)*100,0.1,10)
    n=np.full_like(z,0.5)
    for _ in range(6):
        Qtn=np.clip((qt-sv)/Pa*(Pa/svp)**n,1,1000); Ic=np.sqrt((3.47-np.log10(Qtn))**2+(1.22+np.log10(Fr))**2)
        n=np.clip(0.381*Ic+0.05*(svp/Pa)-0.15,0.5,1.0)
    qc1N=np.clip(qc,1,None)/Pa
    for _ in range(4):
        m=1.338-0.249*np.clip(qc1N,21,254)**0.264; Cn=np.clip((Pa/svp)**m,None,1.7); qc1N=Cn*(np.clip(qc,1,None)/Pa)
    FC=np.clip(80*Ic-137,0,100); dq=(11.9+qc1N/14.6)*np.exp(1.63-9.7/(FC+2)-(15.7/(FC+2))**2)
    return dict(z=z,Ic=Ic,qc1Ncs=qc1N+dq,sv=sv,svp=svp,u0=u0)

def LSN(prof,amax_g,M):
    z=prof["z"]; q=prof["qc1Ncs"]; svp=prof["svp"]
    liq=(prof["Ic"]<=2.6)&(prof["u0"]>0)&(z>=1)&(z<=20)
    if liq.sum()<3: return np.nan
    csr=0.65*amax_g*(prof["sv"]/svp)*rd(z,M)/(MSF(q,M)*Ksig(q,svp))
    FS=np.clip(crr75(q)/np.clip(csr,1e-3,None),0.01,3)
    # post-liquefaction volumetric strain (Zhang/Yoshimine-style, qc1Ncs-based), %
    Dr=np.clip(0.478*np.clip(q,1,254)**0.264-1.063,0.1,1.0)
    ev=np.where(FS>=2,0.0, np.where(FS<=0.5, np.minimum(0.10,1.5*np.exp(-2.5*Dr)*1.0),
        1.5*np.exp(-2.5*Dr)*np.clip(0.08*(2-FS),0,2)))*100   # percent
    dz=np.gradient(z)
    integ=np.where(liq, ev/np.clip(z,0.5,None)*dz, 0.0)
    return float(1000.0*np.nansum(integ)/1000.0)   # van Ballegooy LSN scaling

# spatial helpers
to_nztm=Transformer.from_crs("EPSG:4326","EPSG:2193",always_xy=True)
def shakemap(fp):
    root=ET.fromstring(Path(fp).read_text(encoding="utf-8",errors="replace")); ln=lambda t:t.split('}',1)[-1]
    fld=[];gd=None
    for el in root.iter():
        nm=ln(el.tag)
        if nm=="grid_field": fld.append((int(el.get("index")),el.get("name")))
        elif nm=="grid_data": gd=el.text
    cols=[n for _,n in sorted(fld)]; arr=np.array([float(x) for x in gd.split()]).reshape(-1,len(cols))
    d=pd.DataFrame(arr,columns=cols); x,y=to_nztm.transform(d["LON"].values,d["LAT"].values)
    pc=[c for c in d.columns if c.upper()=="PGA"][0]; tree=cKDTree(np.c_[x,y]); v=np.asarray(d[pc].values)
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

bridge=pd.read_csv(R39/"data/parquet_coord_bridge.csv")
cols=["Original ID","Area ratio (-)","Depth (m)","qc (MPa)","fs (kPa)","u2 (kPa)","u0 (kPa)","σv_tot (kPa)","σv_eff (kPa)"]
tab=pq.read_table(R39/"data/canterbury_cptu.parquet",columns=cols).to_pandas()
tab=tab[tab["Original ID"].isin(set(bridge["oid"]))]
M_S,M_F=7.1,6.2
rows=[]
for oid,df in tab.groupby("Original ID"):
    prof=site_profile(df.sort_values("Depth (m)"))
    br=bridge[bridge["oid"]==oid].iloc[0]; E,N=br["NZTM Easting"],br["NZTM Northing"]
    aS=pgaS(E,N)/100.0; aF=pgaF(E,N)/100.0
    rows.append(dict(oid=oid,year=int(br["year"]),E=E,N=N,
        LSN_S=LSN(prof,aS,M_S),LSN_F=LSN(prof,aF,M_F),pgaF=aF*100,
        mS=inpoly(tS,pS,E,N),mF=inpoly(tF,pF,E,N)))
d=pd.DataFrame(rows).replace([np.inf,-np.inf],np.nan).dropna(subset=["LSN_F"])
d.to_csv(OUT/"mod08_site_table.csv",index=False)
w("# R39 Module 8 -LSN-controlled density-paradox test\n")
w(f"sites: {len(d)} | LSN_Feb median {d['LSN_F'].median():.1f} (IQR {d['LSN_F'].quantile(.25):.1f}-{d['LSN_F'].quantile(.75):.1f})")

# 2. validate LSN discriminates Feb manifestation
auc=roc_auc_score(d["mF"],d["LSN_F"])
w(f"\n## LSN_Feb vs observed Feb manifestation: AUC = {auc:.3f}  (must be >>0.5 for a valid susceptibility proxy)")
w("  Feb manifest rate by LSN_Feb tertile:")
d["lsn_bin"]=pd.qcut(d["LSN_F"],3,labels=["low","mid","high"],duplicates="drop")
w(d.groupby("lsn_bin",observed=True)["mF"].agg(["mean","size"]).to_string())

# 3. the clean test: does Sept liquefaction add risk BEYOND LSN_Feb-
def run(sub,label):
    if sub["mF"].nunique()<2 or sub["mS"].sum()<5: w(f"\n[{label}] insufficient (n={len(sub)})"); return
    m=smf.logit("mF ~ LSN_F + mS",data=sub).fit(disp=0)
    w(f"\n[{label}] n={len(sub)}")
    for t in ["LSN_F","mS"]:
        ci=m.conf_int().loc[t]
        w(f"    {t:6s} coef={m.params[t]:+.3f} OR={np.exp(m.params[t]):.2f} (95%CI {np.exp(ci[0]):.2f}-{np.exp(ci[1]):.2f}) p={m.pvalues[t]:.3g}")
run(d,"ALL bridged CPTs")
run(d[d["year"]<=2010],"pre-Feb (2010) soundings")
run(d[d["year"]<=2011],"<=2011 soundings")
w("\n-> mS OR>1 AFTER controlling LSN_Feb = prior liquefaction adds risk beyond")
w("   CPT-predicted severity = field evidence the EVENT weakened the ground")
w("   (density paradox). mS~1 = it was just persistent susceptibility (LSN captures it).")
(OUT/"mod08_SUMMARY.md").write_text("\n".join(rep),encoding="utf-8")
w(f"\n[saved] {OUT/'mod08_SUMMARY.md'}")
