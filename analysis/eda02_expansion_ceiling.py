"""
R39 EDA 02 -Compute the TRUE data-expansion ceiling from NZGD metadata only.

Decisive question: if we downloaded everything, how many coordinate clusters
have >=2 distinct sounding TIMES (rate-identifying) across the Canterbury event
boundaries, at radii 2/5/10 m- This bounds N_rate and tells us which NZGD IDs
to target. Pure metadata; no downloads.
"""
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.spatial import cKDTree

BASE = Path(__file__).resolve().parents[1]
Q = BASE / "data/processed/nzgd_first_query_2026-06-08"
OUT = BASE / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(Q / "urban_all.csv", encoding="utf-8-sig")
df.columns = [c.strip() for c in df.columns]
print("urban_all rows:", len(df))
print("columns:", list(df.columns))

# Keep CPTs with coordinates and dates
df = df[df["Investigation Type"].str.contains("CPT", na=False)].copy()
df["E"] = pd.to_numeric(df["NZTM Easting"], errors="coerce")
df["N"] = pd.to_numeric(df["NZTM Northing"], errors="coerce")
df["date"] = pd.to_datetime(df["End Date"], errors="coerce", dayfirst=True)
df = df.dropna(subset=["E","N","date"]).copy()
print("CPT with coord+date:", len(df))
print("date range:", df["date"].min().date(), "->", df["date"].max().date())

# Canterbury sequence event boundaries
EV = {
    "Darfield_2010": pd.Timestamp("2010-09-04"),
    "Christchurch_2011-02": pd.Timestamp("2011-02-22"),
    "June_2011": pd.Timestamp("2011-06-13"),
    "December_2011": pd.Timestamp("2011-12-23"),
}
bounds = [pd.Timestamp("1900-01-01")] + list(EV.values()) + [pd.Timestamp("2100-01-01")]
labels = ["pre_Darfield","Darfield_to_Feb","Feb_to_Jun","Jun_to_Dec","post_Dec"]
df["epoch"] = pd.cut(df["date"], bounds, labels=labels, right=False)

print("\nCPT count by epoch:")
print(df["epoch"].value_counts().reindex(labels).to_string())

# Spatial clustering via union-find on KD-tree pairs within radius
def cluster_within(df, radius):
    pts = df[["E","N"]].values
    tree = cKDTree(pts)
    pairs = tree.query_pairs(radius, output_type="ndarray")
    parent = np.arange(len(df))
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for a,b in pairs:
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb
    roots = np.array([find(i) for i in range(len(df))])
    return roots

idx = df.reset_index(drop=True)
print("\n" + "="*70)
print("EXPANSION CEILING by radius")
print("="*70)
summary = []
for r in [2.0, 5.0, 10.0]:
    roots = cluster_within(idx, r)
    idx2 = idx.copy(); idx2["cluster"] = roots
    g = idx2.groupby("cluster")
    n_dates = g["date"].nunique()
    n_epochs = g["epoch"].nunique()
    # rate-identifying: >=2 distinct dates in same cluster
    ri2 = (n_dates >= 2)
    ri3 = (n_dates >= 3)
    # spanning a liquefaction event (>=2 epochs) -true pre/post recovery
    spanning = (n_epochs >= 2) & (n_dates >= 2)
    # clusters that straddle the post-Feb period specifically (Feb_to_Jun or later with an earlier sounding)
    clusters_multi = idx2[idx2["cluster"].isin(n_dates[ri2].index)]
    inv_ids = clusters_multi["NZGD ID"].nunique()
    summary.append(dict(radius_m=r,
                        clusters_total=int(g.ngroups),
                        clusters_ge2dates=int(ri2.sum()),
                        clusters_ge3dates=int(ri3.sum()),
                        clusters_spanning_event=int(spanning.sum()),
                        investigations_in_multidate_clusters=int(inv_ids)))
    print(f"\nradius={r} m:")
    print(f"  total clusters: {g.ngroups}")
    print(f"  clusters with >=2 distinct dates (rate-identifying): {ri2.sum()}")
    print(f"  clusters with >=3 distinct dates: {ri3.sum()}")
    print(f"  clusters spanning >=2 epochs AND >=2 dates (pre/post event): {spanning.sum()}")
    print(f"  NZGD investigations inside the >=2-date clusters: {inv_ids}")

pd.DataFrame(summary).to_csv(OUT/"expansion_ceiling_summary.csv", index=False)
print(f"\n[saved] {OUT/'expansion_ceiling_summary.csv'}")

# For the 5 m radius, dump the multi-date clusters spanning an event with their NZGD IDs
roots = cluster_within(idx, 5.0)
idx["cluster"] = roots
g = idx.groupby("cluster")
keep = g.filter(lambda s: s["date"].nunique() >= 2 and s["epoch"].nunique() >= 2)
keep = keep.sort_values(["cluster","date"])
cl_summary = (keep.groupby("cluster")
              .agg(n_soundings=("NZGD ID","nunique"),
                   n_dates=("date","nunique"),
                   n_epochs=("epoch","nunique"),
                   first_date=("date","min"),
                   last_date=("date","max"),
                   E=("E","mean"), N=("N","mean"),
                   nzgd_ids=("NZGD ID", lambda s: ";".join(map(str,sorted(set(s))))))
              .reset_index())
cl_summary["span_days"] = (cl_summary["last_date"]-cl_summary["first_date"]).dt.days
cl_summary = cl_summary.sort_values(["n_dates","span_days"], ascending=False)
cl_summary.to_csv(OUT/"multidate_event_spanning_clusters_5m.csv", index=False)
print(f"[saved] {OUT/'multidate_event_spanning_clusters_5m.csv'}  ({len(cl_summary)} clusters)")
print("\nTop 15 multi-date event-spanning clusters (5 m):")
print(cl_summary.head(15)[["n_soundings","n_dates","n_epochs","first_date","last_date","span_days"]].to_string())

# How many of these investigations are already downloaded-
have = set()
cov = pd.read_csv(Q/"nzgd_primary_download_coverage_2026-06-08.csv")
allzip = cov.loc[cov["status"]=="all_zip_investigations","ids"].iloc[0]
have = set(int(x) for x in str(allzip).split(";") if x.strip().isdigit())
need_ids = set()
for s in cl_summary["nzgd_ids"]:
    need_ids |= set(int(x) for x in s.split(";"))
print(f"\ninvestigations in 5m event-spanning multi-date clusters: {len(need_ids)}")
print(f"  already downloaded: {len(need_ids & have)}")
print(f"  NOT yet downloaded: {len(need_ids - have)}")
pd.Series(sorted(need_ids - have)).to_csv(OUT/"expansion_target_nzgd_ids.csv", index=False, header=["nzgd_id"])
print(f"[saved] expansion_target_nzgd_ids.csv ({len(need_ids-have)} ids)")
