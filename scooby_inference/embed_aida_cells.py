"""Encode AIDA cells with the substitute scPoli encoder, write the two files Scooby needs, and run the placement
check (as check G for CIMA).

Input:  data/aida/cells/AIDA_DataFreeze_v2.h5ad (CELLxGENE; raw integer counts in raw/X, Ensembl IDs in raw/var,
        cell types in obs['Annotation_Level3']), the encoder from fit_scpoli_encoder.py (scpoli_substitute/).
        The 6 LONZA donors (European controls, not part of AIDA's Asian eQTL cohort) are left out.
Genes:  the encoder's 2,403 genes are OneK1K Ensembl IDs; looked up in AIDA by Ensembl ID, then by GENCODE v32 name
        or OneK1K symbol against AIDA's feature_name. Genes found in none are set to 0 (genes_matched.tsv).
Input to the encoder: log(1 + raw count), exactly as in training and in embed_cima_cells.py.
Placement: each cell -> nearest OneK1K column centre (average of the authors' 10 numbers over the column's cells);
        share nearest the mapped column (aida_mapping.csv) and share nearest a column of the same lineage, next to
        the same numbers for OneK1K's own cells (data/cima_cells/embedded_60donors/onek1k_baseline_per_column.csv).

Run (from scooby_inference/):  python embed_aida_cells.py
  depth test:  python embed_aida_cells.py --thin 0.40 --out data/aida/embedded_thin040 --compare_to data/aida/embedded/embedding.pq
Writes data/aida/embedded/: embedding.pq, celltype.pq (one Scooby column per Level3 type), genes_matched.tsv,
  placement_per_type.csv, summary.json
"""
import argparse
import gzip
import json
import os
import re
import time

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch
from torch import nn

p = argparse.ArgumentParser()
p.add_argument("--h5ad", default="data/aida/cells/AIDA_DataFreeze_v2.h5ad")
p.add_argument("--out", default="data/aida/embedded")
p.add_argument("--encoder_dir", default="scpoli_substitute")
p.add_argument("--gtf", default="data/gencode.v32.annotation.sorted.gtf.gz")
p.add_argument("--mapping", default="aida_mapping.csv")
p.add_argument("--onek1k_embedding", default="data/embedding_no_val_genes.pq")
p.add_argument("--onek1k_celltypes", default="data/celltype_fixed.pq")
p.add_argument("--onek1k_baseline", default="data/cima_cells/embedded_60donors/onek1k_baseline_per_column.csv")
p.add_argument("--block", type=int, default=200_000)
p.add_argument("--thin", type=float, default=0,
               help="depth test: keep each raw count with this probability before encoding (e.g. 0.40), 0 = off")
p.add_argument("--compare_to", default="", help="embedding.pq of the full-depth run: how far each cell moved")
p.add_argument("--seed", type=int, default=0)
args = p.parse_args()
os.makedirs(args.out, exist_ok=True)
dev = "cuda" if torch.cuda.is_available() else "cpu"
t0 = time.time()
log = lambda *a: print(f"[{(time.time() - t0) / 60:5.1f} min]", *a, flush=True)
h = h5py.File(args.h5ad, "r")
dec = lambda xs: [x.decode() if isinstance(x, bytes) else x for x in xs]


def cat(group, key):
    o = group[key]
    return np.array(dec(o["categories"][:]), dtype=object)[o["codes"][:]] if isinstance(o, h5py.Group) else np.array(dec(o[:]), dtype=object)


# ---------- genes ----------
g = pd.read_table(f"{args.encoder_dir}/genes.tsv")
g["symbol"] = g["symbol"].str.strip()
v32 = {}
with gzip.open(args.gtf, "rt") as f:
    for line in f:
        if "\tgene\t" in line:
            v32[re.search(r'gene_id "([^"]+)"', line).group(1).split(".")[0]] = re.search(r'gene_name "([^"]+)"', line).group(1)
g["v32_name"] = g["gene_id"].map(v32)
aida_ids = dec(h["raw/var"][h["raw/var"].attrs["_index"]][:])
aida_sym = cat(h["raw/var"], "feature_name")
pos_id = {s: i for i, s in enumerate(aida_ids)}
sym_ok = pd.Series(aida_sym).value_counts()
pos_sym = {s: i for i, s in enumerate(aida_sym) if sym_ok[s] == 1}
how, col = [], []
for gid, a, b in zip(g["gene_id"], g["v32_name"].fillna(""), g["symbol"]):
    if gid in pos_id:
        how.append("Ensembl ID"); col.append(pos_id[gid])
    elif a in pos_sym:
        how.append("GENCODE v32 name"); col.append(pos_sym[a])
    elif b in pos_sym:
        how.append("OneK1K symbol"); col.append(pos_sym[b])
    else:
        how.append("not found, set to 0"); col.append(-1)
g["how_found"], g["aida_col"] = how, col
g.to_csv(f"{args.out}/genes_matched.tsv", sep="\t", index=False)
found = g["aida_col"].values >= 0
log("genes:\n" + g["how_found"].value_counts().to_string())

# ---------- cells: all Level3-labelled cells except the LONZA donors ----------
donor = cat(h["obs"], "donor_id")
label = cat(h["obs"], "Annotation_Level3")
names = np.array(dec(h["obs"][h["obs"].attrs["_index"]][:]), dtype=object)
rows = np.where(~pd.Series(donor).str.startswith("LONZA").values)[0]
log(f"cells: {len(rows):,} of {len(donor):,} kept ({len(donor) - len(rows):,} LONZA cells left out); "
    f"{len(set(donor[rows]))} donors, {len(set(label[rows]))} Level3 types")

# ---------- encoder ----------
ck = torch.load(f"{args.encoder_dir}/encoder.pt", map_location=dev)
layers, d = [], len(g)
for hd in ck["hidden"]:
    layers += [nn.Linear(d, hd), nn.LayerNorm(hd), nn.ReLU()]
    d = hd
net = nn.Sequential(*layers, nn.Linear(d, 10)).to(dev)
net.load_state_dict(ck["state_dict"])
net.eval()
assert list(ck["genes"]) == list(g["gene_id"]), "genes.tsv does not match the encoder"

# ---------- encode ----------
C = h["raw/X"]
n_genes = int(C.attrs["shape"][1])
indptr = C["indptr"][:]
cols_in, cols_out = g["aida_col"].values[found], np.where(found)[0]
rng = np.random.default_rng(args.seed)
Z = np.zeros((len(rows), 10), np.float32)
depth = np.zeros(len(rows), np.float32)
n_file, done = len(indptr) - 1, 0
for s0 in range(0, n_file, args.block):
    s1 = min(s0 + args.block, n_file)
    i0, i1 = np.searchsorted(rows, s0), np.searchsorted(rows, s1)
    if i0 == i1:
        continue
    a, b = indptr[s0], indptr[s1]
    m = sp.csr_matrix((C["data"][a:b], C["indices"][a:b], indptr[s0:s1 + 1] - a),
                      shape=(s1 - s0, n_genes))[rows[i0:i1] - s0][:, cols_in]
    X = np.zeros((i1 - i0, len(g)), np.float32)
    if args.thin:  # binomial thinning: as if the cell had been sequenced less deeply
        m.data = rng.binomial(np.round(m.data).astype(np.int64), args.thin).astype(np.float32)
    X[:, cols_out] = m.toarray()
    assert np.allclose(X, np.round(X)), "raw/X is not integer counts"
    depth[i0:i1] = X.sum(1)
    with torch.no_grad():
        Z[i0:i1] = net(torch.from_numpy(np.log1p(X)).to(dev)).cpu().numpy()
    done += i1 - i0
    log(f"  encoded {done:,} / {len(rows):,}")

pd.DataFrame({"embedding": list(Z), "obs_names": names[rows]}).to_parquet(f"{args.out}/embedding.pq")
lab = label[rows]
pd.DataFrame([{"celltype": t, "cellindex": np.where(lab == t)[0]} for t in sorted(set(lab))]).to_parquet(
    f"{args.out}/celltype.pq")

# ---------- placement ----------
Y = np.vstack(pd.read_parquet(args.onek1k_embedding)["embedding"].values)
ct = pd.read_parquet(args.onek1k_celltypes)
col_names = np.array([c.replace(" ", "_").replace("/", "_") for c in ct["celltype"]])
cent = np.stack([Y[np.asarray(ix)].mean(0) for ix in ct["cellindex"]])
near, near_d = np.zeros(len(Z), int), np.zeros(len(Z), np.float32)
for s in range(0, len(Z), 100_000):
    dist = np.sqrt(((Z[s:s + 100_000, None, :] - cent[None]) ** 2).sum(-1))
    near[s:s + 100_000], near_d[s:s + 100_000] = dist.argmin(1), dist.min(1)
mp = pd.read_csv(args.mapping).set_index("aida_type")
base = pd.read_csv(args.onek1k_baseline).set_index("onek1k_column")
LIN = {"B_cell": "B", "memory_B_cell": "B", "plasmablast": "B", "monocyte": "Mono", "CD16+_monocyte": "Mono",
       "dendritic_cell": "DC", "plasmacytoid_dendritic_cell": "DC", "NK_cell": "NK", "CD56+_NK_cell": "NK",
       "hematopoietic_precursor_cell": "HSPC", "platelet": "platelet"}
near_lin = np.array([LIN.get(c, "T") for c in col_names[near]])
out = []
for t in sorted(set(lab)):
    k = lab == t
    mapped = [c for c in str(mp.at[t, "pred_columns"]).split(";") if c and c != "nan"]
    top = pd.Series(col_names[near[k]]).value_counts(normalize=True)
    out.append({"aida_type": t, "cells": int(k.sum()), "lineage": mp.at[t, "lineage"], "mapped_to": ";".join(mapped),
                "has_eqtls": bool(mp.at[t, "has_eqtls"]), "nearest_most_often": top.index[0],
                "share_nearest_most_often": round(float(top.iloc[0]), 3),
                "share_nearest_mapped": round(float(np.isin(col_names[near[k]], mapped).mean()), 3) if mapped else np.nan,
                "onek1k_own_cells_nearest_own": base.loc[mapped, "own_cells_nearest_own_column"].mean() if len(mapped) == 1 else np.nan,
                "share_nearest_same_lineage": round(float((near_lin[k] == mp.at[t, "lineage"]).mean()), 3),
                "median_distance": round(float(np.median(near_d[k])), 3)})
pl = pd.DataFrame(out)
pl.to_csv(f"{args.out}/placement_per_type.csv", index=False)
mm = pl["share_nearest_mapped"].notna()
summary = {"cells": int(len(rows)), "donors": int(len(set(donor[rows]))), "types": int(len(set(lab))),
           "genes_found": int(found.sum()), "genes_total": len(g),
           "aida_median_counts_over_encoder_genes": float(np.median(depth)),
           "share_nearest_mapped_column": float((pl.loc[mm, "share_nearest_mapped"] * pl.loc[mm, "cells"]).sum() / pl.loc[mm, "cells"].sum()),
           "share_nearest_same_lineage": float((pl["share_nearest_same_lineage"] * pl["cells"]).sum() / pl["cells"].sum()),
           "median_distance_to_nearest_centre": float(np.median(near_d))}
for lin_, x in pl.groupby("lineage"):
    summary[f"same_lineage_{lin_}"] = float((x["share_nearest_same_lineage"] * x["cells"]).sum() / x["cells"].sum())
summary["thin"] = args.thin
if args.compare_to:
    R = np.vstack(pd.read_parquet(args.compare_to)["embedding"].values)
    assert len(R) == len(Z)
    moved = np.linalg.norm(Z - R, axis=1)
    nearR = np.concatenate([np.sqrt(((R[s:s + 100_000, None, :] - cent[None]) ** 2).sum(-1)).argmin(1)
                            for s in range(0, len(R), 100_000)])
    summary.update({"compare_to": args.compare_to, "median_moved": float(np.median(moved)),
                    "p90_moved": float(np.percentile(moved, 90)), "same_nearest_centre": float((nearR == near).mean())})
    pd.DataFrame({"aida_type": lab, "moved": moved, "same_nearest": nearR == near}).groupby("aida_type").agg(
        cells=("moved", "size"), median_moved=("moved", "median"), same_nearest=("same_nearest", "mean")
    ).round(3).to_csv(f"{args.out}/moved_per_type.csv")
json.dump(summary, open(f"{args.out}/summary.json", "w"), indent=1)
pd.set_option("display.width", 250)
print(pl.to_string(index=False))
print(json.dumps(summary, indent=1))
log(f"written to {args.out}/")
