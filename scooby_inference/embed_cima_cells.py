"""Encode CIMA cells with the substitute scPoli encoder, and write the two files Scooby needs.

Input:  a CIMA .h5ad (raw counts in layers/counts, gene symbols as var names, labels in obs['cell_type_l4'])
        and the encoder from fit_scpoli_encoder.py (scpoli_substitute/encoder.pt, genes.tsv).
Output (--out folder):
  embedding.pq   10 numbers per CIMA cell, same format as data/embedding_no_val_genes.pq
  celltype.pq    cells of each CIMA L4 type, same format as data/celltype_fixed.pq (one Scooby column per type)
  genes_matched.tsv, placement_per_type.csv, summary.json   (the checks below)
Point a config's embedding_path / cell_type_idx_path at the two .pq files to run Scooby with CIMA cells.

Genes: the encoder's 2,403 genes are OneK1K Ensembl IDs. CIMA names genes by symbol, so each ID is looked up
by its GENCODE v32 name first, then by OneK1K's own symbol. Genes found under neither name are set to 0 and
listed in genes_matched.tsv.
Input to the encoder: log(1 + raw count), exactly as in training.

Checks:
  depth      median count per cell over the encoder's genes, CIMA vs OneK1K. The input is not normalised for
             sequencing depth (scPoli's was not either), so a large depth difference shifts the embeddings.
  placement  for each CIMA type: which OneK1K column centre its cells are nearest to (centres = average of the
             authors' 10 numbers over each column's cells), how often that centre is one of the columns the type is
             mapped to in cima_mapping_ontology.csv, and how far cells are from the nearest centre compared with
             OneK1K's own cells.

Run (from scooby_inference/, GPU optional):
  python embed_cima_cells.py --cima data/cima_cells/CIMA_RNA_B_566569cells_36326genes_compressed_20260305.h5ad \
      --out data/cima_cells/embedded_B
  add --max_per_type 5000 to keep at most 5,000 random cells of each type in the Scooby files.
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
p.add_argument("--cima", required=True)
p.add_argument("--out", required=True)
p.add_argument("--encoder_dir", default="scpoli_substitute")
p.add_argument("--gtf", default="data/gencode.v32.annotation.sorted.gtf.gz")
p.add_argument("--onek1k_adata", default="data/scpoli/training_data/onek1k_training_data/OneK1K_only_immune.h5ad",
               help="for the depth check; skipped if missing")
p.add_argument("--onek1k_embedding", default="data/embedding_no_val_genes.pq")
p.add_argument("--onek1k_celltypes", default="data/celltype_fixed.pq")
p.add_argument("--mapping", default="cima_mapping_ontology.csv")
p.add_argument("--max_per_type", type=int, default=0, help="0 = keep all cells")
p.add_argument("--block", type=int, default=200_000, help="cells read and encoded at a time")
p.add_argument("--seed", type=int, default=0)
p.add_argument("--thin", type=float, default=0,
               help="depth test: keep each raw count with this probability before encoding (e.g. 0.76), 0 = off")
p.add_argument("--compare_to", default="",
               help="embedding.pq of another run of the same cells: report how far each cell moved")
args = p.parse_args()
os.makedirs(args.out, exist_ok=True)
dev = "cuda" if torch.cuda.is_available() else "cpu"
t0 = time.time()
log = lambda *a: print(f"[{(time.time() - t0) / 60:5.1f} min]", *a, flush=True)
h = h5py.File(args.cima, "r")
dec = lambda xs: [x.decode() if isinstance(x, bytes) else x for x in xs]

# ---------- genes ----------
g = pd.read_table(f"{args.encoder_dir}/genes.tsv")
g["symbol"] = g["symbol"].str.strip()  # OneK1K symbols carry a leading space
v32 = {}
with gzip.open(args.gtf, "rt") as f:
    for line in f:
        if "\tgene\t" in line:
            gid = re.search(r'gene_id "([^"]+)"', line).group(1).split(".")[0]
            v32[gid] = re.search(r'gene_name "([^"]+)"', line).group(1)
g["v32_name"] = g["gene_id"].map(v32)
cima_genes = pd.Index(dec(h["var"][h["var"].attrs["_index"]][:]))
pos = {s: i for i, s in enumerate(cima_genes)}
g["cima_name"] = [a if a in pos else (b if b in pos else "") for a, b in zip(g["v32_name"].fillna(""), g["symbol"])]
g["cima_col"] = [pos.get(n, -1) for n in g["cima_name"]]
g.to_csv(f"{args.out}/genes_matched.tsv", sep="\t", index=False)
found = g["cima_col"].values >= 0
log(f"genes: {found.sum()} of {len(g)} found in CIMA ({len(g) - found.sum()} set to 0)")

# ---------- cells ----------
obs = h["obs"]
lab = obs["cell_type_l4"]
cats, codes = np.array(dec(lab["categories"][:]), dtype=object), lab["codes"][:]
names = np.array(dec(obs[obs.attrs["_index"]][:]), dtype=object)
rng = np.random.default_rng(args.seed)
keep = []
for c in np.unique(codes[codes >= 0]):
    idx = np.where(codes == c)[0]
    if args.max_per_type and len(idx) > args.max_per_type:
        idx = rng.choice(idx, args.max_per_type, replace=False)
    keep.append(idx)
rows = np.sort(np.concatenate(keep))
log(f"cells: {len(rows):,} of {len(codes):,} kept, {len(keep)} types")

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

# ---------- encode, one block of the file at a time ----------
C = h["layers"]["counts"]
indptr = C["indptr"][:]
cols_in, cols_out = g["cima_col"].values[found], np.where(found)[0]
Z = np.zeros((len(rows), 10), np.float32)
depth = np.zeros(len(rows), np.float32)
n_file, done = len(indptr) - 1, 0
for s0 in range(0, n_file, args.block):  # contiguous stretch of the file, so memory stays bounded
    s1 = min(s0 + args.block, n_file)
    i0, i1 = np.searchsorted(rows, s0), np.searchsorted(rows, s1)
    if i0 == i1:
        continue
    a, b = indptr[s0], indptr[s1]
    m = sp.csr_matrix((C["data"][a:b], C["indices"][a:b], indptr[s0:s1 + 1] - a),
                      shape=(s1 - s0, len(cima_genes)))[rows[i0:i1] - s0][:, cols_in]
    X = np.zeros((i1 - i0, len(g)), np.float32)
    if args.thin:  # binomial thinning: as if the cell had been sequenced less deeply
        m.data = rng.binomial(m.data.astype(np.int64), args.thin).astype(np.float32)
    X[:, cols_out] = m.toarray()
    depth[i0:i1] = X.sum(1)
    with torch.no_grad():
        Z[i0:i1] = net(torch.from_numpy(np.log1p(X)).to(dev)).cpu().numpy()
    done += i1 - i0
    log(f"  encoded {done:,} / {len(rows):,} cells (file rows {s1:,} / {n_file:,})")

# ---------- Scooby files ----------
pd.DataFrame({"embedding": list(Z), "obs_names": names[rows]}).to_parquet(f"{args.out}/embedding.pq")
lab_kept = cats[codes[rows]]
pd.DataFrame([{"celltype": t, "cellindex": np.where(lab_kept == t)[0]} for t in sorted(set(lab_kept))]).to_parquet(
    f"{args.out}/celltype.pq")

# ---------- checks ----------
summary = {"cima_file": args.cima, "thin": args.thin, "cells": int(len(rows)), "types": len(keep),
           "genes_found": int(found.sum()), "genes_total": len(g),
           "cima_median_counts_over_encoder_genes": float(np.median(depth))}
if os.path.exists(args.onek1k_adata):
    o = h5py.File(args.onek1k_adata, "r")
    ovar = pd.Index(dec(o["var"][o["var"].attrs["_index"]][:]))
    oc = ovar.get_indexer(g["gene_id"])
    oip = o["X"]["indptr"][:]
    tot = []
    for s0 in np.sort(rng.choice(len(oip) - 1 - 2500, 20, replace=False)):  # 20 random stretches of 2,500 cells
        a, b = oip[s0], oip[s0 + 2500]
        mo = sp.csr_matrix((o["X"]["data"][a:b], o["X"]["indices"][a:b], oip[s0:s0 + 2501] - a),
                           shape=(2500, len(ovar)))
        tot += list(np.asarray(mo[:, oc].sum(1)).ravel())
    summary["onek1k_median_counts_over_encoder_genes"] = float(np.median(tot))
    log(f"depth: median counts per cell over the encoder genes: CIMA {np.median(depth):.0f}, "
        f"OneK1K {np.median(tot):.0f}")

Y = np.vstack(pd.read_parquet(args.onek1k_embedding)["embedding"].values)
ct = pd.read_parquet(args.onek1k_celltypes)
col_names = [c.replace(" ", "_").replace("/", "_") for c in ct["celltype"]]
cent = np.stack([Y[np.asarray(ix)].mean(0) for ix in ct["cellindex"]])


def nearest(Zs):
    dist = np.sqrt(((Zs[:, None, :] - cent[None]) ** 2).sum(-1))
    return dist.argmin(1), dist.min(1)


ok_lab = np.full(len(Y), -1)
for k, ix in enumerate(ct["cellindex"]):
    ok_lab[np.asarray(ix)] = k
s = rng.choice(np.where(ok_lab >= 0)[0], 100_000, replace=False)
_, base_d = nearest(Y[s])
summary["onek1k_median_distance_to_nearest_centre"] = float(np.median(base_d))
near, near_d = nearest(Z)
mp = pd.read_csv(args.mapping).set_index("tissue_label")["pred_columns"].to_dict()
rows_out = []
for t in sorted(set(lab_kept)):
    m = lab_kept == t
    top = pd.Series(np.array(col_names)[near[m]]).value_counts(normalize=True)
    mapped = [c for c in str(mp.get(t, "")).split(";") if c]
    rows_out.append({"cima_type": t, "cells": int(m.sum()), "mapped_to": ";".join(mapped),
                     "nearest_centre_most_often": top.index[0], "share_of_cells": round(float(top.iloc[0]), 3),
                     "share_nearest_a_mapped_column": round(float(np.isin(np.array(col_names)[near[m]], mapped).mean()), 3)
                     if mapped else None,
                     "median_distance_to_nearest_centre": round(float(np.median(near_d[m])), 3)})
pl = pd.DataFrame(rows_out)
pl.to_csv(f"{args.out}/placement_per_type.csv", index=False)
has_map = pl["share_nearest_a_mapped_column"].notna()
summary["share_nearest_a_mapped_column_all_cells"] = float(
    (pl.loc[has_map, "share_nearest_a_mapped_column"] * pl.loc[has_map, "cells"]).sum() / pl.loc[has_map, "cells"].sum())
summary["cima_median_distance_to_nearest_centre"] = float(np.median(near_d))
if args.compare_to:  # e.g. thinned vs full-depth encoding of the same cells
    ref = pd.read_parquet(args.compare_to)
    assert (ref["obs_names"].values == names[rows]).all(), "--compare_to has different cells or order"
    R = np.vstack(ref["embedding"].values)
    moved = np.linalg.norm(Z - R, axis=1)
    same = float((nearest(R)[0] == near).mean())
    summary.update({"compare_to": args.compare_to, "median_moved": float(np.median(moved)),
                    "p90_moved": float(np.percentile(moved, 90)), "same_nearest_centre": same})
    log(f"compared with {args.compare_to}: cells moved by median {np.median(moved):.4f} "
        f"(90th percentile {np.percentile(moved, 90):.4f}); same nearest centre for {same:.1%} of cells")
json.dump(summary, open(f"{args.out}/summary.json", "w"), indent=1)

pd.set_option("display.width", 250)
print(pl.to_string(index=False))
log(f"distance to the nearest OneK1K centre (median): CIMA {np.median(near_d):.3f}, "
    f"OneK1K's own cells {np.median(base_d):.3f}")
log(f"CIMA cells nearest a column their type is mapped to: {summary['share_nearest_a_mapped_column_all_cells']:.1%}")
log(f"written to {args.out}/: embedding.pq, celltype.pq, genes_matched.tsv, placement_per_type.csv, summary.json")
