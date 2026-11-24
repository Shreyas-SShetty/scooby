"""Train a substitute for the authors' scPoli encoder, which was not released.

The authors turned each OneK1K cell into 10 numbers with a trained scPoli model and saved only the result
(data/embedding_no_val_genes.pq). Scooby was fine-tuned on those 10 numbers, so new cells (e.g. CIMA) must be
encoded the same way. We have, for the same 1,267,768 cells, both scPoli's input (raw counts of 2,403 genes) and
its output (the 10 numbers), so we fit a network that maps the first to the second.

Built as in scooby_reproducibility/notebooks/data_preparation/onek1k1/analyze_onek1k_anndata_new.ipynb:
  genes  = highly variable genes (3,000) minus val/test genes (cells 77-79) -> 2,403
  input  = log(1 + raw count), as scPoli's encoder does with recon_loss='nb' (no library-size normalisation)
  target = embedding_no_val_genes.pq (scPoli get_latent(mean=True), cells 87-94)
Not used: scPoli's batch input (the 75 pools). CIMA cells come from none of those pools.

Checks (written to --out):
  Every cell gets an out-of-fold prediction: pools are split into --folds groups, and each group is predicted by
  a network that never saw its pools. These held-out predictions give the R2 per dimension, the distance to the
  authors' value, and whether each cell lands nearest the right cell type. They are also saved as an embedding
  file, so Scooby can be rerun on OneK1K with them (step 5) and compared with the original predictions.
  A final network trained on all pools is saved for encoding CIMA cells (step 6).

Run (from scooby_inference/, GPU, about 30-60 min with most of it reading the 14 GB file):
  python fit_scpoli_encoder.py
"""
import argparse
import json
import os
import time

import anndata as ad
import numpy as np
import pandas as pd
import torch
from torch import nn

D = "data/scpoli/training_data/onek1k_training_data"
p = argparse.ArgumentParser()
p.add_argument("--adata", default=f"{D}/OneK1K_only_immune.h5ad")
p.add_argument("--val_genes", default=f"{D}/val_gene_sequences.csv")
p.add_argument("--test_genes", default=f"{D}/test_gene_sequences.csv")
p.add_argument("--embedding", default="data/embedding_no_val_genes.pq")
p.add_argument("--celltypes", default="data/celltype_fixed.pq")
p.add_argument("--out", default="scpoli_substitute")
p.add_argument("--folds", type=int, default=5)
p.add_argument("--hidden", default="512,256", help="hidden layer sizes")
p.add_argument("--epochs", type=int, default=40)
p.add_argument("--batch", type=int, default=4096)
p.add_argument("--lr", type=float, default=1e-3)
p.add_argument("--seed", type=int, default=0)
args = p.parse_args()
os.makedirs(args.out, exist_ok=True)
dev = "cuda" if torch.cuda.is_available() else "cpu"
t0 = time.time()
log = lambda *a: print(f"[{(time.time() - t0) / 60:5.1f} min]", *a, flush=True)

# ---------- genes ----------
log("reading", args.adata)
a = ad.read_h5ad(args.adata)
drop = set()
for f in [args.val_genes, args.test_genes]:
    drop |= set(pd.read_csv(f, sep="\t", header=None)[3])  # columns: chrom, start, end, gene_id, strand
genes = [g for g in a.var_names[a.var["highly_variable"].astype(bool)] if g not in drop]
assert len(genes) == 2403, f"expected 2,403 genes (scPoli input size), got {len(genes)}"
pd.DataFrame({"gene_id": genes, "symbol": a.var.loc[genes, "GeneSymbol"].values}).to_csv(
    f"{args.out}/genes.tsv", sep="\t", index=False)
log(f"{len(genes)} genes; {a.n_obs:,} cells")

# ---------- input and target ----------
emb = pd.read_parquet(args.embedding)
assert len(emb) == a.n_obs and (emb["obs_names"].values == a.obs_names.values).all(), \
    "embedding rows are not in the same order as the h5ad cells"
Y = np.vstack(emb["embedding"].values).astype(np.float32)
X = a[:, genes].X
X = np.log1p(X.toarray() if hasattr(X, "toarray") else X).astype(np.float32)
pools = a.obs["sample"].astype(str).values
del a
log(f"input {X.shape}, target {Y.shape}")

# cell -> Scooby column (22 columns; cells in none of them are -1)
ct = pd.read_parquet(args.celltypes)
cols = sorted(ct["celltype"])
label = np.full(len(Y), -1)
for r in ct.itertuples():
    label[np.asarray(r.cellindex)] = cols.index(r.celltype)

Xg, Yg = torch.from_numpy(X).to(dev), torch.from_numpy(Y).to(dev)
del X
hidden = [int(h) for h in args.hidden.split(",")]


def make_net():
    layers, d = [], Xg.shape[1]
    for h in hidden:
        layers += [nn.Linear(d, h), nn.LayerNorm(h), nn.ReLU()]
        d = h
    return nn.Sequential(*layers, nn.Linear(d, Y.shape[1])).to(dev)


def predict(net, idx):
    net.eval()
    with torch.no_grad():
        return torch.cat([net(Xg[idx[i:i + 65536]]) for i in range(0, len(idx), 65536)]).cpu().numpy()


def train(train_idx, name):
    """Train on train_idx, keeping 5% of those cells aside to pick the best epoch."""
    rng = np.random.default_rng(args.seed)
    train_idx = rng.permutation(train_idx)
    n_es = len(train_idx) // 20
    es_idx, fit_idx = train_idx[:n_es], train_idx[n_es:]
    torch.manual_seed(args.seed)
    net = make_net()
    opt = torch.optim.Adam(net.parameters(), lr=args.lr)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, args.epochs)
    best, best_state = np.inf, None
    for ep in range(args.epochs):
        net.train()
        perm = torch.from_numpy(rng.permutation(fit_idx)).to(dev)
        for i in range(0, len(perm), args.batch):
            b = perm[i:i + args.batch]
            loss = ((net(Xg[b]) - Yg[b]) ** 2).mean()
            opt.zero_grad()
            loss.backward()
            opt.step()
        sched.step()
        mse = float(((predict(net, es_idx) - Y[es_idx]) ** 2).mean())
        if mse < best:
            best, best_state = mse, {k: v.detach().clone() for k, v in net.state_dict().items()}
        if ep % 5 == 4 or ep == args.epochs - 1:
            log(f"  {name} epoch {ep + 1}: early-stop-set MSE {mse:.4f} (best {best:.4f})")
    net.load_state_dict(best_state)
    return net


# ---------- out-of-fold predictions: each pool predicted by a network that never saw it ----------
uniq = np.random.default_rng(args.seed).permutation(np.unique(pools))
fold_of_pool = {p_: i % args.folds for i, p_ in enumerate(uniq)}
fold = np.array([fold_of_pool[p_] for p_ in pools])
P = np.zeros_like(Y)
for k in range(args.folds):
    te = np.where(fold == k)[0]
    log(f"fold {k + 1}/{args.folds}: {len(te):,} held-out cells from {sum(v == k for v in fold_of_pool.values())} pools")
    P[te] = predict(train(np.where(fold != k)[0], f"fold {k + 1}"), te)

# ---------- how close are the held-out predictions? ----------
r2 = 1 - ((P - Y) ** 2).sum(0) / ((Y - Y.mean(0)) ** 2).sum(0)
err = np.linalg.norm(P - Y, axis=1)
has = label >= 0
cent = np.stack([Y[label == c].mean(0) for c in range(len(cols))])  # each column's centre in the authors' space
nearest = lambda Z: np.argmin(((Z[:, None, :] - cent[None]) ** 2).sum(-1), 1)
near_true, near_pred = nearest(Y[has]), nearest(P[has])
per_col = []
for c, name in enumerate(cols):
    m = label[has] == c
    per_col.append({"column": name, "cells": int(m.sum()),
                    "median_distance_to_authors": float(np.median(err[has][m])),
                    "nearest_centre_right_authors": float((near_true[m] == c).mean()),
                    "nearest_centre_right_ours": float((near_pred[m] == c).mean()),
                    "nearest_centre_same_as_authors": float((near_pred[m] == near_true[m]).mean())})
per_col = pd.DataFrame(per_col)
per_col.to_csv(f"{args.out}/heldout_per_column.csv", index=False)
gap = np.median([np.linalg.norm(cent[i] - cent[j]) for i in range(len(cols)) for j in range(i + 1, len(cols))])
summary = {"cells": int(len(Y)), "genes": len(genes), "folds": args.folds, "hidden": hidden,
           "r2_per_dimension": [round(float(v), 4) for v in r2],
           "r2_mean": round(float(r2.mean()), 4),
           "median_distance_to_authors": round(float(np.median(err)), 4),
           "median_distance_between_column_centres": round(float(gap), 4),
           "nearest_centre_same_as_authors": round(float((near_pred == near_true).mean()), 4)}
json.dump(summary, open(f"{args.out}/heldout_summary.json", "w"), indent=1)
log("held-out R2 per dimension:", np.round(r2, 3), "| mean", round(float(r2.mean()), 3))
log(f"median distance to the authors' 10 numbers: {summary['median_distance_to_authors']} "
    f"(median distance between two column centres: {summary['median_distance_between_column_centres']})")
log(f"cells whose nearest column centre is the same as with the authors' numbers: "
    f"{summary['nearest_centre_same_as_authors']:.1%}")
print(per_col.to_string(index=False))

# out-of-fold embeddings in the same format as embedding_no_val_genes.pq (for the Scooby check, step 5)
pd.DataFrame({"embedding": list(P), "obs_names": emb["obs_names"].values}).to_parquet(
    f"{args.out}/embedding_heldout_predicted.pq")

# ---------- final network on all pools, for CIMA ----------
log("training the final network on all pools")
net = train(np.arange(len(Y)), "final")
torch.save({"state_dict": net.state_dict(), "hidden": hidden, "genes": genes,
            "input": "log1p(raw counts) of the genes in genes.tsv, in that order"}, f"{args.out}/encoder.pt")
log(f"written to {args.out}/: encoder.pt, genes.tsv, heldout_summary.json, heldout_per_column.csv, "
    f"embedding_heldout_predicted.pq")
