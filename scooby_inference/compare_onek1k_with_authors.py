"""Compare our Scooby predictions with the authors' (Zenodo) predictions for the OneK1K eQTL pairs.

Both are variant effects ln(ALT+1) - ln(REF+1) per cell type (22 OneK1K columns). Only pairs that are OneK1K rows
of all_highpip_susie.tsv and are present in both runs are compared. Three checks:
  1. raw values: are the 22 predictions per pair the same? (Pearson r, |difference|, sign agreement)
  2. eQTL scores per OneK1K cell type (as in evaluate_eqtl.py), computed on the same pairs with both prediction sets
  3. top-1 cell type (largest |effect|) and precision at k = 1, the starting point of Fig. 4f

Run (from scooby_inference/):  python compare_onek1k_with_authors.py
Writes eval_onek1k/eval_out_compare/values_per_celltype.csv, scores_per_celltype.csv, pairs.csv.gz
"""
import argparse
import glob
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

p = argparse.ArgumentParser()
p.add_argument("--ours", default="results/gtex_onek1k_pairs/outs_*.h5ad")
p.add_argument("--theirs", default="data/authors_outs/outs_*.h5ad")
p.add_argument("--eqtl", default="data/all_highpip_susie.tsv")
p.add_argument("--cutoff", type=float, default=0.05, help="non-negligible predicted effect")
p.add_argument("--out", default="eval_onek1k/eval_out_compare")
args = p.parse_args()
out = Path(args.out)
out.mkdir(exist_ok=True)
pd.set_option("display.width", 250)


def load(pattern):
    paths = sorted(glob.glob(pattern))
    assert paths, f"no prediction files match {pattern}"
    a = ad.concat([ad.read_h5ad(f) for f in paths])
    df = pd.DataFrame(np.asarray(a.X), columns=a.var_names)
    df["variant"] = a.obs["variant"].astype(str).values
    df["gene_id"] = a.obs["gene_id"].astype(str).values
    print(f"{pattern}: {len(paths)} files, {len(df):,} rows")
    return df.drop_duplicates(["variant", "gene_id"]).set_index(["variant", "gene_id"])


ours, theirs = load(args.ours), load(args.theirs)
cols = [c for c in theirs.columns if c in ours.columns]
print(f"cell types: ours {ours.shape[1]}, theirs {theirs.shape[1]}, shared {len(cols)}")

# OneK1K eQTL pairs present in both runs
e = pd.read_table(args.eqtl).query("study_label == 'OneK1K'").drop_duplicates(["tissue_label", "variant", "gene_id"])
pairs = pd.MultiIndex.from_frame(e[["variant", "gene_id"]].drop_duplicates())
in_o, in_t = pairs.isin(ours.index), pairs.isin(theirs.index)
print(f"OneK1K eQTL pairs: {len(pairs):,}; with our prediction {in_o.sum():,}; with theirs {in_t.sum():,}; "
      f"both {(in_o & in_t).sum():,}")
pairs = pairs[in_o & in_t]
A, B = ours.loc[pairs, cols], theirs.loc[pairs, cols]

# ---------- 1. raw values ----------
diff = (A - B).abs()
big = B.abs() > args.cutoff
v = pd.DataFrame({
    "pearson_r": [np.corrcoef(A[c], B[c])[0, 1] for c in cols],
    "median_abs_diff": diff.median().values,
    "max_abs_diff": diff.max().values,
    "n_clear_ours": (A.abs() > args.cutoff).sum().values,
    "n_clear_theirs": big.sum().values,
    "sign_agree_if_theirs_clear": [(np.sign(A[c][big[c]]) == np.sign(B[c][big[c]])).mean() for c in cols],
}, index=cols)
v.to_csv(out / "values_per_celltype.csv")
print(f"\n=== 1. raw predictions, {len(pairs):,} pairs x {len(cols)} cell types ===")
print(v.round(4).to_string())
x, y = A.to_numpy().ravel(), B.to_numpy().ravel()
same_class = ((np.abs(x) > args.cutoff) == (np.abs(y) > args.cutoff)).mean()
print(f"All values: Pearson r = {np.corrcoef(x, y)[0, 1]:.5f} | |diff| < 0.01 for {(np.abs(x - y) < 0.01).mean():.1%} | "
      f"same clear/not-clear call for {same_class:.1%}")
worst = diff.stack().sort_values(ascending=False).head(5)
print("largest differences:")
print(pd.DataFrame({"ours": A.stack()[worst.index], "theirs": B.stack()[worst.index], "abs_diff": worst}).round(4))

# ---------- 2. eQTL scores per OneK1K cell type, same pairs, both prediction sets ----------
rows = []
for t, g in e.groupby("tissue_label"):
    g = g.set_index(["variant", "gene_id"])
    g = g[g.index.isin(pairs)]
    if len(g) <= 10:  # as evaluate_eqtl.py --min-n 10
        continue
    row = {"cell_type": t, "n": len(g)}
    beta = g["beta"].to_numpy()
    for name, P in (("ours", A), ("theirs", B)):
        pr = P.loc[g.index, t.replace(" ", "_")].to_numpy()
        clear = np.abs(pr) > args.cutoff
        ok = clear & (np.abs(beta) > args.cutoff) & (np.sign(pr) == np.sign(beta))
        row[f"spearman_{name}"] = stats.spearmanr(pr, beta)[0]
        row[f"n_clear_{name}"] = clear.sum()
        row[f"direction_{name}"] = ok.sum() / clear.sum() if clear.sum() else np.nan
    rows.append(row)
s = pd.DataFrame(rows).sort_values("n", ascending=False)
s.to_csv(out / "scores_per_celltype.csv", index=False)
print("\n=== 2. eQTL scores per OneK1K cell type (prediction = column of the same cell type) ===")
print(s.round(3).to_string(index=False))

# ---------- 3. top-1 cell type and precision at k = 1 ----------
truth = e.groupby(["variant", "gene_id"])["tissue_label"].agg(lambda c: set(c.str.replace(" ", "_"))).loc[pairs]
top_o = np.array(cols)[np.abs(A.to_numpy()).argmax(axis=1)]
top_t = np.array(cols)[np.abs(B.to_numpy()).argmax(axis=1)]
hit_o = np.array([c in tr for c, tr in zip(top_o, truth)])
hit_t = np.array([c in tr for c, tr in zip(top_t, truth)])
print(f"\n=== 3. top-1 cell type (largest |effect|), cutoff {args.cutoff} ===")
for name, keep in (("all pairs", np.ones(len(pairs), bool)),
                   (f"max |effect| > {args.cutoff} (ours)", np.abs(A.to_numpy()).max(axis=1) > args.cutoff),
                   (f"max |effect| > {args.cutoff} (theirs)", np.abs(B.to_numpy()).max(axis=1) > args.cutoff)):
    print(f"{name:32s} n = {keep.sum():5d} | same top-1 in both: {(top_o == top_t)[keep].mean():.1%} | "
          f"precision@1 ours {hit_o[keep].mean():.3f}, theirs {hit_t[keep].mean():.3f}")
pd.DataFrame({"variant": pairs.get_level_values(0), "gene_id": pairs.get_level_values(1),
              "true_cell_types": [";".join(sorted(t)) for t in truth],
              "top1_ours": top_o, "top1_theirs": top_t,
              "max_abs_ours": np.abs(A.to_numpy()).max(axis=1), "max_abs_theirs": np.abs(B.to_numpy()).max(axis=1)}
             ).to_csv(out / "pairs.csv.gz", index=False)
print(f"\nWritten: {out}/values_per_celltype.csv, scores_per_celltype.csv, pairs.csv.gz")
