"""Compare scooby's predicted variant effects with measured eQTL effects (reproduces the paper's Figure 4).

Method copied from the authors' notebooks/figure4/Figure4.ipynb:
  * eQTL rows are taken per (study, tissue); duplicates of the same variant-gene pair are dropped, keeping the
    row from the study with the largest sample size.
  * Predictions are joined on variant + gene_id. Predicted effect = mean of the matching scooby cell-type
    columns: one column for a OneK1K cell type; for GTEx whole blood, the mean over all OneK1K cell types.
  * Distance to TSS uses GENCODE v32 gene starts (gene end for - strand genes).
  * Effects are natural-log fold changes; "non-negligible" means |predicted effect| > 0.05.

Run (reproduce the paper):
  python evaluate_eqtl.py
For other data (e.g. CIMA later), pass --pred, --eqtl and a --mapping CSV with columns
  study_label, tissue_label, pred_columns   (pred_columns = scooby column names separated by ';').
Outputs go to eval_out/: summary.csv (one row per study/tissue), distance_bins.csv, merged.csv.gz.
"""
import argparse
import glob
import gzip
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

p = argparse.ArgumentParser()
p.add_argument("--pred", default="data/authors_outs/outs_*.h5ad", help="glob of scooby output .h5ad files")
p.add_argument("--eqtl", default="data/all_highpip_susie.tsv")
p.add_argument("--gtf", default="data/gencode.v32.annotation.sorted.gtf.gz")
p.add_argument("--mapping", default=None, help="CSV: study_label,tissue_label,pred_columns (default: authors')")
p.add_argument("--cutoff", type=float, default=0.05, help="non-negligible predicted effect")
p.add_argument("--min-n", type=int, default=10, help="report study/tissue groups with more than this many pairs")
p.add_argument("--out", default="eval_out")
args = p.parse_args()
out = Path(args.out)
out.mkdir(exist_ok=True)

# ---------- 1. predictions: all chunks, rows = variant-gene pairs, columns = cell types ----------
paths = sorted(glob.glob(args.pred), key=lambda x: int(re.findall(r"_(\d+)_\d+\.h5ad$", x)[0]))
assert paths, f"no prediction files match {args.pred}"
pred = ad.concat([ad.read_h5ad(f) for f in paths], merge="same")
pred_df = pd.DataFrame(pred.X, columns=pred.var_names)
pred_df["variant"] = pred.obs["variant"].astype(str).values
pred_df["gene_id"] = pred.obs["gene_id"].astype(str).values
n_dup = pred_df.duplicated(["variant", "gene_id"]).sum()
print(f"predictions: {len(paths)} files, {len(pred_df):,} rows ({n_dup} duplicated variant-gene rows, kept as in "
      f"the authors' notebook), {pred.n_vars} cell types")

# ---------- 2. measured eQTLs + distance to TSS ----------
eqtl = pd.read_table(args.eqtl)
genes = []
with gzip.open(args.gtf, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.split("\t", 9)
        if f[2] != "gene":
            continue
        gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1).split(".")[0]
        gname = re.search(r'gene_name "([^"]+)"', f[8]).group(1)
        genes.append((gid, gname, int(f[3]), int(f[4]), f[6]))
genes = pd.DataFrame(genes, columns=["gene_id", "gene_name", "gene_start", "gene_end", "Strand"])
genes = genes.drop_duplicates(["gene_id", "gene_name"])
n_before = len(eqtl)
eqtl = eqtl.merge(genes, on=["gene_id", "gene_name"])
assert len(eqtl) == n_before, "some eQTL genes were not found once in the GTF"
plus = eqtl["Strand"] == "+"
eqtl["dist_to_tss"] = np.where(plus, eqtl["Pos"] - eqtl["gene_start"], eqtl["Pos"] - eqtl["gene_end"])
eqtl["abs_dist_to_tss"] = eqtl["dist_to_tss"].abs()
print(f"eQTLs: {len(eqtl):,} rows; studies: {eqtl['study_label'].nunique()}")

# ---------- 3. which scooby columns stand for which study/tissue ----------
if args.mapping:
    m = pd.read_csv(args.mapping)
    mapping = {(r.study_label, r.tissue_label): r.pred_columns.split(";") for r in m.itertuples()}
else:  # the authors' choice
    onek1k_types = sorted(eqtl.query("study_label == 'OneK1K'")["tissue_label"].unique())
    mapping = {("OneK1K", t): [t.replace(" ", "_")] for t in onek1k_types}
    mapping[("GTEx", "blood")] = [t.replace(" ", "_") for t in onek1k_types]
missing = {c for cols in mapping.values() for c in cols} - set(pred.var_names)
assert not missing, f"mapping names columns not in the predictions: {missing}"


# ---------- 4. metrics ----------
def concordant(x, beta_cut=0.0):
    """Count of non-negligible predictions whose sign matches beta (beta also beyond beta_cut)."""
    big = x["pred_effect"].abs() > args.cutoff
    pos = (x["pred_effect"] > 0) & (x["beta"] > beta_cut)
    neg = (x["pred_effect"] < 0) & (x["beta"] < -beta_cut)
    return int((big & (pos | neg)).sum()), int(big.sum())


rows, merged_all, bins_all = [], [], []
EDGES = [0, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 150000, 262144]
for (study, tissue), cols in mapping.items():
    sub = (eqtl.query("study_label == @study and tissue_label == @tissue")
           .sort_values("sample_size", ascending=False)
           .drop_duplicates(["gene_id", "variant"], keep="first"))
    if sub.empty:
        continue
    pe = pred_df[["variant", "gene_id"]].assign(pred_effect=pred_df[cols].mean(axis=1).values)
    x = sub.merge(pe, on=["variant", "gene_id"])
    if len(x) <= args.min_n:
        continue
    big = x[x["pred_effect"].abs() > args.cutoff]
    c_all, n_big = concordant(x)
    c_strict, _ = concordant(x, beta_cut=args.cutoff)
    rows.append({
        "study": study, "tissue": tissue, "n": len(x),
        "spearman": stats.spearmanr(x["beta"], x["pred_effect"])[0],
        "n_nonneg": n_big,
        "spearman_nonneg": stats.spearmanr(big["beta"], big["pred_effect"])[0] if n_big > 2 else np.nan,
        # paper's scatter-plot quadrants: |pred| > cutoff and |beta| > cutoff, same sign, over all |pred| > cutoff
        "sign_agree_nonneg": c_strict / n_big if n_big else np.nan,
        # paper's distance plot, dashed line: concordant non-negligible predictions over all pairs
        "concordant_over_all": c_all / len(x),
        "spearman_absdist_vs_abseffect": stats.spearmanr(x["abs_dist_to_tss"], x["pred_effect"].abs())[0],
    })
    x = x.assign(study=study, tissue=tissue)
    merged_all.append(x)
    x["distance_bin"] = pd.cut(x["abs_dist_to_tss"], bins=EDGES, include_lowest=True)
    for b, g in x.groupby("distance_bin", observed=False):
        c, nb = concordant(g)
        bins_all.append({"study": study, "tissue": tissue, "distance_bin": str(b), "n": len(g),
                         "n_nonneg": nb, "concordant": c,
                         "binom_p": stats.binomtest(c, nb).pvalue if nb else np.nan})

summary = pd.DataFrame(rows)
summary.to_csv(out / "summary.csv", index=False)
pd.DataFrame(bins_all).to_csv(out / "distance_bins.csv", index=False)
pd.concat(merged_all).to_csv(out / "merged.csv.gz", index=False)

pd.set_option("display.width", 200)
print("\n" + summary.round(3).to_string(index=False))
gtex = summary.query("study == 'GTEx' and tissue == 'blood'")
if len(gtex):
    g = gtex.iloc[0]
    print(f"\nGTEx whole blood: n = {g.n} (paper 908), Spearman = {g.spearman:.3f} (paper 0.449/0.45), "
          f"non-negligible: Spearman = {g.spearman_nonneg:.3f} (paper 0.78), sign agreement = "
          f"{g.sign_agree_nonneg:.1%} (paper 91.6%), concordant/all = {g.concordant_over_all:.3f} (paper 0.228)")
    print("\nGTEx distance bins (paper: 105/113, 26/27, 19/22, 15/16, 22/23, 12/15, 5/7, 2/2, 1/1):")
    d = pd.DataFrame(bins_all).query("study == 'GTEx'")
    print(d[["distance_bin", "n", "concordant", "n_nonneg"]].to_string(index=False))
print(f"\nWritten: {out}/summary.csv, distance_bins.csv, merged.csv.gz")
