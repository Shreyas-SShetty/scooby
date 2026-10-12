"""Step 2e: compare your predictions with the authors' saved predictions for the same rows.

Values are variant effects per cell type: ln(alt+1) - ln(ref+1) of the predicted expression summed over
all cells of the type (natural log, as in the authors' script). Exact equality is not expected (different
GPU, fp16 autocast), but the two runs should agree closely.

Run:  python compare_with_authors.py --ours results/outs_0_199.h5ad
"""
import argparse

import anndata as ad
import numpy as np
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--ours", required=True)
p.add_argument("--theirs", default="data/authors_outs/outs_0_2129.h5ad")
args = p.parse_args()

KEY = ["Chromosome", "Pos", "Ref", "Alt", "gene_name"]


def as_long(path):
    a = ad.read_h5ad(path)
    key = a.obs[KEY].astype(str).agg(":".join, axis=1).to_numpy()
    df = pd.DataFrame(a.X, index=key, columns=a.var_names)
    return df[~df.index.duplicated()]


ours, theirs = as_long(args.ours), as_long(args.theirs)
rows = ours.index.intersection(theirs.index)
cols = ours.columns.intersection(theirs.columns)
print(f"rows: ours {len(ours)}, theirs {len(theirs)}, matched {len(rows)}")
print(f"cell types: ours {len(ours.columns)}, theirs {len(theirs.columns)}, matched {len(cols)}")
if len(rows) < len(ours):
    print("WARNING: some of your rows have no match in the authors' file")
a, b = ours.loc[rows, cols], theirs.loc[rows, cols]

diff = (a - b).abs()
big = b.abs() > 0.05  # the paper's "non-negligible effect" threshold
per_type = pd.DataFrame({
    "pearson_r": [np.corrcoef(a[c], b[c])[0, 1] for c in cols],
    "median_abs_diff": diff.median().values,
    "max_abs_diff": diff.max().values,
    "sign_agree_if_|theirs|>0.05": [(np.sign(a[c][big[c]]) == np.sign(b[c][big[c]])).mean() for c in cols],
}, index=cols)
print("\nPer cell type:")
print(per_type.round(4).to_string())

x, y = a.to_numpy().ravel(), b.to_numpy().ravel()
r = np.corrcoef(x, y)[0, 1]
share_close = (np.abs(x - y) < 0.01).mean()
sign = (np.sign(x[np.abs(y) > 0.05]) == np.sign(y[np.abs(y) > 0.05])).mean()
print(f"\nAll values: Pearson r = {r:.5f} | |diff| < 0.01 for {share_close:.1%} | "
      f"sign agreement (|theirs| > 0.05) = {sign:.1%}")

if r > 0.99 and share_close > 0.95:
    print("PASS: your setup reproduces the authors' predictions.")
else:
    print("CHECK: agreement is lower than expected. Look first at the largest differences:")
    worst = diff.stack().sort_values(ascending=False).head(10)
    print(pd.DataFrame({"ours": a.stack()[worst.index], "theirs": b.stack()[worst.index], "abs_diff": worst}))
