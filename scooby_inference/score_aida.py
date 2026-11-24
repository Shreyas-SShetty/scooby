"""Score the Scooby run on AIDA cells: each AIDA eQTL (Table S6, CellSubtype FDR < 0.05) is compared with the Scooby
column of its own AIDA Level-3 cell type, and set next to OneK1K and CIMA scored the same way.

Run (from scooby_inference/, after run_chunks_shared.sh ... inputs/aida_pairs.tsv step2_config_aida.yaml
results/aida_pairs):  python score_aida.py
  1. checks every row of inputs/aida_pairs.tsv has an output row range in results/aida_pairs/
  2. runs evaluate_eqtl.py with mapping "AIDA type -> its own column" -> eval_out_aida/ (paper's metrics)
  3. adds metrics that do not depend on the size of beta (AIDA betas are ~14x smaller than CIMA's):
       sign_any_beta   = clear predictions (|pred| > 0.05) whose sign matches beta, any |beta|, / clear predictions
       sign_big_beta   = same, only rows with |beta| > 0.05
     and computes the same numbers for OneK1K (eval_onek1k/eval_out_ours) and CIMA (eval_out_cima_cells)
  4. totals by lineage; distance bins (sign_any_beta by distance to TSS)
  5. specificity: Spearman(beta, own column) vs the same pairs scored with every other AIDA column
Writes eval_out_aida/: summary.csv, distance_bins.csv, merged.csv.gz (evaluate_eqtl.py), per_type.csv,
  totals_by_lineage.csv, datasets_compared.csv, distance_bins_any_beta.csv, specificity.csv
"""
import glob
import re
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

here = Path(__file__).resolve().parent
RUN = here / "results/aida_pairs"
OUT = here / "eval_out_aida"
CUT = 0.05
OUT.mkdir(exist_ok=True)

# ---------- 1. complete? ----------
n_rows = sum(1 for _ in open(here / "inputs/aida_pairs.tsv")) - 1
done = set()
for f in glob.glob(f"{RUN}/outs_*.h5ad"):
    s, e = map(int, re.findall(r"outs_(\d+)_(\d+)\.h5ad$", f)[0])
    done.update(range(s, e + 1))
missing = sorted(set(range(n_rows)) - done)
print(f"{RUN.name}: {n_rows:,} input rows, {len(done):,} covered" + (f", MISSING {len(missing)}" if missing else " (complete)"))
if missing:
    sys.exit("run incomplete")

# ---------- 2. own column, evaluate ----------
eq = pd.read_table(here / "data/aida_eqtl.tsv")
types = sorted(eq["tissue_label"].unique())
pd.DataFrame({"study_label": "AIDA", "tissue_label": types, "pred_columns": types}).to_csv(OUT / "mapping_own_column.csv", index=False)
subprocess.run([sys.executable, str(here / "evaluate_eqtl.py"), "--pred", f"{RUN}/outs_*.h5ad",
                "--eqtl", str(here / "data/aida_eqtl.tsv"), "--mapping", str(OUT / "mapping_own_column.csv"),
                "--out", str(OUT)], check=True, stdout=subprocess.DEVNULL)
s = pd.read_csv(OUT / "summary.csv").rename(columns={"tissue": "aida_type"})
mg = pd.read_csv(OUT / "merged.csv.gz")
lin = pd.read_csv(here / "aida_mapping.csv").set_index("aida_type")["lineage"].to_dict()
cells = pd.read_csv(here / "data/aida/embedded/placement_per_type.csv").set_index("aida_type")


# ---------- 3. beta-size-free metrics ----------
def counts(x):
    clear = x["pred_effect"].abs() > CUT
    same = np.sign(x["pred_effect"]) == np.sign(x["beta"])
    bigb = x["beta"].abs() > CUT
    return pd.Series({
        "n": len(x), "n_clear": int(clear.sum()),
        "correct_paper": int((clear & same & bigb).sum()),       # paper: |beta| > 0.05 too, denominator all clear
        "correct_any_beta": int((clear & same).sum()),
        "n_clear_big_beta": int((clear & bigb).sum()), "correct_big_beta": int((clear & same & bigb).sum()),
        "share_beta_over_0.05": bigb.mean(), "median_abs_beta": x["beta"].abs().median(),
        "median_abs_pred": x["pred_effect"].abs().median()})


def pct(g):
    g = g.copy()
    g["pct_clear"] = 100 * g["n_clear"] / g["n"]
    g["paper_sign_agree"] = 100 * g["correct_paper"] / g["n_clear"]
    g["sign_any_beta"] = 100 * g["correct_any_beta"] / g["n_clear"]
    g["sign_big_beta"] = 100 * g["correct_big_beta"] / g["n_clear_big_beta"]
    return g


pt = mg.groupby("tissue").apply(counts, include_groups=False)
pt = pct(pt)
pt["lineage"] = pt.index.map(lin)
pt["cells"] = cells.reindex(pt.index)["cells"].values
pt["share_nearest_mapped"] = cells.reindex(pt.index)["share_nearest_mapped"].values
pt = pt.join(s.set_index("aida_type")[["spearman", "spearman_nonneg"]])
pt.index.name = "aida_type"
pt.sort_values(["lineage", "n"], ascending=[True, False]).to_csv(OUT / "per_type.csv")

sumcols = ["n", "n_clear", "correct_paper", "correct_any_beta", "n_clear_big_beta", "correct_big_beta"]
tl = pt.groupby("lineage")[sumcols].sum()
tl.loc["ALL"] = pt[sumcols].sum()
tl = pct(tl)
tl["median_spearman"] = pt.groupby("lineage")["spearman"].median()
tl.loc["ALL", "median_spearman"] = pt["spearman"].median()
tl.to_csv(OUT / "totals_by_lineage.csv")

# same metrics for OneK1K (our run, OneK1K eQTLs, GTEx excluded) and CIMA (CIMA cells, own column)
ds = []
for name, f, q in [("OneK1K (our run)", here / "eval_onek1k/eval_out_ours/merged.csv.gz", "study == 'OneK1K'"),
                   ("CIMA (CIMA cells)", here / "eval_out_cima_cells/merged.csv.gz", None),
                   ("AIDA (AIDA cells)", OUT / "merged.csv.gz", None)]:
    x = pd.read_csv(f)
    if q:
        x = x.query(q)
    c = counts(x)
    c["dataset"], c["cell_types"], c["spearman_pooled"] = name, x["tissue"].nunique(), stats.spearmanr(x["beta"], x["pred_effect"])[0]
    ds.append(c)
ds = pct(pd.DataFrame(ds).set_index("dataset"))
ds.to_csv(OUT / "datasets_compared.csv")

# ---------- 4. distance bins ----------
EDGES = [0, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 150000, 262144]
mg["distance_bin"] = pd.cut(mg["abs_dist_to_tss"], bins=EDGES, include_lowest=True)
db = pct(mg.groupby("distance_bin", observed=False).apply(counts, include_groups=False))
db["binom_p_any_beta"] = [stats.binomtest(int(a), int(b)).pvalue if b else np.nan
                          for a, b in zip(db["correct_any_beta"], db["n_clear"])]
db.to_csv(OUT / "distance_bins_any_beta.csv")

# ---------- 5. specificity ----------
paths = sorted(glob.glob(f"{RUN}/outs_*.h5ad"), key=lambda x: int(re.findall(r"_(\d+)_\d+\.h5ad$", x)[0]))
P = ad.concat([ad.read_h5ad(f) for f in paths], merge="same")
pdf = pd.DataFrame(P.X, columns=P.var_names)
pdf["variant"], pdf["gene_id"] = P.obs["variant"].astype(str).values, P.obs["gene_id"].astype(str).values
pdf = pdf.drop_duplicates(["variant", "gene_id"])
rows = []
for typ, g in mg.groupby("tissue"):
    x = g[["variant", "gene_id", "beta"]].merge(pdf, on=["variant", "gene_id"])
    r = {c: stats.spearmanr(x["beta"], x[c])[0] for c in P.var_names}
    others = pd.Series({c: v for c, v in r.items() if c != typ})
    other_lin = others[[lin.get(c) != lin.get(typ) for c in others.index]]
    cc = np.corrcoef(x[list(P.var_names)].T.values)
    rows.append({"aida_type": typ, "lineage": lin.get(typ), "n": len(x), "spearman_own": r[typ],
                 "rank_of_own_among_35": int((pd.Series(r) > r[typ]).sum()) + 1,
                 "best_column": max(r, key=r.get), "spearman_best": max(r.values()),
                 "median_other_lineages": other_lin.median(),
                 "median_pearson_between_columns": np.median(cc[np.triu_indices_from(cc, 1)])})
spec = pd.DataFrame(rows)
spec.to_csv(OUT / "specificity.csv", index=False)

# ---------- print ----------
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 40)
r1 = lambda d: d.round(1)
print("\n=== per AIDA type (own column) ===")
print(r1(pt[["lineage", "cells", "share_nearest_mapped", "n", "n_clear", "pct_clear", "paper_sign_agree", "sign_any_beta",
             "n_clear_big_beta", "sign_big_beta", "spearman", "median_abs_beta"]]
         .sort_values(["lineage", "n"], ascending=[True, False]).assign(spearman=pt["spearman"].round(3))).to_string())
print("\n=== totals by lineage ===")
print(r1(tl[["n", "n_clear", "pct_clear", "paper_sign_agree", "sign_any_beta", "n_clear_big_beta", "sign_big_beta"]]
         ).assign(median_spearman=tl["median_spearman"].round(3)).to_string())
print("\n=== three datasets, same metrics ===")
print(ds[["cell_types", "n", "median_abs_beta", "share_beta_over_0.05", "median_abs_pred", "n_clear", "pct_clear",
          "paper_sign_agree", "sign_any_beta", "n_clear_big_beta", "sign_big_beta", "spearman_pooled"]].round(3).to_string())
print("\n=== AIDA by distance to TSS ===")
print(r1(db[["n", "n_clear", "pct_clear", "sign_any_beta"]]).assign(binom_p=db["binom_p_any_beta"].map("{:.1e}".format)).to_string())
print("\n=== specificity ===")
print(spec.round(3).sort_values(["lineage", "aida_type"]).to_string(index=False))
print(f"\nown column best of 35 for {(spec['rank_of_own_among_35'] == 1).sum()} of {len(spec)} types; median rank "
      f"{spec['rank_of_own_among_35'].median():.0f}; own > median other-lineage column for "
      f"{(spec['spearman_own'] > spec['median_other_lineages']).sum()} of {len(spec)}; median Pearson between columns "
      f"{spec['median_pearson_between_columns'].median():.3f}")
print(f"\nWritten to {OUT}/")
