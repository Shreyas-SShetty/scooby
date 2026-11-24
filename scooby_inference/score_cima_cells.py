"""Score the Scooby run that used CIMA's own cells: each CIMA eQTL is compared with the Scooby column of its own
CIMA L4 cell type (no mapping to OneK1K types), then set side by side with the earlier run on OneK1K cells.

Run (from scooby_inference/, after  bash run_all_chunks.sh 0 inputs/cima_pairs.tsv step2_config_cima_cells.yaml
results/cima_pairs_cima_cells):  python score_cima_cells.py
  1. checks every row of inputs/cima_pairs.tsv has an output row range in results/cima_pairs_cima_cells/
  2. runs evaluate_eqtl.py with mapping "CIMA type -> its own column" -> eval_out_cima_cells/
  3. per type: new scores next to the earlier OneK1K-cell scores (eval_out_cima/summary.csv, mapped types only)
  4. totals by lineage and by the earlier match group (types with no OneK1K match are group 'new_unmapped')
  5. same pairs, both runs: agreement of the predicted effects, and sign agreement on the pairs both call clear
  6. specificity: for each type, Spearman(beta, own column) vs the same pairs scored with every other column
Writes eval_out_cima_cells/: summary.csv, distance_bins.csv, merged.csv.gz (from evaluate_eqtl.py),
  per_type_vs_onek1k_cells.csv, totals_by_lineage.csv, totals_by_group.csv, same_pairs.csv, specificity.csv
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
RUN = here / "results/cima_pairs_cima_cells"
OUT = here / "eval_out_cima_cells"
OLD = here / "eval_out_cima"          # earlier run: OneK1K cells, CIMA types mapped by cima_mapping.csv
CUT = 0.05
OUT.mkdir(exist_ok=True)

# ---------- 1. complete? ----------
n_rows = sum(1 for _ in open(here / "inputs/cima_pairs.tsv")) - 1
done = set()
for f in glob.glob(f"{RUN}/outs_*.h5ad"):
    s, e = map(int, re.findall(r"outs_(\d+)_(\d+)\.h5ad$", f)[0])
    done.update(range(s, e + 1))
missing = sorted(set(range(n_rows)) - done)
print(f"{RUN.name}: {n_rows:,} input rows, {len(done):,} covered" + (f", MISSING {len(missing)}" if missing else " (complete)"))
if missing:
    sys.exit("run incomplete")

# ---------- 2. own-column mapping, evaluate ----------
eq = pd.read_table(here / "data/cima_eqtl.tsv")
cols = ad.read_h5ad(sorted(glob.glob(f"{RUN}/outs_*.h5ad"))[0], backed="r").var_names
types = sorted(eq["tissue_label"].unique())
assert set(types) <= set(cols), set(types) - set(cols)
own = pd.DataFrame({"study_label": "CIMA", "tissue_label": types, "pred_columns": types})
own.to_csv(OUT / "mapping_own_column.csv", index=False)
subprocess.run([sys.executable, str(here / "evaluate_eqtl.py"), "--pred", f"{RUN}/outs_*.h5ad",
                "--eqtl", str(here / "data/cima_eqtl.tsv"), "--mapping", str(OUT / "mapping_own_column.csv"),
                "--out", str(OUT)], check=True, stdout=subprocess.DEVNULL)
s = pd.read_csv(OUT / "summary.csv").drop(columns="study").rename(columns={"tissue": "cima_type"})
mg = pd.read_csv(OUT / "merged.csv.gz")

# ---------- 3. per type, next to the OneK1K-cell run ----------
lin = pd.read_csv(here / "data/cima_cells/embedded_60donors/placement_lineage.csv")[["cima_type", "cells", "lineage"]]
m = pd.read_csv(here / "cima_mapping.csv")[["tissue_label", "pred_columns", "group"]] \
    .rename(columns={"tissue_label": "cima_type", "pred_columns": "onek1k_column"})
o = pd.read_csv(OLD / "summary.csv").rename(columns={"tissue": "cima_type"})
keep = ["n", "spearman", "n_nonneg", "spearman_nonneg", "sign_agree_nonneg", "concordant_over_all"]
t = (s[["cima_type"] + keep].merge(lin, on="cima_type", how="left").merge(m, on="cima_type", how="left")
     .merge(o[["cima_type"] + keep].add_prefix("old_").rename(columns={"old_cima_type": "cima_type"}),
            on="cima_type", how="left"))
t["group"] = t["group"].fillna("new_unmapped")
t["lineage"] = t["lineage"].fillna("other")
t = t[["lineage", "cima_type", "cells", "group", "onek1k_column"] + keep + ["old_" + k for k in keep]]
t.sort_values(["lineage", "cima_type"]).to_csv(OUT / "per_type_vs_onek1k_cells.csv", index=False)


# ---------- 4. totals ----------
def totals(df, by):
    df = df.assign(n_correct=(df["sign_agree_nonneg"] * df["n_nonneg"]).round(),
                   old_n_correct=(df["old_sign_agree_nonneg"] * df["old_n_nonneg"]).round())
    g = df.groupby(by).agg(cell_types=("cima_type", "size"), n=("n", "sum"), n_clear=("n_nonneg", "sum"),
                           n_correct=("n_correct", "sum"), old_n_clear=("old_n_nonneg", "sum"),
                           old_n_correct=("old_n_correct", "sum"), median_spearman=("spearman", "median"),
                           old_median_spearman=("old_spearman", "median"))
    g.loc["ALL"] = [len(df), df["n"].sum(), df["n_nonneg"].sum(), df["n_correct"].sum(), df["old_n_nonneg"].sum(),
                    df["old_n_correct"].sum(), df["spearman"].median(), df["old_spearman"].median()]
    g["pct_clear"] = 100 * g["n_clear"] / g["n"]
    g["pct_correct_of_clear"] = 100 * g["n_correct"] / g["n_clear"]
    g["old_pct_correct_of_clear"] = 100 * g["old_n_correct"] / g["old_n_clear"]
    return g.round(3)


tl, tg = totals(t, "lineage"), totals(t, "group")
tl.to_csv(OUT / "totals_by_lineage.csv")
tg.to_csv(OUT / "totals_by_group.csv")

# ---------- 5. same pairs in both runs (mapped types only) ----------
om = pd.read_csv(OLD / "merged.csv.gz")[["tissue", "variant", "gene_id", "beta", "pred_effect"]]
j = mg[["tissue", "variant", "gene_id", "beta", "pred_effect"]].merge(
    om, on=["tissue", "variant", "gene_id"], suffixes=("", "_old"))
assert np.allclose(j["beta"], j["beta_old"])
rows = []
for typ, g in j.groupby("tissue"):
    both = (g["pred_effect"].abs() > CUT) & (g["pred_effect_old"].abs() > CUT)
    rows.append({"cima_type": typ, "pairs": len(g), "pearson_new_vs_old": stats.pearsonr(g["pred_effect"], g["pred_effect_old"])[0],
                 "clear_new": int((g["pred_effect"].abs() > CUT).sum()), "clear_old": int((g["pred_effect_old"].abs() > CUT).sum()),
                 "clear_both": int(both.sum()),
                 "same_sign_new_old_on_clear_both": (np.sign(g.loc[both, "pred_effect"]) == np.sign(g.loc[both, "pred_effect_old"])).mean()})
sp_ = pd.DataFrame(rows).merge(lin[["cima_type", "lineage"]], on="cima_type", how="left")
sp_.to_csv(OUT / "same_pairs.csv", index=False)

# ---------- 6. specificity: own column vs other columns, same pairs ----------
paths = sorted(glob.glob(f"{RUN}/outs_*.h5ad"), key=lambda x: int(re.findall(r"_(\d+)_\d+\.h5ad$", x)[0]))
P = ad.concat([ad.read_h5ad(f) for f in paths], merge="same")
pdf = pd.DataFrame(P.X, columns=P.var_names)
pdf["variant"], pdf["gene_id"] = P.obs["variant"].astype(str).values, P.obs["gene_id"].astype(str).values
pdf = pdf.drop_duplicates(["variant", "gene_id"])
lin_of = dict(zip(lin["cima_type"], lin["lineage"]))
rows = []
for typ, g in mg.groupby("tissue"):
    x = g[["variant", "gene_id", "beta"]].merge(pdf, on=["variant", "gene_id"])
    r = {c: stats.spearmanr(x["beta"], x[c])[0] for c in P.var_names}
    others = pd.Series({c: v for c, v in r.items() if c != typ})
    same_lin = others[[lin_of.get(c) == lin_of.get(typ) for c in others.index]]
    other_lin = others[[lin_of.get(c) != lin_of.get(typ) for c in others.index]]
    rows.append({"cima_type": typ, "lineage": lin_of.get(typ), "n": len(x), "spearman_own": r[typ],
                 "rank_of_own_among_73": int((pd.Series(r) > r[typ]).sum()) + 1,
                 "best_column": max(r, key=r.get), "spearman_best": max(r.values()),
                 "median_same_lineage_others": same_lin.median(), "median_other_lineages": other_lin.median()})
spec = pd.DataFrame(rows)
spec.to_csv(OUT / "specificity.csv", index=False)

# ---------- print ----------
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 30)
pd.set_option("display.max_rows", 200)
print("\n=== per CIMA type (new = CIMA cells, own column; old_ = OneK1K cells, mapped column) ===")
print(t[["lineage", "cima_type", "group", "n", "spearman", "old_spearman", "n_nonneg", "old_n_nonneg",
         "sign_agree_nonneg", "old_sign_agree_nonneg"]].sort_values(["lineage", "cima_type"]).round(3).to_string(index=False))
print("\n=== totals by lineage ===\n" + tl.to_string())
print("\n=== totals by earlier match group ===\n" + tg.to_string())
print("\n=== same pairs, both runs ===")
print(sp_.groupby("lineage")[["pairs", "clear_new", "clear_old", "clear_both"]].sum().assign(
    median_pearson=sp_.groupby("lineage")["pearson_new_vs_old"].median()).round(3).to_string())
print(f"all: pairs {sp_['pairs'].sum():,}, median Pearson {sp_['pearson_new_vs_old'].median():.3f}")
print("\n=== specificity (Spearman with beta, same pairs) ===")
print(spec.round(3).sort_values(["lineage", "cima_type"]).to_string(index=False))
print(f"\nown column best of 73 for {(spec['rank_of_own_among_73'] == 1).sum()} of {len(spec)} types; "
      f"median rank {spec['rank_of_own_among_73'].median():.0f}; own > median other-lineage column for "
      f"{(spec['spearman_own'] > spec['median_other_lineages']).sum()} of {len(spec)}")
print(f"\nWritten to {OUT}/")
