"""Score Scooby's CIMA predictions per CIMA cell type, using the cell-type matches in cima_mapping.csv.

Run (from scooby_inference/, after the GPU run):  python score_cima.py
  1. checks every row of inputs/<run>.tsv has an output chunk in results/<run>/ (run = cima_pairs, or cima_pairs_A/_B)
  2. runs evaluate_eqtl.py with data/cima_eqtl.tsv and cima_mapping.csv -> eval_out_cima/
  3. prints the scores grouped by match strength, next to the OneK1K score for the same Scooby column
     (eval_onek1k/eval_out_ours/summary.csv; note OneK1K uses fine-mapped eQTLs, CIMA lead eQTLs)
  4. prints totals per group, clear predictions by distance to TSS (main groups 1-3), and the unmapped CIMA types
Writes eval_out_cima/cima_scores.csv, group_totals.csv, distance_main.csv.
"""
import glob
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

here = Path(__file__).resolve().parent
OUT = here / "eval_out_cima"

# ---------- 1. is the GPU run complete? ----------
runs = sorted(glob.glob(str(here / "results/cima_pairs*")))
assert runs, "no results/cima_pairs* folder found"
ok = True
for run in runs:
    name = Path(run).name
    n_rows = sum(1 for _ in open(here / "inputs" / f"{name}.tsv")) - 1
    done = set()
    for f in glob.glob(f"{run}/outs_*.h5ad"):
        s, e = map(int, re.findall(r"outs_(\d+)_(\d+)\.h5ad$", f)[0])
        done.update(range(s, e + 1))
    missing = sorted(set(range(n_rows)) - done)
    print(f"{name}: {n_rows:,} input rows, {len(done):,} covered by output files"
          + (f", MISSING rows {missing[0]}-{missing[-1]}" if missing else " (complete)"))
    ok &= not missing
if not ok:
    sys.exit("GPU run incomplete: rerun  bash run_all_chunks.sh <gpu> inputs/<run>.tsv  (done chunks are skipped)")

# ---------- 2. match cell types and score (evaluate_eqtl.py does the join and the metrics) ----------
subprocess.run([sys.executable, str(here / "evaluate_eqtl.py"),
                "--pred", str(here / "results/cima_pairs*/outs_*.h5ad"),
                "--eqtl", str(here / "data/cima_eqtl.tsv"),
                "--mapping", str(here / "cima_mapping.csv"),
                "--out", str(OUT)], check=True)

# ---------- 3. one table: CIMA cell type, Scooby column, match strength, scores, OneK1K reference ----------
m = pd.read_csv(here / "cima_mapping.csv")
s = pd.read_csv(OUT / "summary.csv")
eq = pd.read_table(here / "data/cima_eqtl.tsv")
m["cima_rows"] = m["tissue_label"].map(eq["tissue_label"].value_counts()).fillna(0).astype(int)
t = m.merge(s.drop(columns="study"), left_on="tissue_label", right_on="tissue", how="left").drop(columns="tissue")

ref = here / "eval_onek1k/eval_out_ours/summary.csv"
if ref.exists():
    o = pd.read_csv(ref).query("study == 'OneK1K'")
    o = o.assign(pred_columns=o["tissue"].str.replace(" ", "_"))[
        ["pred_columns", "n", "spearman", "spearman_nonneg", "sign_agree_nonneg"]]
    t = t.merge(o.add_prefix("onek1k_").rename(columns={"onek1k_pred_columns": "pred_columns"}),
                on="pred_columns", how="left")
t.to_csv(OUT / "cima_scores.csv", index=False)

cols = ["tissue_label", "pred_columns", "cima_rows", "n", "spearman", "n_nonneg", "spearman_nonneg",
        "sign_agree_nonneg", "concordant_over_all"] + [c for c in t.columns if c.startswith("onek1k_")]
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 30)
for g, sub in t.groupby("group"):
    print(f"\n=== {g} ===")
    print(sub[cols].round(3).to_string(index=False))
print("\ncima_rows = CIMA lead eQTLs for that cell type; n = those with a Scooby prediction (others were "
      "filtered by ScoobyVariant.py). Blank scores: n <= 10.")

# ---------- 4. totals per group (a pair that is a lead in several cell types counts once per cell type) ----------
t["n_correct"] = (t["sign_agree_nonneg"] * t["n_nonneg"]).round()
tot = t.groupby("group").agg(cell_types=("tissue_label", "size"), n=("n", "sum"),
                             n_nonneg=("n_nonneg", "sum"), n_correct=("n_correct", "sum"))
tot["pct_clear"] = (100 * tot["n_nonneg"] / tot["n"]).round(1)
tot["pct_correct_of_clear"] = (100 * tot["n_correct"] / tot["n_nonneg"]).round(1)
tot.to_csv(OUT / "group_totals.csv")
print("\n=== totals per group (pct_correct_of_clear = sign agreement with |beta| > 0.05, as sign_agree_nonneg) ===")
print(tot.to_string())

# clear predictions by distance to TSS, main groups 1-3 (compare with the GTEx bins in eval_out_ours)
main = t.loc[t["group"] != "4_extra", "tissue_label"]
d = pd.read_csv(OUT / "distance_bins.csv")
d = d[d["tissue"].isin(main)].groupby("distance_bin", sort=False)[["n", "n_nonneg", "concordant"]].sum()
d["pct_clear"] = (100 * d["n_nonneg"] / d["n"]).round(1)
d["pct_concordant_of_clear"] = (100 * d["concordant"] / d["n_nonneg"]).round(1)
d.to_csv(OUT / "distance_main.csv")
print("\n=== main groups 1-3 by distance to TSS (concordant = same sign, any |beta|) ===")
print(d.to_string())

un = eq.loc[~eq["tissue_label"].isin(m["tissue_label"]), "tissue_label"].value_counts()
print(f"\nCIMA cell types not in cima_mapping.csv: {len(un)} ({un.sum():,} of {len(eq):,} eQTL rows)")
print(", ".join(f"{k} ({v})" for k, v in un.items()))
print(f"\nWritten: {OUT}/cima_scores.csv, group_totals.csv, distance_main.csv "
      "(plus summary.csv, distance_bins.csv, merged.csv.gz)")
