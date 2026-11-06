"""Turn CIMA lead cis-eQTLs into (a) a Scooby input list and (b) an eQTL table for evaluate_eqtl.py.

Run:  python prepare_cima.py
Writes ../scooby_inference/inputs/cima_pairs.tsv and ../scooby_inference/data/cima_eqtl.tsv
"""
import gzip
import re
from pathlib import Path

import pandas as pd
from pyfaidx import Fasta

here = Path(__file__).resolve().parent
inf = here.parent / "scooby_inference"

lead = pd.read_csv(here / "xQTL/CIMA_Lead_cis-xQTL.csv", index_col=0)
lead = lead.rename(columns={"A1(ALT/effect allele)": "Alt", "A2(REF)": "Ref", "phenotype_id": "gene_name"})
e = lead[lead["analysis"] == "cis-eQTL"].copy()
n0 = len(e)
e = e[(e["Ref"].str.len() == 1) & (e["Alt"].str.len() == 1)].copy()  # SNVs only, as in the paper
print(f"cis-eQTL rows: {n0:,}; SNVs kept: {len(e):,}")

# position: CIMA chr1_814733 (1-based) -> authors' Pos convention (chr1_841166_A_G has Pos 841165)
cp = e["variant_id"].str.extract(r"^(chr[^_]+)_(\d+)$")
e["Chromosome"], e["pos1"] = cp[0], cp[1].astype(int)
e["Pos"] = e["pos1"] - 1
e["variant"] = e["Chromosome"] + "_" + e["pos1"].astype(str) + "_" + e["Ref"] + "_" + e["Alt"]

# gene names -> Ensembl IDs from the same GENCODE v32 file Scooby uses; keep names that map to one ID
genes = []
with gzip.open(inf / "data/gencode.v32.annotation.sorted.gtf.gz", "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.split("\t", 9)
        if f[2] != "gene":
            continue
        gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1).split(".")[0]
        gname = re.search(r'gene_name "([^"]+)"', f[8]).group(1)
        genes.append((gname, gid))
g = pd.DataFrame(genes, columns=["gene_name", "gene_id"]).drop_duplicates()
g = g.drop_duplicates("gene_name", keep=False)
e = e.merge(g, on="gene_name", how="left")
print(f"rows whose gene maps to one GENCODE v32 ID: {e['gene_id'].notna().mean():.1%}")
e = e.dropna(subset=["gene_id"])

# check REF against the hg38 genome Scooby reads (confirms genome build and position convention)
fa = Fasta(str(inf / "data/genome_human.fa"))
base = [fa[c][p:p + 1].seq.upper() for c, p in zip(e["Chromosome"], e["Pos"])]
e["ref_ok"] = [b == r for b, r in zip(base, e["Ref"])]
alt_hit = sum(b == a for b, a in zip(base, e["Alt"]))
print(f"REF matches hg38: {e['ref_ok'].mean():.2%} | genome shows ALT instead: {alt_hit / len(e):.2%}")
e = e[e["ref_ok"]]

# (a) Scooby input: one row per distinct variant-gene pair, same columns as the authors' input
cols = ["Chromosome", "gene_name", "variant", "Pos", "Ref", "Alt", "gene_id"]
pairs = e[cols].drop_duplicates(["variant", "gene_id"])
pairs.to_csv(inf / "inputs/cima_pairs.tsv", sep="\t", index=False)

# (b) eQTL table in the columns evaluate_eqtl.py reads; beta = CIMA slope (effect of ALT)
eq = e.assign(study_label="CIMA", tissue_label=e["celltype"], beta=e["slope"], sample_size=406)
eq[["study_label", "tissue_label", "variant", "gene_id", "gene_name", "Chromosome", "Pos", "Ref", "Alt",
    "beta", "slope_se", "af", "pval_nominal", "study_wise_qval", "sample_size"]].to_csv(
    inf / "data/cima_eqtl.tsv", sep="\t", index=False)
print(f"wrote {len(pairs):,} pairs to inputs/cima_pairs.tsv and {len(eq):,} rows to data/cima_eqtl.tsv")
