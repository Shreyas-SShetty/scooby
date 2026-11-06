"""Summarise CIMA's lead cis-xQTL table (one top variant per gene/peak per cell type).

Run:  python inspect_lead_variants.py
Writes lead_eqtl_by_celltype.csv (one row per cell type) next to this script.
"""
from pathlib import Path

import numpy as np
import pandas as pd

here = Path(__file__).parent
lead = pd.read_csv(here / "xQTL/CIMA_Lead_cis-xQTL.csv", index_col=0)
lead = lead.rename(columns={"A1(ALT/effect allele)": "alt", "A2(REF)": "ref"})
print(f"{len(lead):,} rows, {lead.shape[1]} columns")
print("columns:", list(lead.columns))
print("\nrows per analysis:\n" + lead["analysis"].value_counts().to_string())

eqtl = lead[lead["analysis"] == "cis-eQTL"].copy()
eqtl["is_snv"] = (eqtl["ref"].str.len() == 1) & (eqtl["alt"].str.len() == 1)
print(f"\ncis-eQTL lead rows: {len(eqtl):,} | cell types: {eqtl['celltype'].nunique()} | "
      f"genes: {eqtl['phenotype_id'].nunique():,} | distinct variants: {eqtl['variant_id'].nunique():,}")
print(f"single-letter variants (SNVs): {eqtl['is_snv'].mean():.1%}; indels: {(~eqtl['is_snv']).sum():,}")
print(f"variant-gene pairs that are lead in more than one cell type: "
      f"{eqtl.duplicated(['variant_id', 'phenotype_id'], keep=False).mean():.1%} of rows")
print(f"distance to TSS (|start_distance|): median {eqtl['start_distance'].abs().median():,.0f} bp, "
      f"within 1 kb {(eqtl['start_distance'].abs() <= 1000).mean():.1%}")
print(f"slope: median |slope| {eqtl['slope'].abs().median():.3f}, positive {(eqtl['slope'] > 0).mean():.1%}")

by_ct = (eqtl.groupby("celltype")
         .agg(n_lead=("phenotype_id", "size"), n_snv=("is_snv", "sum"),
              median_abs_slope=("slope", lambda s: s.abs().median()),
              median_abs_dist=("start_distance", lambda s: s.abs().median()))
         .sort_values("n_lead", ascending=False))
by_ct.to_csv(here / "lead_eqtl_by_celltype.csv")
pd.set_option("display.width", 200)
print("\nlead cis-eQTLs per cell type:\n" + by_ct.round(3).to_string())

cols = ["phenotype_id", "variant_id", "ref", "alt", "af", "slope", "slope_se", "pval_nominal",
        "start_distance", "celltype"]
print("\nstrongest 15 lead cis-eQTLs (by |slope / se|):")
top = eqtl.assign(z=np.abs(eqtl["slope"] / eqtl["slope_se"])).nlargest(15, "z")
print(top[cols].to_string(index=False))
