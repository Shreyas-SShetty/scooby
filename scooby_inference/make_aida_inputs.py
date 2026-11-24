"""Build the AIDA eQTL table and the Scooby input pairs from AIDA Table S6 (Cell 2025, lead eQTLs).

Run (from scooby_inference/):  python make_aida_inputs.py
  1. keeps rows with eQTL_CellSubtype_FDR < 0.05 (the stricter of the two FDRs; 57,875 rows, 11,431 genes)
  2. keeps single-base variants (chr_pos_ref_alt, hg38)
  3. gene symbol -> Ensembl ID: first through AIDA's own gene table (raw/var feature_name of the CELLxGENE h5ad),
     then through a unique GENCODE v32 name. The name Scooby uses is the GENCODE v32 name of that ID (ScoobyVariant
     looks genes up by v32 name in genes_tsscentred.bed and the GTF).
  4. keeps pairs whose variant lies in the gene's Scooby window: the same test ScoobyVariant.py applies
     (bed Start - 163,840 <= Pos < bed End + 163,840), so the GPU does not spend time on rows it would drop
Writes:
  data/aida_eqtl.tsv      one row per kept eQTL (cell type x pair), same columns as data/cima_eqtl.tsv
                          (slope_se = beta / t; af and sample_size not given in Table S6 -> empty)
  inputs/aida_pairs.tsv   distinct variant-gene pairs, ScoobyVariant input format
  data/aida_input_funnel.csv   rows / pairs left after each step
"""
import gzip
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

here = Path(__file__).resolve().parent
S6 = here / "data/aida/cell2025_supplementary/TableS6_AIDA_lead_eQTLs.tsv"
H5 = here / "data/aida/cells/AIDA_DataFreeze_v2.h5ad"
GTF = here / "data/gencode.v32.annotation.sorted.gtf.gz"
BED = here / "scooby_reproducibility/scripts/inference_eqtls/genes_tsscentred.bed"
FDR = 0.05

funnel = []
step = lambda name, df: funnel.append({"step": name, "rows": len(df), "pairs": df[["variant", "Gene"]].drop_duplicates().shape[0],
                                       "genes": df["Gene"].nunique(), "cell_types": df["tissue_label"].nunique()})

d = pd.read_table(S6).rename(columns={"AIDA_CellTypeAnnotation": "tissue_label", "Genetic variant": "variant"})
step("Table S6, all rows (all pass eigenMT FDR < 0.05)", d)
d = d[d["eQTL_CellSubtype_FDR"] < FDR].copy()
step(f"eQTL_CellSubtype_FDR < {FDR}", d)
v = d["variant"].str.split("_", expand=True)
d["Chromosome"], d["pos1"], d["Ref"], d["Alt"] = v[0], v[1].astype(int), v[2], v[3]
d = d[(d["Ref"].str.len() == 1) & (d["Alt"].str.len() == 1)].copy()
step("single-base variants", d)

# ---------- symbol -> Ensembl -> v32 name ----------
v32_name, v32_count = {}, {}
with gzip.open(GTF, "rt") as f:
    for line in f:
        if "\tgene\t" in line:
            gid = re.search(r'gene_id "([^"]+)"', line).group(1).split(".")[0]
            nm = re.search(r'gene_name "([^"]+)"', line).group(1)
            v32_name[gid] = nm
            v32_count[nm] = v32_count.get(nm, []) + [gid]
h = h5py.File(H5, "r")
dec = lambda xs: [x.decode() if isinstance(x, bytes) else x for x in xs]
ids = dec(h["raw/var"][h["raw/var"].attrs["_index"]][:])
fn = h["raw/var/feature_name"]
names = np.array(dec(fn["categories"][:]))[fn["codes"][:]] if isinstance(fn, h5py.Group) else np.array(dec(fn[:]))
aida = pd.DataFrame({"gene_id": ids, "symbol": names})
aida = aida[~aida["symbol"].duplicated(keep=False)].set_index("symbol")["gene_id"].to_dict()
v32_unique = {n: g[0] for n, g in v32_count.items() if len(g) == 1}


def to_id(sym):
    if sym in aida and aida[sym] in v32_name:
        return aida[sym], "AIDA gene table"
    if sym in v32_unique:
        return v32_unique[sym], "GENCODE v32 name"
    return None, "not found"


m = pd.DataFrame([(s, *to_id(s)) for s in d["Gene"].unique()], columns=["Gene", "gene_id", "how"])
print("genes by how matched:\n" + m["how"].value_counts().to_string())
d = d.merge(m, on="Gene")
d = d[d["gene_id"].notna()].copy()
d["gene_name"] = d["gene_id"].map(v32_name)
step("gene matched to a GENCODE v32 ID", d)

# ---------- Scooby window (same test as ScoobyVariant.py) ----------
d["Pos"] = d["pos1"] - 1  # 0-based, as in inputs/cima_pairs.tsv (chr1_814733 -> Pos 814732)
bed = pd.read_table(BED, names=["Chromosome", "Start", "End", "gene_name", "Strand"])
bed["Start"] -= 163840
bed["End"] += 163840
w = d[["Chromosome", "gene_name", "Pos"]].drop_duplicates().merge(bed, on=["Chromosome", "gene_name"])
w = w[(w["Pos"] >= w["Start"]) & (w["Pos"] < w["End"])][["Chromosome", "gene_name", "Pos"]].drop_duplicates()
d = d.merge(w, on=["Chromosome", "gene_name", "Pos"])
step("variant inside the gene's Scooby window", d)

out = pd.DataFrame({"study_label": "AIDA", "tissue_label": d["tissue_label"], "variant": d["variant"],
                    "gene_id": d["gene_id"], "gene_name": d["gene_name"], "Chromosome": d["Chromosome"],
                    "Pos": d["Pos"], "Ref": d["Ref"], "Alt": d["Alt"], "beta": d["eQTL_beta"],
                    "slope_se": d["eQTL_beta"] / d["eQTL_t.stat"], "af": np.nan,
                    "pval_nominal": d["eQTL_p.value"], "study_wise_qval": d["eQTL_CellSubtype_FDR"],
                    "sample_size": np.nan})
out.to_csv(here / "data/aida_eqtl.tsv", sep="\t", index=False)
pairs = (out[["Chromosome", "gene_name", "variant", "Pos", "Ref", "Alt", "gene_id"]].drop_duplicates(["variant", "gene_id"])
         .sort_values(["Chromosome", "Pos", "gene_name"], key=lambda s: s.map(
             lambda c: int(c[3:]) if isinstance(c, str) and c[3:].isdigit() else (99 if isinstance(c, str) else c))))
(here / "inputs").mkdir(exist_ok=True)
pairs.to_csv(here / "inputs/aida_pairs.tsv", sep="\t", index=False)
f = pd.DataFrame(funnel)
f.to_csv(here / "data/aida_input_funnel.csv", index=False)
print(f.to_string(index=False))
print(f"\nwrote data/aida_eqtl.tsv ({len(out):,} rows) and inputs/aida_pairs.tsv ({len(pairs):,} pairs)")
print("rows per cell type:\n" + out["tissue_label"].value_counts().to_string())
