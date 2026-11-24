"""Map each CIMA cell type in celltypes.tsv to the closest of Scooby's 22 OneK1K columns, using the Cell Ontology.

Run (from scooby_inference/):
  curl -L -o cl-basic.obo http://purl.obolibrary.org/obo/cl/cl-basic.obo      # Cell Ontology, ~4 MB
  python map_cima_celltypes.py --obo cl-basic.obo
Writes cima_mapping_ontology.csv (same first 3 columns as cima_mapping.csv, so it can be passed to
evaluate_eqtl.py --mapping), plus how each match was found and whether it agrees with cima_mapping.csv.

How a CIMA cl_id is matched (is_a links only; the first rule that applies wins, then the shortest path):
  1. exact        the cl_id is one of a column's IDs
  2. subtype_of   a column's ID is an ancestor of the cl_id (the CIMA type is a kind of that column's cells)
  3. broader_than the cl_id is an ancestor of a column's ID (CIMA type is broader; often several columns tie)
  4. related_via  only a shared ancestor (reported in 'via'); a weak match, check by hand
Ties (same rule, same distance) give several columns separated by ';' (evaluate_eqtl.py averages them).
Rule-4 picks that do not fit the biology are replaced by OVERRIDES below (relation 'manual'; the ontology's pick
is kept in ontology_pred_columns). Types with no comparable column are left out and listed at the end.
"""
import argparse
from collections import deque

import pandas as pd

# Scooby column -> Cell Ontology IDs of the cells it contains.
# The first ID is the column's own ID: OneK1K's tissue_ontology_id in the eQTL Catalogue table
# (eQTL-Catalogue-resources/ontology_mappings/tissue_ontology_mapping.tsv), which the authors used to build the
# columns (scooby_reproducibility/notebooks/data_preparation/onek1k1/analyze_onek1k_anndata_new.ipynb, cells 16-23).
# Extra "effective" IDs are added where that ID is broader than the cells in the column, so that, for example,
# CIMA's Th17 cells do not become a "subtype" of the naive-only CD4+_T_cell column:
ONEK1K = {
    "CD4+_T_cell": ["CL:0000895"],                  # ID CL:0000624 (all CD4 ab T), but only CD4_Naive cells
    "CD8+_T_cell": ["CL:0000900"],                  # ID CL:0000625 (all CD8 ab T), but only CD8_Naive cells
    "B_cell": ["CL:0000788", "CL:0000818"],         # ID CL:0000236 (all B): B_naive + B_intermediate (+ transitional)
    "NK_cell": ["CL:0000939", "CL:4033071"],        # ID CL:0000623 (all NK): CD56-dim NK + NK_Proliferating
    "dendritic_cell": ["CL:0000990"],               # ID CL:0000451 (all DC): cDC2 (conventional DC); pDC is separate
    "monocyte": ["CL:0002057", "CL:0000860"],       # CD14_Mono = classical monocyte
    "CD16+_monocyte": ["CL:0002396", "CL:0000875"],  # CD16_Mono = non-classical monocyte
    "Treg_memory": ["CL:0002677", "CL:0000815"],    # Treg (cellxgene CL:0000815 was renamed to CL:0002677)
    "memory_B_cell": ["CL:0000787"],
    "CD4+_TCM_cell": ["CL:0000904"],
    "CD4+_TEM_cell": ["CL:0000905"],
    "CD4+_CTL_cell": ["CL:0000934"],
    "CD8+_TCM_cell": ["CL:0000907"],
    "CD8+_TEM_cell": ["CL:0000913"],
    "MAIT_cell": ["CL:0000940"],
    "gdT_cell": ["CL:0000798"],
    "dnT_cell": ["CL:0002489"],
    "CD56+_NK_cell": ["CL:0000938"],
    "plasmablast": ["CL:0000980"],
    "plasmacytoid_dendritic_cell": ["CL:0000784"],
    "hematopoietic_precursor_cell": ["CL:0008001"],
    "platelet": ["CL:0000233"],
}
# The columns' own IDs where they differ from the effective IDs above. They count for rule 1 (exact) only:
# e.g. CIMA DC_CSF2RA has CL:0000451 "dendritic cell", the ID OneK1K itself gave the dendritic_cell column.
OWN_ID = {"CD4+_T_cell": "CL:0000624", "CD8+_T_cell": "CL:0000625", "B_cell": "CL:0000236",
          "NK_cell": "CL:0000623", "dendritic_cell": "CL:0000451"}

# Manual decisions for CIMA types the ontology can only link through a shared ancestor (rule 4), where its pick
# does not fit the biology. "" = leave unmapped (no comparable OneK1K column; the row is left out of the output).
OVERRIDES = {
    **{t: ("CD4+_TEM_cell", "CD4 helper subsets are memory/effector cells; the ontology linked them to naive CD4")
       for t in ["CD4_Th_CCR4", "CD4_Th_CR1", "CD4_Th_LMNA", "CD4_Th_TNFRSF11A", "CD4_Th1-like_GZMK",
                 "CD4_Th17-like_RORC", "CD4_Th22-like_CCR10", "CD4_Tfh-like_CXCR5"]},
    **{t: ("CD8+_TEM_cell", "GZMB/GZMK cytotoxic CD8 cells are effector memory; the ontology linked them to naive CD8")
       for t in ["CD8_CTL_GZMB", "CD8_CTL_GZMK", "CD8_CTL_IFI44L"]},
    "Proliferative_T_MKI67": ("", "cycling T cells of mixed lineage; no matching OneK1K column"),
    "ILC2_IL2RA": ("", "innate lymphoid cells were removed from the OneK1K training data"),
    "NKT_NCR1": ("", "NKT cells have no OneK1K column (MAIT is a different T-cell type)"),
    "NKT_IFNG": ("", "NKT cells have no OneK1K column (MAIT is a different T-cell type)"),
}
RANK = {"exact": 0, "subtype_of": 1, "broader_than": 2, "related_via": 3, "manual": 4}

p = argparse.ArgumentParser()
p.add_argument("--obo", required=True, help="Cell Ontology file (cl-basic.obo)")
p.add_argument("--celltypes", default="celltypes.tsv")
p.add_argument("--current", default="cima_mapping.csv", help="existing mapping, for comparison")
p.add_argument("--out", default="cima_mapping_ontology.csv")
args = p.parse_args()

# ---------- read the ontology: name and is_a parents of every term ----------
names, parents = {}, {}
for stanza in open(args.obo).read().split("\n\n"):
    lines = stanza.strip().split("\n")
    if lines[0] != "[Term]" or "is_obsolete: true" in lines:
        continue
    tags = [line.split(": ", 1) for line in lines[1:] if ": " in line]
    tid = next(v for k, v in tags if k == "id")
    names[tid] = next((v for k, v in tags if k == "name"), "")
    parents[tid] = [v.split(" ")[0] for k, v in tags if k == "is_a"]


def ancestors(x):
    """{ancestor id: number of is_a steps up}, including x itself at 0."""
    dist, q = {x: 0}, deque([x])
    while q:
        c = q.popleft()
        for par in parents.get(c, []):
            if par not in dist:
                dist[par] = dist[c] + 1
                q.append(par)
    return dist


missing = [i for ids in ONEK1K.values() for i in ids if i not in names]
assert not missing, f"OneK1K IDs not in the ontology file: {missing}"
anc = {i: ancestors(i) for ids in ONEK1K.values() for i in ids}


def best_match(x):
    ax = ancestors(x)
    found = []  # (rank, distance, column, via)
    for col, ids in ONEK1K.items():
        options = []
        for m in ids:
            am = anc[m]
            if x == m or x == OWN_ID.get(col):
                options.append((0, 0, ""))
            elif m in ax:
                options.append((1, ax[m], m))
            elif x in am:
                options.append((2, am[x], m))
            else:
                common = set(ax) & set(am)
                if common:
                    c = min(common, key=lambda c: (ax[c] + am[c], c))
                    options.append((3, ax[c] + am[c], c))
        if options:
            r, d, via = min(options)
            found.append((r, d, col, via))
    r, d = min((f[0], f[1]) for f in found)
    cols = sorted(f[2] for f in found if (f[0], f[1]) == (r, d))
    vias = sorted({f[3] for f in found if (f[0], f[1]) == (r, d) and f[3]})
    return cols, [k for k, v in RANK.items() if v == r][0], d, vias


ct = pd.read_table(args.celltypes).query("scheme == 'cima'")
rows = []
for r in ct.itertuples():
    if not isinstance(r.cl_id, str) or r.cl_id not in names:
        rows.append({"tissue_label": r.label, "cl_id": r.cl_id, "relation": "no_cl_id"})
        continue
    cols, rel, d, vias = best_match(r.cl_id)
    rows.append({"tissue_label": r.label, "cl_id": r.cl_id, "cl_name": names[r.cl_id], "pred_columns": ";".join(cols),
                 "relation": rel, "distance": d, "n_columns": len(cols),
                 "via": "; ".join(f"{v} {names[v]}" for v in vias)})
out = pd.DataFrame(rows)
out.insert(0, "study_label", "CIMA")

# manual decisions: keep the ontology's pick for reference
out["ontology_pred_columns"], out["note"] = out["pred_columns"], ""
for t, (cols, why) in OVERRIDES.items():
    i = out["tissue_label"] == t
    out.loc[i, ["pred_columns", "relation", "note"]] = [cols, "manual", why]
unmapped = out[(out["pred_columns"].fillna("") == "")]
out = out.drop(unmapped.index)

cur = pd.read_csv(args.current)[["tissue_label", "pred_columns", "group"]].rename(
    columns={"pred_columns": "current_pred_columns", "group": "current_group"})
out = out.merge(cur, on="tissue_label", how="left")
out["agrees"] = (out["pred_columns"] == out["current_pred_columns"]).astype(object)
out.loc[out["current_pred_columns"].isna(), "agrees"] = pd.NA
out = out.sort_values(["relation", "distance", "tissue_label"], key=lambda s: s.map(RANK) if s.name == "relation" else s)
out.to_csv(args.out, index=False)

pd.set_option("display.width", 250)
pd.set_option("display.max_colwidth", 60)
print(out[["tissue_label", "cl_name", "pred_columns", "relation", "distance", "via", "current_pred_columns", "agrees"]]
      .to_string(index=False))
print(f"\n{len(out)} CIMA types mapped: " + ", ".join(f"{k} {v}" for k, v in out["relation"].value_counts().items()))
print(f"left unmapped ({len(unmapped)}): " + ", ".join(f"{r.tissue_label} ({r.note})" for r in unmapped.itertuples()))
both = out["agrees"].notna()
print(f"in cima_mapping.csv too: {both.sum()}; same column(s): {int(out.loc[both, 'agrees'].sum())}")
print(f"Written: {args.out}")
