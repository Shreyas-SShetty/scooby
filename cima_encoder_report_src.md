---
title: "Using the OneK1K cell encoder for CIMA cells"
subtitle: "Every step and check, and the CIMA eQTL results with CIMA's own cells"
date: "23 November 2026 (update of the 22 November version)"
geometry: margin=1.9cm
fontsize: 10pt
monofont: "JetBrains Mono"
colorlinks: true
toc: true
toc-depth: 2
header-includes:
  - \usepackage{booktabs}
  - \usepackage{float}
  - \renewcommand{\arraystretch}{1.2}
  - \floatplacement{figure}{H}
---

\newpage

# Summary

**Goal.** Run the released OneK1K Scooby model (`lauradmartens/onek1k-scooby`) with **CIMA's own cells** instead of
OneK1K's, so that CIMA eQTLs are predicted in CIMA cell types by cells from the CIMA cohort.

**Problem.** Scooby does not take a cell's gene expression directly. Each cell enters as **10 numbers** (its
embedding), which the authors computed with a trained scPoli encoder. Scooby was fine-tuned on those 10 numbers.
New cells must therefore be encoded by the same encoder, and **the authors did not release it**.

**What we did.** We trained a substitute encoder that reproduces the authors' 10 numbers from each OneK1K cell's
gene counts, tested it on OneK1K pools it never saw, tested whether Scooby's predictions change when the
authors' embeddings are replaced by ours, and then applied it to CIMA cells.

**Results so far.**

| Check | Result |
|:--------------------------------------------------------------|:---------------------------------------------------------------|
| A. Substitute reproduces the authors' 10 numbers on unseen OneK1K pools | R² 0.990–1.000 per dimension (mean 0.997); median error 0.014 vs 1.018 between cell-type centres |
| B. Scooby predictions with substitute vs authors' embeddings (1,585 pairs × 22 columns) | Pearson r = 0.999998; largest difference 0.022; every eQTL score unchanged to 3 decimals |
| C. CIMA B cells (566,569) land on the right OneK1K columns | naive B → B cell 90–93%, plasma → plasmablast 95–99.6%, memory B 45–94% (OneK1K's own memory B: 65%) |
| D. Sequencing depth (CIMA 31% deeper) | thinning CIMA to OneK1K depth moves cells by 0.020 (median); 93% keep the same nearest cell type |
| E. Partial download of the full CIMA file | 2-donor test: 4,075 cells, 6,571,772 values, **0 differences** from the full download |
| F. The 60-donor file is exact | 928,536 cells, 60 donors, 73 types; its 70,394 B-lineage cells equal the B-cell file: 101,375,855 values, **0 differences** |
| G. All 73 CIMA types land in the right place | right **lineage** 96.7% (OneK1K's own cells: 98.2%); right **column** 49.4% (OneK1K's own cells: 63.7%); T-cell subtypes blur (Section 12) |
| H. Scooby runs with CIMA cells (test, 50 pairs) | 73 columns, all values finite; effects vs the earlier run with OneK1K cells: r = 0.996 |
| I. Full run: CIMA eQTLs scored in their own CIMA type | all **69** CIMA eQTL types scored (before: 65); clear predictions with the right sign **87.7%** (2,434 / 2,775) vs **87.2%** (2,463 / 2,825) with OneK1K cells on the same 53,157 eQTLs; predicted effects r = 0.9985 |
| J. Cell-type specificity of the predictions | a type's own column predicts its eQTLs best for only 1 of 69 types (median rank 29 of 73); the 73 columns correlate at median r = 0.958 |
| K. Does Scooby capture measured cell-type differences? | for 8,130 eQTLs measured in several CIMA types, the predicted difference between two types has the right direction in 49.5% of 114,494 cases (chance); typical predicted difference 0.0003 vs measured 0.105. Where Scooby predicts a difference > 0.05 (62 eQTLs) it is right more often than chance (41 of 62, p = 0.015). Same on OneK1K fine-mapped eQTLs with the authors' embeddings (265 eQTLs: 52.2%, not significant per eQTL). Scooby does know the cell types (e.g. LYZ predicted 160x higher per cell in classical monocytes than in naive B cells), but was never trained on variants |
| L. AIDA cells (1,178,372, 619 donors) land in the right place | 2,228 of 2,403 encoder genes found; right **lineage** 97.3%, right **column** 73.3% (Section 19.5) |
| M. AIDA depth (2.5× OneK1K) | thinning to OneK1K depth: right lineage 97.3% → 96.8%; cells move by 0.059 (median); 79.9% keep the same nearest column (Section 19.6) |
| N. AIDA eQTLs scored in their own AIDA type (33,751 rows, 20 types) | clear predictions 3.4%; right sign **84.6%** (957 / 1,131) at any \|beta\|, vs OneK1K 89.2%, CIMA 87.7%; the paper's metric gives only 38.4% because 75% of AIDA betas are below 0.05; own column best for 1 of 20 types, columns correlate at median r = 0.969 (Section 19.8) |

**New since 20 November, 20:00.** The full run finished (07:14 IST, 21 November) and was scored (Section 14). What
changed from the previous CIMA results, and what did not, is listed in Section 15. Section 16 asks whether a missing
clear prediction means a small true effect or a model limit. Section 17 tests whether Scooby captures the measured
differences of an eQTL's effect between cell types, and why it mostly does not.

**New since 21 November.** Section 18 compares how the eQTL inputs were prepared for OneK1K, CIMA and AIDA, how many eQTLs each gives (OneK1K 40–50 times fewer: fine-mapped only), and their effect-size scales (AIDA's betas about 14 times smaller, which matters for the paper's sign metric).
Section 19 covers AIDA so far: which eQTLs (the two FDR columns), the input funnel (27,857 pairs), encoding
AIDA's cells, and the placement and depth checks (L, M).

**New since 22 November.** The AIDA Scooby run finished (17:42 IST, 22 November; all 27,857 pairs) and was scored
(Section 19.8, check N). Section 1.1 has a one-page schematic of how Scooby works in training and in evaluation.

**Main finding.** Using CIMA's own cells removes the cell-type mapping and adds 4 CIMA types that had no OneK1K
counterpart, but changes the eQTL scores very little (+0.5 percentage points). Scooby's *predicted* effects are
nearly the same in every cell type (check J), so which cells fill a column matters little for this evaluation. This
is a property of the model, not a statement about biology: where CIMA measures an eQTL with different strength in
two cell types, Scooby's predictions mostly miss the difference (check K, Section 17). Scooby learned how much of
each gene each cell type makes, but it was trained on the reference genome only and never saw a variant.

\newpage

# 1. Why an encoder is needed at all

## 1.1 What Scooby takes as input

From the OneK1K training script (`scooby_reproducibility/scripts/training_onek1k/train_rna_only.py`) and the
scooby package (`scooby/data/scdata.py`):

| Input per training step | Source |
|:--------------------------------------|:--------------------------------------------------------------------------------|
| DNA sequence (524 kb) | reference genome `genome_human.fa` (no genotypes, no variants) |
| 64 random cells, **each as 10 numbers** | `embedding_no_val_genes.pq` (`scdata.py` lines 528–534) |
| target: those cells' RNA coverage | `snapatac_merged_plus/minus.h5ad` |

Scooby never sees cell-type labels. It learns "a cell with these 10 numbers has this coverage along this DNA".
Cell types enter only (1) when the embedding was made, and (2) at our inference, where `celltype_fixed.pq` lists
which cells to sum into each of the 22 columns.

Figure 1 shows the whole model on one page. The DNA path (Borzoi) is the same for every cell; each cell enters
only through `cell_state_to_conv`, a small network (10 → 128 → 256 → 1,024 → 3,842 numbers; `scooby/modeling/scooby.py`
lines 41–51) that turns the cell's 10 numbers into the weights and bias of that cell's last layer. The two paths
meet only there.

![How Scooby works: DNA goes through Borzoi, each cell embedding goes through the MLP, and the two meet in the last layer. Below: inputs and outputs in training and evaluation.](fig0_schematic.pdf){width=90%}

## 1.2 Where the 10 numbers came from

The authors' notebook `notebooks/data_preparation/onek1k1/analyze_onek1k_anndata_new.ipynb`:

| Cell | What it does |
|:----|:-------------------------------------------------------------------------------------------------|
| 76 | loads `OneK1K_only_immune.h5ad` (1,267,768 cells) |
| 77–79 | keeps highly variable genes minus validation/test genes: **2,403 genes** |
| 82 | `scPoli(adata_train, condition_keys='sample', cell_type_keys='cell_label', recon_loss='nb')`; prints encoder input 2,403, latent **10**, 75 conditions |
| 83 | trains it (50 epochs, 40 pre-training) |
| 85 | saves `onek1k_training_data/scpoli_model_default.pkl` |
| 86 | **loads a different file**, `onek1k_training_data/scpoli_model.pkl` |
| 87–94 | `get_latent(adata_train, mean=True)`, saved as `embedding_no_val_genes.pq` |

## 1.3 Why the 10 numbers must come from the same encoder

Each number is a weighted combination of the cell's gene counts, with weights chosen during training. Example
(OneK1K monocyte `AAACGGGTCTTGTTTG-1`): authors' numbers `−0.07, −0.11, −0.28, 0.04, 0.02, −0.09, 1.68, −0.19,
0.11, −0.70`. Dimension 7 is high (1.68) for monocytes, −0.20 for a CD4 T cell. A newly trained encoder, even on
the same data, would arrange cells differently (for example monocytes high on another dimension), and Scooby,
which learned to read the authors' arrangement, would misread them.

## 1.4 The authors' encoder is not available

* **Not released.** The listing of Zenodo record 15517072 (`training_data.zip` → `onek1k_training_data/`)
  contains `embedding_no_val_genes.pq`, `celltype_fixed.pq`, `OneK1K_only_immune.h5ad.gz`, gene lists and
  pseudobulk tracks, but **no `scpoli_model.pkl` or `scpoli_model_default.pkl`**.
* **Cannot be retrained identically.** No random seed is set in the notebook (`seed`, `manual_seed`,
  `random_state`: no matches); the embedding came from `scpoli_model.pkl`, whose training is not shown; cells
  74–86 were run out of order (execution counts 8, 2, 17, 18, 19, 20, 21, 48, 44, 49, 21, 50, 54, including a
  kernel restart); GPU training is usually not bit-for-bit reproducible.

**Approach chosen.** For the same 1,267,768 OneK1K cells we have both the encoder's **input** (raw counts) and
its **output** (the 10 numbers). A network trained to map one to the other learns the authors' arrangement
directly, without needing their weights or seed.

\newpage

# 2. Data obtained and checked

| File | Source | Size | Checked |
|:----------------------------------------------|:--------------------|:--------------------|:----------------------------------------------------|
| `OneK1K_only_immune.h5ad` | Zenodo 15517072, `training_data.zip` | 2,514,429,404 B (gz), 14.1 GB unzipped | 1,267,768 cells × 32,738 genes; `X` = raw counts (e.g. 1, 1, 1, 3, 1, 18 …; all integers); Ensembl IDs + `GeneSymbol`; `sample` = 75 pools |
| val/test gene lists (`*_gene_sequences.csv`) | same zip | 147,724 B, 130,891 B | tab-separated: chrom, start, end, gene ID, strand; 3,539 and 3,129 genes |
| `embedding_no_val_genes.pq` | same zip (in the repo) | 61 MB | 1,267,768 rows of 10 numbers |
| CIMA B-cell file (`CIMA_RNA_B_…h5ad`) | CNGB, CIMA resource page | 4,327,375,006 B | size equals the server's; read through to the last cell (566,569 cells; 797,445,229 nonzero counts) |

**Input checks before training:**

1. **Genes rebuild exactly.** 3,000 highly variable genes − 324 in the validation list − 273 in the test list =
   **2,403**, the input size scPoli printed.
2. **Same cell order.** The 1,267,768 cell names in `embedding_no_val_genes.pq` equal those in
   `OneK1K_only_immune.h5ad`, in the same order, and are unique. Row *i* of the input and row *i* of the target
   are the same cell.

\newpage

# 3. Training the substitute encoder

Script: `scooby_inference/fit_scpoli_encoder.py`. Output: `scooby_inference/scpoli_substitute/`.
Runtime: 2 minutes on one GPU.

## 3.1 Example cell

`AAACGGGTCTTGTTTG-1` (row 77, pool `pool1`, label `Mono C`), 7,910 counts in total; 183 of the 2,403 genes are
nonzero.

| Gene | Raw count | Input value log(1 + count) |
|:---------------------|---------:|-----------------------------:|
| TMSB4X | 187 | 5.236 |
| FTL | 124 | 4.828 |
| RPLP1 | 72 | 4.290 |
| RPS2 | 70 | 4.263 |
| MT-CO1 | 54 | 4.007 |
| LYZ (monocyte marker) | 52 | 3.970 |
| other 2,220 genes | 0 | 0 |

Target: the authors' `−0.07, −0.11, −0.28, 0.04, 0.02, −0.09, 1.68, −0.19, 0.11, −0.70`.

## 3.2 Settings

| Setting | Value | Why |
|:--------------------|:----------------------------------------------|:-----------------------------------------------|
| input | log(1 + raw count), 2,403 genes, no depth normalisation | what we believe scPoli's encoder uses with `recon_loss='nb'` (not checked in scArches code; supported by R² = 0.997 below) |
| target | authors' 10 numbers | — |
| network | 2,403 → 512 → 256 → 10 (linear, LayerNorm, ReLU) | larger than scPoli's 2,403 → 50 → 10, to imitate it closely |
| loss | mean squared difference from the authors' numbers | — |
| optimiser | Adam, learning rate 0.001, cosine decay, batch 4,096, 40 epochs | — |
| epoch choice | best of 40 on 5% of training cells held aside | avoids memorising |
| pool input | **not used** | CIMA cells belong to none of OneK1K's 75 pools |

## 3.3 Honest test: cross-validation by pool

The 75 pools were shuffled (seed 0) and split into 5 groups of 15 pools (244,556–260,802 cells per group). Five
networks were trained, each without one group, and each predicted the group it did not see. Every cell thus has
a prediction from a network that never saw its pool, the closest test on OneK1K to a new cohort. A sixth network,
trained on all pools, is the one used for CIMA (`encoder.pt`).

\newpage

# 4. Check A: the substitute reproduces the authors' numbers (OneK1K, unseen pools)

From `scpoli_substitute/heldout_summary.json` and `heldout_per_column.csv`.

| Measure | Value |
|:------------------------------------------------------------|:------------------------------------------------|
| R² per dimension (1–10) | 0.9995, 0.9987, 0.9923, 0.9994, 0.9902, 0.9977, 0.9998, 0.9966, 0.9995, 0.9999 |
| mean R² | **0.9974** |
| median distance between our 10 numbers and the authors' | **0.0138** |
| median distance between two cell-type centres | 1.0179 |
| cells nearest the same cell-type centre as with the authors' numbers | **96.1%** |

![Held-out predictions against the authors' numbers. Left: the best dimension (7, which separates monocytes). Middle: the worst dimension (5, which has a small range, −0.03 to 0.10). Right: R² of all 10 dimensions.](fig1_heldout.png){width=100%}

**Per column.** Same nearest centre as with the authors' numbers: 92.7% (Treg memory) to 99.95% (monocyte).
The cells that differ have the same small error (median distance 0.007–0.026 per column); they sit almost exactly
between two similar types (mostly neighbouring T-cell types: CD4+ TCM 93.8%, CD8+ TCM 93.4%, Treg 92.7%), so a
shift of 0.01 changes which centre is nearer.

| Column | Cells | Median distance to authors' | Nearest own centre: authors' / ours | Same nearest centre |
|:-------------------------|-------:|------:|:----------------:|------:|
| monocyte | 35,359 | 0.009 | 99.5% / 99.5% | 99.95% |
| plasmablast | 3,627 | 0.017 | 99.8% / 99.7% | 99.9% |
| CD16+ monocyte | 15,709 | 0.012 | 96.4% / 96.3% | 99.8% |
| NK cell | 164,544 | 0.013 | 88.5% / 88.5% | 99.4% |
| B cell | 95,108 | 0.007 | 83.4% / 83.3% | 98.3% |
| memory B cell | 30,069 | 0.008 | 65.1% / 65.8% | 96.5% |
| CD4+ T cell | 258,664 | 0.017 | 77.7% / 77.9% | 95.9% |
| CD4+ TCM cell | 287,899 | 0.015 | 42.3% / 42.3% | 93.8% |
| CD8+ TCM cell | 16,333 | 0.014 | 23.9% / 23.8% | 93.4% |
| Treg memory | 26,214 | 0.019 | 32.7% / 32.6% | 92.7% |
| platelet | 516 | 0.020 | 54.3% / 54.3% | 96.9% |

(The full 22-column table is in `heldout_per_column.csv`.) The "nearest own centre" share is low for several T-cell
types **even with the authors' numbers** (CD8+ TCM 23.9%): their T-cell subtypes overlap in this space. What
matters is that ours matches theirs (23.8% vs 23.9%).

**Conclusions from check A.** (1) The authors' embedding is almost entirely a function of the cell's own genes:
leaving out scPoli's pool input costs almost nothing. (2) log(1 + raw count) is very likely the right input; with a
wrong transform, R² could not reach 0.997.

\newpage

# 5. Check B: Scooby gives the same predictions with the substitute embeddings

The decisive test: does replacing the authors' embeddings by ours change what Scooby predicts?

**Set-up.** `step2_config_substitute.yaml` is identical to `step2_config.yaml` except `embedding_path` points to
`scpoli_substitute/embedding_heldout_predicted.pq` (every OneK1K cell encoded by the network that did not see its
pool). Same model, same 22 columns, same 1,585 GTEx/OneK1K pairs. Any difference is due to the embeddings
alone. `run_all_chunks.sh` was extended to take a config and an output folder (3rd and 4th arguments), so the
original results were not touched. Run: 31 min 52 s on GPU 0, finished 03:53 IST.

**Raw predictions** (`eval_onek1k/eval_out_substitute_vs_original/`):

| Measure | Value |
|:--------------------------------------------------------|:-------------------------|
| values compared (1,585 pairs × 22 columns) | 34,870 |
| Pearson r | **0.999998** |
| identical to 7 decimals | 2,207 (6.3%) |
| differ by less than 0.0001 | 95.1% |
| differ by more than 0.01 | 6 values |
| largest difference | 0.022 |
| same clear/not-clear call (\|effect\| > 0.05) | 99.99% (3 values cross the cutoff) |
| same sign where the original is clear | 100% |

Largest differences, all in small columns:

| Pair | Column (cells) | Authors' embeddings | Substitute |
|:------------------------------------|:---------------------|------:|------:|
| chr9_37002145_A_G → ENSG00000250850 | platelet (516) | 0.0974 | 0.1192 |
| chr6_31615948_C_T → ENSG00000204472 | dnT cell (1,817) | −0.5775 | −0.5607 |
| chr11_1852991_C_T → ENSG00000130592 | platelet (516) | −0.0070 | −0.0182 |
| chr6_31615948_C_T → ENSG00000204472 | CD4+ CTL cell (17,939) | −0.5723 | −0.5672 |

![Scooby's variant effects with the authors' embeddings vs the substitute. Right: distribution of the absolute differences (log scale; the 2,207 exact zeros are not shown); the red line is the 0.05 cutoff used for "clear" predictions.](fig2_scooby.png){width=92%}

**eQTL scores** (`eval_onek1k/eval_out_substitute/` vs `eval_onek1k/eval_out_ours/`):

| Score | Authors' embeddings | Substitute |
|:----------------------------------------------------------|:---------------:|:--------------:|
| GTEx whole blood, Spearman (n = 908) | 0.450 | 0.450 |
| GTEx, clear pairs: Spearman / sign agreement | 0.783 / 91.5% (224) | 0.782 / 91.6% (225) |
| GTEx, concordant / all | 0.226 | 0.227 |
| OneK1K CD4+ TCM, Spearman (n = 391) | 0.403 | 0.403 |
| OneK1K NK cell, Spearman (n = 230) | 0.416 | 0.416 |
| OneK1K monocyte, Spearman (n = 26) | 0.685 | 0.685 |
| OneK1K Treg memory, Spearman (n = 21) | 0.094 | 0.092 |
| precision@1, 155 pairs with a clear effect | 0.290 | 0.290 |
| same top-1 column, those 155 pairs | — | 100% |

All 11 OneK1K cell types agree to three decimals except Treg memory (0.094 vs 0.092).

**Problem met and fixed.** The first step-5 run crashed in `enformer_pytorch/data.py` line 58
(`np.fromstring … The binary mode of fromstring is removed`). NumPy in `~/venvs/scooby-inf` had been upgraded
from 1.26.4 to 2.4.6 on 6 Nov 14:13, pulled in by the type-checking packages `scipy-stubs`, `optype` and
`numpy-typing-compat` (the last requires NumPy 2.4). They were removed, NumPy 1.26.4 (the version of all earlier
runs) reinstalled, and `numpy==1.26.4` pinned in `requirements-step2.txt`.

\newpage

# 6. Preparing CIMA cells

## 6.1 What the CIMA files contain

CIMA publishes its per-cell RNA data on CNGB (`ftp.cngb.org/pub/SciRAID/trueblood/cima/CIMA_Resource/Cell_Atlas/`):
one file with all 6,484,974 cells (44.4 GB listed; 47.6 GB served) and five lineage files (CD8 T 2,340,449;
CD4 T 2,320,103; Myeloid 754,770; B 566,569; NK 503,083 cells; together exactly 6,484,974). The B-cell file was
downloaded first and inspected (first remotely, then locally):

| Question | Answer |
|:---------------------------|:------------------------------------------------------------------|
| raw counts? | yes, in `layers/counts` (values 1, 1, 1 …; integers) |
| what is `X`? | normalised values (1.9869515 …), not used |
| gene names | **symbols** (`A1BG`, `A2M` …), 36,326 genes, unique |
| cell type | `obs['cell_type_l4']` (the 73 types of the eQTL table) |
| donor / library | `obs['sample']` (420 donors in this file), `obs['library']` (889) |
| cell order | grouped by library (label changes exactly 889 times); cell types mixed within a library (442,783 changes) |

## 6.2 Matching genes

The encoder's genes are OneK1K Ensembl IDs; CIMA uses symbols. Two problems were found and handled:

* OneK1K's symbols carry a **leading space** (`" HES4"`), which made a first direct match find 0 genes; spaces are
  stripped.
* Many genes were **renamed** since OneK1K's annotation. Each Ensembl ID is looked up by its **GENCODE v32 name**
  first (the annotation already used in this project), then by OneK1K's own symbol.

| How found in CIMA | Genes |
|:--------------------------------------------------------------|------:|
| GENCODE v32 name, same as OneK1K's symbol | 1,578 |
| GENCODE v32 name, renamed since OneK1K | 692 |
| OneK1K symbol only | 31 |
| **not found, set to 0** | **102** |
| total | 2,403 |

The 102 missing genes are mostly old clone-based lncRNA names (`RP11-…` 49, `AC…` 9, `LINC…` 5, `CTD-…` 5),
listed in `genes_matched.tsv`. Their effect is tested indirectly by checks C and D; a direct test is listed in
section 10.

## 6.3 Encoding

Script: `scooby_inference/embed_cima_cells.py`. For each CIMA cell: raw counts from `layers/counts` for the 2,301
found genes (the other 102 set to 0), log(1 + count), the final network → 10 numbers. The file is read in
stretches of 200,000 cells so memory stays bounded for the full file. Writes `embedding.pq` (same format as the
authors') and `celltype.pq` (one Scooby column per CIMA L4 type). The B-cell file: 566,569 cells in 0.8 min.

\newpage

# 7. Check C: where CIMA B cells land

For each CIMA cell: which OneK1K column centre (average of the authors' 10 numbers over that column's cells) is
nearest, compared with the column its type is mapped to in `cima_mapping_ontology.csv`
(`data/cima_cells/embedded_B/placement_per_type.csv`).

| CIMA type | Cells | Mapped to | Nearest centre most often | Share nearest the mapped column |
|:---------------------------|-------:|:---------------|:---------------|------:|
| Bn_TCL1A | 208,575 | B_cell | B_cell | 93.0% |
| Bn_IL6 | 14,221 | B_cell | B_cell | 92.0% |
| Bn_IFIT3 | 13,236 | B_cell | B_cell | 90.9% |
| Transitional_B_MME | 57,582 | B_cell | B_cell | 89.6% |
| Plasma_IGHA1 | 13,000 | plasmablast | plasmablast | 99.6% |
| Plasma_IGHG1 | 5,711 | plasmablast | plasmablast | 99.5% |
| Plasmablast_MKI67 | 3,256 | plasmablast | plasmablast | 95.2% |
| Switched_activated_Bm_CD86 | 32,792 | memory_B_cell | memory_B_cell | 93.5% |
| pre-Switched_Bm_IFIT3 | 1,823 | memory_B_cell | memory_B_cell | 78.1% |
| pre-Switched_Bm_JAM3 | 31,377 | memory_B_cell | memory_B_cell | 74.1% |
| Atypical_Bm_ITGAX | 28,714 | memory_B_cell | memory_B_cell | 64.5% |
| Unswitched_Bm_CD1C | 35,408 | memory_B_cell | memory_B_cell | 63.8% |
| Switched_Bm_IGHE | 33,882 | memory_B_cell | memory_B_cell | 62.1% |
| Switched_Bm_IGHDneg | 80,513 | memory_B_cell | memory_B_cell | 54.9% |
| Unswitched_Bm_IL6 | 6,479 | memory_B_cell | B_cell (54.7%) | 45.2% |
| **all B cells** | 566,569 | | | **80.7%** |

![Share of each CIMA type's cells nearest its mapped OneK1K column, at full depth and thinned to OneK1K depth (check D). Black lines: the same share for OneK1K's own cells of that column with the authors' embeddings (B cell 83.4%, memory B cell 65.1%, plasmablast 99.8%).](fig3_placement.png){width=100%}

Memory-B shares of 55–78% are normal for this space: OneK1K's own memory B cells are nearest the memory-B
centre only 65.1% of the time, because memory and naive B overlap.

**Closeness and spread compared with OneK1K's own cells.** The summary line of the run compares CIMA (0.046)
with OneK1K's all-type median (0.087), which is not like-for-like; per type:

| Cells | Median distance to nearest centre | Spread (median distance to own type centre) |
|:-------------------------------|-----:|-----:|
| OneK1K B cell | 0.037 | 0.039 |
| CIMA Bn_TCL1A (naive B) | 0.044 | 0.045 |
| OneK1K memory B cell | 0.042 | 0.049 |
| CIMA Switched_Bm_IGHDneg | 0.052 | 0.059 |
| OneK1K plasmablast | 0.090 | 0.090 |
| CIMA Plasma_IGHA1 | 0.135 | 0.133 |

CIMA naive and memory B cells are about as close to the OneK1K centres as OneK1K's own cells, and not squeezed
together (standard deviation of dimension 7: OneK1K 0.044, CIMA 0.056), which would be the sign of a network
outputting average values for unfamiliar input. CIMA plasma cells sit further out (0.13–0.23 vs 0.09), consistent
with them being more mature (IgA/IgG plasma cells, cycling plasmablasts) than OneK1K's small plasmablast group.

\newpage

# 8. Check D: sequencing depth

The encoder input is not normalised for depth (neither was scPoli's), so deeper sequencing shifts the input.
Median raw count per cell over the encoder's genes: **CIMA 306, OneK1K 233** (CIMA about 31% deeper). In
log(1 + count) terms, an expressed gene is up to about 0.27 higher.

**Test.** Each CIMA count was kept with probability 0.76 (binomial thinning), which brings CIMA to OneK1K's
depth (median 233), and the same 566,569 cells were re-encoded (`--thin 0.76 --compare_to …`;
`data/cima_cells/embedded_B_thin076/`).

| Measure | Full depth | Thinned to OneK1K depth |
|:-------------------------------------------------------------|:---------:|:---------:|
| share nearest a mapped column, all B cells | 80.7% | 78.4% |
| naive B → B_cell (Bn_TCL1A) | 93.0% | 92.9% |
| plasma → plasmablast | 95–99.6% | 95–99.5% |
| memory B → memory_B_cell | 45–94% | 39–90% |
| how far each cell moved: median / 90th percentile | — | 0.020 / 0.056 |
| cells keeping the same nearest centre | — | 93.0% |

The move (0.020) is about 2% of the distance between cell-type centres and about half the spread within a type
(0.045 for naive B): noticeable, but it does not move cells into other types. **Decision:** use full depth for the
main run, and keep thinning as a sensitivity check on the final eQTL scores. (The OneK1K depth median is from 20
random stretches of 2,500 cells and varies by about ±6% between runs: 233 and 247.)

\newpage

# 9. Getting cells of all 73 CIMA types

## 9.1 Why not the whole file, and why whole donors

* The full file (47.6 GB) downloads at about 0.05 MB/s per connection, and Scooby's run time grows with the
  number of cells; OneK1K uses 1.27 million cells, not 6.5 million.
* Scooby's per-type value is a sum over the type's cells, so a representative sample of each type suffices.
* Cells are stored grouped by library, with types mixed. Sampling "N cells per type" would touch almost every
  compressed block (nearly a full download). **All cells of a set of donors** are contiguous, so only their blocks
  need to be read, and each type keeps its natural proportion.
* Only `layers/counts` is needed (not the normalised `X`), roughly halving the download.

## 9.2 The downloader

Scripts: `scooby_inference/fetch_cima_donors.py` and `remote_h5.py`. It reads the file's labels remotely, chooses
donors (random with a seed, or named), looks up where each needed compressed block sits in the file, fetches the
blocks over 8 connections at once (at most 4 requests per second in total), decompresses them, and writes an `.h5ad` that `embed_cima_cells.py` reads
directly. Every fetched block is saved, so a stopped run continues where it left off. Cells are named
`row<n>` after their row in the full file.

**Problems met and fixed.** A first remote read of the file's structure stalled (130 blocks in 9 minutes without
finishing); reading only the needed parts instead (labels, gene names, a few values) took under a minute. The
download-size estimate looked up 50 blocks spread over the whole file and took 9.4 minutes; it now uses 20
consecutive blocks near the first chosen donor.

## 9.3 Check E: the downloader is exact

Two donors (`CIMA_H005`, `CIMA_H006`) were fetched remotely from the B-cell file and compared with the complete
local download of that file:

| Measure | Value |
|:----------------------------------|:------------------|
| cells | 4,075 (in 6 blocks) |
| nonzero counts | 6,571,772 in both |
| values that differ | **0** |
| donors and cell-type labels | identical |

## 9.4 The 60-donor plan (full file)

| | |
|:-------------------------|:----------------------------------------------------------|
| donors | 60 of 421, random (seed 0) |
| cells | **928,536** (14% of 6,484,974) in 114 contiguous blocks |
| cell types | all 73; each gets 9.6–21.1% of its cells (median 14.2%) |
| download | 1,279,168,060 nonzero counts, estimated 2.93 GB |

Types with fewer than 500 cells in the 60 donors (OneK1K's smallest column, platelet, has 516):

| Type | All CIMA cells | In the 60 donors |
|:--------------------------|------:|------:|
| Inflamed_NK_dim_IFIT1 | 204 | 43 |
| pre-Switched_Bm_IFIT3 | 1,823 | 218 |
| AS_DC | 1,662 | 268 |
| CD4_Tem_CCR5 | 2,634 | 403 |
| DC1_CLEC9A | 2,956 | 406 |
| ILC2_IL2RA | 3,042 | 493 (not mapped to OneK1K; dropped anyway) |

Their scores will be noisier and will be flagged; more of their cells can be added from other donors if needed.

\newpage

# 10. The 60-donor download: what happened

| Time (IST, 20 Nov) | Event |
|:------------|:-------------------------------------------------------------------|
| 04:14 | download starts on pasteur; labels, cell boundaries and all 9,574 count-value blocks (725 MB) are fetched |
| 10:07 | ftp.cngb.org starts refusing pasteur's connections (TLS error "wrong version number"), most likely a block after many requests; briefly lifted, then from ~11:00 permanent (still blocked at 19:30) |
| 11:30 | downloader rewritten to be gentler and resumable: looks up only the needed blocks, saves every block location and every block as it arrives, at most 4 requests per second, waits (up to 10 min) instead of failing; re-verified exact against the local B file |
| ~14:10 | download moved to a laptop on another network (scripts and the saved blocks copied there) |
| +57 min | all 19,034 gene-position blocks located |
| +147 to ~+155 min | server errors (HTTP 504 "gateway timeout", the server's own outage); the downloader waited and continued by itself |
| +271 min | 19,033 of 19,034 blocks fetched (~2.2 GB, ~0.16–0.19 MB/s), then the run was stopped by an accidental Ctrl+C |
| evening | rerun fetched the last block; the saved-blocks folder (2.9 GB) copied to pasteur |
| 19:41 | file built on pasteur, offline, in 3.5 min (3.0 GB) |

**Why the build moved to pasteur.** The 1.28 billion nonzero counts need about 10 GB in memory once
decompressed (4 bytes for the count + 4 for the gene number each), and the build briefly holds extra copies
(15–25 GB peak). The downloader was split into `--download_only` (laptop: fetch the blocks, save a small
`meta.npz` with donors, cell types and gene names) and `--build_offline` (pasteur, 251 GB RAM: build the `.h5ad`
from the saved blocks, no network). Tested first on 3 donors of the B-cell file: 3,382 cells, 5,106,974 values,
**0 differences**.

# 11. Check F: the 60-donor file is exact

`data/cima_cells/CIMA_60donors_seed0.h5ad`: **928,536 cells × 36,326 genes**, 60 donors, all 73 types,
1,279,168,060 nonzero counts, exactly as planned (Section 9.4; the six rare types have the planned counts).

All CIMA cell-type files come from the same atlas, so the B-lineage cells of these 60 donors are also in the
B-cell file downloaded separately (Section 6). Compared cell by cell:

| Measure | 60-donor file | B-cell file |
|:-----------------------------|----------:|----------:|
| B-lineage cells of the 60 donors | 70,394 | 70,394 |
| order of donors and types | identical | |
| nonzero counts | 101,375,855 | 101,375,855 |
| values that differ | **0** | |
| gene list (36,326) | identical | |

# 12. Check G: where cells of all 73 types land

All 928,536 cells were encoded with `embed_cima_cells.py` (1.2 min; output `data/cima_cells/embedded_60donors/`).
As in check C, each cell is assigned to the nearest OneK1K column centre (the average of the authors' 10 numbers
over that column's cells). New here: the **same test on OneK1K's own cells** (authors' embeddings), as a
baseline for how well the 22 columns can be told apart at all, and a **lineage-level** version (T / NK / B /
monocyte / DC).

| Measure | CIMA cells | OneK1K's own cells |
|:------------------------------------------|---------------:|---------------:|
| nearest the mapped (own) column | 49.4% | 63.7% |
| nearest a column of the right lineage | **96.7%** | **98.2%** |
| median distance to the nearest centre | 0.121 | 0.087 |

CIMA: the 880,048 cells of the 68 types that have a OneK1K mapping; OneK1K: all 1,267,768 cells.
By lineage (CIMA): B 99.9%, monocytes 99.5%, T 96.8%, NK 92.3%, DC 81.8%. Median counts over the encoder genes:
CIMA 241, OneK1K 233.

![Each CIMA type (cells in brackets). Top: share of cells nearest the OneK1K column the type is mapped to; black
line: the same share for OneK1K's own cells of that column. Bottom: share nearest a column of the same lineage.
Grey: megakaryocytes/platelets and stem cells.](fig4_placement73.png){width=100%}

**Monocytes, NK, B, plasma cells and pDCs** land in their mapped column about as often as OneK1K's own cells
(for example cMono_CD14 99.3% vs 99.5%; Mature_NK_dim 98.0% vs 88.5%; naive B 92.8% vs 83.4%).

**T-cell subtypes blur.** Naive CD4 T cells (CD4_Tn_CCR7) are nearest "CD4+ T cell" for 27% of cells (OneK1K's own
CD4+ T cells: 78%); the rest go to CD8+ T (24%), Treg (16%) and CD8+ TCM (10%). But OneK1K's own T columns are not
well separated either (CD4+ TCM 42%, CD8+ T 43%, CD8+ TCM 24%). The reason is the encoder's gene list. The genes
that classically separate T subtypes are not among the 2,403:

| Gene | Why it is not an encoder input |
|:------------------------------|:----------------------------------------------|
| CD4 | a Scooby **validation** gene (removed on purpose by the authors) |
| CD8A, CD8B, IL7R | Scooby **test** genes (removed on purpose) |
| CCR7, SELL, TCF7, LEF1, CD3E, FOXP3, CD40LG | not among the 3,000 most variable genes the authors selected |

So the authors' 10 numbers carry little CD4-vs-CD8 or naive-vs-memory information, for OneK1K cells as well. This is
a property of the authors' design, not of the substitute encoder (check A) or the download (check F). Consequence:
CIMA eQTL predictions are trustworthy per lineage and for monocyte, NK and B types, but differences between T-cell
subtypes will be small.

The weakest types at lineage level: DC_CSF2RA 45% (1,064 cells), HSPC_CD34 57% (678; OneK1K's own stem cells 55%),
Terminal_NK_dim_CD160neg 67%, gdT2_GZMH 71%.

# 13. Check H: Scooby with CIMA cells, test run

Config `step2_config_cima_cells.yaml` is the substitute-encoder config with two lines changed (the file formats
are identical to the authors' files):

* `embedding_path`: `data/cima_cells/embedded_60donors/embedding.pq`
* `cell_type_idx_path`: `data/cima_cells/embedded_60donors/celltype.pq` (73 CIMA types)

The first 50 CIMA pairs were run (`results/cima_cells_test/`):

* output 45 pairs × **73 columns**, all values finite; the same 45 pairs as in the earlier run (5 rows have no
  usable gene window in both runs);
* speed 1.7 sequences/s, the same as with OneK1K cells (the cell part is cheap; the DNA model dominates);
* each CIMA-cell column against the OneK1K-cell column its type is mapped to: **r = 0.996** over 3,105 values
  (cMono_CD14 vs monocyte 1.000, Bn_TCL1A vs B cell 1.000, Mature_NK_dim vs NK 1.000, MAIT 0.998, naive CD8 T 0.992,
  naive CD4 T 0.990).

![Test run: effect with CIMA cells (own CIMA type) against effect with OneK1K cells (mapped column).](fig5_test_run.png)

Predicted effects are dominated by the DNA sequence; the cells change them a little, most for T cells. So the
full run is expected to be close to the earlier CIMA results, with differences concentrated in T-cell types.

# 14. Check I: the full run, each CIMA eQTL in its own CIMA type

## 14.1 The run

Scooby on all **36,739 CIMA pairs** with CIMA cells (tmux `cima_full`, GPU 0), 19 chunks of 2,000 pairs, started
19:57 IST on 20 November, finished **07:14 IST on 21 November** (11.3 h). Every input row is covered by an output file.

```
bash run_all_chunks.sh 0 inputs/cima_pairs.tsv \
     step2_config_cima_cells.yaml results/cima_pairs_cima_cells
python score_cima_cells.py        # -> eval_out_cima_cells/
```

Scoring follows the paper's method (`evaluate_eqtl.py`, unchanged): a prediction is **clear** when
\|predicted effect\| > 0.05 (natural-log fold change); it is **right** when its sign matches the measured CIMA beta
(and \|beta\| > 0.05). The only difference from before: the mapping file now says *each CIMA type -> the Scooby
column of the same name* (`eval_out_cima_cells/mapping_own_column.csv`).

## 14.2 From CIMA's eQTL table to the scored rows

| Step | Count | Why it drops |
|:----------------------------------------------------------|----------:|:-------------------------------------------|
| CIMA lead eQTLs (`data/cima_eqtl.tsv`), 69 cell types | 59,961 | one row per (cell type, gene, lead variant) |
| distinct variant–gene pairs sent to Scooby | 36,739 | a pair can be the lead in several types |
| pairs Scooby could predict | 32,607 | variant outside the gene's ±262 kb model window |
| eQTL rows with a prediction (pairs back in their types) | **55,209** | |
| clear predictions (\|effect\| > 0.05) | 2,863 (5.2%) | Scooby predicts near-zero effects for most |
| clear and right sign | **2,512 (87.7%)** | |

The number of cells (928,536) does not enter these counts: all cells of one CIMA type are summed into that type's
column, and each column gives one prediction per pair.

## 14.3 Results by lineage (all 69 types)

| Lineage | Types | eQTL rows | Clear | Right sign of clear | Median Spearman (beta, prediction) |
|:----------------|------:|----------:|------:|--------------------:|-----:|
| T | 32 | 33,362 | 1,754 | 1,553 (88.5%) | 0.121 |
| NK | 5 | 3,678 | 167 | 147 (88.0%) | 0.082 |
| B | 14 | 7,374 | 362 | 310 (85.6%) | 0.143 |
| Monocyte | 8 | 6,515 | 401 | 347 (86.5%) | 0.173 |
| DC | 4 | 1,850 | 80 | 69 (86.3%) | 0.186 |
| ILC2, NKT (2), proliferating T | 4 | 2,052 | 88 | 78 (88.6%) | 0.142 |
| HSPC, megakaryocyte | 2 | 378 | 11 | 8 | |
| **All** | **69** | **55,209** | **2,863** | **2,512 (87.7%)** | **0.132** |

Per-type numbers: `eval_out_cima_cells/per_type_vs_onek1k_cells.csv`.

## 14.4 Check J: is a type's own column the best one for its eQTLs?

If the cells mattered for the prediction, the eQTLs of, say, cMono_CD14 should be predicted best by the cMono_CD14
column. For each type, its eQTLs were scored against **all 73 columns** (Spearman with beta, same pairs;
`specificity.csv`):

* the own column is the best of 73 for **1 of 69** types (pDC_IRF4); median rank **29**;
* the own column beats the median column of *other lineages* for 39 of 69 types (about half, i.e. chance);
* the 73 columns are very similar: median pairwise Pearson **0.958** over the 32,607 pairs (lowest 0.714).
  Of the 32,607 pairs, 1,111 are clear in at least one column, 220 in all 73.

So Scooby's *predicted* effects for these eQTLs are largely the same in every cell type; in the model, the
cell-type signal is small compared with the DNA-sequence signal. This describes the model's output only; whether the
true effects differ between cell types is tested against CIMA's measurements in Section 17. This agrees with the T-subtype blur in check G, and it means a
comparison of models (e.g. the planned CIMA fine-tune, tested on AIDA) should also look at specificity, not only at
overall sign agreement.

![Left: predicted effects with CIMA cells (own type) vs OneK1K cells (ontology-mapped column), same 53,157 eQTLs. Middle: per CIMA type, percent of clear predictions with the right sign, both runs (dot size = clear predictions). Right: rank of each type's own column among the 73.](fig6_full_run.png)

# 15. What changed from the previous CIMA results

The previous CIMA results (`scooby_gtex_onek1k_cima.pdf`) used Scooby with **OneK1K's cells** and mapped each CIMA
type onto OneK1K columns. Two mappings exist: the first one (`cima_mapping.csv`, 48 types, `eval_out_cima/`) and the
later ontology mapping (`cima_mapping_ontology.csv`, 69 types, 65 scored, `eval_out_cima_ontology/`). The comparison
below is against the ontology mapping unless stated.

## 15.1 What is different in the setup

| | Previous run | This run |
|:-------------------------|:-----------------------------------------|:-----------------------------------------|
| Cells filling the columns | OneK1K, 1,267,768 cells | CIMA, 928,536 cells from 60 donors |
| Columns | 22 OneK1K cell types | 73 CIMA L4 types |
| Cell embeddings | authors' scPoli | substitute encoder (checks A, B) |
| CIMA type -> column | ontology mapping: 69 types onto the 22 columns; 39 exact, 13 subtype-of, 11 manual, 5 related, 1 broader; 10 types averaged over several columns | none: each type its own column |
| Types sharing one column | many (e.g. all CD4 naive subtypes -> CD4+ T cell) | none |
| CIMA types scored | 65 | **69** (+ ILC2_IL2RA, NKT_IFNG, NKT_NCR1, Proliferative_T_MKI67) |
| eQTL rows scored | 53,157 | **55,209** (+2,052) |

Unchanged: the same Scooby model and weights, the same 36,739 input pairs, the same 32,607 predicted pairs, the same
scoring script and cutoffs.

## 15.2 What changed in the results

Same 65 types and the same 53,157 eQTL rows in both runs (`compare_ontology.csv`):

| Lineage | Clear, previous | Right sign, previous | Clear, now | Right sign, now | Pearson of effects |
|:---------|------:|--------------:|------:|--------------:|------:|
| T | 1,788 | 1,575 (88.1%) | 1,754 | 1,553 (88.5%) | 0.9985 |
| NK | 173 | 152 (87.9%) | 167 | 147 (88.0%) | 0.9991 |
| B | 372 | 313 (84.1%) | 362 | 310 (85.6%) | 0.9993 |
| Monocyte | 405 | 350 (86.4%) | 401 | 347 (86.5%) | 0.9988 |
| DC | 79 | 68 (86.1%) | 80 | 69 (86.2%) | 0.9983 |
| HSPC, MK | 8 | 5 | 11 | 8 | |
| **All 65** | **2,825** | **2,463 (87.2%)** | **2,775** | **2,434 (87.7%)** | **0.9985** |

Median Spearman (beta, prediction) per type: 0.128 previous, 0.132 now. Against the first 48-type mapping: 86.6%
(1,982 / 2,288) previous, 87.4% (1,959 / 2,241) now, on the same 42,713 rows (`compare_mapping48.csv`).

Two real rows (from `merged.csv.gz` in `eval_out_cima_cells/` and `eval_out_cima/`):

| CIMA type | Pair | Measured beta | Now (CIMA cells) | Previous (OneK1K cells) |
|:-----------|:--------------------------------|------:|------:|------:|
| cMono_CD14 | C19orf33, chr19_38304146_T_G | −0.768 | −1.1675 | −1.1674 |
| CD4_Tn_CCR7 | AC136475.9, chr11_327334_C_T | +1.251 | 0.0510 | 0.0479 |

The first is typical: the effect is the same to three decimals. The second shows where most of the small count
differences come from: an effect near the 0.05 cutoff lands on different sides in the two runs.

## 15.3 Reading

* The gain is mainly in **coverage and cleanliness**: no mapping decisions (11 were manual), 4 more types, and
  each type scored with its own cells.
* The scores move by +0.5 percentage points overall (largest in B: +1.5); with about 2,800 clear predictions this is
  within noise (one standard error of a proportion near 87% is about 0.6 points).
* The predicted effects are nearly identical (r = 0.9985), because they are dominated by the DNA sequence and are
  largely shared across cell types (check J).

## 15.4 Checks not yet done

| Check | Why |
|:-----------------------------------------------------------------|:-------------------------------------------------|
| effect of the 102 missing genes: set them to 0 in OneK1K cells and measure how far the embeddings move | direct test of the gene-matching loss |
| depth: rerun the CIMA scoring with thinned counts | does the 0.020 shift change any eQTL score? (given check J, unlikely) |
| cell numbers: rerun with, say, 2,000 cells per column | confirms a sample of cells per type gives the same effects |

## 15.5 Next: AIDA

AIDA's lead eQTLs (Cell 2025, Table S6: 143,918 rows, 57,875 significant, 20 cell types) and its cells
(CELLxGENE, 1,265,624 cells, raw counts) are downloaded (`data/aida/`). The same pipeline (encode AIDA cells,
Scooby, score in own type) can be run with the authors' model now. Fine-tuning on CIMA needs CIMA's raw reads
(GSA-Human HRA009968), which are controlled access.

\newpage

# 16. No clear prediction: small true effect, or the model?

Most eQTLs get no clear prediction (\|predicted effect\| ≤ 0.05), even close to the gene. In the earlier results
(`scooby_gtex_onek1k_cima.pdf`, pages 2, 4, 8), within 0–1 kb of the TSS a clear prediction is made for 50.5% of
GTEx pairs, 44.5% of OneK1K pairs (40.0% at 1–2 kb) and 25.9% of CIMA pairs. Note that these percentages are the
share of pairs that **get** a clear prediction; of the clear predictions, 83–93% have the right sign.

A missing prediction could mean the model fails, or that the true effect is small. To separate the two, each eQTL
set was sorted by the measured effect size \|beta\| and cut into **four equal-sized groups (quarters)**: quarter 1
holds the 25% of eQTLs with the smallest measured \|beta\|, quarter 4 the 25% with the largest. Each data set (and
its 0–1 kb subset) is split on its own, so the \|beta\| ranges differ between rows.

| Data | Distance | Quarter 1 (smallest \|beta\|) | Quarter 2 | Quarter 3 | Quarter 4 (largest \|beta\|) |
|:-----------------|:-------|:------------------|:------------------|:------------------|:------------------|
| GTEx whole blood | all | 16.7% (0.01–0.34) | 23.3% (0.34–0.49) | 26.9% (0.50–0.73) | 31.7% (0.74–2.38) |
| GTEx whole blood | 0–1 kb | 33.9% (0.06–0.35) | 47.3% (0.35–0.53) | 60.0% (0.54–0.74) | 60.7% (0.74–2.13) |
| OneK1K | all | 12.7% (0.00–0.38) | 21.0% (0.38–0.58) | 29.4% (0.58–0.85) | 28.0% (0.85–2.53) |
| OneK1K | 0–1 kb | 37.0% (0.17–0.45) | 60.3% (0.46–0.70) | 43.8% (0.71–1.00) | 37.0% (1.00–1.77) |
| CIMA (own cells) | all | 2.4% (0.08–0.32) | 3.4% (0.32–0.40) | 5.0% (0.40–0.55) | 9.9% (0.55–1.59) |
| CIMA (own cells) | 0–1 kb | 16.2% (0.11–0.34) | 20.9% (0.34–0.45) | 27.6% (0.45–0.64) | 35.7% (0.64–1.51) |

Cells: percent of eQTLs in that quarter with a clear prediction (in brackets: the \|beta\| range of the quarter).
Group sizes: GTEx 227 (0–1 kb: 55–56), OneK1K 313–314 (0–1 kb: 73), CIMA 13,802–13,803 (0–1 kb: 1,138). Of the
clear predictions, the share with the same sign as beta also rises with \|beta\|, e.g. GTEx 0–1 kb 78.9% in
quarter 1 and 100% in quarter 4 (all numbers: `eval_out_cima_cells/clear_by_beta_quarter.csv`).

**Reading.**

* **Partly a small true effect.** Larger measured effects get clear predictions more often (GTEx 16.7% → 31.7%,
  CIMA 2.4% → 9.9%), and those predictions are more often right. So some missing predictions reflect weaker effects.
* **Not only that.** Every eQTL here is significant, so none has a zero effect; OneK1K's smallest quarter near the
  gene still has \|beta\| 0.17–0.45. And large, nearby effects are missed as well: in GTEx's largest quarter within
  1 kb (\|beta\| 0.74–2.13) 39% get no clear prediction; in OneK1K's largest quarter (all distances, \|beta\|
  0.85–2.53) 72% get none. These cannot be called negligible, so model limits are part of the explanation.
* **Different scales.** beta is on the eQTL study's normalised-expression scale, the 0.05 cutoff on Scooby's
  log-coverage scale; "negligible" cannot be defined the same way on both sides, only the ordering can be compared.
* **Lead vs fine-mapped variants.** CIMA's eQTLs are **lead** variants (the most significant SNP per gene), not
  necessarily causal. If a lead SNP only tags the causal variant through linkage, its own DNA effect may truly be
  near zero, and Scooby's near-zero prediction would be right. This likely contributes to CIMA's much lower rate of
  clear predictions (about 5% vs 20–30% for the fine-mapped GTEx and OneK1K variants). AIDA also releases only lead
  variants, so the same applies there.

\newpage

# 17. Check K: does Scooby capture how an eQTL's effect differs between cell types?

Check J showed that Scooby's 73 columns give nearly the same predicted effect for each pair. That is a statement
about the model's output only; it cannot say whether the *true* effects differ between cell types. Check K compares
the predictions with CIMA's **measured** effects.

## 17.1 The question and how the test cases were built

**Question.** When CIMA measures the same eQTL with a different strength in two cell types, does Scooby's
prediction differ between those two cell types in the same direction?

**Cases.** A variant–gene pair appears in CIMA's table once for every cell type in which it is the lead eQTL.
**8,130** pairs are leads in two or more types. For each such pair, every combination of two of its types (A, B)
gives one case, **114,494 cases** in total. Each case has:

* measured: beta_A and beta_B, and their standard errors (CIMA's `slope_se`);
* predicted: Scooby's effect in column A and in column B (full run with CIMA cells, Section 14).

Measured difference = beta_A − beta_B; predicted difference = pred_A − pred_B. The direction is **right** when the
two differences have the same sign. A measured difference is called **significant** when it exceeds 3 standard
errors, \|beta_A − beta_B\| / sqrt(se_A² + se_B²) > 3, so that it is unlikely to be measurement noise (the median
standard error of a difference is 0.084).

**A real case** (CHI3L2, variant chr1_111225954_G_A):

| CIMA type | Measured beta | Scooby prediction |
|:---------------------------|------:|------:|
| Bn_TCL1A (naive B) | **−1.23** | −0.0148 |
| CD8_Tn_CCR7 (naive CD8 T) | **+1.05** | −0.0163 |

The measured effect changes sign between B and T cells; the difference is 27 standard errors. Scooby predicts
practically the same small value in both.

## 17.2 Results

Each row keeps fewer cases, with more reliable measured differences (`eval_out_cima_cells/check_k_filters.csv`). Right dir. = predicted difference has the same sign as the measured one; Spearman = between measured and predicted
differences; Med. meas. / Med. pred. = median absolute measured / predicted difference.

| Cases kept | Cases | eQTLs | Right dir. | Spearman | Med. meas. | Med. pred. |
|:--------------------------------------|-------:|------:|------:|------:|------:|------:|
| all | 114,494 | 8,130 | 49.5% | +0.005 | 0.105 | 0.0003 |
| \|measured diff\| > 0.3 | 15,081 | 1,366 | 50.4% | +0.059 | 0.397 | 0.0007 |
| measured diff significant | 19,589 | 1,797 | 51.1% | +0.064 | 0.353 | 0.0007 |
| … different lineages (B vs T, …) | 12,197 | 1,281 | 51.0% | +0.091 | 0.365 | 0.0011 |
| … same lineage (CD4 vs CD8 T, …) | 7,392 | 1,125 | 51.2% | +0.018 | 0.334 | 0.0003 |
| … measured signs opposite | 174 | 50 | 48.3% | −0.012 | 0.877 | 0.0016 |
| … predicted diff > 0.01 | 3,003 | 221 | 61.0% | +0.337 | 0.404 | 0.027 |
| … predicted diff > 0.05 | 827 | 62 | **82.2%** | +0.513 | 0.482 | 0.073 |

**Reading.**

* **Mostly no predicted difference.** The typical predicted difference (0.0003–0.001) is hundreds of times smaller
  than the measured one (about 0.35). The direction of a near-zero difference is essentially random, so most rows
  are at about 50%, including the 174 cases where the measured effect changes sign.
* **Lineages separate a little better than subtypes.** Spearman 0.091 between lineages vs 0.018 within one
  lineage, in line with check G (the embeddings separate lineages, but hardly T-cell subtypes).
* **When Scooby does predict a real difference, it is more often right.** For predicted differences > 0.05 the
  direction is right in 82.2% of cases.

## 17.3 The 82%, counted carefully

The 827 cases come from only **62 distinct eQTLs**, because an eQTL measured in many types gives many cases (LYZ
alone gives 170; `check_k_large_predicted.csv`). Counting each eQTL once, the majority of its cases is right for
**41 of 62 eQTLs (66%)**, binomial p = 0.015: better than chance, but much weaker than 82% suggests.

A right direction of the difference does not mean each prediction is right. LYZ (chr12_69338325_T_C):

| CIMA type | Measured beta | Scooby prediction |
|:-----------|------:|------:|
| cMono_CD14 | −0.875 | −0.011 |
| Bn_TCL1A | −0.328 | +0.063 |

The measured effect is stronger (more negative) in monocytes, and Scooby's monocyte value is lower than its B-cell
value, so the difference points the right way; but Scooby's B-cell value has the wrong sign.

## 17.4 Caveats

* **Selection.** Only pairs that are leads in several types are included, which favours eQTLs shared across types:
  measured beta_A and beta_B correlate at r = 0.958 and have the same sign in 99.8% of cases.
* **Measurement noise.** Handled by the significance filter (rows 3–8 of the table).
* **Lead variants.** If a lead SNP only tags the causal variant through linkage, part of a measured difference may
  come from linkage patterns, which Scooby cannot see in the DNA sequence.
* **Scales.** beta is on CIMA's normalised-expression scale, Scooby's effect on log coverage; only directions and
  rankings are compared, never sizes.
* Non-independence: one eQTL contributes many cases (handled in 17.3 by counting eQTLs).

## 17.5 Why: Scooby knows the cell types, but was never shown a variant

**Is it because Scooby never learned to tell cell types apart?** No. To check, a copy of `ScoobyVariant.py` that
also saves the predicted expression for both alleles was run on 6 of the eQTLs above, with CIMA cells
(`eval_out_cima_cells/levels_6pairs/`; the copy differs from the original by two lines,
`changes_vs_ScoobyVariant.diff`; its effects equal the full run, e.g. C19orf33 in cMono_CD14 −1.1675 in both).

**How to read the numbers.** For one column, e.g. cMono_CD14 (49,729 CIMA cells), Scooby predicts each cell's
LYZ RNA (read coverage over the gene) and adds them up over all cells, once with the reference allele and once
with the variant allele:

| cMono_CD14, LYZ | Total over 49,729 cells | Per cell |
|:----------------------------------|----------:|------:|
| DNA with the reference allele (T) | 2,334,947 | **46.95** |
| DNA with the variant allele (C) | 2,309,865 | 46.45 |

*Per cell* = total / number of cells: how much LYZ Scooby expects in a typical cell of this type. *Variant effect*
= ln(2,309,865 + 1) − ln(2,334,947 + 1) = **−0.011** (about 1% less). For naive B cells (Bn_TCL1A, 24,669 cells)
the same gives 0.29 per cell (REF) and 0.31 (ALT), effect +0.063.

| CIMA type | Predicted LYZ per cell (reference allele) | Variant effect (chr12_69338325_T_C) |
|:---------------|------:|------:|
| cMono_CD14 | **46.95** | −0.011 |
| ncMono_FCGR3A | 21.40 | +0.009 |
| pDC_IRF4 | 8.22 | +0.031 |
| Bn_TCL1A | **0.29** | +0.063 |
| CD4_Tn_CCR7 | 0.15 | +0.048 |

* **Expression level: strongly cell-type specific.** Scooby expects about **160 times** more LYZ per cell in
  classical monocytes than in naive B cells (over all 73 columns: 738-fold range). That is correct biology: LYZ is
  a monocyte gene. The cell embeddings are therefore used, and the cell types are distinguished.
* **Variant effect: near zero everywhere**, while CIMA measures −0.875 in cMono_CD14.

**The model can give different effects per cell type**, for some eQTLs:

| Variant effect | cMono_CD14 | ncMono | Bn_TCL1A | CD4_Tn | NK_dim |
|:-------------------------------|------:|------:|------:|------:|------:|
| HKDC1 chr10_69232960_G_A | −0.86 | −0.98 | −1.14 | **−1.53** | −1.51 |
| C19orf33 chr19_38304146_T_G | **−1.17** | −1.15 | −0.67 | −0.69 | −0.83 |

So different effects are possible; they are just mostly small or do not match the measurements.

![Left: measured vs predicted difference between two CIMA types for the same eQTL (significant measured differences only); most predicted differences are near zero. Right: 6 eQTLs, each dot one CIMA column: predicted expression per cell (reference allele, log scale) against the predicted variant effect. LYZ spans a 738-fold range of expression with near-zero effects; HKDC1 and C19orf33 have large effects that vary between columns.](fig7_check_k.png)

**What training never showed it.** Scooby was fine-tuned on the reads of all OneK1K donors aligned to **one
reference genome** (`genome_human.fa` in `scripts/training_onek1k/config_rna_only.yaml`); the donors' own genotypes
were not used. It therefore learned how DNA sequence relates to expression in each cell type by comparing
**different genes** (e.g. "this promoter is active in monocytes"), but never saw "this base changed, and expression
changed in this cell type". Every variant effect is an extrapolation.

**Why the extrapolation tends to be shared.** In the model code (`scooby/modeling/scooby.py`,
`forward_convs_on_emb`) the DNA part (Borzoi, with LoRA) is the same for all cells; each cell only adds its own
final layer: a weighted sum of the same 1,920 DNA features plus a bias, followed by softplus. The weights come from
the cell's 10-number embedding (`forward_cell_embs_only`). A SNP changes the shared DNA features, and every cell
type reads that same change with different weights. The effect is a log ratio, ln(ALT) − ln(REF): where a cell
type's weights mainly scale a gene up or down, the scale cancels in the ratio. LYZ is this case: a 160-fold level
difference between monocytes and B cells, about zero effect in both.

## 17.6 What this means for the planned CIMA fine-tune and AIDA test

* Fine-tuning on CIMA in the same way (reads aligned to the reference genome) would teach the model CIMA's cell
  types and expression levels, and may sharpen T-cell subtypes, but it still would not show the model any variant
  effect. A large improvement of cell-type-specific eQTL predictions on AIDA is therefore not expected from that alone.
* Learning cell-type-specific variant effects directly would need training on each donor's own DNA sequence
  (personal genomes), i.e. the donors' genotypes, which are controlled access for both CIMA and AIDA.
* Metric for the comparison of models: besides overall sign agreement, the two numbers of this check, i.e. in how
  many cases a model predicts a real cell-type difference, and how often its direction is right.

CIMA numbers in this section: `python check_k_celltype_differences.py` (log `eval_out_cima_cells/check_k.log`).

## 17.7 The same test on OneK1K eQTLs

To check that the CIMA result is not caused by the substitute encoder, by CIMA's cells or by lead variants, the
test was repeated on **OneK1K**, with the authors' setup.

**Subset.** OneK1K's fine-mapped eQTLs (`data/all_highpip_susie.tsv`, study OneK1K, PIP ≥ 0.9, i.e. each variant
very likely the causal one): 762 variant–gene pairs have a Scooby prediction, and **265** of them are eQTLs in two
or more OneK1K cell types. Every combination of two of their cell types gives **1,507 cases**. Predictions: the
existing run with the **authors' embeddings and OneK1K's own cells** (`results/gtex_onek1k_pairs/`, 22 columns); no
new GPU run was needed. A pair listed twice for one cell type (two credible sets) keeps the higher-PIP row. The run
with the substitute embeddings gives the same numbers to 3 decimals (`eval_onek1k/check_k_substitute/`).

So this subset removes three possible explanations at once: fine-mapped instead of lead variants, the authors'
embeddings instead of the substitute encoder, and the very cells the model was trained on instead of new cells.

**A real case** (RPL28, variant chr19_55385959_A_G):

| OneK1K type | Measured beta | Scooby prediction |
|:-------------|------:|------:|
| CD4+ T cell | **−1.09** | −0.0081 |
| plasmablast | **−0.25** | −0.0064 |

The values are the effects **in each cell type**. The differences are: measured −1.09 − (−0.25) = **−0.84**
(23 standard errors, the effect is four times stronger in CD4+ T cells); predicted −0.0081 − (−0.0064) = **−0.0017**.
Both are negative, so this case counts as *right direction*, but the predicted difference is about 500 times smaller
than the measured one: Scooby predicts practically the same small value in both.

**Results, next to CIMA** (same definitions as in 17.1–17.2; `eval_onek1k/check_k/filters.csv`):

| Cases kept | OneK1K cases (eQTLs) | OneK1K right dir. | OneK1K Spearman | CIMA right dir. | CIMA Spearman |
|:----------------------------|------------:|------:|------:|------:|------:|
| all | 1,507 (265) | 52.2% | 0.124 | 49.5% | 0.005 |
| measured diff significant | 624 (94) | 57.1% | 0.201 | 51.1% | 0.064 |
| … different lineages | 377 (60) | 63.4% | 0.253 | 51.0% | 0.091 |
| … same lineage | 247 (65) | 47.4% | 0.016 | 51.2% | 0.018 |
| … measured signs opposite | 0 | – | – | 48.3% | −0.012 |
| … predicted diff > 0.05 | 24 (9) | 58.3% | −0.046 | 82.2% | 0.513 |

Typical size of the difference (all cases): measured 0.163 in OneK1K (CIMA 0.105), predicted **0.0012** (CIMA
0.0003), i.e. more than 100 times smaller than measured. Measured effects of the same eQTL in two OneK1K types
correlate at r = 0.965 (same sign 99.9%); no OneK1K case has a significant sign flip.

**Counted per eQTL.** One eQTL measured in many cell types gives many cases, so the OneK1K percentages above rest on
few eQTLs. Counting each eQTL once (is the majority of its cases right?), and giving the Spearman correlation a 95%
interval by resampling eQTLs (300 bootstrap samples; `eval_onek1k/check_k/per_eqtl_summary.csv`):

| Cases kept | Data | eQTLs with majority right | p (binomial) | Spearman, 95% interval |
|:----------------------------|:-------|---------------:|------:|:----------------|
| all | OneK1K | 135 / 257 (52.5%) | 0.45 | 0.124 (−0.015 to 0.230) |
| | CIMA | 3,834 / 7,896 (48.6%) | 0.01 | 0.005 (−0.037 to 0.046) |
| significant | OneK1K | 48 / 88 (54.5%) | 0.46 | 0.201 (−0.037 to 0.362) |
| | CIMA | 852 / 1,681 (50.7%) | 0.59 | 0.064 (−0.046 to 0.187) |
| significant, different lineages | OneK1K | 33 / 57 (57.9%) | 0.29 | 0.253 (−0.035 to 0.437) |
| | CIMA | 603 / 1,220 (49.4%) | 0.71 | 0.091 (−0.034 to 0.223) |
| significant, same lineage | OneK1K | 28 / 59 (47.5%) | 0.80 | 0.016 (−0.191 to 0.251) |
| | CIMA | 546 / 1,061 (51.5%) | 0.36 | 0.018 (−0.083 to 0.114) |

(eQTLs whose cases split exactly half right, half wrong are left out of the count.)

**Reading.**

* OneK1K looks slightly better than CIMA, mainly **between lineages** (63.4% of cases right). Counted per eQTL this is
  33 of 57 (p = 0.29), and every 95% interval includes 0: with 60–94 eQTLs, the difference from chance is not
  significant. **Within a lineage** (e.g. CD4 vs CD8 T) there is no signal at all, as in CIMA.
* Where Scooby predicts a sizeable difference, the rate of right directions is the same as in CIMA, but OneK1K has
  too few such eQTLs to be sure (table below).
* **Conclusion:** the CIMA result holds for OneK1K. Even with fine-mapped variants, the authors' embeddings and the
  cells the model was trained on, Scooby's predicted differences between cell types are about 100 times smaller than
  the measured ones, and their direction is not reliably better than chance. This is a property of the model (trained
  on the reference genome, never shown a variant; 17.5), not of the substitute encoder, the CIMA cells or the lead
  variants.

**When Scooby does predict a sizeable difference.** Two different numbers (cases with a significant measured
difference and \|pred_A − pred_B\| > 0.05; each eQTL counted once, as right when most of its cases are right):

| | OneK1K | CIMA |
|:------------------------------------------|:-------------------------|:-------------------------|
| (a) how often Scooby predicts a sizeable difference between two cell types | 9 of 265 eQTLs = **3.4%** | 62 of 8,130 eQTLs = **0.8%** |
| (b) when it does, how often the direction is right | 6 of 9 = **67%** | 41 of 62 = **66%** |
| chance of doing at least this well by guessing (coin flip) | about 1 in 4 | about 1 in 132 |
| p-value (two-sided binomial test) | 0.51 | 0.015 |

The 9 OneK1K eQTLs: right for SLC12A7 (4 of 5 cases), ITGB2 (3/3), PPFIBP2 (2/2), SDK2 (2/2), MTUS1 (2/2) and
EML2-AS1 (1/1); wrong for CALHM6 (0/6), AC024028.1 (0/2) and ACCS (0/1).

**How to read the p-value.** Suppose Scooby had no knowledge of cell-type differences, so that each direction were a
pure guess, like tossing a coin. The p-value is how often a coin would do at least as well (counting the equally
unusual opposite outcome too). A coin tossed 9 times gives 6 or more heads about 1 time in 4 (p = 0.51): 6 of 9 is
ordinary for a coin and does not show that Scooby knows anything. A coin tossed 62 times gives 41 or more heads about
1 time in 132 (p = 0.015): rare for a coin, so on these 62 CIMA eQTLs Scooby probably does better than guessing. The
rate is the same in both data sets (about two-thirds); only the number of eQTLs differs, and with few eQTLs a lucky
streak easily looks like skill. The usual threshold is p < 0.05 (less than 1 chance in 20).

**In short:** Scooby rarely predicts a sizeable difference between cell types (3.4% of OneK1K eQTLs, 0.8% of CIMA
eQTLs); when it does, the direction is right about two-thirds of the time in both data sets, which is significantly
better than guessing in CIMA (62 eQTLs) but cannot be told apart from guessing in OneK1K (9 eQTLs).

OneK1K numbers: `python check_k_onek1k.py` (log `eval_onek1k/check_k_authors.log`; with
`--pred results/gtex_onek1k_pairs_substitute/outs_0_1584.h5ad --out eval_onek1k/check_k_substitute` for the
substitute embeddings).


\newpage

# 18. Pre-processing of the eQTL inputs: OneK1K, CIMA and AIDA

Before Scooby can score an eQTL, the published table has to become a list of variant–gene pairs in the format
`ScoobyVariant.py` reads (`Chromosome, gene_name, variant, Pos, Ref, Alt, gene_id`). Four steps decide which eQTLs
enter. The three datasets needed different amounts of work, because they were published in different forms.

| Step | OneK1K (+ GTEx) | CIMA | AIDA |
|:--------------------|:--------------------------------|:--------------------------------|:--------------------------------|
| 1. significance filter | different kind: fine-mapped variants, PIP $\geq$ 0.9 (1,412 OneK1K rows; lowest PIP 0.9004) | q < 0.05 (largest `study_wise_qval` in `data/cima_eqtl.tsv`: 0.049999) | CellSubtype FDR < 0.05 (Section 19.2) |
| 2. single-base variants only | yes: 0 indels among the 1,585 input pairs | yes: 0 indels among the 36,739 | yes: 3,827 indel rows removed, 54,048 kept |
| 3. gene symbol → Ensembl ID | not needed: the eQTL Catalogue gives Ensembl IDs | not needed: `cima_eqtl.tsv` already has `gene_id` | **new**: Table S6 gives only symbols (Section 19.3) |
| 4. Scooby window | done beforehand by the authors: their input list (1,585 pairs) leaves out the 56 of 818 OneK1K pairs outside the window | not done beforehand; 36,739 in, 32,607 predicted | done beforehand, same rule: 27,857 pairs in |

**1. Significance filter.** This decides how sure we are that a row is a real eQTL with the right variant.

- **OneK1K and GTEx:** these come from the eQTL Catalogue's fine-mapping (SuSiE). A row is kept only if the
  variant has a posterior inclusion probability of at least 0.9, i.e. it is very likely the causal variant itself.
  This is the authors' own selection (`data/all_highpip_susie.tsv`).
- **CIMA and AIDA:** they release only **lead** variants (the most significant variant per gene and cell type),
  not fine-mapped ones. The filter is a false discovery rate. A lead variant can be a neighbour in linkage with
  the causal one, so these sets are noisier than the PIP $\geq$ 0.9 set however strict the FDR.

**2. Single-base variants.** `ScoobyVariant.py` builds the ALT sequence by swapping one base, so insertions and
deletions are left out. Only AIDA's table contained any.

**3. Gene IDs.** Scooby's gene table (`genes_tsscentred.bed`) and GTF use GENCODE v32 names and IDs. For AIDA,
each symbol is converted to an Ensembl ID, then to the GENCODE v32 name Scooby looks up.

**4. Scooby window.** Scooby predicts the gene's expression from a 524 kb sequence. `ScoobyVariant.py` keeps a pair
only if the variant lies between the gene region's start − 163,840 and end + 163,840 (`genes_tsscentred.bed`
coordinates), and drops the rest itself.

- **OneK1K (+ GTEx):** the authors' input list `inputs/gtex_onek1k_pairs.tsv` was already filtered. 762 of the 818
  OneK1K pairs are in it; the other 56 are outside the window (55) or have no gene region (1). All 1,585 listed
  pairs were predicted.
- **CIMA:** Scooby dropped 4,132 of 36,739 pairs after loading them.
- **AIDA:** I applied the same rule first (copied from `ScoobyVariant.py`), because about 20,000 of the 47,953
  pairs fall outside. The result is identical; only GPU time is saved.

## 18.1 How many eQTLs each dataset gives

| | OneK1K | CIMA | AIDA |
|:------------------------------------------|------:|------:|------:|
| eQTL rows after steps 1–3 (cell type × pair) | 1,412 | 59,961 | 54,027 |
| distinct variant–gene pairs | 818 | 36,739 | 47,953 |
| … inside Scooby's window | 762 | 32,607 | 27,857 |
| cell types with eQTLs | 21 | 69 | 20 |
| eQTL rows with a prediction | 1,255 | 55,209 | 33,751 |

OneK1K has **40–50 times fewer** eQTLs than CIMA or AIDA. The window is not the reason (it removes 7% of OneK1K
pairs, 11% of CIMA's, 42% of AIDA's). The reason is step 1:

- **OneK1K:** keeps only variants that fine-mapping calls almost certainly causal (PIP $\geq$ 0.9). Most
  significant OneK1K eQTLs have no single variant that certain, because the causal variant cannot be separated
  from its neighbours in linkage, so they are not in the set at all.
- **CIMA and AIDA:** keep one lead variant for every significant gene in every cell type, whether or not that
  variant is the causal one.

The trade-off is size against certainty. OneK1K's 1,255 rows are few but each variant is very likely causal. CIMA's
and AIDA's tens of thousands of rows include many lead variants that only tag the causal variant. Scooby's
prediction for a tagging variant can be near zero even when the gene has a strong eQTL. This is one reason why
clear predictions are rarer on CIMA than on OneK1K.

## 18.2 Effect sizes are on different scales

The measured effect (beta) is the change in expression per copy of the ALT allele, in whatever expression unit
each study used.

| | OneK1K | CIMA | AIDA (scored rows) |
|:------------------------------------|------:|------:|------:|
| median \|beta\| | 0.592 | 0.397 | **0.029** |
| middle half of \|beta\| (25th–75th percentile) | 0.388–0.868 | 0.314–0.538 | 0.018–0.051 |
| smallest \|beta\| | 0.0037 | 0.0754 | 0.0033 |
| share with \|beta\| > 0.05 | 99.4% | **100%** | **25.4%** |

(OneK1K: `data/all_highpip_susie.tsv`; CIMA: `data/cima_eqtl.tsv`; AIDA: `data/aida_eqtl.tsv`.)

AIDA's betas are about 14 times smaller than CIMA's. They are significant, so the effects are not weak; they are
measured in a different unit. OneK1K's and CIMA's expression values are most likely standardised before testing,
so beta is in standard deviations. AIDA's are most likely not. This should be confirmed in AIDA's Methods.

**Why it matters for the score.** The paper's sign metric (`sign_agree_nonneg` in `evaluate_eqtl.py`, used for
every "clear predictions with the right sign" number in this report) counts a clear prediction (|pred| > 0.05) as
right only if the measured |beta| is **also** > 0.05 and the signs agree. The denominator is all clear
predictions.

| | OneK1K | CIMA | AIDA |
|:--------------------------------------------------------------|------:|------:|------:|
| rows with \|beta\| > 0.05, i.e. rows where a prediction can count as right | 99.4% | 100% | 25.4% |

For OneK1K and CIMA the beta cut-off changes almost nothing. For AIDA it marks about 3 in 4 clear predictions wrong
whatever their sign.

Two real AIDA rows:

| Cell type | Variant | Gene | beta | p-value | Can a right-sign prediction count? |
|:----------------|:--------------------|:------|------:|------:|:------------------|
| atypical_B | chr1_3125817_G_A | PEX10 | 0.056 | 1.5 × 10⁻⁵ | yes (\|beta\| > 0.05) |
| CD14+_Monocyte | chr7_76404639_T_C | ZP3 | 0.029 | 1.5 × 10⁻³⁸ | **no**: \|beta\| < 0.05 |

ZP3 is one of the strongest eQTLs in the table (p = 1.5 × 10⁻³⁸), yet its beta is the AIDA median. Under the
paper's metric even a correct clear prediction for it would count as wrong.

**How AIDA will be scored.** Both numbers will be reported:

1. the paper's metric unchanged, for comparability, with the caveat above;
2. **plain sign agreement**: among clear predictions, the share whose sign equals beta's sign, at any |beta|. This
   does not depend on the unit. For CIMA it is almost identical to the paper's metric, since every CIMA |beta| is
   above 0.05.

The results are in Section 19.8: on AIDA the two numbers are 38.4% and 84.6%.

\newpage

# 19. AIDA: steps and results (authors' model, AIDA's own cells)

**Aim.** Score AIDA's eQTLs (Asian Immune Diversity Atlas; Cell 2025, doi 10.1016/j.cell.2025.02.017) with the
authors' OneK1K Scooby model, using **AIDA's own cells** in each AIDA cell type, exactly as Section 14 did for CIMA.
This is the baseline for the planned comparison with a CIMA fine-tuned model. That fine-tune still needs CIMA's raw
reads (Section 15.5).

## 19.1 Data

| Item | File (under `scooby_inference/data/aida/`) | Content |
|:----------------|:----------------------------------------|:------------------------------------------------|
| eQTLs | `cell2025_supplementary/` `TableS6_AIDA_lead_eQTLs.tsv` (from `mmc4.xlsx`) | 143,918 rows: one lead variant per gene per cell type; 20 cell types (Annotation_Level3 names); variants `chr_pos_ref_alt`, hg38; gene **symbols** |
| cells | `cells/AIDA_DataFreeze_v2.h5ad` (CELLxGENE, 14.24 GB) | 1,265,624 cells × 35,477 genes (Ensembl IDs); raw integer counts in `raw/X`; 625 donors; 10x **5′** v2 (OneK1K: 3′) |
| genotypes | not downloaded (managed access) | not needed: Scooby never uses genotypes |

## 19.2 Which eQTLs: the two FDR columns

Table S6 has two FDR columns. **Every** one of the 143,918 rows already has `eigenMT_FDR` < 0.05.
`eQTL_CellSubtype_FDR` is a stricter cut.

| | eigenMT FDR | CellSubtype FDR |
|:------------------|:-------------------------------------------|:-------------------------------------------|
| test unit | one test per gene: the lead variant's p-value × number of independent variants near the gene (`eigenMT_pvalue` = `eQTL_p.value` × `eigenMT_tests`, exact for all rows) | one test per variant–gene pair: ranks rows exactly as the raw p-value (Spearman 1.0 in every type) |
| question | does this gene have an eQTL in this type? | is this variant–gene association significant? |
| rows < 0.05 | 143,918 (all) | **57,875** (11,431 genes, the paper's eGene count) |

Two real rows (both atypical_B):

| Variant | Gene | beta | p-value | eigenMT FDR | CellSubtype FDR | Kept |
|:-------------------|:----------|------:|--------:|-------:|-------:|:----|
| chr1_3125817_G_A | PEX10 | 0.056 | 0.000015 | 0.0029 | **0.040** | yes |
| chr1_975014_C_T | LINC01128 | −0.102 | 0.0013 | 0.032 | **0.566** | no |

| eigenMT FDR of a row | ≤ 0.001 | 0.001–0.01 | 0.01–0.03 | 0.03–0.05 |
|:----------------------------------|----:|----:|----:|----:|
| share also passing CellSubtype FDR < 0.05 | 100% | 89% | 32% | 1.5% |

**Chosen: CellSubtype FDR < 0.05.** Strong eQTLs pass both; the difference is weak eQTLs whose lead variant is an
uncertain pick. Neither column could be recomputed from the table, because genes that failed were not released
(for CD14+ monocytes, eigenMT FDR implies about 17,600 genes tested; 8,021 are in the table). The exact definition
of CellSubtype FDR should be confirmed in the paper's Methods.

## 19.3 From Table S6 to Scooby's input

Script: `scooby_inference/make_aida_inputs.py`; counts in `data/aida_input_funnel.csv`.

| Step | eQTL rows | Pairs | Genes |
|:------------------------------------------------|------:|------:|------:|
| Table S6, all rows (all pass eigenMT FDR < 0.05) | 143,918 | 135,718 | 12,117 |
| CellSubtype FDR < 0.05 | 57,875 | 51,373 | 11,431 |
| single-base variants | 54,048 | 47,972 | 11,347 |
| gene matched to a GENCODE v32 ID | 54,027 | 47,953 | 11,343 |
| **variant inside the gene's Scooby window** | **33,751** | **27,857** | **9,034** |

**Gene matching.** Each symbol is looked up first in AIDA's own gene table (`raw/var/feature_name` → Ensembl ID
in the h5ad), then as a unique GENCODE v32 name. 10,504 genes were found the first way, 839 the second, and 4 not
at all. Scooby then uses the GENCODE v32 name of that ID.

**Window.** 42% of the matched pairs (20,096 of 47,953) fall outside Scooby's window: AIDA tested variants up to 1 Mb away, and
the median lead variant is 158 kb from the gene start.

**Format.**

- `Pos` = position − 1 (0-based), as in `inputs/cima_pairs.tsv` (`chr1_814733_A_T` → `Pos` 814732).
- `data/aida_eqtl.tsv` has the same columns as `data/cima_eqtl.tsv`: `beta` = `eQTL_beta`, `slope_se` = beta / t,
  `af` and `sample_size` empty.
- Outputs: `data/aida_eqtl.tsv` (33,751 rows) and `inputs/aida_pairs.tsv` (27,857 pairs).

Rows per cell type after all steps: CD14+_Monocyte 3,993, CD4+_T_cm 3,665, CD4+_T_naive 3,543, CD16+_NK 3,474,
CD16+_Monocyte 2,851, CD8+_T_naive 2,264, CD8+_T_GZMBhi 2,016, CD8+_T_GZMKhi 1,775, cDC2 1,609, CD4+_T_em 1,359,
naive_B 1,255, memory_B_IGHMlo 1,187, MAIT 996, Treg 825, gdT_GZMKhi 802, memory_B_IGHMhi 789, CD4+_T_cyt 610,
gdT_GZMBhi 343, CD56+_NK 224, atypical_B 171.

## 19.4 Encoding AIDA's cells

Script: `scooby_inference/embed_aida_cells.py`; output `data/aida/embedded/`; 2.5 min for all cells.

**Cells.** The 6 LONZA donors are left out: they are European controls (87,252 cells), not part of AIDA's Asian eQTL
cohort. That leaves **1,178,372 cells from 619 donors**, in 35 Annotation_Level3 types, all 20 eQTL types among
them. As for CIMA, all cells of one type are summed into that type's Scooby column (`celltype.pq`).

**Genes.** The encoder reads 2,403 genes, given as OneK1K Ensembl IDs. AIDA has Ensembl IDs too, so most match
directly.

| How found in AIDA | Genes |
|:------------------------------------------|------:|
| Ensembl ID | 2,206 |
| GENCODE v32 name | 3 |
| OneK1K symbol | 19 |
| **not found, set to 0** | **175** (7.3%; CIMA: 102) |
| total | 2,403 |

The 175 are mostly old clone-based names: `RP…-` 81, `AC/AL/AP…` 30, `CT…-` 16, `LINC…` 6, other 42
(`genes_matched.tsv`). Encoder input: log(1 + raw count) from `raw/X`, as in training. The script asserts that the
counts are integers.

## 19.5 Placement check (as check G)

Each AIDA cell goes to the nearest OneK1K column centre (the average of the authors' 10 numbers over that column's
cells). The check counts how often that centre is the column the AIDA type corresponds to (`aida_mapping.csv`)
and how often it is at least the right lineage. OneK1K's own cells are the baseline.

| Measure | AIDA | CIMA (check G) | OneK1K's own cells |
|:-------------------------------------------|------:|------:|------:|
| nearest the mapped column | 73.3% | 49.4% | 63.7% |
| nearest a column of the right lineage | **97.3%** | 96.7% | 98.2% |
| median distance to the nearest centre | 0.107 | 0.121 | 0.087 |

By lineage (AIDA): B 99.98%, monocytes 99.5%, T 98.2%, NK 90.5%, DC 89.7%.

| AIDA eQTL type | Cells | Mapped to | Nearest that column | OneK1K's own cells | Right lineage |
|:----------------|------:|:----------------|------:|------:|------:|
| CD14+_Monocyte | 183,568 | monocyte | 99.5% | 99.5% | 99.9% |
| CD16+_Monocyte | 42,726 | CD16+ monocyte | 95.4% | 96.4% | 99.9% |
| naive_B | 31,744 | B cell | 97.6% | 83.4% | 100% |
| CD16+_NK | 146,608 | NK cell | 89.8% | 88.5% | 94.0% |
| MAIT | 21,251 | MAIT cell | 88.7% | 83.3% | 99.9% |
| CD4+_T_naive | 145,215 | CD4+ T cell | 65.7% | 77.7% | 100% |
| CD8+_T_GZMKhi | 27,444 | CD8+ TEM / TCM | 20.8% | – | 99.3% |

The mapped-column share is higher than CIMA's partly because AIDA's Level3 types are coarser (35 vs 73), so each
matches one OneK1K column more directly. The weakest eQTL type, CD8+_T_GZMKhi, most often lands nearest CD4+ CTL,
but its cells are still T cells. **Conclusion:** the 175 missing genes do not spoil the embeddings, so the AIDA
columns can be used. T subtypes blur as before (Section 12).

## 19.6 Depth check (as check D)

"Encoder-gene counts" means a cell's raw counts summed over the 2,403 encoder genes, i.e. how much signal the
encoder sees. The encoder input is not normalised for depth.

| | Median encoder-gene counts per cell |
|:---------------------|------:|
| OneK1K (training data) | 233 |
| CIMA | 241 |
| AIDA | **586** (2.5× OneK1K) |

**Test.** Each AIDA count was kept with probability 0.40 (binomial thinning, as if sequenced less deeply), and the
cells were re-encoded and placed again (`embed_aida_cells.py --thin 0.40`; output
`data/aida/embedded_thin040/`, per type in `moved_per_type.csv`).

| | Full depth | Thinned to OneK1K depth |
|:----------------------------------|------:|------:|
| median encoder-gene counts | 586 | 234 |
| nearest the mapped column | 73.3% | 69.9% |
| nearest a column of the right lineage | 97.3% | 96.8% |
| cells moved (median / 90th percentile) | – | 0.059 / 0.128 |
| same nearest column as at full depth | – | 79.9% |

Types with one clear column barely change: CD14+_Monocyte 99.5% → 99.3%, naive_B 97.6% → 99.1%, CD16+_NK
89.8% → 87.3%. The memory B subtypes drift toward the naive B column (memory_B_IGHMhi 64.6% → 29.1% nearest memory
B). Cells move more than CIMA's did in check D (median 0.059 vs 0.020), but most of the 20% that change column stay
in their lineage.

**Reading.** The 2.5× depth does not break lineage placement, and full depth places cells slightly better, so the
full-depth embeddings are used. Whether a move of about 0.06 changes Scooby's predicted *effects* is not tested
yet. It would need a second Scooby run on a sample of pairs with the thinned embeddings.

## 19.7 The Scooby run

| Item | Value |
|:----------------------|:---------------------------------------------------------------------|
| model | authors' `lauradmartens/onek1k-scooby`, unchanged |
| config | `scooby_inference/step2_config_aida.yaml` (as `step2_config_cima_cells.yaml`, with the AIDA `embedding.pq` and `celltype.pq`) |
| input | `inputs/aida_pairs.tsv`, 27,857 pairs × 35 columns |
| runner | `run_chunks_shared.sh`: 1,000-pair chunks; each GPU claims the next free chunk, so faster GPUs take more and no chunk runs twice |
| GPUs | 0 and 2 (RTX PRO 6000 Blackwell), 1 (RTX A6000; 35 of 49 GB used) |
| started | 22 November, 12:13 IST; about 30 min per chunk per GPU, 28 chunks |
| finished | 22 November, 17:42 IST; all 27,857 pairs have predictions, no chunk failed |
| output | `results/aida_pairs/outs_<first>_<last>.h5ad` |

## 19.8 Check N: results

Each AIDA eQTL is compared with the Scooby column of its own AIDA Level-3 type (`score_aida.py`, which calls
`evaluate_eqtl.py` as Section 14 did for CIMA). "Clear" means |predicted effect| > 0.05. Two sign metrics are
reported (Section 18.2):

- **paper's metric:** clear predictions whose sign matches beta *and* whose |beta| > 0.05, over all clear predictions;
- **sign, any |beta|:** clear predictions whose sign matches beta, over all clear predictions.

**The three datasets, scored the same way** (`eval_out_aida/datasets_compared.csv`; OneK1K: our run, OneK1K eQTLs
only; CIMA: CIMA cells, own column):

| | OneK1K | CIMA | AIDA |
|:-----------------------------------------------|------:|------:|------:|
| cell types | 11 | 69 | 20 |
| eQTL rows with a prediction | 1,255 | 55,209 | 33,751 |
| median \|beta\| (scored rows) | 0.583 | 0.403 | **0.029** |
| share with \|beta\| > 0.05 | 99.4% | 100% | **25.4%** |
| clear predictions | 286 (22.8%) | 2,863 (5.2%) | 1,131 (3.4%) |
| **paper's metric** | 89.2% (255) | 87.7% (2,512) | **38.4%** (434) |
| **sign, any \|beta\|** | 89.2% (255) | 87.7% (2,512) | **84.6%** (957) |
| sign, rows with \|beta\| > 0.05 only | 89.2% (255 / 286) | 87.7% (2,512 / 2,863) | 85.9% (434 / 505) |
| Spearman(beta, predicted effect), all rows | 0.385 | 0.132 | 0.076 |

For OneK1K and CIMA the two sign metrics are the same number, because nearly every beta is above 0.05. For AIDA they
differ by 46 points. Of the 1,131 clear AIDA predictions, 626 have |beta| < 0.05; the paper's metric marks all of
them wrong, although 523 of them have the right sign. Restricted to the rows where |beta| > 0.05, AIDA scores 85.9%,
close to the other two.

Two real rows (CD14+ monocytes, `eval_out_aida/merged.csv.gz`):

| Variant | Gene | Distance to TSS | beta | p-value | Predicted effect | Paper's metric | Sign, any \|beta\| |
|:--------------------|:--------|------:|------:|------:|------:|:------:|:------:|
| chr13_100588735_C_A | GGACT | 794 bp | +0.049 | 6.5 × 10⁻⁸¹ | +0.481 | wrong | right |
| chr17_50547162_A_C | SPATA20 | 4,103 bp | −0.081 | 5.5 × 10⁻⁵⁷ | −0.546 | right | right |

GGACT's eQTL is extremely significant and Scooby predicts it in the right direction, but its beta (0.049) is just
under 0.05.

**By lineage** (`eval_out_aida/totals_by_lineage.csv`):

| Lineage | Types | eQTL rows | Clear | Paper's metric | Sign, any \|beta\| | Median Spearman |
|:------|---:|------:|------:|------:|------:|------:|
| B | 4 | 3,402 | 99 (2.9%) | 52.5% | 89.9% (89) | 0.077 |
| DC | 1 | 1,609 | 47 (2.9%) | 46.8% | 93.6% (44) | 0.088 |
| Monocytes | 2 | 6,844 | 261 (3.8%) | 29.1% | 85.4% (223) | 0.093 |
| NK | 2 | 3,698 | 99 (2.7%) | 35.4% | 81.8% (81) | 0.121 |
| T | 11 | 18,198 | 625 (3.4%) | 39.8% | 83.2% (520) | 0.074 |
| **All** | **20** | **33,751** | **1,131 (3.4%)** | **38.4%** | **84.6% (957)** | **0.075** |

Per type (`eval_out_aida/per_type.csv`), sign at any |beta| ranges from 76.9% (CD8+ T naive, 78 clear) to 93.6%
(cDC2, 47 clear); the two types with fewer than 5 clear predictions (atypical B 3, CD56+ NK 4) are 100%. Placement
quality does not predict the score: CD8+ T GZMKhi, whose cells are nearest their mapped OneK1K column only 20.8% of
the time (Section 19.5), scores 83.3% (55 of 66).

**By distance to the TSS** (`eval_out_aida/distance_bins_any_beta.csv`):

| Distance | eQTL rows | Clear | Sign, any \|beta\| | Binomial p |
|:------------------|------:|------:|------:|------:|
| ≤ 1 kb | 1,962 | 458 (23.3%) | 81.9% | 2 × 10⁻⁴⁵ |
| 1–2 kb | 1,063 | 84 (7.9%) | 81.0% | 9 × 10⁻⁹ |
| 2–5 kb | 2,415 | 146 (6.0%) | 91.8% | 3 × 10⁻²⁷ |
| 5–10 kb | 3,005 | 166 (5.5%) | 80.1% | 2 × 10⁻¹⁵ |
| 10–20 kb | 4,300 | 94 (2.2%) | 91.5% | 1 × 10⁻¹⁷ |
| 20–50 kb | 7,450 | 124 (1.7%) | 87.9% | 9 × 10⁻¹⁹ |
| 50–100 kb | 5,859 | 49 (0.8%) | 85.7% | 4 × 10⁻⁷ |
| 100–150 kb | 3,296 | 7 (0.2%) | 100% | 0.016 |
| > 150 kb | 4,401 | 3 (0.1%) | 100% | 0.25 |

As for OneK1K and CIMA, Scooby makes clear predictions mostly near the gene, and when it does, the sign is right
at every distance. AIDA has relatively many distant eQTLs (21,006 of 33,751 rows beyond 20 kb), which is one
reason its share of clear predictions (3.4%) is lower than CIMA's (5.2%).

**Cell-type specificity (as check J)** (`eval_out_aida/specificity.csv`). For each type, Spearman(beta, predicted
effect) was computed with its own column and with each of the other 34 AIDA columns on the same eQTLs:

- the own column is the best of 35 for **1 of 20** types (CD16+ monocytes); median rank 20 of 35;
- the own column beats the median column of *other* lineages for 10 of 20 types;
- the 35 columns correlate with each other at median Pearson r = **0.969** over the 27,857 pairs (lowest 0.806; CIMA: 0.958 over 32,607 pairs, lowest 0.714).

**Reading.** On AIDA, with AIDA's own cells, Scooby gets the direction of an eQTL right as often as on OneK1K and
CIMA (84.6% vs 89.2% and 87.7%), once the score does not depend on the unit of beta. The paper's metric (38.4%)
understates this because of the beta cut-off, not because of the predictions. As with CIMA, the predictions are
almost the same in every cell type, so this is evidence that Scooby predicts the direction of eQTL effects across
cohorts, not that it predicts cell-type-specific effects.

## 19.9 Next

1. Check K on AIDA: eQTLs measured in several AIDA types (measured vs predicted differences between types).
2. Optional: thinned embeddings on a sample of pairs, to see whether depth changes the predicted effects.
3. Confirm the CellSubtype FDR definition and the unit of beta in AIDA's Methods.

Caveats:

- Only 3.4% of AIDA rows get a clear prediction; the 84.6% rests on 1,131 rows.
- 5′ chemistry (OneK1K 3′); lead variants, not fine-mapped ones.

\newpage

# 20. Limitations

* The encoder is a **substitute**, not the authors' scPoli. Checks A and B show it reproduces the authors'
  embeddings and Scooby's predictions on OneK1K; on CIMA there is no ground truth, only checks C and D.
* scPoli's pool (batch) correction is not reproduced; check A shows this costs little within OneK1K, but CIMA's
  technical differences (lab, protocol, 31% more depth) are larger than those between OneK1K pools.
* 102 of 2,403 genes (4%) are missing in CIMA and set to 0.
* CIMA's L4 types are finer than OneK1K's 22 columns; types with no OneK1K counterpart (e.g. NKT cells, ILC2)
  are placed in the space Scooby learned from OneK1K, where no such cells existed.
* 60 of 421 donors; six types have fewer than 500 cells.
* The authors' embedding does not separate T-cell subtypes well (key marker genes are Scooby validation/test genes
  or not among the selected variable genes; Section 12), so T-subtype-specific effects are weak, with any cells.
* AIDA (Section 19): 175 of 2,403 genes (7.3%) set to 0; cells 2.5× deeper than OneK1K over the encoder genes;
  10x 5′ chemistry (OneK1K 3′); betas on a smaller scale than CIMA's; lead variants only.

# 21. Files

| File | Contents |
|:------------------------------------------------|:-----------------------------------------------------------|
| `scooby_inference/fit_scpoli_encoder.py` | trains the substitute encoder, runs check A |
| `scooby_inference/scpoli_substitute/` | `encoder.pt`, `genes.tsv`, `heldout_summary.json`, `heldout_per_column.csv`, `embedding_heldout_predicted.pq` |
| `scooby_inference/step2_config_substitute.yaml` | Scooby config with the substitute embeddings (check B) |
| `scooby_inference/run_all_chunks.sh` | now takes a config and output folder (3rd, 4th arguments) |
| `results/gtex_onek1k_pairs_substitute/` | Scooby predictions with the substitute embeddings |
| `eval_onek1k/eval_out_substitute_vs_original/`, `eval_onek1k/eval_out_substitute/` | check B |
| `scooby_inference/embed_cima_cells.py` | encodes CIMA cells, checks C and D |
| `data/cima_cells/embedded_B/`, `data/cima_cells/embedded_B_thin076/` | checks C and D |
| `scooby_inference/fetch_cima_donors.py`, `remote_h5.py` | partial download of CIMA donors, check E |
| `data/cima_cells/test_B_2donors.h5ad`, `plan_60donors.log`, `fetch_60donors.log` | check E, plan, first download attempt |
| `data/cima_cells/CIMA_60donors_seed0.h5ad`, `build_60donors.log` | the 60-donor file (check F) |
| `data/cima_cells/embedded_60donors/` | `embedding.pq`, `celltype.pq`, `placement_per_type.csv`, `placement_lineage.csv`, `onek1k_baseline_per_column.csv` (check G) |
| `scooby_inference/step2_config_cima_cells.yaml` | Scooby config with CIMA cells |
| `results/cima_cells_test/`, `results/cima_pairs_cima_cells/` | test run (check H) and the full run |
| `scooby_inference/score_cima_cells.py` | scores the full run in own types, compares with the 48-type run, check J |
| `eval_out_cima_cells/` | `summary.csv`, `merged.csv.gz`, `per_type_vs_onek1k_cells.csv`, `totals_by_lineage.csv`, `same_pairs.csv`, `specificity.csv`, `compare_ontology.csv`, `compare_mapping48.csv` (checks I, J), `clear_by_beta_quarter.csv` (Section 16) |
| `scooby_inference/` `check_k_celltype_differences.py` | check K (Section 17): measured vs predicted differences between cell types |
| `eval_out_cima_cells/check_k_cases.csv`, `check_k_filters.csv`, `check_k_large_predicted.csv`, `check_k.log` | check K results |
| `scooby_inference/check_k_onek1k.py`, `eval_onek1k/check_k/`, `eval_onek1k/check_k_substitute/`, `eval_onek1k/check_k_authors.log` | check K on OneK1K (Section 17.7) |
| `eval_out_cima_cells/levels_6pairs/` | 6 eQTLs with predicted REF/ALT expression per column (`levels_table.csv`), the script copy and its diff, `README.txt` |
| `scooby_inference/make_aida_inputs.py` | AIDA Table S6 → `data/aida_eqtl.tsv`, `inputs/aida_pairs.tsv`, `data/aida_input_funnel.csv` (Section 19.3) |
| `scooby_inference/embed_aida_cells.py`, `aida_mapping.csv` | encodes AIDA cells, placement and depth checks (Sections 19.4–19.6) |
| `data/aida/embedded/`, `data/aida/embedded_thin040/` | AIDA `embedding.pq`, `celltype.pq`, `genes_matched.tsv`, `placement_per_type.csv`, `summary.json`; thinned run with `moved_per_type.csv` |
| `scooby_inference/step2_config_aida.yaml`, `run_chunks_shared.sh` | AIDA Scooby config; multi-GPU chunk runner |
| `results/aida_pairs/` | AIDA Scooby run (finished 22 November) |
| `scooby_inference/score_aida.py`, `eval_out_aida.log` | scores the AIDA run (Section 19.8, check N) |
| `eval_out_aida/` | `summary.csv`, `merged.csv.gz`, `per_type.csv`, `totals_by_lineage.csv`, `datasets_compared.csv`, `distance_bins_any_beta.csv`, `specificity.csv` |
| `figures/make_scooby_schematic.py`, `figures/scooby_schematic.pdf/.png` | Figure 1 (Section 1.1) |
| `requirements-step2.txt` | NumPy pinned to 1.26.4 |
| `fit_scpoli_encoder_explained.pdf` (repo root) | step-by-step explanation of the encoder script |
