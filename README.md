# Scooby eQTL inference: OneK1K and CIMA

This repository runs the released OneK1K Scooby model (`lauradmartens/onek1k-scooby`) on variant-gene pairs and scores its predicted effects against measured eQTLs:
- **OneK1K / GTEx whole blood:** a reproduction of the paper's Fig. 4.
- **CIMA:** lead cis-eQTLs from an independent cohort.

For each pair, the prediction is `ln(ALT+1) − ln(REF+1)` of the expression predicted for each of the 22 OneK1K cell types, summed over that type's cells (natural log).

Layout: `scooby_inference/` holds the model run and scoring; `cima/` holds the CIMA eQTL table and the script that turns it into Scooby inputs.

All commands are run from `scooby_inference/` (`cd scooby_inference`) unless stated otherwise.

## 1. Environment (once)
Requires a CUDA GPU. Tested on an RTX Blackwell (sm_120) card.
```bash
bash setup_env.sh                   # creates ~/venvs/scooby-inf (Python 3.11), installs torch 2.7.1 (CUDA 12.8),
                                    # requirements-step2.txt, scooby (--no-deps), and clones scooby_reproducibility/
source ~/venvs/scooby-inf/bin/activate
```
The model weights are downloaded from Hugging Face on the first run (~765 MB).

## 2. Data (once)
These files go into `scooby_inference/data/`. Files marked "(in git)" come with the repo; download the rest from the authors' Zenodo records:

| File in `scooby_inference/data/` | Source |
|---|---|
| `embedding_no_val_genes.pq`, `celltype_fixed.pq` (in git) | Zenodo 15517072, `training_data.zip` → `onek1k_training_data/` |
| `genome_human.fa` | same zip → `scooby_training_data/genome_human.fa.gz` (unzip it) |
| `gencode.v32.annotation.sorted.gtf.gz` (+ `.tbi`) | Zenodo 15517764, `reproducibility_data.zip` |
| `all_highpip_susie.tsv` (in git) | same zip → `eQTL_fig4/` |
| `authors_outs/` (11 `.h5ad`, in git) | same zip → `eQTL_fig4/scooby_variants_tsscentred/` (the authors' predictions, only needed for the comparison) |
| `cima_eqtl.tsv` (in git) | built by `cima/prepare_cima.py` (see 4a) |

`step2_config.yaml` contains **absolute paths**. Edit them if your checkout is somewhere other than `/home/shreyas/scooby`.

## 3. OneK1K / GTEx whole blood
```bash
# a. predict (GPU; run inside tmux). Writes results/gtex_onek1k_pairs/outs_*.h5ad
bash run_all_chunks.sh 0 inputs/gtex_onek1k_pairs.tsv        # 1,585 pairs; 0 = GPU id

# b. score against measured eQTLs
python evaluate_eqtl.py --pred "results/gtex_onek1k_pairs/outs_*.h5ad" --out eval_onek1k/eval_out_ours
python evaluate_eqtl.py                                      # authors' predictions -> eval_onek1k/eval_out_authors_eqtl

# c. check that our predictions reproduce the authors'
python compare_onek1k_with_authors.py                        # -> eval_onek1k/eval_out_compare
```
`inputs/gtex_onek1k_pairs.tsv` holds the GTEx-blood and OneK1K pairs the authors scored, rebuilt from the `.obs` of their output files.

## 4. CIMA
```bash
# a. build the inputs (already in git; only rerun if needed). Needs data/genome_human.fa and the GTF
cd ../cima && python prepare_cima.py && cd ../scooby_inference
#    xQTL/CIMA_Lead_cis-xQTL.csv -> inputs/cima_pairs.tsv (36,739 pairs) and data/cima_eqtl.tsv
#    keeps SNV cis-eQTLs whose REF matches hg38; beta = CIMA slope

# b. predict (GPU, several hours; rerunning skips finished chunks)
bash run_all_chunks.sh 0 inputs/cima_pairs.tsv               # -> results/cima_pairs/

# c. score (run step 3b first: it adds the OneK1K scores for comparison)
python score_cima.py                                         # -> eval_out_cima/
```
CIMA cell types are matched to Scooby's 22 OneK1K columns through `cima_mapping.csv`. The predictions use **OneK1K cell embeddings**; only the variants and measured betas come from CIMA.

## 5. Outputs
| Folder | Main files |
|---|---|
| `eval_onek1k/eval_out_ours/`, `eval_onek1k/eval_out_authors_eqtl/` | `summary.csv` (n, Spearman, sign agreement for \|pred\| > 0.05, per study/cell type), `distance_bins.csv`, `merged.csv.gz` |
| `eval_onek1k/eval_out_compare/` | `values_per_celltype.csv` (ours vs authors), `scores_per_celltype.csv`, `pairs.csv.gz` |
| `eval_out_cima/` | the same three files, plus `cima_scores.csv` (per CIMA cell type, with OneK1K reference), `group_totals.csv`, `distance_main.csv` |

`results/` (the raw prediction `.h5ad` files) is git-ignored.
