# Scooby eQTL inference: OneK1K, GTEx, CIMA and AIDA

Runs the released OneK1K Scooby model (`lauradmartens/onek1k-scooby`) on variant-gene pairs and scores the predicted
effect, `ln(ALT+1) − ln(REF+1)` of the expression summed over a cell type's cells, against measured eQTLs.

All commands are run from `scooby_inference/` (`cd scooby_inference`) unless stated otherwise. Long GPU runs belong
in tmux; `run_all_chunks.sh` / `run_chunks_shared.sh` skip chunks that are already done, so a rerun continues.
The `step2_config*.yaml` files contain absolute paths: edit them if the checkout is not at `/home/shreyas/scooby`.

## Setup (once)
```bash
bash setup_env.sh                      # ~/venvs/scooby-inf (Python 3.11, torch 2.7.1 / CUDA 12.8) + clones scooby_reproducibility/
source ~/venvs/scooby-inf/bin/activate
```
Data in `scooby_inference/data/` (from the authors' Zenodo records 15517072 and 15517764 unless in git):
`embedding_no_val_genes.pq`, `celltype_fixed.pq`, `all_highpip_susie.tsv` (in git);
`genome_human.fa`, `gencode.v32.annotation.sorted.gtf.gz` (+ `.tbi`); `scpoli/training_data/` (OneK1K counts, for the encoder).

## Other experiments 

**1. GTEx whole blood and OneK1K eQTLs with OneK1K cells: reproduce the paper's Fig. 4.**
```bash
bash run_all_chunks.sh 0 inputs/gtex_onek1k_pairs.tsv                       # -> results/gtex_onek1k_pairs/
python evaluate_eqtl.py --pred "results/gtex_onek1k_pairs/outs_*.h5ad" --out eval_onek1k/eval_out_ours
```

**2. CIMA eQTLs with OneK1K cells: CIMA cell types mapped to OneK1K columns through the Cell Ontology.**
```bash
cd ../cima && python prepare_cima.py && cd ../scooby_inference               # -> inputs/cima_pairs.tsv, data/cima_eqtl.tsv
bash run_all_chunks.sh 0 inputs/cima_pairs.tsv                               # -> results/cima_pairs/
python map_cima_celltypes.py --obo cl-basic.obo                              # -> cima_mapping_ontology.csv
python evaluate_eqtl.py --pred "results/cima_pairs/outs_*.h5ad" --eqtl data/cima_eqtl.tsv \
    --mapping cima_mapping_ontology.csv --out eval_out_cima_ontology
```

**3. Substitute scPoli encoder: learn the authors' (unreleased) encoder from OneK1K counts → embeddings.**
```bash
python fit_scpoli_encoder.py                                                 # -> scpoli_substitute/encoder.pt
```

**4. CIMA eQTLs with CIMA cells: download 60 CIMA donors, encode them (incl. the lineage check), score each eQTL in its own CIMA type.**
```bash
python fetch_cima_donors.py --n_donors 60                                    # partial download over HTTP -> data/cima_cells/CIMA_60donors_seed0.h5ad
python embed_cima_cells.py --cima data/cima_cells/CIMA_60donors_seed0.h5ad --out data/cima_cells/embedded_60donors
bash run_all_chunks.sh 0 inputs/cima_pairs.tsv step2_config_cima_cells.yaml results/cima_pairs_cima_cells
python score_cima_cells.py                                                   # -> eval_out_cima_cells/
```

**5. Does the cell type change the variant effect? REF and ALT predictions saved separately for all 755 OneK1K fine-mapped eQTLs.**
```bash
D=eval_onek1k/levels_755; i=0
for r in 0_251 252_503 504_754; do
  CUDA_VISIBLE_DEVICES=$i python $D/ScoobyVariant_levels.py -in_path $D/pairs.tsv -out_path $D/outs_$r.h5ad \
      -config_path step2_config.yaml & i=$((i+1))
done; wait                                                                   # -> outs_*.h5ad with layers log_ref / log_alt
```

**6. AIDA eQTLs with AIDA cells: AIDA Table S6 lead eQTLs, AIDA cells (CELLxGENE) encoded with the substitute encoder.**
```bash
python make_aida_inputs.py                                                   # -> data/aida_eqtl.tsv, inputs/aida_pairs.tsv
python embed_aida_cells.py                                                   # -> data/aida/embedded/
bash run_chunks_shared.sh 0 inputs/aida_pairs.tsv step2_config_aida.yaml results/aida_pairs   # one per GPU: 0, 1, 2
python score_aida.py                                                         # -> eval_out_aida/
```

Not in git: `results/` (raw predictions), raw cell files (`data/cima_cells/`, `data/aida/cells/`, `data/scpoli/`),
the AIDA paper's supplementary files, large embedding files, PDFs, decks and `figures/`.
