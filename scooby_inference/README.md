# scooby step 2: run the OneK1K model on the Blackwell GPU and check it against the authors' numbers

Put this folder next to the `scooby_assets/` folder from step 1 (for example `~/scooby/scooby_step2/`), then:

```bash
cd ~/scooby/scooby_step2
nvidia-smi                                       # pick a free Blackwell card; GPU 0 below
bash setup_env.sh                                # once, ~10 min
source ~/venvs/scooby-inf/bin/activate

CUDA_VISIBLE_DEVICES=0 python smoke_test.py      # GPU + model check, ~5 min (downloads the 765 MB weights once)
python fetch_step2_files.py                      # ~1.1 GB into ./data (+3 GB unzipped genome)
python make_chunk_input.py --n 200               # first 200 rows of the authors' chunk 0
# run the command it prints, e.g.:
CUDA_VISIBLE_DEVICES=0 python scooby_reproducibility/scripts/inference_eqtls/ScoobyVariant.py \
    -in_path inputs/chunk0.tsv -out_path results/outs_0_199.h5ad -config_path step2_config.yaml
python compare_with_authors.py --ours results/outs_0_199.h5ad
```

What each step checks:

| Script | Checks | Pass looks like |
|---|---|---|
| `smoke_test.py` | PyTorch has Blackwell kernels; the released model loads and runs | `sm_120` listed, 191.2 M parameters, decoder input (128, 10), three timed runs |
| `fetch_step2_files.py` | the needed Zenodo files exist and download | 11 prediction chunks, genome unzipped |
| `make_chunk_input.py` | the authors' scored rows can be rebuilt | ~2,100 rows in chunk 0, columns include Chromosome, gene_name, Pos, Ref, Alt |
| `ScoobyVariant.py` (authors' code) | your run on the same rows | finishes without the REF assertion error |
| `compare_with_authors.py` | your numbers vs theirs | `PASS`: Pearson r > 0.99 and > 95% of values within 0.01 |

Notes:
- Long runs: start them inside `tmux` so they survive a dropped SSH connection.
- Run time: about 2 model passes per row and allele; the smoke test's timing × 4 × rows gives an estimate.
- If the comparison says CHECK, keep the output and the printed largest differences; they show whether the
  problem is a few rows (position or strand handling) or everything (environment or precision).
- Effects are natural-log fold changes (as in the authors' script), not log2.
