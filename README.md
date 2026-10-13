```
cd ~/scooby/scooby_inference
nvidia-smi                                      
bash setup_env.sh                               
source ~/venvs/scooby-inf/bin/activate

CUDA_VISIBLE_DEVICES=0 python smoke_test.py      
python fetch_step2_files.py                      # ~1.1 GB into ./data (+3 GB unzipped genome)
python make_chunk_input.py --n 200               # first 200 rows of the authors' chunk 0
# run the command it prints, e.g.:
CUDA_VISIBLE_DEVICES=0 python scooby_reproducibility/scripts/inference_eqtls/ScoobyVariant.py \
    -in_path inputs/chunk0.tsv -out_path results/outs_0_199.h5ad -config_path step2_config.yaml
python compare_with_authors.py --ours results/outs_0_199.h5ad
```
