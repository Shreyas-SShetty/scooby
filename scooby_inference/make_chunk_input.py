"""Step 2d: rebuild the authors' input rows for chunk 0 and write a config for their script.

The authors' input table (all_highpip_vars.tsv) is not on Zenodo, but their script stores the rows it
scored in each output file's .obs, together with 3 columns it added from the gene-window file
(Start, End, Strand). Dropping those 3 gives back the input rows, in the authors' own position convention.

Run:  python make_chunk_input.py --n 200      # quick test: first 200 rows (~5-10 min on the GPU)
      python make_chunk_input.py              # whole chunk (~2,100 rows, about 1 h)
"""
import argparse
from pathlib import Path

import anndata as ad
import yaml

p = argparse.ArgumentParser()
p.add_argument("--n", type=int, default=None, help="use only the first n rows")
p.add_argument("--chunk", default="data/authors_outs/outs_0_2129.h5ad")
args = p.parse_args()

theirs = ad.read_h5ad(args.chunk)
obs = theirs.obs.drop(columns=["Start", "End", "Strand"])
for col in obs.columns:  # categorical columns from AnnData -> plain values
    if str(obs[col].dtype) == "category":
        obs[col] = obs[col].astype(str)
if args.n:
    obs = obs.iloc[: args.n]
print(f"{len(obs)} rows; columns: {list(obs.columns)}")
print(f"authors' output: {theirs.n_obs} rows x {theirs.n_vars} cell types")

Path("inputs").mkdir(exist_ok=True)
Path("results").mkdir(exist_ok=True)
in_path = Path("inputs/chunk0.tsv").resolve()
obs.to_csv(in_path, sep="\t", index=False)

data = Path("data").resolve()
config = {
    "pretrained_path": "lauradmartens/onek1k-scooby",
    "embedding_path": str(data / "embedding_no_val_genes.pq"),
    "cell_type_idx_path": str(data / "celltype_fixed.pq"),
    "fasta_file": str(data / "genome_human.fa"),
    "gtf_file": str(data / "gencode.v32.annotation.sorted.gtf.gz"),
    "bed_file": str(Path("scooby_reproducibility/scripts/inference_eqtls/genes_tsscentred.bed").resolve()),
    "data_path": str(in_path),
    "span": True,
    "return_center_bins_only": True,
    "replicate": None,
    "disable_autocast": False,  # authors ran with fp16 autocast
}
Path("step2_config.yaml").write_text(yaml.safe_dump(config, sort_keys=False))

# The script takes the row range from the output file name: outs_<first>_<last>.h5ad
out_path = f"results/outs_0_{len(obs) - 1}.h5ad"
print("\nRun the authors' script on GPU 0:")
print(f"  CUDA_VISIBLE_DEVICES=0 python scooby_reproducibility/scripts/inference_eqtls/ScoobyVariant.py "
      f"-in_path {in_path} -out_path {out_path} -config_path step2_config.yaml")
print(f"Then compare:\n  python compare_with_authors.py --ours {out_path}")
