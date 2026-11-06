#!/usr/bin/env bash
# Run the authors' ScoobyVariant.py on all rows of an input table, in chunks of 2,000 pairs.
# Each chunk is saved as results/<input name>/outs_<first>_<last>.h5ad; rerunning skips chunks already done.
# Run inside tmux:  bash run_all_chunks.sh [gpu] [input]
#   bash run_all_chunks.sh 0 inputs/gtex_onek1k_pairs.tsv     # GTEx blood + OneK1K pairs (1,585)
#   bash run_all_chunks.sh 0 inputs/all_pairs.tsv             # all pairs the authors scored (19,613)
set -uo pipefail
cd "$(dirname "$0")"
source ~/venvs/scooby-inf/bin/activate

GPU=${1:-0}
IN=${2:-inputs/gtex_onek1k_pairs.tsv}
OUT=results/$(basename "$IN" .tsv)
mkdir -p "$OUT"
LAST=$(( $(wc -l < "$IN") - 2 ))   # last row index (file has a header line)
STEP=2000

for s in $(seq 0 $STEP $LAST); do
  e=$(( s + STEP - 1 < LAST ? s + STEP - 1 : LAST ))
  out=$OUT/outs_${s}_${e}.h5ad
  if [ -f "$out" ]; then echo "skip $out (done)"; continue; fi
  echo "$(date '+%F %T') start rows $s-$e of $IN on GPU $GPU"
  CUDA_VISIBLE_DEVICES=$GPU python scooby_reproducibility/scripts/inference_eqtls/ScoobyVariant.py -in_path "$IN" -out_path "$out" -config_path step2_config.yaml 2>&1 | tee "$OUT/run_${s}_${e}.log"
  if [ "${PIPESTATUS[0]}" -ne 0 ] || [ ! -f "$out" ]; then
    echo "$(date '+%F %T') FAILED rows $s-$e; see $OUT/run_${s}_${e}.log. Stopping."
    exit 1
  fi
  echo "$(date '+%F %T') done rows $s-$e"
done
echo "$(date '+%F %T') all chunks finished; results in $OUT/"
