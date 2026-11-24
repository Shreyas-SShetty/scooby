#!/usr/bin/env bash
# Like run_all_chunks.sh, but several GPUs can work through the same input at once: each chunk is claimed with an
# atomic mkdir (OUT/claims/<s>_<e>) before it is run, so no chunk runs twice and faster GPUs take more chunks.
# Output rows keep the row numbers of the single input file, as with run_all_chunks.sh.
# Run one per GPU, each in its own tmux window:
#   bash run_chunks_shared.sh 0 inputs/aida_pairs.tsv step2_config_aida.yaml results/aida_pairs
#   bash run_chunks_shared.sh 1 ...   bash run_chunks_shared.sh 2 ...
# A failed chunk leaves its claim and a FAILED file; delete OUT/claims/<s>_<e> to let it run again.
set -uo pipefail
cd "$(dirname "$0")"
source ~/venvs/scooby-inf/bin/activate
GPU=$1; IN=$2; CONFIG=$3; OUT=$4; STEP=${5:-1000}
mkdir -p "$OUT/claims"
LAST=$(( $(wc -l < "$IN") - 2 ))
for s in $(seq 0 $STEP $LAST); do
  e=$(( s + STEP - 1 < LAST ? s + STEP - 1 : LAST ))
  out=$OUT/outs_${s}_${e}.h5ad
  [ -f "$out" ] && continue
  mkdir "$OUT/claims/${s}_${e}" 2>/dev/null || continue   # another GPU has it
  echo "$(date '+%F %T') GPU $GPU start rows $s-$e"
  CUDA_VISIBLE_DEVICES=$GPU python scooby_reproducibility/scripts/inference_eqtls/ScoobyVariant.py -in_path "$IN" -out_path "$out" -config_path "$CONFIG" > "$OUT/run_${s}_${e}.log" 2>&1
  if [ $? -ne 0 ] || [ ! -f "$out" ]; then
    echo "$(date '+%F %T') GPU $GPU FAILED rows $s-$e (see $OUT/run_${s}_${e}.log)" | tee "$OUT/FAILED_${s}_${e}"
    continue
  fi
  echo "$(date '+%F %T') GPU $GPU done rows $s-$e"
done
echo "$(date '+%F %T') GPU $GPU: no unclaimed chunks left"
