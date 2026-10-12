"""Step 2b: check that PyTorch sees the Blackwell GPU and that the released OneK1K model runs on it.

Run:  CUDA_VISIBLE_DEVICES=0 python smoke_test.py
Uses a random DNA sequence and random cell embeddings, so it needs no data files.
"""
import time

import torch
from scooby.modeling import Scooby
from scooby.utils.utils import fix_rev_comp_rna, get_pseudobulk_count_pred

# 1. GPU and PyTorch
assert torch.cuda.is_available(), "PyTorch does not see a GPU"
print("torch", torch.__version__, "| CUDA", torch.version.cuda, "|", torch.cuda.get_device_name(0))
arch = torch.cuda.get_arch_list()
print("kernels built for:", arch)
assert "sm_120" in arch, "this PyTorch build has no Blackwell (sm_120) kernels; reinstall from the cu128 index"
x = torch.randn(4096, 4096, device="cuda", dtype=torch.bfloat16)
print("matmul OK:", float((x @ x).float().abs().mean()) > 0)

# 2. Model, loaded exactly as in the authors' ScoobyVariant.py
model = Scooby.from_pretrained(
    "lauradmartens/onek1k-scooby", cell_emb_dim=10, embedding_dim=1920, n_tracks=2,
    return_center_bins_only=True, disable_cache=False, use_transform_borzoi_emb=True,
).eval().to("cuda")
print("parameters (M):", round(sum(p.numel() for p in model.parameters()) / 1e6, 1), "(expected 191.2)")
print("decoder input:", tuple(model.cell_state_to_conv[0].weight.shape), "(expected (128, 10))")


# 3. The authors' predict(): forward + reverse-complement pass, averaged (fp16 autocast, as they ran it)
def predict(model, seqs, seqs_rev_comp, conv_weights, conv_biases, bins_to_predict=None):
    with torch.no_grad(), torch.autocast("cuda"):
        out = model.forward_sequence_w_convs(seqs, conv_weights, conv_biases, bins_to_predict=bins_to_predict)
        out_rc = model.forward_sequence_w_convs(seqs_rev_comp, conv_weights, conv_biases,
                                                bins_to_predict=(6143 - bins_to_predict))
    return (out + fix_rev_comp_rna(torch.flip(out_rc, (1, -3)))) / 2


seq = torch.nn.functional.one_hot(torch.randint(0, 4, (524288,)), 4).float()
seq = seq.unsqueeze(0).permute(0, 2, 1).cuda()          # (1, 4, 524288), as in the script
cells = torch.randn(1, 50_000, 10, device="cuda")        # 50k random "cells"
with torch.no_grad():
    w, b = model.forward_cell_embs_only(cells)
weights = [(w.half(), b.half())]
bins = torch.arange(3000, 3040, device="cuda")           # 40 output bins, like a gene's exons

for i in range(3):  # first call includes CUDA warm-up
    torch.cuda.synchronize(); t = time.time()
    total = get_pseudobulk_count_pred(model, seq, cell_emb_conv_weights_and_biases=weights, gene_slice=bins,
                                      strand="+", predict=predict, clip_soft=5, model_type="rna")
    torch.cuda.synchronize()
    print(f"run {i}: summed prediction {float(total[0]):.2f}, {time.time() - t:.2f} s")
    model.sequences, model.last_embs = [], []            # clear the cache so each run is a full pass

print(f"peak GPU memory: {torch.cuda.max_memory_allocated() / 1e9:.1f} GB")
print("Smoke test passed. Next: python fetch_step2_files.py")
