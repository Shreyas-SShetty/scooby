"""Download the raw counts of all cells from a random set of CIMA donors, without downloading the whole 44 GB file.

The full CIMA file (CIMA_RNA_6484974cells_36326genes_compressed.h5ad on ftp.cngb.org) stores cells grouped by
library, and each donor's cells sit in one or a few contiguous blocks. So all cells of a set of donors can be read
with HTTP range requests: only layers/counts (raw counts) for those blocks, not the normalised X.
Taking whole donors keeps every cell type in its natural proportion.

Output: an .h5ad with layers/counts (raw counts), obs = cell_type_l4, sample, library, file_row (the cell's row
in the full file; cell names are 'row<file_row>'), and var = CIMA gene symbols. embed_cima_cells.py reads it as is.

Run (from scooby_inference/):
  python fetch_cima_donors.py --n_donors 60 --plan     # which donors, cells per type, estimated size and time
  python fetch_cima_donors.py --n_donors 60            # download (inside tmux; rerun continues after a stop)
  python fetch_cima_donors.py --n_donors 60 --download_only   # download only; build later, e.g. on another machine:
  python fetch_cima_donors.py --n_donors 60 --build_offline   # build the .h5ad from <out>.chunks/ (no network)
The build needs ~15-25 GB of RAM for 60 donors; --download_only + --build_offline lets a small machine download
and a big one build (copy the whole <out>.chunks/ folder, which includes meta.npz).
"""
import argparse
import os
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from remote_h5 import (chunk_bytes_per_element, dset_info, open_remote, read_cached, read_ranges,  # noqa: E402
                       read_rows)

URL = ("https://ftp.cngb.org/pub/SciRAID/trueblood/cima/CIMA_Resource/Cell_Atlas/"
       "CIMA_RNA_6484974cells_36326genes_compressed.h5ad")
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--url", default=URL)
    p.add_argument("--n_donors", type=int, default=60)
    p.add_argument("--donors", default="", help="comma-separated sample IDs instead of a random choice")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", default="")
    p.add_argument("--workers", type=int, default=8, help="parallel connections (the server blocks clients that ask too much)")
    p.add_argument("--plan", action="store_true", help="only show what would be downloaded")
    p.add_argument("--download_only", action="store_true", help="download into <out>.chunks/, do not build")
    p.add_argument("--build_offline", action="store_true", help="build from <out>.chunks/ only (no network)")
    args = p.parse_args()
    out = args.out or f"data/cima_cells/CIMA_{args.donors.count(',') + 1 if args.donors else args.n_donors}donors_seed{args.seed}.h5ad"
    cache = out + ".chunks"
    t0 = time.time()
    log = lambda *a: print(f"[{(time.time() - t0) / 60:5.1f} min]", *a, flush=True)
    dec = lambda xs: [x.decode() if isinstance(x, bytes) else x for x in xs]

    if args.build_offline:
        build(out, cache, log)
        return

    f, h = open_remote(args.url)
    log(f"opened {args.url.split('/')[-1]} ({f.raw.size / 1e9:.1f} GB)")
    o, C = h["obs"], h["layers"]["counts"]
    cat = {}
    for k in ["sample", "cell_type_l4", "library"]:
        cat[k] = np.array(dec(o[k]["categories"][:]), dtype=object)
        cat[k + "_codes"] = read_rows(o[k]["codes"], args.url, 0, o[k]["codes"].shape[0], workers=args.workers,
                                      cache_dir=cache, label=f"obs {k}")
    log(f"obs read: {len(cat['sample_codes']):,} cells, {len(cat['sample'])} donors, {len(cat['cell_type_l4'])} types")

    # ---------- choose donors ----------
    if args.donors:
        chosen = [cat["sample"].tolist().index(d) for d in args.donors.split(",")]
    else:
        present = np.unique(cat["sample_codes"][cat["sample_codes"] >= 0])
        chosen = sorted(np.random.default_rng(args.seed).choice(present, args.n_donors, replace=False).tolist())
    rows = np.where(np.isin(cat["sample_codes"], chosen))[0]
    brk = np.where(np.diff(rows) != 1)[0]
    runs = list(zip(np.r_[rows[0], rows[brk + 1]], np.r_[rows[brk], rows[-1]] + 1))  # [start, end) in file rows
    log(f"{len(chosen)} donors: {len(rows):,} cells in {len(runs)} contiguous blocks")

    types = pd.Series(cat["cell_type_l4"][cat["cell_type_l4_codes"]])
    per_type = pd.DataFrame({"all_donors": types.value_counts(), "chosen_donors": types.iloc[rows].value_counts()}) \
        .fillna(0).astype(int).sort_values("all_donors")
    pd.set_option("display.max_rows", 200)
    print(per_type.to_string())

    ip = read_ranges(C["indptr"], args.url, [(a, b + 1) for a, b in runs], workers=args.workers, cache_dir=cache,
                     label="indptr")
    n_el = sum(int(x[-1] - x[0]) for x in ip)
    e0 = int(ip[0][0])  # estimate from the blocks at the first chosen donor
    bpe = chunk_bytes_per_element(C["data"], e0) + chunk_bytes_per_element(C["indices"], e0)
    log(f"nonzero values to fetch: {n_el:,}; estimated download {n_el * bpe / 1e9:.2f} GB "
        f"(about {n_el * bpe / 0.45e6 / 3600:.1f} h at 0.45 MB/s)")
    if args.plan:
        sys.exit(0)

    # everything the build needs besides the chunk files, so it can run offline (plain arrays, no pickle)
    meta = {"rows": rows, "indptr": np.concatenate([[0], np.cumsum(np.concatenate([np.diff(x) for x in ip]))]),
            "el_ranges": np.array([(int(x[0]), int(x[-1])) for x in ip], dtype=np.int64),
            "n_genes": np.array(int(C.attrs["shape"][1])),
            "genes": np.array(dec(h["var"][h["var"].attrs["_index"]][:]), dtype=str)}
    for k in ["sample", "cell_type_l4", "library"]:
        meta[k] = np.array(cat[k], dtype=str)
        meta[k + "_codes"] = cat[k + "_codes"][rows]
    for k in ["data", "indices"]:
        i = dset_info(C[k])
        meta[k + "_name"], meta[k + "_chunk"] = np.array(i["name"]), np.array(i["chunk"])
        meta[k + "_filters"], meta[k + "_dtype"] = np.array(i["filters"], dtype=np.int64), np.array(i["dtype"])
    np.savez(os.path.join(cache, "meta.npz"), **meta)

    # ---------- download counts ----------
    el_ranges = [(int(x[0]), int(x[-1])) for x in ip]
    log("downloading layers/counts/data (count values)")
    read_ranges(C["data"], args.url, el_ranges, workers=args.workers, cache_dir=cache, label="data", load=False)
    log("downloading layers/counts/indices (gene of each value)")
    read_ranges(C["indices"], args.url, el_ranges, workers=args.workers, cache_dir=cache, label="indices", load=False)
    log("all chunks saved")

    if args.download_only:
        log(f"done; copy the whole {cache}/ folder to the machine that builds, then run with --build_offline")
        return
    build(out, cache, log)


def build(out, cache, log):
    m = dict(np.load(os.path.join(cache, "meta.npz")))
    el_ranges = [tuple(map(int, r)) for r in m["el_ranges"]]
    info = {k: {"name": str(m[k + "_name"]), "chunk": int(m[k + "_chunk"]),
                "filters": m[k + "_filters"].tolist(), "dtype": str(m[k + "_dtype"])} for k in ["data", "indices"]}
    log("reading count values from saved chunks")
    data = read_cached(info["data"], cache, el_ranges, np.float32)
    log("reading gene positions from saved chunks")
    ind = read_cached(info["indices"], cache, el_ranges, np.int32)
    n_genes, rows = int(m["n_genes"]), m["rows"]
    assert ind.min() >= 0 and ind.max() < n_genes, "gene position out of range: chunks damaged?"
    X = sp.csr_matrix((data, ind, m["indptr"]), shape=(len(rows), n_genes))
    obs = pd.DataFrame({
        "cell_type_l4": pd.Categorical.from_codes(m["cell_type_l4_codes"], m["cell_type_l4"]),
        "sample": pd.Categorical.from_codes(m["sample_codes"], m["sample"]),
        "library": pd.Categorical.from_codes(m["library_codes"], m["library"]),
        "file_row": rows}, index=[f"row{r}" for r in rows])
    a = ad.AnnData(obs=obs, var=pd.DataFrame(index=m["genes"]), layers={"counts": X})
    log(f"writing {out}")
    a.write_h5ad(out, compression="gzip")
    log(f"written {out}: {a.n_obs:,} cells x {a.n_vars:,} genes, {X.nnz:,} nonzero counts")
    log(f"the downloaded chunks in {cache}/ can be deleted now")


# Run only when called as a script (so the functions can be imported, e.g. by worker processes).
if __name__ == "__main__":
    main()
