"""Step 2c: fetch only the files needed to rerun the authors' variant scoring (~1.1 GB, not 78 GB).

Run:  python fetch_step2_files.py
If Zenodo refuses partial reads, download the zips with wget and pass them:
      python fetch_step2_files.py --training-zip training_data.zip --repro-zip reproducibility_data.zip
Everything goes into ./data. Files already present (e.g. from step 1 in ../scooby_assets) are reused.
"""
import argparse
import gzip
import shutil
import zipfile
from pathlib import Path

from remotezip import RemoteZip

p = argparse.ArgumentParser()
p.add_argument("--training-zip", default="https://zenodo.org/records/15517072/files/training_data.zip?download=1")
p.add_argument("--repro-zip", default="https://zenodo.org/records/15517764/files/reproducibility_data.zip?download=1")
p.add_argument("--step1-dir", default="../scooby_assets", help="where step 1 saved the two OneK1K files")
args = p.parse_args()

DATA = Path("data")
DATA.mkdir(exist_ok=True)

# (zip, path ending inside the zip, local name). Path endings come from the Zenodo file listings.
WANTED = [
    ("training", "onek1k_training_data/embedding_no_val_genes.pq", "embedding_no_val_genes.pq"),
    ("training", "onek1k_training_data/celltype_fixed.pq", "celltype_fixed.pq"),
    ("training", "scooby_training_data/genome_human.fa.gz", "genome_human.fa.gz"),
    ("repro", "gencode.v32.annotation.sorted.gtf.gz", "gencode.v32.annotation.sorted.gtf.gz"),
    ("repro", "gencode.v32.annotation.sorted.gtf.gz.tbi", "gencode.v32.annotation.sorted.gtf.gz.tbi"),
    ("repro", "eQTL_fig4/all_highpip_susie.tsv", "all_highpip_susie.tsv"),
]
# The authors' predictions: all 11 chunks (~5 MB in total) go into data/authors_outs/.
AUTHORS_OUTS = "eQTL_fig4/scooby_variants_tsscentred/outs_"


def open_zip(source):
    if Path(source).exists():
        return zipfile.ZipFile(source)
    try:
        return RemoteZip(source)
    except Exception:
        return RemoteZip(source, support_suffix_range=False)


def copy_member(z, info, dest):
    with z.open(info) as src, open(dest, "wb") as out:
        shutil.copyfileobj(src, out, length=16 * 1024 * 1024)
    print(f"  {info.file_size / 1e6:8.1f} MB  {dest}")


# Reuse the step-1 copies of the two OneK1K files if they exist.
for name in ("embedding_no_val_genes.pq", "celltype_fixed.pq"):
    found = list(Path(args.step1_dir).rglob(f"onek1k_training_data/{name}"))
    if found and not (DATA / name).exists():
        shutil.copy(found[0], DATA / name)
        print(f"  reused step-1 file -> {DATA / name}")

for key, source in (("training", args.training_zip), ("repro", args.repro_zip)):
    have = lambda local: (DATA / local).exists() or (local.endswith(".fa.gz") and (DATA / local[:-3]).exists())
    todo = [(end, local) for k, end, local in WANTED if k == key and not have(local)]
    want_outs = key == "repro" and not (DATA / "authors_outs").exists()
    if not todo and not want_outs:
        continue
    print(f"[{key}] reading table of contents ...")
    with open_zip(source) as z:
        infos = [i for i in z.infolist() if not i.filename.endswith("/")]
        for end, local in todo:
            hits = sorted((i for i in infos if i.filename.endswith(end)), key=lambda i: len(i.filename))
            if not hits:
                raise SystemExit(f"not found in {key} zip: {end}")
            copy_member(z, hits[0], DATA / local)
        if want_outs:
            outs = [i for i in infos if AUTHORS_OUTS in i.filename and i.filename.endswith(".h5ad")]
            (DATA / "authors_outs").mkdir()
            for i in outs:
                copy_member(z, i, DATA / "authors_outs" / Path(i.filename).name)
            print(f"  {len(outs)} prediction chunks (expected 11)")

# The script reads the genome with pyfaidx, which needs the plain FASTA (the shipped .fai is not reused;
# pyfaidx builds its own index on first use).
fa = DATA / "genome_human.fa"
if not fa.exists():
    print("unzipping genome (~3 GB) ...")
    with gzip.open(DATA / "genome_human.fa.gz", "rb") as src, open(fa, "wb") as out:
        shutil.copyfileobj(src, out, length=16 * 1024 * 1024)
    (DATA / "genome_human.fa.gz").unlink()

print("\nDone. Contents of ./data:")
for f in sorted(DATA.rglob("*")):
    if f.is_file():
        print(f"  {f.stat().st_size / 1e6:9.1f} MB  {f}")
print("Next: python make_chunk_input.py --n 200")
