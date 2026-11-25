"""Evaluate AlphaGenome 2 on the scooby Figure 4 GTEx whole-blood pairs."""

import argparse
from importlib.metadata import version
import os
from pathlib import Path

import numpy as np
import pandas as pd

from evaluate_alphagenome2 import (
    MODEL, ROOT, add_windows, metrics, paired_intervals, read_scooby, score_rows,
)

CONDITION = "UBERON:0013756"
TRACK = f"{CONDITION}_polyA plus RNA-seq-gtex-paired-ended"
SCORER = "RNA_SEQ"
CLEAR_CUTOFF = 0.05


def prepare_manifest(out):
    saved = pd.read_csv(ROOT / "eval_onek1k/eval_out_authors_eqtl/merged.csv.gz")
    rows = saved.loc[saved.study.eq("GTEx") & saved.tissue.eq("blood")].rename(
        columns={"pred_effect": "scooby_effect"},
    ).sort_values(["variant", "gene_id"]).reset_index(drop=True)
    if rows.empty or rows.duplicated(["variant", "gene_id"]).any():
        raise ValueError("Expected unique GTEx whole-blood variant–gene pairs")
    rows["position_1based"] = rows.Pos + 1
    identifiers = (rows.Chromosome + "_" + rows.position_1based.astype(str)
                   + "_" + rows.Ref + "_" + rows.Alt)
    if not rows.variant.eq(identifiers).all():
        raise ValueError("Variant IDs disagree with the REF, ALT, or one-based positions")
    if not np.isfinite(rows[["scooby_effect", "beta"]].to_numpy()).all():
        raise ValueError("GTEx benchmark includes non-finite effects")

    scooby = read_scooby(sorted((ROOT / "data/authors_outs").glob("outs_*.h5ad")))
    scooby["tss_center"] = (scooby.Start + scooby.End) // 2
    centers = scooby[["variant", "gene_id", "tss_center"]].drop_duplicates()
    rows = rows.merge(centers, on=["variant", "gene_id"], how="left", validate="one_to_one")
    if rows.tss_center.isna().any():
        raise ValueError("GTEx overlap is missing a scooby gene window")
    rows = add_windows(rows)
    rows["ag2_condition"] = CONDITION
    rows["scooby_clear"] = rows.scooby_effect.abs().gt(CLEAR_CUTOFF)
    rows.to_csv(out / "manifest.csv", index=False)
    print(f"{len(rows)} GTEx whole-blood associations, {rows.scooby_clear.sum()} scooby-clear",
          flush=True)
    return rows


def select_gtex_score(tracks, gene):
    target = tracks.loc[tracks.gene_id.eq(gene) & tracks.track_name.eq(TRACK)
                        & tracks.imputation_mode.eq("observed")]
    if len(target) > 1:
        raise ValueError(f"Multiple scores for {gene} on the fixed GTEx track")
    if target.empty or not np.isfinite(target.raw_score.iloc[0]):
        return np.nan, "", "", "target_or_track_missing"
    return target.raw_score.iloc[0], "observed", TRACK, "scored"


def write_comparison(rows, out, phase):
    summaries = []
    for name, cohort in [
        ("all matchable pairs", rows),
        ("scooby clear pairs", rows.loc[rows.scooby_clear]),
    ]:
        common = cohort.dropna(subset=["scooby_effect", "ag2_effect"])
        summary = {"subset": name, **metrics(cohort), **paired_intervals(cohort)}
        if not common.empty:
            clear_beta = common.beta.abs().gt(CLEAR_CUTOFF)
            for model in ("scooby", "ag2"):
                hits = clear_beta & np.sign(common[f"{model}_effect"]).eq(np.sign(common.beta))
                summary[f"{model}_correct_with_beta_cutoff"] = int(hits.sum())
                summary[f"{model}_direction_with_beta_cutoff"] = hits.mean()
        summaries.append(summary)
    result = pd.DataFrame(summaries)
    result.to_csv(out / f"{phase}_comparison.csv", index=False)
    print(result[["subset", "n_common", "scooby_spearman", "ag2_spearman",
                  "scooby_sign_accuracy", "ag2_sign_accuracy"]].round(3).to_string(index=False),
          flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "ag2_out/gtex")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    args.out.mkdir(parents=True, exist_ok=True)
    rows = prepare_manifest(args.out)
    if args.prepare_only:
        return

    from alphagenome.models.v2 import dna_client

    key = os.environ.get("ALPHAGENOME_API_KEY")
    if not key:
        parser.error("set ALPHAGENOME_API_KEY in the current environment")
    model = dna_client.create(key, model_version=dna_client.ModelVersion[MODEL], timeout=20)
    rna = model.observed_metadata(dna_client.Organism.HOMO_SAPIENS).rna_seq
    track = rna.loc[rna.name.eq(TRACK) & rna.biosample_curie.eq(CONDITION)
                    & rna.strand.eq(".") & rna.data_source.eq("gtex")
                    & rna.imputation_mode.eq("observed")]
    if len(track) != 1:
        raise ValueError("The exact observed AlphaGenome 2 GTEx whole-blood track is unavailable")
    track.assign(model_version=MODEL, client_version=version("alphagenome"), scorer=SCORER).to_csv(
        args.out / "track.csv", index=False,
    )
    pilot = rows.sample(n=min(32, len(rows)), random_state=42).sort_index()
    for phase, selected in [("pilot", pilot), ("full", rows)]:
        print(f"Scoring {phase}: {len(selected)} associations with {TRACK}", flush=True)
        scores = score_rows(model, selected, args.out / "cache", args.workers,
                            scorer=SCORER, select_effect=select_gtex_score, retry_timeouts=True)
        scores.to_csv(args.out / f"{phase}_scores.csv", index=False)
        result = write_comparison(scores, args.out, phase)
        if phase == "pilot":
            if result.iloc[0].n_common < len(pilot) * 0.8 or pd.isna(result.iloc[0].ag2_spearman):
                raise SystemExit("Pilot failed the coverage or correlation check")
            print("Pilot passed. Starting the full panel.", flush=True)


if __name__ == "__main__":
    main()
