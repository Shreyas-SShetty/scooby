"""Run a paired CIMA lead cis-eQTL pilot using the archived Scooby scores."""

import argparse
from importlib.metadata import version
import os
from pathlib import Path

import numpy as np
import pandas as pd

from evaluate_alphagenome2 import (
    KEYS, MODEL, ROOT, SCORER, add_windows, metrics, paired_intervals, score_rows, select_score,
)

PILOT_CELLS = [
    "CD4_Tcm_CXCR5", "CD8_Tem_CCR7neg", "NK_bright_XCL1", "MAIT_SLC4A10",
    "CD4_Tn_CCR7", "DC2_CD1C", "cMono_CD14", "Bn_TCL1A",
]
PILOT_PER_CELL = 12
TRACK_SUFFIX = "_scRNA-seq 10X 5-cellranger-cell"


def read_inputs():
    measured = pd.read_csv(ROOT / "data/cima_eqtl.tsv", sep="\t")
    scooby = pd.read_csv(ROOT / "eval_out_cima/merged.csv.gz").rename(
        columns={"pred_effect": "scooby_effect"},
    )
    if not measured.study_label.eq("CIMA").all() or not scooby.study_label.eq("CIMA").all():
        raise ValueError("Unexpected study in the CIMA sources")
    if measured.duplicated(KEYS).any() or scooby.duplicated(KEYS).any():
        raise ValueError("Expected exactly one row per variant, gene and CIMA cell type")
    joined = scooby.merge(measured[KEYS + ["beta"]], on=KEYS, how="left",
                          validate="one_to_one", suffixes=("", "_source"))
    if not np.isfinite(joined[["scooby_effect", "beta", "beta_source"]].to_numpy()).all():
        raise ValueError("Missing or non-finite values in the archived Scooby evaluation")
    if not np.allclose(joined.beta, joined.beta_source, rtol=1e-12, atol=1e-12):
        raise ValueError("Archived Scooby rows disagree with the measured CIMA effects")
    scored = joined.drop(columns="beta_source")
    scored["position_1based"] = scored.Pos + 1
    identifiers = (scored.Chromosome + "_" + scored.position_1based.astype(str)
                   + "_" + scored.Ref + "_" + scored.Alt)
    if not scored.variant.eq(identifiers).all():
        raise ValueError("Variant IDs disagree with the alleles or one-based positions")
    if not scored.Strand.isin(["+", "-"]).all():
        raise ValueError("Missing or invalid gene strand in the archived Scooby evaluation")
    # GENCODE starts/ends in the archived evaluation are one-based and inclusive.
    scored["tss_center"] = scored.gene_start.where(scored.Strand.eq("+"), scored.gene_end) - 1
    return measured, scored


def resolve_mapping(measured, scored, catalog):
    scooby_map = pd.read_csv(ROOT / "cima_mapping.csv").rename(
        columns={"group": "scooby_match_group", "note": "scooby_match_note"},
    )
    ag2_map = pd.read_csv(ROOT / "cima_ag2_reference_map.csv")
    mapped_labels = set(scooby_map.tissue_label)
    if (mapped_labels != set(ag2_map.tissue_label)
            or mapped_labels != set(scored.tissue_label)):
        raise ValueError("CIMA mappings must cover the cell types with archived Scooby scores")
    coverage = pd.concat([
        measured.groupby("tissue_label").size().rename("n_measured"),
        scored.groupby("tissue_label").size().rename("n_scooby"),
    ], axis=1).fillna(0).astype(int).reset_index()
    mapping = coverage.merge(scooby_map, on="tissue_label", how="left", validate="one_to_one")
    mapping = mapping.merge(ag2_map, on="tissue_label", how="left", validate="one_to_one")
    mapping = mapping.merge(
        catalog[["condition", "biosample_name", "sc_gene_expression"]].rename(
            columns={"condition": "ag2_condition"},
        ),
        on="ag2_condition", how="left", validate="many_to_one",
    )
    modes = ["observed", "imputed_supported"]
    if not mapping.loc[mapping.ag2_condition.notna(), "sc_gene_expression"].isin(modes).all():
        raise ValueError("A configured CIMA reference condition is unavailable for expression")
    mapping["in_pilot"] = mapping.tissue_label.isin(PILOT_CELLS)
    mapping.loc[mapping.n_scooby.eq(0), "ag2_match_quality"] = "no_archived_scooby_score"
    return mapping


def pick_pilot(scored, mapping):
    chosen = mapping.loc[mapping.in_pilot]
    if set(chosen.tissue_label) != set(PILOT_CELLS) or not chosen.sc_gene_expression.eq("observed").all():
        raise ValueError("Every pilot cell type requires an observed single-cell expression condition")
    rows = scored.merge(chosen[[
        "tissue_label", "pred_columns", "scooby_match_group", "ag2_condition",
        "ag2_match_quality", "ag2_match_note",
    ]], on="tissue_label", validate="many_to_one")
    # Stable ordering and seed fix the sample before scoring. Include small Scooby effects.
    rows = rows.sort_values(KEYS)
    pilot = rows.groupby("tissue_label", sort=True).sample(n=PILOT_PER_CELL, random_state=42)
    return add_windows(pilot.sort_values(KEYS).reset_index(drop=True))


def observed_tracks(model, pilot, organism):
    expected = set(pilot.ag2_condition + TRACK_SUFFIX)
    tracks = model.observed_metadata(organism).sc_gene_expression
    tracks = tracks.loc[tracks.name.isin(expected) & tracks.imputation_mode.eq("observed")
                        & tracks.suspension_type.eq("cell")]
    if set(tracks.name) != expected:
        raise ValueError("The common observed 10x 5′ expression assay is missing for a pilot condition")
    return tracks


def select_cima_score(tracks, gene):
    matching = tracks.loc[tracks.track_name.str.endswith(TRACK_SUFFIX, na=False)
                          & tracks.imputation_mode.eq("observed")
                          & tracks.suspension_type.eq("cell")]
    return select_score(matching, gene)


def summarize(scores, out):
    pooled = {"group": "all (pooled)", **metrics(scores), **paired_intervals(scores)}
    rows = [pooled]
    rows.extend({"group": cell, **metrics(group)}
                for cell, group in scores.groupby("tissue_label"))
    summary = pd.DataFrame(rows)
    summary.to_csv(out / "pilot_comparison.csv", index=False)
    print(summary[["group", "n_common", "scooby_spearman", "ag2_spearman",
                   "scooby_sign_accuracy", "ag2_sign_accuracy"]].round(3).to_string(index=False),
          flush=True)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "ag2_out/cima")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--prepare-only", action="store_true",
                        help="Build the mapping and pilot manifest using the saved condition catalog")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    args.out.mkdir(parents=True, exist_ok=True)
    catalog_path = args.out / "ag2_human_conditions.csv"

    if args.prepare_only:
        catalog = pd.read_csv(catalog_path)
    else:
        from alphagenome.models.v2 import dna_client

        key = os.environ.get("ALPHAGENOME_API_KEY")
        if not key:
            parser.error("set ALPHAGENOME_API_KEY in the current environment")
        model = dna_client.create(key, model_version=dna_client.ModelVersion[MODEL], timeout=20)
        organism = dna_client.Organism.HOMO_SAPIENS
        catalog = model.supported_conditions(organism)
        catalog.to_csv(catalog_path, index=False)

    measured, scooby = read_inputs()
    mapping = resolve_mapping(measured, scooby, catalog)
    mapping.to_csv(args.out / "condition_mapping.csv", index=False)
    pilot = pick_pilot(scooby, mapping)
    pilot.to_csv(args.out / "pilot_manifest.csv", index=False)
    print(f"CIMA has {len(measured):,} measured rows and {len(scooby):,} archived Scooby scores. "
          f"Pilot: {len(pilot)} rows across {pilot.tissue_label.nunique()} cell types.", flush=True)
    if args.prepare_only:
        return

    tracks = observed_tracks(model, pilot, organism)
    tracks.assign(model_version=MODEL, client_version=version("alphagenome"), scorer=SCORER).to_csv(
        args.out / "pilot_tracks.csv", index=False,
    )
    scores = score_rows(model, pilot, args.out / "cache", args.workers,
                        select_effect=select_cima_score, retry_timeouts=True)
    scores.to_csv(args.out / "pilot_scores.csv", index=False)
    summary = summarize(scores, args.out)
    if (summary.iloc[0].n_common < 0.8 * len(pilot)
            or summary.iloc[1:].n_common.lt(8).any()
            or summary.iloc[1:].ag2_spearman.isna().any()):
        raise SystemExit("CIMA pilot failed the coverage or correlation check")
    print("CIMA pilot passed the coverage and correlation checks.", flush=True)


if __name__ == "__main__":
    main()
