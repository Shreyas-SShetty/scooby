"""Score all archived CIMA Scooby associations with AlphaGenome 2."""

import argparse
import fcntl
from functools import partial
from hashlib import sha256
from importlib.metadata import version
import json
import os
from pathlib import Path
import tempfile
from datetime import datetime, timezone

import pandas as pd

from evaluate_alphagenome2 import (
    KEYS, MODEL, REQUEST, ROOT, SCORER, add_windows, metrics, paired_intervals,
    score_rows, select_score, write_cache,
)
from evaluate_alphagenome2_cima import TRACK_SUFFIX, read_inputs, resolve_mapping

BATCH_REQUESTS = 200


def atomic_text(path, contents):
    with tempfile.NamedTemporaryFile(
        mode="w", dir=path.parent, delete=False, encoding="utf-8"
    ) as handle:
        tmp = Path(handle.name)
    try:
        tmp.write_text(contents, encoding="utf-8")
        tmp.replace(path)
    finally:
        tmp.unlink(missing_ok=True)


def prepare(out):
    measured, scooby = read_inputs()
    mapping = resolve_mapping(measured, scooby, pd.read_csv(out / "ag2_human_conditions.csv"))
    mapped = mapping.loc[mapping.n_scooby.gt(0)]
    if mapped.ag2_condition.isna().any():
        raise ValueError("Every CIMA cell type with Scooby scores needs a reference in the full run")
    observed = pd.read_csv(out / "full_observed_tracks.csv")
    chosen = {}
    for condition in mapped.loc[mapped.sc_gene_expression.eq("observed"), "ag2_condition"].unique():
        names = observed.loc[observed.biosample_curie.eq(condition)
                             & observed.imputation_mode.eq("observed")]
        if names.empty:
            raise ValueError(f"No observed single-cell expression tracks for {condition}")
        preferred = condition + TRACK_SUFFIX
        whole_cell = names.loc[names.suspension_type.eq("cell")]
        if names.name.eq(preferred).any():
            chosen[condition] = preferred
        elif not whole_cell.empty:
            chosen[condition] = ";".join(sorted(set(whole_cell.name)))
        else:
            chosen[condition] = ";".join(sorted(set(names.name)))
    mapping["ag2_expected_tracks"] = mapping.ag2_condition.map(chosen)
    mapping["ag2_scope"] = mapping.sc_gene_expression.map({
        "observed": "observed_reference", "imputed_supported": "imputed_reference",
    })
    mapping.loc[mapping.ag2_match_quality.eq("shared_scooby_proxy"), "ag2_scope"] = "shared_proxy"
    mapping.to_csv(out / "full_condition_mapping.csv", index=False)

    manifest = scooby.merge(mapping[[
        "tissue_label", "pred_columns", "scooby_match_group", "ag2_condition",
        "ag2_match_quality", "ag2_match_note", "sc_gene_expression", "ag2_scope",
        "ag2_expected_tracks",
    ]], on="tissue_label", validate="many_to_one")
    manifest = add_windows(manifest).sort_values(REQUEST + KEYS).reset_index(drop=True)
    manifest_path = out / "full_manifest.csv"
    atomic_text(manifest_path, manifest.to_csv(index=False))
    fingerprint = sha256()
    fingerprint.update(f"cima-full-v1:{version('alphagenome')}:{MODEL}:{SCORER}:".encode())
    fingerprint.update(manifest_path.read_bytes())
    return manifest, fingerprint.hexdigest()


def select_full_score(tracks, gene, *, mode, expected):
    matching = tracks.loc[tracks.imputation_mode.eq(mode)]
    if expected:
        matching = matching.loc[matching.track_name.isin(expected)]
    else:
        preferred = matching.track_name.str.endswith(TRACK_SUFFIX, na=False)
        if preferred.any():
            matching = matching.loc[preferred]
    return select_score(matching, gene)


def batches(manifest):
    for _, condition_rows in manifest.groupby("ag2_condition", sort=True):
        groups = list(condition_rows.groupby(REQUEST, sort=False).indices.values())
        for start in range(0, len(groups), BATCH_REQUESTS):
            block = groups[start:start + BATCH_REQUESTS]
            indices = sorted(i for group in block for i in group)
            yield condition_rows.iloc[indices], len(block)


def select_for_rows(rows):
    first = rows.iloc[0]
    expected = first.ag2_expected_tracks
    return partial(select_full_score, mode=first.sc_gene_expression,
                   expected=expected.split(";") if isinstance(expected, str) else [])


def write_results(scores, out):
    pooled = [
        ("all 48 CIMA labels (pooled)", scores),
        ("without shared proxies (pooled)", scores.loc[scores.ag2_scope.ne("shared_proxy")]),
    ]
    rows = [{"group": label, **metrics(subset), **paired_intervals(subset)}
            for label, subset in pooled]
    for column in ("ag2_scope", "scooby_match_group"):
        rows.extend({"group": label, **metrics(group)}
                    for label, group in scores.groupby(column))
    rows.extend({"group": label, **metrics(group)}
                for label, group in scores.groupby("tissue_label"))
    summary = pd.DataFrame(rows)
    atomic_text(out / "full_scores.csv", scores.to_csv(index=False))
    atomic_text(out / "full_comparison.csv", summary.to_csv(index=False))
    print(summary.iloc[:len(pooled)][[
        "group", "n_common", "scooby_spearman", "ag2_spearman",
        "scooby_sign_accuracy", "ag2_sign_accuracy",
    ]].round(3).to_string(index=False), flush=True)


def run_full(model, manifest, out, workers, run_id):
    progress_file = out / "full_progress.json"
    plan = list(batches(manifest))
    expected_requests = sum(n_requests for _, n_requests in plan)
    if expected_requests != len(manifest[REQUEST].drop_duplicates()):
        raise ValueError("Run plan must process each variant window and condition exactly once")
    checkpoint_dir = out / "full_shards" / run_id
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    parts = []
    progress = {"run_id": run_id, "state": "running", "n_rows": len(manifest),
                "n_cell_types": manifest.tissue_label.nunique(), "n_requests": expected_requests,
                "requests_processed": 0, "rows_processed": 0, "rows_with_ag2_score": 0,
                "completed_batches": 0, "n_batches": len(plan)}

    def update_progress():
        progress["updated_at_utc"] = datetime.now(timezone.utc).isoformat()
        atomic_text(progress_file, json.dumps(progress, indent=2) + "\n")

    update_progress()
    try:
        for index, (rows, n_requests) in enumerate(plan):
            path = checkpoint_dir / f"{index:04d}.csv.gz"
            if path.exists():
                batch = pd.read_csv(path)
                if not batch[KEYS].equals(rows[KEYS].reset_index(drop=True)):
                    raise ValueError(f"Checkpoint does not match the prepared rows: {path.name}")
            else:
                batch = score_rows(model, rows, out / "cache", workers,
                                   select_effect=select_for_rows(rows), retry_timeouts=True)
                write_cache(batch, path)
            parts.append(batch)
            progress["completed_batches"] += 1
            progress["requests_processed"] += n_requests
            progress["rows_processed"] += len(rows)
            progress["rows_with_ag2_score"] += int(batch.ag2_effect.notna().sum())
            update_progress()
            print(f"Full CIMA progress: {progress['requests_processed']}/{expected_requests} requests, "
                  f"{progress['rows_processed']}/{len(manifest)} rows", flush=True)

        scores = pd.concat(parts, ignore_index=True)
        if len(scores) != len(manifest) or scores.duplicated(KEYS).any():
            raise ValueError("Checkpoints do not cover the full benchmark")
        write_results(scores.sort_values(KEYS), out)
    except Exception as error:
        progress["state"] = "interrupted"
        progress["error_type"] = type(error).__name__
        if callable(getattr(error, "code", None)):
            progress["api_status"] = error.code().name
        update_progress()
        raise

    progress["state"] = "complete"
    progress["score_status"] = scores.ag2_status.value_counts().to_dict()
    update_progress()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "ag2_out/cima")
    parser.add_argument("--workers", type=int, default=5)
    parser.add_argument("--prepare-only", action="store_true",
                        help="Build the full plan from saved AlphaGenome catalog metadata")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / "full.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            parser.error(f"A full CIMA run already owns {args.out / 'full.lock'}")
        catalog_path = args.out / "ag2_human_conditions.csv"
        track_path = args.out / "full_observed_tracks.csv"
        have_metadata = catalog_path.is_file() and track_path.is_file()
        if args.prepare_only and not have_metadata:
            parser.error("the saved catalog and observed tracks are required for offline preparation")
        if not args.prepare_only:
            from alphagenome.models.v2 import dna_client

            key = os.environ.get("ALPHAGENOME_API_KEY")
            if not key:
                parser.error("set ALPHAGENOME_API_KEY in the current environment")
            model = dna_client.create(key, model_version=dna_client.ModelVersion[MODEL], timeout=20)
            if not have_metadata:
                organism = dna_client.Organism.HOMO_SAPIENS
                catalog = model.supported_conditions(organism)
                tracks = model.observed_metadata(organism).sc_gene_expression
                atomic_text(catalog_path, catalog.to_csv(index=False))
                atomic_text(track_path, tracks.to_csv(index=False))
        manifest, run_id = prepare(args.out)
        print(f"Full CIMA plan: {len(manifest):,} rows across {manifest.tissue_label.nunique()} labels, "
              f"{len(manifest[REQUEST].drop_duplicates()):,} unique requests.", flush=True)
        if args.prepare_only:
            return

        run_full(model, manifest, args.out, args.workers, run_id)


if __name__ == "__main__":
    main()
