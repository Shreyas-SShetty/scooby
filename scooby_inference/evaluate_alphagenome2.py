"""Run a matched OneK1K eQTL comparison using the released scooby predictions."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
from importlib.metadata import version
import json
import os
from pathlib import Path
import tempfile
import time

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent
KEYS = ["variant", "gene_id", "tissue_label"]
REQUEST = ["Chromosome", "position_1based", "Ref", "Alt", "interval_start", "ag2_condition"]
PILOT_CELLS = ["B cell", "CD4+ TCM cell", "NK cell", "monocyte"]
SCORER = "SC_GENE_EXPRESSION"
MODEL = "ALL_FOLDS"
WIDTH = 1_048_576
BLOOD = "UBERON:0000178"
MARROW = "UBERON:0002371"
# OneK1K labels are from the eQTL Catalogue. CD4+ T cell and CD8+ T cell
# are its naive populations. The other proxies are recorded in every output.
CELL_MAP = {
    "B cell": (BLOOD, "0000236", "direct"),
    "CD16+ monocyte": (BLOOD, "0002396", "direct"),
    "CD4+ CTL cell": (BLOOD, "0000934", "direct"),
    "CD4+ TCM cell": (BLOOD, "0000904", "direct"),
    "CD4+ TEM cell": (BLOOD, "0000905", "direct"),
    "CD4+ T cell": (BLOOD, "0000895", "direct"),
    "CD56+ NK cell": (BLOOD, "0000938", "direct"),
    "CD8+ TCM cell": (BLOOD, "0000907", "direct"),
    "CD8+ TEM cell": (BLOOD, "0000913", "direct"),
    "CD8+ T cell": (BLOOD, "0000900", "direct"),
    "MAIT cell": (BLOOD, "0000940", "direct"),
    "NK cell": (BLOOD, "0000623", "direct"),
    "Treg memory": (BLOOD, "0000792", "broader regulatory T cell"),
    "dendritic cell": (BLOOD, "0000990", "conventional dendritic cell proxy"),
    "gdT cell": (BLOOD, "0000798", "direct"),
    "hematopoietic precursor cell": (MARROW, "0000037", "bone marrow stem cell proxy"),
    "memory B cell": (BLOOD, "0000787", "direct"),
    "monocyte": (BLOOD, "0001054", "direct"),
    "plasmablast": (BLOOD, "0000980", "direct"),
    "plasmacytoid dendritic cell": (BLOOD, "0000784", "direct"),
    "platelet": (BLOOD, "0000233", "direct"),
}
MODES = {"observed": 0, "imputed_supported": 1, "imputed_exploratory": 2}
TRACK_COLS = ["gene_id", "raw_score", "track_name", "imputation_mode", "suspension_type"]


def read_scooby(paths):
    chunks = []
    for path in paths:
        pred = ad.read_h5ad(path)
        meta = pred.obs[["variant", "gene_id", "Chromosome", "Start", "End"]].copy()
        meta = meta.rename(columns={"Chromosome": "scooby_chromosome"})
        effects = pd.DataFrame(pred.X, index=meta.index, columns=pred.var_names)
        chunks.append(pd.concat([meta, effects], axis=1))
    if not chunks:
        raise ValueError("No scooby prediction chunks found")
    return pd.concat(chunks, ignore_index=True)


def add_windows(rows):
    rows = rows.copy()
    rows["tss_distance"] = (rows.Pos - rows.tss_center).abs().astype("Int64")
    rows["window_anchor"] = np.where(rows.tss_center.notna(), "scooby_gene_window", "variant")
    lengths = pd.read_csv(ROOT / "data/genome_human.fa.fai", sep="\t", header=None)
    chr_length = rows.Chromosome.map(lengths.set_index(0)[1])
    if chr_length.isna().any() or (chr_length < WIDTH).any():
        raise ValueError("Missing or too-short benchmark chromosome in the reference index")
    anchor = rows.tss_center.fillna(rows.Pos).astype(int)
    distal = (anchor - rows.Pos).abs().ge(WIDTH // 2)
    anchor = anchor.where(~distal, (anchor + rows.Pos) // 2)
    rows.loc[distal, "window_anchor"] = "gene_variant_midpoint"
    center_start = (anchor - WIDTH // 2).clip(lower=0, upper=chr_length - WIDTH)
    if not (rows.Pos.ge(center_start) & rows.Pos.lt(center_start + WIDTH)).all():
        raise ValueError("A variant falls outside its AlphaGenome interval")
    rows["interval_start"] = center_start
    rows["interval_end"] = center_start + WIDTH
    return rows


def prepare_manifest(out):
    eqtl = pd.read_csv(ROOT / "data/all_highpip_susie.tsv", sep="\t")
    eqtl = eqtl.loc[eqtl.study_label.eq("OneK1K")].sort_values(
        "sample_size", ascending=False, kind="stable"
    ).drop_duplicates(KEYS).sort_values(KEYS).reset_index(drop=True)
    eqtl["position_1based"] = eqtl.Pos + 1
    identifiers = (eqtl.Chromosome + "_" + eqtl.position_1based.astype(str)
                   + "_" + eqtl.Ref + "_" + eqtl.Alt)
    if not eqtl.variant.eq(identifiers).all():
        raise ValueError("Variant IDs disagree with the REF, ALT, or one-based positions")
    if not np.isfinite(eqtl.beta).all():
        raise ValueError("OneK1K contains non-finite measured effects")
    if set(eqtl.tissue_label) != set(CELL_MAP):
        raise ValueError("Cell mapping does not cover exactly the OneK1K labels")

    scooby = read_scooby(sorted((ROOT / "data/authors_outs").glob("outs_*.h5ad")))
    scooby["tss_center"] = (scooby.Start + scooby.End) // 2
    gene_keys = ["gene_id", "scooby_chromosome"]
    centers = scooby[gene_keys + ["tss_center"]].drop_duplicates()
    if centers.duplicated(gene_keys).any():
        raise ValueError("Conflicting scooby gene window centers")
    eqtl = eqtl.merge(
        centers.rename(columns={"scooby_chromosome": "Chromosome"}),
        on=["gene_id", "Chromosome"], how="left", validate="many_to_one",
    )
    score_groups = []
    for tissue in CELL_MAP:
        column = tissue.replace(" ", "_")
        label = eqtl.loc[eqtl.tissue_label.eq(tissue)]
        scores = scooby[["variant", "gene_id", column]].drop_duplicates().rename(
            columns={column: "scooby_effect"}
        )
        score_groups.append(label.merge(scores, on=KEYS[:2], how="left", validate="one_to_one"))
    eqtl = pd.concat(score_groups, ignore_index=True).sort_values(KEYS).reset_index(drop=True)
    eqtl = add_windows(eqtl)
    eqtl["ag2_condition"] = eqtl.tissue_label.map({
        cell: f"{tissue}_CL:{term}" for cell, (tissue, term, _) in CELL_MAP.items()
    })
    eqtl["cell_match"] = eqtl.tissue_label.map({
        cell: quality for cell, (_, _, quality) in CELL_MAP.items()
    })
    eqtl.to_csv(out / "manifest.csv", index=False)
    print(f"{len(eqtl)} OneK1K rows, {eqtl.scooby_effect.notna().sum()} with scooby scores",
          flush=True)
    return eqtl


def add_condition_metadata(rows, catalog, out):
    mapping = rows[["tissue_label", "ag2_condition", "cell_match"]].drop_duplicates()
    mapping = mapping.merge(
        catalog[["condition", "biosample_name", "sc_gene_expression"]],
        left_on="ag2_condition", right_on="condition", how="left", validate="many_to_one",
    ).drop(columns="condition")
    mapping.to_csv(out / "conditions.csv", index=False)
    missing = mapping.loc[~mapping.sc_gene_expression.isin(MODES), "tissue_label"].tolist()
    if missing:
        raise ValueError(f"AlphaGenome catalog lacks expression conditions for: {missing}")
    return rows.merge(mapping, on=["tissue_label", "ag2_condition", "cell_match"],
                      validate="many_to_one")


def pilot_rows(rows):
    samples = []
    for cell in PILOT_CELLS:
        candidates = rows.loc[rows.tissue_label.eq(cell) & rows.scooby_effect.notna()]
        samples.append(candidates.sample(n=min(16, len(candidates)), random_state=42))
    return pd.concat(samples).sort_values(KEYS).reset_index(drop=True)


def write_cache(table, path):
    with tempfile.NamedTemporaryFile(dir=path.parent, suffix=".csv.gz", delete=False) as handle:
        tmp = Path(handle.name)
    try:
        table.to_csv(tmp, index=False)
        tmp.replace(path)
    finally:
        tmp.unlink(missing_ok=True)


def get_tracks(model, key, cache, scorer=SCORER):
    from alphagenome.data import genome
    from alphagenome.models.v2 import variant_scorers

    request_id = sha256(json.dumps([version("alphagenome"), MODEL, scorer, WIDTH, *key])
                        .encode()).hexdigest()
    path = cache / f"{request_id}.csv.gz"
    if not path.exists():
        chrom, position, ref, alt, start, condition = key
        variant = genome.Variant(chrom, position, ref, alt)
        interval = genome.Interval(chrom, start, start + WIDTH)
        scores = model.score_variant(
            interval=interval, variant=variant, condition=condition,
            variant_scorers=[variant_scorers.RECOMMENDED_VARIANT_SCORERS[scorer]],
        )
        table = variant_scorers.tidy_scores(
            scores, match_gene_strand=True, filter_protein_coding=False,
        )
        tracks = (table.reindex(columns=TRACK_COLS) if table is not None and not table.empty
                  else pd.DataFrame(columns=TRACK_COLS))
        write_cache(tracks, path)
    return pd.read_csv(path)


def select_score(tracks, gene):
    target = tracks.loc[tracks.gene_id.eq(gene) & np.isfinite(tracks.raw_score.astype(float))]
    if target.empty:
        return np.nan, "", "", "target_missing"
    target = target.assign(mode_priority=target.imputation_mode.map(MODES))
    if target.mode_priority.isna().any():
        raise ValueError(f"Unrecognized AlphaGenome track metadata for gene {gene}")
    target = target.loc[target.mode_priority.eq(target.mode_priority.min())]
    whole_cell = target.suspension_type.eq("cell")
    if whole_cell.any():
        target = target.loc[whole_cell]
    return (target.raw_score.mean(), target.imputation_mode.iloc[0],
            ";".join(sorted(set(target.track_name))), "scored")


def score_rows(model, rows, cache, workers, *, scorer=SCORER, select_effect=select_score,
               retry_timeouts=False):
    import grpc

    cache.mkdir(exist_ok=True)
    groups = list(rows.groupby(REQUEST, sort=False))

    def score_group(group):
        key, subset = group
        attempts = 2 if retry_timeouts else 1
        for attempt in range(attempts):
            try:
                tracks = get_tracks(model, key, cache, scorer)
            except grpc.RpcError as error:
                if error.code() == grpc.StatusCode.INVALID_ARGUMENT:
                    print(f"Rejected request {key}: {error.details()}", flush=True)
                    return [(index, np.nan, "", "", "invalid_request") for index in subset.index]
                if retry_timeouts and error.code() == grpc.StatusCode.DEADLINE_EXCEEDED:
                    if attempt + 1 < attempts:
                        print(f"Retrying timed-out request {key}", flush=True)
                        time.sleep(2)
                        continue
                    print(f"Request timed out again {key}", flush=True)
                    return [(index, np.nan, "", "", "request_timeout") for index in subset.index]
                raise
            return [(index, *select_effect(tracks, gene)) for index, gene in subset.gene_id.items()]

    results = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        # Submit one small batch at a time so a quota/auth failure stops new requests.
        batch_size = workers * 5
        for start in range(0, len(groups), batch_size):
            for entries in pool.map(score_group, groups[start:start + batch_size]):
                results.extend(entries)
            print(f"Completed {min(start + batch_size, len(groups))}/{len(groups)} requests",
                  flush=True)
    scores = pd.DataFrame(results, columns=[
        "index", "ag2_effect", "ag2_track_mode", "ag2_tracks", "ag2_status"
    ]).set_index("index")
    return rows.join(scores, validate="one_to_one")


def metrics(group):
    common = group.dropna(subset=["scooby_effect", "ag2_effect"])
    result = {
        "n_manifest": len(group),
        "n_scooby": group.scooby_effect.notna().sum(),
        "n_ag2": group.ag2_effect.notna().sum(),
        "n_common": len(common),
    }
    directional = common.loc[common.beta.ne(0)]
    result["n_directional"] = len(directional)
    for name, column in [("scooby", "scooby_effect"), ("ag2", "ag2_effect")]:
        effect = common[column]
        result[f"{name}_spearman"] = (
            spearmanr(common.beta, effect).statistic
            if len(common) >= 3 and effect.nunique() > 1 and common.beta.nunique() > 1
            else np.nan
        )
        result[f"{name}_sign_accuracy"] = (
            np.sign(directional[column]).eq(np.sign(directional.beta)).mean()
        )
        # Stable variant-key order breaks ties at the fraction cutoff.
        top = directional.loc[directional[column].abs().nlargest(
            (len(directional) + 3) // 4
        ).index]
        result[f"{name}_sign_accuracy_top_quarter"] = (
            np.sign(top[column]).eq(np.sign(top.beta)).mean()
        )
    result["delta_spearman"] = result["ag2_spearman"] - result["scooby_spearman"]
    result["delta_sign_accuracy"] = result["ag2_sign_accuracy"] - result["scooby_sign_accuracy"]
    return result


def paired_intervals(scores, n_resamples=1000):
    """Resample genes with replacement, retaining all their paired observations."""
    common = scores.dropna(subset=["scooby_effect", "ag2_effect"]).reset_index(drop=True)
    blocks = list(common.groupby("gene_id").indices.values())
    if len(blocks) < 2:
        return {}
    beta, scooby, ag2 = common[["beta", "scooby_effect", "ag2_effect"]].to_numpy().T
    sign_delta = ((np.sign(ag2) == np.sign(beta)).astype(float)
                  - (np.sign(scooby) == np.sign(beta)))
    rng = np.random.default_rng(42)
    rho_deltas, sign_deltas = [], []
    for _ in range(n_resamples):
        indices = np.concatenate([blocks[i] for i in rng.integers(len(blocks), size=len(blocks))])
        b, s, a = beta[indices], scooby[indices], ag2[indices]
        if len(indices) >= 3 and all(np.ptp(values) > 0 for values in (b, s, a)):
            rho_deltas.append(spearmanr(b, a).statistic - spearmanr(b, s).statistic)
        nonzero = indices[b != 0]
        if len(nonzero):
            sign_deltas.append(sign_delta[nonzero].mean())
    result = {}
    for metric, samples in [("delta_spearman", rho_deltas), ("delta_sign_accuracy", sign_deltas)]:
        if samples:
            result[f"{metric}_ci_low"], result[f"{metric}_ci_high"] = np.quantile(
                samples, [0.025, 0.975],
            )
    return result


def write_comparison(scores, out, phase):
    pooled = {"group": "all (pooled)", **metrics(scores), **paired_intervals(scores)}
    summaries = [pooled]
    summaries.extend({"group": cell, **metrics(group)}
                     for cell, group in scores.groupby("tissue_label"))
    direct = scores.loc[scores.cell_match.eq("direct") & scores.sc_gene_expression.eq("observed")]
    summaries.append({"group": "direct observed conditions (pooled)", **metrics(direct),
                      **paired_intervals(direct)})
    summary = pd.DataFrame(summaries)
    summary.to_csv(out / f"{phase}_comparison.csv", index=False)
    bins = pd.cut(scores.tss_distance, [-1, 1000, 10_000, 35_000, 100_000, np.inf])
    distance = [{"distance_bp": str(label), **metrics(group)}
                for label, group in scores.groupby(bins, observed=True)]
    pd.DataFrame(distance).to_csv(out / f"{phase}_distance.csv", index=False)
    print(summary[["group", "n_common", "scooby_spearman", "ag2_spearman",
                   "scooby_sign_accuracy", "ag2_sign_accuracy"]].round(3).to_string(index=False),
          flush=True)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=["prepare", "pilot", "full", "run", "summarize"],
                        help="'run' scores and checks the pilot before running the full set")
    parser.add_argument("--out", type=Path, default=ROOT / "ag2_out")
    parser.add_argument("--workers", type=int, default=3)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    args.out.mkdir(parents=True, exist_ok=True)
    if args.phase == "summarize":
        saved = [p for p in ("pilot", "full") if (args.out / f"{p}_scores.csv").is_file()]
        if not saved:
            parser.error(f"no completed score tables in {args.out}")
        for phase in saved:
            write_comparison(pd.read_csv(args.out / f"{phase}_scores.csv"), args.out, phase)
        return
    rows = prepare_manifest(args.out)
    if args.phase == "prepare":
        return

    from alphagenome.models.v2 import dna_client

    key = os.environ.get("ALPHAGENOME_API_KEY")
    if not key:
        parser.error("set ALPHAGENOME_API_KEY in the current environment")
    model = dna_client.create(key, model_version=dna_client.ModelVersion[MODEL], timeout=20)
    catalog = model.supported_conditions(dna_client.Organism.HOMO_SAPIENS)
    rows = add_condition_metadata(rows, catalog, args.out)

    phases = ["pilot", "full"] if args.phase == "run" else [args.phase]
    for phase in phases:
        selected = pilot_rows(rows) if phase == "pilot" else rows
        print(f"Scoring {phase}: {len(selected)} rows using {version('alphagenome')} / {MODEL}",
              flush=True)
        scores = score_rows(model, selected, args.out / "cache", args.workers)
        scores.to_csv(args.out / f"{phase}_scores.csv", index=False)
        summary = write_comparison(scores, args.out, phase)
        if phase == "pilot":
            min_scored = (len(selected) + 1) // 2
            per_cell = summary.loc[summary.group.isin(PILOT_CELLS)]
            if (summary.iloc[0].n_common < min_scored
                    or per_cell.ag2_spearman.notna().sum() < 2):
                raise SystemExit("Pilot did not produce enough usable matched predictions")
            print("Pilot coverage and metrics checks passed.", flush=True)


if __name__ == "__main__":
    main()
