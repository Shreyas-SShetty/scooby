# AlphaGenome 2 eQTL benchmarks

## AlphaGenome 2 versus Scooby on OneK1K

From the repository root, use a Python 3.12 environment and install `alphagenome==2.0.0`.
With your research API key set as `ALPHAGENOME_API_KEY` in the running environment:

```bash
python scooby_inference/evaluate_alphagenome2.py run
```

`run` scores 64 pilot rows across four cell types, checks that at least half have
paired scores and that at least two cell types yield usable correlations, then
scores the full set. `prepare` builds the manifest offline. `pilot` and `full`
run those phases individually. `summarize` rebuilds the comparison from saved
scores offline. Reuse the same `--out` directory to resume with cached scores.
The default is `scooby_inference/ag2_out`.

The benchmark uses the OneK1K rows in the high-PIP SuSiE table and the 11 released
scooby prediction chunks used by the existing evaluator. It keeps the largest
study sample for each variant, Ensembl gene and cell type. All three keys must
match for a comparison. The local 200-row GPU rerun only covers a subset.

AlphaGenome 2 uses the `ALL_FOLDS` model and the signed `SC_GENE_EXPRESSION` scorer.
Its 1 MiB input uses the released scooby gene-window center when available. It
uses the variant as the center for genes absent from the scooby annotation and
the midpoint for variants outside a gene-centered window. The manifest records
the window for each row. Coordinates are GRCh38. The manifest checks that REF,
ALT and the one-based position agree with the input variant identifier.

`conditions.csv` records the live cell condition catalog matches and the three
proxies. For each target gene, scores use strand-compatible tracks and prefer
observed, then supported imputed, then exploratory imputed tracks. Within that
priority, whole-cell tracks take precedence over other suspension types. The
benchmark averages remaining tracks and records their names per row.

The `pilot_` and `full_` outputs include row-level scores, per-cell and pooled
Spearman correlation and sign accuracy against measured beta, and metrics by
distance to the window center. Each comparison uses the same scored rows for
both models and reports each model's coverage. Top-quarter sign accuracy uses
an equal number of predictions selected by each model's absolute effect size.
The pooled comparisons give 95% intervals for the paired differences using
1,000 bootstrap resamples of genes, retaining every matched observation for a
selected gene. The pooled row for direct observed conditions checks sensitivity
to the mapping proxies.

## AlphaGenome 2 versus Scooby on GTEx whole blood

Use the same Python environment and `ALPHAGENOME_API_KEY` to run the 32-row
pilot and then the full GTEx whole-blood comparison:

```bash
python scooby_inference/evaluate_alphagenome2_gtex.py
```

The input is the 908 variant–gene associations from the saved scooby Figure 4
benchmark. AlphaGenome 2 scores the observed GTEx bulk RNA-seq track for venous
blood using the `RNA_SEQ` variant scorer and the same window-centering policy
as the OneK1K comparison. The Scooby baseline comes from the Figure 4 archive
`eval_onek1k/eval_out_authors_eqtl/merged.csv.gz`. The output in `ag2_out/gtex` includes the exact track,
row-level scores, and paired performance for all associations with both scores.
The comparison also reports the subset with `|scooby effect| > 0.05`, using the
same rows for both models. The stricter direction metric also requires `|beta| > 0.05` for a correct prediction.

## AlphaGenome 2 versus Scooby on CIMA

Use the same environment and `ALPHAGENOME_API_KEY` to run a paired pilot:

```bash
python scooby_inference/evaluate_alphagenome2_cima.py
```

The pilot draws 12 rows per cell type from eight cell types with archived Scooby
scores and an observed AlphaGenome 2 reference for the same 10x 5′ assay. The
sample uses seed 42 and includes small predicted effects. Both models are compared
with the same measured CIMA lead cis-eQTL slopes on the shared scored rows. The
AlphaGenome 2 input is centered on the GENCODE v32 gene TSS recorded in the
archived evaluation, with the usual window safeguards.

The reference map in `cima_ag2_reference_map.csv` records lineage matches, proxies
and exclusions. These references come from the AlphaGenome 2 catalog and may
cover broader cell populations than the marker-defined CIMA subtypes. Scooby's
baseline still uses the OneK1K columns in `cima_mapping.csv`. The resolved
coverage and mapping, chosen tracks, pilot scores and comparison are saved in
`ag2_out/cima`. Rerunning reuses cached requests. `--prepare-only` rebuilds the
pilot manifest and mapping offline once the live catalog has been saved.

The full CIMA runner scores the 42,713 associations with saved Scooby effects
across all 48 labels. It uses the observed 10x 5′ assay where available, other
observed whole-cell tracks when needed, and supported imputed tracks for
conditions without an observed single-cell expression track. For the three
CD8 cytotoxic labels and the megakaryocyte label, both models use the same
broader Scooby reference. These four labels are marked `shared_proxy` and have
a separate pooled sensitivity result.

```bash
python scooby_inference/evaluate_alphagenome2_cima_full.py
```

The full runner saves the live condition catalog and observed track metadata
in `ag2_out/cima` if they are missing. It writes status to `full_progress.json` and
checkpoints completed blocks in `full_shards/`. Rerun the same command to
resume after an API interruption. It writes `full_scores.csv` and
`full_comparison.csv` after processing the complete panel.
