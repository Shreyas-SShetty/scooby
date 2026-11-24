Predicted expression level (REF and ALT) and variant effect in all 73 CIMA columns, for 6 CIMA eQTLs (check K).
ScoobyVariant_levels.py = scooby_reproducibility/scripts/inference_eqtls/ScoobyVariant.py plus 2 lines that also
save the log(total+1) predictions for REF and ALT as layers log_ref / log_alt (changes_vs_ScoobyVariant.diff).
Run (from scooby_inference/, 2026-11-21, GPU 1):
  python eval_out_cima_cells/levels_6pairs/ScoobyVariant_levels.py -in_path eval_out_cima_cells/levels_6pairs/pairs.tsv \
     -out_path eval_out_cima_cells/levels_6pairs/outs_0_5.h5ad -config_path step2_config_cima_cells.yaml
(The output name must be outs_<first>_<last>.h5ad: the script reads the row range from it.)
Effects equal those of the full run (e.g. C19orf33 in cMono_CD14 -1.1675 in both).
