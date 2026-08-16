#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(Cairo)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("Usage: finalize_gate12v2_axis1_interpretation.R <outdir> <gate12r_dir> <model.rds> <scores.tsv> <pseudobulk.rds>")
}
outdir <- normalizePath(args[[1L]], mustWork = TRUE)
gate12r_dir <- normalizePath(args[[2L]], mustWork = TRUE)
model_file <- normalizePath(args[[3L]], mustWork = TRUE)
score_file <- normalizePath(args[[4L]], mustWork = TRUE)
pseudobulk_file <- normalizePath(args[[5L]], mustWork = TRUE)

needed <- c(
  "axis1_feature_loading_stability.tsv", "axis1_bootstrap_summary.tsv",
  "hallmark_gene_set_audit.tsv", "hallmark_activity_scores.tsv.gz",
  "hallmark_primary_associations.tsv", "hallmark_cancer_specific_associations.tsv",
  "hallmark_lopo_associations.tsv.gz", "hallmark_interpretability_summary.tsv"
)
if (any(!file.exists(file.path(outdir, needed)))) stop("Core Gate12V2 outputs are incomplete")

loading_stability <- fread(file.path(outdir, "axis1_feature_loading_stability.tsv"))
bootstrap_summary <- fread(file.path(outdir, "axis1_bootstrap_summary.tsv"))
hallmark_audit <- fread(file.path(outdir, "hallmark_gene_set_audit.tsv"))
activity_scores <- fread(file.path(outdir, "hallmark_activity_scores.tsv.gz"), select = c(
  "pseudobulk_id", "sample_id", "patient_id", "cancer", "meta_state"
))
primary <- fread(file.path(outdir, "hallmark_primary_associations.tsv"))
cancer_specific <- fread(file.path(outdir, "hallmark_cancer_specific_associations.tsv"))
lopo <- fread(file.path(outdir, "hallmark_lopo_associations.tsv.gz"))
interpretability <- fread(file.path(outdir, "hallmark_interpretability_summary.tsv"))
model <- readRDS(model_file)
scores <- fread(score_file)
pb <- readRDS(pseudobulk_file)

stopifnot(
  nrow(loading_stability) == 26L,
  unique(loading_stability$valid_bootstraps) == 1999L,
  nrow(hallmark_audit) == 50L,
  nrow(primary) == 50L,
  nrow(cancer_specific) == 100L,
  nrow(lopo) == 50L * 18L,
  nrow(interpretability) == 50L
)

promoted_files <- c(
  axis1_block_contributions = "axis1_block_contributions.tsv",
  axis1_block_loading_mass = "axis1_block_loading_mass.tsv",
  axis1_categorical_paired_contrasts = "axis1_categorical_paired_contrasts.tsv",
  axis1_simple_mixture_lopo_predictions = "axis1_simple_mixture_lopo_predictions.tsv",
  axis1_simple_mixture_lopo_summary = "axis1_simple_mixture_lopo_summary.tsv",
  axis1_residual_anatomical_contrasts = "axis1_residual_anatomical_contrasts.tsv"
)
promoted_paths <- setNames(file.path(gate12r_dir, unname(promoted_files)), names(promoted_files))
if (any(!file.exists(promoted_paths))) stop("One or more promoted Gate12R files are missing")
promoted_manifest <- data.table(
  artifact = names(promoted_files),
  source_path = normalizePath(promoted_paths),
  bytes = file.info(promoted_paths)$size,
  md5 = unname(tools::md5sum(promoted_paths)),
  model_status = "promoted_unchanged"
)
fwrite(promoted_manifest, file.path(outdir, "promoted_gate12r_artifact_manifest.tsv"), sep = "\t")

block_mass <- fread(promoted_paths[["axis1_block_loading_mass"]])
categorical <- fread(promoted_paths[["axis1_categorical_paired_contrasts"]])
simple_mixture <- fread(promoted_paths[["axis1_simple_mixture_lopo_summary"]])
residual <- fread(promoted_paths[["axis1_residual_anatomical_contrasts"]])
block_contrib <- fread(promoted_paths[["axis1_block_contributions"]])

decision_status <- if (any(interpretability$interpretable_association)) {
  "INTERPRETABLE_ASSOCIATION"
} else {
  "COMPOSITION_ONLY_INTERPRETATION"
}
decision <- data.table(
  gate = "Gate12V2",
  execution_status = "COMPLETE",
  scientific_decision = decision_status,
  frozen_features = nrow(loading_stability),
  valid_bootstraps = unique(loading_stability$valid_bootstraps),
  hallmarks_tested = nrow(interpretability),
  hallmarks_fdr = sum(interpretability$fdr_pass),
  hallmarks_cross_cancer = sum(interpretability$fdr_pass & interpretability$cross_cancer_same_direction),
  hallmarks_interpretable = sum(interpretability$interpretable_association),
  primary_singular_models = sum(primary$singular, na.rm = TRUE),
  simple_mixture_cv_r2 = simple_mixture$cross_validated_R2,
  max_exact_reconstruction_error = max(block_contrib$absolute_reconstruction_error),
  interpretation_scope = "within_discovery_not_independent_validation"
)
fwrite(decision, file.path(outdir, "gate12v2_decision.tsv"), sep = "\t")

input_audit <- data.table(
  metric = c(
    "axis_samples", "axis_patients", "axis_cancers", "frozen_features", "broad_features",
    "myeloid_features", "t_nk_features", "bootstrap_target", "bootstrap_attempts",
    "pseudobulks_input", "pseudobulks_ncells_ge20", "pseudobulks_analyzed",
    "pseudobulk_states", "hallmarks_available", "hallmarks_tested", "seed"
  ),
  value = as.character(c(
    nrow(scores), uniqueN(scores$patient_id), uniqueN(scores$cancer), nrow(loading_stability),
    sum(model$blocks == "Broad"), sum(model$blocks == "Myeloid"), sum(model$blocks == "T_NK"),
    1999L, bootstrap_summary$total_attempts, nrow(pb$metadata), sum(pb$metadata$n_cells >= 20L),
    nrow(activity_scores), uniqueN(activity_scores$meta_state), nrow(hallmark_audit),
    nrow(interpretability), 12022026L
  ))
)
fwrite(input_audit, file.path(outdir, "input_audit.tsv"), sep = "\t")

loading_plot <- copy(loading_stability)
loading_plot[, feature_label := sub("^(Broad|Myeloid|T_NK)__", "", feature)]
loading_plot[, feature_label := factor(feature_label, levels = rev(feature_label))]
p1 <- ggplot(loading_plot, aes(frozen_loading, feature_label, colour = block)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.18, linewidth = 0.35) +
  geom_point(aes(size = sign_stability), alpha = 0.9) +
  facet_grid(block ~ ., scales = "free_y", space = "free_y") +
  scale_size_continuous(range = c(1.6, 3.4), limits = c(0, 1)) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom") +
  labs(title = "A  Frozen Axis1 loadings with patient-bootstrap intervals",
       x = "Frozen loading (95% aligned-bootstrap interval)", y = NULL,
       size = "Sign stability", colour = "Block")

p2 <- ggplot(block_mass, aes(block, squared_fraction, fill = block)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * squared_fraction)), vjust = -0.35, size = 3.1) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"),
                     limits = c(0, max(block_mass$squared_fraction) * 1.22)) +
  theme_bw(base_size = 9) +
  theme(legend.position = "none", panel.grid.minor = element_blank()) +
  labs(title = "B  Exact squared-loading mass by block",
       subtitle = sprintf("LOPO simple-mixture R² = %.3f; max reconstruction error = %.2e",
                          simple_mixture$cross_validated_R2, max(block_contrib$absolute_reconstruction_error)),
       x = NULL, y = "Fraction of Axis1 squared loading")

hallmark_plot <- copy(interpretability)
hallmark_plot[, pathway_label := gsub("_", " ", sub("^HALLMARK_", "", pathway))]
hallmark_plot[, pathway_label := factor(pathway_label, levels = rev(pathway_label[order(beta)]))]
p3 <- ggplot(hallmark_plot, aes(beta, pathway_label, colour = interpretable_association)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.12, linewidth = 0.25) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = c(`FALSE` = "grey35", `TRUE` = "#D55E00")) +
  theme_bw(base_size = 7.6) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom") +
  labs(title = "C  All 50 Hallmark associations after stage and cancer adjustment",
       x = "Activity SD per Axis1 SD (95% CI)", y = NULL, colour = "Passes all criteria")

robust <- interpretability[, .(
  pathway, Pooled = beta, Prostate = beta_prostate, Renal = beta_renal,
  LOPO = sign(beta) * lopo_direction_fraction
)]
robust <- melt(robust, id.vars = "pathway", variable.name = "analysis", value.name = "signed_metric")
robust[, pathway_label := gsub("_", " ", sub("^HALLMARK_", "", pathway))]
robust[, pathway_label := factor(pathway_label, levels = rev(unique(hallmark_plot$pathway_label)))]
p4 <- ggplot(robust, aes(analysis, pathway_label, fill = signed_metric)) +
  geom_tile(colour = "white", linewidth = 0.12) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  theme_bw(base_size = 7.6) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(title = "D  Cross-cancer direction and leave-one-patient-out stability",
       subtitle = "LOPO cells show signed direction-retention fraction",
       x = NULL, y = NULL, fill = "Signed metric")

fig <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title = "Axis1 biological-interpretability audit",
    subtitle = paste0("Decision: ", decision_status,
                      "; associations are internal to the discovery cohort and are not independent validation")
  )
CairoPDF(file.path(outdir, "Figure_gate12v2_axis1_interpretation.pdf"), width = 14, height = 17)
print(fig)
dev.off()
CairoPNG(file.path(outdir, "Figure_gate12v2_axis1_interpretation.png"), width = 4200, height = 5100, res = 300)
print(fig)
dev.off()

passing <- interpretability[interpretable_association == TRUE]
passing_lines <- if (nrow(passing)) {
  paste0("- ", passing$pathway, ": beta=", sprintf("%.3f", passing$beta),
         ", 95% CI ", sprintf("%.3f", passing$ci_low), " to ", sprintf("%.3f", passing$ci_high),
         ", q=", format(passing$q_value, digits = 3),
         ", prostate=", sprintf("%.3f", passing$beta_prostate),
         ", renal=", sprintf("%.3f", passing$beta_renal),
         ", LOPO=", sprintf("%.1f%%", 100 * passing$lopo_direction_fraction))
} else {
  "- None."
}
writeLines(c(
  "# Gate12V2 Axis1 biological-interpretability checkpoint", "",
  "## Execution status", "",
  "- Computational run: COMPLETE (independent audit pending).",
  paste0("- Scientific decision: **", decision_status, "**."),
  paste0("- Feature-loading bootstrap: 1,999 valid patient-cluster replicates; ",
         bootstrap_summary$rejected_attempts, " rejected attempts."),
  paste0("- Hallmark scope: ", nrow(interpretability), "/50 pathways tested; ",
         sum(interpretability$fdr_pass), " primary FDR hits; ", nrow(passing),
         " passed all frozen criteria."),
  paste0("- Primary mixed models with a boundary (singular) random-effect fit: ",
         sum(primary$singular, na.rm = TRUE), "/50. Fixed Axis1 estimates remain finite; this is retained as a limitation."),
  "", "## Hallmarks passing every frozen criterion", "", passing_lines,
  "", "## Promoted existing Axis1 interpretation", "",
  paste0("- Exact block squared-loading fractions: Broad ",
         sprintf("%.1f%%", 100 * block_mass[block == "Broad", squared_fraction]),
         ", Myeloid ", sprintf("%.1f%%", 100 * block_mass[block == "Myeloid", squared_fraction]),
         ", T/NK ", sprintf("%.1f%%", 100 * block_mass[block == "T_NK", squared_fraction]), "."),
  paste0("- Patient-held-out simple-mixture cross-validated R2: ",
         sprintf("%.3f", simple_mixture$cross_validated_R2), "."),
  paste0("- Pooled categorical tumour-minus-distal median: ",
         sprintf("%.3f", categorical[group == "Pooled" & contrast == "tumor_minus_distal", median_difference]), "."),
  paste0("- Mixture/depth/cancer-adjusted residual tumour-minus-distal median: ",
         sprintf("%.3f", residual[contrast == "tumor_minus_distal", median_difference]), "."),
  "", "## Interpretation boundary", "",
  "Hallmark activity and Axis1 were calculated from the same discovery samples. Surviving associations provide biological interpretation only; they are not independent validation, causal pathway evidence, prediction or a therapeutic target screen. Hallmark sets are overlapping and the count of qualifying pathways must not be presented as 33 independent mechanisms."
), file.path(outdir, "GATE12V2_CHECKPOINT.md"))

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
cat("GATE12V2_FINALIZATION_STATUS=COMPLETE\n")
print(decision)
