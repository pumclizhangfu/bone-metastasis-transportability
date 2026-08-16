#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(edgeR)
  library(msigdbr)
  library(lme4)
  library(lmerTest)
  library(matrixStats)
  library(ggplot2)
  library(patchwork)
  library(Cairo)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop(paste(
    "Usage: run_gate12v2_axis1_interpretation.R",
    "<gate12g_model.rds> <axis_scores.tsv> <state_pseudobulk_counts.rds>",
    "<gate12r_statistical_dir> <outdir> <n_boot> <seed>"
  ))
}

model_file <- normalizePath(args[[1L]], mustWork = TRUE)
score_file <- normalizePath(args[[2L]], mustWork = TRUE)
pseudobulk_file <- normalizePath(args[[3L]], mustWork = TRUE)
gate12r_dir <- normalizePath(args[[4L]], mustWork = TRUE)
outdir <- args[[5L]]
n_boot <- as.integer(args[[6L]])
seed <- as.integer(args[[7L]])
if (!is.finite(n_boot) || n_boot != 1999L) stop("Gate12V2 requires exactly 1,999 valid bootstraps")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

required <- function(x, cols, label) {
  missing <- setdiff(cols, names(x))
  if (length(missing)) stop(label, " missing columns: ", paste(missing, collapse = ", "))
}

zscore <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

block_scale <- function(x, blocks) {
  out <- x
  scales <- setNames(numeric(length(unique(blocks))), unique(blocks))
  for (block in unique(blocks)) {
    j <- which(blocks == block)
    total_variance <- sum(apply(x[, j, drop = FALSE], 2L, var))
    if (!is.finite(total_variance) || total_variance <= 0) stop("Non-positive variance in block ", block)
    scales[[block]] <- sqrt(total_variance)
    out[, j] <- out[, j, drop = FALSE] / scales[[block]]
  }
  list(x = out, scales = scales)
}

fit_axis_term <- function(dat, include_cancer = TRUE) {
  fixed <- if (include_cancer) {
    "Axis1_z + compartment_f + cancer_f + meta_state_f"
  } else {
    "Axis1_z + compartment_f + meta_state_f"
  }
  form <- as.formula(paste0("activity_z ~ ", fixed, " + (1 | patient_id/sample_id)"))
  warnings_seen <- character()
  fit <- tryCatch(
    withCallingHandlers(
      lmerTest::lmer(
        form, data = dat, REML = FALSE,
        control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000),
                              check.conv.singular = "ignore")
      ),
      warning = function(w) {
        warnings_seen <<- c(warnings_seen, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(data.table(
      beta = NA_real_, se = NA_real_, df = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      statistic = NA_real_, p_value = NA_real_, singular = NA, converged = FALSE,
      warning = conditionMessage(fit)
    ))
  }
  tab <- coef(summary(fit))
  if (!"Axis1_z" %in% rownames(tab)) {
    return(data.table(
      beta = NA_real_, se = NA_real_, df = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      statistic = NA_real_, p_value = NA_real_, singular = lme4::isSingular(fit), converged = FALSE,
      warning = paste(c(warnings_seen, "Axis1_z coefficient unavailable"), collapse = " | ")
    ))
  }
  row <- tab["Axis1_z", ]
  df_value <- unname(row[["df"]])
  critical <- if (is.finite(df_value)) qt(0.975, df = df_value) else 1.96
  estimate <- unname(row[["Estimate"]])
  standard_error <- unname(row[["Std. Error"]])
  data.table(
    beta = estimate,
    se = standard_error,
    df = df_value,
    ci_low = estimate - critical * standard_error,
    ci_high = estimate + critical * standard_error,
    statistic = unname(row[["t value"]]),
    p_value = unname(row[["Pr(>|t|)"]]),
    singular = lme4::isSingular(fit, tol = 1e-5),
    converged = is.finite(estimate) && is.finite(standard_error),
    warning = paste(unique(warnings_seen), collapse = " | ")
  )
}

model <- readRDS(model_file)
scores <- fread(score_file)
pb <- readRDS(pseudobulk_file)
required(scores, c("accession", "cancer", "patient_id", "sample_id", "compartment", "Axis1"), "Axis scores")
required(as.data.table(pb$metadata), c("pseudobulk_id", "cancer", "patient_id", "sample_id",
                                      "compartment", "meta_state", "n_cells"), "Pseudobulk metadata")
if (!identical(rownames(model$x_raw_clr), scores$sample_id)) stop("Model and score sample orders differ")
if (!identical(names(model$blocks), colnames(model$x_raw_clr))) stop("Model block and feature orders differ")
if (ncol(model$x_raw_clr) != 26L) stop("Expected 26 frozen features")
if (nrow(scores) != 41L || uniqueN(scores$patient_id) != 18L) stop("Unexpected Axis1 discovery scope")
if (!identical(colnames(pb$counts), pb$metadata$pseudobulk_id)) stop("Pseudobulk matrix and metadata differ")

# -----------------------------------------------------------------------------
# 1. Full feature-level patient-cluster bootstrap for the frozen Axis1.
# -----------------------------------------------------------------------------
set.seed(seed)
reference_loading <- model$pca$rotation[, 1L]
patients <- unique(model$meta$patient_id)
bootstrap_rows <- vector("list", n_boot)
match_rows <- vector("list", n_boot)
accepted <- 0L
attempted <- 0L
max_attempts <- n_boot * 10L
while (accepted < n_boot && attempted < max_attempts) {
  attempted <- attempted + 1L
  sampled_patients <- sample(patients, length(patients), replace = TRUE)
  idx <- unlist(lapply(sampled_patients, function(pid) which(model$meta$patient_id == pid)), use.names = FALSE)
  fit <- try({
    scaled <- block_scale(model$x_raw_clr[idx, , drop = FALSE], model$blocks)$x
    prcomp(scaled, center = TRUE, scale. = FALSE, rank. = 4L)
  }, silent = TRUE)
  if (inherits(fit, "try-error") || ncol(fit$rotation) < 4L) next
  candidate_cor <- vapply(seq_len(4L), function(k) cor(reference_loading, fit$rotation[, k]), numeric(1))
  if (!any(is.finite(candidate_cor))) next
  matched_axis <- which.max(abs(candidate_cor))
  alignment_cor <- candidate_cor[[matched_axis]]
  aligned <- fit$rotation[, matched_axis] * sign(alignment_cor)
  if (any(!is.finite(aligned))) next
  accepted <- accepted + 1L
  bootstrap_rows[[accepted]] <- data.table(
    iteration = accepted,
    feature = names(reference_loading),
    block = unname(model$blocks[names(reference_loading)]),
    frozen_loading = as.numeric(reference_loading),
    aligned_loading = as.numeric(aligned)
  )
  match_rows[[accepted]] <- data.table(
    iteration = accepted,
    draw_attempt = attempted,
    matched_axis = matched_axis,
    preflip_loading_correlation = alignment_cor,
    aligned_loading_correlation = abs(alignment_cor),
    unique_patients_in_draw = uniqueN(sampled_patients)
  )
}
if (accepted != n_boot) stop("Unable to obtain 1,999 valid patient-cluster bootstrap replicates")
bootstrap_loadings <- rbindlist(bootstrap_rows)
bootstrap_match <- rbindlist(match_rows)
loading_stability <- bootstrap_loadings[, .(
  frozen_loading = unique(frozen_loading),
  bootstrap_median = median(aligned_loading),
  ci_low = quantile(aligned_loading, 0.025),
  ci_high = quantile(aligned_loading, 0.975),
  sign_stability = mean(sign(aligned_loading) == sign(unique(frozen_loading))),
  interval_excludes_zero = quantile(aligned_loading, 0.025) * quantile(aligned_loading, 0.975) > 0,
  valid_bootstraps = .N
), by = .(feature, block)]
loading_stability[, abs_frozen_loading := abs(frozen_loading)]
setorder(loading_stability, block, -abs_frozen_loading)
bootstrap_match_summary <- bootstrap_match[, .(
  n_bootstraps = .N,
  median_aligned_correlation = median(aligned_loading_correlation),
  q025_aligned_correlation = quantile(aligned_loading_correlation, 0.025),
  q975_aligned_correlation = quantile(aligned_loading_correlation, 0.975),
  median_unique_patients = median(unique_patients_in_draw),
  total_attempts = max(draw_attempt),
  rejected_attempts = max(draw_attempt) - .N
)]
matched_axis_frequency <- bootstrap_match[, .N, by = matched_axis][, fraction := N / sum(N)]

fwrite(bootstrap_loadings, file.path(outdir, "axis1_bootstrap_feature_loadings.tsv.gz"), sep = "\t")
fwrite(bootstrap_match, file.path(outdir, "axis1_bootstrap_axis_matches.tsv.gz"), sep = "\t")
fwrite(loading_stability, file.path(outdir, "axis1_feature_loading_stability.tsv"), sep = "\t")
fwrite(bootstrap_match_summary, file.path(outdir, "axis1_bootstrap_summary.tsv"), sep = "\t")
fwrite(matched_axis_frequency, file.path(outdir, "axis1_matched_axis_frequency.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# 2. Hallmark activities from exact patient-sample-state pseudobulk counts.
# -----------------------------------------------------------------------------
pb_meta <- as.data.table(pb$metadata)
keep_pb <- pb_meta$n_cells >= 20L & pb_meta$sample_id %chin% scores$sample_id
pb_meta <- copy(pb_meta[keep_pb])
counts <- pb$counts[, keep_pb, drop = FALSE]
if (ncol(counts) != nrow(pb_meta) || nrow(pb_meta) < 200L) stop("Unexpected retained pseudobulk scope")

dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge, method = "TMM")
log_cpm <- cpm(dge, log = TRUE, prior.count = 2)
gene_mean <- rowMeans(log_cpm)
gene_sd <- matrixStats::rowSds(log_cpm)
valid_gene <- is.finite(gene_mean) & is.finite(gene_sd) & gene_sd > 0
log_cpm <- log_cpm[valid_gene, , drop = FALSE]
gene_mean <- gene_mean[valid_gene]
gene_sd <- gene_sd[valid_gene]
gene_z <- sweep(sweep(log_cpm, 1L, gene_mean, "-"), 1L, gene_sd, "/")
upper_gene <- toupper(rownames(gene_z))

msig <- as.data.table(msigdbr(db_species = "HS", species = "human", collection = "H"))
hallmark_sets <- lapply(split(msig$gene_symbol, msig$gs_name), function(x) unique(toupper(x)))
if (length(hallmark_sets) != 50L) stop("Expected all 50 Hallmark gene sets")

hallmark_audit_rows <- vector("list", length(hallmark_sets))
activity_matrix <- matrix(NA_real_, nrow = ncol(gene_z), ncol = length(hallmark_sets),
                          dimnames = list(colnames(gene_z), names(hallmark_sets)))
for (i in seq_along(hallmark_sets)) {
  pathway <- names(hallmark_sets)[[i]]
  requested <- hallmark_sets[[i]]
  idx <- match(requested, upper_gene)
  idx <- unique(idx[is.finite(idx)])
  hallmark_audit_rows[[i]] <- data.table(
    pathway = pathway,
    genes_in_set = length(requested),
    measured_nonconstant_genes = length(idx),
    coverage = length(idx) / length(requested),
    eligible = length(idx) >= 15L
  )
  if (length(idx) >= 15L) activity_matrix[, i] <- colMeans(gene_z[idx, , drop = FALSE])
}
hallmark_audit <- rbindlist(hallmark_audit_rows)
if (sum(hallmark_audit$eligible) != 50L) stop("Not all 50 Hallmarks have adequate measured-gene coverage")

activity <- as.data.table(activity_matrix, keep.rownames = "pseudobulk_id")
activity_data <- merge(pb_meta, activity, by = "pseudobulk_id", sort = FALSE)
activity_data <- merge(activity_data, scores[, .(sample_id, frozen_Axis1 = Axis1)], by = "sample_id", sort = FALSE)
if (nrow(activity_data) != nrow(pb_meta)) stop("Axis1 merge lost pseudobulks")
activity_data[, Axis1_z := zscore(frozen_Axis1)]
activity_data[, `:=`(
  compartment_f = factor(compartment, levels = c("distal", "involved", "tumor")),
  cancer_f = factor(cancer),
  meta_state_f = factor(meta_state),
  patient_id = factor(patient_id),
  sample_id = factor(sample_id)
)]

pathways <- names(hallmark_sets)
for (pathway in pathways) activity_data[, paste0(pathway, "__z") := zscore(get(pathway))]

score_output <- activity_data[, c(
  "pseudobulk_id", "accession", "cancer", "sample_id", "patient_id", "compartment",
  "lineage", "meta_state", "n_cells", "raw_umi_sum", "frozen_Axis1", "Axis1_z", pathways
), with = FALSE]
fwrite(score_output, file.path(outdir, "hallmark_activity_scores.tsv.gz"), sep = "\t")
fwrite(hallmark_audit, file.path(outdir, "hallmark_gene_set_audit.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# 3. Primary, cancer-specific and leave-one-patient-out association models.
# -----------------------------------------------------------------------------
primary_rows <- vector("list", length(pathways))
cancer_rows <- list()
lopo_rows <- list()
all_patients <- levels(activity_data$patient_id)
for (i in seq_along(pathways)) {
  pathway <- pathways[[i]]
  dat <- copy(activity_data)
  dat[, activity_z := get(paste0(pathway, "__z"))]
  dat <- dat[is.finite(activity_z) & is.finite(Axis1_z)]
  ans <- fit_axis_term(dat, include_cancer = TRUE)
  ans[, `:=`(
    pathway = pathway,
    n_pseudobulks = nrow(dat),
    n_samples = uniqueN(dat$sample_id),
    n_patients = uniqueN(dat$patient_id),
    model = "primary_adjusted_nested_random_intercept"
  )]
  primary_rows[[i]] <- ans

  for (cancer_value in levels(activity_data$cancer_f)) {
    sub <- droplevels(dat[cancer == cancer_value])
    ca <- fit_axis_term(sub, include_cancer = FALSE)
    ca[, `:=`(
      pathway = pathway,
      cancer = cancer_value,
      n_pseudobulks = nrow(sub),
      n_samples = uniqueN(sub$sample_id),
      n_patients = uniqueN(sub$patient_id),
      model = "cancer_specific_adjusted_nested_random_intercept"
    )]
    cancer_rows[[length(cancer_rows) + 1L]] <- ca
  }

  for (patient_value in all_patients) {
    sub <- droplevels(dat[as.character(patient_id) != patient_value])
    lp <- fit_axis_term(sub, include_cancer = TRUE)
    lp[, `:=`(
      pathway = pathway,
      deleted_patient = patient_value,
      n_pseudobulks = nrow(sub),
      n_samples = uniqueN(sub$sample_id),
      n_patients = uniqueN(sub$patient_id),
      model = "leave_one_patient_out"
    )]
    lopo_rows[[length(lopo_rows) + 1L]] <- lp
  }
}

primary <- rbindlist(primary_rows, fill = TRUE)
setcolorder(primary, c("pathway", setdiff(names(primary), "pathway")))
primary[, q_value := p.adjust(p_value, method = "BH")]
cancer_specific <- rbindlist(cancer_rows, fill = TRUE)
lopo <- rbindlist(lopo_rows, fill = TRUE)
primary_sign <- setNames(sign(primary$beta), primary$pathway)
lopo[, primary_direction := unname(primary_sign[pathway])]

cancer_wide <- dcast(cancer_specific, pathway ~ cancer, value.var = "beta")
setnames(cancer_wide, setdiff(names(cancer_wide), "pathway"),
         paste0("beta_", setdiff(names(cancer_wide), "pathway")))
lopo_summary <- lopo[, .(
  n_lopo_evaluable = sum(converged & is.finite(beta)),
  n_lopo_same_direction = sum(converged & is.finite(beta) & sign(beta) == primary_direction),
  lopo_direction_fraction = mean(sign(beta[converged & is.finite(beta)]) ==
                                   primary_direction[converged & is.finite(beta)]),
  lopo_beta_min = min(beta[converged & is.finite(beta)]),
  lopo_beta_max = max(beta[converged & is.finite(beta)]),
  lopo_singular_fraction = mean(singular[converged & !is.na(singular)])
), by = pathway]

interpretability <- merge(primary, cancer_wide, by = "pathway", all.x = TRUE)
interpretability <- merge(interpretability, lopo_summary, by = "pathway", all.x = TRUE)
if (!all(c("beta_prostate", "beta_renal") %in% names(interpretability))) stop("Missing cancer-specific coefficients")
interpretability[, `:=`(
  fdr_pass = is.finite(q_value) & q_value < 0.05,
  cross_cancer_same_direction = is.finite(beta_prostate) & is.finite(beta_renal) &
    sign(beta_prostate) == sign(beta) & sign(beta_renal) == sign(beta) & sign(beta) != 0,
  lopo_pass = n_lopo_evaluable == length(all_patients) & lopo_direction_fraction >= 0.90
)]
interpretability[, interpretable_association := fdr_pass & cross_cancer_same_direction & lopo_pass]
setorder(interpretability, q_value, pathway)

fwrite(primary, file.path(outdir, "hallmark_primary_associations.tsv"), sep = "\t")
fwrite(cancer_specific, file.path(outdir, "hallmark_cancer_specific_associations.tsv"), sep = "\t")
fwrite(lopo, file.path(outdir, "hallmark_lopo_associations.tsv.gz"), sep = "\t")
fwrite(interpretability, file.path(outdir, "hallmark_interpretability_summary.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# 4. Promote existing exact interpretation outputs without modifying their models.
# -----------------------------------------------------------------------------
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
    n_boot, attempted, nrow(pb$metadata), sum(pb$metadata$n_cells >= 20L), nrow(activity_data),
    uniqueN(activity_data$meta_state), length(hallmark_sets), nrow(interpretability), seed
  ))
)
fwrite(input_audit, file.path(outdir, "input_audit.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# 5. Diagnostic figure for scientific review. Final publication layout is Gate12V4.
# -----------------------------------------------------------------------------
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
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, max(block_mass$squared_fraction) * 1.22)) +
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
  pathway,
  Pooled = beta,
  Prostate = beta_prostate,
  Renal = beta_renal,
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
  paste0("- Feature-loading bootstrap: ", n_boot, " valid patient-cluster replicates; ",
         attempted - n_boot, " rejected attempts."),
  paste0("- Hallmark scope: ", nrow(interpretability), "/50 pathways tested; ",
         sum(interpretability$fdr_pass), " primary FDR hits; ", nrow(passing), " passed all frozen criteria."),
  "", "## Hallmarks passing every frozen criterion", "", passing_lines,
  "", "## Promoted existing Axis1 interpretation", "",
  paste0("- Exact block squared-loading fractions: Broad ",
         sprintf("%.1f%%", 100 * block_mass[block == "Broad", squared_fraction]),
         ", Myeloid ", sprintf("%.1f%%", 100 * block_mass[block == "Myeloid", squared_fraction]),
         ", T/NK ", sprintf("%.1f%%", 100 * block_mass[block == "T_NK", squared_fraction]), "."),
  paste0("- Patient-held-out simple-mixture cross-validated R2: ", sprintf("%.3f", simple_mixture$cross_validated_R2), "."),
  paste0("- Pooled categorical tumour-minus-distal median: ",
         sprintf("%.3f", categorical[group == "Pooled" & contrast == "tumor_minus_distal", median_difference]), "."),
  paste0("- Mixture/depth/cancer-adjusted residual tumour-minus-distal median: ",
         sprintf("%.3f", residual[contrast == "tumor_minus_distal", median_difference]), "."),
  "", "## Interpretation boundary", "",
  "Hallmark activity and Axis1 were calculated from the same discovery samples. Surviving associations provide biological interpretation only; they are not independent validation, causal pathway evidence, prediction or a therapeutic target screen."
), file.path(outdir, "GATE12V2_CHECKPOINT.md"))

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
cat("GATE12V2_EXECUTION_STATUS=COMPLETE\n")
print(decision)
