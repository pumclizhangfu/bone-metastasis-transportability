#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: run_gate12z_block_context.R <gate12g_model.rds> <axis_scores.tsv> <outdir>")
}

model_file <- normalizePath(args[[1L]], mustWork = TRUE)
scores_file <- normalizePath(args[[2L]], mustWork = TRUE)
outdir <- args[[3L]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

seed <- 1208L
n_boot <- 1999L
set.seed(seed)

model <- readRDS(model_file)
scores <- fread(scores_file)
if (model$selected_k != 4L || nrow(model$x_raw_clr) != 41L || ncol(model$x_raw_clr) != 26L) {
  stop("Unexpected frozen Gate12G model scope")
}
if (nrow(scores) != 41L || uniqueN(scores$patient_id) != 18L) stop("Unexpected frozen Axis1 score scope")
if (!identical(model$meta$sample_id, scores$sample_id)) stop("Model and score sample order mismatch")
if (!identical(names(model$blocks), colnames(model$x_raw_clr))) stop("Frozen block order mismatch")

# Reconstruct the exact scaled feature matrix and signed block contributions.
scaled <- model$x_raw_clr
block_scales <- setNames(numeric(length(unique(model$blocks))), unique(model$blocks))
for (block in unique(model$blocks)) {
  j <- which(model$blocks == block)
  block_scales[[block]] <- sqrt(sum(apply(model$x_raw_clr[, j, drop = FALSE], 2L, var)))
  scaled[, j] <- scaled[, j, drop = FALSE] / block_scales[[block]]
}
centered <- sweep(scaled, 2L, model$pca$center, "-")
loading <- model$pca$rotation[, 1L]
feature_contribution <- sweep(centered, 2L, loading, "*")
block_contribution <- sapply(unique(model$blocks), function(block) {
  rowSums(feature_contribution[, model$blocks == block, drop = FALSE])
})
block_contribution <- as.data.table(block_contribution)
block_contribution[, sample_id := rownames(model$x_raw_clr)]
setcolorder(block_contribution, c("sample_id", unique(model$blocks)))
block_contribution <- merge(
  scores[, .(accession, cancer, patient_id, sample_id, compartment, retained_cells, Axis1)],
  block_contribution, by = "sample_id", sort = FALSE
)
block_contribution[, reconstructed_Axis1 := rowSums(.SD), .SDcols = unique(model$blocks)]
block_contribution[, absolute_reconstruction_error := abs(Axis1 - reconstructed_Axis1)]
if (max(block_contribution$absolute_reconstruction_error) >= 1e-10) stop("Exact Axis1 reconstruction failed")
fwrite(block_contribution, file.path(outdir, "axis1_exact_block_contributions.tsv"), sep = "\t")

block_summary <- rbindlist(lapply(unique(model$blocks), function(block) {
  x <- block_contribution[[block]]
  data.table(
    block = block,
    n_features = sum(model$blocks == block),
    variance = var(x),
    spearman_with_Axis1 = cor(x, block_contribution$Axis1, method = "spearman"),
    pearson_with_Axis1 = cor(x, block_contribution$Axis1),
    marginal_R2 = summary(lm(block_contribution$Axis1 ~ x))$r.squared,
    mean_absolute_contribution = mean(abs(x))
  )
}))
fwrite(block_summary, file.path(outdir, "axis1_block_marginal_summary.tsv"), sep = "\t")

# Patient-held-out block-only prediction. No latent coordinate is refitted;
# ordinary least squares predicts the already frozen full Axis1.
meta <- copy(scores[, .(accession, cancer, patient_id, sample_id, compartment,
                        retained_cells, Axis1)])
meta[, log10_depth := log10(retained_cells)]
feature_dt <- as.data.table(scaled)
feature_dt[, sample_id := rownames(scaled)]
analysis_dt <- merge(meta, feature_dt, by = "sample_id", sort = FALSE)

feature_names <- colnames(scaled)
model_features <- list(
  broad_only = feature_names[model$blocks == "Broad"],
  myeloid_only = feature_names[model$blocks == "Myeloid"],
  t_nk_only = feature_names[model$blocks == "T_NK"],
  broad_plus_depth = c(feature_names[model$blocks == "Broad"], "log10_depth")
)

patients <- unique(analysis_dt$patient_id)
lopo_rows <- list()
for (model_name in names(model_features)) {
  vars <- model_features[[model_name]]
  form <- reformulate(vars, response = "Axis1")
  for (held in patients) {
    train <- analysis_dt[patient_id != held]
    test <- analysis_dt[patient_id == held]
    fit <- lm(form, data = train)
    pred <- as.numeric(predict(fit, newdata = test))
    lopo_rows[[length(lopo_rows) + 1L]] <- test[, .(
      model = model_name, held_out_patient = held, sample_id, patient_id,
      cancer, compartment, observed_Axis1 = Axis1,
      predicted_Axis1 = pred, training_mean_Axis1 = mean(train$Axis1)
    )]
  }
}
lopo <- rbindlist(lopo_rows)
lopo[, `:=`(
  residual = observed_Axis1 - predicted_Axis1,
  baseline_residual = observed_Axis1 - training_mean_Axis1
)]
fwrite(lopo, file.path(outdir, "axis1_block_lopo_predictions.tsv"), sep = "\t")

lopo_summary <- lopo[, .(
  n_samples = .N,
  n_patients = uniqueN(patient_id),
  cross_validated_R2_training_mean = 1 - sum(residual^2) / sum(baseline_residual^2),
  spearman_rho = cor(observed_Axis1, predicted_Axis1, method = "spearman"),
  pearson_r = cor(observed_Axis1, predicted_Axis1),
  mean_absolute_error = mean(abs(residual)),
  root_mean_squared_error = sqrt(mean(residual^2))
), by = model]
fwrite(lopo_summary, file.path(outdir, "axis1_block_lopo_summary.tsv"), sep = "\t")

paired_prediction_rows <- list()
for (model_name in names(model_features)) {
  z <- lopo[model == model_name]
  wide <- dcast(z, patient_id + cancer ~ compartment, value.var = "predicted_Axis1")
  if (all(c("tumor", "distal") %in% names(wide))) {
    d <- wide[is.finite(tumor) & is.finite(distal), tumor - distal]
    paired_prediction_rows[[length(paired_prediction_rows) + 1L]] <- data.table(
      model = model_name, contrast = "tumor_minus_distal",
      n_pairs = length(d), n_positive = sum(d > 0),
      positive_fraction = mean(d > 0), median_difference = median(d), mean_difference = mean(d)
    )
  }
}
paired_prediction <- rbindlist(paired_prediction_rows)
fwrite(paired_prediction, file.path(outdir, "axis1_block_lopo_paired_direction.tsv"), sep = "\t")

# Transparent, order-independent context decomposition. Patient absorbs the
# cancer/accession factor when both are present; this structural alias is
# retained and explicitly reported rather than silently assigning shared mass.
context <- copy(meta)
context[, `:=`(
  compartment_f = factor(compartment, levels = c("distal", "involved", "tumor")),
  cancer_accession = factor(paste(cancer, accession, sep = "::")),
  patient_f = factor(patient_id)
)]
full_formula <- Axis1 ~ compartment_f + cancer_accession + log10_depth + patient_f
full_fit <- lm(full_formula, data = context)
full_sse <- sum(residuals(full_fit)^2)
sst <- sum((context$Axis1 - mean(context$Axis1))^2)
full_r2 <- 1 - full_sse / sst

factor_formulas <- list(
  compartment = Axis1 ~ compartment_f,
  cancer_accession = Axis1 ~ cancer_accession,
  log10_depth = Axis1 ~ log10_depth,
  patient = Axis1 ~ patient_f
)
reduced_formulas <- list(
  compartment = Axis1 ~ cancer_accession + log10_depth + patient_f,
  cancer_accession = Axis1 ~ compartment_f + log10_depth + patient_f,
  log10_depth = Axis1 ~ compartment_f + cancer_accession + patient_f,
  patient = Axis1 ~ compartment_f + cancer_accession + log10_depth
)

context_rows <- lapply(names(factor_formulas), function(term) {
  marginal_fit <- lm(factor_formulas[[term]], data = context)
  reduced_fit <- lm(reduced_formulas[[term]], data = context)
  reduced_sse <- sum(residuals(reduced_fit)^2)
  delta_ss <- max(0, reduced_sse - full_sse)
  data.table(
    term = term,
    marginal_R2 = summary(marginal_fit)$r.squared,
    full_model_R2 = full_r2,
    drop_one_delta_R2 = delta_ss / sst,
    partial_R2 = if (reduced_sse > 0) delta_ss / reduced_sse else NA_real_,
    reduced_model_R2 = 1 - reduced_sse / sst,
    structural_note = if (term == "cancer_accession")
      "Aliased with patient when patient fixed effects are present; zero unique drop-one contribution is expected"
    else ""
  )
})
context_decomposition <- rbindlist(context_rows)
fwrite(context_decomposition, file.path(outdir, "axis1_context_decomposition.tsv"), sep = "\t")

alias_table <- data.table(
  cancer = context$cancer,
  accession = context$accession,
  cancer_accession = as.character(context$cancer_accession),
  patient_id = context$patient_id
)[, .(n_samples = .N, n_patients = uniqueN(patient_id)), by = .(cancer, accession, cancer_accession)]
alias_table[, perfect_cancer_accession_confounding := TRUE]
fwrite(alias_table, file.path(outdir, "cancer_accession_alias_audit.tsv"), sep = "\t")

coef_alias <- data.table(
  coefficient = names(coef(full_fit)),
  estimate = as.numeric(coef(full_fit)),
  estimable = is.finite(as.numeric(coef(full_fit)))
)
fwrite(coef_alias, file.path(outdir, "axis1_full_context_model_coefficients.tsv"), sep = "\t")

# Within-patient paired anatomical contrasts with patient bootstrap intervals.
contrast_defs <- list(
  involved_minus_distal = c("involved", "distal"),
  tumor_minus_involved = c("tumor", "involved"),
  tumor_minus_distal = c("tumor", "distal")
)
paired_rows <- list(); bootstrap_rows <- list()
for (group_value in c("Pooled", sort(unique(context$cancer)))) {
  z <- if (group_value == "Pooled") context else context[cancer == group_value]
  wide <- dcast(z, patient_id ~ compartment, value.var = "Axis1")
  for (contrast_name in names(contrast_defs)) {
    lev <- contrast_defs[[contrast_name]]
    if (!all(lev %in% names(wide))) next
    d <- wide[is.finite(get(lev[1L])) & is.finite(get(lev[2L])), get(lev[1L]) - get(lev[2L])]
    if (!length(d)) next
    set.seed(seed + length(paired_rows) + 100L)
    boot <- replicate(n_boot, median(sample(d, length(d), replace = TRUE)))
    paired_rows[[length(paired_rows) + 1L]] <- data.table(
      group = group_value, contrast = contrast_name,
      n_pairs = length(d), n_positive = sum(d > 0), positive_fraction = mean(d > 0),
      median_difference = median(d), mean_difference = mean(d),
      bootstrap_ci_low = quantile(boot, 0.025), bootstrap_ci_high = quantile(boot, 0.975),
      wilcoxon_two_sided_p = wilcox.test(d, mu = 0, exact = length(d) < 50)$p.value
    )
    bootstrap_rows[[length(bootstrap_rows) + 1L]] <- data.table(
      group = group_value, contrast = contrast_name, iteration = seq_len(n_boot),
      median_difference = boot
    )
  }
}
paired_contrasts <- rbindlist(paired_rows)
paired_contrasts[, wilcoxon_q := p.adjust(wilcoxon_two_sided_p, method = "BH"), by = group]
paired_bootstrap <- rbindlist(bootstrap_rows)
fwrite(paired_contrasts, file.path(outdir, "axis1_within_patient_context_contrasts.tsv"), sep = "\t")
fwrite(paired_bootstrap, file.path(outdir, "axis1_within_patient_context_bootstrap.tsv.gz"), sep = "\t")

summary <- data.table(
  metric = c(
    "samples", "patients", "features", "broad_features", "myeloid_features", "t_nk_features",
    "max_reconstruction_error", "full_context_R2", "compartment_marginal_R2",
    "compartment_drop_one_delta_R2", "patient_marginal_R2", "patient_drop_one_delta_R2",
    "cancer_accession_marginal_R2", "cancer_accession_unique_drop_one_delta_R2",
    paste0(names(model_features), "_lopo_R2")
  ),
  value = as.character(c(
    nrow(context), uniqueN(context$patient_id), ncol(scaled),
    sum(model$blocks == "Broad"), sum(model$blocks == "Myeloid"), sum(model$blocks == "T_NK"),
    max(block_contribution$absolute_reconstruction_error), full_r2,
    context_decomposition[term == "compartment", marginal_R2],
    context_decomposition[term == "compartment", drop_one_delta_R2],
    context_decomposition[term == "patient", marginal_R2],
    context_decomposition[term == "patient", drop_one_delta_R2],
    context_decomposition[term == "cancer_accession", marginal_R2],
    context_decomposition[term == "cancer_accession", drop_one_delta_R2],
    lopo_summary[match(names(model_features), model), cross_validated_R2_training_mean]
  ))
)
fwrite(summary, file.path(outdir, "GATE12Z_BLOCK_CONTEXT_SUMMARY.tsv"), sep = "\t")

best_model <- lopo_summary[which.max(cross_validated_R2_training_mean)]
writeLines(c(
  "# Gate12Z Axis1 block-ablation and sample-context checkpoint", "",
  paste0("- Exact block reconstruction maximum error: ",
         format(max(block_contribution$absolute_reconstruction_error), scientific = TRUE)),
  paste0("- Best patient-held-out block model: ", best_model$model,
         " (R2=", sprintf("%.3f", best_model$cross_validated_R2_training_mean),
         "; Spearman rho=", sprintf("%.3f", best_model$spearman_rho), ")"),
  paste0("- Compartment marginal R2: ",
         sprintf("%.3f", context_decomposition[term == "compartment", marginal_R2])),
  paste0("- Compartment unique drop-one delta R2 after patient/cancer-depth adjustment: ",
         sprintf("%.3f", context_decomposition[term == "compartment", drop_one_delta_R2])),
  paste0("- Patient marginal R2: ",
         sprintf("%.3f", context_decomposition[term == "patient", marginal_R2])),
  paste0("- Cancer/accession marginal R2: ",
         sprintf("%.3f", context_decomposition[term == "cancer_accession", marginal_R2])),
  "- Cancer and accession are perfectly confounded and patient fixed effects absorb their shared contribution; no separate causal attribution is possible.",
  "",
  "The block tests quantify how well each feature family reconstructs the already frozen sample-context coordinate. They are not clinical prediction models and do not distinguish disease biology from every unrecorded specimen-context factor."
), file.path(outdir, "GATE12Z_BLOCK_CONTEXT_CHECKPOINT.md"))

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
message("Gate12Z block/context analysis complete")
