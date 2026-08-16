#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(metafor)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) {
  stop(paste(
    "Usage: run_gate12r_statistical_sensitivities.R",
    "<broad.tsv> <states.tsv> <state_cohort_effects.tsv> <model.rds>",
    "<axis_scores.tsv> <external_scores.tsv> <external_assignments.tsv.gz> <outdir>"
  ))
}

broad_file <- args[[1L]]
states_file <- args[[2L]]
cohort_effect_file <- args[[3L]]
model_file <- args[[4L]]
score_file <- args[[5L]]
external_score_file <- args[[6L]]
external_assignment_file <- args[[7L]]
outdir <- args[[8L]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

seed <- 1208L
set.seed(seed)

broad <- fread(broad_file)
states <- fread(states_file)
cohort_effects <- fread(cohort_effect_file)
model <- readRDS(model_file)
scores <- fread(score_file)
external_scores <- fread(external_score_file)
external_assignments <- fread(external_assignment_file)

stopifnot(
  all(c("sample_id", "patient_id", "cancer", "compartment", "Axis1") %in% names(scores)),
  identical(rownames(model$x_raw_clr), scores$sample_id),
  identical(names(model$blocks), colnames(model$x_raw_clr))
)

bootstrap_ci <- function(x, statistic, n_boot = 9999L, seed_offset = 0L) {
  x <- as.data.table(x)
  if (!nrow(x)) return(c(low = NA_real_, high = NA_real_))
  set.seed(seed + seed_offset)
  value <- replicate(n_boot, statistic(x[sample.int(.N, .N, replace = TRUE)]))
  quantile(value[is.finite(value)], c(0.025, 0.975), names = FALSE, na.rm = TRUE)
}

# 1. Patient-paired categorical anatomical contrasts for the frozen coordinate.
contrast_def <- list(
  involved_minus_distal = c("involved", "distal"),
  tumor_minus_involved = c("tumor", "involved"),
  tumor_minus_distal = c("tumor", "distal")
)

axis_contrast_rows <- list()
group_specs <- c(unique(scores$cancer), "Pooled")
for (group_name in group_specs) {
  z <- if (group_name == "Pooled") copy(scores) else scores[cancer == group_name]
  wide <- dcast(z, patient_id + cancer ~ compartment, value.var = "Axis1")
  for (contrast_name in names(contrast_def)) {
    lev <- contrast_def[[contrast_name]]
    if (!all(lev %in% names(wide))) next
    d <- wide[[lev[[1L]]]] - wide[[lev[[2L]]]]
    d <- d[is.finite(d)]
    ci <- if (length(d)) bootstrap_ci(data.table(d = d), function(q) median(q$d),
                                      seed_offset = length(axis_contrast_rows) + 1L) else c(NA, NA)
    wt <- if (length(d) >= 2L && any(d != 0)) suppressWarnings(wilcox.test(d, mu = 0, exact = FALSE)) else NULL
    axis_contrast_rows[[length(axis_contrast_rows) + 1L]] <- data.table(
      group = group_name,
      contrast = contrast_name,
      n_pairs = length(d),
      n_positive = sum(d > 0),
      median_difference = if (length(d)) median(d) else NA_real_,
      bootstrap_ci_low = ci[[1L]],
      bootstrap_ci_high = ci[[2L]],
      wilcoxon_p = if (is.null(wt)) NA_real_ else wt$p.value
    )
  }
}
axis_contrasts <- rbindlist(axis_contrast_rows)
axis_contrasts[, wilcoxon_q := p.adjust(wilcoxon_p, method = "BH")]
fwrite(axis_contrasts, file.path(outdir, "axis1_categorical_paired_contrasts.tsv"), sep = "\t")

# 2. Count-aware patient-fixed-effect state models, using categorical compartments.
state_model_rows <- list()
for (ca in unique(states$cancer)) {
  ca_states <- states[cancer == ca]
  complete_patients <- ca_states[, .(n_compartments = uniqueN(compartment)),
                                 by = patient_id][n_compartments == 3L, patient_id]
  ca_states <- ca_states[patient_id %chin% complete_patients]
  for (lin in unique(ca_states$lineage)) {
    for (st in unique(ca_states[lineage == lin, state])) {
      z <- ca_states[lineage == lin & state == st]
      z[, compartment_f := factor(compartment, levels = c("distal", "involved", "tumor"))]
      fit <- try(glm(cbind(n_state, lineage_total - n_state) ~ factor(patient_id) + compartment_f,
                     data = z, family = quasibinomial()), silent = TRUE)
      for (term in c("compartment_finvolved", "compartment_ftumor")) {
        label <- if (term == "compartment_finvolved") "involved_vs_distal" else "tumor_vs_distal"
        if (inherits(fit, "try-error") || !term %in% rownames(coef(summary(fit)))) {
          state_model_rows[[length(state_model_rows) + 1L]] <- data.table(
            cancer = ca, lineage = lin, state = st, contrast = label,
            n_patients = uniqueN(z$patient_id), beta = NA_real_, se = NA_real_,
            odds_ratio = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_
          )
        } else {
          tab <- coef(summary(fit))[term, ]
          critical <- qt(0.975, df = max(1, fit$df.residual))
          state_model_rows[[length(state_model_rows) + 1L]] <- data.table(
            cancer = ca, lineage = lin, state = st, contrast = label,
            n_patients = uniqueN(z$patient_id), beta = unname(tab[["Estimate"]]),
            se = unname(tab[["Std. Error"]]), odds_ratio = exp(unname(tab[["Estimate"]])),
            ci_low = exp(unname(tab[["Estimate"]]) - critical * unname(tab[["Std. Error"]])),
            ci_high = exp(unname(tab[["Estimate"]]) + critical * unname(tab[["Std. Error"]])),
            p_value = unname(tab[["Pr(>|t|)"]])
          )
        }
      }
    }
  }
}
state_models <- rbindlist(state_model_rows)
state_models[, q_value := p.adjust(p_value, method = "BH")]
fwrite(state_models, file.path(outdir, "state_patient_fixed_categorical_quasibinomial.tsv"), sep = "\t")

# 3. REML with Hartung-Knapp uncertainty. With two cohorts these remain exploratory.
meta_rows <- list()
for (lin in unique(cohort_effects$lineage)) {
  for (st in unique(cohort_effects[lineage == lin, state])) {
    z <- cohort_effects[lineage == lin & state == st & is.finite(beta_per_step) & is.finite(se) & se > 0]
    if (nrow(z) < 2L) next
    fit <- try(rma(yi = beta_per_step, sei = se, data = z, method = "REML", test = "knha"), silent = TRUE)
    if (inherits(fit, "try-error")) next
    meta_rows[[length(meta_rows) + 1L]] <- data.table(
      lineage = lin, state = st, k_cohorts = nrow(z),
      beta_reml_knha = as.numeric(fit$b), se_reml_knha = fit$se,
      ci_low = fit$ci.lb, ci_high = fit$ci.ub, p_value = fit$pval,
      tau2 = fit$tau2, I2 = fit$I2,
      direction_concordant = uniqueN(sign(z$beta_per_step)) == 1L,
      interpretation = "exploratory_k_equals_2"
    )
  }
}
meta_sensitivity <- rbindlist(meta_rows)
meta_sensitivity[, q_value := p.adjust(p_value, method = "BH")]
fwrite(meta_sensitivity, file.path(outdir, "state_reml_hartung_knapp_sensitivity.tsv"), sep = "\t")

# 4. Exact Axis1 feature and block decomposition.
block_scales <- setNames(numeric(length(unique(model$blocks))), unique(model$blocks))
scaled <- model$x_raw_clr
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
block_contribution <- merge(scores[, .(sample_id, patient_id, cancer, compartment, Axis1)],
                            block_contribution, by = "sample_id", sort = FALSE)
block_contribution[, reconstructed_Axis1 := rowSums(.SD), .SDcols = unique(model$blocks)]
block_contribution[, absolute_reconstruction_error := abs(Axis1 - reconstructed_Axis1)]
fwrite(block_contribution, file.path(outdir, "axis1_block_contributions.tsv"), sep = "\t")

loading_audit <- data.table(
  feature = names(loading), block = unname(model$blocks[names(loading)]), loading = as.numeric(loading)
)
loading_audit[, `:=`(abs_loading = abs(loading), squared_loading = loading^2)]
loading_mass <- loading_audit[, .(
  n_features = .N,
  absolute_loading_mass = sum(abs_loading),
  squared_loading_mass = sum(squared_loading)
), by = block]
loading_mass[, `:=`(
  absolute_fraction = absolute_loading_mass / sum(absolute_loading_mass),
  squared_fraction = squared_loading_mass / sum(squared_loading_mass)
)]
fwrite(loading_audit, file.path(outdir, "axis1_feature_loading_audit.tsv"), sep = "\t")
fwrite(loading_mass, file.path(outdir, "axis1_block_loading_mass.tsv"), sep = "\t")

# 5. Can a simple broad-mixture/depth model explain Axis1?
broad_wide <- dcast(broad[sample_id %chin% scores$sample_id],
                    sample_id ~ broad_class, value.var = "fraction", fill = 0)
baseline <- merge(scores, broad_wide, by = "sample_id", sort = FALSE)
for (nm in c("Malignant", "T_NK", "Myeloid")) if (!nm %in% names(baseline)) baseline[, (nm) := 0]
baseline[, log10_cells := log10(retained_cells)]

lopo_rows <- list()
for (pid in unique(baseline$patient_id)) {
  train <- baseline[patient_id != pid]
  test <- baseline[patient_id == pid]
  fit <- lm(Axis1 ~ Malignant + T_NK + Myeloid + log10_cells, data = train)
  lopo_rows[[length(lopo_rows) + 1L]] <- test[, .(
    sample_id, patient_id, cancer, compartment, Axis1,
    predicted_simple_mixture = as.numeric(predict(fit, newdata = test))
  )]
}
lopo <- rbindlist(lopo_rows)
lopo[, residual_lopo := Axis1 - predicted_simple_mixture]
fwrite(lopo, file.path(outdir, "axis1_simple_mixture_lopo_predictions.tsv"), sep = "\t")

lopo_summary <- data.table(
  n_samples = nrow(lopo), n_patients = uniqueN(lopo$patient_id),
  spearman_rho = cor(lopo$Axis1, lopo$predicted_simple_mixture, method = "spearman"),
  pearson_r = cor(lopo$Axis1, lopo$predicted_simple_mixture),
  cross_validated_R2 = 1 - sum(lopo$residual_lopo^2) / sum((lopo$Axis1 - mean(lopo$Axis1))^2),
  RMSE = sqrt(mean(lopo$residual_lopo^2))
)
fwrite(lopo_summary, file.path(outdir, "axis1_simple_mixture_lopo_summary.tsv"), sep = "\t")

full_baseline_fit <- lm(Axis1 ~ Malignant + T_NK + Myeloid + log10_cells + factor(cancer), data = baseline)
baseline[, residual_adjusted := residuals(full_baseline_fit)]
residual_wide <- dcast(baseline, patient_id + cancer ~ compartment, value.var = "residual_adjusted")
residual_rows <- list()
for (contrast_name in names(contrast_def)) {
  lev <- contrast_def[[contrast_name]]
  if (!all(lev %in% names(residual_wide))) next
  d <- residual_wide[[lev[[1L]]]] - residual_wide[[lev[[2L]]]]
  d <- d[is.finite(d)]
  ci <- bootstrap_ci(data.table(d = d), function(q) median(q$d), seed_offset = 500L + length(residual_rows))
  residual_rows[[length(residual_rows) + 1L]] <- data.table(
    contrast = contrast_name, n_pairs = length(d), n_positive = sum(d > 0),
    median_difference = median(d), bootstrap_ci_low = ci[[1L]], bootstrap_ci_high = ci[[2L]]
  )
}
residual_contrasts <- rbindlist(residual_rows)
fwrite(residual_contrasts, file.path(outdir, "axis1_residual_anatomical_contrasts.tsv"), sep = "\t")

# 6. External technical-dependence uncertainty, sensitivity to decision boundary,
#    and a frozen discovery broad-mixture projection diagnostic.
ext <- external_scores[projectable == TRUE]
boot_rho <- function(x, y, n_boot = 9999L, offset = 0L) {
  set.seed(seed + offset)
  v <- replicate(n_boot, {
    idx <- sample.int(length(x), length(x), replace = TRUE)
    suppressWarnings(cor(x[idx], y[idx], method = "spearman"))
  })
  quantile(v[is.finite(v)], c(0.025, 0.975), names = FALSE)
}

tech <- rbindlist(list(
  data.table(metric = "log10_qc_cells", rho = cor(ext$Axis1, log10(ext$qc_cells), method = "spearman"),
             ci_low = boot_rho(ext$Axis1, log10(ext$qc_cells), offset = 601L)[1L],
             ci_high = boot_rho(ext$Axis1, log10(ext$qc_cells), offset = 602L)[2L]),
  data.table(metric = "broad_assigned_fraction", rho = cor(ext$Axis1, ext$broad_assigned_fraction, method = "spearman"),
             ci_low = boot_rho(ext$Axis1, ext$broad_assigned_fraction, offset = 603L)[1L],
             ci_high = boot_rho(ext$Axis1, ext$broad_assigned_fraction, offset = 604L)[2L])
))
fwrite(tech, file.path(outdir, "external_axis1_technical_bootstrap.tsv"), sep = "\t")

threshold_sweep <- CJ(metric = tech$metric, absolute_rho_threshold = seq(0.30, 0.60, by = 0.025))
threshold_sweep <- merge(threshold_sweep, tech[, .(metric, rho)], by = "metric")
threshold_sweep[, passes := abs(rho) < absolute_rho_threshold]
fwrite(threshold_sweep, file.path(outdir, "external_axis1_technical_threshold_sweep.tsv"), sep = "\t")

ext_broad <- external_assignments[sample_id %chin% ext$sample_id, .N,
                                  by = .(sample_id, gate12g_broad)]
ext_total <- ext_broad[, .(assignment_cells = sum(N)), by = sample_id]
ext_broad[, fraction := N / sum(N), by = sample_id]
ext_broad_wide <- dcast(ext_broad, sample_id ~ gate12g_broad, value.var = "fraction", fill = 0)
ext_baseline <- merge(ext, ext_broad_wide, by = "sample_id", sort = FALSE)
ext_baseline <- merge(ext_baseline, ext_total, by = "sample_id", sort = FALSE)
for (nm in c("Malignant", "T_NK", "Myeloid")) if (!nm %in% names(ext_baseline)) ext_baseline[, (nm) := 0]
ext_baseline[, log10_cells := log10(qc_cells)]
discovery_simple_fit <- lm(Axis1 ~ Malignant + T_NK + Myeloid + log10_cells, data = baseline)
ext_baseline[, predicted_discovery_simple_mixture := as.numeric(predict(discovery_simple_fit, newdata = ext_baseline))]
ext_baseline[, residual_from_discovery_simple_mixture := Axis1 - predicted_discovery_simple_mixture]
fwrite(ext_baseline[, .(sample_id, patient_id, cancer_code, Axis1, qc_cells,
                        broad_assigned_fraction, Malignant, T_NK, Myeloid,
                        predicted_discovery_simple_mixture, residual_from_discovery_simple_mixture)],
       file.path(outdir, "external_axis1_simple_mixture_diagnostic.tsv"), sep = "\t")

external_mixture_summary <- data.table(
  n_patients = nrow(ext_baseline),
  rho_axis1_vs_predicted = cor(ext_baseline$Axis1, ext_baseline$predicted_discovery_simple_mixture,
                              method = "spearman"),
  residual_origin_eta2 = {
    fit <- lm(residual_from_discovery_simple_mixture ~ factor(cancer_code), data = ext_baseline)
    an <- anova(fit)
    an$`Sum Sq`[1L] / sum(an$`Sum Sq`)
  }
)
fwrite(external_mixture_summary, file.path(outdir, "external_axis1_simple_mixture_summary.tsv"), sep = "\t")

# Publication-facing sensitivity figure.
p1 <- ggplot(axis_contrasts[group == "Pooled"], aes(contrast, median_difference)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey45") +
  geom_errorbar(aes(ymin = bootstrap_ci_low, ymax = bootstrap_ci_high), width = 0.12) +
  geom_point(size = 2.4, colour = "#0072B2") +
  coord_flip() + theme_bw(base_size = 10) +
  labs(title = "A  Categorical paired Axis1 contrasts", x = NULL, y = "Median paired difference (95% bootstrap CI)")

p2 <- ggplot(lopo, aes(predicted_simple_mixture, Axis1, colour = compartment)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
  geom_point(size = 1.8, alpha = 0.85) + theme_bw(base_size = 10) +
  labs(title = "B  Patient-held-out simple-mixture model",
       subtitle = sprintf("Cross-validated R2 = %.3f", lopo_summary$cross_validated_R2),
       x = "Held-out prediction", y = "Frozen Axis1")

block_long <- melt(block_contribution, id.vars = c("sample_id", "patient_id", "cancer", "compartment", "Axis1"),
                   measure.vars = unique(model$blocks), variable.name = "block", value.name = "contribution")
p3 <- ggplot(block_long, aes(compartment, contribution, fill = block)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  geom_jitter(aes(colour = block), width = 0.12, size = 0.9, alpha = 0.65, show.legend = FALSE) +
  facet_wrap(~block, scales = "free_y") + theme_bw(base_size = 10) +
  labs(title = "C  Axis1 block decomposition", x = NULL, y = "Exact block contribution")

p4 <- ggplot(tech, aes(metric, rho)) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = 2, colour = "#D55E00") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.12) +
  geom_point(size = 2.5, colour = "#009E73") +
  coord_flip() + theme_bw(base_size = 10) +
  labs(title = "D  External technical-dependence uncertainty", x = NULL, y = "Spearman rho (95% patient bootstrap CI)")

fig <- (p1 | p2) / (p3 | p4) +
  plot_annotation(title = "Gate12R statistical and compositional sensitivity analyses")
ggsave(file.path(outdir, "FigureR2_statistical_sensitivities.pdf"), fig, width = 12, height = 9)
ggsave(file.path(outdir, "FigureR2_statistical_sensitivities.png"), fig, width = 12, height = 9, dpi = 300, bg = "white")

decision <- data.table(
  categorical_tumor_vs_distal_positive = axis_contrasts[group == "Pooled" & contrast == "tumor_minus_distal", median_difference] > 0,
  residual_tumor_vs_distal_positive = residual_contrasts[contrast == "tumor_minus_distal", median_difference] > 0,
  simple_mixture_cv_r2 = lopo_summary$cross_validated_R2,
  simple_mixture_not_sufficient = lopo_summary$cross_validated_R2 < 0.80 &&
    residual_contrasts[contrast == "tumor_minus_distal", median_difference] > 0,
  max_axis_reconstruction_error = max(block_contribution$absolute_reconstruction_error),
  status = "COMPLETE"
)
fwrite(decision, file.path(outdir, "gate12r_statistical_sensitivity_decision.tsv"), sep = "\t")

writeLines(c(
  "# Gate12R statistical sensitivity checkpoint", "",
  paste0("- Categorical pooled tumour-minus-distal Axis1 median: ",
         sprintf("%.4f", axis_contrasts[group == "Pooled" & contrast == "tumor_minus_distal", median_difference])),
  paste0("- Residual tumour-minus-distal Axis1 median after broad-mixture/depth/cancer adjustment: ",
         sprintf("%.4f", residual_contrasts[contrast == "tumor_minus_distal", median_difference])),
  paste0("- Patient-held-out simple-mixture cross-validated R2: ", sprintf("%.4f", lopo_summary$cross_validated_R2)),
  paste0("- Exact decomposition maximum reconstruction error: ",
         format(max(block_contribution$absolute_reconstruction_error), scientific = TRUE)),
  "- REML/Hartung-Knapp estimates use only two cancer/accession cohorts and are explicitly exploratory.",
  "- Technical-dependence intervals are patient bootstrap intervals; the original frozen 0.50 gate is unchanged."
), file.path(outdir, "GATE12R_STATISTICAL_SENSITIVITY_CHECKPOINT.md"))

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
cat("GATE12R_STATISTICAL_SENSITIVITY_STATUS=COMPLETE\n")
print(decision)
