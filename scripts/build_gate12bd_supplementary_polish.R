#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

options(stringsAsFactors = FALSE, warn = 1)
set.seed(20260812L)

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
out <- file.path(root, "results/gate12bd_aligned_submission_package")
fig_dir <- file.path(out, "figures", "supplementary")
src_dir <- file.path(out, "source_data", "supplementary")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)

af <- file.path(root, "results/gate12af_major_redesign")
afs <- file.path(af, "supplementary_redesign", "source_data")
ai <- file.path(root, "results/gate12ai_submission_ready", "source_data")
tech_path <- file.path(root, "results/gate12ab_minor_revision_closure", "analysis", "gse266330_axis1", "gse266330_technical_correlations.tsv")

pal <- list(
  prostate = "#0072B2", renal = "#D55E00", broad = "#009E73", myeloid = "#E69F00",
  tnk = "#0072B2", purple = "#6F4C9B", dark = "#222222", blue = "#2166AC", red = "#B2182B"
)
cancer_cols <- c(Prostate = pal$prostate, Renal = pal$renal)

theme_pub <- theme_classic(base_size = 9, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 9.2, face = "bold", hjust = 0, margin = margin(b = 2.0)),
    plot.subtitle = element_text(size = 7.4, colour = "#4A4A4A", lineheight = .96, margin = margin(b = 4.0)),
    axis.title = element_text(size = 8.2), axis.text = element_text(size = 7.6, colour = "#303030"),
    strip.text = element_text(size = 7.7, face = "bold"),
    strip.background = element_rect(fill = "#F2F2F2", colour = "#D0D0D0", linewidth = .25),
    legend.title = element_text(size = 7.5), legend.text = element_text(size = 7.4),
    plot.tag = element_text(size = 11, face = "bold"), plot.margin = margin(7, 8, 6, 6)
  )
theme_heat <- theme_minimal(base_size = 9, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 9.2, face = "bold", hjust = 0, margin = margin(b = 2.0)),
    plot.subtitle = element_text(size = 7.4, colour = "#4A4A4A", lineheight = .96, margin = margin(b = 4.0)),
    axis.title = element_text(size = 8.2), axis.text = element_text(size = 7.6, colour = "#303030"),
    strip.text = element_text(size = 7.7, face = "bold"), panel.grid = element_blank(),
    legend.title = element_text(size = 7.5), legend.text = element_text(size = 7.4),
    plot.tag = element_text(size = 11, face = "bold"), plot.margin = margin(7, 8, 6, 6)
  )
tag_theme <- theme(plot.tag = element_text(family = "Arial", size = 11, face = "bold"))

cairo_png_device <- function(filename, width, height, bg = "white", ...) {
  grDevices::png(filename, width = width, height = height, units = "in", res = 450,
                 type = "cairo", bg = bg, ...)
}
save_figure <- function(plot, stem, width, height) {
  ggsave(file.path(fig_dir, paste0(stem, ".png")), plot, width = width, height = height,
         dpi = 450, device = cairo_png_device, bg = "white", limitsize = FALSE)
  ggsave(file.path(fig_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         device = cairo_pdf, bg = "white", limitsize = FALSE)
}
write_source <- function(x, stem) fwrite(x, file.path(src_dir, stem), sep = "\t", quote = FALSE, na = "NA",
                                         compress = if (endsWith(stem, ".gz")) "gzip" else "none")

## Figure S3: conventional bar, dot-heatmap, and coefficient forest grammar.
persist3 <- fread(file.path(afs, "FigureS3A_persistence_n_over_N.tsv"))
hall3 <- fread(file.path(afs, "FigureS3B_robust_hallmarks.tsv"))
cat3 <- fread(file.path(af, "source_data", "Figure3E_categorical_anatomical_contrasts.tsv"))
persist3[, state := factor(state, levels = c("CD4 T", "CD8/CTL", "NK/NKT"))]
pS3a <- ggplot(persist3, aes(state, fraction, fill = domain)) +
  geom_col(position = position_dodge(width = .68), width = .58, colour = "#333333", linewidth = .2) +
  geom_text(aes(label = label), position = position_dodge(width = .68), vjust = -.34, size = 2.45) +
  scale_fill_manual(values = c(Genes = "#A7C7E7", Hallmarks = "#F2C894")) +
  scale_y_continuous(labels = percent, limits = c(0, 1.13), expand = expansion(mult = c(0, .02))) +
  labs(title = "Sensitivity persistence", subtitle = "Labels show persistent / eligible features",
       x = NULL, y = "Persistent / eligible", fill = NULL) + theme_pub +
  theme(legend.position = "bottom", legend.box.margin = margin(t = -3))
hall3[, state := factor(state, levels = c("CD4 T", "CD8/CTL", "NK/NKT"))]
hall3[, pathway_label := factor(pathway_label, levels = rev(unique(hall3[order(mean_nes)]$pathway_label)))]
pS3b <- ggplot(hall3, aes(state, pathway_label, colour = mean_nes, size = abs(mean_nes))) +
  geom_point(alpha = .92) +
  scale_colour_gradient(low = "#D9D9D9", high = pal$purple, name = "Mean NES") +
  scale_size_continuous(range = c(1.2, 3.1), name = "|Mean NES|") +
  labs(title = "Robust Hallmark effects", subtitle = "Internal functional annotation",
       x = "T/NK state", y = NULL) + theme_pub +
  theme(axis.text.y = element_text(size = 7.1), axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = "right")
cat3[, contrast_label := factor(contrast_label, levels = unique(contrast_label))]
pS3c <- ggplot(cat3, aes(logFC, gene_state, colour = contrast_label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#777777", linewidth = .35) +
  geom_errorbar(aes(xmin = low, xmax = high), orientation = "y",
                position = position_dodge(width = .5), width = .12, linewidth = .35) +
  geom_point(position = position_dodge(width = .5), size = 1.25) +
  facet_wrap(~cancer_label, nrow = 1) +
  labs(title = "Categorical anatomical contrasts", subtitle = "Anatomical comparisons; no temporal interpretation",
       x = "Patient-paired pseudobulk contrast", y = NULL, colour = NULL) + theme_pub +
  theme(legend.position = "bottom")
figS3 <- (pS3a | pS3b) / pS3c +
  plot_layout(heights = c(1.10, 1.0), widths = c(.78, 1.22)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)
save_figure(figS3, "FigureS3_tnk_transcription_robustness", 7.2, 8.2)
write_source(persist3, "FigureS3A_persistence_n_over_N.tsv")
write_source(hall3, "FigureS3B_robust_hallmarks.tsv")
write_source(cat3, "FigureS3C_categorical_anatomical_contrasts.tsv")

## Figure S5: no table-like dashboard, no connecting lines, no overlapping headers.
contr5 <- fread(file.path(af, "source_data", "Figure5B_exact_block_contributions.tsv"))
context5 <- fread(file.path(afs, "FigureS5B_context_attribution.tsv"))
oep5 <- fread(file.path(afs, "FigureS5C_oep_scores.tsv"))
contr5[, compartment := factor(compartment, levels = c("Distal", "Involved", "Tumor"), labels = c("D", "I", "T"))]
sum5 <- contr5[, .(median = median(contribution), q1 = quantile(contribution, .25), q3 = quantile(contribution, .75)),
               by = .(cancer, compartment, block_display)]
pS5a <- ggplot(contr5, aes(compartment, contribution, colour = cancer, shape = cancer)) +
  geom_hline(yintercept = 0, colour = "#D0D0D0") +
  geom_point(position = position_jitterdodge(jitter.width = .11, dodge.width = .42), alpha = .28, size = .82) +
  geom_errorbar(data = sum5, aes(y = median, ymin = q1, ymax = q3),
                position = position_dodge(width = .42), width = .11, linewidth = .48) +
  geom_point(data = sum5, aes(y = median), position = position_dodge(width = .42), size = 1.65) +
  facet_wrap(~block_display, nrow = 1) +
  scale_colour_manual(values = cancer_cols) + scale_shape_manual(values = c(Prostate = 16, Renal = 17)) +
  labs(title = "Exact block contributions", subtitle = "Patient-samples with median and IQR; D/I/T are anatomical categories",
       x = NULL, y = "Signed contribution to Axis1", colour = NULL, shape = NULL) + theme_pub +
  theme(legend.position = "bottom")
context5[, attribution := factor(attribution, levels = c("Marginal R2", "Drop-one delta R2"))]
context5[, term_label := factor(term_label, levels = rev(c("Anatomical compartment", "Patient", "Cancer/accession", "Cell depth")))]
pS5b <- ggplot(context5, aes(R2, term_label, colour = attribution, shape = attribution)) +
  geom_point(position = position_dodge(width = .45), size = 2.15) +
  scale_colour_manual(values = c(`Marginal R2` = pal$tnk, `Drop-one delta R2` = pal$myeloid), name = NULL) +
  scale_shape_manual(values = c(`Marginal R2` = 16, `Drop-one delta R2` = 17), name = NULL) +
  labs(title = "Context attribution", subtitle = "Descriptive; cancer/accession is aliased with patient",
       x = "R2 or drop-one delta R2", y = NULL) + theme_pub +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank(), legend.box.margin = margin(t = -2))
stats5 <- oep5[, .(n = .N, median = median(Axis1)), by = cancer_code][order(median)]
oep5[, cancer_code := factor(cancer_code, levels = stats5$cancer_code)]
stats5[, cancer_code := factor(cancer_code, levels = levels(oep5$cancer_code))]
pS5c <- ggplot(oep5, aes(cancer_code, Axis1)) +
  geom_hline(yintercept = 0, linetype = 3, colour = "#777777") +
  geom_jitter(width = .10, height = 0, colour = pal$purple, size = 1.15, alpha = .78) +
  stat_summary(data = oep5[cancer_code %chin% stats5[n >= 3, cancer_code]], fun = median,
               geom = "crossbar", width = .46, linewidth = .4) +
  geom_text(data = stats5, aes(y = max(oep5$Axis1) + .12, label = paste0("n=", n)), size = 2.4) +
  labs(title = "OEP005136 technical projection", subtitle = "Source categories; medians shown for n >= 3",
       x = "Source-defined category", y = "Frozen Axis1") + theme_pub +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(.04, .16)))
figS5 <- pS5a / (pS5b | pS5c) +
  plot_layout(heights = c(.92, 1.08), widths = c(1.02, .98)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)
save_figure(figS5, "FigureS5_axis1_context_diagnostics", 7.2, 7.8)
write_source(contr5, "FigureS5A_exact_block_contributions.tsv")
write_source(context5, "FigureS5B_context_attribution.tsv")
write_source(oep5, "FigureS5C_oep_scores.tsv")

## Figure S8: remove edge labels from the deterministic deletion panel.
gse8 <- fread(file.path(afs, "FigureS8A_all_projectable_scores.tsv"))
loco8 <- fread(file.path(afs, "FigureS8B_leave_one_origin_out.tsv"))
tech8 <- fread(file.path(afs, "FigureS8C_technical_diagnostics.tsv"))
techsum8 <- fread(tech_path)
full_hl8 <- median(outer(gse8[condition == "bone_metastasis", Axis1], gse8[condition != "bone_metastasis", Axis1], "-"))
stats8 <- gse8[, .(n = .N, median = median(Axis1)), by = origin][order(median)]
gse8[, origin := factor(origin, levels = stats8$origin)]
stats8[, origin := factor(origin, levels = levels(gse8$origin))]
pS8a <- ggplot(gse8, aes(origin, Axis1)) +
  geom_jitter(width = .10, height = 0, colour = pal$purple, size = 1.12, alpha = .78) +
  stat_summary(data = gse8[origin %chin% stats8[n >= 3, origin]], fun = median,
               geom = "crossbar", width = .45, linewidth = .38) +
  geom_text(data = stats8, aes(y = max(gse8$Axis1) + .12, label = paste0("n=", n)), size = 2.4) +
  labs(title = "Every projectable donor or patient", subtitle = "Descriptive origins; medians shown for n >= 3",
       x = NULL, y = "Frozen Axis1") + theme_pub +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) + scale_y_continuous(expand = expansion(mult = c(.04, .16)))
loco8[, delta_from_full := hodges_lehmann_shift - full_hl8]
loco8[, dropped_origin := factor(dropped_origin, levels = rev(dropped_origin))]
pS8b <- ggplot(loco8, aes(delta_from_full, dropped_origin)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#777777") +
  geom_segment(aes(x = 0, xend = delta_from_full, yend = dropped_origin), colour = "#B9B9B9", linewidth = .45) +
  geom_point(colour = pal$broad, size = 2.05) +
  labs(title = "Leave-one-origin-out sensitivity", subtitle = sprintf("Deviation from full Hodges-Lehmann shift %.3f", full_hl8),
       x = "Deleted-origin estimate - full estimate", y = "Deleted origin") + theme_pub
metric_map <- c(`log10 QC cells` = "rho_log10_qc_cells", `Broad assignment fraction` = "rho_broad_assigned_fraction")
ann8 <- data.table(metric = names(metric_map), stat_key = unname(metric_map))
ann8 <- merge(ann8, techsum8, by.x = "stat_key", by.y = "metric", all.x = TRUE)
ann8[, label := sprintf("Spearman rho %.2f\n95%% bootstrap %.2f to %.2f\nfrozen |rho| < 0.50", estimate, ci_low, ci_high)]
loc8 <- tech8[, .(x = min(value), y = max(Axis1)), by = metric]
ann8 <- merge(ann8, loc8, by = "metric", all.x = TRUE)
pS8c <- ggplot(tech8, aes(value, Axis1)) +
  geom_point(size = 1.12, alpha = .75, colour = pal$purple) +
  geom_text(data = ann8, aes(x = x, y = y, label = label), hjust = 0, vjust = 1, size = 2.35, inherit.aes = FALSE) +
  facet_wrap(~metric, scales = "free_x", nrow = 1) +
  labs(title = "Technical rank associations", subtitle = "Donor/patient points; no fitted line",
       x = NULL, y = "Frozen Axis1") + theme_pub
figS8 <- pS8a / (pS8b | pS8c) + plot_layout(heights = c(.80, 1.0), widths = c(.82, 1.18)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)
save_figure(figS8, "FigureS8_external_heterogeneity_sensitivity", 7.2, 8.1)
write_source(gse8, "FigureS8A_all_projectable_scores.tsv")
write_source(loco8, "FigureS8B_leave_one_origin_out.tsv")
write_source(tech8, "FigureS8C_technical_diagnostics.tsv")
write_source(techsum8, "FigureS8C_technical_correlations_with_bootstrap.tsv")

## Figure S10: a standard sensitivity heatmap plus the section-stratified effect matrix.
summary10 <- fread(file.path(afs, "FigureS9D_grid_summary.tsv"))
neigh10 <- fread(file.path(ai, "Figure6E_nine_class_neighborhood_matrix.tsv"))
summary10[, label := sprintf("%d/%d", supporting_sections, decision_evaluable_sections)]
pS10a <- ggplot(summary10, aes(factor(block_width), factor(k), fill = supporting_sections / decision_evaluable_sections)) +
  geom_tile(colour = "white", linewidth = .6) +
  geom_tile(data = summary10[k == 6 & block_width == 20], fill = NA, colour = "black", linewidth = .9) +
  geom_text(aes(label = label), size = 2.8, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "#BFDAD1", limits = c(0, 1), name = "Supporting /\nevaluable") +
  labs(title = "Linked-section sensitivity grid", subtitle = "Outline marks k=6 and width=20",
       x = "Spatial-block width", y = "Symmetric graph k") + theme_heat
samples10 <- unique(neigh10$sample)
neigh10[, sample := factor(sample, levels = samples10, labels = paste0("S", seq_along(samples10)))]
neigh10[, class_label := factor(class_label, levels = rev(sort(unique(class_label))))]
pS10b <- ggplot(neigh10, aes(sample, class_label)) +
  geom_tile(aes(fill = ifelse(decision_evaluable, estimate, NA_real_)), colour = "white", linewidth = .5) +
  geom_point(data = neigh10[decision_evaluable == TRUE & interval_excludes_zero == TRUE],
             shape = 21, fill = NA, colour = "#111111", size = 2.3, stroke = .5) +
  geom_text(data = neigh10[decision_evaluable == FALSE], label = "NE", colour = "#777777", size = 2.45) +
  scale_fill_gradient2(low = pal$blue, mid = "white", high = pal$red, midpoint = 0,
                       na.value = "#EEEEEE", name = "Boundary - distal") +
  labs(title = "Nine-class section effects", subtitle = "S1-S4 shown separately; no pooled inference",
       x = "Linked section", y = NULL) + theme_heat +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "right")
figS10 <- (pS10a | pS10b) + plot_layout(widths = c(.68, 1.32)) +
  plot_annotation(tag_levels = "a", theme = tag_theme)
save_figure(figS10, "FigureS10_spatial_sensitivity", 7.2, 4.8)
write_source(summary10, "FigureS10A_grid_support_summary.tsv")
write_source(neigh10, "FigureS10B_nine_class_neighborhood_matrix.tsv")

## Complete panel-level source-data archive for the retained supplementary pages.
source_copies <- setNames(
  c(
    "FigureS4A_consensus_diagnostics.tsv", "FigureS4B_ari_null.tsv.gz", "FigureS4B_ari_receipt.tsv",
    "TableS4_representation_evidence_boundary.tsv", "FigureS6A_rule_retention.tsv",
    "FigureS6B_resampling_intervals.tsv", "FigureS6B_resampling_values.tsv",
    "TableS6_claim_evidence_boundary.tsv", "FigureS7A_all_hallmark_effects.tsv",
    "FigureS7B_support_matrix.tsv", "FigureS7C_recurrent_member_genes.tsv",
    "FigureS9A_C_spatial_scores.tsv.gz"
  ),
  c(
    file.path(afs, "FigureS4A_consensus_diagnostics.tsv"),
    file.path(afs, "FigureS4B_ari_null.tsv.gz"),
    file.path(afs, "FigureS4B_ari_receipt.tsv"),
    file.path(afs, "FigureS4C_evidence_table.tsv"),
    file.path(afs, "FigureS6A_rule_retention.tsv"),
    file.path(afs, "FigureS6B_resampling_intervals.tsv"),
    file.path(afs, "FigureS6B_resampling_values.tsv"),
    file.path(afs, "FigureS6C_evidence_table.tsv"),
    file.path(afs, "FigureS7A_all_hallmark_effects.tsv"),
    file.path(afs, "FigureS7B_support_matrix.tsv"),
    file.path(afs, "FigureS7C_recurrent_member_genes.tsv"),
    file.path(afs, "FigureS9A_B_spatial_scores.tsv.gz")
  )
)
for (from in names(source_copies)) {
  if (!file.exists(from)) stop("Missing supplementary source input: ", from)
  to <- file.path(src_dir, unname(source_copies[[from]]))
  if (!file.copy(from, to, overwrite = TRUE)) stop("Could not copy supplementary source data: ", from)
}

message("Gate12BD supplementary polish completed: S3, S5, S8 and S10.")
