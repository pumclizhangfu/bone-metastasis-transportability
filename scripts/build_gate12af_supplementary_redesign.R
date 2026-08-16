suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(ggplot2)
  library(patchwork)
  library(png)
  library(scales)
})

options(stringsAsFactors = FALSE, warn = 1)
set.seed(20260810L)

root <- normalizePath(".", mustWork = TRUE)
out_dir <- file.path(root, "results/gate12af_major_redesign/supplementary_redesign")
fig_dir <- file.path(out_dir, "figures/supplementary")
source_dir <- file.path(out_dir, "source_data")
admin_dir <- file.path(out_dir, "admin")
prov_dir <- file.path(out_dir, "provenance")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(admin_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(prov_dir, recursive = TRUE, showWarnings = FALSE)

parent <- "results/gate12ab_minor_revision_closure/source_data"
phase_a <- "results/gate12ad_figure_restructure/phase_a_source_provenance"
old_supp <- "results/gate12ad_figure_restructure/phase_c_supplement_manuscript/source_data"

used_inputs <- character()
track <- function(rel) {
  if (!file.exists(rel)) stop("Missing input: ", rel, call. = FALSE)
  used_inputs <<- unique(c(used_inputs, rel))
  rel
}
read_dt <- function(rel) fread(track(rel))
sha256_file <- function(x) digest(x, file = TRUE, algo = "sha256")
rel_path <- function(x) substring(normalizePath(x, mustWork = TRUE), nchar(root) + 2L)
assert <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
write_source <- function(x, filename) {
  target <- file.path(source_dir, filename)
  fwrite(x, target, sep = "\t", quote = TRUE, na = "NA", compress = if (endsWith(filename, ".gz")) "gzip" else "none")
  invisible(target)
}
cairo_png_device <- function(filename, width, height, bg = "white", ...) {
  grDevices::png(filename, width = width, height = height, units = "in", res = 400, type = "cairo", bg = bg, ...)
}

pal <- list(blue = "#0077BB", cyan = "#33BBEE", teal = "#009988", orange = "#EE7733",
            red = "#CC3311", purple = "#6F4C9B", grey = "#A9A9A9", dark = "#252525")
theme_pub <- theme_classic(base_size = 8.2, base_family = "Arial") +
  theme(plot.title = element_text(size = 9.2, face = "bold", hjust = 0, margin = margin(b = 1.5)),
        plot.subtitle = element_text(size = 7.0, colour = "#4A4A4A", margin = margin(b = 2.2)),
        axis.title = element_text(size = 7.8), axis.text = element_text(size = 7.0, colour = "#303030"),
        strip.text = element_text(size = 7.2, face = "bold"),
        strip.background = element_rect(fill = "#F1F1F1", colour = "#D0D0D0", linewidth = 0.25),
        legend.title = element_text(size = 7.0), legend.text = element_text(size = 6.8),
        panel.grid.minor = element_blank(), plot.tag = element_text(size = 11, face = "bold"),
        plot.margin = margin(3, 3, 3, 3))
theme_heat <- theme_minimal(base_size = 8.0, base_family = "Arial") +
  theme(plot.title = element_text(size = 9.2, face = "bold", hjust = 0, margin = margin(b = 1.5)),
        plot.subtitle = element_text(size = 7.0, colour = "#4A4A4A", margin = margin(b = 2.2)),
        axis.title = element_text(size = 7.8), axis.text = element_text(size = 7.0, colour = "#303030"),
        strip.text = element_text(size = 7.2, face = "bold"), panel.grid = element_blank(),
        legend.title = element_text(size = 7.0), legend.text = element_text(size = 6.8),
        plot.tag = element_text(size = 11, face = "bold"), plot.margin = margin(3, 3, 3, 3))

receipts <- list()
save_figure <- function(plot, filename, width, height) {
  png_path <- file.path(fig_dir, paste0(filename, ".png"))
  pdf_path <- file.path(fig_dir, paste0(filename, ".pdf"))
  ggsave(png_path, plot, width = width, height = height, dpi = 400, device = cairo_png_device, bg = "white", limitsize = FALSE)
  ggsave(pdf_path, plot, width = width, height = height, device = cairo_pdf, bg = "white", limitsize = FALSE)
  data.table(figure = filename, width_in = width, height_in = height, dpi_png = 400,
             png_relative_path = rel_path(png_path), png_bytes = file.info(png_path)$size, png_sha256 = sha256_file(png_path),
             pdf_relative_path = rel_path(pdf_path), pdf_bytes = file.info(pdf_path)$size, pdf_sha256 = sha256_file(pdf_path))
}

## S1: cohort, QC and annotation ------------------------------------------------
roles <- read_dt(file.path(parent, "Figure1_dataset_role_matrix.tsv"))
avail <- read_dt(file.path(parent, "Figure1_patient_sample_availability.tsv"))
markers <- read_dt(file.path(phase_a, "source_data/Figure1C_canonical_marker_dotplot.tsv"))
composition <- read_dt(file.path(phase_a, "source_data/Figure1D_all_sample_broad_composition.tsv"))

unit_map <- c(GSE143791 = "patient", GSE202813 = "patient", GSE266330 = "donor/patient",
              OEP005136 = "patient", GSE323357 = "linked section/block")
roles[, dataset_label := paste0(dataset, "  [", unit_map[dataset], "]")]
roles[, role := factor(role, levels = c("Discovery", "Fixed program", "Ecotype", "Axis transfer", "Paired biology", "Spatial"))]
roles[, dataset_label := factor(dataset_label, levels = rev(unique(dataset_label)))]
pS1a <- ggplot(roles, aes(role, dataset_label)) +
  geom_tile(fill = "white", colour = "#D9D9D9", linewidth = 0.45) +
  geom_text(aes(label = ifelse(used, "YES", "-"), colour = used), size = 2.4, fontface = "bold") +
  scale_colour_manual(values = c(`TRUE` = pal$blue, `FALSE` = "#B8B8B8"), guide = "none") +
  labs(title = "Dataset roles and biological units", subtitle = "Frozen roles; brackets show the inferential unit",
       x = NULL, y = NULL) + theme_heat + theme(axis.text.x = element_text(angle = 25, hjust = 1))

avail[, cancer := factor(cancer, levels = c("prostate", "renal"), labels = c("Prostate", "Renal"))]
avail[, compartment := factor(compartment, levels = c("Distal", "Involved", "Tumor"))]
avail[, patient_id := factor(patient_id, levels = rev(unique(patient_id)))]
pS1b <- ggplot(avail, aes(compartment, patient_id, fill = present)) +
  geom_tile(colour = "white", linewidth = 0.55) + facet_grid(cancer ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c(`TRUE` = pal$teal, `FALSE` = "white"), guide = "none") +
  labs(title = "Discovery specimen availability", subtitle = "18 patients | 42 samples | 11 complete triplets",
       x = NULL, y = NULL) + theme_heat

markers[, broad_class_label := fifelse(broad_class_label == "Unassigned", "Unassigned (QC)", broad_class_label)]
markers[, broad_class_label := factor(broad_class_label, levels = rev(unique(markers[order(broad_class_order)]$broad_class_label)))]
markers[, gene := factor(gene, levels = unique(markers[order(marker_gene_order)]$gene))]
pS1c <- ggplot(markers, aes(gene, broad_class_label, size = detected_percent, colour = mean_scaled_expression)) +
  geom_point(alpha = 0.90) +
  scale_size_continuous(range = c(0.25, 4.2), breaks = c(10, 40, 70), name = "Detected (%)") +
  scale_colour_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0, name = "Gene-wise\nscaled mean") +
  labs(title = "Canonical marker annotation audit", subtitle = "Unassigned is a QC category, not a biological lineage",
       x = NULL, y = NULL) + theme_heat + theme(axis.text.x = element_text(angle = 38, hjust = 1), legend.position = "bottom")

composition[, cancer_label := factor(cancer, levels = c("prostate", "renal"), labels = c("Prostate", "Renal"))]
composition[, compartment_label := factor(compartment, levels = c("distal", "involved", "tumor"), labels = c("D", "I", "T"))]
composition[, sample_id := factor(sample_id, levels = unique(composition[order(sample_panel_order)]$sample_id))]
composition[, broad_class_label := factor(broad_class_label, levels = unique(composition[order(broad_class_order)]$broad_class_label))]
broad_cols <- c("T/NK" = pal$blue, Myeloid = pal$orange, "B cell" = "#33BBEE", Progenitor = "#D8B365",
                Stromal = pal$teal, Endothelial = "#44AA99", Osteoclast = "#AA4499", Osteoblast = "#CC6677",
                Malignant = pal$purple, Unassigned = "#BBBBBB")
pS1d <- ggplot(composition, aes(sample_id, fraction, fill = broad_class_label)) +
  geom_col(width = 0.93, colour = "white", linewidth = 0.08) +
  facet_grid(~cancer_label + compartment_label, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = broad_cols, name = NULL) + scale_y_continuous(labels = percent, expand = c(0, 0)) +
  labs(title = "Complete 42-sample broad composition", subtitle = "Descriptive only; no cell-level inference is drawn from stacked bars",
       x = "Sample", y = "Cell fraction") + theme_pub +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4.7), legend.position = "bottom")

figS1 <- (pS1a | pS1b) / pS1c / pS1d + plot_layout(heights = c(0.82, 0.90, 1.05), widths = c(1.15, 0.85)) + plot_annotation(tag_levels = "A")
receipts[[1]] <- save_figure(figS1, "FigureS1_cohort_qc_annotation", 7.2, 10.0)
write_source(roles, "FigureS1A_dataset_role_table.tsv"); write_source(avail, "FigureS1B_availability.tsv")
write_source(markers, "FigureS1C_marker_audit.tsv"); write_source(composition, "FigureS1D_full_composition.tsv")

## S2: compositional robustness ------------------------------------------------
clr <- read_dt(file.path(parent, "Figure2_clr_sensitivity_concordance.tsv"))
clr[, state_label := gsub("_", " ", state)]
clr[, state_label := gsub("T NK", "T/NK", state_label, fixed = TRUE)]
clr[, state_label := gsub("CD8 CTL", "CD8/CTL", state_label, fixed = TRUE)]
clr[, lineage := factor(lineage, levels = c("Myeloid", "T_NK"))]
rho_clr <- cor(clr$original_beta, clr$clr_beta, method = "spearman")
lim2 <- max(abs(c(clr$original_beta, clr$clr_beta))) * 1.08
pS2a <- ggplot(clr, aes(original_beta, clr_beta, colour = lineage, shape = cross_cancer_stable)) +
  geom_hline(yintercept = 0, colour = "#D0D0D0") + geom_vline(xintercept = 0, colour = "#D0D0D0") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#777777") + geom_point(size = 2.0) +
  coord_equal(xlim = c(-lim2, lim2), ylim = c(-lim2, lim2)) +
  scale_colour_manual(values = c(Myeloid = pal$orange, T_NK = pal$blue), labels = c(T_NK = "T/NK")) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.25, label = sprintf("Spearman rho %.2f", rho_clr), size = 2.1) +
  labs(title = "Empirical-logit versus paired CLR", subtitle = "Same scale and 1:1 line; filled symbols were originally stable",
       x = "Empirical-logit effect", y = "Paired CLR effect", colour = "Lineage", shape = "Originally stable") +
  theme_pub + theme(legend.position = "none")
clr[, state_factor := factor(state_label, levels = rev(clr[order(lineage, clr_beta)]$state_label))]
pS2b <- ggplot(clr, aes(clr_beta, state_factor, colour = lineage)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#666666") +
  geom_errorbar(aes(xmin = clr_low, xmax = clr_high), orientation = "y", width = 0.12, linewidth = 0.5) +
  geom_point(aes(fill = cross_cancer_stable), shape = 21, size = 1.9, stroke = 0.45) +
  scale_colour_manual(values = c(Myeloid = pal$orange, T_NK = pal$blue), guide = "none") +
  scale_fill_manual(values = c(`TRUE` = pal$dark, `FALSE` = "white"), guide = "none") +
  labs(title = "Paired CLR sensitivity for all 16 states", subtitle = "16/16 directions retained across 6/6 prespecified sensitivities",
       x = "CLR change per anatomical step (95% patient-bootstrap CI)", y = NULL) + theme_pub + theme(panel.grid.major.y = element_blank())
figS2 <- (pS2a | pS2b) + plot_annotation(tag_levels = "A")
receipts[[2]] <- save_figure(figS2, "FigureS2_composition_robustness", 7.2, 4.8)
write_source(clr, "FigureS2_composition_robustness.tsv")

## S3: T/NK transcription robustness ------------------------------------------
persist <- read_dt(file.path(parent, "Figure3_robustness_persistence.tsv"))
hall3 <- read_dt(file.path(parent, "Figure3_version_frozen_hallmarks.tsv"))
persist[, label := sprintf("%d/%d", persistent, candidate)]
persist[, label_x := pmin(fraction + 0.075, 1.10)]
persist[, state := factor(state, levels = c("CD4 T", "CD8/CTL", "NK/NKT"))]
persist[, domain := factor(domain, levels = c("Genes", "Hallmarks"))]
pS3a <- ggplot(persist, aes(fraction, state, colour = domain, shape = domain)) +
  geom_segment(aes(x = 0, xend = fraction, yend = state), position = position_dodge(width = 0.42), colour = "#D0D0D0", linewidth = 0.55) +
  geom_point(position = position_dodge(width = 0.42), size = 2.4) +
  geom_text(aes(x = label_x, label = label), position = position_dodge(width = 0.42), hjust = 0,
            size = 1.9, show.legend = FALSE) +
  scale_colour_manual(values = c(Genes = pal$blue, Hallmarks = pal$teal)) +
  scale_shape_manual(values = c(Genes = 16, Hallmarks = 17)) + scale_x_continuous(labels = percent, limits = c(0, 1.18)) +
  labs(title = "Sensitivity-persistent T/NK results", subtitle = "Labels are persistent/eligible results, not cell counts",
       x = "Persistent fraction", y = NULL, colour = NULL, shape = NULL) + theme_pub + theme(legend.position = "bottom")
hall3[, state := factor(state, levels = c("CD4 T", "CD8/CTL", "NK/NKT"))]
hall3[, pathway_label := factor(pathway_label, levels = rev(unique(hall3[order(mean_nes)]$pathway_label)))]
pS3b <- ggplot(hall3, aes(mean_nes, pathway_label, colour = state, shape = state)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#777777") + geom_point(position = position_dodge(width = 0.48), size = 1.8) +
  scale_colour_manual(values = c(`CD4 T` = pal$blue, `CD8/CTL` = pal$orange, `NK/NKT` = pal$teal)) +
  scale_shape_manual(values = c(`CD4 T` = 16, `CD8/CTL` = 17, `NK/NKT` = 15)) +
  labs(title = "Version-frozen robust Hallmark effects", subtitle = "Cohort-direction and patient-deletion support required",
       x = "Mean prostate/renal normalized enrichment score", y = NULL, colour = "State", shape = "State") +
  theme_pub + theme(legend.position = "bottom", panel.grid.major.y = element_blank())
figS3 <- (pS3a | pS3b) + plot_layout(widths = c(0.78, 1.22)) + plot_annotation(tag_levels = "A")
receipts[[3]] <- save_figure(figS3, "FigureS3_tnk_transcription_robustness", 7.2, 5.6)
write_source(persist, "FigureS3A_persistence_n_over_N.tsv"); write_source(hall3, "FigureS3B_robust_hallmarks.tsv")

## S4: fixed-program/ecotype diagnostics --------------------------------------
diag4 <- read_dt(file.path(parent, "Figure4_consensus_diagnostics.tsv"))
null4 <- read_dt(file.path(phase_a, "source_data/Figure4C_ari_permutation_null.tsv.gz"))
receipt4 <- read_dt(file.path(phase_a, "provenance/Figure4C_ari_permutation_receipt.tsv"))
rep4 <- read_dt(file.path(parent, "Figure4_representation_boundary.tsv"))
thresholds <- data.table(metric = c("pac", "mean_silhouette", "min_cluster_size"), threshold = c(0.20, 0.25, 5))
diag4 <- merge(diag4, thresholds, by = "metric", all.x = TRUE)
diag4[, metric_label := factor(metric, levels = c("mean_silhouette", "min_cluster_size", "pac"),
                               labels = c("Silhouette (higher better)", "Minimum size", "PAC (lower better)"))]
pS4a <- ggplot(diag4, aes(k, value)) + geom_line(colour = "#888888", linewidth = 0.45) +
  geom_point(aes(fill = k == 3), shape = 21, size = 2.1) +
  geom_hline(aes(yintercept = threshold), linetype = 2, colour = pal$red, linewidth = 0.4) +
  facet_wrap(~metric_label, scales = "free_y", nrow = 1) + scale_fill_manual(values = c(`TRUE` = pal$dark, `FALSE` = "white"), guide = "none") +
  scale_x_continuous(breaks = 2:5) + labs(title = "Consensus diagnostics", subtitle = "Black: frozen k=3 | red: prespecified threshold",
       x = "Candidate k", y = NULL) + theme_pub
obs4 <- receipt4$observed_ari[[1]]
pS4b <- ggplot(null4, aes(null_ari)) + geom_histogram(bins = 44, fill = "#B7C9D6", colour = "white") +
  geom_vline(xintercept = obs4, colour = pal$purple, linewidth = 0.75) +
  annotate("text", x = obs4, y = Inf, vjust = 1.2, hjust = -0.05,
           label = sprintf("Observed ARI %.3f\nExact P %.3f", obs4, receipt4$empirical_p[[1]]), colour = pal$purple, size = 2.0) +
  labs(title = "Exact 10,000-permutation replay", x = "Permuted adjusted Rand index", y = "Permutations") + theme_pub
rep4[, representation := factor(representation, levels = rev(unique(representation)))]
rep4[, evidence_layer := factor(evidence_layer, levels = c("Discovery", "External direction", "External structure/agreement"))]
rep_cols <- c(`Stable intersection` = "#DCEFEA", `Internally specified` = "#DCEAF5", `Not supported` = "#E5E5E5", `Not tested` = "white")
pS4c <- ggplot(rep4, aes(evidence_layer, representation, fill = status)) +
  geom_tile(colour = "#BFBFBF", linewidth = 0.45) + geom_text(aes(label = status), size = 1.8, lineheight = 0.9) +
  scale_fill_manual(values = rep_cols, guide = "none") +
  labs(title = "Evidence table", subtitle = "Discovery stability cannot replace an external criterion", x = NULL, y = NULL) +
  theme_heat + theme(axis.text.x = element_text(angle = 20, hjust = 1))
figS4 <- pS4a / (pS4b | pS4c) + plot_layout(heights = c(0.78, 1.0), widths = c(0.88, 1.12)) + plot_annotation(tag_levels = "A")
receipts[[4]] <- save_figure(figS4, "FigureS4_fixed_representation_diagnostics", 7.2, 6.2)
write_source(diag4, "FigureS4A_consensus_diagnostics.tsv"); write_source(null4, "FigureS4B_ari_null.tsv.gz")
write_source(receipt4, "FigureS4B_ari_receipt.tsv"); write_source(rep4, "FigureS4C_evidence_table.tsv")

## S5: Axis1 construction/context diagnostics ---------------------------------
contrib5 <- read_dt(file.path(parent, "Figure5_block_contributions.tsv"))
context5 <- read_dt(file.path(parent, "Figure5_context_variance_decomposition.tsv"))
oep5 <- read_dt(file.path(parent, "FigureS1_OEP_projectable_scores.tsv"))
contrib5[, compartment := factor(tolower(compartment), levels = c("distal", "involved", "tumor"), labels = c("D", "I", "T"))]
contrib5[, block_display := factor(block_display, levels = c("Broad", "Myeloid", "T/NK"))]
contrib_sum <- contrib5[, .(median = median(contribution), q1 = quantile(contribution, .25), q3 = quantile(contribution, .75)),
                        by = .(cancer, compartment, block_display)]
pS5a <- ggplot(contrib5, aes(compartment, contribution, colour = cancer, shape = cancer)) +
  geom_hline(yintercept = 0, colour = "#D0D0D0") +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.42), alpha = 0.30, size = 0.9) +
  geom_errorbar(data = contrib_sum, aes(y = median, ymin = q1, ymax = q3), position = position_dodge(width = 0.42), width = 0.12, linewidth = 0.52) +
  geom_point(data = contrib_sum, aes(y = median), position = position_dodge(width = 0.42), size = 1.9) +
  facet_wrap(~block_display, nrow = 1) + scale_colour_manual(values = c(prostate = pal$blue, renal = pal$orange)) +
  scale_shape_manual(values = c(prostate = 16, renal = 17)) +
  labs(title = "Complete exact block contributions", subtitle = "Raw patient-samples with median and IQR; D/I/T are anatomical categories",
       x = NULL, y = "Signed contribution to Axis1", colour = "Cohort", shape = "Cohort") + theme_pub + theme(legend.position = "bottom")
context5[attribution == "Unique drop-one delta R2", attribution := "Drop-one delta R2"]
context5[, attribution := factor(attribution, levels = c("Marginal R2", "Drop-one delta R2"))]
context5[, term_label := factor(term_label, levels = rev(c("Anatomical compartment", "Patient", "Cancer/accession", "Cell depth")))]
pS5b <- ggplot(context5, aes(R2, term_label, colour = attribution, shape = attribution)) +
  geom_segment(aes(x = 0, xend = R2, yend = term_label), position = position_dodge(width = 0.44), colour = "#C0C0C0", linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.44), size = 2.1) +
  scale_colour_manual(values = c(`Marginal R2` = pal$blue, `Drop-one delta R2` = pal$orange)) +
  scale_shape_manual(values = c(`Marginal R2` = 16, `Drop-one delta R2` = 17)) +
  labs(title = "Context attribution", subtitle = "Cancer/accession aliased with patient; descriptive only",
       x = "R-squared", y = NULL, colour = NULL, shape = NULL) + theme_pub +
  theme(legend.position = "bottom", plot.title = element_text(size = 8.5, face = "bold"), plot.subtitle = element_text(size = 6.3))
oep_stats <- oep5[, .(n = .N, median = median(Axis1)), by = cancer_code][order(median)]
oep5[, cancer_code := factor(cancer_code, levels = oep_stats$cancer_code)]
oep_stats[, cancer_code := factor(cancer_code, levels = levels(oep5$cancer_code))]
pS5c <- ggplot(oep5, aes(cancer_code, Axis1)) + geom_hline(yintercept = 0, linetype = 3, colour = "#777777") +
  geom_jitter(width = 0.10, height = 0, colour = pal$purple, size = 1.25, alpha = 0.82) +
  stat_summary(data = oep5[cancer_code %chin% oep_stats[n >= 3, cancer_code]], fun = median, geom = "crossbar", width = 0.48, linewidth = 0.4) +
  geom_text(data = oep_stats, aes(y = max(oep5$Axis1) + 0.12, label = paste0("n=", n)), size = 1.7) +
  labs(title = "OEP005136 projection", subtitle = "Raw source categories; median for n>=3",
       x = "Source-defined category", y = "Frozen Axis1") + theme_pub +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), plot.title = element_text(size = 8.5, face = "bold"),
        plot.subtitle = element_text(size = 6.3)) + scale_y_continuous(expand = expansion(mult = c(.04, .14)))
figS5 <- pS5a / (pS5b | pS5c) + plot_layout(heights = c(0.92, 1.08), widths = c(1.05, 0.95)) + plot_annotation(tag_levels = "A")
receipts[[5]] <- save_figure(figS5, "FigureS5_axis1_context_diagnostics", 7.2, 6.8)
write_source(contrib5, "FigureS5A_all_block_contributions.tsv"); write_source(context5, "FigureS5B_context_attribution.tsv")
write_source(oep5, "FigureS5C_oep_scores.tsv")

## S6: paired endpoint robustness ----------------------------------------------
rates6 <- read_dt(file.path(parent, "Figure7_subsampling_rule_rates.tsv"))
stability6 <- read_dt(file.path(parent, "Figure7_subsampling_stability.tsv"))
claims6 <- read_dt(file.path(parent, "Figure7_claim_boundary.tsv"))
rates6[, endpoint_short := factor(endpoint, levels = c("BM - matched normal bone", "BM - matched primary tumor"),
                                  labels = c("Matched normal bone", "Matched primary tumor"))]
pS6a <- ggplot(rates6, aes(rule_fraction, endpoint_short)) +
  geom_vline(xintercept = .80, linetype = 2, colour = pal$red, linewidth = .45) +
  geom_segment(aes(x = 0, xend = rule_fraction, yend = endpoint_short), colour = "#B0B0B0", linewidth = .65) +
  geom_point(aes(fill = rule_fraction >= .80), shape = 21, size = 2.7) +
  geom_text(aes(label = percent(rule_fraction, accuracy = .1)), nudge_y = .16, size = 1.9) +
  scale_fill_manual(values = c(`TRUE` = pal$teal, `FALSE` = pal$red), guide = "none") + scale_x_continuous(labels = percent, limits = c(0, 1.02)) +
  labs(title = "Frozen-rule retention", subtitle = "499 within-specimen resamples; dashed line is the 80% gate",
       x = "Resamples meeting sign-count rule", y = NULL) + theme_pub
stability6[, endpoint_short := factor(endpoint, levels = c("BM - normal bone", "BM - primary tumor"),
                                      labels = c("Matched normal bone", "Matched primary tumor"))]
q6 <- stability6[, .(median = median(median_difference), low95 = quantile(median_difference, .025), high95 = quantile(median_difference, .975),
                     low90 = quantile(median_difference, .05), high90 = quantile(median_difference, .95)), by = endpoint_short]
pS6b <- ggplot(q6, aes(median, endpoint_short)) + geom_vline(xintercept = 0, linetype = 2, colour = "#666666") +
  geom_errorbar(aes(xmin = low95, xmax = high95), orientation = "y", width = .10, linewidth = .55, colour = "#777777") +
  geom_errorbar(aes(xmin = low90, xmax = high90), orientation = "y", width = .18, linewidth = 2.0, colour = pal$purple) +
  geom_point(size = 2.4, colour = "black") +
  labs(title = "Resampling intervals", subtitle = "90% thick | 95% thin; same specimens",
       x = "Median paired Axis1 difference", y = NULL) + theme_pub
claims6[claim == "Mechanistic driver", claim := "Causal mechanism"]
claims6[, claim := factor(claim, levels = rev(claim))]
claim_cols <- c(Supported = "#DCEFEA", `Supported with boundary` = "#DCEAF5", `Not established` = "#E5E5E5", `Not tested` = "white")
pS6c <- ggplot(claims6, aes("Evidence", claim, fill = status)) +
  geom_tile(colour = "#BEBEBE", linewidth = .5) + geom_text(aes(label = paste(status, basis, sep = " | ")), size = 1.75) +
  scale_fill_manual(values = claim_cols, guide = "none") +
  labs(title = "Final evidence table", subtitle = "A favorable unpaired contrast cannot replace the failed matched-normal endpoint",
       x = NULL, y = NULL) + theme_heat + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
figS6 <- (pS6a | pS6b) / pS6c + plot_layout(heights = c(.90, 1.10)) + plot_annotation(tag_levels = "A")
receipts[[6]] <- save_figure(figS6, "FigureS6_paired_endpoint_robustness", 7.2, 5.8)
write_source(rates6, "FigureS6A_rule_retention.tsv"); write_source(stability6, "FigureS6B_resampling_values.tsv")
write_source(q6, "FigureS6B_resampling_intervals.tsv"); write_source(claims6, "FigureS6C_evidence_table.tsv")

## S7: internal functional annotation -----------------------------------------
hall7 <- read_dt(file.path(parent, "FigureS2_all_hallmark_effects.tsv"))
support7 <- read_dt(file.path(parent, "FigureS2_hallmark_layer_support.tsv"))
members7 <- read_dt(file.path(parent, "FigureS2_driver_recurrence_top15.tsv"))
component_map <- unique(support7[, .(pathway_label, redundancy_component)])
hall7 <- merge(hall7, component_map, by = "pathway_label", all.x = TRUE)
hall7[, row_label := sprintf("C%02d | %s", redundancy_component, pathway_label)]
hall_order <- unique(hall7[layer == "Patient-cluster"][order(estimate)]$row_label)
hall7[, row_label := factor(row_label, levels = hall_order)]
pS7a <- ggplot(hall7, aes(estimate, row_label, colour = layer)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#777777") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = .08, position = position_dodge(width = .40), linewidth = .34) +
  geom_point(position = position_dodge(width = .40), size = 1.05) +
  scale_colour_manual(values = c(`Patient-cluster` = pal$orange, `Sample-aggregated` = pal$blue)) +
  labs(title = "All 50 calibrated Hallmark estimates", subtitle = "C-prefixes denote redundancy groups",
       x = "Activity SD per Axis1 SD", y = NULL, colour = NULL) + theme_pub +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 4.8), plot.title = element_text(size = 8.5, face = "bold"))
support7 <- merge(support7, component_map, by = c("pathway_label", "redundancy_component"), all.x = TRUE)
support7[, row_label := sprintf("C%02d | %s", redundancy_component, pathway_label)]
support7[, row_label := factor(row_label, levels = hall_order)]
support7[, layer := factor(layer, levels = c("Patient-cluster", "Sample aggregate", "Matched null"))]
pS7b <- ggplot(support7, aes(layer, row_label, fill = supported)) +
  geom_tile(colour = "#D8D8D8", linewidth = .20) + geom_text(aes(label = ifelse(supported, "YES", "-")), size = 1.25) +
  scale_fill_manual(values = c(`TRUE` = "#DCEFEA", `FALSE` = "white"), guide = "none") +
  labs(title = "Layer support", subtitle = "Distinct from effect size", x = NULL, y = NULL) +
  theme_heat + theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.text.y = element_blank(),
                     plot.title = element_text(size = 8.5, face = "bold"), plot.subtitle = element_text(size = 6.3))
members7[, gene := factor(gene, levels = rev(gene))]
pS7c <- ggplot(members7, aes(n_hallmarks, gene)) +
  geom_segment(aes(x = 0, xend = n_hallmarks, yend = gene), colour = "#A8A8A8", linewidth = .55) +
  geom_point(aes(size = median_abs_adjusted_beta), colour = pal$red) +
  scale_size_continuous(range = c(1.5, 3.4), name = "Median |beta|") +
  labs(title = "Recurrently represented genes across calibrated Hallmark sets",
       subtitle = "Membership recurrence is descriptive and does not identify regulatory drivers",
       x = "Calibrated Hallmark sets containing gene", y = NULL) + theme_pub + theme(legend.position = "bottom")
figS7 <- (pS7a | pS7b) / pS7c + plot_layout(heights = c(1.72, .58), widths = c(1.38, .62)) + plot_annotation(tag_levels = "A")
receipts[[7]] <- save_figure(figS7, "FigureS7_internal_functional_annotation", 7.2, 10.2)
write_source(hall7, "FigureS7A_all_hallmark_effects.tsv"); write_source(support7, "FigureS7B_support_matrix.tsv")
write_source(members7, "FigureS7C_recurrent_member_genes.tsv")

## S8: GSE266330 directional-test audit ---------------------------------------
gse8 <- read_dt(file.path(parent, "FigureS3_GSE266330_projectable_scores.tsv"))
loco8 <- read_dt(file.path(parent, "FigureS3_GSE266330_leave_one_origin_out.tsv"))
tech8 <- read_dt(file.path(parent, "FigureS3_GSE266330_technical_diagnostics.tsv"))
stats8 <- gse8[, .(n = .N, median = median(Axis1)), by = origin][order(median)]
gse8[, origin := factor(origin, levels = stats8$origin)]
stats8[, origin := factor(origin, levels = levels(gse8$origin))]
pS8a <- ggplot(gse8, aes(origin, Axis1)) + geom_jitter(width = .10, height = 0, colour = pal$purple, size = 1.3, alpha = .82) +
  stat_summary(data = gse8[origin %chin% stats8[n >= 3, origin]], fun = median, geom = "crossbar", width = .5, linewidth = .4) +
  geom_text(data = stats8, aes(y = max(gse8$Axis1) + .12, label = paste0("n=", n)), size = 1.7) +
  labs(title = "Every projectable donor/patient", subtitle = "Source categories are descriptive; medians only for n>=3",
       x = NULL, y = "Frozen Axis1") + theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(.04, .14)))
loco8[, dropped_origin := factor(dropped_origin, levels = rev(dropped_origin))]
pS8b <- ggplot(loco8, aes(hodges_lehmann_shift, dropped_origin)) + geom_vline(xintercept = 0, linetype = 2, colour = "#777777") +
  geom_segment(aes(x = 0, xend = hodges_lehmann_shift, yend = dropped_origin), colour = "#B0B0B0", linewidth = .55) +
  geom_point(colour = pal$teal, size = 2.0) +
  labs(title = "Leave-one-origin-out test", x = "Hodges-Lehmann metastasis-minus-control shift", y = "Dropped origin") + theme_pub
pS8c <- ggplot(tech8, aes(value, Axis1)) + geom_point(size = 1.1, alpha = .75, colour = pal$purple) +
  geom_smooth(method = "lm", se = FALSE, colour = "#555555", linewidth = .45) + facet_wrap(~metric, scales = "free_x", nrow = 1) +
  labs(title = "Continuous technical diagnostics", subtitle = "Single non-projectable case remains in source data",
       x = NULL, y = "Axis1") + theme_pub
figS8 <- (pS8a | pS8b) / pS8c + plot_layout(heights = c(1.05, .80), widths = c(1.30, .70)) + plot_annotation(tag_levels = "A")
receipts[[8]] <- save_figure(figS8, "FigureS8_gse266330_directional_audit", 7.2, 6.7)
write_source(gse8, "FigureS8A_all_projectable_scores.tsv"); write_source(loco8, "FigureS8B_leave_one_origin_out.tsv")
write_source(tech8, "FigureS8C_technical_diagnostics.tsv")

## S9: complete spatial sensitivity archive -----------------------------------
overlay9 <- read_dt(file.path(phase_a, "source_data/Figure6A_C_spot_overlay_coordinates_and_scores.tsv.gz"))
assets9 <- read_dt(file.path(phase_a, "provenance/Figure6_histology_asset_index.tsv"))
bars9 <- read_dt(file.path(phase_a, "provenance/Figure6_scale_bar_audit.tsv"))
scale9 <- read_dt(file.path(phase_a, "provenance/Figure6_shared_color_scale_audit.tsv"))
neigh9 <- read_dt(file.path(phase_a, "source_data/Figure6F_neighborhood_class_effects.tsv"))
effects9 <- read_dt(file.path(parent, "FigureS4_spatial_grid_all_effects.tsv"))
summary9 <- read_dt(file.path(parent, "FigureS4_spatial_grid_summary.tsv"))
shared_limit9 <- unique(scale9$shared_abs_q98_limit)
overlay9[, score_difference := full_axis1_display - malignant_excluded_axis1_display]
diff_limit9 <- unname(quantile(abs(overlay9$score_difference), .98, na.rm = TRUE))

theme_spatial <- theme_void(base_size = 7, base_family = "Arial") +
  theme(plot.title = element_text(face = "bold", size = 6.8, hjust = .5), legend.position = "bottom",
        legend.title = element_text(size = 6.4), legend.text = element_text(size = 6.0), plot.margin = margin(1, 1, 1, 1))
spatial_assets <- list()
for (i in seq_len(nrow(assets9))) {
  s <- assets9$sample[[i]]; img <- readPNG(track(assets9$image_relative_path[[i]])); d <- overlay9[sample == s]
  xr <- range(d$plot_x_lowres); yr <- range(d$plot_y_lowres); px <- diff(xr) * .035; py <- diff(yr) * .035
  spatial_assets[[s]] <- list(image = img * .88 + .12, width = dim(img)[2], height = dim(img)[1],
                              crop = c(xmin = max(0, xr[1] - px), xmax = min(dim(img)[2], xr[2] + px),
                                       ymin = max(0, yr[1] - py), ymax = min(dim(img)[1], yr[2] + py)),
                              bar = bars9[sample == s, pixels_per_500um_lowres])
}
spatial_map <- function(sample_value, difference = FALSE) {
  a <- spatial_assets[[sample_value]]; d <- copy(overlay9[sample == sample_value]); cr <- a$crop
  d[, score := if (difference) score_difference else full_axis1_display]
  lim <- if (difference) diff_limit9 else shared_limit9
  dx <- cr[["xmax"]] - cr[["xmin"]]; dy <- cr[["ymax"]] - cr[["ymin"]]
  bx0 <- cr[["xmin"]] + .055 * dx; bx1 <- bx0 + a$bar; by <- cr[["ymin"]] + .045 * dy
  ggplot() + annotation_raster(a$image, xmin = 0, xmax = a$width, ymin = 0, ymax = a$height) +
    geom_point(data = d, aes(plot_x_lowres, plot_y_lowres, fill = score), shape = 21, colour = alpha("white", 0), stroke = 0, size = .54, alpha = .70) +
    geom_point(data = d[source_annotated_tumor == TRUE], aes(plot_x_lowres, plot_y_lowres), shape = 21, fill = NA, colour = "#333333", stroke = .10, size = .63) +
    annotate("segment", x = bx0, xend = bx1, y = by, yend = by, linewidth = .55, colour = "black") +
    annotate("text", x = (bx0 + bx1) / 2, y = by + .025 * dy, label = "500 um", size = 1.35, fontface = "bold") +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
                         limits = c(-lim, lim), oob = squish, name = if (difference) "Full - malignant-excluded" else "Full Axis1") +
    coord_fixed(xlim = cr[c("xmin", "xmax")], ylim = cr[c("ymin", "ymax")], expand = FALSE) +
    labs(title = sample_value) + theme_spatial
}
maps_full <- wrap_plots(lapply(assets9$sample, spatial_map, difference = FALSE), ncol = 2, guides = "collect") & theme(legend.position = "bottom")
maps_diff <- wrap_plots(lapply(assets9$sample, spatial_map, difference = TRUE), ncol = 2, guides = "collect") & theme(legend.position = "bottom")
neigh9[, sample := factor(sample, levels = assets9$sample)]
neigh9[, class_label := factor(class_label, levels = rev(sort(unique(class_label))))]
pS9c <- ggplot(neigh9, aes(sample, class_label, fill = estimate)) +
  geom_tile(colour = "white", linewidth = .45) + geom_text(aes(label = sprintf("%.2f", estimate)), size = 1.65) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "Boundary - distal") +
  labs(title = "Complete nine-class neighborhood-associated composition", subtitle = "Section-specific estimates only; no interaction or pooled-section claim",
       x = "Linked section", y = NULL) + theme_heat + theme(axis.text.x = element_text(angle = 25, hjust = 1))
effects9[, layer_label := factor(layer_label, levels = c("Full Axis1", "Malignant-excluded Axis1", "One-ring RCTD proxy"))]
effects9[, sample := factor(sample, levels = assets9$sample)]
effects9[, width_label := factor(block_width, levels = c(15, 20, 25))]
effects9[, k_label := factor(k, levels = c(4, 6, 8))]
lim_grid <- unname(quantile(abs(effects9$estimate), .98, na.rm = TRUE))
pS9d <- ggplot(effects9, aes(width_label, k_label, fill = ifelse(decision_evaluable, estimate, NA_real_))) +
  geom_tile(colour = "white", linewidth = .42) +
  geom_tile(data = effects9[k == 6 & block_width == 20], fill = NA, colour = "black", linewidth = .65) +
  geom_text(data = effects9[decision_evaluable == FALSE], label = "NE", size = 1.7, colour = "#666666") +
  facet_grid(layer_label ~ sample) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-lim_grid, lim_grid), oob = squish, na.value = "#F0F0F0", name = "Effect") +
  labs(title = "Complete k-by-block sensitivity matrix", subtitle = "Black outline: primary k=6, width=20 | NE: insufficient spatial blocks",
       x = "Spatial-block width", y = "Symmetric graph k") + theme_heat

pS9a <- wrap_elements(full = maps_full + plot_annotation(title = "Full Axis1 on all four linked sections",
  theme = theme(plot.title = element_text(family = "Arial", size = 9.2, face = "bold"))))
pS9b <- wrap_elements(full = maps_diff + plot_annotation(title = "Full minus malignant-excluded Axis1 on all four sections",
  theme = theme(plot.title = element_text(family = "Arial", size = 9.2, face = "bold"))))
figS9 <- pS9a / pS9b / pS9c / pS9d + plot_layout(heights = c(1.12, 1.12, .72, 1.35)) + plot_annotation(tag_levels = "A")
receipts[[9]] <- save_figure(figS9, "FigureS9_spatial_full_sensitivity", 7.2, 14.0)
write_source(overlay9, "FigureS9A_B_spatial_scores.tsv.gz"); write_source(neigh9, "FigureS9C_full_neighborhood_matrix.tsv")
write_source(summary9, "FigureS9D_grid_summary.tsv"); write_source(effects9, "FigureS9D_all_grid_effects.tsv")

## Contracts, legends and provenance ------------------------------------------
contract <- data.table(
  figure = paste0("Figure S", 1:9),
  purpose = c("Cohort/QC/annotation", "Composition robustness", "T/NK transcription robustness",
              "Fixed-program/ecotype diagnostics", "Axis1 construction/context diagnostics",
              "Independent paired endpoint robustness", "Internal functional annotation",
              "GSE266330 directional audit", "Spatial full sensitivity"),
  main_link = c("Figure 1", "Figure 2", "Figure 3", "Figure 4", "Figure 5", "Figure 5", "Figure 5", "Figure 5", "Figure 6"),
  prohibited_inference = c("Marker QC is not discovery evidence", "CLR effects are relative", "Persistence is not external validation",
                           "Stable computation is not transportability", "Context attribution is descriptive", "499 resamples are not biological replicates",
                           "Hallmarks are not mechanisms", "Origin groups are not powered comparisons", "Linked sections are not independent animals")
)
fwrite(contract, file.path(admin_dir, "GATE12AF_SUPPLEMENT_CONTRACT.tsv"), sep = "\t")

legends <- c(
  "# Gate12AF supplementary figure legends", "",
  "## Supplementary Figure S1. Cohort, quality-control and annotation archive", "",
  "A, Frozen dataset roles with biological units. B, Discovery specimen availability. C, Canonical marker audit; Unassigned is explicitly a QC category. D, Full ten-class composition across 42 samples; stacked bars are descriptive.", "",
  "## Supplementary Figure S2. Paired compositional robustness", "",
  "A, Empirical-logit versus paired CLR effects on the same scale. B, CLR estimates for all 16 states; all 16 directions were retained in all six prespecified sensitivity settings.", "",
  "## Supplementary Figure S3. T/NK transcriptional robustness", "",
  "A, Persistent-over-eligible gene and Hallmark results after frozen sensitivity filters. B, Version-frozen Hallmark effects that required cohort-direction and patient-deletion support. Persistence is internal discovery robustness, not external validation.", "",
  "## Supplementary Figure S4. Fixed-representation failure diagnostics", "",
  "A, Consensus diagnostics with k=3 in black and frozen thresholds in red. B, Exact 10,000-permutation ARI null. C, Discrete evidence table separating discovery stability from external direction and structure/agreement.", "",
  "## Supplementary Figure S5. Axis1 context diagnostics", "",
  "A, Complete exact block contributions. B, Descriptive context attribution; cancer/accession is aliased with patient. C, Raw OEP005136 technical-projection scores by descriptive source category; medians are shown only for categories with at least three patients.", "",
  "## Supplementary Figure S6. Independent paired endpoint robustness", "",
  "A, Frozen-rule retention in 499 within-specimen resamples. B, Median and 90%/95% resampling intervals; the iterations reuse the same specimens and are not biological replicates. C, Final evidence table retaining the failed matched-normal endpoint.", "",
  "## Supplementary Figure S7. Internal functional annotation", "",
  "A, All 50 Hallmark estimates with redundancy-component labels. B, Calibration support shown separately from effect magnitude. C, Recurrently represented member genes across calibrated Hallmark sets. These panels do not establish independent pathways, causal targets or regulatory drivers.", "",
  "## Supplementary Figure S8. GSE266330 directional-test audit", "",
  "A, Every projectable donor/patient by descriptive source category. B, Leave-one-origin-out Hodges-Lehmann shifts. C, Continuous technical diagnostics. The comparison remains unpaired and source-confounded.", "",
  "## Supplementary Figure S9. Complete spatial sensitivity archive", "",
  "A, Full Axis1 maps in all four linked sections. B, Full-minus-malignant-excluded maps. C, Complete nine-class neighborhood-associated composition matrix. D, Complete graph-k by block-width effect matrix; NE marks insufficient spatial blocks and the primary setting is outlined. Sections are not pooled as independent animals."
)
writeLines(legends, file.path(out_dir, "GATE12AF_SUPPLEMENTARY_FIGURE_LEGENDS.md"))
writeLines(capture.output(sessionInfo()), file.path(prov_dir, "sessionInfo.txt"))

manifest <- rbindlist(receipts)
fwrite(manifest, file.path(prov_dir, "GATE12AF_SUPPLEMENT_FIGURE_MANIFEST.tsv"), sep = "\t")
used_inputs <- unique(c(used_inputs, "scripts/build_gate12af_supplementary_redesign.R"))
input_manifest <- data.table(relative_path = used_inputs, bytes = file.info(used_inputs)$size,
                             sha256 = vapply(used_inputs, sha256_file, character(1)))
setorder(input_manifest, relative_path)
fwrite(input_manifest, file.path(prov_dir, "GATE12AF_SUPPLEMENT_INPUT_MANIFEST.tsv"), sep = "\t")
cat("GATE12AF_SUPPLEMENT_STATUS=GENERATED_PENDING_AUDIT\n")
cat("GATE12AF_SUPPLEMENT_FIGURES=9\n")
cat("GATE12AF_SUPPLEMENT_OUTPUT=", rel_path(out_dir), "\n", sep = "")
