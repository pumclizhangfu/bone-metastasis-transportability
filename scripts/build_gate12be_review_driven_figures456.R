#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(png)
  library(scales)
  library(grid)
})

options(stringsAsFactors = FALSE, warn = 1)
set.seed(20260812L)
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
out <- file.path(root, "results", "gate12be_review_driven_redesign")
fig_out <- file.path(out, "figures", "main")
src_out <- file.path(out, "source_data", "main")
admin_out <- file.path(out, "admin")
for (d in c(fig_out, src_out, admin_out)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

read_rel <- function(path, ...) {
  full <- file.path(root, path)
  if (!file.exists(full)) stop("Missing input: ", path)
  fread(full, ...)
}
write_source <- function(x, stem) fwrite(x, file.path(src_out, stem), sep = "\t", quote = FALSE, na = "NA",
                                         compress = if (endsWith(stem, ".gz")) "gzip" else "none")
cairo_png <- function(filename, width, height, bg = "white", ...) {
  grDevices::png(filename, width = width, height = height, units = "in", res = 450,
                 type = "cairo", bg = bg, ...)
}
save_figure <- function(plot, stem, height) {
  ggsave(file.path(fig_out, paste0(stem, ".png")), plot, width = 7.2, height = height,
         units = "in", dpi = 450, device = cairo_png, bg = "white", limitsize = FALSE)
  ggsave(file.path(fig_out, paste0(stem, ".pdf")), plot, width = 7.2, height = height,
         units = "in", device = cairo_pdf, bg = "white", limitsize = FALSE)
}

blue <- "#3973B9"; orange <- "#E17C2B"; green <- "#009E73"; purple <- "#7047A3"
red <- "#C44E52"; grey <- "#BFC3C7"; dark <- "#30343B"
theme_pub <- theme_classic(base_size = 8.7, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 9.4, face = "bold", hjust = 0, margin = margin(b = 2)),
    plot.subtitle = element_text(size = 7.3, colour = "#525A61", hjust = 0, lineheight = .96, margin = margin(b = 3)),
    axis.title = element_text(size = 8.0), axis.text = element_text(size = 7.2, colour = "#30343B"),
    strip.text = element_text(size = 7.3, face = "bold"),
    strip.background = element_rect(fill = "#F3F4F5", colour = NA),
    legend.title = element_text(size = 7.2), legend.text = element_text(size = 6.9),
    legend.key.height = unit(2.8, "mm"), legend.key.width = unit(2.8, "mm"),
    panel.grid = element_blank(), plot.tag = element_text(size = 11.5, face = "bold"),
    plot.margin = margin(4, 5, 4, 5)
  )
title_theme <- theme(
  plot.title = element_text(family = "Arial", size = 13.2, face = "bold", hjust = 0),
  plot.subtitle = element_text(family = "Arial", size = 8.0, colour = "#4B5563", hjust = 0),
  plot.tag = element_text(family = "Arial", size = 11.5, face = "bold"),
  plot.margin = margin(7, 7, 6, 7)
)

## Figure 4 -----------------------------------------------------------------
shortlist <- read_rel("results/gate12c_sender_receiver/gate12c_provisional_shortlist.tsv")
ligand_activity <- read_rel("results/gate12c_sender_receiver/nichenet_ligand_activity.tsv")
ligand_targets <- read_rel("results/gate12c_sender_receiver/nichenet_top_ligand_target_links.tsv")
cellchat <- read_rel("results/gate12c_sender_receiver/cellchat_condition_communications.tsv.gz")
sensitivity <- read_rel("results/gate12c_external_validation/external_threshold_sensitivity.tsv")
decision <- read_rel("results/gate12c_external_validation/gate12c_external_axis_decision.tsv")

top_ligands <- unique(rbindlist(list(
  ligand_activity[order(nichenet_rank)][1:min(12L, .N)],
  ligand_activity[ligand %in% shortlist$ligand]
), use.names = TRUE))
setorder(top_ligands, nichenet_rank)
top_ligands[, shortlisted := ligand %in% shortlist$ligand]
top_ligands[, ligand := factor(ligand, levels = rev(as.character(ligand)))]
p4a <- ggplot(top_ligands, aes(pearson, ligand)) +
  geom_col(aes(fill = shortlisted), width = .70) +
  geom_text(aes(label = paste0("rank ", nichenet_rank)), hjust = -0.12,
            size = 2.05, colour = dark) +
  scale_fill_manual(values = c(`TRUE` = orange, `FALSE` = "#AAB2BA"),
                    labels = c(`TRUE` = "Shortlisted pair", `FALSE` = "Top-rank context"),
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, .18))) +
  labs(title = "NicheNet-predicted ligand activity",
       subtitle = "Top 12 plus all shortlisted ligands; labels give global rank",
       x = "Ligand-target Pearson correlation", y = NULL) + theme_pub +
  theme(legend.position = "bottom", legend.justification = "left")

short_chat <- merge(cellchat, shortlist[, .(ligand, receptor, best_source, best_target, axis)],
                    by = c("ligand", "receptor"), all = FALSE, sort = FALSE)
short_chat <- short_chat[source == best_source & target == best_target,
                         .(prob = max(prob), pval = min(pval)), by = .(axis, accession, compartment)]
short_grid <- CJ(axis = shortlist$axis, accession = c("GSE143791", "GSE202813"),
                 compartment = c("distal", "involved", "tumor"), unique = TRUE)
short_chat <- merge(short_grid, short_chat, by = c("axis", "accession", "compartment"), all.x = TRUE)
short_chat[, compartment := factor(compartment, c("distal", "involved", "tumor"), c("Distal", "Involved", "Tumor"))]
short_chat[, cohort := factor(accession, c("GSE143791", "GSE202813"), c("Prostate", "Renal"))]
short_chat[, axis := factor(axis, levels = rev(shortlist$axis))]
p4b <- ggplot(short_chat, aes(compartment, axis)) +
  geom_point(data = short_chat[is.na(prob)], shape = 4, size = 2.6,
             colour = "#AAB2BA", stroke = .65) +
  geom_point(data = short_chat[!is.na(prob)], aes(fill = prob),
             shape = 21, size = 3.8, colour = dark, stroke = .22) +
  facet_wrap(~cohort, nrow = 1, drop = FALSE) +
  scale_fill_viridis_c(option = "C", direction = -1, trans = "sqrt",
                       limits = c(0, .05), oob = squish,
                       name = "Inferred CellChat\nprobability") +
  scale_x_discrete(drop = FALSE) +
  labs(title = "Expression-supported ligand-receptor pairs",
       subtitle = "State-matched CellChat support; x = not retained", x = NULL, y = NULL) + theme_pub +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "right") +
  guides(fill = guide_colourbar(title.position = "top", barheight = unit(17, "mm")))

candidate_targets <- ligand_targets[ligand %in% shortlist$ligand]
target_grid <- CJ(ligand = shortlist$ligand, target = unique(candidate_targets$target), unique = TRUE)
candidate_targets <- merge(target_grid, candidate_targets, by = c("ligand", "target"), all.x = TRUE)
candidate_targets[is.na(weight), weight := 0]
candidate_targets[, ligand := factor(ligand, levels = shortlist$ligand)]
candidate_targets[, target := factor(target, levels = rev(unique(target)))]
p4c <- ggplot(candidate_targets, aes(ligand, target, fill = weight)) +
  geom_tile(colour = "white", linewidth = .45) +
  scale_fill_gradient(low = "#F7F7F7", high = "#B2182B", name = "Regulatory\npotential") +
  labs(title = "NicheNet regulatory-potential screen", subtitle = "Targets in the frozen receiver programme",
       x = "Ligand", y = "Target") + theme_pub +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "right")

ligand_effects <- rbindlist(list(
  shortlist[, .(ligand, cohort = "Prostate", log2fc = ligand_log2fc_GSE143791, q_value = ligand_q_GSE143791)],
  shortlist[, .(ligand, cohort = "Renal", log2fc = ligand_log2fc_GSE202813, q_value = ligand_q_GSE202813)]
))
ligand_effects[, cohort := factor(cohort, levels = c("Prostate", "Renal"))]
ligand_effects[, ligand := factor(ligand, levels = rev(shortlist$ligand))]
ligand_effects[, fdr_strength := -log10(pmax(q_value, 1e-12))]
p4d <- ggplot(ligand_effects, aes(cohort, ligand)) +
  geom_point(aes(size = fdr_strength, fill = log2fc), shape = 21, colour = dark, stroke = .18) +
  scale_fill_gradient(low = "#FCEBEA", high = "#B94450", name = "Ligand\nlog2FC") +
  scale_size_continuous(range = c(1.6, 4.6), name = expression(-log[10](FDR))) +
  labs(title = "Patient-pseudobulk ligand effects", subtitle = "Discovery-cohort patient-level expression evidence",
       x = NULL, y = NULL) + theme_pub +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.box = "horizontal", legend.justification = "left",
        legend.spacing.x = unit(2, "mm"), legend.margin = margin(0),
        plot.subtitle = element_text(size = 7.0)) +
  guides(
    size = guide_legend(order = 1, title.position = "top", nrow = 1),
    fill = guide_colourbar(order = 2, title.position = "top",
                           barwidth = unit(18, "mm"), barheight = unit(2.4, "mm"))
  )

sensitivity_plot <- copy(sensitivity)
sensitivity_plot[, detection_label := factor(
  paste0(round(100 * min_detect), "%"),
  levels = c("2%", "5%", "10%")
)]
sensitivity_plot[, min_cells_label := factor(min_cells, levels = c(10, 20, 30))]
sensitivity_plot[, pair_text := sub(" -> ", "-", axis, fixed = TRUE)]
sensitivity_plot[axis == "CXCL16 -> CXCR6", pair_text := "CXCL16-CXCR6 (primary)"]
sensitivity_plot[, pair_label := factor(
  pair_text,
  levels = c("SPP1-CD44", "CXCL16-CXCR6 (primary)", "CCL4-CCR5")
)]
sensitivity_plot[, selected := axis == scenario_winner]
sensitivity_plot[, support_label := sprintf("%.0f%%", 100 * min_core)]
sensitivity_plot[, text_colour := fifelse(min_core >= .80, "white", dark)]

p4e_grid <- ggplot(sensitivity_plot, aes(min_cells_label, pair_label, fill = min_core)) +
  geom_tile(colour = "white", linewidth = .58) +
  geom_tile(data = sensitivity_plot[selected == TRUE], fill = NA,
            colour = dark, linewidth = .88) +
  geom_text(aes(label = support_label, colour = text_colour), size = 2.05) +
  facet_grid(. ~ detection_label, scales = "free_x", space = "free_x") +
  scale_fill_viridis_c(option = "C", direction = -1, limits = c(.5, 1), oob = squish,
                       labels = percent_format(accuracy = 1), name = "Minimum core\nsupport") +
  scale_colour_identity() +
  labs(title = "Threshold sensitivity of external support",
       subtitle = "Text: minimum support; outline: selected pair",
       x = "Minimum sender/receiver cells", y = NULL) + theme_pub +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        panel.spacing.x = unit(1.4, "mm"), legend.position = "right",
        legend.key.height = unit(13, "mm"), plot.margin = margin(4, 5, 4, 4)) +
  guides(fill = guide_colourbar(title.position = "top", barheight = unit(15, "mm")))

selection_stability <- sensitivity_plot[, .(
  scenarios_selected = sum(selected),
  total_scenarios = .N,
  retention = mean(selected)
), by = .(axis, pair_label)]
p4e_frequency <- ggplot(selection_stability, aes(retention, pair_label)) +
  geom_vline(xintercept = .75, linetype = 2, linewidth = .36, colour = "#777777") +
  geom_col(aes(fill = axis == "CXCL16 -> CXCR6"), width = .58) +
  geom_text(aes(label = paste0(scenarios_selected, "/", total_scenarios)), hjust = -0.18,
            size = 2.55, fontface = "bold", colour = dark) +
  annotate("text", x = .75, y = 3.45, label = "75% retention rule",
           hjust = 1.05, size = 2.15, colour = "#666666") +
  scale_fill_manual(values = c(`TRUE` = purple, `FALSE` = "#AAB2BA"), guide = "none") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, .25, .50, .75, 1),
                     labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, .02))) +
  labs(title = "Selection stability", subtitle = "No pair met the prespecified rule",
       x = "Scenarios selected", y = NULL) +
  theme_pub +
  theme(axis.text.y = element_text(size = 6.9), plot.margin = margin(4, 4, 4, 8))
p4e <- wrap_elements(full = p4e_grid | p4e_frequency + plot_layout(widths = c(1.55, .45)))

figure4 <- (p4a | p4b) / (p4c | p4d) / p4e +
  plot_layout(heights = c(.90, .78, 1.25), widths = c(.92, 1.08)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(family = "Arial", size = 11.5, face = "bold"),
      plot.margin = margin(4, 7, 5, 7)
    )
  )
save_figure(figure4, "Figure4_review_driven_communication_hypotheses", 8.4)

## Figure 5 -----------------------------------------------------------------
gse_scores <- read_rel("results/gate8a_gse266330_server_v2/sample_lineage_scores.tsv")
gse_threshold <- read_rel("results/gate8a_gse266330_server_v2/cell_threshold_sensitivity.tsv")
oep_pairs <- read_rel("results/gate8b_oep005136/paired_differences.tsv")
oep_threshold <- read_rel("results/gate8b_oep005136/cell_threshold_sensitivity.tsv")
axis_src <- file.path(root, "results", "gate12al_figure1_integration", "source_data")
axis_recon <- read_rel(paste0("results/gate12ad_figure_restructure/phase_b_main_figures/",
                              "source_data/Figure5B_block_lopo_predictions.tsv"))
axis_metrics <- fread(file.path(axis_src, "Figure5C_reconstruction_metrics.tsv"))
axis_external <- fread(file.path(axis_src, "Figure5D_unpaired_external_scores.tsv"))
axis_paired <- fread(file.path(axis_src, "Figure5E_paired_endpoints.tsv"))
axis_model <- readRDS(file.path(root, "results/gate12g_ecological_full/gate12g_ecological_program_full.rds"))
axis_cells <- read_rel(paste0("results/gate12z_major_revision/analysis/gse266330_axis1/",
                             "gse266330_gate12z_cell_assignments.tsv.gz"))

lineage_labels <- c(CD4_T = "CD4 T", CD8_CTL = "CD8/CTL", NK_NKT = "NK/NKT")
gse_common <- copy(gse_scores[signature == "common_81" & eligible_20 == TRUE])
gse_common[, lineage_label := factor(lineage_labels[lineage], levels = lineage_labels)]
gse_common[, condition_label := factor(condition, c("healthy_bm", "bone_metastasis"),
                                        c("Healthy marrow", "Bone metastasis"))]
gse_counts <- gse_common[, .(n = .N), by = .(lineage_label, condition_label)]
gse_counts[, label_y := max(gse_common$score, na.rm = TRUE) + .008]
p5a <- ggplot(gse_common, aes(condition_label, score, fill = condition_label)) +
  geom_boxplot(width = .58, outlier.shape = NA, linewidth = .38, alpha = .34) +
  geom_point(position = position_jitter(width = .10), size = .95, alpha = .76, stroke = 0) +
  geom_text(data = gse_counts, aes(condition_label, label_y, label = paste0("n=", n)),
            inherit.aes = FALSE, size = 2.0, colour = dark) +
  facet_wrap(~lineage_label, nrow = 1) +
  scale_fill_manual(values = c(`Healthy marrow` = grey, `Bone metastasis` = red), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(.03, .15))) +
  labs(title = "External 81-gene programme", x = NULL, y = "Programme score") + theme_pub +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

oep_common <- copy(oep_pairs[signature == "common_81"])
oep_common[, lineage_label := factor(lineage_labels[lineage], levels = lineage_labels)]
oep_common[, contrast_label := factor(contrast, c("bm_vs_normal_bone", "bm_vs_primary"),
                                      c("Normal bone", "Primary tumour"))]
oep_common[, pair_id := paste0(patient_id, fifelse(contrast == "bm_vs_normal_bone", "-N", "-P"))]
pair_order <- unique(oep_common[order(contrast, patient_id), pair_id])
oep_pair_grid <- CJ(pair_id = pair_order, lineage = names(lineage_labels), unique = TRUE)
oep_pair_heat <- merge(
  oep_pair_grid,
  oep_common[, .(pair_id, patient_id, cancer, contrast, contrast_label, lineage,
                 lineage_label, difference, control_score, bone_metastasis_score)],
  by = c("pair_id", "lineage"), all.x = TRUE, sort = FALSE
)
oep_pair_heat[, pair_id := factor(pair_id, levels = pair_order)]
oep_pair_heat[, lineage_label := factor(lineage_labels[lineage], levels = rev(lineage_labels))]
pair_meta <- unique(oep_common[, .(pair_id, patient_id, cancer, contrast, contrast_label)])
pair_meta[, pair_id := factor(pair_id, levels = pair_order)]
pair_meta[, contrast_bar := fifelse(contrast == "bm_vs_normal_bone", "Normal bone", "Primary tumour")]
pair_meta[, cancer_bar := cancer]
pair_annotation <- melt(
  pair_meta[, .(pair_id, Comparison = contrast_bar, Cancer = cancer_bar)],
  id.vars = "pair_id", variable.name = "annotation", value.name = "annotation_value"
)
pair_annotation[, annotation := factor(annotation, c("Cancer", "Comparison"))]
pair_annotation[, annotation_colour := fcase(
  annotation_value == "LUCA", green,
  annotation_value == "PRAD", purple,
  annotation_value == "Normal bone", "#8FA3B8",
  annotation_value == "Primary tumour", "#D7A26D",
  default = "#DDDDDD"
)]
pair_annotation[, annotation_label := fcase(
  annotation_value == "Normal bone", "Normal",
  annotation_value == "Primary tumour", "Primary",
  default = annotation_value
)]
pair_annotation[, annotation_text_colour := fifelse(
  annotation_value %in% c("LUCA", "PRAD"), "white", dark
)]
pair_effect_limit <- max(abs(oep_pair_heat$difference), na.rm = TRUE)
p5b_annotation <- ggplot(pair_annotation, aes(pair_id, annotation, fill = annotation_colour)) +
  geom_tile(colour = "white", linewidth = .35) +
  geom_vline(xintercept = 2.5, colour = "white", linewidth = 1.0) +
  geom_text(aes(label = annotation_label, colour = annotation_text_colour),
            size = 1.72, fontface = "bold") +
  scale_fill_identity() +
  scale_colour_identity() +
  labs(x = NULL, y = NULL) +
  theme_void(base_family = "Arial", base_size = 7.2) +
  theme(axis.text.y = element_text(colour = dark, size = 6.2, hjust = 1),
        axis.text.x = element_blank(), plot.margin = margin(0, 5, 0, 4))
p5b_heat <- ggplot(oep_pair_heat, aes(pair_id, lineage_label, fill = difference)) +
  geom_tile(colour = "white", linewidth = .58) +
  geom_vline(xintercept = 2.5, colour = "white", linewidth = 1.0) +
  geom_text(aes(label = fifelse(is.na(difference), "NE", sprintf("%+.4f", difference))),
            size = 1.58, colour = dark) +
  scale_fill_gradient2(low = blue, mid = "white", high = red, midpoint = 0,
                       limits = c(-pair_effect_limit, pair_effect_limit), oob = squish,
                       na.value = "#EFEFEF", name = expression(Delta*" programme score")) +
  scale_x_discrete(labels = setNames(pair_meta$patient_id, as.character(pair_meta$pair_id)), drop = FALSE) +
  labs(x = NULL, y = NULL) + theme_pub +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(size = 6.3), axis.text.y = element_text(size = 6.4),
        legend.position = "right", legend.key.height = unit(10, "mm"),
        plot.margin = margin(0, 5, 3, 4)) +
  guides(fill = guide_colourbar(title.position = "top", barheight = unit(11, "mm")))
p5b <- wrap_elements(full =
  (p5b_annotation / p5b_heat + plot_layout(heights = c(.37, 1.00))) +
    plot_annotation(title = "Paired 81-gene programme") &
    theme(plot.title = element_text(family = "Arial", size = 9.4, face = "bold", hjust = 0,
                                    margin = margin(l = 12, b = 2)),
          plot.margin = margin(4, 5, 3, 5))
)

threshold_plot <- rbindlist(list(
  gse_threshold[, .(panel_text = "GSE266330", min_cells, lineage,
                    effect_value = effect, unit_label = sprintf("%d/%d", n_bm, n_control))],
  oep_threshold[contrast == "bm_vs_normal_bone",
                .(panel_text = "OEP: normal bone", min_cells, lineage,
                  effect_value = median_difference, unit_label = sprintf("%d/%d", n_positive, n_pairs))],
  oep_threshold[contrast == "bm_vs_primary",
                .(panel_text = "OEP: primary tumour", min_cells, lineage,
                  effect_value = median_difference, unit_label = sprintf("%d/%d", n_positive, n_pairs))]
))
threshold_plot[, panel_label := factor(
  panel_text, c("GSE266330", "OEP: normal bone", "OEP: primary tumour")
)]
threshold_plot[, lineage_label := factor(lineage_labels[lineage], levels = rev(lineage_labels))]
threshold_plot[, threshold_label := factor(min_cells, levels = sort(unique(min_cells)))]
threshold_plot[, tile_label := sprintf("%.4f\n%s", effect_value, unit_label)]
sensitivity_limit <- ceiling(max(abs(threshold_plot$effect_value), na.rm = TRUE) * 100) / 100
p5c <- ggplot(threshold_plot, aes(threshold_label, lineage_label, fill = effect_value)) +
  geom_tile(colour = "white", linewidth = .42) +
  geom_text(aes(label = tile_label), size = 1.82, lineheight = .90) +
  facet_grid(. ~ panel_label, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = blue, mid = "white", high = red, midpoint = 0,
                       limits = c(-sensitivity_limit, sensitivity_limit), oob = squish,
                       name = "Programme-score\ndifference") +
  labs(title = "Cell-threshold sensitivity", x = "Minimum cells", y = NULL) + theme_pub +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "right",
        panel.spacing.x = unit(1.8, "mm"), legend.key.height = unit(13, "mm")) +
  guides(fill = guide_colourbar(title.position = "top", barheight = unit(15, "mm")))

axis_recon <- merge(
  axis_recon,
  axis_metrics[, .(model, spearman_rho, root_mean_squared_error, n_samples)],
  by = "model", all.x = TRUE, sort = FALSE
)
model_short <- c(`Broad block` = "Broad block", `Myeloid block` = "Myeloid block",
                 `T/NK block` = "T/NK block", `Broad block + depth` = "Broad + depth")
axis_recon[, model_label := factor(model_short[model_label], levels = unname(model_short))]
axis_recon[, facet_label := as.character(model_label)]
axis_recon[, facet_label := factor(facet_label, levels = unique(facet_label[order(model_order)]))]
recon_limit <- max(abs(c(axis_recon$observed_Axis1, axis_recon$predicted_Axis1)), na.rm = TRUE)
recon_annotations <- unique(axis_recon[, .(facet_label, spearman_rho, root_mean_squared_error)])
recon_annotations[, `:=`(
  x = -recon_limit * .92, y = recon_limit * .92,
  label = sprintf("rho %.2f\nRMSE %.2f", spearman_rho, root_mean_squared_error)
)]
p5d <- ggplot(axis_recon, aes(observed_Axis1, predicted_Axis1, colour = cancer)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = .31, colour = "#8B8F94") +
  geom_point(size = .82, alpha = .78, stroke = 0) +
  geom_text(data = recon_annotations, aes(x = x, y = y, label = label), inherit.aes = FALSE,
            hjust = 0, vjust = 1, size = 1.52, colour = dark, lineheight = .88) +
  facet_wrap(~facet_label, nrow = 2) +
  scale_colour_manual(values = c(prostate = blue, renal = orange), name = "Discovery cancer") +
  coord_equal(xlim = c(-recon_limit, recon_limit), ylim = c(-recon_limit, recon_limit)) +
  labs(title = "Held-out Axis1 reconstruction",
       x = "Observed Axis1", y = "Reconstructed Axis1") + theme_pub +
  theme(legend.position = "bottom", legend.title = element_blank(),
        strip.text = element_text(size = 6.15, lineheight = .88),
        panel.spacing = unit(1.4, "mm"), axis.text = element_text(size = 6.3))

axis_external <- axis_external[projectable == TRUE]
origin_labels <- c(BC = "Breast", BDC = "Bile duct", CC = "Colon", EC = "Oesophageal",
                   KC = "Kidney", LC = "Lung", PC = "Prostate", TC = "Thyroid")
if (axis_model$selected_k != 4L || ncol(axis_model$x_raw_clr) != 26L ||
    !identical(names(axis_model$blocks), colnames(axis_model$x_raw_clr))) {
  stop("Frozen Gate12G model contract is invalid")
}
axis_features <- colnames(axis_model$x_raw_clr)
broad_levels <- sub("^Broad__", "", axis_features[axis_model$blocks == "Broad"])
myeloid_levels <- sub("^Myeloid__", "", axis_features[axis_model$blocks == "Myeloid"])
tnk_levels <- sub("^T_NK__", "", axis_features[axis_model$blocks == "T_NK"])
axis_base <- unique(axis_external[, .(patient_id, cancer, analysis_role, condition, projectable)])
make_axis_count_wide <- function(dt, label_col, levels) {
  observed <- dt[!is.na(get(label_col)), .N, by = .(patient_id, label = get(label_col))]
  grid <- axis_base[, .(label = levels), by = .(patient_id, cancer, analysis_role, condition, projectable)]
  z <- merge(grid, observed, by = c("patient_id", "label"), all.x = TRUE)
  z[is.na(N), N := 0L]
  wide <- dcast(z, patient_id + cancer + analysis_role + condition + projectable ~ label,
                value.var = "N", fill = 0)
  for (lev in levels) if (!lev %chin% names(wide)) wide[, (lev) := 0L]
  setcolorder(wide, c("patient_id", "cancer", "analysis_role", "condition", "projectable", levels))
  setkey(wide, patient_id)
  wide
}
axis_broad_w <- make_axis_count_wide(axis_cells, "gate12z_broad", broad_levels)
axis_myeloid_w <- make_axis_count_wide(axis_cells[gate12z_broad == "Myeloid"], "gate12z_state", myeloid_levels)
axis_tnk_w <- make_axis_count_wide(axis_cells[gate12z_broad == "T_NK"], "gate12z_state", tnk_levels)
if (!identical(axis_broad_w$patient_id, axis_myeloid_w$patient_id) ||
    !identical(axis_broad_w$patient_id, axis_tnk_w$patient_id)) stop("External feature blocks are misaligned")
close_clr <- function(m) {
  p <- as.matrix(m) + .5
  p <- p / rowSums(p)
  log(p) - rowMeans(log(p))
}
axis_xb <- close_clr(axis_broad_w[, ..broad_levels])
axis_xm <- close_clr(axis_myeloid_w[, ..myeloid_levels])
axis_xt <- close_clr(axis_tnk_w[, ..tnk_levels])
colnames(axis_xb) <- paste0("Broad__", broad_levels)
colnames(axis_xm) <- paste0("Myeloid__", myeloid_levels)
colnames(axis_xt) <- paste0("T_NK__", tnk_levels)
axis_x_transformed <- cbind(axis_xb, axis_xm, axis_xt)[, axis_features, drop = FALSE]
rownames(axis_x_transformed) <- axis_broad_w$patient_id
axis_block_scales <- vapply(unique(axis_model$blocks), function(block) {
  j <- which(axis_model$blocks == block)
  sqrt(sum(apply(axis_model$x_raw_clr[, j, drop = FALSE], 2L, var)))
}, numeric(1))
for (block in unique(axis_model$blocks)) {
  j <- which(axis_model$blocks == block)
  axis_x_transformed[, j] <- axis_x_transformed[, j, drop = FALSE] / axis_block_scales[[block]]
}
axis_check <- sweep(axis_x_transformed, 2L, axis_model$pca$center, "-") %*%
  axis_model$pca$rotation[, seq_len(axis_model$selected_k), drop = FALSE]
saved_axes <- as.matrix(axis_external[match(rownames(axis_check), patient_id), .(Axis1, Axis2, Axis3, Axis4)])
if (max(abs(axis_check - saved_axes)) >= 1e-10) stop("External frozen-feature reconstruction does not reproduce saved axes")
axis_external[, origin_text := fifelse(condition == "healthy_marrow", "Healthy marrow",
                                       unname(origin_labels[cancer]))]
axis_external[, condition_text := fifelse(condition == "healthy_marrow", "Healthy", "Bone metastasis")]
axis_external[, origin_rank := match(origin_text, c("Healthy marrow", "Breast", "Bile duct", "Colon",
                                                   "Oesophageal", "Kidney", "Lung", "Prostate", "Thyroid"))]
setorder(axis_external, origin_rank, Axis1, patient_id)
axis_external[, sample_order := seq_len(.N)]
axis_external[, patient_factor := factor(patient_id, levels = patient_id)]
axis_display_matrix <- t(scale(axis_x_transformed[axis_external$patient_id, , drop = FALSE],
                               center = TRUE, scale = TRUE))
axis_feature_heat <- as.data.table(as.table(axis_display_matrix))
setnames(axis_feature_heat, c("feature", "patient_id", "display_z"))
axis_feature_meta <- data.table(
  feature = axis_features,
  block = factor(unname(axis_model$blocks), c("Broad", "Myeloid", "T_NK")),
  feature_order = seq_along(axis_features),
  feature_label = gsub("_", " ", sub("^(Broad|Myeloid|T_NK)__", "", axis_features))
)
axis_feature_heat <- merge(axis_feature_heat, axis_feature_meta, by = "feature", all.x = TRUE, sort = FALSE)
axis_feature_heat <- merge(axis_feature_heat,
                           axis_external[, .(patient_id, patient_factor, sample_order, origin_text,
                                             condition_text, Axis1)],
                           by = "patient_id", all.x = TRUE, sort = FALSE)
axis_feature_heat[, feature_label := factor(feature_label, levels = rev(axis_feature_meta$feature_label))]
origin_palette <- c(`Healthy marrow` = "#9AA0A6", Breast = "#D55E00", `Bile duct` = "#E69F00",
                    Colon = "#56B4E9", Oesophageal = "#CC79A7", Kidney = "#0072B2",
                    Lung = "#009E73", Prostate = "#6A3D9A", Thyroid = "#A6761D")
axis_annotation <- rbindlist(list(
  axis_external[, .(patient_factor, sample_order, annotation = "Condition",
                    annotation_colour = fifelse(condition_text == "Healthy", "#A8ADB3", red))],
  axis_external[, .(patient_factor, sample_order, annotation = "Origin",
                    annotation_colour = unname(origin_palette[origin_text]))]
))
axis_annotation[, annotation := factor(annotation, c("Origin", "Condition"))]
origin_bounds <- axis_external[, .(xmin = min(sample_order) - .5, xmax = max(sample_order) + .5,
                                   midpoint = mean(range(sample_order)), n = .N), by = origin_text]
origin_bounds[, short_label := c(`Healthy marrow` = "Healthy", Breast = "BC", `Bile duct` = "BDC",
                                 Colon = "CC", Oesophageal = "EC", Kidney = "KC", Lung = "LC",
                                 Prostate = "PC", Thyroid = "TC")[origin_text]]
p5e_annotation <- ggplot(axis_annotation, aes(sample_order, annotation, fill = annotation_colour)) +
  geom_tile(colour = "white", linewidth = .12) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(.5, nrow(axis_external) + .5), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) + theme_void(base_family = "Arial", base_size = 6.5) +
  theme(axis.text.y = element_text(colour = dark, size = 5.7, hjust = 1),
        axis.text.x = element_blank(), plot.margin = margin(0, 3, 0, 2))
p5e_heat <- ggplot(axis_feature_heat, aes(sample_order, feature_label, fill = display_z)) +
  geom_tile() +
  geom_vline(data = origin_bounds[-.N], aes(xintercept = xmax), inherit.aes = FALSE,
             colour = "white", linewidth = .48) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-2.5, 2.5), oob = squish,
                       name = "Row z-score") +
  scale_x_continuous(limits = c(.5, nrow(axis_external) + .5),
                     breaks = origin_bounds$midpoint, labels = origin_bounds$short_label,
                     expand = c(0, 0)) +
  labs(x = NULL, y = NULL) + theme_pub +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(size = 5.5, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 5.7), legend.position = "right",
        legend.key.height = unit(12, "mm"), plot.margin = margin(0, 3, 3, 2)) +
  guides(fill = guide_colourbar(title.position = "top", barheight = unit(13, "mm")))
p5e <- wrap_elements(full =
  p5e_heat +
    plot_annotation(title = "External composition landscape") &
    theme(plot.title = element_text(family = "Arial", size = 9.4, face = "bold", hjust = 0,
                                    margin = margin(l = 12, b = 2)),
          plot.margin = margin(4, 4, 3, 4))
)

axis_paired <- axis_paired[contrast %in% c("bm_vs_normal_bone", "bm_vs_primary")]
axis_paired[, comparison_label := factor(
  contrast, c("bm_vs_normal_bone", "bm_vs_primary"),
  c("Normal bone", "Primary tumour")
)]
axis_pair_effect <- copy(axis_paired)
axis_pair_effect[, direction := factor(fifelse(difference >= 0, "BM higher", "Control higher"),
                                       c("Control higher", "BM higher"))]
axis_pair_effect[, patient_key := paste(contrast, patient_id, sep = "::")]
patient_key_order <- axis_pair_effect[order(contrast, difference), patient_key]
axis_pair_effect[, patient_factor := factor(patient_key, levels = rev(patient_key_order))]
effect_limit <- max(abs(axis_pair_effect$difference)) * 1.16
p5f <- ggplot(axis_pair_effect, aes(difference, patient_factor, fill = direction)) +
  geom_vline(xintercept = 0, linewidth = .35, colour = dark) +
  geom_col(width = .60, colour = "white", linewidth = .25) +
  geom_text(aes(label = sprintf("%+.2f", difference),
                hjust = fifelse(difference >= 0, -.12, 1.12)),
            size = 1.75, colour = dark) +
  facet_wrap(~comparison_label, nrow = 1, scales = "free_y",
             labeller = as_labeller(c(
               `Normal bone` = "Normal bone\n1/2 positive; required 2/2",
               `Primary tumour` = "Primary tumour\n2/3 positive; required 2/3"
             ))) +
  scale_fill_manual(values = c(`Control higher` = blue, `BM higher` = red), name = NULL) +
  scale_x_continuous(limits = c(-effect_limit, effect_limit),
                     breaks = pretty_breaks(4), expand = expansion(mult = c(.02, .02))) +
  scale_y_discrete(labels = function(x) sub("^.*::", "", x)) +
  labs(title = "Paired Axis1 effects",
       x = expression(Delta*" frozen Axis1"), y = NULL) + theme_pub +
  theme(legend.position = "bottom", strip.text = element_text(size = 6.8),
        axis.text.y = element_text(size = 6.3), panel.spacing.x = unit(3.0, "mm"))

figure5 <- (p5a | p5b) / p5c / (p5d | p5e) / p5f +
  plot_layout(heights = c(1.00, .79, 1.32, .84), widths = c(.82, 1.18)) +
  plot_annotation(tag_levels = "A", theme = theme(plot.margin = margin(4, 7, 5, 7)))
save_figure(figure5, "Figure5_review_driven_representation_stress_test", 10.6)

## Figure 6 -----------------------------------------------------------------
overlay <- read_rel("results/gate12ai_submission_ready/source_data/Figure6A_B_spatial_scores.tsv.gz")
effects <- read_rel("results/gate12ai_submission_ready/source_data/Figure6D_section_effects.tsv")
assets <- read_rel(paste0("results/gate12ad_figure_restructure/phase_a_source_provenance/",
                         "provenance/Figure6_histology_asset_index.tsv"))
scale_bars <- read_rel(paste0("results/gate12ad_figure_restructure/phase_a_source_provenance/",
                             "provenance/Figure6_scale_bar_audit.tsv"))
scale_audit <- read_rel(paste0("results/gate12ad_figure_restructure/phase_a_source_provenance/",
                              "provenance/Figure6_shared_color_scale_audit.tsv"))
shared_limit <- unique(scale_audit$shared_abs_q98_limit)
if (length(shared_limit) != 1L || !is.finite(shared_limit)) stop("Invalid shared spatial scale")
samples <- assets$sample
theme_spatial <- theme_void(base_size = 8, base_family = "Arial") +
  theme(plot.title = element_text(face = "bold", size = 7.6, hjust = .5, margin = margin(b = 1)),
        legend.title = element_text(size = 7.1), legend.text = element_text(size = 6.8),
        legend.key.width = unit(13, "mm"), legend.key.height = unit(2.6, "mm"),
        plot.margin = margin(1, 1, 1, 1))
lighten_histology <- function(image, amount = .12) image * (1 - amount) + amount
spatial_assets <- list()
for (i in seq_len(nrow(assets))) {
  sid <- assets$sample[[i]]
  image <- readPNG(file.path(root, assets$image_relative_path[[i]]))
  d <- overlay[sample == sid]
  xr <- range(d$plot_x_lowres); yr <- range(d$plot_y_lowres)
  px <- diff(xr) * .035; py <- diff(yr) * .035
  spatial_assets[[sid]] <- list(image = image, light_image = lighten_histology(image),
    width = dim(image)[2], height = dim(image)[1],
    crop = c(xmin = max(0, xr[1] - px), xmax = min(dim(image)[2], xr[2] + px),
             ymin = max(0, yr[1] - py), ymax = min(dim(image)[1], yr[2] + py)),
    bar = scale_bars[sample == sid, pixels_per_500um_lowres])
}
spatial_panel <- function(sid, layer = c("histology", "full", "excluded"), show_column_title = FALSE) {
  layer <- match.arg(layer); asset <- spatial_assets[[sid]]; d <- copy(overlay[sample == sid]); crop <- asset$crop
  dx <- crop[["xmax"]] - crop[["xmin"]]; dy <- crop[["ymax"]] - crop[["ymin"]]
  x0 <- crop[["xmin"]] + .055 * dx; x1 <- x0 + asset$bar; y0 <- crop[["ymin"]] + .045 * dy
  p <- ggplot() + annotation_raster(if (layer == "histology") asset$image else asset$light_image,
                                    xmin = 0, xmax = asset$width, ymin = 0, ymax = asset$height)
  if (layer == "histology") {
    p <- p + geom_point(data = d[source_annotated_tumor == TRUE], aes(plot_x_lowres, plot_y_lowres),
                        shape = 21, fill = NA, colour = "#202020", stroke = .16, size = .78)
  } else {
    d[, score := if (layer == "full") full_axis1_display else malignant_excluded_axis1_display]
    p <- p + geom_point(data = d, aes(plot_x_lowres, plot_y_lowres, colour = score),
                        shape = 16, size = .55, alpha = .70) +
      scale_colour_gradient2(low = blue, mid = "#F7F7F7", high = red, midpoint = 0,
                             limits = c(-shared_limit, shared_limit), oob = squish,
                             name = "Frozen Axis1", guide = guide_colourbar(title.position = "top"))
  }
  p + annotate("segment", x = x0, xend = x1, y = y0, yend = y0, linewidth = .60, colour = "black") +
    annotate("text", x = x0, y = y0 + .028 * dy, label = "500~mu*m", parse = TRUE,
             hjust = 0, size = 2.25, fontface = "bold") +
    coord_fixed(xlim = crop[c("xmin", "xmax")], ylim = crop[c("ymin", "ymax")], expand = FALSE) +
    labs(title = if (show_column_title) paste0("S", match(sid, samples)) else NULL) + theme_spatial
}
hist_row <- wrap_plots(lapply(samples, spatial_panel, layer = "histology", show_column_title = TRUE), nrow = 1)
full_row <- wrap_plots(lapply(samples, spatial_panel, layer = "full"), nrow = 1, guides = "collect") & theme(legend.position = "none")
excl_row <- wrap_plots(lapply(samples, spatial_panel, layer = "excluded"), nrow = 1, guides = "collect") & theme(legend.position = "none")
row_header <- function(plot, title) wrap_elements(full = plot + plot_annotation(
  title = title,
  theme = theme(plot.title = element_text(family = "Arial", size = 9.1, face = "bold", margin = margin(b = 2)),
                plot.margin = margin(1, 2, 1, 2))))
p6a <- row_header(hist_row, "H&E and tumour-labelled spots")
p6b_maps <- row_header(full_row, "Full Axis1")
p6c_maps <- row_header(excl_row, "Malignant-excluded Axis1")
scale_data <- data.table(score = seq(-shared_limit, shared_limit, length.out = 240L), y = 1)
shared_scale_plot <- ggplot(scale_data, aes(score, y, fill = score)) +
  geom_tile() +
  scale_fill_gradient2(low = blue, mid = "#F7F7F7", high = red, midpoint = 0,
                       limits = c(-shared_limit, shared_limit), guide = "none") +
  scale_x_continuous(breaks = c(-.2, 0, .2)) +
  labs(x = "Axis1", y = NULL) +
  theme_void(base_family = "Arial") +
  theme(axis.title.x = element_text(size = 7.2, margin = margin(t = 2)),
        axis.text.x = element_text(size = 6.8, colour = dark),
        axis.ticks.x = element_line(colour = dark), axis.ticks.length = unit(1.5, "mm"),
        plot.margin = margin(0, 125, 0, 125))
p6c <- wrap_elements(full = p6c_maps / shared_scale_plot + plot_layout(heights = c(1, .11)))
p6b <- p6b_maps

block_scores <- overlay[distance_ring != "unreachable", .(block_score = median(full_axis1, na.rm = TRUE), spots = .N),
                        by = .(sample, spatial_block, distance_ring)]
block_scores[, region := fcase(distance_ring == "tumor", "Tumour", distance_ring == "boundary", "Boundary",
  distance_ring %in% c("proximal_2_3", "intermediate_4_6"), "Intermediate",
  distance_ring == "distal_7plus", "Distal", default = NA_character_)]
block_scores <- block_scores[!is.na(region)]
block_scores[, region := factor(region, c("Tumour", "Boundary", "Intermediate", "Distal"))]
block_scores[, section := factor(sample, samples, paste0("S", seq_along(samples)))]
p6d <- ggplot(block_scores, aes(region, block_score, fill = region)) +
  geom_boxplot(width = .62, linewidth = .36, outlier.shape = NA, alpha = .32) +
  geom_point(position = position_jitter(width = .10), size = .78, alpha = .70, stroke = 0) +
  facet_wrap(~section, nrow = 2) +
  scale_fill_manual(values = c(Tumour = red, Boundary = orange, Intermediate = "#E6B43C", Distal = blue), guide = "none") +
  labs(title = "Full Axis1 by tumour distance",
       x = NULL, y = "Spatial-block median") + theme_pub +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), strip.text = element_text(size = 7.1))
effects[, sample := factor(sample, levels = samples, labels = paste0("S", seq_along(samples)))]
effects[, layer_label := factor(layer_label, c("Full Axis1", "Malignant-excluded Axis1", "One-ring RCTD proxy"))]
effects[, layer_short := factor(layer_label,
  levels = c("Full Axis1", "Malignant-excluded Axis1", "One-ring RCTD proxy"),
  labels = c("Full", "Malignant-excluded", "One-ring RCTD"))]
p6e <- ggplot(effects, aes(estimate, sample)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = .32, colour = "#777777") +
  geom_errorbar(data = effects[decision_evaluable == TRUE], aes(xmin = ci_low, xmax = ci_high),
                orientation = "y", width = .13, linewidth = .48, colour = "#356A9A") +
  geom_point(data = effects[decision_evaluable == TRUE], size = 1.72, colour = "#356A9A") +
  geom_text(data = effects[decision_evaluable == FALSE], aes(x = .13, label = "NE"),
            colour = "#777777", size = 2.25) +
  facet_grid(cols = vars(layer_short)) +
  scale_y_discrete(limits = rev(paste0("S", seq_along(samples))), drop = FALSE) +
  scale_x_continuous(limits = c(-.03, .26), breaks = c(0, .1, .2), expand = expansion(mult = c(.01, .01))) +
  labs(title = "Boundary - distal effects",
       x = "Score difference", y = NULL) + theme_pub +
  theme(strip.background = element_blank(), strip.text = element_text(size = 6.6, face = "bold"),
        panel.spacing.x = unit(2.0, "mm"), axis.text.x = element_text(size = 6.5),
        axis.line.x = element_line(colour = dark, linewidth = .3), axis.line.y = element_blank(),
        axis.ticks.y = element_blank())

figure6 <- p6a / p6b / p6c / (p6d | p6e) +
  plot_layout(heights = c(1.10, 1.08, 1.50, 1.30), widths = c(1.03, .97)) +
  plot_annotation(tag_levels = "A", theme = theme(plot.margin = margin(4, 7, 5, 7)))
save_figure(figure6, "Figure6_review_driven_spatial_organization", 10.3)

## Source-data outputs and figure contract ----------------------------------
write_source(top_ligands, "Figure4A_nichenet_predicted_ligand_activity.tsv")
write_source(short_chat, "Figure4B_inferred_cellchat_pairs.tsv")
write_source(candidate_targets, "Figure4C_nichenet_regulatory_potential.tsv")
write_source(ligand_effects, "Figure4D_patient_pseudobulk_ligand_effects.tsv")
write_source(sensitivity_plot, "Figure4E_quantitative_threshold_sensitivity.tsv")
write_source(selection_stability, "Figure4E_selection_stability.tsv")
write_source(decision, "Figure4E_external_decision.tsv")
write_source(gse_common, "Figure5A_GSE266330_program_failure.tsv")
write_source(oep_pair_heat, "Figure5B_OEP_paired_programme_differences.tsv")
write_source(pair_annotation, "Figure5B_OEP_pair_annotations.tsv")
write_source(threshold_plot, "Figure5C_combined_programme_sensitivity.tsv")
write_source(axis_recon, "Figure5D_Axis1_four_model_reconstruction.tsv")
write_source(axis_metrics, "Figure5D_Axis1_reconstruction_metrics.tsv")
write_source(axis_feature_heat, "Figure5E_external_frozen_feature_heatmap.tsv.gz")
write_source(axis_external, "Figure5E_external_patient_annotations.tsv")
write_source(data.table(feature = rep(axis_features, each = nrow(axis_x_transformed)),
                        block = rep(unname(axis_model$blocks), each = nrow(axis_x_transformed)),
                        frozen_block_scaled_clr = as.vector(axis_x_transformed),
                        patient_id = rep(rownames(axis_x_transformed), times = ncol(axis_x_transformed))),
             "Figure5E_external_frozen_transformed_matrix.tsv.gz")
write_source(axis_pair_effect, "Figure5F_Axis1_paired_effects.tsv")
write_source(overlay, "Figure6A_C_spatial_feature_maps.tsv.gz")
write_source(block_scores, "Figure6D_block_region_scores.tsv")
write_source(effects, "Figure6E_section_effects.tsv")

contract <- data.table(
  figure = c(rep("Figure4", 5), rep("Figure5", 6), rep("Figure6", 5)),
  panel = c(LETTERS[1:5], LETTERS[1:6], LETTERS[1:5]),
  disposition = c("retain", "retain_polish", "retain_polish", "retain", "promoted_negative_boundary",
                  "enlarged_failure", "published_paired_difference_heatmap", "merged_sensitivity",
                  "published_four_model_parity", "published_annotated_patient_heatmap",
                  "published_zero_centered_patient_effect", "enlarged_map", "enlarged_map",
                  "enlarged_map", "retain_polish", "retain_enlarge"),
  refit_or_reanalysis = FALSE
)
fwrite(contract, file.path(admin_out, "GATE12BE_FIGURE_CONTRACT.tsv"), sep = "\t")
cat("GATE12BE_FIGURES456_STATUS=COMPLETE\n")
