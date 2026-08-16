#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(trailingOnly = TRUE)
project <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else normalizePath(".")
out <- file.path(project, "results", "gate12bd_aligned_submission_package")
figure_out <- file.path(out, "figures", "supplementary")
source_out <- file.path(out, "source_data", "supplementary")
admin_out <- file.path(out, "admin")
dir.create(figure_out, recursive = TRUE, showWarnings = FALSE)
dir.create(source_out, recursive = TRUE, showWarnings = FALSE)
dir.create(admin_out, recursive = TRUE, showWarnings = FALSE)
set.seed(20260812L)

read_rel <- function(path) fread(file.path(project, path))
blue <- "#3973B9"
orange <- "#E17C2B"
green <- "#009E73"
purple <- "#7047A3"
red <- "#C44E52"
grey <- "#BFC3C7"

theme_pub <- theme_classic(base_size = 8.2, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 9.2, face = "bold", margin = margin(b = 1)),
    plot.subtitle = element_text(size = 6.8, colour = "#55606E", lineheight = 0.95,
                                 margin = margin(b = 3)),
    axis.title = element_text(size = 7.4),
    axis.text = element_text(size = 6.6, colour = "#30343B"),
    strip.text = element_text(size = 6.8, face = "bold"),
    strip.background = element_rect(fill = "#F3F4F5", colour = NA),
    legend.title = element_text(size = 6.8),
    legend.text = element_text(size = 6.3),
    panel.grid = element_blank(),
    plot.tag = element_text(size = 11.0, face = "bold"),
    plot.margin = margin(4, 4, 4, 4)
  )

save_figure <- function(plot, filename, width = 7.2, height) {
  cairo_png <- function(filename, width, height, bg = "white", ...) {
    grDevices::png(filename, width = width, height = height, units = "in", res = 450,
                   type = "cairo", bg = bg, ...)
  }
  ggsave(file.path(figure_out, paste0(filename, ".png")), plot, width = width,
         height = height, units = "in", device = cairo_png, bg = "white", limitsize = FALSE)
  ggsave(file.path(figure_out, paste0(filename, ".pdf")), plot, width = width,
         height = height, units = "in", device = grDevices::cairo_pdf,
         bg = "white", limitsize = FALSE)
}

## Supplementary Figure S11: the diagnostic panels moved out of main Figure 5.
axis_source <- file.path(project, "results", "gate12al_figure1_integration", "source_data")
loadings <- fread(file.path(axis_source, "Figure5A_axis1_loadings.tsv"))
contributions <- fread(file.path(axis_source, "Figure5B_block_contributions.tsv"))
recon <- fread(file.path(axis_source, "Figure5C_held_out_reconstructions.tsv"))
recon_metrics <- fread(file.path(axis_source, "Figure5C_reconstruction_metrics.tsv"))
external_scores <- fread(file.path(axis_source, "Figure5D_unpaired_external_scores.tsv"))
paired <- fread(file.path(axis_source, "Figure5E_paired_endpoints.tsv"))

block_cols <- c(Broad = green, Myeloid = "#E69F00", `T/NK` = "#0072B2")
loadings[, feature_label := factor(feature_label, levels = feature_label[order(frozen_loading)])]
p11a <- ggplot(loadings, aes(frozen_loading, feature_label, colour = block_display)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.30, colour = "#777777") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y",
                width = 0.10, linewidth = 0.42) +
  geom_point(aes(fill = interval_excludes_zero), shape = 21, size = 1.65,
             colour = "#30343B", stroke = 0.28) +
  scale_colour_manual(values = block_cols, name = "Feature block") +
  scale_fill_manual(values = c(`TRUE` = "#30343B", `FALSE` = "white"), guide = "none") +
  labs(title = "Frozen Axis1 loadings", subtitle = "Patient-cluster bootstrap intervals; relative compositional features",
       x = "Loading (95% bootstrap interval)", y = NULL) +
  theme_pub + theme(legend.position = "bottom")

contributions[, compartment_label := factor(tolower(compartment),
                                             levels = c("distal", "involved", "tumor"),
                                             labels = c("D", "I", "T"))]
contributions[, cancer_label := factor(tolower(cancer), levels = c("prostate", "renal"),
                                       labels = c("Prostate", "Renal"))]
contributions[, block_label := factor(block_display, levels = c("Broad", "Myeloid", "T/NK"))]
p11b <- ggplot(contributions, aes(compartment_label, contribution, colour = cancer_label,
                                  shape = cancer_label)) +
  geom_hline(yintercept = 0, linewidth = 0.28, colour = "#BBBBBB") +
  geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.38),
             size = 0.75, alpha = 0.32, stroke = 0) +
  stat_summary(fun = median, geom = "point", position = position_dodge(width = 0.38), size = 1.9) +
  stat_summary(fun.data = function(x) data.frame(y = median(x), ymin = quantile(x, .25),
                                                 ymax = quantile(x, .75)),
               geom = "errorbar", position = position_dodge(width = 0.38),
               width = 0.12, linewidth = 0.42) +
  facet_wrap(~block_label, nrow = 1) +
  scale_colour_manual(values = c(Prostate = blue, Renal = orange), name = NULL) +
  scale_shape_manual(values = c(Prostate = 16, Renal = 17), name = NULL) +
  labs(title = "Exact block decomposition", subtitle = "The three blocks decompose the same coordinate",
       x = "Anatomical compartment", y = "Signed contribution") +
  theme_pub + theme(legend.position = "bottom")

recon_metrics[, model_label := factor(model_label, levels = rev(model_label[order(model_order)]))]
p11c_left <- ggplot(recon_metrics, aes(cross_validated_R2_training_mean, model_label)) +
  geom_vline(xintercept = 0.90, linewidth = 0.30, colour = "#BBBBBB") +
  geom_point(size = 2.2, colour = purple) +
  scale_x_continuous(limits = c(0.85, 0.965)) +
  labs(title = "Patient-held-out block reconstruction",
       subtitle = "Cross-validated against the training-mean baseline",
       x = expression("Cross-validated reconstruction " * R^2), y = NULL) + theme_pub
recon_one <- recon[model == "myeloid_only"]
p11c_right <- ggplot(recon_one, aes(predicted_Axis1, observed_Axis1, colour = cancer)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.30, colour = "#888888") +
  geom_point(size = 1.0, alpha = 0.78) +
  scale_colour_manual(values = c(prostate = blue, renal = orange), guide = "none") +
  coord_equal() +
  labs(title = "Observed versus reconstructed Axis1",
       subtitle = "Myeloid-only model; 41 held-out samples",
       x = "Reconstructed Axis1", y = "Observed Axis1") + theme_pub
p11c <- p11c_left | p11c_right + plot_layout(widths = c(1.05, 0.95))

external_scores <- external_scores[projectable == TRUE]
external_scores[, condition_label := factor(display_group,
                                             levels = c("Healthy marrow", "Bone metastasis"))]
p11d <- ggplot(external_scores, aes(condition_label, Axis1, fill = condition_label)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.32, linewidth = 0.35) +
  geom_point(position = position_jitter(width = 0.10), size = 0.72, alpha = 0.65, stroke = 0) +
  scale_fill_manual(values = c(`Healthy marrow` = grey, `Bone metastasis` = red), guide = "none") +
  labs(title = "Unpaired GSE266330 directional support",
       subtitle = "Source-confounded; this comparison cannot replace matched normal bone",
       x = NULL, y = "Frozen Axis1") + theme_pub +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))

paired_primary <- paired[contrast %in% c("bm_vs_normal_bone", "bm_vs_primary")]
paired_primary[, comparison_label := factor(contrast,
  levels = c("bm_vs_normal_bone", "bm_vs_primary"),
  labels = c("Matched normal bone", "Matched primary tumour"))]
paired_primary[, patient_label := factor(patient_id, levels = unique(patient_id))]
p11e <- ggplot(paired_primary, aes(difference, patient_label, colour = cancer_code)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.32, colour = "#777777") +
  geom_point(size = 1.8) +
  facet_wrap(~comparison_label, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = c(LUCA = green, PRAD = purple), name = "Cancer") +
  labs(title = "Independent paired Axis1 endpoint",
       subtitle = "Metastasis minus control; patients are the biological units",
       x = "Paired Axis1 difference", y = NULL) +
  theme_pub + theme(legend.position = "bottom")

figure11 <- (p11a | p11b) +
  plot_layout(widths = c(1.03, 0.97)) +
  plot_annotation(
    title = "Supplementary Figure S11 | Axis1 construction and block redundancy",
    subtitle = "Frozen loadings and exact patient-sample block contributions",
    tag_levels = "A",
    theme = theme(plot.title = element_text(family = "Arial", size = 12.0, face = "bold"),
                  plot.subtitle = element_text(family = "Arial", size = 7.6, colour = "#4B5563"),
                  plot.tag = element_text(family = "Arial", size = 11.0, face = "bold"))
  )
save_figure(figure11, "FigureS11_axis1_construction_redundancy_endpoint", height = 5.2)

p11c_left <- p11c_left + labs(tag = "C")
p11c_right <- p11c_right + labs(tag = "D")
p11d <- p11d + labs(tag = "E")
p11e <- p11e + labs(tag = "F")
figure11_continued <- (p11c_left | p11c_right) / (p11d | p11e) +
  plot_layout(heights = c(1.0, 1.0), widths = c(1.02, 0.98)) +
  plot_annotation(
    title = "Supplementary Figure S11 continued | Reconstruction and endpoint boundary",
    subtitle = "Technical reconstruction, source-confounded projection and matched human comparisons",
    theme = theme(plot.title = element_text(family = "Arial", size = 12.0, face = "bold"),
                  plot.subtitle = element_text(family = "Arial", size = 7.6, colour = "#4B5563"),
                  plot.tag = element_text(family = "Arial", size = 11.0, face = "bold"))
  )
save_figure(figure11_continued, "FigureS11_axis1_construction_redundancy_endpoint_continued", height = 7.6)

fwrite(loadings, file.path(source_out, "FigureS11A_axis1_loadings.tsv"), sep = "\t")
fwrite(contributions, file.path(source_out, "FigureS11B_block_contributions.tsv"), sep = "\t")
fwrite(recon_metrics, file.path(source_out, "FigureS11C_reconstruction_metrics.tsv"), sep = "\t")
fwrite(recon, file.path(source_out, "FigureS11D_held_out_reconstructions.tsv"), sep = "\t")
fwrite(external_scores, file.path(source_out, "FigureS11E_unpaired_external_scores.tsv"), sep = "\t")
fwrite(paired_primary, file.path(source_out, "FigureS11F_paired_axis1_endpoints.tsv"), sep = "\t")

## Supplementary Figure S12: independent co-detection and unique-target sensitivity.
external <- read_rel("results/gate12c_external_validation/external_dataset_axis_summary.tsv")
patients <- read_rel("results/gate12c_external_validation/external_patient_axis_support.tsv")
sensitivity <- read_rel("results/gate12c_external_validation/external_threshold_sensitivity.tsv")
decision <- read_rel("results/gate12c_external_validation/gate12c_external_axis_decision.tsv")

core <- external[dataset %in% c("GSE266330", "OEP005136")]
core[, dataset := factor(dataset, levels = c("GSE266330", "OEP005136"))]
axis_levels <- c("SPP1 -> CD44", "CXCL16 -> CXCR6", "CCL4 -> CCR5")
core[, axis := factor(axis, levels = axis_levels)]
p12a <- ggplot(core, aes(support_fraction, axis, fill = dataset)) +
  geom_vline(xintercept = 0.50, linetype = 2, linewidth = 0.32, colour = "#777777") +
  geom_point(aes(size = eligible_patients), shape = 21, colour = "#30343B", stroke = 0.25,
             position = position_dodge(width = 0.38)) +
  scale_fill_manual(values = c(GSE266330 = blue, OEP005136 = orange), name = NULL) +
  scale_size_continuous(range = c(2.4, 5.0), name = "Eligible patients") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0.45, 1.01)) +
  labs(title = "Independent patient-level co-detection",
       subtitle = "Blue, GSE266330; orange, OEP005136; expression plausibility only",
       x = "Supporting eligible patients", y = NULL) + theme_pub +
  theme(legend.position = "none")

origin <- patients[dataset %in% c("GSE266330", "OEP005136") & condition == "bone_metastasis" &
                     eligible == TRUE,
                   .(support_fraction = mean(patient_support), eligible_patients = .N),
                   by = .(axis, dataset, cancer)]
origin[, axis := factor(axis, levels = axis_levels)]
origin[, cancer := factor(cancer, levels = sort(unique(cancer)))]
p12b <- ggplot(origin, aes(cancer, axis, fill = support_fraction)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = eligible_patients), size = 2.0) +
  facet_wrap(~dataset, ncol = 1, scales = "free_x") +
  scale_fill_gradient(low = "white", high = "#3973B9", limits = c(0, 1),
                      name = "Support fraction", labels = scales::percent) +
  labs(title = "Cross-origin expression support", subtitle = "Numbers are eligible patients",
       x = "Source-defined cancer origin", y = NULL) +
  theme_pub + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                    axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")

winner <- unique(sensitivity[, .(min_cells, min_detect, scenario_winner)])
winner[, min_detect_label := factor(paste0(round(100 * min_detect), "%"),
                                    levels = c("10%", "5%", "2%"))]
p12c <- ggplot(winner, aes(factor(min_cells), min_detect_label, fill = scenario_winner)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = sub(" -> ", "\n", scenario_winner, fixed = TRUE)), size = 2.25) +
  scale_fill_manual(values = c("CCL4 -> CCR5" = blue, "CXCL16 -> CXCR6" = "#B2182B",
                               "SPP1 -> CD44" = "#E69F00"), name = "Scenario winner") +
  labs(title = "Unique-target sensitivity",
       subtitle = "CXCL16-CXCR6 retained the primary rank in only 2/9 scenarios",
       x = "Minimum cells in sender and receiver state", y = "Detection threshold") +
  theme_pub + theme(legend.position = "none")

decision_plot <- decision[, .(axis, supporting_dataset_count, min_core_fraction,
                              final_freeze_eligible)]
decision_plot[, axis := factor(axis, levels = axis_levels)]
p12d <- ggplot(decision_plot, aes(min_core_fraction, axis)) +
  geom_vline(xintercept = 0.50, linetype = 2, linewidth = 0.32, colour = "#777777") +
  geom_point(aes(size = supporting_dataset_count, fill = final_freeze_eligible),
             shape = 21, colour = "#30343B", stroke = 0.30) +
  scale_fill_manual(values = c(`TRUE` = red, `FALSE` = "white"),
                    labels = c(`TRUE` = "Primary-rule selected before sensitivity",
                               `FALSE` = "External support"), name = NULL) +
  scale_size_continuous(range = c(2.5, 4.5), name = "Supporting datasets") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0.45, 1.01)) +
  labs(title = "Prespecified external decision",
       subtitle = "Red fill marks the pre-sensitivity primary-rule selection",
       x = "Minimum core-cohort support", y = NULL) +
  theme_pub + theme(legend.position = "none")

figure12 <- (p12a | p12b) / (p12c | p12d) +
  plot_layout(heights = c(1.02, 0.98), widths = c(0.92, 1.08)) +
  plot_annotation(
    title = "Supplementary Figure S12 | Communication plausibility and target-selection sensitivity",
    subtitle = "Co-detection supports a shortlist, but sensitivity prevents a unique target claim",
    tag_levels = "A",
    theme = theme(plot.title = element_text(family = "Arial", size = 12.0, face = "bold"),
                  plot.subtitle = element_text(family = "Arial", size = 7.6, colour = "#4B5563"),
                  plot.tag = element_text(family = "Arial", size = 11.0, face = "bold"))
  )
save_figure(figure12, "FigureS12_sender_receiver_external_sensitivity", height = 8.4)

fwrite(core, file.path(source_out, "FigureS12A_external_dataset_support.tsv"), sep = "\t")
fwrite(origin, file.path(source_out, "FigureS12B_external_origin_support.tsv"), sep = "\t")
fwrite(winner, file.path(source_out, "FigureS12C_target_sensitivity_winners.tsv"), sep = "\t")
fwrite(decision_plot, file.path(source_out, "FigureS12D_external_decision.tsv"), sep = "\t")

## The complete CellChat overview was deliberately removed from main Figure 4.
## Retain it as a descriptive supplementary archive rather than a mechanistic result.
communication_overview <- fread(file.path(
  project, "results", "gate12bd_aligned_submission_package", "source_data", "main",
  "Figure4C_tumor_communication_matrix.tsv"
))
sender_order <- communication_overview[, .(total = sum(probability)), by = source_label][
  order(total), source_label
]
receiver_order <- communication_overview[, .(total = sum(probability)), by = target_label][
  order(-total), target_label
]
communication_overview[, source_label := factor(source_label, levels = sender_order)]
communication_overview[, target_label := factor(target_label, levels = receiver_order)]
communication_overview[, cohort := factor(cohort, levels = c("Prostate", "Renal"))]
p12e <- ggplot(communication_overview, aes(target_label, source_label, fill = probability)) +
  geom_tile(colour = "white", linewidth = 0.42) +
  facet_wrap(~cohort, nrow = 1) +
  scale_fill_viridis_c(option = "C", trans = "sqrt", name = "Summed inferred\nprobability") +
  labs(
    tag = "E",
    title = "Descriptive tumour-compartment communication overview",
    subtitle = "Expression-derived CellChat screen; top sender and receiver states only",
    x = "Receiver state", y = "Sender state"
  ) +
  theme_pub +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 38, hjust = 1),
    legend.position = "right"
  )
figure12_continued <- p12e +
  plot_annotation(
    title = "Supplementary Figure S12 continued | Global descriptive communication screen",
    subtitle = "This overview is hypothesis generating and does not measure binding or signal transmission",
    theme = theme(
      plot.title = element_text(family = "Arial", size = 12.0, face = "bold"),
      plot.subtitle = element_text(family = "Arial", size = 7.6, colour = "#4B5563"),
      plot.tag = element_text(family = "Arial", size = 11.0, face = "bold")
    )
  )
save_figure(
  figure12_continued,
  "FigureS12_sender_receiver_external_sensitivity_continued",
  height = 5.7
)
fwrite(
  communication_overview,
  file.path(source_out, "FigureS12E_tumour_communication_overview.tsv"),
  sep = "\t"
)

audit <- data.table(
  figure = c(rep("FigureS11", 6), rep("FigureS12", 5)),
  panel = c(LETTERS[1:6], LETTERS[1:5]),
  visual_form = c("loading_forest", "faceted_contribution_plot", "reconstruction_metric_scatter",
                  "observed_reconstructed_scatter", "boxplot_with_patient_points", "paired_difference_dotplot",
                  "cohort_support_dotplot", "origin_support_heatmap",
                  "threshold_sensitivity_heatmap", "external_decision_dotplot",
                  "descriptive_sender_receiver_heatmap"),
  claim_level = c(rep("technical_or_endpoint_diagnostic", 6),
                  rep("hypothesis_generating_or_negative", 5)),
  main_text = FALSE
)
fwrite(audit, file.path(admin_out, "GATE12BD_SUPPLEMENTARY_VISUAL_AUDIT.tsv"), sep = "\t")
cat("GATE12BD_SUPPLEMENT_STATUS=COMPLETE\n")
