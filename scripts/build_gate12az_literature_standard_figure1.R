#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrastr)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(grid)
})

options(stringsAsFactors = FALSE, warn = 1)
set.seed(20260811L)

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) normalizePath(args[[1L]], mustWork = TRUE) else
  normalizePath(".", mustWork = TRUE)
output_gate <- if (length(args) >= 2L) args[[2L]] else
  "gate12az_literature_calibrated_redesign"
use_upgraded_umap <- output_gate %in% c(
  "gate12bb_umap_integrated_figures",
  "gate12be_review_driven_redesign"
)
out <- file.path(root, "results", output_gate)
figure_out <- file.path(out, "figures", "main")
source_out <- if (identical(output_gate, "gate12be_review_driven_redesign")) {
  file.path(out, "source_data", "main")
} else {
  file.path(out, "source_data")
}
admin_out <- file.path(out, "admin")
dir.create(figure_out, recursive = TRUE, showWarnings = FALSE)
dir.create(source_out, recursive = TRUE, showWarnings = FALSE)
dir.create(admin_out, recursive = TRUE, showWarnings = FALSE)

phase_a <- file.path(root, "results", "gate12ad_figure_restructure",
                     "phase_a_source_provenance", "source_data")
umap_file <- file.path(root, "results", "gate12an_umap_reembedding",
                       "K_compact_spca_coordinates.tsv.gz")
availability_file <- file.path(root, "results", "gate12ai_submission_ready",
                               "source_data", "Figure1B_patient_availability.tsv")
marker_file <- file.path(phase_a, "Figure1C_canonical_marker_dotplot.tsv")
composition_file <- file.path(phase_a, "Figure1D_all_sample_broad_composition.tsv")
required <- c(umap_file, availability_file, marker_file, composition_file)
if (use_upgraded_umap) {
  required <- c(
    required,
    file.path(root, "results", "gate12ba_umap_upgrade", "source_data",
              "Figure1B_global_UMAP_coordinates.tsv.gz"),
    file.path(root, "results", "gate12ba_umap_upgrade", "source_data",
              "Figure1B_global_label_anchors.tsv")
  )
}
if (any(!file.exists(required))) {
  stop("Missing Figure 1 input: ", paste(required[!file.exists(required)], collapse = ", "))
}

umap <- fread(umap_file)
availability <- fread(availability_file)
markers <- fread(marker_file)
composition <- fread(composition_file)
if (use_upgraded_umap) {
  upgraded_coords <- fread(
    file.path(root, "results", "gate12ba_umap_upgrade", "source_data",
              "Figure1B_global_UMAP_coordinates.tsv.gz"),
    select = c("cell_id", "plot_x", "plot_y")
  )
  upgraded_anchors <- fread(
    file.path(root, "results", "gate12ba_umap_upgrade", "source_data",
              "Figure1B_global_label_anchors.tsv")
  )
  umap <- merge(umap, upgraded_coords, by = "cell_id", all.x = TRUE, sort = FALSE)
  if (nrow(umap) != 107886L || any(!is.finite(umap$plot_x)) ||
      any(!is.finite(umap$plot_y))) stop("Upgraded global UMAP join failed")
}

broad_order <- c("T_NK", "Myeloid", "B", "Progenitor", "Stromal", "Endothelial",
                 "Osteoclast", "Osteoblast", "Malignant", "Erythroid", "Unassigned")
broad_labels <- c(
  T_NK = "T/NK", Myeloid = "Myeloid", B = "B cell", Progenitor = "Progenitor",
  Stromal = "Stromal", Endothelial = "Endothelial", Osteoclast = "Osteoclast",
  Osteoblast = "Osteoblast", Malignant = "Malignant", Erythroid = "Erythroid",
  Unassigned = "Unassigned"
)
# Okabe-Ito-derived and muted extensions; stable across every atlas panel.
broad_cols <- c(
  T_NK = "#3973B9", Myeloid = "#E68435", B = "#4FA3BE",
  Progenitor = "#C8A42A", Stromal = "#35976D", Endothelial = "#69B3A2",
  Osteoclast = "#C75B56", Osteoblast = "#8B79B9", Malignant = "#7047A3",
  Erythroid = "#D95B8A", Unassigned = "#B9BDC2"
)
cancer_cols <- c(Prostate = "#3973B9", Renal = "#E68435")
compartment_cols <- c(Distal = "#5AA9E6", Involved = "#E6B43C", Tumor = "#C44E52")

theme_pub <- theme_classic(base_size = 8.4, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 9.2, face = "bold", hjust = 0, margin = margin(b = 2)),
    plot.subtitle = element_text(size = 7.1, colour = "#525A61", hjust = 0,
                                 lineheight = 0.95, margin = margin(b = 3)),
    axis.title = element_text(size = 7.8),
    axis.text = element_text(size = 7.0, colour = "#30343B"),
    strip.text = element_text(size = 7.2, face = "bold"),
    strip.background = element_rect(fill = "#F3F4F5", colour = NA),
    legend.title = element_text(size = 7.0),
    legend.text = element_text(size = 6.7),
    legend.key.height = unit(2.6, "mm"),
    legend.key.width = unit(2.6, "mm"),
    panel.grid = element_blank(),
    plot.tag = element_text(size = 11.5, face = "bold"),
    plot.margin = margin(4, 4, 4, 4)
  )

## A. Cohort matrix: standard tile overview, no workflow arrows.
availability[, cancer := factor(cancer, levels = c("Prostate", "Renal"))]
availability[, compartment := factor(compartment, levels = c("Distal", "Involved", "Tumor"))]
availability[, patient_id := factor(patient_id, levels = rev(unique(availability[order(cancer, patient_id)]$patient_id)))]
p1a <- ggplot(availability, aes(compartment, patient_id, fill = present)) +
  geom_tile(width = 0.84, height = 0.84, colour = "white", linewidth = 0.25) +
  facet_grid(cancer ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c(`TRUE` = "#4C78A8", `FALSE` = "#E8EBED"),
                    breaks = c(TRUE, FALSE), labels = c("Available", "Not available"),
                    name = NULL) +
  labs(title = "Discovery cohort",
       subtitle = "18 patients | 42 bone specimens | three anatomical compartments",
       x = NULL, y = NULL) +
  theme_pub +
  theme(axis.ticks = element_blank(), axis.line = element_blank(),
        legend.position = "bottom", panel.spacing.y = unit(1.2, "mm"))

## B. Global UMAP: the standard atlas overview used in published single-cell studies.
umap[, broad_class := factor(broad_class, levels = broad_order)]
class_counts <- umap[, .N, by = broad_class][order(-N)]
draw_levels <- class_counts$broad_class
umap[, draw_group := match(broad_class, draw_levels)]
umap[, draw_random := runif(.N)]
setorder(umap, draw_group, draw_random)
if (use_upgraded_umap) {
  p1b <- ggplot(umap, aes(plot_x, plot_y, colour = broad_class, alpha = broad_class)) +
    ggrastr::geom_point_rast(size = 0.105, raster.dpi = 900, stroke = 0) +
    ggrepel::geom_text_repel(
      data = upgraded_anchors,
      aes(plot_x, plot_y, label = label), inherit.aes = FALSE,
      family = "Arial", fontface = "bold", size = 1.80,
      colour = "#202124", box.padding = 0.08,
      point.padding = 0.02, min.segment.length = Inf, segment.colour = NA,
      max.overlaps = Inf, force = 0.20, max.time = 2, max.iter = 10000,
      seed = 20260812L, show.legend = FALSE
    ) +
    scale_colour_manual(values = broad_cols, breaks = broad_order,
                        labels = broad_labels[broad_order], drop = FALSE, name = NULL) +
    scale_alpha_manual(values = c(setNames(rep(0.82, length(broad_order) - 1L),
                                           broad_order[-length(broad_order)]),
                                  Unassigned = 0.22), guide = "none") +
    scale_x_continuous(expand = expansion(mult = 0.012)) +
    scale_y_continuous(expand = expansion(mult = 0.012)) +
    coord_equal(clip = "off") +
    labs(title = "Integrated single-cell atlas",
         subtitle = sprintf("%s cells | fixed descriptive UMAP", comma(nrow(umap)))) +
    theme_void(base_size = 8.4, base_family = "Arial") +
    theme(plot.title = element_text(size = 9.2, face = "bold", hjust = 0),
          plot.subtitle = element_text(size = 7.1, colour = "#525A61", hjust = 0,
                                       margin = margin(b = 3)),
          legend.position = "none", plot.margin = margin(4, 4, 4, 4))
} else {
  p1b <- ggplot(umap, aes(umap_1, umap_2, colour = broad_class, alpha = broad_class)) +
    ggrastr::geom_point_rast(size = 0.028, raster.dpi = 800) +
    scale_colour_manual(values = broad_cols, breaks = broad_order,
                        labels = broad_labels[broad_order], drop = FALSE, name = NULL) +
    scale_alpha_manual(values = c(setNames(rep(0.66, length(broad_order) - 1L),
                                           broad_order[-length(broad_order)]),
                                  Unassigned = 0.16), guide = "none") +
    coord_equal() +
    labs(title = "Integrated single-cell atlas",
         subtitle = sprintf("%s cells | frozen integrated embedding", comma(nrow(umap)))) +
    theme_void(base_size = 8.4, base_family = "Arial") +
    theme(plot.title = element_text(size = 9.2, face = "bold", hjust = 0),
          plot.subtitle = element_text(size = 7.1, colour = "#525A61", hjust = 0,
                                       margin = margin(b = 3)),
          legend.position = "bottom", legend.justification = "left",
          legend.text = element_text(size = 6.5),
          legend.key.height = unit(2.4, "mm"), legend.key.width = unit(2.4, "mm"),
          plot.margin = margin(4, 4, 4, 4)) +
    guides(colour = guide_legend(ncol = 3, byrow = TRUE,
                                 override.aes = list(size = 2.4, alpha = 1)))
}

## C. Marker DotPlot: expression by colour, detection by point area.
markers[, broad_class_label := factor(
  broad_class_label,
  levels = rev(unique(markers[order(broad_class_order)]$broad_class_label))
)]
markers[, gene := factor(gene, levels = unique(markers[order(marker_gene_order)]$gene))]
p1c <- ggplot(markers, aes(gene, broad_class_label)) +
  geom_point(aes(size = detected_percent, fill = mean_scaled_expression),
             shape = 21, colour = "#4A4A4A", stroke = 0.18) +
  scale_fill_gradient2(low = "#3B6FB6", mid = "#F4F4F4", high = "#B94450",
                       midpoint = 0, limits = c(-2.5, 2.5), oob = squish,
                       name = "Scaled\nexpression") +
  scale_size_continuous(range = c(0.3, 3.4), limits = c(0, 100),
                        breaks = c(25, 50, 75), name = "Detected (%)") +
  labs(title = "Canonical lineage markers",
       subtitle = "Colour: gene-wise scaled mean | Size: detected-cell fraction",
       x = NULL, y = NULL) +
  theme_pub +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 50, hjust = 1), legend.position = "right")

## D. Per-sample stacked composition: a standard atlas composition panel.
composition[, broad_class := factor(broad_class, levels = broad_order)]
composition[, cancer_label := factor(cancer, levels = c("prostate", "renal"),
                                     labels = c("Prostate", "Renal"))]
composition[, compartment_label := factor(
  compartment, levels = c("distal", "involved", "tumor"),
  labels = c("Distal", "Involved", "Tumor")
)]
composition[, sample_axis := factor(sample_panel_order, levels = sort(unique(sample_panel_order)))]
sample_labels <- unique(composition[, .(sample_panel_order, patient_id)])[order(sample_panel_order)]
p1d <- ggplot(composition, aes(sample_axis, fraction, fill = broad_class)) +
  geom_col(width = 0.9, colour = "white", linewidth = 0.05) +
  facet_grid(~cancer_label + compartment_label, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = broad_cols, breaks = broad_order,
                    labels = broad_labels[broad_order], drop = FALSE, name = NULL) +
  scale_x_discrete(labels = setNames(sample_labels$patient_id, sample_labels$sample_panel_order)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(title = "Cellular composition across discovery samples",
       subtitle = "Every bar is one specimen; all broad classes are retained",
       x = NULL, y = "Cell fraction") +
  theme_pub +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom",
        panel.spacing.x = unit(0.7, "mm")) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

## E. Sample-level distributions: replaces all patient spaghetti lines.
major_classes <- c("T_NK", "Myeloid", "B", "Progenitor", "Stromal", "Malignant")
composition_box <- copy(composition[broad_class %in% major_classes])
composition_box[, broad_class_label := factor(
  broad_labels[as.character(broad_class)],
  levels = broad_labels[major_classes]
)]
p1e <- ggplot(composition_box,
              aes(compartment_label, fraction, fill = cancer_label, colour = cancer_label)) +
  geom_boxplot(position = position_dodge(width = 0.72), width = 0.60,
               linewidth = 0.35, outlier.shape = NA, alpha = 0.28) +
  geom_point(position = position_jitterdodge(jitter.width = 0.10, dodge.width = 0.72),
             size = 0.85, alpha = 0.78, stroke = 0) +
  facet_wrap(~broad_class_label, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = cancer_cols, name = NULL) +
  scale_colour_manual(values = cancer_cols, name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0.03, 0.10))) +
  labs(title = "Major-lineage fractions by anatomical compartment",
       subtitle = "Boxplots show medians and interquartile ranges; points are individual specimens",
       x = NULL, y = "Cell fraction") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "bottom", strip.text = element_text(size = 6.9))

figure1 <- (p1a | p1b) / p1c / p1d / p1e +
  plot_layout(heights = if (use_upgraded_umap) c(1.34, 0.88, 0.86, 1.12) else
                c(1.10, 0.92, 0.90, 1.16),
              widths = if (use_upgraded_umap) c(0.72, 1.28) else c(0.82, 1.18)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(family = "Arial", size = 11.5, face = "bold"),
      plot.margin = margin(4, 7, 5, 7)
    )
  )

cairo_png <- function(filename, width, height, bg = "white", ...) {
  grDevices::png(filename, width = width, height = height, units = "in", res = 450,
                 type = "cairo", bg = bg, ...)
}
png_path <- file.path(figure_out, "Figure1_literature_standard_atlas.png")
pdf_path <- file.path(figure_out, "Figure1_literature_standard_atlas.pdf")
ggsave(png_path, figure1, width = 7.2, height = 9.1, units = "in",
       device = cairo_png, bg = "white", limitsize = FALSE)
ggsave(pdf_path, figure1, width = 7.2, height = 9.1, units = "in",
       device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE)

fwrite(availability, file.path(source_out, "Figure1A_cohort_matrix.tsv"), sep = "\t")
umap_source <- if (use_upgraded_umap) {
  umap[, .(cell_id, accession, cancer, sample_id, patient_id, compartment,
           broad_class, broad_class_label, umap_1, umap_2, plot_x, plot_y)]
} else {
  umap[, .(cell_id, accession, cancer, sample_id, patient_id, compartment,
           broad_class, broad_class_label, umap_1, umap_2)]
}
fwrite(umap_source, file.path(source_out, "Figure1B_umap_coordinates.tsv.gz"),
       sep = "\t", compress = "gzip")
fwrite(markers, file.path(source_out, "Figure1C_marker_dotplot.tsv"), sep = "\t")
fwrite(composition, file.path(source_out, "Figure1D_sample_composition.tsv"), sep = "\t")
fwrite(composition_box, file.path(source_out, "Figure1E_major_lineage_distributions.tsv"), sep = "\t")

audit <- data.table(
  panel = c("A", "B", "C", "D", "E"),
  visual_form = c("cohort_matrix", "UMAP", "marker_DotPlot", "stacked_composition_bar",
                  "boxplot_with_sample_points"),
  bespoke_form = FALSE,
  connected_patient_lines = FALSE,
  main_text_eligible = TRUE
)
fwrite(audit, file.path(admin_out, "GATE12AZ_FIGURE1_VISUAL_AUDIT.tsv"), sep = "\t")
writeLines(c(
  "# Gate12AZ Figure 1 build record",
  "",
  "- Every panel uses a conventional single-cell atlas visual form.",
  "- No arrow-box workflow, patient trajectory, rule-retention, or bespoke diagnostic graphic is present.",
  "- The global panel uses frozen descriptive display coordinates and contains every integrated cell.",
  "- Boxplots are descriptive sample-level displays; patient-aware inference remains in the dedicated statistical figure/table.",
  "- PNG was rendered at 450 dpi; a vector PDF was produced from the same code."
), file.path(admin_out, "GATE12AZ_FIGURE1_BUILD_RECORD.md"))

cat(if (use_upgraded_umap) "GATE12BB_FIGURE1_STATUS=COMPLETE\n" else
      "GATE12AZ_FIGURE1_STATUS=COMPLETE\n")
cat("PNG=", png_path, "\n", sep = "")
cat("PDF=", pdf_path, "\n", sep = "")
