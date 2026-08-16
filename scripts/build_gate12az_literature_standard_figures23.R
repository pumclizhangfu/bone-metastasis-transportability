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

required <- c(
  "results/gate12b_cell_states/gate12b_cell_state_coordinates.tsv.gz",
  "results/gate12b_cell_states/patient_state_composition.tsv",
  "results/gate12b_cell_states/state_meta_effects.tsv",
  "results/gate12b_cell_states/cluster_module_scores.tsv",
  "results/gate12b_cell_states/cluster_annotation_audit.tsv",
  "results/gate12ay_multiscale_atlas/source_data/Gate12AY_myeloid_local_UMAP.tsv.gz",
  "results/gate12ay_multiscale_atlas/source_data/Gate12AY_t_nk_local_UMAP.tsv.gz",
  "results/gate12ai_submission_ready/source_data/Figure3B_gene_concordance.tsv.gz",
  "results/gate12ad_figure_restructure/phase_b_main_figures/source_data/Figure3B_hallmark_bubble_heatmap.tsv"
)
missing <- required[!file.exists(file.path(root, required))]
if (length(missing)) stop("Missing Figure 2/3 input: ", paste(missing, collapse = ", "))
if (use_upgraded_umap) {
  upgraded_required <- c(
    "results/gate12ba_umap_upgrade/source_data/Figure2A_myeloid_UMAP_coordinates.tsv.gz",
    "results/gate12ba_umap_upgrade/source_data/Figure2A_myeloid_label_anchors.tsv",
    "results/gate12ba_umap_upgrade/source_data/Figure3A_tnk_UMAP_coordinates.tsv.gz",
    "results/gate12ba_umap_upgrade/source_data/Figure3A_tnk_label_anchors.tsv"
  )
  upgraded_missing <- upgraded_required[!file.exists(file.path(root, upgraded_required))]
  if (length(upgraded_missing)) {
    stop("Missing upgraded UMAP input: ", paste(upgraded_missing, collapse = ", "))
  }
}

read_rel <- function(x, ...) fread(file.path(root, x), ...)
state_coords <- read_rel(required[[1]],
                         select = c("barcode", "accession", "cancer", "sample_id", "patient_id",
                                    "compartment", "lineage", "gate12b_state"))
composition <- read_rel(required[[2]])
effects <- read_rel(required[[3]])
module_scores <- read_rel(required[[4]])
cluster_audit <- read_rel(required[[5]], select = c("cluster", "state", "lineage", "cells"))
myeloid_local <- read_rel(required[[6]],
                          select = c("cell_id", "local_umap_1", "local_umap_2"))
tnk_local <- read_rel(required[[7]],
                      select = c("cell_id", "local_umap_1", "local_umap_2"))
gene_meta <- read_rel(required[[8]])
hallmarks <- read_rel(required[[9]])

state_cols <- c(
  Classical_monocyte = "#E69F00", Inflammatory_monocyte = "#D55E00",
  C1QC_macrophage = "#0072B2", Resident_macrophage = "#009E73",
  cDC = "#56B4E9", pDC = "#7B3294", Proliferating_myeloid = "#111111",
  Unresolved_myeloid = "#B8B8B8", CD4_naive = "#56B4E9",
  CD4_memory = "#0072B2", Treg = "#CC79A7", CD8_effector = "#E69F00",
  CD8_exhausted = "#A50F15", NK_adaptive = "#009E73",
  Proliferating_T_NK = "#111111", Unresolved_T_NK = "#B8B8B8"
)
cancer_cols <- c(prostate = "#3973B9", renal = "#E68435")

myeloid_order <- c("Classical_monocyte", "Inflammatory_monocyte", "C1QC_macrophage",
                   "Resident_macrophage", "cDC", "pDC", "Proliferating_myeloid",
                   "Unresolved_myeloid")
tnk_order <- c("CD4_naive", "CD4_memory", "Treg", "CD8_effector", "CD8_exhausted",
               "NK_adaptive", "Proliferating_T_NK", "Unresolved_T_NK")
clean_state <- function(x) {
  out <- gsub("_", " ", x, fixed = TRUE)
  out <- gsub("T NK", "T/NK", out, fixed = TRUE)
  out
}
clean_pathway <- function(x) {
  out <- sub("^HALLMARK_", "", x)
  out <- gsub("_", " ", out)
  out <- tools::toTitleCase(tolower(out))
  out <- gsub("Tnfa Signaling Via Nfkb", "TNF/NF-kB", out, fixed = TRUE)
  out <- gsub("Il2 Stat5 Signaling", "IL2/STAT5", out, fixed = TRUE)
  out <- gsub("Mtorc1 Signaling", "mTORC1 signaling", out, fixed = TRUE)
  out <- gsub("Interferon Gamma Response", "IFN-gamma response", out, fixed = TRUE)
  out
}

theme_pub <- theme_classic(base_size = 8.4, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 9.2, face = "bold", hjust = 0, margin = margin(b = 2)),
    plot.subtitle = element_text(size = 7.1, colour = "#525A61", hjust = 0,
                                 lineheight = 0.95, margin = margin(b = 3)),
    axis.title = element_text(size = 7.8),
    axis.text = element_text(size = 6.9, colour = "#30343B"),
    strip.text = element_text(size = 7.0, face = "bold"),
    strip.background = element_rect(fill = "#F3F4F5", colour = NA),
    legend.title = element_text(size = 7.0), legend.text = element_text(size = 6.6),
    legend.key.height = unit(2.6, "mm"), legend.key.width = unit(2.6, "mm"),
    panel.grid = element_blank(), plot.tag = element_text(size = 11.5, face = "bold"),
    plot.margin = margin(4, 4, 4, 4)
  )

merge_local <- function(lineage_value, local_dt) {
  x <- merge(
    state_coords[lineage == lineage_value], local_dt,
    by.x = "barcode", by.y = "cell_id", all = FALSE, sort = FALSE
  )
  expected <- state_coords[lineage == lineage_value, uniqueN(barcode)]
  if (nrow(x) != expected || uniqueN(x$barcode) != expected) {
    stop("Local UMAP/state join failed for ", lineage_value, ": ", nrow(x), " / ", expected)
  }
  x
}
myeloid_umap <- merge_local("Myeloid", myeloid_local)
tnk_umap <- merge_local("T_NK", tnk_local)
if (use_upgraded_umap) {
  myeloid_upgrade <- read_rel(
    "results/gate12ba_umap_upgrade/source_data/Figure2A_myeloid_UMAP_coordinates.tsv.gz",
    select = c("barcode", "plot_x", "plot_y")
  )
  tnk_upgrade <- read_rel(
    "results/gate12ba_umap_upgrade/source_data/Figure3A_tnk_UMAP_coordinates.tsv.gz",
    select = c("barcode", "plot_x", "plot_y")
  )
  myeloid_anchors <- read_rel(
    "results/gate12ba_umap_upgrade/source_data/Figure2A_myeloid_label_anchors.tsv"
  )
  tnk_anchors <- read_rel(
    "results/gate12ba_umap_upgrade/source_data/Figure3A_tnk_label_anchors.tsv"
  )
  myeloid_umap <- merge(myeloid_umap, myeloid_upgrade, by = "barcode",
                        all.x = TRUE, sort = FALSE)
  tnk_umap <- merge(tnk_umap, tnk_upgrade, by = "barcode",
                    all.x = TRUE, sort = FALSE)
  if (any(!is.finite(myeloid_umap$plot_x)) || any(!is.finite(myeloid_umap$plot_y)) ||
      any(!is.finite(tnk_umap$plot_x)) || any(!is.finite(tnk_umap$plot_y))) {
    stop("Upgraded lineage UMAP join failed")
  }
  ## Separate two adjacent labels without changing any UMAP coordinate or the
  ## observed density-mode anchors to which the labels refer.
  myeloid_anchors[, `:=`(
    anchor_x = plot_x, anchor_y = plot_y,
    label_x = plot_x, label_y = plot_y,
    manual_label = FALSE
  )]
  myeloid_anchors[group == "Classical_monocyte",
                  `:=`(label_x = -5.05, label_y = 3.05, manual_label = TRUE)]
  myeloid_anchors[group == "Resident_macrophage",
                  `:=`(label_x = -0.05, label_y = 0.15, manual_label = TRUE)]
  tnk_anchors[, `:=`(
    anchor_x = plot_x, anchor_y = plot_y,
    label_x = plot_x, label_y = plot_y,
    manual_label = FALSE
  )]
  tnk_anchors[group == "NK_adaptive",
              `:=`(label_x = -5.70, label_y = 2.55, manual_label = TRUE)]
  tnk_anchors[group == "CD8_effector",
              `:=`(label_x = -2.45, label_y = 1.00, manual_label = TRUE)]
  tnk_anchors[group == "CD8_exhausted",
              `:=`(label_x = -3.35, label_y = -5.40, manual_label = TRUE)]
  tnk_anchors[group == "Proliferating_T_NK",
              `:=`(label_x = -6.05, label_y = -6.78, manual_label = TRUE)]
}

make_umap <- function(dt, state_order, title, upgraded_anchors = NULL) {
  dt <- copy(dt)
  dt[, gate12b_state := factor(gate12b_state, levels = state_order)]
  counts <- dt[, .N, by = gate12b_state][order(-N)]
  dt[, draw_group := match(gate12b_state, counts$gate12b_state)]
  dt[, draw_random := runif(.N)]
  setorder(dt, draw_group, draw_random)
  alpha_vals <- setNames(rep(0.70, length(state_order)), state_order)
  alpha_vals[grepl("^Unresolved", names(alpha_vals))] <- 0.16

  if (!is.null(upgraded_anchors)) {
    point_size <- if (nrow(dt) > 50000L) 0.080 else 0.115
    has_manual <- "manual_label" %in% names(upgraded_anchors)
    direct_anchors <- if (has_manual) upgraded_anchors[manual_label == FALSE] else upgraded_anchors
    manual_anchors <- if (has_manual) upgraded_anchors[manual_label == TRUE] else NULL
    p <- ggplot(dt, aes(plot_x, plot_y, colour = gate12b_state,
                     alpha = gate12b_state)) +
        ggrastr::geom_point_rast(size = point_size, raster.dpi = 900, stroke = 0) +
        ggrepel::geom_text_repel(
          data = direct_anchors,
          aes(plot_x, plot_y, label = label), inherit.aes = FALSE,
          family = "Arial", fontface = "bold",
          size = 2.08,
          colour = "#202124", box.padding = 0.08,
          point.padding = 0.02, min.segment.length = Inf, segment.colour = NA,
          max.overlaps = Inf, force = 0.20, max.time = 2, max.iter = 10000,
          seed = 20260812L, show.legend = FALSE
        ) +
        scale_colour_manual(values = state_cols[state_order], breaks = state_order,
                            labels = clean_state(state_order), drop = FALSE, name = NULL) +
        scale_alpha_manual(values = alpha_vals, guide = "none") +
        scale_x_continuous(expand = expansion(mult = 0.012)) +
        scale_y_continuous(expand = expansion(mult = 0.012)) +
        coord_equal(clip = "off") +
        labs(title = title,
             subtitle = sprintf("%s cells | lineage-local UMAP", comma(nrow(dt)))) +
        theme_void(base_size = 8.4, base_family = "Arial") +
        theme(plot.title = element_text(size = 9.2, face = "bold", hjust = 0),
              plot.subtitle = element_text(size = 7.1, colour = "#525A61", hjust = 0,
                                           margin = margin(b = 3)),
              legend.position = "none", plot.margin = margin(4, 4, 4, 4))
    if (has_manual && nrow(manual_anchors)) {
      p <- p +
        geom_segment(
          data = manual_anchors,
          aes(x = anchor_x, y = anchor_y, xend = label_x, yend = label_y,
              colour = group),
          inherit.aes = FALSE, linewidth = 0.30, alpha = 0.88,
          lineend = "round", show.legend = FALSE
        ) +
        geom_label(
          data = manual_anchors,
          aes(x = label_x, y = label_y, label = label),
          inherit.aes = FALSE, family = "Arial", fontface = "bold",
          size = 2.08, colour = "#202124", fill = scales::alpha("white", 0.90),
          linewidth = 0, label.padding = unit(0.055, "lines"),
          show.legend = FALSE
        )
    }
    return(p)
  }

  ggplot(dt, aes(local_umap_1, local_umap_2, colour = gate12b_state,
                 alpha = gate12b_state)) +
      ggrastr::geom_point_rast(size = 0.035, raster.dpi = 800) +
      scale_colour_manual(values = state_cols[state_order], breaks = state_order,
                          labels = clean_state(state_order), drop = FALSE, name = NULL) +
      scale_alpha_manual(values = alpha_vals, guide = "none") +
      coord_equal() +
      labs(title = title,
           subtitle = sprintf("%s cells | lineage-local UMAP", comma(nrow(dt)))) +
      theme_void(base_size = 8.4, base_family = "Arial") +
      theme(plot.title = element_text(size = 9.2, face = "bold", hjust = 0),
            plot.subtitle = element_text(size = 7.1, colour = "#525A61", hjust = 0,
                                         margin = margin(b = 3)),
            legend.position = "bottom", legend.justification = "left",
            legend.text = element_text(size = 6.4),
            legend.key.height = unit(2.4, "mm"), legend.key.width = unit(2.4, "mm"),
            plot.margin = margin(4, 4, 4, 4)) +
      guides(colour = guide_legend(ncol = 2, byrow = TRUE,
                                   override.aes = list(size = 2.4, alpha = 1)))
}

module_matrix <- merge(module_scores, cluster_audit,
                       by = c("lineage", "cluster"), all.x = TRUE, sort = FALSE)
module_matrix <- module_matrix[!is.na(state)]
module_matrix <- module_matrix[, .(
  module_z = weighted.mean(module_z, w = cells),
  raw_module_score = weighted.mean(raw_module_score, w = cells)
), by = .(lineage, state, module)]

make_module_heatmap <- function(lineage_value, state_order, title) {
  x <- copy(module_matrix[lineage == lineage_value & !grepl("^Unresolved", state)])
  shown_states <- state_order[!grepl("^Unresolved", state_order)]
  shown_modules <- shown_states[shown_states %in% unique(x$module)]
  x <- x[state %in% shown_states & module %in% shown_modules]
  x[, state_label := factor(clean_state(state), levels = rev(clean_state(shown_states)))]
  x[, module_label := factor(clean_state(module), levels = clean_state(shown_modules))]
  x[, module_z_plot := pmax(-2.5, pmin(2.5, module_z))]

  ggplot(x, aes(module_label, state_label, fill = module_z_plot)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#B94450",
                         midpoint = 0, limits = c(-2.5, 2.5), name = "Module z") +
    labs(title = title,
         subtitle = "Rows: assigned states | Columns: frozen annotation modules",
         x = NULL, y = NULL) +
    theme_pub +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
}

make_composition_boxplots <- function(lineage_value, state_order, title) {
  shown_states <- state_order[!grepl("^Unresolved", state_order)]
  x <- copy(composition[lineage == lineage_value & state %in% shown_states])
  x[, compartment := factor(compartment, levels = c("distal", "involved", "tumor"),
                            labels = c("Distal", "Involved", "Tumor"))]
  x[, state_label := factor(clean_state(state), levels = clean_state(shown_states))]

  ggplot(x, aes(compartment, fraction, fill = cancer, colour = cancer)) +
    geom_boxplot(position = position_dodge(width = 0.72), width = 0.60,
                 linewidth = 0.34, outlier.shape = NA, alpha = 0.28) +
    geom_point(position = position_jitterdodge(jitter.width = 0.10, dodge.width = 0.72),
               size = 0.78, alpha = 0.78, stroke = 0) +
    facet_wrap(~state_label, scales = "free_y", nrow = 2) +
    scale_fill_manual(values = cancer_cols, labels = c(prostate = "Prostate", renal = "Renal"),
                      name = NULL) +
    scale_colour_manual(values = cancer_cols, labels = c(prostate = "Prostate", renal = "Renal"),
                        name = NULL) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0.03, 0.10))) +
    labs(title = title,
         subtitle = "Boxplots summarize specimens; points are individual samples",
         x = NULL, y = "Within-lineage fraction") +
    theme_pub +
    theme(axis.text.x = element_text(angle = 25, hjust = 1),
          strip.text = element_text(size = 6.3), legend.position = "bottom")
}

save_figure <- function(plot, filename, height) {
  cairo_png <- function(filename, width, height, bg = "white", ...) {
    grDevices::png(filename, width = width, height = height, units = "in", res = 450,
                   type = "cairo", bg = bg, ...)
  }
  ggsave(file.path(figure_out, paste0(filename, ".png")), plot, width = 7.2,
         height = height, units = "in", device = cairo_png, bg = "white", limitsize = FALSE)
  ggsave(file.path(figure_out, paste0(filename, ".pdf")), plot, width = 7.2,
         height = height, units = "in", device = grDevices::cairo_pdf,
         bg = "white", limitsize = FALSE)
}

## ----------------------------------------------------------------------
## Figure 2: myeloid state remodelling
## ----------------------------------------------------------------------
p2a <- make_umap(myeloid_umap, myeloid_order, "Myeloid state atlas",
                 if (use_upgraded_umap) myeloid_anchors else NULL)
p2b <- make_module_heatmap("Myeloid", myeloid_order, "Marker-module annotation")
p2c <- make_composition_boxplots("Myeloid", myeloid_order,
                                 "Myeloid-state abundance across compartments")

myeloid_effects <- copy(effects[lineage == "Myeloid" & !grepl("^Unresolved", state)])
myeloid_effects[, state_label := factor(clean_state(state),
                                        levels = rev(clean_state(myeloid_order[!grepl("^Unresolved", myeloid_order)])))]
p2d <- ggplot(myeloid_effects, aes(odds_ratio_per_step, state_label)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.35, colour = "#777777") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y",
                width = 0.13, linewidth = 0.48, colour = "#3F4A54") +
  geom_point(aes(fill = cross_cancer_stable), shape = 21, size = 2.0,
             stroke = 0.45, colour = "#2F3439") +
  scale_fill_manual(values = c(`TRUE` = "#3973B9", `FALSE` = "white"),
                    labels = c(`TRUE` = "Stable", `FALSE` = "Not stable"), name = NULL) +
  scale_x_log10() +
  labs(title = "Patient-level composition effects",
       subtitle = "Meta-analytic odds ratio per anatomical step",
       x = "Odds ratio per step (95% CI, log scale)", y = NULL) +
  theme_pub + theme(legend.position = "bottom")

cohort_effects <- read_rel("results/gate12b_cell_states/state_cohort_effects.tsv")
myeloid_concordance <- dcast(
  cohort_effects[lineage == "Myeloid"], state ~ cancer, value.var = "beta_per_step"
)
myeloid_concordance <- merge(
  myeloid_concordance,
  effects[lineage == "Myeloid", .(state, cross_cancer_stable)], by = "state", all.x = TRUE
)
myeloid_concordance[, state_label := fcase(
  state == "Classical_monocyte", "Classical mono.",
  state == "Inflammatory_monocyte", "Inflammatory mono.",
  state == "C1QC_macrophage", "C1QC macro.",
  state == "Resident_macrophage", "Resident macro.",
  state == "Proliferating_myeloid", "Proliferating",
  default = clean_state(state)
)]
p2e <- ggplot(myeloid_concordance, aes(prostate, renal)) +
  geom_hline(yintercept = 0, linewidth = 0.28, colour = "#B3B3B3") +
  geom_vline(xintercept = 0, linewidth = 0.28, colour = "#B3B3B3") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.35,
              colour = "#777777") +
  geom_point(aes(fill = cross_cancer_stable), shape = 21, size = 2.2,
             stroke = 0.45, colour = "#2F3439") +
  ggrepel::geom_text_repel(data = myeloid_concordance,
                           aes(label = state_label), size = 1.90, seed = 20260811L,
                           max.overlaps = Inf, min.segment.length = 0,
                           segment.size = 0.22, box.padding = 0.30,
                           point.padding = 0.12, force = 1.2) +
  scale_fill_manual(values = c(`TRUE` = "#3973B9", `FALSE` = "white"),
                    breaks = c(FALSE, TRUE), labels = c("Not stable", "Stable"),
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = 0.16)) +
  scale_y_continuous(expand = expansion(mult = 0.16)) +
  coord_cartesian(clip = "off") +
  labs(title = "Cross-cancer concordance",
       subtitle = "Cohort-specific log-odds effects",
       x = "Prostate", y = "Renal") + theme_pub +
  theme(legend.position = "bottom", plot.margin = margin(4, 9, 4, 5)) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 2.4)))

figure2 <- (p2a | p2b) / p2c / (p2d | p2e) +
  plot_layout(heights = c(1.12, 1.18, 0.90), widths = c(1.02, 0.98)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(family = "Arial", size = 11.5, face = "bold"),
      plot.margin = margin(4, 7, 5, 7)
    )
  )
save_figure(figure2, "Figure2_literature_standard_myeloid", 8.3)

## ----------------------------------------------------------------------
## Figure 3: T/NK state and transcriptional remodelling
## ----------------------------------------------------------------------
p3a <- make_umap(tnk_umap, tnk_order, "T/NK state atlas",
                 if (use_upgraded_umap) tnk_anchors else NULL)
p3b <- make_module_heatmap("T_NK", tnk_order, "Marker-module annotation")
p3c <- make_composition_boxplots("T_NK", tnk_order,
                                 "T/NK-state abundance across compartments")

gene_plot <- copy(gene_meta[is.finite(meta_logFC) & is.finite(meta_FDR)])
gene_plot[, state_label := factor(state, levels = unique(state))]
gene_plot[, neglog10_fdr := pmin(60, -log10(pmax(meta_FDR, 1e-60)))]
gene_plot[, category := fifelse(meta_FDR < 0.05 & meta_logFC >= 1, "Up",
                                fifelse(meta_FDR < 0.05 & meta_logFC <= -1, "Down", "Not significant"))]
label_genes <- gene_plot[meta_FDR < 0.05][order(state, -abs(meta_z)), head(.SD, 3), by = state]
p3d <- ggplot(gene_plot, aes(meta_logFC, neglog10_fdr)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.30, colour = "#888888") +
  geom_vline(xintercept = c(-1, 1), linetype = 2, linewidth = 0.30, colour = "#888888") +
  ggrastr::geom_point_rast(aes(colour = category), size = 0.34, alpha = 0.52,
                           raster.dpi = 600) +
  ggrepel::geom_text_repel(data = label_genes, aes(label = gene), size = 2.0,
                           seed = 20260811L, max.overlaps = Inf, min.segment.length = 0,
                           segment.size = 0.20, box.padding = 0.20) +
  facet_wrap(~state_label, nrow = 1) +
  scale_colour_manual(values = c(Down = "#3973B9", `Not significant` = "#BFC3C7", Up = "#C44E52"),
                      name = NULL) +
  labs(title = "Cross-cancer pseudobulk transcriptional changes",
       subtitle = "Meta-analysis by T/NK lineage state; labels mark the strongest genes",
       x = "Meta-analytic log2 fold change", y = expression(-log[10](FDR))) +
  theme_pub + theme(legend.position = "bottom", strip.text = element_text(size = 6.7))

hallmark_plot <- copy(hallmarks[robust_conserved == TRUE])
hallmark_plot[, state := factor(state, levels = c("CD4 T", "CD8/CTL", "NK/NKT"))]
if (!"pathway_label" %in% names(hallmark_plot)) {
  hallmark_plot[, pathway_label := clean_pathway(pathway)]
}
path_order <- hallmark_plot[, .(mean_abs = mean(abs(mean_nes))), by = pathway_label][
  order(mean_abs), pathway_label]
hallmark_plot[, pathway_label := factor(pathway_label, levels = path_order)]
hallmark_plot[, fdr_strength_plot := pmin(35, -log10(pmax(combined_FDR, 1e-35)))]
p3e <- ggplot(hallmark_plot, aes(state, pathway_label)) +
  geom_point(aes(size = fdr_strength_plot, fill = mean_nes), shape = 21,
             colour = "#444444", stroke = 0.18) +
  scale_fill_gradient(low = "#F4F7FB", high = "#B94450", name = "Mean NES") +
  scale_size_continuous(range = c(1.2, 4.2), name = expression(-log[10](FDR))) +
  labs(title = "Conserved Hallmark programmes",
       subtitle = "Concordant in prostate and renal cohorts",
       x = NULL, y = NULL) +
  theme_pub + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                    legend.position = "right")

figure3 <- (p3a | p3b) / p3c / (p3d | p3e) +
  plot_layout(heights = c(1.10, 1.16, 1.02), widths = c(1.12, 0.88)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(family = "Arial", size = 11.5, face = "bold"),
      plot.margin = margin(4, 7, 5, 7)
    )
  )
save_figure(figure3, "Figure3_literature_standard_tnk", 8.7)

fwrite(myeloid_umap, file.path(source_out, "Figure2A_myeloid_umap.tsv.gz"), sep = "\t", compress = "gzip")
if (use_upgraded_umap) {
  fwrite(myeloid_anchors, file.path(source_out, "Figure2A_myeloid_label_positions.tsv"), sep = "\t")
}
fwrite(module_matrix[lineage == "Myeloid"], file.path(source_out, "Figure2B_myeloid_module_heatmap.tsv"), sep = "\t")
fwrite(composition[lineage == "Myeloid"], file.path(source_out, "Figure2C_myeloid_composition.tsv"), sep = "\t")
fwrite(myeloid_effects, file.path(source_out, "Figure2D_myeloid_effects.tsv"), sep = "\t")
fwrite(myeloid_concordance, file.path(source_out, "Figure2E_myeloid_concordance.tsv"), sep = "\t")
fwrite(tnk_umap, file.path(source_out, "Figure3A_tnk_umap.tsv.gz"), sep = "\t", compress = "gzip")
if (use_upgraded_umap) {
  fwrite(tnk_anchors, file.path(source_out, "Figure3A_tnk_label_positions.tsv"), sep = "\t")
}
fwrite(module_matrix[lineage == "T_NK"], file.path(source_out, "Figure3B_tnk_module_heatmap.tsv"), sep = "\t")
fwrite(composition[lineage == "T_NK"], file.path(source_out, "Figure3C_tnk_composition.tsv"), sep = "\t")
fwrite(gene_plot, file.path(source_out, "Figure3D_tnk_gene_meta.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(hallmark_plot, file.path(source_out, "Figure3E_tnk_hallmarks.tsv"), sep = "\t")

audit <- data.table(
  figure = rep(c("Figure2", "Figure3"), each = 5),
  panel = rep(LETTERS[1:5], 2),
  visual_form = c("UMAP", "module_heatmap", "boxplot_with_sample_points",
                  "forest_plot", "concordance_scatter",
                  "UMAP", "module_heatmap", "boxplot_with_sample_points",
                  "volcano_plot", "GSEA_bubble_plot"),
  bespoke_form = FALSE,
  connected_patient_lines = FALSE,
  main_text_eligible = TRUE
)
fwrite(audit, file.path(admin_out, "GATE12AZ_FIGURES23_VISUAL_AUDIT.tsv"), sep = "\t")

cat(if (use_upgraded_umap) "GATE12BB_FIGURES23_STATUS=COMPLETE\n" else
      "GATE12AZ_FIGURES23_STATUS=COMPLETE\n")
