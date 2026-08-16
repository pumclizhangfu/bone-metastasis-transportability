#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrastr)
  library(patchwork)
  library(grid)
})

options(stringsAsFactors = FALSE, warn = 1)
set.seed(20260811L)
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) normalizePath(args[[1L]], mustWork = TRUE) else
  normalizePath(".", mustWork = TRUE)
out <- file.path(root, "results", "gate12ay_multiscale_atlas")
source_out <- file.path(out, "source_data")
admin_out <- file.path(out, "admin")
dir.create(source_out, recursive = TRUE, showWarnings = FALSE)
dir.create(admin_out, recursive = TRUE, showWarnings = FALSE)

global_file <- file.path(root, "results", "gate12an_umap_reembedding",
                         "K_compact_spca_coordinates.tsv.gz")
local_files <- setNames(
  file.path(source_out, sprintf("Gate12AY_%s_local_UMAP.tsv.gz",
                               c("t_nk", "b_cell", "myeloid", "niche"))),
  c("t_nk", "b_cell", "myeloid", "niche")
)
required <- c(global_file, local_files)
if (any(!file.exists(required))) stop("Missing multiscale-atlas input: ",
                                      paste(required[!file.exists(required)], collapse = ", "))

broad_order <- c("T_NK", "Myeloid", "B", "Progenitor", "Stromal", "Endothelial",
                 "Osteoclast", "Osteoblast", "Malignant", "Erythroid", "Unassigned")
broad_labels <- c(
  T_NK = "T/NK", Myeloid = "Myeloid", B = "B cell", Progenitor = "Progenitor",
  Stromal = "Stromal", Endothelial = "Endothelial", Osteoclast = "Osteoclast",
  Osteoblast = "Osteoblast", Malignant = "Malignant", Erythroid = "Erythroid",
  Unassigned = "Unassigned"
)
broad_cols <- c(
  T_NK = "#3973B9", Myeloid = "#E68435", B = "#4FA3BE",
  Progenitor = "#C8A42A", Stromal = "#35976D", Endothelial = "#6AB7AE",
  Osteoclast = "#C75B56", Osteoblast = "#8B79B9", Malignant = "#7047A3",
  Erythroid = "#D95B8A", Unassigned = "#B9BDC2"
)

state_cols <- c(
  T_NK = "#4C78A8", CD4_T = "#6C8EBF", CD4_naive = "#9BB7D4",
  CD8_naive = "#4E9FD1", CTL_1 = "#2A6FBB", CTL_2 = "#174A7E",
  CTL_3 = "#0D3B66", NKT = "#7569C9", NK_1 = "#5B5BD6",
  NK_2 = "#918AE1", Treg = "#B279A2", Proliferating_T = "#D95F9C",
  B = "#4FA3BE", Mature_B = "#2C7FB8", Memory_B = "#7FCDBB",
  Pro_B = "#41B6C4", Immature_B = "#A1DAB4",
  Myeloid = "#F28E2B", Monocyte_1 = "#E07B39", Monocyte_2 = "#F4A259",
  Monocyte_3 = "#C66B2E", Macrophage = "#B84A3A", mDC = "#C44E52",
  pDC = "#9C755F", Monocyte_progenitor = "#EDC948", Osteoclast = "#D55E50",
  Progenitor = "#C8A42A", HSPC = "#E3C565", Stromal = "#35976D",
  MSC_1 = "#4DAF7C", MSC_2 = "#238B6B", MSC_3 = "#78C69A",
  Pericyte_1 = "#2A9D8F", Pericyte_2 = "#5AB4AC", Pericyte_3 = "#8DD3C7",
  Endothelial = "#62B6A8", Malignant = "#7047A3", Erythroid = "#D95B8A",
  Osteoblast = "#8B79B9", Unassigned = "#B9BDC2"
)

clean_state_label <- function(x) gsub("_", " ", x, fixed = TRUE)

draw_ordered <- function(dt, group_col) {
  counts <- dt[, .N, by = c(group_col)]
  ordered_groups <- counts[order(-N), get(group_col)]
  dt[, draw_group := match(get(group_col), ordered_groups)]
  dt[, draw_random := runif(.N)]
  setorderv(dt, c("draw_group", "draw_random"))
  dt
}

panel_theme <- function() {
  theme_void(base_family = "Arial", base_size = 8.5) +
    theme(
      plot.title = element_text(size = 10.2, face = "bold", hjust = 0,
                                margin = margin(b = 1.2)),
      plot.subtitle = element_text(size = 7.2, colour = "#59636E", hjust = 0,
                                   margin = margin(b = 2.0)),
      plot.margin = margin(4, 5, 4, 5),
      legend.position = "bottom",
      legend.justification = "left",
      legend.text = element_text(size = 6.4, colour = "#30343B"),
      legend.key.height = unit(2.4, "mm"),
      legend.key.width = unit(2.4, "mm"),
      legend.spacing.x = unit(0.6, "mm"),
      legend.margin = margin(t = 1.5),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

global <- fread(global_file)
global[, broad_class := factor(broad_class, levels = broad_order)]
global <- draw_ordered(global, "broad_class")
p_global <- ggplot(global, aes(umap_1, umap_2, colour = broad_class,
                               alpha = broad_class)) +
  ggrastr::geom_point_rast(size = 0.024, raster.dpi = 800) +
  scale_colour_manual(values = broad_cols, breaks = broad_order,
                      labels = broad_labels[broad_order], drop = FALSE) +
  scale_alpha_manual(values = c(setNames(rep(0.62, length(broad_order) - 1L),
                                         broad_order[-length(broad_order)]),
                                Unassigned = 0.15), guide = "none") +
  coord_equal() +
  labs(title = "Integrated atlas",
       subtitle = sprintf("%s cells | all 42 bone specimens",
                          format(nrow(global), big.mark = ",")), colour = NULL) +
  panel_theme() +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE,
                               override.aes = list(size = 2.6, alpha = 1)))

local_titles <- c(
  t_nk = "T and NK cells",
  b_cell = "B-cell compartment",
  myeloid = "Myeloid compartment",
  niche = "Bone-niche and tumor compartment"
)

local_plot <- function(panel_name) {
  dt <- fread(local_files[[panel_name]])
  present_states <- dt[, .N, by = harmonized_state][order(-N), harmonized_state]
  missing_colours <- setdiff(present_states, names(state_cols))
  if (length(missing_colours)) stop("Missing state colours: ", paste(missing_colours, collapse = ", "))
  dt[, state_factor := factor(harmonized_state, levels = present_states)]
  dt <- draw_ordered(dt, "harmonized_state")
  legend_rows <- if (length(present_states) >= 13L) 3L else 2L
  point_size <- if (nrow(dt) < 10000L) 0.050 else if (nrow(dt) < 20000L) 0.040 else 0.030

  ggplot(dt, aes(local_umap_1, local_umap_2, colour = state_factor)) +
    ggrastr::geom_point_rast(size = point_size, alpha = 0.70, raster.dpi = 800) +
    scale_colour_manual(
      values = state_cols[present_states], breaks = present_states,
      labels = clean_state_label(present_states), drop = FALSE
    ) +
    coord_equal() +
    labs(
      title = local_titles[[panel_name]],
      subtitle = sprintf("%s cells | lineage-local UMAP",
                         format(nrow(dt), big.mark = ",")),
      colour = NULL
    ) +
    panel_theme() +
    guides(colour = guide_legend(nrow = legend_rows, byrow = TRUE,
                                 override.aes = list(size = 2.5, alpha = 1)))
}

p_t_nk <- local_plot("t_nk")
p_b_cell <- local_plot("b_cell")
p_myeloid <- local_plot("myeloid")
p_niche <- local_plot("niche")

atlas_design <- "
ABC
ADE
"
atlas <- (p_global + p_t_nk + p_b_cell + p_myeloid + p_niche) +
  plot_layout(design = atlas_design, widths = c(1.32, 1, 1)) +
  plot_annotation(
    title = "Multiscale single-cell atlas of the bone-metastatic ecosystem",
    subtitle = "A global integrated view is paired with lineage-local embeddings to resolve cellular states without altering global coordinates",
    tag_levels = "a",
    theme = theme(
      plot.title = element_text(family = "Arial", size = 14.0, face = "bold", hjust = 0),
      plot.subtitle = element_text(family = "Arial", size = 8.2, colour = "#4B5563", hjust = 0),
      plot.tag = element_text(family = "Arial", size = 12.0, face = "bold"),
      plot.margin = margin(7, 7, 6, 7)
    )
  )

cairo_png <- function(filename, width, height, bg = "white", ...) {
  grDevices::png(filename, width = width, height = height, units = "in", res = 450,
                 type = "cairo", bg = bg, ...)
}
ggsave(file.path(out, "Gate12AY_multiscale_cell_atlas_v2.pdf"), atlas,
       width = 14.2, height = 8.8, units = "in", device = grDevices::cairo_pdf,
       bg = "white")
ggsave(file.path(out, "Gate12AY_multiscale_cell_atlas_v2.png"), atlas,
       width = 14.2, height = 8.8, units = "in", device = cairo_png,
       bg = "white")

state_counts <- rbindlist(lapply(names(local_files), function(panel_name) {
  fread(local_files[[panel_name]], select = c("harmonized_state"))[
    , .N, by = harmonized_state][, panel := panel_name]
}), use.names = TRUE)
setcolorder(state_counts, c("panel", "harmonized_state", "N"))
fwrite(state_counts, file.path(source_out, "Gate12AY_state_counts.tsv"), sep = "\t")

writeLines(c(
  "# Gate12AY multiscale atlas",
  "",
  "- Global panel: K_compact_spca from the frozen 30-dimensional Harmony representation; all 107,886 cells.",
  "- Local panels: de novo UMAPs recomputed within T/NK, B-cell, myeloid and niche/tumor subsets from the same frozen Harmony representation.",
  "- Local UMAPs are visual magnifications, not evidence that distances are comparable across panels.",
  "- No global UMAP island was translated, packed, clipped or removed.",
  "- All panels use deterministic seeds and render every cell in the corresponding subset.",
  "- The local panels address visual crowding while retaining the audited global atlas as the overview."
), file.path(admin_out, "GATE12AY_MULTISCALE_ATLAS_RECORD.md"))

cat("GATE12AY_STATUS=COMPLETE\n")
