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
set.seed(20260812L)

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) normalizePath(args[[1L]], mustWork = TRUE) else
  normalizePath(".", mustWork = TRUE)
out <- file.path(root, "results", "gate12ba_umap_upgrade")
panel_out <- file.path(out, "panels")
source_out <- file.path(out, "source_data")
admin_out <- file.path(out, "admin")
dir.create(panel_out, recursive = TRUE, showWarnings = FALSE)
dir.create(source_out, recursive = TRUE, showWarnings = FALSE)
dir.create(admin_out, recursive = TRUE, showWarnings = FALSE)

global_file <- file.path(
  root, "results", "gate12ap_selected_umap", "source_data",
  "Gate12AP_selected_UMAP_coordinates.tsv.gz"
)
state_file <- file.path(
  root, "results", "gate12b_cell_states", "gate12b_cell_state_coordinates.tsv.gz"
)
myeloid_file <- file.path(
  root, "results", "gate12ay_multiscale_atlas", "source_data",
  "Gate12AY_myeloid_local_UMAP.tsv.gz"
)
tnk_file <- file.path(
  root, "results", "gate12ay_multiscale_atlas", "source_data",
  "Gate12AY_t_nk_local_UMAP.tsv.gz"
)
required <- c(global_file, state_file, myeloid_file, tnk_file)
if (any(!file.exists(required))) {
  stop("Missing UMAP input: ", paste(required[!file.exists(required)], collapse = ", "))
}

global <- fread(global_file)
state <- fread(
  state_file,
  select = c("barcode", "lineage", "gate12b_state")
)
myeloid_local <- fread(
  myeloid_file,
  select = c("cell_id", "local_umap_1", "local_umap_2")
)
tnk_local <- fread(
  tnk_file,
  select = c("cell_id", "local_umap_1", "local_umap_2")
)

join_state <- function(lineage_value, local_dt) {
  x <- merge(
    state[lineage == lineage_value], local_dt,
    by.x = "barcode", by.y = "cell_id", all = FALSE, sort = FALSE
  )
  expected <- state[lineage == lineage_value, uniqueN(barcode)]
  if (nrow(x) != expected || uniqueN(x$barcode) != expected) {
    stop("State/local-UMAP join failed for ", lineage_value)
  }
  x
}

myeloid <- join_state("Myeloid", myeloid_local)
tnk <- join_state("T_NK", tnk_local)

broad_order <- c(
  "T_NK", "Myeloid", "B", "Progenitor", "Stromal", "Endothelial",
  "Osteoclast", "Osteoblast", "Malignant", "Erythroid", "Unassigned"
)
broad_labels <- c(
  T_NK = "T/NK", Myeloid = "Myeloid", B = "B cell", Progenitor = "Progenitor",
  Stromal = "Stromal", Endothelial = "Endothelial", Osteoclast = "Osteoclast",
  Osteoblast = "Osteoblast", Malignant = "Malignant", Erythroid = "Erythroid",
  Unassigned = "Unassigned"
)
broad_cols <- c(
  T_NK = "#3B73B9", Myeloid = "#E68632", B = "#42A4BE",
  Progenitor = "#C9A227", Stromal = "#2F9971", Endothelial = "#67B9AE",
  Osteoclast = "#C95B55", Osteoblast = "#8D79BD", Malignant = "#7046A5",
  Erythroid = "#D85787", Unassigned = "#C5C9CD"
)

myeloid_order <- c(
  "Classical_monocyte", "Inflammatory_monocyte", "C1QC_macrophage",
  "Resident_macrophage", "cDC", "pDC", "Proliferating_myeloid",
  "Unresolved_myeloid"
)
myeloid_labels <- c(
  Classical_monocyte = "Classical mono.",
  Inflammatory_monocyte = "Inflammatory mono.",
  C1QC_macrophage = "C1QC macrophage",
  Resident_macrophage = "Resident macrophage",
  cDC = "cDC", pDC = "pDC",
  Proliferating_myeloid = "Proliferating",
  Unresolved_myeloid = "Unresolved"
)
myeloid_cols <- c(
  Classical_monocyte = "#E5A000", Inflammatory_monocyte = "#D55E00",
  C1QC_macrophage = "#1479B8", Resident_macrophage = "#169B75",
  cDC = "#58AFE0", pDC = "#8A3AA8", Proliferating_myeloid = "#3A3A3A",
  Unresolved_myeloid = "#C8C8C8"
)

tnk_order <- c(
  "CD4_naive", "CD4_memory", "Treg", "CD8_effector", "CD8_exhausted",
  "NK_adaptive", "Proliferating_T_NK", "Unresolved_T_NK"
)
tnk_labels <- c(
  CD4_naive = "CD4 naive", CD4_memory = "CD4 memory", Treg = "Treg",
  CD8_effector = "CD8 effector", CD8_exhausted = "CD8 exhausted",
  NK_adaptive = "NK adaptive", Proliferating_T_NK = "Proliferating",
  Unresolved_T_NK = "Unresolved"
)
tnk_cols <- c(
  CD4_naive = "#56B4E9", CD4_memory = "#2979B9", Treg = "#CC79A7",
  CD8_effector = "#E69F00", CD8_exhausted = "#B31B21",
  NK_adaptive = "#009E73", Proliferating_T_NK = "#333333",
  Unresolved_T_NK = "#C8C8C8"
)

## UMAP axes are arbitrary. This function applies only one global rigid rotation
## and optional reflections, preserving all pairwise distances exactly.
orient_embedding <- function(dt, x_col, y_col, target_aspect = 1.12,
                             left_group = NULL, right_group = NULL,
                             upper_group = NULL, lower_group = NULL,
                             group_col) {
  ans <- copy(dt)
  xy <- as.matrix(ans[, .SD, .SDcols = c(x_col, y_col)])
  centre <- colMeans(xy)
  xy0 <- sweep(xy, 2L, centre, "-")
  angle_grid <- 0:179
  scores <- vapply(angle_grid, function(deg) {
    th <- deg * pi / 180
    rot <- matrix(c(cos(th), sin(th), -sin(th), cos(th)), nrow = 2L, byrow = TRUE)
    trial <- xy0 %*% rot
    spans <- apply(trial, 2L, function(z) diff(range(z, finite = TRUE)))
    max(spans[1L] / target_aspect, spans[2L])
  }, numeric(1L))
  angle <- angle_grid[which.min(scores)]
  th <- angle * pi / 180
  rot <- matrix(c(cos(th), sin(th), -sin(th), cos(th)), nrow = 2L, byrow = TRUE)
  plot_xy <- xy0 %*% rot
  ans[, `:=`(plot_x = plot_xy[, 1L], plot_y = plot_xy[, 2L])]
  group_values <- ans[[group_col]]
  if (!is.null(left_group) && !is.null(right_group)) {
    if (median(ans[group_values == left_group]$plot_x) >
        median(ans[group_values == right_group]$plot_x)) ans[, plot_x := -plot_x]
  }
  if (!is.null(upper_group) && !is.null(lower_group)) {
    if (median(ans[group_values == upper_group]$plot_y) <
        median(ans[group_values == lower_group]$plot_y)) ans[, plot_y := -plot_y]
  }
  attr(ans, "rotation_degrees") <- angle
  ans
}

## An observed high-density cell is used as a reproducible label anchor. The
## displayed label contains no connector line, keeping the standard atlas look.
dense_anchors <- function(dt, group_col, label_map, group_order,
                          unresolved_pattern = "Unassigned|Unresolved") {
  group_counts <- dt[, .N, by = group_col]
  setnames(group_counts, group_col, "anchor_group")
  eligible_groups <- group_counts[
    N >= 50L & !grepl(unresolved_pattern, anchor_group), anchor_group
  ]
  group_order <- group_order[group_order %in% eligible_groups]
  rbindlist(lapply(group_order, function(g) {
    x <- dt[get(group_col) == g]
    if (!nrow(x)) return(NULL)
    take <- min(3500L, nrow(x))
    idx <- if (nrow(x) > take) sample.int(nrow(x), take) else seq_len(nrow(x))
    sample_xy <- x[idx, .(plot_x, plot_y)]
    gx <- 38L
    xr <- range(sample_xy$plot_x, finite = TRUE)
    yr <- range(sample_xy$plot_y, finite = TRUE)
    xb <- cut(sample_xy$plot_x, breaks = seq(xr[1], xr[2], length.out = gx + 1L),
              include.lowest = TRUE, labels = FALSE)
    yb <- cut(sample_xy$plot_y, breaks = seq(yr[1], yr[2], length.out = gx + 1L),
              include.lowest = TRUE, labels = FALSE)
    dens <- data.table(i = seq_len(nrow(sample_xy)), xb = xb, yb = yb)[
      , n_bin := .N, by = .(xb, yb)][order(-n_bin)]
    best <- dens$i[1L]
    data.table(group = g, label = unname(label_map[g]),
               plot_x = sample_xy$plot_x[best], plot_y = sample_xy$plot_y[best])
  }))
}

make_umap <- function(dt, group_col, group_order, label_map, colour_map,
                      title, subtitle, point_size, label_size,
                      unresolved_pattern = "Unassigned|Unresolved") {
  x <- copy(dt)
  x[, plot_group := factor(get(group_col), levels = group_order)]
  counts <- x[, .N, by = plot_group]
  background <- group_order[grepl(unresolved_pattern, group_order)]
  resolved <- group_order[!grepl(unresolved_pattern, group_order)]
  resolved_draw <- as.character(counts[plot_group %in% resolved][order(-N)]$plot_group)
  draw_order <- c(background, resolved_draw)
  x[, draw_group := match(as.character(plot_group), draw_order)]
  x[, draw_random := runif(.N)]
  setorder(x, draw_group, draw_random)
  anchors <- dense_anchors(
    x, "plot_group", label_map, group_order,
    unresolved_pattern = unresolved_pattern
  )
  anchors[, plot_group := factor(group, levels = group_order)]
  alpha_values <- setNames(rep(0.82, length(group_order)), group_order)
  alpha_values[grepl(unresolved_pattern, group_order)] <- 0.22

  p <- ggplot(x, aes(plot_x, plot_y, colour = plot_group, alpha = plot_group)) +
    ggrastr::geom_point_rast(size = point_size, raster.dpi = 900, stroke = 0) +
    ggrepel::geom_text_repel(
      data = anchors,
      aes(plot_x, plot_y, label = label),
      inherit.aes = FALSE, family = "Arial", fontface = "bold", size = label_size,
      colour = "#202124", box.padding = 0.08,
      point.padding = 0.02, min.segment.length = Inf, segment.colour = NA,
      max.overlaps = Inf, force = 0.20, max.time = 2, max.iter = 10000,
      seed = 20260812L, show.legend = FALSE
    ) +
    scale_colour_manual(values = colour_map, breaks = group_order, drop = FALSE) +
    scale_alpha_manual(values = alpha_values, guide = "none") +
    scale_x_continuous(expand = expansion(mult = 0.012)) +
    scale_y_continuous(expand = expansion(mult = 0.012)) +
    coord_equal(clip = "off") +
    labs(title = title, subtitle = subtitle) +
    theme_void(base_size = 8.5, base_family = "Arial") +
    theme(
      plot.title = element_text(size = 9.4, face = "bold", hjust = 0,
                                margin = margin(b = 1.5)),
      plot.subtitle = element_text(size = 7.2, colour = "#59616A", hjust = 0,
                                   margin = margin(b = 2.5)),
      legend.position = "none",
      plot.margin = margin(4, 4, 4, 4),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
  list(plot = p, coordinates = x, anchors = anchors)
}

## Use the vetted K_compact_spca coordinates and rotate the already-vetted
## display coordinates by 90 degrees for a landscape main-figure viewport.
## This is a rigid transform and preserves every pairwise distance.
global[, `:=`(plot_x = plot_umap_2, plot_y = -plot_umap_1)]
if (global[broad_class == "T_NK", median(plot_x)] >
    global[broad_class == "Myeloid", median(plot_x)]) global[, plot_x := -plot_x]
if (global[broad_class == "B", median(plot_y)] <
    global[broad_class == "Stromal", median(plot_y)]) global[, plot_y := -plot_y]
global_result <- make_umap(
  global, "broad_class", broad_order, broad_labels, broad_cols,
  "Integrated cell atlas", sprintf("%s cells | fixed K-compact embedding", comma(nrow(global))),
  point_size = 0.105, label_size = 1.80
)

myeloid <- orient_embedding(
  myeloid, "local_umap_1", "local_umap_2", target_aspect = 1.06,
  left_group = "Classical_monocyte", right_group = "C1QC_macrophage",
  upper_group = "Proliferating_myeloid", lower_group = "Inflammatory_monocyte",
  group_col = "gate12b_state"
)
myeloid_result <- make_umap(
  myeloid, "gate12b_state", myeloid_order, myeloid_labels, myeloid_cols,
  "Myeloid state atlas", sprintf("%s cells | lineage-local UMAP", comma(nrow(myeloid))),
  point_size = 0.115, label_size = 1.72
)

tnk <- orient_embedding(
  tnk, "local_umap_1", "local_umap_2", target_aspect = 1.06,
  left_group = "NK_adaptive", right_group = "CD4_naive",
  upper_group = "CD4_naive", lower_group = "CD8_exhausted",
  group_col = "gate12b_state"
)
tnk_result <- make_umap(
  tnk, "gate12b_state", tnk_order, tnk_labels, tnk_cols,
  "T/NK state atlas", sprintf("%s cells | lineage-local UMAP", comma(nrow(tnk))),
  point_size = 0.080, label_size = 1.68
)

cairo_png <- function(filename, width, height, bg = "white", ...) {
  grDevices::png(filename, width = width, height = height, units = "in", res = 450,
                 type = "cairo", bg = bg, ...)
}

save_panel <- function(result, stem, width = 3.35, height = 3.10) {
  ggsave(file.path(panel_out, paste0(stem, ".png")), result$plot,
         width = width, height = height, units = "in", device = cairo_png,
         bg = "white", limitsize = FALSE)
  ggsave(file.path(panel_out, paste0(stem, ".pdf")), result$plot,
         width = width, height = height, units = "in", device = grDevices::cairo_pdf,
         bg = "white", limitsize = FALSE)
  invisible(NULL)
}

save_panel(global_result, "Figure1B_global_UMAP_upgraded", 3.65, 3.15)
save_panel(myeloid_result, "Figure2A_myeloid_UMAP_upgraded", 3.35, 3.10)
save_panel(tnk_result, "Figure3A_tnk_UMAP_upgraded", 3.35, 3.10)

review <- global_result$plot | myeloid_result$plot | tnk_result$plot
review <- review +
  plot_layout(widths = c(1.10, 1, 1)) +
  plot_annotation(
    title = "Publication-scale UMAP upgrade",
    subtitle = "Fixed embeddings; rigid orientation only; direct labels; no connector lines",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(family = "Arial", size = 13.0, face = "bold", hjust = 0),
      plot.subtitle = element_text(family = "Arial", size = 8.0, colour = "#59616A", hjust = 0),
      plot.tag = element_text(family = "Arial", size = 10.5, face = "bold"),
      plot.margin = margin(7, 7, 6, 7)
    )
  )
ggsave(file.path(out, "GATE12BA_UMAP_upgrade_review.png"), review,
       width = 10.4, height = 3.75, units = "in", device = cairo_png,
       bg = "white", limitsize = FALSE)
ggsave(file.path(out, "GATE12BA_UMAP_upgrade_review.pdf"), review,
       width = 10.4, height = 3.75, units = "in", device = grDevices::cairo_pdf,
       bg = "white", limitsize = FALSE)

fwrite(
  global_result$coordinates[, .(cell_id, broad_class, plot_x, plot_y)],
  file.path(source_out, "Figure1B_global_UMAP_coordinates.tsv.gz"),
  sep = "\t", compress = "gzip"
)
fwrite(
  myeloid_result$coordinates[, .(barcode, gate12b_state, plot_x, plot_y)],
  file.path(source_out, "Figure2A_myeloid_UMAP_coordinates.tsv.gz"),
  sep = "\t", compress = "gzip"
)
fwrite(
  tnk_result$coordinates[, .(barcode, gate12b_state, plot_x, plot_y)],
  file.path(source_out, "Figure3A_tnk_UMAP_coordinates.tsv.gz"),
  sep = "\t", compress = "gzip"
)
fwrite(global_result$anchors, file.path(source_out, "Figure1B_global_label_anchors.tsv"), sep = "\t")
fwrite(myeloid_result$anchors, file.path(source_out, "Figure2A_myeloid_label_anchors.tsv"), sep = "\t")
fwrite(tnk_result$anchors, file.path(source_out, "Figure3A_tnk_label_anchors.tsv"), sep = "\t")

audit <- data.table(
  panel = c("Figure1B", "Figure2A", "Figure3A"),
  cells = c(nrow(global), nrow(myeloid), nrow(tnk)),
  all_cells_rendered = TRUE,
  embedding_recomputed = FALSE,
  rigid_orientation_only = TRUE,
  connector_lines = FALSE,
  direct_labels = TRUE,
  raster_dpi = 900L,
  output_dpi = 450L,
  vector_pdf = TRUE
)
fwrite(audit, file.path(admin_out, "GATE12BA_UMAP_TECHNICAL_AUDIT.tsv"), sep = "\t")

writeLines(c(
  "# Gate12BA UMAP upgrade",
  "",
  "- The three selected embeddings are unchanged; no cluster was translated, scaled or reshaped.",
  "- Global UMAP uses the vetted K_compact_spca embedding (45 neighbours, min_dist 0.30, cosine, SPCA initialization, seed 20260811).",
  "- Myeloid and T/NK panels use the frozen Gate12AY lineage-local embeddings (35/40 neighbours, min_dist 0.25, cosine, SPCA initialization).",
  "- Only a single rigid rotation/reflection was allowed for panel fit; pairwise distances are preserved.",
  "- Every cell is shown. Rare states are drawn last; unresolved cells are visually de-emphasized.",
  "- Labels are placed near observed density modes and contain no connector segments.",
  "- UMAP is descriptive; distances between clusters and apparent cluster shapes are not interpreted biologically.",
  "- Standalone panels were rendered at their approximate final physical size, 450 dpi, with vector PDFs."
), file.path(admin_out, "GATE12BA_UMAP_UPGRADE_RECORD.md"))

cat("GATE12BA_UMAP_STATUS=COMPLETE\n")
cat("GLOBAL_CELLS=", nrow(global), "\n", sep = "")
cat("MYELOID_CELLS=", nrow(myeloid), "\n", sep = "")
cat("TNK_CELLS=", nrow(tnk), "\n", sep = "")
