#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrastr)
  library(ggrepel)
  library(RANN)
  library(grid)
})

options(stringsAsFactors = FALSE, warn = 1)
set.seed(20260811L)

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) normalizePath(args[[1L]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
selected_candidate <- if (length(args) >= 2L) args[[2L]] else "G_default_open"
out_name <- if (length(args) >= 3L) args[[3L]] else "gate12ap_selected_umap"
source_dir <- file.path(root, "results", "gate12an_umap_reembedding")
out <- file.path(root, "results", out_name)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "source_data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "admin"), recursive = TRUE, showWarnings = FALSE)

coord_file <- file.path(source_dir, paste0(selected_candidate, "_coordinates.tsv.gz"))
parameter_file <- file.path(source_dir, "candidate_parameters.tsv")
metric_file <- file.path(source_dir, "candidate_metrics.tsv")
required <- c(coord_file, parameter_file, metric_file)
if (any(!file.exists(required))) stop("Missing selected-UMAP source file(s)")

coords <- fread(coord_file)
all_parameters <- fread(parameter_file)
all_metrics <- fread(metric_file)
parameters <- all_parameters[candidate == selected_candidate]
metrics <- all_metrics[candidate == selected_candidate]
baseline_metrics <- all_metrics[candidate == "G_default_open"]
if (nrow(parameters) != 1L || nrow(metrics) != 1L) stop("Selected candidate is not uniquely defined")
if (nrow(coords) != 107886L) stop("Expected 107,886 cells; found ", nrow(coords))

class_order <- c("T_NK", "Myeloid", "B", "Progenitor", "Stromal", "Endothelial",
                 "Osteoclast", "Osteoblast", "Malignant", "Erythroid", "Unassigned")
class_labels <- c(
  T_NK = "T/NK", Myeloid = "Myeloid", B = "B cell", Progenitor = "Progenitor",
  Stromal = "Stromal", Endothelial = "Endothelial", Osteoclast = "Osteoclast",
  Osteoblast = "Osteoblast", Malignant = "Malignant", Erythroid = "Erythroid",
  Unassigned = "Unassigned"
)
class_cols <- c(
  T_NK = "#3973B9", Myeloid = "#E68435", B = "#4FA3BE",
  Progenitor = "#C8A42A", Stromal = "#35976D", Endothelial = "#6AB7AE",
  Osteoclast = "#C75B56", Osteoblast = "#8B79B9", Malignant = "#7047A3",
  Erythroid = "#D95B8A", Unassigned = "#B9BDC2"
)

## UMAP orientation is arbitrary. Apply one distance-preserving rigid transform.
## Search integer rotation angles for the orientation that gives the largest
## possible plotting scale in the target viewport, using the complete coordinate
## range (not class labels). Then reflect axes to make T/NK left of myeloid and
## B above stromal. No scaling, translation by class, or nonlinear manipulation.
raw_xy <- as.matrix(coords[, .(umap_1, umap_2)])
centre <- colMeans(raw_xy)
centered_xy <- sweep(raw_xy, 2L, centre, "-")
target_aspect <- 1.25
angle_grid <- 0:179
orientation_scores <- vapply(angle_grid, function(angle_degrees) {
  theta <- angle_degrees * pi / 180
  trial_rotation <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)),
                           nrow = 2L, byrow = TRUE)
  trial_xy <- centered_xy %*% trial_rotation
  trial_ranges <- apply(trial_xy, 2L, function(z) diff(range(z, finite = TRUE)))
  max(trial_ranges[1L] / target_aspect, trial_ranges[2L])
}, numeric(1L))
orientation_angle <- angle_grid[which.min(orientation_scores)]
theta <- orientation_angle * pi / 180
rotation <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)),
                   nrow = 2L, byrow = TRUE)
plot_xy <- centered_xy %*% rotation
colnames(plot_xy) <- c("plot_umap_1", "plot_umap_2")
coords[, `:=`(plot_umap_1 = plot_xy[, 1L], plot_umap_2 = plot_xy[, 2L])]

if (coords[broad_class == "T_NK", median(plot_umap_1)] >
    coords[broad_class == "Myeloid", median(plot_umap_1)]) {
  coords[, plot_umap_1 := -plot_umap_1]
  rotation[, 1L] <- -rotation[, 1L]
}
if (coords[broad_class == "B", median(plot_umap_2)] <
    coords[broad_class == "Stromal", median(plot_umap_2)]) {
  coords[, plot_umap_2 := -plot_umap_2]
  rotation[, 2L] <- -rotation[, 2L]
}

coords[, broad_class := factor(broad_class, levels = class_order)]
class_counts <- coords[, .N, by = broad_class]
draw_order <- c("Unassigned", as.character(class_counts[broad_class != "Unassigned"][order(-N), broad_class]))
coords[, draw_index := match(as.character(broad_class), draw_order)]
coords[, random_index := runif(.N)]
setorder(coords, draw_index, random_index)

## Choose an observed cell at the local density maximum as each label anchor.
## This avoids placing labels at a median that may lie inside a UMAP hole.
anchor_list <- lapply(class_order, function(cl) {
  dt <- coords[as.character(broad_class) == cl]
  if (!nrow(dt)) return(NULL)
  take <- min(nrow(dt), 5000L)
  idx <- if (nrow(dt) > take) sample.int(nrow(dt), take) else seq_len(nrow(dt))
  xy <- as.matrix(dt[idx, .(plot_umap_1, plot_umap_2)])
  k <- min(40L, nrow(xy))
  if (k <= 1L) {
    best <- 1L
  } else {
    nn <- RANN::nn2(xy, xy, k = k)
    best <- which.min(nn$nn.dists[, k])
  }
  data.table(
    broad_class = factor(cl, levels = class_order),
    label = unname(class_labels[cl]),
    plot_umap_1 = xy[best, 1L],
    plot_umap_2 = xy[best, 2L]
  )
})
anchors <- rbindlist(anchor_list)

alpha_values <- setNames(rep(0.68, length(class_order)), class_order)
alpha_values[c("T_NK", "Myeloid")] <- 0.52
alpha_values["Unassigned"] <- 0.16

base_umap <- ggplot(coords, aes(plot_umap_1, plot_umap_2,
                                colour = broad_class, alpha = broad_class)) +
  ggrastr::geom_point_rast(size = 0.028, raster.dpi = 800) +
  scale_colour_manual(values = class_cols, breaks = class_order,
                      labels = class_labels[class_order], drop = FALSE) +
  scale_alpha_manual(values = alpha_values, guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.018)) +
  scale_y_continuous(expand = expansion(mult = 0.018)) +
  coord_equal() +
  theme_void(base_family = "Arial", base_size = 8.5) +
  theme(
    plot.title = element_text(size = 10.2, face = "bold", hjust = 0, margin = margin(b = 2)),
    plot.subtitle = element_text(size = 7.6, colour = "#555555", hjust = 0, margin = margin(b = 2)),
    legend.position = "bottom",
    legend.text = element_text(size = 7.3),
    legend.key.height = unit(3.0, "mm"),
    legend.key.width = unit(3.0, "mm"),
    legend.spacing.x = unit(1.0, "mm"),
    plot.margin = margin(4, 4, 4, 4),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.box.background = element_rect(fill = "white", colour = NA)
  ) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE,
                               override.aes = list(size = 2.8, alpha = 1)))

subtitle <- sprintf("%s cells from 42 bone samples", format(nrow(coords), big.mark = ","))

p_clean <- base_umap +
  labs(title = "Integrated cell atlas", subtitle = subtitle, colour = NULL)

p_labelled <- base_umap +
  ggrepel::geom_label_repel(
    data = anchors,
    aes(plot_umap_1, plot_umap_2, label = label, colour = broad_class),
    inherit.aes = FALSE,
    family = "Arial", fontface = "bold", size = 2.45,
    fill = scales::alpha("white", 0.88), label.size = 0.12,
    label.padding = unit(0.12, "lines"), box.padding = 0.26,
    point.padding = 0.10, min.segment.length = 0,
    segment.color = "#777777", segment.size = 0.24,
    max.overlaps = Inf, seed = 20260811L, show.legend = FALSE
  ) +
  labs(title = "Integrated cell atlas", subtitle = subtitle, colour = NULL) +
  theme(legend.position = "none")

cairo_png <- function(filename, width, height, bg = "white", ...) {
  grDevices::png(filename, width = width, height = height, units = "in", res = 450,
                 type = "cairo", bg = bg, ...)
}

for (nm in c("clean", "labelled")) {
  p <- get(paste0("p_", nm))
  ggsave(file.path(out, paste0("Gate12AP_selected_UMAP_", nm, ".pdf")), p,
         width = 6.4, height = 5.55, units = "in", device = grDevices::cairo_pdf,
         bg = "white")
  ggsave(file.path(out, paste0("Gate12AP_selected_UMAP_", nm, ".png")), p,
         width = 6.4, height = 5.55, units = "in", device = cairo_png,
         bg = "white")
}

fwrite(coords[, .(cell_id, accession, cancer, sample_id, patient_id, compartment,
                  broad_class, broad_class_label, harmonized_state, label_source,
                  raw_umap_1 = umap_1, raw_umap_2 = umap_2,
                  plot_umap_1, plot_umap_2)],
       file.path(out, "source_data", "Gate12AP_selected_UMAP_coordinates.tsv.gz"),
       sep = "\t", compress = "gzip")
fwrite(anchors, file.path(out, "source_data", "Gate12AP_label_anchors.tsv"), sep = "\t")
fwrite(metrics, file.path(out, "source_data", "Gate12AP_selected_metrics.tsv"), sep = "\t")
fwrite(parameters, file.path(out, "source_data", "Gate12AP_selected_parameters.tsv"), sep = "\t")
fwrite(data.table(
  item = c("raw_center_umap_1", "raw_center_umap_2", "orientation_angle_degrees",
           "target_plot_aspect",
           "rotation_11", "rotation_21", "rotation_12", "rotation_22"),
  value = c(centre[1L], centre[2L], orientation_angle, target_aspect,
            rotation[1L, 1L], rotation[2L, 1L],
            rotation[1L, 2L], rotation[2L, 2L])
), file.path(out, "source_data", "Gate12AP_rigid_transform.tsv"), sep = "\t")

writeLines(c(
  "# Gate12AP selected UMAP",
  "",
  paste0("- Selected candidate: `", selected_candidate, "`."),
  sprintf("- Parameters: n_neighbors=%d; min_dist=%.2f; init=%s; repulsion=%.2f; cosine metric; 400 epochs; seed=20260811.",
          parameters$n_neighbors, parameters$min_dist, parameters$init_method,
          parameters$repulsion_strength),
  sprintf("- High-dimensional 30-neighbour overlap: %.4f.", metrics$neighbor_overlap_at_k),
  sprintf("- Default-open baseline overlap: %.4f; compact-layout delta: %.4f.",
          baseline_metrics$neighbor_overlap_at_k,
          metrics$neighbor_overlap_at_k - baseline_metrics$neighbor_overlap_at_k),
  "- Selection basis: highest neighbourhood retention among the compact SPCA/PCA sensitivity candidates, with visibly reduced inter-island whitespace.",
  sprintf("- Display coordinates use a full-range, viewport-aware rigid rotation (%d degrees before reflection); distances are unchanged.", orientation_angle),
  "- Every one of the 107,886 cells is rendered; no display downsampling or balanced-density sampling is used.",
  "- UMAP is descriptive. Inter-cluster distances and 2D shapes are not treated as inferential evidence."
), file.path(out, "admin", "GATE12AP_SELECTION_RECORD.md"))

cat("GATE12AP_STATUS=COMPLETE\n")
