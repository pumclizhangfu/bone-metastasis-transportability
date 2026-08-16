#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(igraph)
  library(matrixStats)
  library(ggplot2)
  library(patchwork)
  library(png)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) {
  stop(paste(
    "Usage: run_gate12v3_spatial_geometry.R",
    "<section_objects.rds> <full_axis1_spots.tsv.gz> <excluded_axis1_spots.tsv.gz>",
    "<histology_dir> <frozen_graph_receipt.tsv> <outdir> <n_boot> <seed> <block_width>"
  ))
}

section_file <- normalizePath(args[[1L]], mustWork = TRUE)
full_file <- normalizePath(args[[2L]], mustWork = TRUE)
excluded_file <- normalizePath(args[[3L]], mustWork = TRUE)
histology_dir <- normalizePath(args[[4L]], mustWork = TRUE)
graph_receipt_file <- normalizePath(args[[5L]], mustWork = TRUE)
outdir <- args[[6L]]
n_boot <- as.integer(args[[7L]])
seed <- as.integer(args[[8L]])
block_width <- as.numeric(args[[9L]])
if (!is.finite(n_boot) || n_boot != 1999L) stop("Gate12V3 requires exactly 1,999 valid effect bootstraps")
if (!is.finite(block_width) || block_width != 20) stop("Gate12V3 frozen spatial block width is 20 array units")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

required <- function(x, cols, label) {
  missing <- setdiff(cols, names(x))
  if (length(missing)) stop(label, " missing columns: ", paste(missing, collapse = ", "))
}

sym_knn_graph <- function(coords, k = 6L) {
  d <- as.matrix(dist(coords))
  diag(d) <- Inf
  n <- nrow(d)
  adjacency <- matrix(FALSE, nrow = n, ncol = n)
  for (i in seq_len(n)) adjacency[i, order(d[i, ])[seq_len(k)]] <- TRUE
  adjacency <- adjacency | t(adjacency)
  diag(adjacency) <- FALSE
  Matrix(adjacency * 1, sparse = TRUE)
}

distance_to_tumour <- function(adjacency, tumour) {
  graph <- graph_from_adjacency_matrix(adjacency, mode = "undirected", weighted = NULL, diag = FALSE)
  if (!any(tumour)) stop("Section has no submitted tumour-labelled spot")
  n_original <- vcount(graph)
  graph_augmented <- add_vertices(graph, 1L)
  virtual <- vcount(graph_augmented)
  graph_augmented <- add_edges(
    graph_augmented,
    as.vector(rbind(rep(virtual, sum(tumour)), which(tumour)))
  )
  distance <- as.numeric(distances(graph_augmented, v = virtual, to = seq_len(n_original))) - 1
  list(distance = distance, components = components(graph)$no)
}

distance_ring <- function(distance) {
  out <- rep("unreachable", length(distance))
  out[is.finite(distance) & distance == 0] <- "tumor"
  out[is.finite(distance) & distance == 1] <- "boundary"
  out[is.finite(distance) & distance >= 2 & distance <= 3] <- "proximal_2_3"
  out[is.finite(distance) & distance >= 4 & distance <= 6] <- "intermediate_4_6"
  out[is.finite(distance) & distance >= 7] <- "distal_7plus"
  factor(out, levels = c("tumor", "boundary", "proximal_2_3", "intermediate_4_6",
                         "distal_7plus", "unreachable"))
}

row_median_na <- function(x) {
  out <- matrixStats::rowMedians(x, na.rm = TRUE)
  out[!is.finite(out)] <- NA_real_
  out
}

bootstrap_group_summary <- function(block_table, group_cols, value_col, block_ids, boot_index) {
  groups <- unique(block_table[, ..group_cols])
  rows <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    selector <- rep(TRUE, nrow(block_table))
    for (column in group_cols) selector <- selector & block_table[[column]] == groups[[column]][[i]]
    z <- block_table[selector]
    values <- setNames(rep(NA_real_, length(block_ids)), block_ids)
    values[as.character(z$spatial_block)] <- z[[value_col]]
    boot_matrix <- matrix(values[boot_index], nrow = nrow(boot_index), ncol = ncol(boot_index))
    boot_stat <- row_median_na(boot_matrix)
    valid <- boot_stat[is.finite(boot_stat)]
    row <- copy(groups[i])
    row[, `:=`(
      estimate = median(z[[value_col]], na.rm = TRUE),
      ci_low = if (length(valid)) quantile(valid, 0.025) else NA_real_,
      ci_high = if (length(valid)) quantile(valid, 0.975) else NA_real_,
      occupied_blocks = sum(is.finite(values)),
      valid_bootstraps = length(valid),
      interval_type = "spatial_block_bootstrap"
    )]
    rows[[i]] <- row
  }
  rbindlist(rows, fill = TRUE)
}

bootstrap_effect <- function(block_table, layer_value, ring_a, ring_b, block_ids,
                             n_valid = 1999L, seed_value = 1L) {
  z <- block_table[layer == layer_value & ring %chin% c(ring_a, ring_b)]
  va <- setNames(rep(NA_real_, length(block_ids)), block_ids)
  vb <- setNames(rep(NA_real_, length(block_ids)), block_ids)
  za <- z[ring == ring_a]
  zb <- z[ring == ring_b]
  va[as.character(za$spatial_block)] <- za$value
  vb[as.character(zb$spatial_block)] <- zb$value
  point <- median(va, na.rm = TRUE) - median(vb, na.rm = TRUE)
  set.seed(seed_value)
  accepted <- numeric()
  attempts <- 0L
  max_attempts <- n_valid * 20L
  batch <- max(500L, n_valid)
  while (length(accepted) < n_valid && attempts < max_attempts) {
    draws <- matrix(sample.int(length(block_ids), batch * length(block_ids), replace = TRUE),
                    nrow = batch, ncol = length(block_ids))
    ea <- row_median_na(matrix(va[draws], nrow = batch, ncol = length(block_ids)))
    eb <- row_median_na(matrix(vb[draws], nrow = batch, ncol = length(block_ids)))
    effect <- ea - eb
    effect <- effect[is.finite(effect)]
    accepted <- c(accepted, effect)
    attempts <- attempts + batch
  }
  if (length(accepted) < n_valid) stop("Could not obtain required valid block-bootstrap effects")
  accepted <- accepted[seq_len(n_valid)]
  data.table(
    estimate = point,
    ci_low = quantile(accepted, 0.025),
    ci_high = quantile(accepted, 0.975),
    positive_fraction = mean(accepted > 0),
    negative_fraction = mean(accepted < 0),
    valid_bootstraps = n_valid,
    bootstrap_draws_attempted = attempts,
    ring_a_blocks = sum(is.finite(va)),
    ring_b_blocks = sum(is.finite(vb))
  )
}

sections <- readRDS(section_file)
full_spots <- fread(full_file)
excluded_spots <- fread(excluded_file)
frozen_graph <- fread(graph_receipt_file)
required(full_spots, c("sample", "barcode", "x", "y", "submitted_label", "axis1_weighted_score",
                       "rctd_axis1_proxy"), "Full Axis1 spots")
required(excluded_spots, c("sample", "barcode", "axis1_full_recalculated", "axis1_malignant_excluded",
                           "malignant_excluded_rctd_proxy"), "Malignant-excluded spots")
required(frozen_graph, c("sample", "spots", "k", "symmetric_edges", "min_degree", "median_degree",
                         "max_degree"), "Frozen graph receipt")
if (length(sections) != 4L) stop("Expected four spatial sections")
if (nrow(full_spots) != 8190L || nrow(excluded_spots) != 8190L) stop("Expected 8,190 spots in both score layers")
setkey(full_spots, sample, barcode)
setkey(excluded_spots, sample, barcode)
if (!identical(full_spots[, .(sample, barcode)], excluded_spots[, .(sample, barcode)])) {
  stop("Full and malignant-excluded spot keys do not match")
}

sample_manifest <- data.table(
  sample = c("GSM9564255", "GSM9564256", "GSM9564257", "GSM9564258"),
  library = c("24664", "24665", "24666", "24667")
)
rctd_classes <- c("B cell", "Endothelial", "Erythroid", "Macrophage", "MSCs", "Myeloid",
                  "plasma B cell", "pre-B cell", "T cell")

spot_rows <- list()
neighbor_rows <- list()
graph_rows <- list()
curve_rows <- list()
ring_rows <- list()
effect_rows <- list()
neighborhood_summary_rows <- list()
neighborhood_effect_rows <- list()
trend_rows <- list()
section_support_rows <- list()
block_manifest_rows <- list()

for (section_index in seq_len(nrow(sample_manifest))) {
  sample_value <- sample_manifest$sample[[section_index]]
  section <- sections[[sample_value]]
  barcodes <- rownames(section$weights)
  if (!identical(colnames(section$weights), rctd_classes)) stop(sample_value, ": unexpected RCTD class order")
  coords <- as.matrix(section$coords[barcodes, c("x", "y"), drop = FALSE])
  adjacency <- sym_knn_graph(coords, k = 6L)
  degree <- as.numeric(rowSums(adjacency))
  graph_distance <- distance_to_tumour(
    adjacency,
    as.character(section$submitted_label[barcodes]) == "Tumor"
  )
  distance <- graph_distance$distance
  ring <- distance_ring(distance)
  spatial_block <- paste(floor(coords[, "x"] / block_width), floor(coords[, "y"] / block_width), sep = "_")
  block_ids <- sort(unique(spatial_block))
  set.seed(seed + section_index * 1000L)
  boot_index <- matrix(sample.int(length(block_ids), n_boot * length(block_ids), replace = TRUE),
                       nrow = n_boot, ncol = length(block_ids))

  full <- full_spots[.(sample_value, barcodes)]
  excluded <- excluded_spots[.(sample_value, barcodes)]
  if (anyNA(full$barcode) || anyNA(excluded$barcode)) stop(sample_value, ": score layer alignment failed")
  if (max(abs(full$axis1_weighted_score - excluded$axis1_full_recalculated)) > 1e-12) {
    stop(sample_value, ": stored full Axis1 layers do not agree")
  }
  proxy_difference <- max(abs(full$rctd_axis1_proxy - excluded$malignant_excluded_rctd_proxy))
  if (!is.finite(proxy_difference) || proxy_difference > 1e-12) stop(sample_value, ": RCTD proxy layers differ")

  weights <- section$weights[barcodes, rctd_classes, drop = FALSE]
  neighbor_weights <- as.matrix(adjacency %*% weights) / degree
  colnames(neighbor_weights) <- rctd_classes
  rownames(neighbor_weights) <- barcodes
  neighbor_proxy_full <- as.numeric(adjacency %*% full$rctd_axis1_proxy) / degree
  neighbor_proxy_excluded <- as.numeric(adjacency %*% excluded$malignant_excluded_rctd_proxy) / degree

  positions <- section$positions[barcodes, , drop = FALSE]
  section_spots <- data.table(
    sample = sample_value,
    barcode = barcodes,
    array_col = coords[, "x"],
    array_row = coords[, "y"],
    pixel_row = positions$pixel_row,
    pixel_col = positions$pixel_col,
    submitted_label = as.character(section$submitted_label[barcodes]),
    graph_distance = distance,
    distance_finite = is.finite(distance),
    distance_ring = ring,
    spatial_block = spatial_block,
    full_axis1 = full$axis1_weighted_score,
    malignant_excluded_axis1 = excluded$axis1_malignant_excluded,
    rctd_proxy_full = full$rctd_axis1_proxy,
    rctd_proxy_malignant_excluded = excluded$malignant_excluded_rctd_proxy,
    one_ring_rctd_proxy_full = neighbor_proxy_full,
    one_ring_rctd_proxy_malignant_excluded = neighbor_proxy_excluded
  )
  spot_rows[[section_index]] <- section_spots

  neighbor_long <- as.data.table(neighbor_weights, keep.rownames = "barcode")
  neighbor_long[, `:=`(
    sample = sample_value,
    graph_distance = distance,
    distance_ring = ring,
    spatial_block = spatial_block
  )]
  neighbor_long <- melt(
    neighbor_long,
    id.vars = c("sample", "barcode", "graph_distance", "distance_ring", "spatial_block"),
    measure.vars = rctd_classes,
    variable.name = "rctd_class",
    value.name = "one_ring_weight"
  )
  neighbor_rows[[section_index]] <- neighbor_long

  graph_row <- data.table(
    sample = sample_value,
    spots = nrow(adjacency),
    k = 6L,
    symmetric_edges = sum(adjacency) / 2,
    min_degree = min(degree),
    median_degree = median(degree),
    max_degree = max(degree),
    graph_components = graph_distance$components,
    tumour_spots = sum(ring == "tumor"),
    boundary_spots = sum(ring == "boundary"),
    unreachable_spots = sum(ring == "unreachable"),
    max_finite_distance = max(distance[is.finite(distance)]),
    occupied_spatial_blocks = length(block_ids),
    boundary_blocks = uniqueN(spatial_block[ring == "boundary"]),
    finite_distal_blocks = uniqueN(spatial_block[ring == "distal_7plus"]),
    proxy_max_absolute_difference = proxy_difference
  )
  frozen_row <- frozen_graph[sample == sample_value]
  graph_row[, frozen_receipt_match :=
              spots == frozen_row$spots & k == frozen_row$k &
              symmetric_edges == frozen_row$symmetric_edges &
              min_degree == frozen_row$min_degree & median_degree == frozen_row$median_degree &
              max_degree == frozen_row$max_degree]
  graph_rows[[section_index]] <- graph_row
  block_manifest_rows[[section_index]] <- section_spots[, .(
    spots = .N,
    finite_distance_spots = sum(distance_finite),
    tumour_spots = sum(distance_ring == "tumor"),
    boundary_spots = sum(distance_ring == "boundary"),
    finite_distal_spots = sum(distance_ring == "distal_7plus"),
    unreachable_spots = sum(distance_ring == "unreachable")
  ), by = .(sample, spatial_block)]

  score_long <- melt(
    section_spots,
    id.vars = c("sample", "barcode", "graph_distance", "distance_finite", "distance_ring", "spatial_block"),
    measure.vars = c("full_axis1", "malignant_excluded_axis1", "one_ring_rctd_proxy_full",
                     "one_ring_rctd_proxy_malignant_excluded"),
    variable.name = "layer",
    value.name = "score"
  )
  score_long[, ring := as.character(distance_ring)]
  block_ring <- score_long[, .(value = median(score)), by = .(spatial_block, ring, layer)]
  block_exact <- score_long[distance_finite == TRUE,
                            .(value = median(score)), by = .(spatial_block, graph_distance, layer)]

  ring_summary <- bootstrap_group_summary(
    block_ring, c("ring", "layer"), "value", block_ids, boot_index
  )
  ring_summary[, sample := sample_value]
  ring_rows[[section_index]] <- ring_summary
  curve_summary <- bootstrap_group_summary(
    block_exact, c("graph_distance", "layer"), "value", block_ids, boot_index
  )
  curve_summary[, sample := sample_value]
  curve_rows[[section_index]] <- curve_summary

  contrast_definitions <- list(
      boundary_minus_distal = c("boundary", "distal_7plus"),
      tumor_minus_distal = c("tumor", "distal_7plus")
    )
  for (layer_value in unique(block_ring$layer)) {
    for (comparison_name in names(contrast_definitions)) {
      rings <- contrast_definitions[[comparison_name]]
      result <- bootstrap_effect(
        block_ring, layer_value, rings[[1L]], rings[[2L]], block_ids,
        n_valid = n_boot,
        seed_value = seed + section_index * 10000L +
          match(layer_value, unique(block_ring$layer)) * 100L +
          ifelse(comparison_name == "boundary_minus_distal", 1L, 2L)
      )
      result[, `:=`(
        sample = sample_value,
        layer = layer_value,
        comparison = comparison_name,
        ring_a = rings[[1L]],
        ring_b = rings[[2L]]
      )]
      effect_rows[[length(effect_rows) + 1L]] <- result
    }
  }

  finite_curve <- curve_summary[is.finite(graph_distance)]
  trend_rows[[section_index]] <- finite_curve[, .(
    n_distances = .N,
    spearman_rho_distance = suppressWarnings(cor(graph_distance, estimate, method = "spearman")),
    first_distance = min(graph_distance),
    last_distance = max(graph_distance)
  ), by = .(sample, layer)]

  neighbor_block_ring <- neighbor_long[, .(value = median(one_ring_weight)),
                                        by = .(spatial_block, ring = as.character(distance_ring), rctd_class)]
  neighborhood_summary <- bootstrap_group_summary(
    neighbor_block_ring, c("ring", "rctd_class"), "value", block_ids, boot_index
  )
  neighborhood_summary[, sample := sample_value]
  neighborhood_summary_rows[[section_index]] <- neighborhood_summary
  for (class_value in rctd_classes) {
    z <- copy(neighbor_block_ring[rctd_class == class_value])
    z[, layer := class_value]
    result <- bootstrap_effect(
      z, class_value, "boundary", "distal_7plus", block_ids,
      n_valid = n_boot,
      seed_value = seed + section_index * 20000L + match(class_value, rctd_classes)
    )
    result[, `:=`(
      sample = sample_value,
      rctd_class = class_value,
      comparison = "boundary_minus_distal"
    )]
    neighborhood_effect_rows[[length(neighborhood_effect_rows) + 1L]] <- result
  }

  primary <- rbindlist(effect_rows)[sample == sample_value & comparison == "boundary_minus_distal"]
  primary_wide <- dcast(primary, sample ~ layer, value.var = c("estimate", "positive_fraction"))
  evaluable <- graph_row$boundary_blocks >= 3L & graph_row$finite_distal_blocks >= 5L
  supports <- evaluable &&
    primary_wide$estimate_full_axis1 > 0 & primary_wide$positive_fraction_full_axis1 >= 0.75 &
    primary_wide$estimate_malignant_excluded_axis1 > 0 &
      primary_wide$positive_fraction_malignant_excluded_axis1 >= 0.75 &
    primary_wide$estimate_one_ring_rctd_proxy_full > 0 &
      primary_wide$positive_fraction_one_ring_rctd_proxy_full >= 0.75 &
    primary_wide$estimate_one_ring_rctd_proxy_malignant_excluded > 0 &
      primary_wide$positive_fraction_one_ring_rctd_proxy_malignant_excluded >= 0.75
  section_support_rows[[section_index]] <- data.table(
    sample = sample_value,
    decision_evaluable = evaluable,
    boundary_blocks = graph_row$boundary_blocks,
    finite_distal_blocks = graph_row$finite_distal_blocks,
    full_axis1_effect = primary_wide$estimate_full_axis1,
    full_axis1_positive_fraction = primary_wide$positive_fraction_full_axis1,
    malignant_excluded_effect = primary_wide$estimate_malignant_excluded_axis1,
    malignant_excluded_positive_fraction = primary_wide$positive_fraction_malignant_excluded_axis1,
    one_ring_rctd_proxy_effect = primary_wide$estimate_one_ring_rctd_proxy_full,
    one_ring_rctd_proxy_positive_fraction = primary_wide$positive_fraction_one_ring_rctd_proxy_full,
    proxy_layers_identical = proxy_difference <= 1e-12,
    supports_tumour_proximal_direction = supports,
    non_evaluable_reason = if (evaluable) "" else "fewer_than_3_boundary_blocks_or_5_distal_blocks"
  )
}

spots <- rbindlist(spot_rows)
neighbors <- rbindlist(neighbor_rows)
graph_receipt <- rbindlist(graph_rows)
block_manifest <- rbindlist(block_manifest_rows)
distance_curves <- rbindlist(curve_rows)
distance_ring_summary <- rbindlist(ring_rows)
distance_effects <- rbindlist(effect_rows)
neighborhood_summary <- rbindlist(neighborhood_summary_rows)
neighborhood_effects <- rbindlist(neighborhood_effect_rows)
distance_trends <- rbindlist(trend_rows)
section_support <- rbindlist(section_support_rows)

if (!all(graph_receipt$frozen_receipt_match)) stop("Reconstructed six-neighbor graph does not match frozen receipt")
if (nrow(spots) != 8190L || spots[, anyDuplicated(barcode), by = sample][, any(V1 > 0)]) {
  stop("Spot output scope mismatch")
}
if (!all(section_support$proxy_layers_identical)) stop("Full and malignant-excluded RCTD proxy layers differ")

supporting_sections <- sum(section_support$supports_tumour_proximal_direction)
decision_status <- if (supporting_sections >= 3L) "CONSISTENT_ORGANIZATION" else "HETEROGENEOUS_ORGANIZATION"
decision <- data.table(
  gate = "Gate12V3",
  execution_status = "COMPLETE",
  scientific_decision = decision_status,
  measured_spots = nrow(spots),
  sections = uniqueN(spots$sample),
  decision_evaluable_sections = sum(section_support$decision_evaluable),
  supporting_sections = supporting_sections,
  unreachable_spots = sum(!spots$distance_finite),
  all_graph_receipts_match = all(graph_receipt$frozen_receipt_match),
  all_proxy_layers_identical = all(section_support$proxy_layers_identical),
  independent_animal_inference = FALSE,
  interpretation_scope = "linked_cross_species_descriptive"
)

fwrite(spots, file.path(outdir, "spot_spatial_geometry.tsv.gz"), sep = "\t")
fwrite(neighbors, file.path(outdir, "one_ring_neighborhood_composition.tsv.gz"), sep = "\t")
fwrite(graph_receipt, file.path(outdir, "spatial_graph_distance_receipt.tsv"), sep = "\t")
fwrite(block_manifest, file.path(outdir, "spatial_block_manifest.tsv"), sep = "\t")
fwrite(distance_curves, file.path(outdir, "axis1_graph_distance_curves.tsv"), sep = "\t")
fwrite(distance_ring_summary, file.path(outdir, "axis1_distance_ring_summary.tsv"), sep = "\t")
fwrite(distance_effects, file.path(outdir, "axis1_distance_effects.tsv"), sep = "\t")
fwrite(distance_trends, file.path(outdir, "axis1_distance_trends.tsv"), sep = "\t")
fwrite(neighborhood_summary, file.path(outdir, "neighborhood_distance_summary.tsv"), sep = "\t")
fwrite(neighborhood_effects, file.path(outdir, "neighborhood_class_effects.tsv"), sep = "\t")
fwrite(section_support, file.path(outdir, "section_organization_support.tsv"), sep = "\t")
fwrite(decision, file.path(outdir, "gate12v3_decision.tsv"), sep = "\t")

input_audit <- data.table(
  metric = c(
    "measured_spots", "sections", "full_score_spots", "malignant_excluded_score_spots",
    "rctd_classes", "knn_k", "spatial_block_width", "bootstrap_iterations",
    "submitted_tumour_spots", "unreachable_spots", "decision_evaluable_sections", "seed"
  ),
  value = as.character(c(
    nrow(spots), uniqueN(spots$sample), nrow(full_spots), nrow(excluded_spots),
    length(rctd_classes), 6L, block_width, n_boot,
    sum(spots$submitted_label == "Tumor"), sum(!spots$distance_finite),
    sum(section_support$decision_evaluable), seed
  ))
)
fwrite(input_audit, file.path(outdir, "input_audit.tsv"), sep = "\t")

# Native H&E overlays.
lighten_histology <- function(image, amount = 0.30) {
  image[, , seq_len(3L)] <- image[, , seq_len(3L)] * (1 - amount) + amount
  image
}
robust_limit <- function(x, probability = 0.98) unname(quantile(abs(x), probability, na.rm = TRUE))
assets <- list()
plot_data <- list()
for (i in seq_len(nrow(sample_manifest))) {
  sample_value <- sample_manifest$sample[[i]]
  library <- sample_manifest$library[[i]]
  prefix <- paste0(sample_value, "_Bone_ST_", library)
  image_path <- file.path(histology_dir, paste0(prefix, "_tissue_lowres_image.png"))
  scale_path <- file.path(histology_dir, paste0(prefix, "_scalefactors_json.json"))
  if (!file.exists(image_path) || !file.exists(scale_path)) stop(sample_value, ": missing H&E assets")
  image <- readPNG(image_path)
  scale <- fromJSON(scale_path)$tissue_lowres_scalef
  data <- copy(spots[sample == sample_value])
  data[, `:=`(
    plot_x = pixel_col * scale,
    plot_y = dim(image)[1L] - pixel_row * scale
  )]
  if (min(data$plot_x) < 0 || max(data$plot_x) > dim(image)[2L] ||
      min(data$plot_y) < 0 || max(data$plot_y) > dim(image)[1L]) {
    stop(sample_value, ": spot coordinate outside H&E image")
  }
  assets[[sample_value]] <- list(
    image = image,
    light = lighten_histology(image),
    width = dim(image)[2L],
    height = dim(image)[1L],
    scale = scale
  )
  plot_data[[sample_value]] <- data
}

score_limit_full <- robust_limit(spots$full_axis1)
score_limit_excluded <- robust_limit(spots$malignant_excluded_axis1)
ring_colours <- c(
  tumor = "#D73027", boundary = "#FC8D59", proximal_2_3 = "#FEE08B",
  intermediate_4_6 = "#91CF60", distal_7plus = "#4575B4", unreachable = "#7F7F7F"
)
theme_spatial <- theme_void(base_size = 8.5) +
  theme(plot.title = element_text(face = "bold", size = 8.8, hjust = 0.5),
        plot.margin = margin(1, 1, 1, 1), legend.position = "none")

spatial_panel <- function(sample_value, layer, title) {
  asset <- assets[[sample_value]]
  data <- copy(plot_data[[sample_value]])
  background <- if (layer == "rings") asset$image else asset$light
  p <- ggplot() +
    annotation_raster(background, xmin = 0, xmax = asset$width, ymin = 0, ymax = asset$height) +
    coord_fixed(xlim = c(0, asset$width), ylim = c(0, asset$height), expand = FALSE) +
    theme_spatial + labs(title = title)
  if (layer == "rings") {
    p <- p +
      geom_point(data = data, aes(plot_x, plot_y, fill = distance_ring), shape = 21,
                 colour = "white", stroke = 0.08, size = 0.95, alpha = 0.78) +
      scale_fill_manual(values = ring_colours, drop = FALSE)
  } else {
    value_col <- if (layer == "full") "full_axis1" else "malignant_excluded_axis1"
    limit <- if (layer == "full") score_limit_full else score_limit_excluded
    data[, display_value := pmax(pmin(get(value_col), limit), -limit)]
    p <- p +
      geom_point(data = data, aes(plot_x, plot_y, fill = display_value), shape = 21,
                 colour = "#252525", stroke = 0.06, size = 1.05, alpha = 0.76) +
      geom_point(data = data[submitted_label == "Tumor"], aes(plot_x, plot_y), shape = 21,
                 fill = NA, colour = "black", stroke = 0.18, size = 1.25) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                           limits = c(-limit, limit))
  }
  p
}

spatial_row <- function(layer, row_label) {
  panels <- lapply(seq_len(nrow(sample_manifest)), function(i) {
    sample_value <- sample_manifest$sample[[i]]
    spatial_panel(sample_value, layer, paste0(if (i == 1L) paste0(row_label, "  ") else "", sample_value))
  })
  wrap_plots(panels, nrow = 1)
}

p_rings <- spatial_row("rings", "A")
p_full <- spatial_row("full", "B")
p_excluded <- spatial_row("excluded", "C")

curve_plot_data <- distance_curves[layer %chin% c("full_axis1", "malignant_excluded_axis1")]
curve_plot_data[, layer_label := factor(layer,
  levels = c("full_axis1", "malignant_excluded_axis1"),
  labels = c("Full Axis1", "Malignant-excluded Axis1"))]
p_curves <- ggplot(curve_plot_data, aes(graph_distance, estimate, colour = layer_label, fill = layer_label)) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.13, colour = NA) +
  geom_line(linewidth = 0.55) + geom_point(size = 0.65) +
  facet_wrap(~sample, scales = "free_x", nrow = 1) +
  theme_bw(base_size = 8.3) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom") +
  labs(title = "D  Axis1 versus finite graph distance from submitted tumour labels",
       subtitle = "Curves are medians of spatial-block medians; bands are 95% block-bootstrap intervals",
       x = "Six-neighbor graph distance", y = "Spatial-block median score", colour = NULL, fill = NULL)

primary_effect_plot <- distance_effects[
  comparison == "boundary_minus_distal" &
    layer %chin% c("full_axis1", "malignant_excluded_axis1", "one_ring_rctd_proxy_full")
]
primary_effect_plot[, layer_label := factor(layer,
  levels = c("full_axis1", "malignant_excluded_axis1", "one_ring_rctd_proxy_full"),
  labels = c("Full Axis1", "Malignant-excluded Axis1", "One-ring RCTD proxy"))]
p_effect <- ggplot(primary_effect_plot, aes(estimate, sample, colour = layer_label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.14,
                 position = position_dodge(width = 0.45)) +
  geom_point(position = position_dodge(width = 0.45), size = 1.8) +
  theme_bw(base_size = 8.3) + theme(panel.grid.minor = element_blank(), legend.position = "bottom") +
  labs(title = "E  Boundary-minus-finite-distal effects",
       subtitle = "GSM9564257 is displayed but not decision-evaluable (one boundary block)",
       x = "Difference of spatial-block medians (95% block-bootstrap CI)", y = NULL, colour = NULL)

figure <- p_rings / p_full / p_excluded / (p_curves | p_effect) +
  plot_layout(heights = c(1, 1, 1, 1.15)) +
  plot_annotation(
    title = "Frozen Axis1 is spatially organized along submitted tumour-boundary distance",
    subtitle = paste0(
      "All 8,190 measured spots are shown on native H&E; unreachable graph components are grey. Decision: ",
      decision_status, ". Sections are not asserted to be independent animals."
    ),
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 9, colour = "#333333"))
  )

png_tmp <- tempfile(pattern = "gate12v3_histology_", fileext = ".png")
ggsave(png_tmp, figure, width = 18, height = 14.8, dpi = 300, device = ragg::agg_png, bg = "white")
if (!file.copy(png_tmp, file.path(outdir, "Figure_gate12v3_native_histology_geometry.png"), overwrite = TRUE)) {
  stop("Failed to copy Gate12V3 PNG")
}
unlink(png_tmp)
ggsave(file.path(outdir, "Figure_gate12v3_native_histology_geometry.pdf"), figure,
       width = 18, height = 14.8, device = cairo_pdf, bg = "white")

# Complete nine-class neighborhood display.
neighborhood_plot <- copy(neighborhood_summary)
neighborhood_plot[, ring := factor(ring, levels = levels(spots$distance_ring))]
neighborhood_plot[, rctd_class := factor(rctd_class, levels = rev(rctd_classes))]
p_neighbor_heat <- ggplot(neighborhood_plot, aes(ring, rctd_class, fill = estimate)) +
  geom_tile(colour = "white", linewidth = 0.15) +
  facet_wrap(~sample, nrow = 2) +
  scale_fill_viridis_c(option = "C", limits = c(0, max(neighborhood_plot$estimate, na.rm = TRUE))) +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(title = "A  One-ring RCTD composition across graph-distance rings",
       x = NULL, y = NULL, fill = "Median neighbor weight")

neighborhood_effects[, rctd_class := factor(rctd_class, levels = rev(rctd_classes))]
p_neighbor_effect <- ggplot(neighborhood_effects, aes(estimate, rctd_class)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.14) +
  geom_point(size = 1.5, colour = "#0072B2") +
  facet_wrap(~sample, nrow = 2, scales = "free_x") +
  theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank()) +
  labs(title = "B  Class-specific boundary-minus-finite-distal neighborhood effects",
       x = "Difference in spatial-block median neighbor weight (95% CI)", y = NULL)

neighborhood_figure <- p_neighbor_heat / p_neighbor_effect +
  plot_annotation(
    title = "Nine-class neighborhood ecology without cell-class selection",
    subtitle = "Every submitted RCTD class is retained; intervals use within-section spatial-block resampling."
  )
png_tmp2 <- tempfile(pattern = "gate12v3_neighborhood_", fileext = ".png")
ggsave(png_tmp2, neighborhood_figure, width = 14.5, height = 13.5, dpi = 300,
       device = ragg::agg_png, bg = "white")
if (!file.copy(png_tmp2, file.path(outdir, "FigureS_gate12v3_neighborhood_ecology.png"), overwrite = TRUE)) {
  stop("Failed to copy Gate12V3 neighborhood PNG")
}
unlink(png_tmp2)
ggsave(file.path(outdir, "FigureS_gate12v3_neighborhood_ecology.pdf"), neighborhood_figure,
       width = 14.5, height = 13.5, device = cairo_pdf, bg = "white")

support_lines <- apply(section_support, 1L, function(row) {
  paste0(
    "- ", row[["sample"]], ": evaluable=", row[["decision_evaluable"]],
    ", full effect=", sprintf("%.4f", as.numeric(row[["full_axis1_effect"]])),
    " (positive bootstrap ", sprintf("%.1f%%", 100 * as.numeric(row[["full_axis1_positive_fraction"]])), ")",
    ", excluded effect=", sprintf("%.4f", as.numeric(row[["malignant_excluded_effect"]])),
    " (", sprintf("%.1f%%", 100 * as.numeric(row[["malignant_excluded_positive_fraction"]])), ")",
    ", neighbor proxy effect=", sprintf("%.4f", as.numeric(row[["one_ring_rctd_proxy_effect"]])),
    " (", sprintf("%.1f%%", 100 * as.numeric(row[["one_ring_rctd_proxy_positive_fraction"]])), ")",
    ", supports=", row[["supports_tumour_proximal_direction"]], "."
  )
})
writeLines(c(
  "# Gate12V3 spatial geometry and neighborhood ecology checkpoint", "",
  "## Execution status", "",
  "- Computational run: COMPLETE (independent audit pending).",
  paste0("- Scientific decision: **", decision_status, "**."),
  paste0("- Measured spots: ", nrow(spots), " across four sections."),
  paste0("- Spots in graph components without a submitted tumour label: ", sum(!spots$distance_finite),
         "; retained as `unreachable` and excluded only from ordered distance effects."),
  paste0("- Decision-evaluable sections: ", sum(section_support$decision_evaluable), "/4; supporting sections: ",
         supporting_sections, "/4."),
  "", "## Section-specific decision effects", "", support_lines,
  "", "## Interpretation boundary", "",
  "All effects are section-specific spatial-block summaries. The four sections are not asserted to be four independent animals. Submitted tumour labels are source annotations rather than new pathologist polygons. Full and malignant-excluded RCTD proxy layers are identical because the nine-class RCTD reference contains no malignant class; this is a linked sensitivity, not independent validation."
), file.path(outdir, "GATE12V3_CHECKPOINT.md"))

saveRDS(list(
  decision = decision,
  section_support = section_support,
  graph_receipt = graph_receipt,
  seed = seed,
  n_boot = n_boot,
  block_width = block_width,
  ring_levels = levels(spots$distance_ring),
  independent_animal_inference = FALSE
), file.path(outdir, "gate12v3_results.rds"), compress = "xz")
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
cat("GATE12V3_EXECUTION_STATUS=COMPLETE\n")
print(decision)
