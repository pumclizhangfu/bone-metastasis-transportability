#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(igraph)
  library(matrixStats)
  library(ggplot2)
  library(patchwork)
})

if (!requireNamespace("digest", quietly = TRUE)) stop("Package digest is required")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: run_gate12ab_spatial_sensitivity_grid.R <project_root> <outdir> <n_boot> <seed>")
}
root <- normalizePath(args[[1L]], mustWork = TRUE)
outdir <- args[[2L]]
n_boot <- as.integer(args[[3L]])
seed <- as.integer(args[[4L]])
if (!identical(n_boot, 999L)) stop("Gate12AB frozen grid requires exactly 999 valid bootstrap effects")
if (!is.finite(seed)) stop("Seed must be finite")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

revision_root <- normalizePath(file.path(outdir, "..", ".."), mustWork = TRUE)
figure_dir <- file.path(revision_root, "figures", "supplementary")
source_dir <- file.path(revision_root, "source_data")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

spot_file <- file.path(root, "results", "gate12v3_spatial_geometry", "spot_spatial_geometry.tsv.gz")
primary_effect_file <- file.path(root, "results", "gate12v3_spatial_geometry", "axis1_distance_effects.tsv")
spots <- fread(spot_file)
required_columns <- c("sample", "barcode", "array_col", "array_row", "submitted_label",
                      "full_axis1", "malignant_excluded_axis1", "rctd_proxy_full")
missing_columns <- setdiff(required_columns, names(spots))
if (length(missing_columns)) stop("Missing spot columns: ", paste(missing_columns, collapse = ", "))
if (nrow(spots) != 8190L || uniqueN(spots$sample) != 4L) stop("Expected 8,190 spots in four sections")
if (spots[, anyDuplicated(barcode), by = sample][, any(V1 > 0L)]) stop("Duplicated spot barcode within section")

k_grid <- c(4L, 6L, 8L)
block_width_grid <- c(15, 20, 25)

sym_knn_graph <- function(coords, k) {
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
  if (!any(tumour)) stop("Section has no submitted tumor-labelled spot")
  n_original <- vcount(graph)
  graph_augmented <- add_vertices(graph, 1L)
  virtual <- vcount(graph_augmented)
  graph_augmented <- add_edges(graph_augmented,
                               as.vector(rbind(rep(virtual, sum(tumour)), which(tumour))))
  as.numeric(distances(graph_augmented, v = virtual, to = seq_len(n_original))) - 1
}

distance_ring <- function(distance) {
  out <- rep("unreachable", length(distance))
  out[is.finite(distance) & distance == 0] <- "tumor"
  out[is.finite(distance) & distance == 1] <- "boundary"
  out[is.finite(distance) & distance >= 2 & distance <= 3] <- "proximal_2_3"
  out[is.finite(distance) & distance >= 4 & distance <= 6] <- "intermediate_4_6"
  out[is.finite(distance) & distance >= 7] <- "distal_7plus"
  out
}

row_median_na <- function(x) {
  out <- matrixStats::rowMedians(x, na.rm = TRUE)
  out[!is.finite(out)] <- NA_real_
  out
}

bootstrap_effect <- function(block_table, block_ids, n_valid, seed_value) {
  va <- setNames(rep(NA_real_, length(block_ids)), block_ids)
  vb <- setNames(rep(NA_real_, length(block_ids)), block_ids)
  a <- block_table[ring == "boundary"]
  b <- block_table[ring == "distal_7plus"]
  va[as.character(a$spatial_block)] <- a$value
  vb[as.character(b$spatial_block)] <- b$value
  point <- median(va, na.rm = TRUE) - median(vb, na.rm = TRUE)
  set.seed(seed_value)
  accepted <- numeric()
  attempts <- 0L
  batch <- max(500L, n_valid)
  while (length(accepted) < n_valid && attempts < n_valid * 20L) {
    draws <- matrix(sample.int(length(block_ids), batch * length(block_ids), replace = TRUE),
                    nrow = batch, ncol = length(block_ids))
    ea <- row_median_na(matrix(va[draws], nrow = batch, ncol = length(block_ids)))
    eb <- row_median_na(matrix(vb[draws], nrow = batch, ncol = length(block_ids)))
    effect <- ea - eb
    accepted <- c(accepted, effect[is.finite(effect)])
    attempts <- attempts + batch
  }
  if (length(accepted) < n_valid) stop("Could not obtain required valid block-bootstrap effects")
  accepted <- accepted[seq_len(n_valid)]
  data.table(
    estimate = point,
    ci_low = unname(quantile(accepted, 0.025)),
    ci_high = unname(quantile(accepted, 0.975)),
    positive_fraction = mean(accepted > 0),
    valid_bootstraps = n_valid,
    bootstrap_draws_attempted = attempts,
    boundary_blocks = sum(is.finite(va)),
    finite_distal_blocks = sum(is.finite(vb))
  )
}

effect_rows <- list()
section_rows <- list()
graph_rows <- list()
row_index <- 0L
section_index <- 0L

for (sample_value in sort(unique(spots$sample))) {
  section_index <- section_index + 1L
  z <- copy(spots[sample == sample_value])
  coords <- as.matrix(z[, .(array_col, array_row)])
  for (k_value in k_grid) {
    adjacency <- sym_knn_graph(coords, k_value)
    degree <- as.numeric(rowSums(adjacency))
    graph_distance <- distance_to_tumour(adjacency, z$submitted_label == "Tumor")
    ring <- distance_ring(graph_distance)
    one_ring_proxy <- as.numeric(adjacency %*% z$rctd_proxy_full) / degree
    graph_rows[[length(graph_rows) + 1L]] <- data.table(
      sample = sample_value, k = k_value, spots = nrow(z),
      symmetric_edges = sum(adjacency) / 2, min_degree = min(degree),
      median_degree = median(degree), max_degree = max(degree),
      graph_components = components(graph_from_adjacency_matrix(adjacency, mode = "undirected"))$no,
      tumour_spots = sum(ring == "tumor"), boundary_spots = sum(ring == "boundary"),
      unreachable_spots = sum(ring == "unreachable"),
      max_finite_distance = max(graph_distance[is.finite(graph_distance)])
    )
    score_table <- data.table(
      full_axis1 = z$full_axis1,
      malignant_excluded_axis1 = z$malignant_excluded_axis1,
      one_ring_rctd_proxy = one_ring_proxy
    )
    for (block_width in block_width_grid) {
      spatial_block <- paste(floor(z$array_col / block_width),
                             floor(z$array_row / block_width), sep = "_")
      block_ids <- sort(unique(spatial_block))
      boundary_blocks <- uniqueN(spatial_block[ring == "boundary"])
      distal_blocks <- uniqueN(spatial_block[ring == "distal_7plus"])
      evaluable <- boundary_blocks >= 3L && distal_blocks >= 5L
      layer_results <- list()
      for (layer_index in seq_along(score_table)) {
        layer_value <- names(score_table)[[layer_index]]
        block_ring <- data.table(spatial_block = spatial_block, ring = ring,
                                 value = score_table[[layer_index]])[
          ring %chin% c("boundary", "distal_7plus"),
          .(value = median(value)), by = .(spatial_block, ring)
        ]
        result <- bootstrap_effect(
          block_ring, block_ids, n_boot,
          seed + section_index * 100000L + k_value * 1000L +
            as.integer(block_width) * 10L + layer_index
        )
        result[, `:=`(
          sample = sample_value, k = k_value, block_width = block_width,
          layer = layer_value, decision_evaluable = evaluable,
          layer_support = evaluable && estimate > 0 && positive_fraction >= 0.75
        )]
        layer_results[[layer_index]] <- result
      }
      layer_results <- rbindlist(layer_results)
      effect_rows[[length(effect_rows) + 1L]] <- layer_results
      row_index <- row_index + 1L
      section_rows[[row_index]] <- data.table(
        sample = sample_value, k = k_value, block_width = block_width,
        occupied_spatial_blocks = length(block_ids), boundary_blocks = boundary_blocks,
        finite_distal_blocks = distal_blocks, decision_evaluable = evaluable,
        all_three_layers_support = all(layer_results$layer_support),
        full_axis1_effect = layer_results[layer == "full_axis1", estimate],
        full_axis1_positive_fraction = layer_results[layer == "full_axis1", positive_fraction],
        malignant_excluded_effect = layer_results[layer == "malignant_excluded_axis1", estimate],
        malignant_excluded_positive_fraction = layer_results[layer == "malignant_excluded_axis1", positive_fraction],
        one_ring_rctd_proxy_effect = layer_results[layer == "one_ring_rctd_proxy", estimate],
        one_ring_rctd_proxy_positive_fraction = layer_results[layer == "one_ring_rctd_proxy", positive_fraction]
      )
    }
  }
}

effects <- rbindlist(effect_rows)
section_support <- rbindlist(section_rows)
graphs <- unique(rbindlist(graph_rows), by = c("sample", "k"))
grid_summary <- section_support[, .(
  sections = .N,
  decision_evaluable_sections = sum(decision_evaluable),
  supporting_sections = sum(all_three_layers_support),
  combination_supports = sum(all_three_layers_support) >= 3L,
  all_section_effects_positive = all(full_axis1_effect > 0 & malignant_excluded_effect > 0 &
                                       one_ring_rctd_proxy_effect > 0)
), by = .(k, block_width)][order(k, block_width)]

if (nrow(grid_summary) != 9L || nrow(section_support) != 36L || nrow(effects) != 108L) {
  stop(sprintf("Incomplete frozen grid: summary=%d section=%d effects=%d",
               nrow(grid_summary), nrow(section_support), nrow(effects)))
}

primary_old <- fread(primary_effect_file)[
  comparison == "boundary_minus_distal" &
    layer %chin% c("full_axis1", "malignant_excluded_axis1", "one_ring_rctd_proxy_full"),
  .(sample, old_layer = layer, old_estimate = estimate)
]
primary_old[, layer := fifelse(old_layer == "one_ring_rctd_proxy_full", "one_ring_rctd_proxy", old_layer)]
primary_check <- merge(
  effects[k == 6L & block_width == 20, .(sample, layer, gate12ab_estimate = estimate)],
  primary_old[, .(sample, layer, gate12v3_estimate = old_estimate)],
  by = c("sample", "layer"), all = TRUE
)
primary_check[, absolute_difference := abs(gate12ab_estimate - gate12v3_estimate)]
primary_check[, matches := is.finite(absolute_difference) & absolute_difference <= 1e-12]
if (nrow(primary_check) != 12L || !all(primary_check$matches)) stop("Primary k=6, width=20 point estimates do not reproduce Gate12V3")

fwrite(effects, file.path(outdir, "spatial_sensitivity_grid_effects.tsv"), sep = "\t")
fwrite(section_support, file.path(outdir, "spatial_sensitivity_grid_section_support.tsv"), sep = "\t")
fwrite(grid_summary, file.path(outdir, "spatial_sensitivity_grid_summary.tsv"), sep = "\t")
fwrite(graphs, file.path(outdir, "spatial_sensitivity_graph_receipt.tsv"), sep = "\t")
fwrite(primary_check, file.path(outdir, "primary_setting_reproduction.tsv"), sep = "\t")

layer_labels <- c(full_axis1 = "Full Axis1", malignant_excluded_axis1 = "Malignant-excluded Axis1",
                  one_ring_rctd_proxy = "One-ring RCTD proxy")
effects[, layer_label := factor(layer_labels[layer], levels = unname(layer_labels))]
effects[, sample := factor(sample, levels = sort(unique(as.character(sample))))]
effect_limit <- unname(quantile(abs(effects$estimate), 0.98, na.rm = TRUE))

theme_grid <- theme_minimal(base_size = 8.2, base_family = "Arial") +
  theme(plot.title = element_text(size = 9.4, face = "bold"),
        plot.subtitle = element_text(size = 7.2, colour = "#444444"),
        plot.tag = element_text(size = 11.5, face = "bold"),
        axis.title = element_text(size = 8.0), axis.text = element_text(size = 7.2),
        strip.text = element_text(size = 7.2, face = "bold"), panel.grid = element_blank(),
        legend.title = element_text(size = 7.2), legend.text = element_text(size = 6.8),
        plot.margin = margin(5, 6, 5, 6))

p_summary <- ggplot(grid_summary, aes(factor(block_width), factor(k), fill = supporting_sections)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = paste0(supporting_sections, "/4")), size = 3.0, fontface = "bold") +
  scale_fill_gradientn(colours = c("#F2F2F2", "#D9EAF4", "#0077BB"), limits = c(0, 4), breaks = 0:4) +
  labs(title = "Grid-level support", subtitle = "Supporting sections under the frozen three-layer rule",
       x = "Spatial-block width (array units)", y = "Symmetric kNN k", fill = "Sections") + theme_grid

p_effect <- ggplot(effects, aes(factor(block_width), factor(k), fill = estimate)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_point(data = effects[decision_evaluable == TRUE & layer_support == TRUE],
             shape = 21, fill = "black", colour = "black", size = 1.7, stroke = 0.4) +
  geom_point(data = effects[decision_evaluable == TRUE & layer_support == FALSE],
             shape = 21, fill = "white", colour = "black", size = 1.7, stroke = 0.4) +
  geom_text(data = effects[decision_evaluable == FALSE], label = "x", size = 2.7, fontface = "bold") +
  facet_grid(layer_label ~ sample) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-effect_limit, effect_limit), oob = scales::squish) +
  labs(title = "Section- and layer-specific boundary-minus-finite-distal effects",
       subtitle = "Filled dot: evaluable and >=75% positive bootstrap draws; open dot: evaluable without support; x: not evaluable",
       x = "Spatial-block width (array units)", y = "Symmetric kNN k", fill = "Effect") + theme_grid

figure_s4 <- p_summary / p_effect + patchwork::plot_layout(heights = c(0.72, 2.28)) +
  patchwork::plot_annotation(tag_levels = "A")
png_path <- file.path(figure_dir, "FigureS4_spatial_parameter_sensitivity.png")
pdf_path <- file.path(figure_dir, "FigureS4_spatial_parameter_sensitivity.pdf")
ggsave(png_path, figure_s4, width = 7.2, height = 9.0, units = "in", dpi = 400,
       device = grDevices::png, bg = "white", limitsize = FALSE)
ggsave(pdf_path, figure_s4, width = 7.2, height = 9.0, units = "in", device = cairo_pdf,
       bg = "white", limitsize = FALSE)
fwrite(effects, file.path(source_dir, "FigureS4_spatial_grid_all_effects.tsv"), sep = "\t")
fwrite(grid_summary, file.path(source_dir, "FigureS4_spatial_grid_summary.tsv"), sep = "\t")

decision <- data.table(
  gate = "Gate12AB_spatial_sensitivity",
  execution_status = "COMPLETE",
  parameter_combinations = nrow(grid_summary),
  combinations_supporting_frozen_rule = sum(grid_summary$combination_supports),
  combinations_with_all_section_effects_positive = sum(grid_summary$all_section_effects_positive),
  primary_point_estimates_reproduced = all(primary_check$matches),
  independent_animal_inference = FALSE,
  interpretation_scope = "linked_cross_species_parameter_sensitivity"
)
fwrite(decision, file.path(outdir, "gate12ab_spatial_sensitivity_decision.tsv"), sep = "\t")

checkpoint <- c(
  "# Gate12AB spatial parameter-sensitivity checkpoint", "",
  "- Status: COMPLETE", "- Frozen grid: k = {4, 6, 8}; block width = {15, 20, 25}",
  paste0("- Valid block-bootstrap effects per section/layer/combination: ", n_boot),
  paste0("- Supporting combinations: ", sum(grid_summary$combination_supports), "/9"),
  paste0("- Combinations with positive point estimates in all 12 section-layer cells: ",
         sum(grid_summary$all_section_effects_positive), "/9"),
  paste0("- Primary k=6, width=20 point estimates reproduced: ", all(primary_check$matches)),
  "- Reporting scope: linked cross-species descriptive sensitivity; no independent-animal inference"
)
writeLines(checkpoint, file.path(outdir, "GATE12AB_SPATIAL_GRID_CHECKPOINT.md"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))

receipt_files <- list.files(outdir, full.names = TRUE)
receipt_files <- receipt_files[basename(receipt_files) != "SHA256SUMS.txt"]
figure_files <- c(png_path, pdf_path,
                  file.path(source_dir, "FigureS4_spatial_grid_all_effects.tsv"),
                  file.path(source_dir, "FigureS4_spatial_grid_summary.tsv"))
receipt_files <- c(receipt_files, figure_files)
receipt <- data.table(
  path = sub(paste0("^", root, "/"), "", normalizePath(receipt_files, mustWork = TRUE)),
  bytes = file.info(receipt_files)$size,
  sha256 = vapply(receipt_files, digest::digest, character(1), file = TRUE, algo = "sha256")
)
fwrite(receipt, file.path(outdir, "SHA256SUMS.txt"), sep = "\t", col.names = FALSE)
message("Gate12AB spatial sensitivity complete: ", sum(grid_summary$combination_supports), "/9 combinations support the frozen rule")
