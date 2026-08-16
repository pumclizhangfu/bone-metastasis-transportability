#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(uwot)
})

options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Usage: run_gate12ay_lineage_umap.R PROJECT PANEL")
project <- normalizePath(args[[1L]], mustWork = TRUE)
panel <- tolower(args[[2L]])
panel_name <- panel

specs <- data.table(
  panel = c("lymphoid", "t_nk", "b_cell", "myeloid", "niche"),
  n_neighbors = c(40L, 40L, 25L, 35L, 30L),
  min_dist = c(0.25, 0.25, 0.25, 0.25, 0.25),
  repulsion_strength = c(0.80, 0.80, 0.80, 0.80, 0.80),
  seed = c(20260821L, 20260824L, 20260825L, 20260822L, 20260823L)
)
spec <- specs[panel == panel_name]
if (nrow(spec) != 1L) stop("Unknown lineage panel: ", panel)

source_dir <- file.path(project, "results", "gate12an_umap_reembedding")
out <- file.path(project, "results", "gate12ay_multiscale_atlas", "source_data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
checkpoint <- readRDS(file.path(source_dir, "harmony30_checkpoint.rds"))
harmony30 <- checkpoint$harmony30
meta <- as.data.table(checkpoint$meta)
if (!identical(rownames(harmony30), meta$cell_id)) stop("Harmony/metadata order mismatch")
if (!"broad_class_label" %in% names(meta)) {
  broad_label_map <- c(
    T_NK = "T/NK", Myeloid = "Myeloid", B = "B cell", Progenitor = "Progenitor",
    Stromal = "Stromal", Endothelial = "Endothelial", Osteoclast = "Osteoclast",
    Osteoblast = "Osteoblast", Malignant = "Malignant", Erythroid = "Erythroid",
    Unassigned = "Unassigned"
  )
  meta[, broad_class_label := unname(broad_label_map[as.character(broad_class)])]
}

meta[, lineage_panel := fcase(
  broad_class %chin% c("T_NK", "B"), "lymphoid",
  broad_class %chin% c("Myeloid", "Osteoclast"), "myeloid",
  default = "niche"
)]
idx <- switch(
  panel,
  t_nk = which(meta$broad_class == "T_NK"),
  b_cell = which(meta$broad_class == "B"),
  which(meta$lineage_panel == panel)
)
if (length(idx) < 100L) stop("Too few cells for panel: ", panel)

set.seed(spec$seed)
message(sprintf(
  "Running %s local UMAP: cells=%d n_neighbors=%d min_dist=%.2f repulsion=%.2f",
  panel, length(idx), spec$n_neighbors, spec$min_dist, spec$repulsion_strength
))
embedding <- uwot::umap(
  harmony30[idx, , drop = FALSE],
  n_neighbors = spec$n_neighbors,
  n_components = 2L,
  metric = "cosine",
  n_epochs = 400L,
  min_dist = spec$min_dist,
  spread = 1.0,
  init = "spca",
  repulsion_strength = spec$repulsion_strength,
  negative_sample_rate = 5L,
  scale = FALSE,
  n_threads = 4L,
  n_sgd_threads = 1L,
  fast_sgd = FALSE,
  verbose = TRUE,
  seed = spec$seed
)

result <- cbind(
  meta[idx, .(cell_id, accession, cancer, sample_id, patient_id, compartment,
              broad_class, broad_class_label, harmonized_state, label_source,
              lineage_panel)],
  data.table(local_umap_1 = embedding[, 1L], local_umap_2 = embedding[, 2L])
)
result[, lineage_panel := panel]
output_file <- file.path(out, sprintf("Gate12AY_%s_local_UMAP.tsv.gz", panel))
temporary_file <- paste0(output_file, ".incomplete")
fwrite(result, temporary_file, sep = "\t", compress = "gzip")
if (!file.rename(temporary_file, output_file)) stop("Could not publish ", output_file)
fwrite(spec, file.path(out, sprintf("Gate12AY_%s_parameters.tsv", panel)), sep = "\t")
cat(sprintf("GATE12AY_LINEAGE_STATUS=COMPLETE PANEL=%s CELLS=%d\n", panel, nrow(result)))
