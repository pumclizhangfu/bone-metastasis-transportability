#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(spacexr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: run_gate12e_reconstruction_rctd.R RAW_DIR OUT_DIR [smoke|full]")
}
raw_dir <- normalizePath(args[[1]], mustWork = TRUE)
out_dir <- normalizePath(args[[2]], mustWork = FALSE)
mode <- if (length(args) >= 3L) args[[3]] else "full"
if (!mode %in% c("smoke", "full")) stop("Mode must be smoke or full")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260808)

candidate_genes <- c("Ccl4", "Ccr5", "Cxcl16", "Cxcr6", "Spp1", "Cd44")
spatial_samples <- data.frame(
  sample = c("GSM9564255", "GSM9564256", "GSM9564257", "GSM9564258"),
  library = c("24664", "24665", "24666", "24667"),
  stringsAsFactors = FALSE
)
if (mode == "smoke") spatial_samples <- spatial_samples[1, , drop = FALSE]

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

read_gz_lines <- function(path) {
  connection <- gzfile(path, "rt")
  on.exit(close(connection))
  readLines(connection, warn = FALSE)
}

read_10x_triplet <- function(matrix_path, feature_path, barcode_path) {
  connection <- gzfile(matrix_path, "rb")
  on.exit(close(connection))
  counts <- as(readMM(connection), "dgCMatrix")
  feature_fields <- strsplit(read_gz_lines(feature_path), "\t", fixed = TRUE)
  feature_column <- if (all(lengths(feature_fields) >= 2L)) 2L else 1L
  features <- vapply(feature_fields, `[[`, character(1), feature_column)
  barcodes <- vapply(strsplit(read_gz_lines(barcode_path), "\t", fixed = TRUE), `[[`, character(1), 1L)
  if (nrow(counts) != length(features) || ncol(counts) != length(barcodes)) {
    stop("10x dimension mismatch: ", basename(matrix_path))
  }
  duplicate_features <- sum(duplicated(features))
  rownames(counts) <- make.unique(features)
  colnames(counts) <- barcodes
  list(
    counts = counts,
    features = features,
    barcodes = barcodes,
    feature_column = feature_column,
    duplicate_features = duplicate_features
  )
}

read_scrna <- function() {
  prefix <- file.path(raw_dir, "GSM9564259_Bone_scRNA_filtered_counts")
  data <- read_10x_triplet(
    paste0(prefix, "_matrix.mtx.gz"),
    paste0(prefix, "_features.tsv.gz"),
    paste0(prefix, "_barcodes.tsv.gz")
  )
  metadata <- read.delim(
    gzfile(file.path(raw_dir, "GSM9564259_Bone_scRNA_cell_metadata.tsv.gz")),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  index <- match(data$barcodes, metadata$barcode)
  if (anyNA(index)) stop("scRNA metadata does not cover every barcode")
  metadata <- metadata[index, , drop = FALSE]
  libsize <- Matrix::colSums(data$counts)
  detected <- Matrix::colSums(data$counts > 0)
  mt_index <- grepl("^mt-", rownames(data$counts), ignore.case = TRUE)
  mt_fraction <- if (any(mt_index)) Matrix::colSums(data$counts[mt_index, , drop = FALSE]) / pmax(libsize, 1) else rep(0, length(libsize))
  keep <- detected >= 200 & detected < 7500 & libsize >= 500 & mt_fraction < 0.20 & !is.na(metadata$anno) & metadata$anno != ""
  qc <- data.frame(
    barcode = data$barcodes,
    library_size = as.numeric(libsize),
    detected_genes = as.numeric(detected),
    mitochondrial_fraction = as.numeric(mt_fraction),
    submitted_anno = metadata$anno,
    retained = keep,
    stringsAsFactors = FALSE
  )
  write_tsv(qc, file.path(out_dir, "scrna_cell_qc.tsv.gz"))
  counts <- data$counts[, keep, drop = FALSE]
  metadata <- metadata[keep, , drop = FALSE]
  libsize <- libsize[keep]
  labels <- factor(metadata$anno)
  label_counts <- as.data.frame(table(labels), stringsAsFactors = FALSE)
  names(label_counts) <- c("submitted_anno", "cells")
  label_counts$retained_for_reference <- label_counts$cells >= 25
  write_tsv(label_counts, file.path(out_dir, "scrna_reference_cell_types.tsv"))
  eligible_labels <- label_counts$submitted_anno[label_counts$retained_for_reference]
  eligible <- as.character(labels) %in% eligible_labels
  counts <- counts[, eligible, drop = FALSE]
  metadata <- metadata[eligible, , drop = FALSE]
  libsize <- libsize[eligible]
  labels <- droplevels(labels[eligible])
  names(labels) <- colnames(counts)
  names(libsize) <- colnames(counts)

  detection_rows <- lapply(candidate_genes, function(gene) {
    if (!gene %in% rownames(counts)) {
      return(data.frame(gene = gene, cell_type = c("Macrophage", "T cell"), cells = NA, detected_cells = NA, detection_fraction = NA, status = "GENE_ABSENT"))
    }
    do.call(rbind, lapply(c("Macrophage", "T cell"), function(cell_type) {
      selected <- which(as.character(labels) == cell_type)
      detected_cells <- if (length(selected)) sum(counts[gene, selected, drop = TRUE] > 0) else 0
      data.frame(
        gene = gene,
        cell_type = cell_type,
        cells = length(selected),
        detected_cells = detected_cells,
        detection_fraction = if (length(selected)) detected_cells / length(selected) else NA_real_,
        status = if (length(selected)) "OK" else "CELL_TYPE_ABSENT",
        stringsAsFactors = FALSE
      )
    }))
  })
  detection <- do.call(rbind, detection_rows)
  write_tsv(detection, file.path(out_dir, "scrna_candidate_gene_detection.tsv"))
  list(
    counts = counts,
    labels = labels,
    libsize = libsize,
    metadata = metadata,
    detection = detection,
    total_cells = ncol(data$counts),
    qc_cells = sum(keep),
    reference_cells = ncol(counts),
    duplicate_features = data$duplicate_features
  )
}

read_spatial <- function(sample, library) {
  prefix <- file.path(raw_dir, paste0(sample, "_Bone_ST_filtered_", library))
  data <- read_10x_triplet(
    paste0(prefix, "_matrix.mtx.gz"),
    paste0(prefix, "_features.tsv.gz"),
    paste0(prefix, "_barcodes.tsv.gz")
  )
  metadata <- read.delim(
    gzfile(file.path(raw_dir, paste0(sample, "_Bone_ST_", library, "_spot_metadata.tsv.gz"))),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  metadata$raw_barcode <- sub("_[0-9]+$", "", metadata$barcode)
  metadata_index <- match(data$barcodes, metadata$raw_barcode)
  if (anyNA(metadata_index)) stop(sample, ": metadata does not cover every filtered spot")
  metadata <- metadata[metadata_index, , drop = FALSE]
  rownames(metadata) <- data$barcodes
  positions <- read.csv(
    gzfile(file.path(raw_dir, paste0(sample, "_Bone_ST_", library, "_tissue_positions_list.csv.gz"))),
    header = FALSE, stringsAsFactors = FALSE
  )
  if (ncol(positions) != 6L) stop(sample, ": unexpected tissue-position width")
  names(positions) <- c("barcode", "in_tissue", "array_row", "array_col", "pixel_row", "pixel_col")
  position_index <- match(data$barcodes, positions$barcode)
  if (anyNA(position_index)) stop(sample, ": positions do not cover every filtered spot")
  positions <- positions[position_index, , drop = FALSE]
  libsize <- Matrix::colSums(data$counts)
  detected <- Matrix::colSums(data$counts > 0)
  keep <- detected >= 200 & libsize >= 500 & positions$in_tissue == 1
  qc <- data.frame(
    sample = sample,
    barcode = data$barcodes,
    library_size = as.numeric(libsize),
    detected_genes = as.numeric(detected),
    in_tissue = positions$in_tissue,
    submitted_label = metadata$S2C2_label,
    retained = keep,
    stringsAsFactors = FALSE
  )
  counts <- data$counts[, keep, drop = FALSE]
  metadata <- metadata[keep, , drop = FALSE]
  positions <- positions[keep, , drop = FALSE]
  libsize <- libsize[keep]
  rownames(positions) <- colnames(counts)
  list(
    counts = counts,
    metadata = metadata,
    positions = positions,
    libsize = libsize,
    qc = qc,
    duplicate_features = data$duplicate_features,
    total_spots = length(keep),
    retained_spots = sum(keep)
  )
}

message("Reading and QC-filtering the scRNA reference")
scrna <- read_scrna()
reference <- Reference(
  scrna$counts,
  scrna$labels,
  nUMI = scrna$libsize,
  n_max_cells = 20000,
  min_UMI = 100
)

section_objects <- list()
section_qc <- list()
rctd_receipts <- list()
for (i in seq_len(nrow(spatial_samples))) {
  sample <- spatial_samples$sample[[i]]
  library <- spatial_samples$library[[i]]
  message("Reconstructing ", sample, " / ", library)
  spatial <- read_spatial(sample, library)
  common_genes <- intersect(rownames(spatial$counts), rownames(scrna$counts))
  if (length(common_genes) < 5000L) {
    stop(sample, ": only ", length(common_genes), " common genes between spatial and scRNA matrices")
  }
  message("Common spatial/reference genes for ", sample, ": ", length(common_genes))
  write_tsv(spatial$qc, file.path(out_dir, paste0(sample, "_spot_qc.tsv.gz")))
  coords <- data.frame(
    x = spatial$positions$array_col,
    y = spatial$positions$array_row,
    row.names = rownames(spatial$positions)
  )
  spatial_rna <- SpatialRNA(coords, spatial$counts, nUMI = spatial$libsize)
  message("Running RCTD for ", sample, " (mode=", mode, ")")
  rctd <- create.RCTD(
    spatial_rna,
    reference,
    max_cores = if (mode == "smoke") 2 else 8,
    test_mode = mode == "smoke",
    CELL_MIN_INSTANCE = 25,
    keep_reference = FALSE
  )
  rctd <- run.RCTD(rctd, doublet_mode = "full")
  weights <- normalize_weights(rctd@results$weights)
  weights <- as.matrix(weights)
  if (mode == "full") {
    weight_index <- match(colnames(spatial$counts), rownames(weights))
    if (anyNA(weight_index)) stop(sample, ": full-mode RCTD weights do not cover every retained spot")
    weights <- weights[weight_index, , drop = FALSE]
  } else {
    if (!all(rownames(weights) %in% colnames(spatial$counts))) {
      stop(sample, ": smoke-mode RCTD returned unknown spot barcodes")
    }
    smoke_barcodes <- rownames(weights)
    spatial$counts <- spatial$counts[, smoke_barcodes, drop = FALSE]
    spatial$metadata <- spatial$metadata[smoke_barcodes, , drop = FALSE]
    spatial$positions <- spatial$positions[smoke_barcodes, , drop = FALSE]
    spatial$libsize <- spatial$libsize[smoke_barcodes]
    message(
      "Smoke-mode RCTD returned ", nrow(weights), " of ",
      spatial$retained_spots, " retained spots"
    )
  }
  write_tsv(
    data.frame(barcode = rownames(weights), weights, check.names = FALSE),
    file.path(out_dir, paste0(sample, "_RCTD_weights.tsv.gz"))
  )
  genes_present <- candidate_genes[candidate_genes %in% rownames(spatial$counts)]
  expression <- matrix(0, nrow = length(candidate_genes), ncol = ncol(spatial$counts), dimnames = list(candidate_genes, colnames(spatial$counts)))
  if (length(genes_present)) {
    expression[genes_present, ] <- log1p(
      1e4 * as.matrix(spatial$counts[genes_present, , drop = FALSE]) /
        matrix(as.numeric(spatial$libsize), nrow = length(genes_present), ncol = ncol(spatial$counts), byrow = TRUE)
    )
  }
  section_objects[[sample]] <- list(
    sample = sample,
    library = library,
    coords = coords,
    positions = spatial$positions,
    submitted_label = setNames(spatial$metadata$S2C2_label, colnames(spatial$counts)),
    weights = weights,
    candidate_expression = expression,
    candidate_detection = spatial$counts[candidate_genes[candidate_genes %in% rownames(spatial$counts)], , drop = FALSE] > 0,
    retained_spots = spatial$retained_spots,
    total_spots = spatial$total_spots
  )
  section_qc[[sample]] <- data.frame(
    sample = sample,
    library = library,
    total_filtered_spots = spatial$total_spots,
    retained_spots = spatial$retained_spots,
    duplicate_features = spatial$duplicate_features,
    submitted_tumor_spots = sum(spatial$metadata$S2C2_label == "Tumor"),
    submitted_msc_spots = sum(spatial$metadata$S2C2_label == "MSCs"),
    stringsAsFactors = FALSE
  )
  rctd_receipts[[sample]] <- data.frame(
    sample = sample,
    mode = mode,
    spots = nrow(weights),
    retained_input_spots = spatial$retained_spots,
    rctd_coverage_fraction = nrow(weights) / spatial$retained_spots,
    cell_types = ncol(weights),
    macrophage_column = "Macrophage" %in% colnames(weights),
    t_cell_column = "T cell" %in% colnames(weights),
    finite_weights = all(is.finite(weights)),
    row_sum_min = min(rowSums(weights)),
    row_sum_max = max(rowSums(weights)),
    status = if (all(is.finite(weights)) && "Macrophage" %in% colnames(weights) && "T cell" %in% colnames(weights)) "PASS" else "FAIL",
    stringsAsFactors = FALSE
  )
  rm(rctd, spatial_rna, spatial)
  gc()
}

section_qc_table <- do.call(rbind, section_qc)
rctd_receipt_table <- do.call(rbind, rctd_receipts)
write_tsv(section_qc_table, file.path(out_dir, "spatial_section_qc.tsv"))
write_tsv(rctd_receipt_table, file.path(out_dir, "rctd_receipt.tsv"))
saveRDS(section_objects, file.path(out_dir, "gate12e_section_objects.rds"), compress = "xz")
saveRDS(
  list(
    mode = mode,
    seed = 20260808,
    scrna_total_cells = scrna$total_cells,
    scrna_qc_cells = scrna$qc_cells,
    scrna_reference_cells = scrna$reference_cells,
    scrna_duplicate_features = scrna$duplicate_features,
    candidate_detection = scrna$detection,
    section_qc = section_qc_table,
    rctd_receipt = rctd_receipt_table
  ),
  file.path(out_dir, "gate12e_reconstruction_receipt.rds")
)
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo_reconstruction.txt"))

overall <- all(rctd_receipt_table$status == "PASS") && nrow(rctd_receipt_table) == nrow(spatial_samples)
cat("GATE12E_RECONSTRUCTION=", if (overall) "PASS" else "FAIL", "\n", sep = "")
cat("MODE=", mode, "\n", sep = "")
cat("SCRNA_REFERENCE_CELLS=", scrna$reference_cells, "\n", sep = "")
cat("SPATIAL_SECTIONS=", nrow(section_qc_table), "\n", sep = "")
