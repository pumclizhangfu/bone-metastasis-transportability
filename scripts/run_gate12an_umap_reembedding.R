#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(Seurat)
  library(harmony)
  library(data.table)
  library(uwot)
  library(RANN)
})

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
project <- if (length(args) >= 1L) normalizePath(args[[1L]], mustWork = TRUE) else
  "."
out <- if (length(args) >= 2L) args[[2L]] else
  file.path(project, "results", "gate12an_umap_reembedding")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

seed <- 20260811L
set.seed(seed)
checkpoint_file <- file.path(out, "harmony30_checkpoint.rds")

if (file.exists(checkpoint_file)) {
  message("Loading frozen Harmony checkpoint: ", checkpoint_file)
  checkpoint <- readRDS(checkpoint_file)
  harmony30 <- checkpoint$harmony30
  meta <- checkpoint$meta
} else {
  sce_dir <- file.path(project, "data", "gate3b_work", "annotated_sce")
  files <- sort(list.files(sce_dir, pattern = "\\.rds$", full.names = TRUE))
  if (length(files) != 42L) stop("Expected 42 annotated SCE files; found ", length(files))

  message("Reading 42 annotated SCE objects")
  sces <- lapply(files, readRDS)
  gene_order <- rownames(sces[[1L]])
  if (!all(vapply(sces, function(x) identical(rownames(x), gene_order), logical(1)))) {
    stop("Gene order differs among SCE objects")
  }
  mats <- lapply(sces, function(x) as(assay(x, "counts"), "dgCMatrix"))
  counts_mat <- do.call(cbind, mats)
  meta <- rbindlist(lapply(sces, function(x) as.data.table(as.data.frame(colData(x)))), fill = TRUE)
  if (ncol(counts_mat) != nrow(meta)) stop("Count/metadata cell mismatch")
  if (anyDuplicated(colnames(counts_mat))) stop("Cell identifiers are not globally unique")
  rownames(meta) <- colnames(counts_mat)
  rm(mats, sces); invisible(gc())

  meta[, compartment := factor(compartment, levels = c("distal", "involved", "tumor"))]
  meta[, cancer := factor(cancer, levels = c("prostate", "renal"))]
  meta[, broad_class := factor(
    broad_class,
    levels = c("T_NK", "Myeloid", "B", "Progenitor", "Stromal", "Endothelial",
               "Osteoclast", "Osteoblast", "Malignant", "Erythroid", "Unassigned")
  )]

  obj <- CreateSeuratObject(
    counts = counts_mat, meta.data = as.data.frame(meta),
    project = "paired_bone_metastasis", min.cells = 1, min.features = 0
  )
  rm(counts_mat, meta); invisible(gc())

  message("Normalizing and selecting 3,000 variable features")
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,
                       verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000,
                              verbose = FALSE)
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 40,
                seed.use = seed, verbose = FALSE)

  message("Running Harmony by accession; downstream candidates share this exact representation")
  obj <- RunHarmony(
    obj, group.by.vars = "accession", reduction.use = "pca", dims.use = 1:35,
    assay.use = "RNA", max_iter = 30, plot_convergence = FALSE, verbose = FALSE
  )
  harmony30 <- Embeddings(obj, reduction = "harmony")[, 1:30, drop = FALSE]
  keep_meta <- c("accession", "cancer", "sample_id", "patient_id", "compartment",
                 "broad_class", "broad_class_label", "harmonized_state", "label_source")
  keep_meta <- intersect(keep_meta, colnames(obj@meta.data))
  meta <- as.data.table(obj@meta.data[rownames(harmony30), keep_meta, drop = FALSE],
                        keep.rownames = "cell_id")
  checkpoint <- list(
    harmony30 = harmony30,
    meta = meta,
    variable_features = VariableFeatures(obj),
    pca_stdev = Stdev(obj, reduction = "pca"),
    harmony_stdev = Stdev(obj, reduction = "harmony"),
    seed = seed,
    normalization = "LogNormalize(scale.factor=10000); VST 3000; ScaleData",
    integration = "Harmony(accession; PCA1:35); UMAP input Harmony1:30"
  )
  saveRDS(checkpoint, checkpoint_file, compress = "xz")
  rm(obj, checkpoint); invisible(gc())
}

if (nrow(harmony30) != nrow(meta)) stop("Harmony/metadata row mismatch")
if (!identical(rownames(harmony30), meta$cell_id)) stop("Harmony/metadata cell order mismatch")
if (!"broad_class_label" %in% names(meta)) {
  broad_label_map <- c(
    T_NK = "T/NK", Myeloid = "Myeloid", B = "B cell", Progenitor = "Progenitor",
    Stromal = "Stromal", Endothelial = "Endothelial", Osteoclast = "Osteoclast",
    Osteoblast = "Osteoblast", Malignant = "Malignant", Erythroid = "Erythroid",
    Unassigned = "Unassigned"
  )
  meta[, broad_class_label := unname(broad_label_map[as.character(broad_class)])]
}

candidates <- data.table(
  candidate = c("A_local_crisp", "B_local_balanced", "C_balanced",
                "D_global_crisp", "E_global_balanced", "F_global_open",
                "G_default_open", "H_local_open", "I_compact_spectral",
                "J_compact_pca", "K_compact_spca", "L_global_compact",
                "M_mid_compact", "N_compact_agspectral", "O_compact_lvrandom",
                "P_superglobal_spca"),
  n_neighbors = c(15L, 30L, 45L, 60L, 75L, 100L, 30L, 20L,
                  45L, 45L, 45L, 60L, 55L, 55L, 55L, 100L),
  min_dist = c(0.05, 0.10, 0.20, 0.10, 0.25, 0.40, 0.30, 0.35,
               0.30, 0.30, 0.30, 0.35, 0.32, 0.30, 0.30, 0.35),
  spread = 1.0,
  metric = "cosine",
  n_epochs = 400L,
  init_method = c(rep("spectral", 9L), "pca", "spca", "spca", "spca",
                  "agspectral", "lvrandom", "spca"),
  repulsion_strength = c(rep(1.0, 8L), 0.65, 0.65, 0.65, 0.55, 0.60, 0.65, 0.65, 0.50),
  negative_sample_rate = c(rep(5L, 11L), 4L, 5L, 5L, 5L, 5L),
  seed = seed
)
fwrite(candidates, file.path(out, "candidate_parameters.tsv"), sep = "\t")

## One deterministic, class-stratified sample is used only for neighborhood metrics.
meta[, row_index := seq_len(.N)]
metric_idx <- meta[, {
  take <- min(.N, 1800L)
  .(row_index = sample(row_index, take))
}, by = broad_class]$row_index
metric_idx <- sort(metric_idx)
metric_harmony <- harmony30[metric_idx, , drop = FALSE]
metric_meta <- meta[metric_idx]
metric_k <- 30L
message("Computing cosine-equivalent reference high-dimensional neighbors for ", length(metric_idx), " cells")
## UMAP uses cosine distance. For non-zero vectors, Euclidean neighbours after
## row-wise L2 normalization have the same ordering as cosine-distance neighbours.
metric_norm <- sqrt(rowSums(metric_harmony^2))
if (any(!is.finite(metric_norm) | metric_norm == 0)) stop("Invalid Harmony row norm in metric sample")
metric_harmony_unit <- metric_harmony / metric_norm
hd_nn <- RANN::nn2(metric_harmony_unit, metric_harmony_unit,
                    k = metric_k + 1L)$nn.idx[, -1L, drop = FALSE]
label_vec <- as.character(metric_meta$broad_class)
cancer_vec <- as.character(metric_meta$cancer)
compartment_vec <- as.character(metric_meta$compartment)
row_purity <- function(nn, labels) {
  mean(vapply(seq_len(nrow(nn)), function(i) mean(labels[nn[i, ]] == labels[i]), numeric(1)))
}
hd_class_purity <- row_purity(hd_nn, label_vec)
hd_cancer_same <- row_purity(hd_nn, cancer_vec)
hd_compartment_same <- row_purity(hd_nn, compartment_vec)

metrics <- vector("list", nrow(candidates))
for (i in seq_len(nrow(candidates))) {
  spec <- candidates[i]
  coord_file <- file.path(out, paste0(spec$candidate, "_coordinates.tsv.gz"))
  if (file.exists(coord_file)) {
    message("Reusing completed coordinates: ", spec$candidate)
    coords <- fread(coord_file, select = c("cell_id", "umap_1", "umap_2"))
    if (!identical(coords$cell_id, meta$cell_id)) stop("Cell order mismatch in ", coord_file)
    embedding <- as.matrix(coords[, .(umap_1, umap_2)])
  } else {
    message(sprintf(
      "Running %s: n_neighbors=%d min_dist=%.2f init=%s repulsion=%.2f",
      spec$candidate, spec$n_neighbors, spec$min_dist,
      spec$init_method, spec$repulsion_strength
    ))
    set.seed(spec$seed)
    embedding <- uwot::umap(
      harmony30,
      n_neighbors = spec$n_neighbors,
      n_components = 2,
      metric = spec$metric,
      n_epochs = spec$n_epochs,
      min_dist = spec$min_dist,
      spread = spec$spread,
      init = spec$init_method,
      repulsion_strength = spec$repulsion_strength,
      negative_sample_rate = spec$negative_sample_rate,
      scale = FALSE,
      n_threads = 8L,
      n_sgd_threads = 1L,
      fast_sgd = FALSE,
      verbose = TRUE,
      seed = spec$seed
    )
    coord_dt <- cbind(
      meta[, .(cell_id, accession, cancer, sample_id, patient_id, compartment,
               broad_class, broad_class_label, harmonized_state, label_source)],
      data.table(umap_1 = embedding[, 1L], umap_2 = embedding[, 2L])
    )
    fwrite(coord_dt, coord_file, sep = "\t", compress = "gzip")
  }

  metric_embedding <- embedding[metric_idx, , drop = FALSE]
  umap_nn <- RANN::nn2(metric_embedding, metric_embedding, k = metric_k + 1L)$nn.idx[, -1L, drop = FALSE]
  overlap <- mean(vapply(seq_len(nrow(umap_nn)), function(j) {
    length(intersect(hd_nn[j, ], umap_nn[j, ])) / metric_k
  }, numeric(1)))
  metrics[[i]] <- data.table(
    candidate = spec$candidate,
    metric_cells = length(metric_idx),
    k = metric_k,
    neighbor_overlap_at_k = overlap,
    hd_class_purity = hd_class_purity,
    umap_class_purity = row_purity(umap_nn, label_vec),
    hd_cancer_same = hd_cancer_same,
    umap_cancer_same = row_purity(umap_nn, cancer_vec),
    hd_compartment_same = hd_compartment_same,
    umap_compartment_same = row_purity(umap_nn, compartment_vec)
  )
  fwrite(rbindlist(metrics[seq_len(i)]), file.path(out, "candidate_metrics.tsv"), sep = "\t")
  rm(embedding, metric_embedding, umap_nn); invisible(gc())
}

metrics_dt <- rbindlist(metrics)
fwrite(metrics_dt, file.path(out, "candidate_metrics.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
writeLines(c(
  "# Gate12AN UMAP re-embedding record", "",
  "- Input: 42 annotated discovery SCE objects.",
  "- Shared high-dimensional representation: LogNormalize, 3,000 VST features, PCA, Harmony by accession.",
  "- Every candidate uses Harmony dimensions 1-30 and cosine distance.",
  "- High-dimensional neighbourhood retention uses cosine-equivalent kNN after row-wise L2 normalization.",
  "- Candidate sensitivity analysis varies n_neighbors, min_dist, initialization and repulsion strength.",
  "- Seed: 20260811.",
  "- UMAP is descriptive; 2D distances and cluster shapes are not biological endpoints.",
  "- Candidate choice requires both neighborhood-retention review and blinded visual review."
), file.path(out, "GATE12AN_REEMBEDDING_RECORD.md"))

cat("GATE12AN_STATUS=COMPLETE\n")
