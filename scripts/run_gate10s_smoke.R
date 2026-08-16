#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(SingleCellExperiment)
  library(scuttle)
  library(SingleR)
  library(UCell)
  library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: run_gate10s_smoke.R PROJECT GSE225_EXTRACTED BREAST_SELECTED_DIR OUT_DIR", call. = FALSE)
}
project <- normalizePath(args[[1]], mustWork = TRUE)
lung_dir <- normalizePath(args[[2]], mustWork = TRUE)
breast_selected_dir <- normalizePath(args[[3]], mustWork = TRUE)
out_dir <- normalizePath(args[[4]], mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260807)

read_tsv <- function(path) fread(path, sep = "\t", header = TRUE)
split_markers <- function(x) unique(strsplit(x, "\\|")[[1L]])
detected_count <- function(mat, markers) {
  idx <- match(markers, rownames(mat), nomatch = 0L)
  idx <- idx[idx > 0L]
  if (!length(idx)) return(integer(ncol(mat)))
  as.integer(Matrix::colSums(mat[idx, , drop = FALSE] > 0))
}
aggregate_symbols <- function(mat, symbols) {
  keep <- !is.na(symbols) & nzchar(symbols)
  mat <- mat[keep, , drop = FALSE]; symbols <- symbols[keep]
  u <- unique(symbols)
  if (length(u) == length(symbols)) {
    rownames(mat) <- u
    return(as(mat, "dgCMatrix"))
  }
  map <- sparseMatrix(i = match(symbols, u), j = seq_along(symbols), x = 1,
                      dims = c(length(u), length(symbols)))
  ans <- map %*% mat
  rownames(ans) <- u; colnames(ans) <- colnames(mat)
  as(ans, "dgCMatrix")
}

signature_table <- read_tsv(file.path(project, "config/gate9_frozen_ecotype_definition.tsv"))
gate_table <- read_tsv(file.path(project, "config/gate9a_raw_state_gates.tsv"))
monaco <- readRDS(file.path(project, "metadata/monaco_immune_data.rds"))
thresholds <- read_tsv(file.path(project, "results/gate9a_reference/raw_reconstruction_v1/cell_margin_thresholds.tsv"))
primary_states <- signature_table[primary_role == "primary", state_id]
all_states <- signature_table$state_id
signatures <- setNames(lapply(signature_table$analysis_markers, split_markers), signature_table$state_id)
gate_markers <- setNames(lapply(gate_table$core_markers, split_markers), gate_table$state_id)
gate_min <- setNames(as.integer(gate_table$min_core_detected), gate_table$state_id)
myeloid_fallback <- c("LST1", "TYROBP", "FCER1G", "CTSS", "AIF1")
t_fallback <- c("CD3D", "CD3E", "TRAC", "CD247")
cd4_support <- c("CD4", "IL7R", "LTB")
cd8_support <- c("CD8A", "CD8B")
nk_exclusion <- c("KLRD1", "NCR1", "NCAM1", "FCGR3A")
myeloid_labels <- c("Classical monocytes", "Intermediate monocytes", "Non classical monocytes", "Myeloid dendritic cells")
t_labels <- c("Naive CD4 T cells", "Naive CD8 T cells", "Central memory CD8 T cells", "Effector memory CD8 T cells",
              "T regulatory cells", "Follicular helper T cells", "Terminal effector CD4 T cells",
              "Terminal effector CD8 T cells", "Th1 cells", "Th1/Th17 cells", "Th17 cells", "Th2 cells")
myeloid_margin <- thresholds[lineage == "myeloid", threshold][[1L]]
t_margin <- thresholds[lineage == "t_cell", threshold][[1L]]

make_core_eligibility <- function(mat) {
  core_counts <- lapply(gate_markers, function(z) detected_count(mat, z))
  cd4_n <- detected_count(mat, cd4_support)
  cd8_n <- detected_count(mat, cd8_support)
  raw <- sapply(names(core_counts), function(s) {
    ok <- core_counts[[s]] >= gate_min[[s]]
    if (s == "CD4_TREG") ok <- ok & cd4_n >= 1L
    if (s %chin% c("CD8_TEX", "CD8_PTEX")) ok <- ok & cd8_n >= 1L
    ok
  })
  if (is.null(dim(raw))) raw <- matrix(raw, ncol = 1L)
  colnames(raw) <- names(core_counts)
  list(raw = raw, cd4_n = cd4_n, cd8_n = cd8_n)
}

lung_prefix <- "GSM7041480_sg1"
lung_files <- c(
  matrix = file.path(lung_dir, paste0(lung_prefix, "-matrix.mtx.gz")),
  features = file.path(lung_dir, paste0(lung_prefix, "-features.tsv.gz")),
  barcodes = file.path(lung_dir, paste0(lung_prefix, "-barcodes.tsv.gz"))
)
if (!all(file.exists(lung_files))) stop("Missing GSE225209 smoke triplet", call. = FALSE)
feat <- fread(lung_files[["features"]], header = FALSE)
if (ncol(feat) == 2L) feat[, V3 := "Gene Expression"]
bc <- fread(lung_files[["barcodes"]], header = FALSE)[[1L]]
raw <- as(Matrix::readMM(gzfile(lung_files[["matrix"]])), "dgCMatrix")
gex_idx <- which(feat[[3L]] == "Gene Expression")
gex <- aggregate_symbols(raw[gex_idx, , drop = FALSE], feat[[2L]][gex_idx])
colnames(gex) <- bc
n_count <- as.numeric(Matrix::colSums(gex))
n_feature <- as.integer(Matrix::colSums(gex > 0))
mito <- startsWith(rownames(gex), "MT-")
pct_mito <- if (any(mito)) 100 * as.numeric(Matrix::colSums(gex[mito, , drop = FALSE])) / pmax(n_count, 1) else rep(0, ncol(gex))
qc <- n_feature >= 200 & n_feature <= 8000 & n_count >= 500 & pct_mito <= 20
lung_mat <- gex[, qc, drop = FALSE]
colnames(lung_mat) <- paste0("LUCA_BM_01::", colnames(lung_mat))
lung_core <- make_core_eligibility(lung_mat)
candidate <- rowSums(lung_core$raw[, primary_states, drop = FALSE]) > 0
sr_label <- sr_pruned <- rep(NA_character_, ncol(lung_mat))
if (any(candidate)) {
  sce <- SingleCellExperiment(list(counts = lung_mat[, candidate, drop = FALSE]))
  sce <- scuttle::logNormCounts(sce)
  pred <- SingleR(test = sce, ref = monaco, labels = SummarizedExperiment::colData(monaco)$label.fine,
                  fine.tune = TRUE, prune = TRUE, assay.type.test = "logcounts", assay.type.ref = "logcounts",
                  BPPARAM = MulticoreParam(8L))
  sr_label[candidate] <- as.character(pred$labels)
  sr_pruned[candidate] <- as.character(pred$pruned.labels)
}
myeloid_n <- detected_count(lung_mat, myeloid_fallback)
t_n <- detected_count(lung_mat, t_fallback)
nk_n <- detected_count(lung_mat, nk_exclusion)
myeloid_ok <- (!is.na(sr_pruned) & sr_pruned %chin% myeloid_labels) | myeloid_n >= 2L
t_ok <- (!is.na(sr_pruned) & sr_pruned %chin% t_labels) | (t_n >= 2L & nk_n < 2L)
dual <- myeloid_ok & t_ok
myeloid_ok[dual] <- FALSE; t_ok[dual] <- FALSE
lung_scores <- UCell::ScoreSignatures_UCell(matrix = lung_mat, features = signatures, maxRank = 1500,
                                            BPPARAM = MulticoreParam(8L), ncores = 8L, force.gc = TRUE)
lung_scores <- as.matrix(lung_scores)
colnames(lung_scores) <- sub("_UCell$", "", colnames(lung_scores))
eligible <- lung_core$raw
eligible[, c("CD14HI_MONO", "CD16HI_MONO", "MACROPHAGE", "OSTEOCLAST")] <-
  eligible[, c("CD14HI_MONO", "CD16HI_MONO", "MACROPHAGE", "OSTEOCLAST"), drop = FALSE] & myeloid_ok
eligible[, c("CD4_TREG", "CD8_TEX", "CD8_PTEX")] <-
  eligible[, c("CD4_TREG", "CD8_TEX", "CD8_PTEX"), drop = FALSE] & t_ok
lung_assignment <- rep("unassigned", ncol(lung_mat))
for (lineage_states in list(c("CD14HI_MONO", "CD16HI_MONO", "MACROPHAGE", "OSTEOCLAST"), c("CD4_TREG", "CD8_TEX"))) {
  idx <- which(rowSums(eligible[, lineage_states, drop = FALSE]) > 0)
  for (j in idx) {
    elig <- lineage_states[eligible[j, lineage_states]]
    win <- elig[which.max(lung_scores[j, elig])]
    competitors <- setdiff(lineage_states, win)
    margin <- lung_scores[j, win] - max(lung_scores[j, competitors])
    threshold <- if (win %chin% c("CD4_TREG", "CD8_TEX")) t_margin else myeloid_margin
    if (is.finite(margin) && margin >= threshold) lung_assignment[j] <- win
  }
}
lung_counts <- data.table(dataset = "GSE225209", sample_or_lesion = "sz", patient_id = "LUCA_BM_01",
                          final_assignment = lung_assignment)[, .N, by = .(dataset, sample_or_lesion, patient_id, final_assignment)]

breast_matrix_path <- file.path(breast_selected_dir, "GSE190772_BoM_logCounts.selected_markers.tsv.gz")
breast_meta_path <- file.path(project, "data/raw/GSE190772/GSE190772_BoM_MetaData.txt.gz")
bm <- fread(breast_matrix_path)
genes <- bm[[1L]]
bm_mat <- as.matrix(bm[, -1L, with = FALSE])
storage.mode(bm_mat) <- "double"
rownames(bm_mat) <- genes
meta <- fread(breast_meta_path)
setkey(meta, Cell)
lesion <- "BoM1"
lesion_cells <- names(bm)[-1L][startsWith(names(bm)[-1L], paste0(lesion, "_"))]
meta_lesion <- meta[J(lesion_cells)]
if (anyNA(meta_lesion$Cell)) stop("BoM1 metadata alignment failed", call. = FALSE)
qc_b <- meta_lesion$nFeature_RNA >= 200 & meta_lesion$nFeature_RNA <= 8000 &
  meta_lesion$nCount_RNA >= 500 & meta_lesion$percent.mt <= 20
bm_mat <- bm_mat[, lesion_cells[qc_b], drop = FALSE]
meta_lesion <- meta_lesion[qc_b]
bcore <- make_core_eligibility(bm_mat)
myeloid_author <- meta_lesion$cellType %chin% c("Macrophage", "Osteoclast", "Monocyte")
t_author <- meta_lesion$cellType == "T"
bdual <- myeloid_author & t_author
myeloid_author[bdual] <- FALSE; t_author[bdual] <- FALSE
beligible <- bcore$raw[, primary_states, drop = FALSE]
beligible[, c("CD14HI_MONO", "CD16HI_MONO", "MACROPHAGE", "OSTEOCLAST")] <-
  beligible[, c("CD14HI_MONO", "CD16HI_MONO", "MACROPHAGE", "OSTEOCLAST"), drop = FALSE] & myeloid_author
beligible[, c("CD4_TREG", "CD8_TEX")] <-
  beligible[, c("CD4_TREG", "CD8_TEX"), drop = FALSE] & t_author
breast_assignment <- rep("unassigned", ncol(bm_mat))
for (lineage_states in list(c("CD14HI_MONO", "CD16HI_MONO", "MACROPHAGE", "OSTEOCLAST"), c("CD4_TREG", "CD8_TEX"))) {
  idx <- which(rowSums(beligible[, lineage_states, drop = FALSE]) == 1L)
  if (length(idx)) breast_assignment[idx] <- lineage_states[max.col(beligible[idx, lineage_states, drop = FALSE], ties.method = "first")]
}
breast_counts <- data.table(dataset = "GSE190772", sample_or_lesion = lesion, patient_id = "BRCA_BM_01",
                            final_assignment = breast_assignment)[, .N, by = .(dataset, sample_or_lesion, patient_id, final_assignment)]

counts <- rbind(lung_counts, breast_counts)
fwrite(counts, file.path(out_dir, "smoke_state_counts.tsv"), sep = "\t")
qc_summary <- data.table(
  dataset = c("GSE225209", "GSE190772"),
  sample_or_lesion = c("sz", "BoM1"),
  input_type = c("raw_10x_counts", "author_log_normalized"),
  raw_or_author_cells = c(ncol(gex), length(lesion_cells)),
  qc_cells = c(ncol(lung_mat), ncol(bm_mat)),
  assigned_cells = c(sum(lung_assignment != "unassigned"), sum(breast_assignment != "unassigned")),
  assigned_fraction = c(mean(lung_assignment != "unassigned"), mean(breast_assignment != "unassigned"))
)
qc_summary[, S3_ge_500_qc := qc_cells >= 500]
qc_summary[, S4_ge_50_assigned := assigned_cells >= 50]
fwrite(qc_summary, file.path(out_dir, "smoke_qc_summary.tsv"), sep = "\t")
smoke_pass <- all(qc_summary$S3_ge_500_qc) && all(qc_summary$S4_ge_50_assigned)
writeLines(c(
  "# Gate10S two-path smoke audit", "",
  paste0("- Decision: **", if (smoke_pass) "PASS" else "FAIL", "**"),
  "- GSE225209 path: frozen raw-count gates, Monaco, UCell and source margins",
  "- GSE190772 path: author-lineage-constrained, unique-core normalized projection",
  "- Biological confirmation claimed: **NO**", "- Full run authorized: **only if PASS**", ""
), file.path(out_dir, "GATE10S_SMOKE_AUDIT.md"))
files <- list.files(out_dir, full.names = TRUE)
files <- files[basename(files) != "SHA256SUMS"]
status <- system(sprintf("cd %s && sha256sum %s > SHA256SUMS", shQuote(out_dir), paste(shQuote(basename(files)), collapse = " ")))
if (status != 0L) stop("SHA256 manifest failed", call. = FALSE)
cat("GATE10S_SMOKE_DECISION=", if (smoke_pass) "PASS" else "FAIL", "\n", sep = "")
