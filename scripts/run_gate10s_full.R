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

options(stringsAsFactors = FALSE, warn = 1)
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) != 4L) {
  stop("Usage: run_gate10s_full.R PROJECT GSE225_EXTRACTED BREAST_SELECTED_DIR OUT_DIR", call. = FALSE)
}
project <- normalizePath(argv[[1]], mustWork = TRUE)
lung_dir <- normalizePath(argv[[2]], mustWork = TRUE)
breast_selected_dir <- normalizePath(argv[[3]], mustWork = TRUE)
out_dir <- normalizePath(argv[[4]], mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260807)
workers <- 8L

log_path <- file.path(out_dir, "gate10s_full.log")
log_con <- file(log_path, open = "wt")
on.exit(try(close(log_con), silent = TRUE), add = TRUE)
log_msg <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(line, "\n")
  writeLines(line, log_con); flush(log_con)
}
read_tsv <- function(path) fread(path, sep = "\t", header = TRUE)
split_markers <- function(x) unique(strsplit(x, "\\|")[[1L]])
bool <- function(x) toupper(as.character(x)) == "TRUE"

signature_table <- read_tsv(file.path(project, "config/gate9_frozen_ecotype_definition.tsv"))
gate_table <- read_tsv(file.path(project, "config/gate9a_raw_state_gates.tsv"))
source_scaling <- read_tsv(file.path(project, "results/gate9a_reference/raw_reconstruction_v1/feature_scaling.tsv"))
monaco <- readRDS(file.path(project, "metadata/monaco_immune_data.rds"))
thresholds <- read_tsv(file.path(project, "results/gate9a_reference/raw_reconstruction_v1/cell_margin_thresholds.tsv"))
crosswalk <- read_tsv(file.path(project, "results/gate10s_supportive_projection/input_audit_v1/sample_lesion_crosswalk.tsv"))

primary_states <- signature_table[primary_role == "primary", state_id]
expected_states <- c("CD14HI_MONO", "CD16HI_MONO", "MACROPHAGE", "OSTEOCLAST", "CD4_TREG", "CD8_TEX")
if (!identical(primary_states, expected_states)) stop("Frozen state order mismatch", call. = FALSE)
all_states <- signature_table$state_id
signatures <- setNames(lapply(signature_table$analysis_markers, split_markers), signature_table$state_id)
gate_markers <- setNames(lapply(gate_table$core_markers, split_markers), gate_table$state_id)
gate_min <- setNames(as.integer(gate_table$min_core_detected), gate_table$state_id)
source_mean <- setNames(source_scaling$training_mean, source_scaling$state_id)[primary_states]
source_sd <- setNames(source_scaling$training_sd, source_scaling$state_id)[primary_states]

myeloid_fallback <- c("LST1", "TYROBP", "FCER1G", "CTSS", "AIF1")
t_fallback <- c("CD3D", "CD3E", "TRAC", "CD247")
cd4_support <- c("CD4", "IL7R", "LTB")
cd8_support <- c("CD8A", "CD8B")
nk_exclusion <- c("KLRD1", "NCR1", "NCAM1", "FCGR3A")
myeloid_labels <- c("Classical monocytes", "Intermediate monocytes", "Non classical monocytes", "Myeloid dendritic cells")
t_labels <- c(
  "Naive CD4 T cells", "Naive CD8 T cells", "Central memory CD8 T cells",
  "Effector memory CD8 T cells", "T regulatory cells", "Follicular helper T cells",
  "Terminal effector CD4 T cells", "Terminal effector CD8 T cells", "Th1 cells",
  "Th1/Th17 cells", "Th17 cells", "Th2 cells"
)
myeloid_margin <- thresholds[lineage == "myeloid", threshold][[1L]]
t_margin <- thresholds[lineage == "t_cell", threshold][[1L]]

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

make_core_eligibility <- function(mat) {
  core_counts <- lapply(gate_markers, function(z) detected_count(mat, z))
  cd4_n <- detected_count(mat, cd4_support)
  cd8_n <- detected_count(mat, cd8_support)
  raw <- sapply(names(core_counts), function(state) {
    ok <- core_counts[[state]] >= gate_min[[state]]
    if (state == "CD4_TREG") ok <- ok & cd4_n >= 1L
    if (state %chin% c("CD8_TEX", "CD8_PTEX")) ok <- ok & cd8_n >= 1L
    ok
  })
  if (is.null(dim(raw))) raw <- matrix(raw, ncol = 1L)
  colnames(raw) <- names(core_counts)
  list(raw = raw, cd4_n = cd4_n, cd8_n = cd8_n)
}

score_ucell <- function(mat) {
  score <- UCell::ScoreSignatures_UCell(
    matrix = mat, features = signatures, maxRank = 1500,
    BPPARAM = MulticoreParam(workers), ncores = workers, force.gc = TRUE
  )
  score <- as.matrix(score)
  colnames(score) <- sub("_UCell$", "", colnames(score))
  if (!all(all_states %chin% colnames(score))) stop("UCell output columns mismatch", call. = FALSE)
  if (!is.null(rownames(score)) && all(colnames(mat) %chin% rownames(score))) {
    score <- score[colnames(mat), , drop = FALSE]
  }
  score
}

process_lung <- function(gsm, stem, sample_name, patient_id) {
  prefix <- paste(gsm, stem, sep = "_")
  files <- c(
    matrix = file.path(lung_dir, paste0(prefix, "-matrix.mtx.gz")),
    features = file.path(lung_dir, paste0(prefix, "-features.tsv.gz")),
    barcodes = file.path(lung_dir, paste0(prefix, "-barcodes.tsv.gz"))
  )
  if (!all(file.exists(files))) stop("Missing 10X triplet for ", sample_name, call. = FALSE)
  log_msg("START lung sample=", sample_name, " patient=", patient_id)
  feat <- fread(files[["features"]], header = FALSE)
  if (ncol(feat) == 2L) feat[, V3 := "Gene Expression"]
  barcodes <- fread(files[["barcodes"]], header = FALSE)[[1L]]
  raw <- as(Matrix::readMM(gzfile(files[["matrix"]])), "dgCMatrix")
  if (nrow(raw) != nrow(feat) || ncol(raw) != length(barcodes)) stop("10X dimension mismatch", call. = FALSE)
  gex_idx <- which(feat[[3L]] == "Gene Expression")
  gex <- aggregate_symbols(raw[gex_idx, , drop = FALSE], feat[[2L]][gex_idx])
  colnames(gex) <- barcodes
  n_count <- as.numeric(Matrix::colSums(gex))
  n_feature <- as.integer(Matrix::colSums(gex > 0))
  mito <- startsWith(rownames(gex), "MT-")
  pct_mito <- if (any(mito)) 100 * as.numeric(Matrix::colSums(gex[mito, , drop = FALSE])) / pmax(n_count, 1) else rep(0, ncol(gex))
  qc <- n_feature >= 200 & n_feature <= 8000 & n_count >= 500 & pct_mito <= 20
  mat <- gex[, qc, drop = FALSE]
  cell_ids <- paste0(patient_id, "::", colnames(mat))
  colnames(mat) <- cell_ids
  core <- make_core_eligibility(mat)
  candidate <- rowSums(core$raw[, primary_states, drop = FALSE]) > 0
  sr_label <- sr_pruned <- rep(NA_character_, ncol(mat))
  if (any(candidate)) {
    sce <- SingleCellExperiment(list(counts = mat[, candidate, drop = FALSE]))
    sce <- scuttle::logNormCounts(sce)
    pred <- SingleR(
      test = sce, ref = monaco, labels = SummarizedExperiment::colData(monaco)$label.fine,
      fine.tune = TRUE, prune = TRUE, assay.type.test = "logcounts", assay.type.ref = "logcounts",
      BPPARAM = MulticoreParam(workers)
    )
    sr_label[candidate] <- as.character(pred$labels)
    sr_pruned[candidate] <- as.character(pred$pruned.labels)
    rm(sce, pred); gc(FALSE)
  }
  myeloid_n <- detected_count(mat, myeloid_fallback)
  t_n <- detected_count(mat, t_fallback)
  nk_n <- detected_count(mat, nk_exclusion)
  myeloid_ok <- (!is.na(sr_pruned) & sr_pruned %chin% myeloid_labels) | myeloid_n >= 2L
  t_ok <- (!is.na(sr_pruned) & sr_pruned %chin% t_labels) | (t_n >= 2L & nk_n < 2L)
  dual <- myeloid_ok & t_ok
  myeloid_ok[dual] <- FALSE; t_ok[dual] <- FALSE
  eligible <- core$raw
  eligible[, c("CD14HI_MONO", "CD16HI_MONO", "MACROPHAGE", "OSTEOCLAST")] <-
    eligible[, c("CD14HI_MONO", "CD16HI_MONO", "MACROPHAGE", "OSTEOCLAST"), drop = FALSE] & myeloid_ok
  eligible[, c("CD4_TREG", "CD8_TEX", "CD8_PTEX")] <-
    eligible[, c("CD4_TREG", "CD8_TEX", "CD8_PTEX"), drop = FALSE] & t_ok
  scores <- score_ucell(mat)
  assignment <- rep("unassigned", ncol(mat))
  for (lineage_states in list(c("CD14HI_MONO", "CD16HI_MONO", "MACROPHAGE", "OSTEOCLAST"), c("CD4_TREG", "CD8_TEX"))) {
    idx <- which(rowSums(eligible[, lineage_states, drop = FALSE]) > 0)
    for (j in idx) {
      elig <- lineage_states[eligible[j, lineage_states]]
      winner <- elig[which.max(scores[j, elig])]
      competitor <- setdiff(lineage_states, winner)
      margin <- scores[j, winner] - max(scores[j, competitor])
      cutoff <- if (winner %chin% c("CD4_TREG", "CD8_TEX")) t_margin else myeloid_margin
      if (is.finite(margin) && margin >= cutoff) assignment[j] <- winner
    }
  }
  cells <- data.table(
    cell_id = cell_ids, dataset = "GSE225209", sample_or_lesion = sample_name,
    patient_id = patient_id, cancer_code = "LUCA", input_type = "raw_10x_counts",
    nCount_RNA = n_count[qc], nFeature_RNA = n_feature[qc], pct_mito = pct_mito[qc],
    broad_lineage = fifelse(dual, "ambiguous", fifelse(myeloid_ok, "myeloid", fifelse(t_ok, "t_cell", "other_unassigned"))),
    final_assignment = assignment
  )
  audit <- data.table(
    dataset = "GSE225209", sample_or_lesion = sample_name, patient_id = patient_id,
    input_type = "raw_10x_counts", raw_or_author_cells = ncol(gex), qc_cells = ncol(mat),
    assigned_cells = sum(assignment != "unassigned"), assigned_fraction = mean(assignment != "unassigned"),
    dual_lineage_cells = sum(dual)
  )
  log_msg("DONE lung sample=", sample_name, " qc=", ncol(mat), " assigned=", audit$assigned_cells)
  rm(raw, gex, mat, scores); gc(FALSE)
  list(cells = cells, audit = audit)
}

read_selected_matrix <- function(name) {
  path <- file.path(breast_selected_dir, name)
  x <- fread(path)
  genes <- x[[1L]]
  mat <- as.matrix(x[, -1L, with = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- genes
  mat
}

process_breast_lesion <- function(mat_all, metadata_path, lesion, patient_id, author_mode) {
  meta <- fread(metadata_path)
  id_col <- if ("Cell" %chin% names(meta)) "Cell" else "cell.ID"
  setkeyv(meta, id_col)
  all_cells <- colnames(mat_all)
  lesion_cells <- all_cells[startsWith(all_cells, paste0(lesion, "_"))]
  meta_lesion <- meta[J(lesion_cells)]
  if (nrow(meta_lesion) != length(lesion_cells) || anyNA(meta_lesion[[id_col]])) stop("Breast metadata alignment failed: ", lesion, call. = FALSE)
  log_msg("START breast lesion=", lesion, " patient=", patient_id)
  qc <- meta_lesion$nFeature_RNA >= 200 & meta_lesion$nFeature_RNA <= 8000 &
    meta_lesion$nCount_RNA >= 500 & meta_lesion$percent.mt <= 20
  meta_lesion <- meta_lesion[qc]
  mat <- mat_all[, lesion_cells[qc], drop = FALSE]
  core <- make_core_eligibility(mat)
  t_n <- detected_count(mat, t_fallback)
  nk_n <- detected_count(mat, nk_exclusion)
  if (author_mode == "BoM12") {
    author_label <- meta_lesion$cellType
    myeloid_ok <- author_label %chin% c("Monocyte", "Macrophage", "Osteoclast")
    t_ok <- author_label == "T"
  } else {
    author_label <- meta_lesion$cell_type3
    myeloid_ok <- author_label %chin% c("Monocyte", "Macrophage", "Osteoclast")
    t_ok <- author_label == "Other.immune" & t_n >= 2L & nk_n < 2L
  }
  dual <- myeloid_ok & t_ok
  myeloid_ok[dual] <- FALSE; t_ok[dual] <- FALSE
  eligible <- core$raw[, primary_states, drop = FALSE]
  eligible[, c("CD14HI_MONO", "CD16HI_MONO", "MACROPHAGE", "OSTEOCLAST")] <-
    eligible[, c("CD14HI_MONO", "CD16HI_MONO", "MACROPHAGE", "OSTEOCLAST"), drop = FALSE] & myeloid_ok
  eligible[, c("CD4_TREG", "CD8_TEX")] <-
    eligible[, c("CD4_TREG", "CD8_TEX"), drop = FALSE] & t_ok
  assignment <- rep("unassigned", ncol(mat))
  for (lineage_states in list(c("CD14HI_MONO", "CD16HI_MONO", "MACROPHAGE", "OSTEOCLAST"), c("CD4_TREG", "CD8_TEX"))) {
    idx <- which(rowSums(eligible[, lineage_states, drop = FALSE]) == 1L)
    if (length(idx)) assignment[idx] <- lineage_states[max.col(eligible[idx, lineage_states, drop = FALSE], ties.method = "first")]
  }
  cells <- data.table(
    cell_id = colnames(mat), dataset = "GSE190772", sample_or_lesion = lesion,
    patient_id = patient_id, cancer_code = "BRCA", input_type = "author_log_normalized",
    nCount_RNA = meta_lesion$nCount_RNA, nFeature_RNA = meta_lesion$nFeature_RNA,
    pct_mito = meta_lesion$percent.mt, broad_lineage = fifelse(dual, "ambiguous", fifelse(myeloid_ok, "myeloid", fifelse(t_ok, "t_cell", "other_unassigned"))),
    final_assignment = assignment
  )
  audit <- data.table(
    dataset = "GSE190772", sample_or_lesion = lesion, patient_id = patient_id,
    input_type = "author_log_normalized", raw_or_author_cells = length(lesion_cells), qc_cells = ncol(mat),
    assigned_cells = sum(assignment != "unassigned"), assigned_fraction = mean(assignment != "unassigned"),
    dual_lineage_cells = sum(dual)
  )
  log_msg("DONE breast lesion=", lesion, " qc=", ncol(mat), " assigned=", audit$assigned_cells)
  list(cells = cells, audit = audit)
}

log_msg("Gate10S full start seed=20260807 workers=8")
lung_spec <- data.table(
  gsm = c("GSM7041480", "GSM7041481", "GSM7041482"),
  stem = c("sg1", "sg2", "sg3"), sample = c("sz", "s13", "s14"),
  patient = c("LUCA_BM_01", "LUCA_BM_02", "LUCA_BM_03")
)
lung_results <- lapply(seq_len(nrow(lung_spec)), function(i) {
  process_lung(lung_spec$gsm[[i]], lung_spec$stem[[i]], lung_spec$sample[[i]], lung_spec$patient[[i]])
})

bom12 <- read_selected_matrix("GSE190772_BoM_logCounts.selected_markers.tsv.gz")
bom7 <- read_selected_matrix("GSM6870693_BoM7_scRNA_LogCounts.selected_markers.tsv.gz")
bom8 <- read_selected_matrix("GSM6870694_BoM8_scRNA_LogCounts.selected_markers.tsv.gz")
breast_results <- list(
  process_breast_lesion(bom12, file.path(project, "data/raw/GSE190772/GSE190772_BoM_MetaData.txt.gz"), "BoM1", "BRCA_BM_01", "BoM12"),
  process_breast_lesion(bom12, file.path(project, "data/raw/GSE190772/GSE190772_BoM_MetaData.txt.gz"), "BoM2", "BRCA_BM_01", "BoM12"),
  process_breast_lesion(bom7, file.path(project, "data/raw/GSE190772/GSM6870693_BoM7_scRNA_MetaData.txt.gz"), "BoM7", "BRCA_BM_02", "BoM78"),
  process_breast_lesion(bom8, file.path(project, "data/raw/GSE190772/GSM6870694_BoM8_scRNA_MetaData.txt.gz"), "BoM8", "BRCA_BM_02", "BoM78")
)
all_results <- c(lung_results, breast_results)
cells <- rbindlist(lapply(all_results, `[[`, "cells"), use.names = TRUE, fill = TRUE)
audit <- rbindlist(lapply(all_results, `[[`, "audit"), use.names = TRUE, fill = TRUE)
audit[, S3_ge_500_qc := qc_cells >= 500]
audit[, S4_ge_50_assigned := assigned_cells >= 50]
audit[, projection_status := fifelse(!S3_ge_500_qc, "NOT_PROJECTABLE", fifelse(!S4_ge_50_assigned, "LOW_COVERAGE", "PROJECTABLE"))]

sample_counts_long <- cells[final_assignment %chin% primary_states, .N,
  by = .(dataset, sample_or_lesion, patient_id, cancer_code, input_type, state_id = final_assignment)]
sample_base <- unique(cells[, .(dataset, sample_or_lesion, patient_id, cancer_code, input_type)])
grid <- sample_base[, .(state_id = primary_states),
  by = .(dataset, sample_or_lesion, patient_id, cancer_code, input_type)]
sample_counts_long <- merge(grid, sample_counts_long,
  by = c("dataset", "sample_or_lesion", "patient_id", "cancer_code", "input_type", "state_id"), all.x = TRUE)
sample_counts_long[is.na(N), N := 0L]
sample_counts <- dcast(sample_counts_long, dataset + sample_or_lesion + patient_id + cancer_code + input_type ~ state_id, value.var = "N")
sample_counts <- merge(sample_counts, audit[, .(dataset, sample_or_lesion, qc_cells, assigned_cells, assigned_fraction, S3_ge_500_qc, S4_ge_50_assigned, projection_status)],
                       by = c("dataset", "sample_or_lesion"), all.x = TRUE)

transform_rows <- function(dt) {
  z <- copy(dt)
  for (state in primary_states) {
    p <- (z[[state]] + 0.5) / (z$qc_cells + 1)
    logit <- log(p / (1 - p))
    z[[state]] <- (logit - source_mean[[state]]) / source_sd[[state]]
  }
  z
}
add_axes <- function(dt) {
  z <- copy(dt)
  z[, Mphi_OC_axis := (MACROPHAGE + OSTEOCLAST) / 2]
  z[, Mono_axis := (CD14HI_MONO + CD16HI_MONO) / 2]
  z[, Treg_Tex_axis := (CD4_TREG + CD8_TEX) / 2]
  z
}

sample_z <- transform_rows(sample_counts)
sample_axes <- add_axes(sample_z)

lung_patient <- sample_z[dataset == "GSE225209", c("dataset", "patient_id", "cancer_code", "input_type", primary_states), with = FALSE]
lung_patient[, aggregation := "single_lesion"]
breast_equal <- sample_z[dataset == "GSE190772", lapply(.SD, mean), by = .(dataset, patient_id, cancer_code, input_type), .SDcols = primary_states]
breast_equal[, aggregation := "equal_lesion_mean_z"]
patient_z_primary <- rbindlist(list(lung_patient, breast_equal), use.names = TRUE, fill = TRUE)
patient_axes_primary <- add_axes(patient_z_primary)

breast_pooled_counts <- sample_counts[dataset == "GSE190772", c(lapply(.SD, sum), list(qc_cells = sum(qc_cells))),
  by = .(dataset, patient_id, cancer_code, input_type), .SDcols = primary_states]
breast_pooled_z <- transform_rows(breast_pooled_counts)
breast_pooled_z[, aggregation := "pooled_lesion_cell_counts_sensitivity"]
breast_pooled_axes <- add_axes(breast_pooled_z)

smoke_counts <- read_tsv(file.path(project, "results/gate10s_supportive_projection/smoke_v1/smoke_state_counts.tsv"))
full_smoke_subset <- sample_counts_long[
  (dataset == "GSE225209" & sample_or_lesion == "sz") |
  (dataset == "GSE190772" & sample_or_lesion == "BoM1"),
  .(dataset, sample_or_lesion, patient_id, final_assignment = state_id, N)
][N > 0L]
full_unassigned <- audit[
  (dataset == "GSE225209" & sample_or_lesion == "sz") |
  (dataset == "GSE190772" & sample_or_lesion == "BoM1"),
  .(dataset, sample_or_lesion, patient_id, final_assignment = "unassigned", N = qc_cells - assigned_cells)
]
full_smoke_subset <- rbind(full_smoke_subset, full_unassigned)
compare_keys <- c("dataset", "sample_or_lesion", "patient_id", "final_assignment")
smoke_compare <- merge(smoke_counts, full_smoke_subset, by = compare_keys, suffixes = c("_smoke", "_full"), all = TRUE)
smoke_reproduced <- nrow(smoke_compare) == nrow(smoke_counts) &&
  nrow(smoke_compare) == nrow(full_smoke_subset) &&
  !anyNA(smoke_compare) && all(smoke_compare$N_smoke == smoke_compare$N_full)

fwrite(audit, file.path(out_dir, "sample_projection_qc.tsv"), sep = "\t")
fwrite(sample_counts, file.path(out_dir, "sample_state_counts.tsv"), sep = "\t")
fwrite(sample_z, file.path(out_dir, "sample_state_z.tsv"), sep = "\t")
fwrite(sample_axes[, c("dataset", "sample_or_lesion", "patient_id", "cancer_code", "input_type", "Mphi_OC_axis", "Mono_axis", "Treg_Tex_axis", "projection_status"), with = FALSE],
       file.path(out_dir, "sample_axes.tsv"), sep = "\t")
fwrite(patient_axes_primary, file.path(out_dir, "patient_axes_primary.tsv"), sep = "\t")
fwrite(breast_pooled_axes, file.path(out_dir, "breast_patient_axes_pooled_sensitivity.tsv"), sep = "\t")
fwrite(cells, file.path(out_dir, "gate10s_cell_assignments.tsv.gz"), sep = "\t")
fwrite(data.table(check = "smoke_sample_counts_exact_match", pass = smoke_reproduced), file.path(out_dir, "smoke_reproduction_audit.tsv"), sep = "\t")

five_unique <- uniqueN(patient_axes_primary$patient_id) == 5L && nrow(patient_axes_primary) == 5L
seven_samples <- nrow(audit) == 7L && uniqueN(paste(audit$dataset, audit$sample_or_lesion)) == 7L
if (!all(audit$S3_ge_500_qc)) decision <- "NOT_PROJECTABLE" else if (!all(audit$S4_ge_50_assigned)) decision <- "PROJECTABLE_WITH_LOW_COVERAGE" else decision <- "PROJECTABLE"
hard_pass <- five_unique && seven_samples && smoke_reproduced
if (!hard_pass) decision <- "HARD_AUDIT_FAILURE"

writeLines(c(
  "# Gate10S full supportive projection final audit", "",
  "## Material Passport", "",
  "- Origin Skill: academic-research-suite / experiment-agent",
  "- Version: `gate10s_supportive_projection_v1`",
  "- Verification Status: ANALYZED; deterministic smoke subset reproduced",
  "", "## Decision", "",
  paste0("- Projection decision: **", decision, "**"),
  paste0("- Lung patients: ", uniqueN(audit[dataset == "GSE225209", patient_id])),
  paste0("- Breast lesions/patients: ", audit[dataset == "GSE190772", .N], "/", uniqueN(audit[dataset == "GSE190772", patient_id])),
  paste0("- Unique patient outputs: ", nrow(patient_axes_primary)),
  paste0("- S3 samples with >=500 QC cells: ", sum(audit$S3_ge_500_qc), "/", nrow(audit)),
  paste0("- S4 samples with >=50 assigned cells: ", sum(audit$S4_ge_50_assigned), "/", nrow(audit)),
  paste0("- Smoke subset exact reproduction: ", if (smoke_reproduced) "PASS" else "FAIL"),
  "- GSE190772 output label: `NORMALIZED_SUPPORTIVE_ONLY`",
  "- P values calculated: **NO**",
  "- Datasets pooled for inference: **NO**",
  "- Gate10A decision changed: **NO**",
  "- Gate10B authorized: **NO**",
  "", "## Interpretation boundary", "",
  "These outputs test scoring portability in small independent datasets. They cannot confirm a pan-cancer association, estimate cancer-origin effects, or rescue Gate10A/Gate10B. Breast primary values are equal-lesion means; cell-pooled values are sensitivity analyses only.", ""
), file.path(out_dir, "GATE10S_FINAL_AUDIT.md"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
pkg <- data.table(package = c("R", "Matrix", "data.table", "SingleCellExperiment", "scuttle", "SingleR", "UCell"),
                  version = c(as.character(getRversion()), vapply(c("Matrix", "data.table", "SingleCellExperiment", "scuttle", "SingleR", "UCell"), function(p) as.character(packageVersion(p)), character(1))))
fwrite(pkg, file.path(out_dir, "package_versions.tsv"), sep = "\t")
log_msg("Gate10S full completed decision=", decision, " patients=", nrow(patient_axes_primary), " smoke_reproduced=", smoke_reproduced)
close(log_con)

manifest_files <- list.files(out_dir, full.names = TRUE)
manifest_files <- manifest_files[basename(manifest_files) != "SHA256SUMS"]
status <- system(sprintf("cd %s && sha256sum %s > SHA256SUMS", shQuote(out_dir), paste(shQuote(basename(manifest_files)), collapse = " ")))
if (status != 0L) stop("SHA256 manifest failed", call. = FALSE)
cat("GATE10S_FULL_DECISION=", decision, "\n", sep = "")
