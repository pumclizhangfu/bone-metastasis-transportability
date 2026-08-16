#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(SingleR)
  library(BiocParallel)
  library(parallel)
})

options(error = function() {
  traceback(3L)
  quit(status = 1L, save = "no")
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) {
  stop(paste(
    "Usage: run_gate12z_gse266330_axis1.R <raw_flat_dir>",
    "<library_crosswalk.tsv> <patient_crosswalk.tsv> <gate9a_cells.tsv.gz>",
    "<transfer_reference.rds> <gate12g_model.rds> <outdir> <workers> <n_subsample>"
  ))
}

raw_dir <- normalizePath(args[[1L]], mustWork = TRUE)
library_crosswalk_file <- normalizePath(args[[2L]], mustWork = TRUE)
patient_crosswalk_file <- normalizePath(args[[3L]], mustWork = TRUE)
gate9a_cells_file <- normalizePath(args[[4L]], mustWork = TRUE)
transfer_file <- normalizePath(args[[5L]], mustWork = TRUE)
model_file <- normalizePath(args[[6L]], mustWork = TRUE)
outdir <- args[[7L]]
workers <- as.integer(args[[8L]])
n_subsample <- as.integer(args[[9L]])
if (!is.finite(workers) || workers < 1L) stop("workers must be >=1")
if (!is.finite(n_subsample) || n_subsample < 1L) stop("n_subsample must be >=1")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
cache_dir <- file.path(outdir, "library_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

seed <- 1208L
n_permutation <- 9999L
n_boot <- 1999L
subsample_fraction <- 0.80
set.seed(seed)

message("Gate12Z GSE266330 Axis1 start; workers=", workers,
        "; n_subsample=", n_subsample, "; seed=", seed)

transfer <- readRDS(transfer_file)
model <- readRDS(model_file)
if (!identical(transfer$status, "PASS")) stop("Frozen Gate12G transfer reference must have PASS status")
if (model$selected_k != 4L) stop("Frozen Gate12G model must select K=4")
if (ncol(model$x_raw_clr) != 26L) stop("Frozen Gate12G model must contain 26 features")
if (!identical(names(model$blocks), colnames(model$x_raw_clr))) stop("Frozen feature and block order mismatch")

library_crosswalk <- fread(library_crosswalk_file)
patient_crosswalk <- fread(patient_crosswalk_file)
cells_qc <- fread(gate9a_cells_file,
                  select = c("cell_id", "cell_barcode", "library_id", "patient_id",
                             "cancer", "analysis_role"))
if (nrow(library_crosswalk) != 63L || uniqueN(patient_crosswalk$patient_id) != 47L) {
  stop("Frozen GSE266330 scope must be 63 libraries and 47 patients/donors")
}
if (nrow(cells_qc) != 195934L) stop("Frozen Gate9A QC cell universe must contain 195,934 cells")
if (uniqueN(cells_qc$patient_id) != 47L) stop("Gate9A QC cells must contain 47 patients/donors")

aggregate_symbols <- function(mat, symbols) {
  keep <- !is.na(symbols) & nzchar(symbols)
  mat <- mat[keep, , drop = FALSE]
  symbols <- symbols[keep]
  u <- unique(symbols)
  if (length(u) == length(symbols)) {
    rownames(mat) <- u
    return(as(mat, "dgCMatrix"))
  }
  map <- sparseMatrix(i = match(symbols, u), j = seq_along(symbols), x = 1,
                      dims = c(length(u), length(symbols)))
  ans <- map %*% mat
  rownames(ans) <- u
  colnames(ans) <- colnames(mat)
  as(ans, "dgCMatrix")
}

make_test_logcpm <- function(counts, genes) {
  lib <- Matrix::colSums(counts)
  if (any(lib <= 0)) stop("QC matrix contains zero-count cells")
  idx <- match(genes, rownames(counts), nomatch = 0L)
  out <- matrix(0, nrow = length(genes), ncol = ncol(counts),
                dimnames = list(genes, colnames(counts)))
  present <- which(idx > 0L)
  if (length(present)) out[present, ] <- as.matrix(counts[idx[present], , drop = FALSE])
  log1p(t(t(out) / lib) * 1e6)
}

classify_taxonomy <- function(counts, ref_obj) {
  test <- make_test_logcpm(counts, ref_obj$marker_genes)
  pred <- SingleR(
    test = test, ref = ref_obj$logcpm, labels = ref_obj$meta$label,
    de.method = "classic", fine.tune = TRUE, prune = TRUE,
    BPPARAM = SerialParam()
  )
  data.table(
    cell_id = colnames(counts),
    label = as.character(pred$labels),
    pruned = as.character(pred$pruned.labels),
    delta_next = as.numeric(pred$delta.next)
  )
}

process_library <- function(i) {
  row <- library_crosswalk[i]
  lib_id <- row$geo_prefix
  cache_file <- file.path(cache_dir, paste0(lib_id, ".rds"))
  if (file.exists(cache_file)) {
    ans <- readRDS(cache_file)
    message("CACHED ", lib_id, " cells=", nrow(ans))
    return(ans)
  }
  files <- c(
    matrix = file.path(raw_dir, paste0(lib_id, "_matrix.mtx.gz")),
    features = file.path(raw_dir, paste0(lib_id, "_features.tsv.gz")),
    barcodes = file.path(raw_dir, paste0(lib_id, "_barcodes.tsv.gz"))
  )
  if (!all(file.exists(files))) stop("Incomplete raw matrix triplet for ", lib_id)
  feat <- fread(files[["features"]], header = FALSE)
  if (ncol(feat) == 2L) feat[, V3 := "Gene Expression"]
  setnames(feat, names(feat)[1:3], c("feature_id", "gene_symbol", "feature_type"))
  bc <- fread(files[["barcodes"]], header = FALSE)[[1L]]
  raw <- as(Matrix::readMM(gzfile(files[["matrix"]])), "dgCMatrix")
  colnames(raw) <- bc
  if (nrow(raw) != nrow(feat) || ncol(raw) != length(bc)) stop("Raw dimensions mismatch for ", lib_id)
  gex_idx <- which(feat$feature_type == "Gene Expression")
  gex <- aggregate_symbols(raw[gex_idx, , drop = FALSE], feat$gene_symbol[gex_idx])

  old <- cells_qc[library_id == lib_id]
  if (!nrow(old)) stop("No frozen QC cells for ", lib_id)
  idx <- match(old$cell_barcode, colnames(gex))
  if (anyNA(idx)) stop("Frozen Gate9A QC barcode missing from raw matrix for ", lib_id)
  gex <- gex[, idx, drop = FALSE]
  colnames(gex) <- old$cell_id

  broad_pred <- classify_taxonomy(gex, transfer$references$Broad)
  broad_pred[, broad := fifelse(is.na(pruned) | !nzchar(pruned), "Unassigned", pruned)]
  state <- rep(NA_character_, ncol(gex))
  state_delta <- rep(NA_real_, ncol(gex))
  for (lineage in c("Myeloid", "T_NK")) {
    idx_lineage <- which(broad_pred$broad == lineage)
    if (!length(idx_lineage)) next
    state_pred <- classify_taxonomy(gex[, idx_lineage, drop = FALSE], transfer$references[[lineage]])
    unresolved <- if (lineage == "Myeloid") "Unresolved_myeloid" else "Unresolved_T_NK"
    state_value <- fifelse(is.na(state_pred$pruned) | !nzchar(state_pred$pruned),
                           unresolved, state_pred$pruned)
    state[idx_lineage] <- state_value
    state_delta[idx_lineage] <- state_pred$delta_next
  }
  ans <- data.table(
    cell_id = old$cell_id, cell_barcode = old$cell_barcode, library_id = lib_id,
    patient_id = old$patient_id, cancer = old$cancer, analysis_role = old$analysis_role,
    gate12z_broad = broad_pred$broad, gate12z_broad_delta = broad_pred$delta_next,
    gate12z_state = state, gate12z_state_delta = state_delta
  )
  saveRDS(ans, cache_file, compress = "xz")
  message("COMPLETED ", lib_id, " cells=", nrow(ans),
          " broad_assigned=", sprintf("%.3f", mean(ans$gate12z_broad != "Unassigned")))
  ans
}

idx <- seq_len(nrow(library_crosswalk))
assignments <- if (.Platform$OS.type == "unix" && workers > 1L) {
  mclapply(idx, process_library, mc.cores = workers, mc.preschedule = FALSE)
} else {
  lapply(idx, process_library)
}
if (any(vapply(assignments, inherits, logical(1), what = "try-error"))) {
  stop("At least one library classification failed")
}
cells <- rbindlist(assignments, use.names = TRUE)
if (nrow(cells) != nrow(cells_qc)) stop("Classified cell count differs from frozen QC universe")
fwrite(cells, file.path(outdir, "gse266330_gate12z_cell_assignments.tsv.gz"), sep = "\t")

broad_levels <- sub("^Broad__", "", colnames(model$x_raw_clr)[model$blocks == "Broad"])
my_levels <- sub("^Myeloid__", "", colnames(model$x_raw_clr)[model$blocks == "Myeloid"])
tn_levels <- sub("^T_NK__", "", colnames(model$x_raw_clr)[model$blocks == "T_NK"])

make_count_wide <- function(dt, label_col, levels) {
  base <- unique(patient_crosswalk[, .(
    patient_id, cancer = cancer_harmonized, analysis_role,
    condition = fifelse(analysis_role == "negative_control", "healthy_marrow", "bone_metastasis")
  )])
  obs <- dt[!is.na(get(label_col)), .N, by = .(patient_id, label = get(label_col))]
  grid <- base[, .(label = levels), by = .(patient_id, cancer, analysis_role, condition)]
  z <- merge(grid, obs, by = c("patient_id", "label"), all.x = TRUE)
  z[is.na(N), N := 0L]
  wide <- dcast(z, patient_id + cancer + analysis_role + condition ~ label,
                value.var = "N", fill = 0)
  for (lev in levels) if (!lev %chin% names(wide)) wide[, (lev) := 0L]
  setcolorder(wide, c("patient_id", "cancer", "analysis_role", "condition", levels))
  wide
}

broad_w <- make_count_wide(cells, "gate12z_broad", broad_levels)
my_w <- make_count_wide(cells[gate12z_broad == "Myeloid"], "gate12z_state", my_levels)
tn_w <- make_count_wide(cells[gate12z_broad == "T_NK"], "gate12z_state", tn_levels)
setkey(broad_w, patient_id); setkey(my_w, patient_id); setkey(tn_w, patient_id)
if (!identical(broad_w$patient_id, my_w$patient_id) || !identical(broad_w$patient_id, tn_w$patient_id)) {
  stop("Patient composition blocks are misaligned")
}

close_clr <- function(m) {
  m <- as.matrix(m)
  p <- m + 0.5
  p <- p / rowSums(p)
  lp <- log(p)
  lp - rowMeans(lp)
}

project_counts <- function(bmat, mmat, tmat) {
  xb <- close_clr(bmat); xm <- close_clr(mmat); xt <- close_clr(tmat)
  colnames(xb) <- paste0("Broad__", colnames(bmat))
  colnames(xm) <- paste0("Myeloid__", colnames(mmat))
  colnames(xt) <- paste0("T_NK__", colnames(tmat))
  x <- cbind(xb, xm, xt)
  x <- x[, colnames(model$x_raw_clr), drop = FALSE]
  for (block in unique(model$blocks)) {
    j <- which(model$blocks == block)
    block_scale <- sqrt(sum(apply(model$x_raw_clr[, j, drop = FALSE], 2L, var)))
    x[, j] <- x[, j, drop = FALSE] / block_scale
  }
  centered <- sweep(x, 2L, model$pca$center, "-")
  score <- centered %*% model$pca$rotation[, seq_len(model$selected_k), drop = FALSE]
  colnames(score) <- paste0("Axis", seq_len(ncol(score)))
  list(score = score, transformed = x)
}

bmat <- as.matrix(broad_w[, ..broad_levels])
mmat <- as.matrix(my_w[, ..my_levels])
tmat <- as.matrix(tn_w[, ..tn_levels])
projection <- project_counts(bmat, mmat, tmat)
full_scores <- projection$score
all_finite <- apply(projection$transformed, 1L, function(z) all(is.finite(z)))

patient_qc <- broad_w[, .(patient_id, cancer, analysis_role, condition)]
patient_qc[, `:=`(
  qc_cells = rowSums(bmat),
  broad_assigned_fraction = 1 - bmat[, "Unassigned"] / rowSums(bmat),
  myeloid_cells = rowSums(mmat),
  t_nk_cells = rowSums(tmat),
  all_26_features_representable = TRUE,
  all_transforms_finite = all_finite
)]
patient_qc[, projectable := qc_cells >= 500L & myeloid_cells >= 50L &
             t_nk_cells >= 50L & all_26_features_representable & all_transforms_finite]
score_dt <- cbind(patient_qc, as.data.table(full_scores))
fwrite(score_dt, file.path(outdir, "gse266330_axis_scores.tsv"), sep = "\t")

class_metrics <- as.data.table(transfer$class_metrics)
feature_reliability <- data.table(feature = colnames(model$x_raw_clr), block = unname(model$blocks))
feature_reliability[, label := sub("^(Broad|Myeloid|T_NK)__", "", feature)]
feature_reliability[, taxonomy := block]
feature_reliability[, rejection_bin := label %chin% c("Unassigned", "Unresolved_myeloid", "Unresolved_T_NK")]
feature_reliability <- merge(
  feature_reliability,
  class_metrics[, .(taxonomy, label = truth, evaluable_patients, recall)],
  by = c("taxonomy", "label"), all.x = TRUE, sort = FALSE
)
feature_reliability[, reliable := rejection_bin |
                      (!is.na(recall) & evaluable_patients >= 3L & recall >= 0.50)]
feature_reliability <- feature_reliability[match(colnames(model$x_raw_clr), feature)]
loading <- model$pca$rotation[, seq_len(model$selected_k), drop = FALSE]
loading_coverage <- rbindlist(lapply(seq_len(model$selected_k), function(a) {
  data.table(
    axis = paste0("Axis", a),
    reliable_loading_coverage = sum(loading[feature_reliability$reliable, a]^2) / sum(loading[, a]^2)
  )
}))
fwrite(feature_reliability, file.path(outdir, "gse266330_feature_reliability.tsv"), sep = "\t")
fwrite(loading_coverage, file.path(outdir, "gse266330_loading_coverage.tsv"), sep = "\t")

# Frozen sensitivity: collapse unreliable resolved labels into lineage rejection bins.
bmat_s <- bmat; mmat_s <- mmat; tmat_s <- tmat
unreliable <- feature_reliability[reliable == FALSE]
for (i in seq_len(nrow(unreliable))) {
  lab <- unreliable$label[i]; tax <- unreliable$taxonomy[i]
  if (tax == "Myeloid" && lab %chin% colnames(mmat_s)) {
    mmat_s[, "Unresolved_myeloid"] <- mmat_s[, "Unresolved_myeloid"] + mmat_s[, lab]
    mmat_s[, lab] <- 0
  } else if (tax == "T_NK" && lab %chin% colnames(tmat_s)) {
    tmat_s[, "Unresolved_T_NK"] <- tmat_s[, "Unresolved_T_NK"] + tmat_s[, lab]
    tmat_s[, lab] <- 0
  } else if (tax == "Broad" && lab %chin% colnames(bmat_s)) {
    bmat_s[, "Unassigned"] <- bmat_s[, "Unassigned"] + bmat_s[, lab]
    bmat_s[, lab] <- 0
  }
}
sensitivity_scores <- project_counts(bmat_s, mmat_s, tmat_s)$score
sensitivity_dt <- data.table(
  axis = paste0("Axis", seq_len(model$selected_k)),
  score_spearman = vapply(seq_len(model$selected_k), function(a) {
    cor(full_scores[, a], sensitivity_scores[, a], method = "spearman")
  }, numeric(1))
)
fwrite(sensitivity_dt, file.path(outdir, "gse266330_unreliable_state_sensitivity.tsv"), sep = "\t")

hl_shift <- function(x, y) median(as.vector(outer(x, y, "-")))
rank_biserial <- function(x, y) {
  n1 <- length(x); n0 <- length(y)
  ranks <- rank(c(x, y), ties.method = "average")
  u <- sum(ranks[seq_len(n1)]) - n1 * (n1 + 1) / 2
  2 * u / (n1 * n0) - 1
}

primary <- score_dt[projectable == TRUE]
bm <- primary[condition == "bone_metastasis", Axis1]
ctrl <- primary[condition == "healthy_marrow", Axis1]
origins_ge2 <- primary[condition == "bone_metastasis", .N, by = cancer][N >= 2L]
dataset_evaluable <- length(bm) >= 20L && length(ctrl) >= 3L && nrow(origins_ge2) >= 4L
observed_hl <- if (length(bm) && length(ctrl)) hl_shift(bm, ctrl) else NA_real_
wilcox_p <- if (length(bm) && length(ctrl)) wilcox.test(bm, ctrl, exact = FALSE)$p.value else NA_real_
set.seed(seed + 100L)
all_axis <- c(bm, ctrl)
n_bm <- length(bm)
perm_hl <- if (length(bm) && length(ctrl)) replicate(n_permutation, {
  z <- sample(all_axis, replace = FALSE)
  hl_shift(z[seq_len(n_bm)], z[-seq_len(n_bm)])
}) else rep(NA_real_, n_permutation)
permutation_p <- if (all(is.finite(perm_hl))) (1 + sum(abs(perm_hl) >= abs(observed_hl))) /
  (n_permutation + 1) else NA_real_
primary_contrast <- data.table(
  endpoint = "Axis1_bone_metastasis_minus_healthy_marrow",
  n_bone_metastasis = length(bm), n_healthy_marrow = length(ctrl),
  median_bone_metastasis = if (length(bm)) median(bm) else NA_real_,
  median_healthy_marrow = if (length(ctrl)) median(ctrl) else NA_real_,
  hodges_lehmann_shift = observed_hl,
  rank_biserial = if (length(bm) && length(ctrl)) rank_biserial(bm, ctrl) else NA_real_,
  wilcoxon_two_sided_p = wilcox_p,
  label_permutation_two_sided_p = permutation_p,
  directional_expectation = "positive",
  direction_positive = is.finite(observed_hl) && observed_hl > 0,
  dataset_evaluable = dataset_evaluable,
  axis1_loading_coverage_ge_0_70 = loading_coverage[axis == "Axis1", reliable_loading_coverage] >= 0.70,
  directional_support = dataset_evaluable && is.finite(observed_hl) && observed_hl > 0 &&
    is.finite(permutation_p) && permutation_p < 0.05
)
fwrite(primary_contrast, file.path(outdir, "gse266330_primary_contrast.tsv"), sep = "\t")

cancer_summary <- primary[condition == "bone_metastasis", .(
  n_patients = .N, median_Axis1 = median(Axis1), q1_Axis1 = quantile(Axis1, 0.25),
  q3_Axis1 = quantile(Axis1, 0.75)
), by = cancer][order(cancer)]
fwrite(cancer_summary, file.path(outdir, "gse266330_cancer_origin_summary.tsv"), sep = "\t")

loco <- rbindlist(lapply(origins_ge2$cancer, function(drop_origin) {
  b <- primary[condition == "bone_metastasis" & cancer != drop_origin, Axis1]
  data.table(dropped_origin = drop_origin, n_bm = length(b), n_control = length(ctrl),
             hodges_lehmann_shift = hl_shift(b, ctrl), direction_positive = hl_shift(b, ctrl) > 0)
}))
fwrite(loco, file.path(outdir, "gse266330_leave_one_origin_out.tsv"), sep = "\t")

draw_multihyper <- function(counts, fraction = 0.8) {
  counts <- as.integer(counts)
  draw_total <- floor(sum(counts) * fraction)
  ans <- integer(length(counts))
  remaining_total <- sum(counts); remaining_draw <- draw_total
  if (length(counts) > 1L) {
    for (i in seq_len(length(counts) - 1L)) {
      ans[i] <- rhyper(1L, counts[i], remaining_total - counts[i], remaining_draw)
      remaining_total <- remaining_total - counts[i]
      remaining_draw <- remaining_draw - ans[i]
    }
  }
  ans[length(counts)] <- remaining_draw
  ans
}

projectable_idx <- which(patient_qc$projectable)
subsample_rows <- vector("list", n_subsample)
set.seed(seed + 200L)
for (b in seq_len(n_subsample)) {
  bs_b <- t(apply(bmat, 1L, draw_multihyper, fraction = subsample_fraction))
  bs_m <- t(apply(mmat, 1L, draw_multihyper, fraction = subsample_fraction))
  bs_t <- t(apply(tmat, 1L, draw_multihyper, fraction = subsample_fraction))
  colnames(bs_b) <- colnames(bmat)
  colnames(bs_m) <- colnames(mmat)
  colnames(bs_t) <- colnames(tmat)
  bs_score <- project_counts(bs_b, bs_m, bs_t)$score[, 1L]
  bmv <- bs_score[patient_qc$condition == "bone_metastasis" & patient_qc$projectable]
  cv <- bs_score[patient_qc$condition == "healthy_marrow" & patient_qc$projectable]
  subsample_rows[[b]] <- data.table(
    iteration = b,
    patient_rank_spearman = cor(full_scores[projectable_idx, 1L], bs_score[projectable_idx],
                                method = "spearman"),
    hodges_lehmann_shift = hl_shift(bmv, cv),
    direction_positive = hl_shift(bmv, cv) > 0
  )
}
subsampling <- rbindlist(subsample_rows)
fwrite(subsampling, file.path(outdir, "gse266330_cell_subsampling.tsv.gz"), sep = "\t")

set.seed(seed + 300L)
tech_boot <- rbindlist(lapply(seq_len(n_boot), function(b) {
  idx <- sample(seq_len(nrow(primary)), replace = TRUE)
  z <- primary[idx]
  data.table(
    iteration = b,
    rho_log10_qc_cells = suppressWarnings(cor(z$Axis1, log10(z$qc_cells), method = "spearman")),
    rho_broad_assigned_fraction = suppressWarnings(cor(z$Axis1, z$broad_assigned_fraction,
                                                       method = "spearman"))
  )
}))
tech_point <- data.table(
  metric = c("rho_log10_qc_cells", "rho_broad_assigned_fraction"),
  estimate = c(cor(primary$Axis1, log10(primary$qc_cells), method = "spearman"),
               cor(primary$Axis1, primary$broad_assigned_fraction, method = "spearman")),
  ci_low = c(quantile(tech_boot$rho_log10_qc_cells, 0.025, na.rm = TRUE),
             quantile(tech_boot$rho_broad_assigned_fraction, 0.025, na.rm = TRUE)),
  ci_high = c(quantile(tech_boot$rho_log10_qc_cells, 0.975, na.rm = TRUE),
              quantile(tech_boot$rho_broad_assigned_fraction, 0.975, na.rm = TRUE))
)
fwrite(tech_boot, file.path(outdir, "gse266330_technical_bootstrap.tsv.gz"), sep = "\t")
fwrite(tech_point, file.path(outdir, "gse266330_technical_correlations.tsv"), sep = "\t")

summary <- data.table(
  metric = c(
    "libraries", "qc_cells", "patients_total", "bone_metastasis_total", "healthy_controls_total",
    "projectable_total", "projectable_bone_metastasis", "projectable_healthy_controls",
    "origins_ge2_projectable", "axis1_reliable_loading_coverage", "primary_hl_shift",
    "primary_permutation_p", "subsample_median_rank_spearman", "subsample_q10_rank_spearman",
    "subsample_positive_shift_fraction", "collapsed_state_axis1_spearman", "dataset_evaluable",
    "directional_support"
  ),
  value = as.character(c(
    nrow(library_crosswalk), nrow(cells), nrow(score_dt),
    sum(score_dt$condition == "bone_metastasis"), sum(score_dt$condition == "healthy_marrow"),
    sum(score_dt$projectable), sum(score_dt$projectable & score_dt$condition == "bone_metastasis"),
    sum(score_dt$projectable & score_dt$condition == "healthy_marrow"), nrow(origins_ge2),
    loading_coverage[axis == "Axis1", reliable_loading_coverage], observed_hl, permutation_p,
    median(subsampling$patient_rank_spearman), quantile(subsampling$patient_rank_spearman, 0.10),
    mean(subsampling$direction_positive), sensitivity_dt[axis == "Axis1", score_spearman],
    dataset_evaluable, primary_contrast$directional_support
  ))
)
fwrite(summary, file.path(outdir, "GATE12Z_GSE266330_SUMMARY.tsv"), sep = "\t")

decision <- if (!dataset_evaluable || loading_coverage[axis == "Axis1", reliable_loading_coverage] < 0.70) {
  "NON_PROJECTABLE"
} else if (primary_contrast$directional_support) {
  "PROJECTABLE_DIRECTIONAL_SUPPORT"
} else {
  "PROJECTABLE_NO_DIRECTIONAL_SUPPORT"
}
writeLines(c(
  "# Gate12Z GSE266330 frozen Axis1 checkpoint", "",
  paste0("- Decision: **", decision, "**"),
  paste0("- Frozen QC cells classified: ", format(nrow(cells), big.mark = ",")),
  paste0("- Projectable patients/donors: ", sum(score_dt$projectable), "/", nrow(score_dt),
         " (bone metastasis ", sum(score_dt$projectable & score_dt$condition == "bone_metastasis"),
         "; healthy marrow ", sum(score_dt$projectable & score_dt$condition == "healthy_marrow"), ")"),
  paste0("- Axis1 reliable-loading coverage: ", sprintf("%.4f", loading_coverage[axis == "Axis1", reliable_loading_coverage])),
  paste0("- Hodges-Lehmann metastasis-minus-control shift: ", sprintf("%.4f", observed_hl)),
  paste0("- Two-sided label-permutation P: ", format(permutation_p, digits = 4)),
  paste0("- 80% cell-subsampling positive-shift fraction: ", sprintf("%.4f", mean(subsampling$direction_positive))),
  paste0("- Collapsed-unreliable-state Axis1 correlation: ", sprintf("%.4f", sensitivity_dt[axis == "Axis1", score_spearman])),
  "",
  "The GSE266330 classifier, feature order, rejection bins, transforms, centers, loadings and orientation were reused without refitting. Technical projectability and a group contrast do not establish a universal, bone-specific, mechanistic or clinical state."
), file.path(outdir, "GATE12Z_GSE266330_CHECKPOINT.md"))

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
message("Gate12Z GSE266330 Axis1 complete; decision=", decision)
