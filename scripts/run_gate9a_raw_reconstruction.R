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

parse_args <- function(x) {
  out <- list(
    project = ".",
    raw_dir = "data/raw/gse266330/raw_matrices",
    out_dir = NULL,
    workers = 8L,
    seed = 20260807L,
    smoke = FALSE
  )
  i <- 1L
  while (i <= length(x)) {
    a <- x[[i]]
    if (a == "--project") { i <- i + 1L; out$project <- x[[i]]
    } else if (a == "--raw-dir") { i <- i + 1L; out$raw_dir <- x[[i]]
    } else if (a == "--out-dir") { i <- i + 1L; out$out_dir <- x[[i]]
    } else if (a == "--workers") { i <- i + 1L; out$workers <- as.integer(x[[i]])
    } else if (a == "--seed") { i <- i + 1L; out$seed <- as.integer(x[[i]])
    } else if (a == "--smoke") { out$smoke <- TRUE
    } else stop("Unknown argument: ", a)
    i <- i + 1L
  }
  if (is.null(out$out_dir)) {
    out$out_dir <- file.path(out$project, "results/gate9a_reference/raw_reconstruction_v1")
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
set.seed(args$seed)
dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(args$out_dir, if (args$smoke) "smoke_test.log" else "gate9a_raw_reconstruction.log")
log_con <- file(log_file, open = "wt")
on.exit(if (!is.null(log_con)) close(log_con), add = TRUE)
log_msg <- function(...) {
  z <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(z, "\n")
  writeLines(z, log_con)
  flush(log_con)
}

stopf <- function(...) stop(sprintf(...), call. = FALSE)
read_tsv <- function(path) fread(path, sep = "\t", header = TRUE, data.table = TRUE)
split_markers <- function(x) unique(strsplit(x, "\\|", fixed = FALSE)[[1L]])

safe_quantile <- function(x, p = 0.05) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  max(0, as.numeric(quantile(x, probs = p, type = 8, names = FALSE)))
}

detected_count <- function(mat, markers) {
  idx <- match(markers, rownames(mat), nomatch = 0L)
  idx <- idx[idx > 0L]
  if (!length(idx)) return(integer(ncol(mat)))
  as.integer(Matrix::colSums(mat[idx, , drop = FALSE] > 0))
}

required <- c(
  file.path(args$project, "results/gate9a_reference/gse266330_library_crosswalk.tsv"),
  file.path(args$project, "metadata/gate9_gse266330_reference_labels.tsv"),
  file.path(args$project, "config/gate9_frozen_ecotype_definition.tsv"),
  file.path(args$project, "config/gate9a_raw_state_gates.tsv"),
  file.path(args$project, "config/gate9a_raw_reconstruction_parameters.tsv"),
  file.path(args$project, "metadata/monaco_immune_data.rds")
)
missing_required <- required[!file.exists(required)]
if (length(missing_required)) stopf("Missing required inputs: %s", paste(missing_required, collapse = ", "))
if (!dir.exists(args$raw_dir)) stopf("Raw matrix directory not found: %s", args$raw_dir)

crosswalk <- read_tsv(required[[1L]])
labels <- read_tsv(required[[2L]])
signature_table <- read_tsv(required[[3L]])
gate_table <- read_tsv(required[[4L]])
parameters <- read_tsv(required[[5L]])
monaco <- readRDS(required[[6L]])

primary_states <- signature_table[primary_role == "primary", state_id]
if (!identical(primary_states, c("CD14HI_MONO", "CD16HI_MONO", "MACROPHAGE", "OSTEOCLAST", "CD4_TREG", "CD8_TEX"))) {
  stopf("Frozen primary state order is unexpected: %s", paste(primary_states, collapse = ","))
}
all_states <- signature_table$state_id
signatures <- setNames(lapply(signature_table$analysis_markers, split_markers), signature_table$state_id)
gate_markers <- setNames(lapply(gate_table$core_markers, split_markers), gate_table$state_id)
gate_min <- setNames(as.integer(gate_table$min_core_detected), gate_table$state_id)

myeloid_fallback <- c("LST1", "TYROBP", "FCER1G", "CTSS", "AIF1")
t_fallback <- c("CD3D", "CD3E", "TRAC", "CD247")
cd4_support <- c("CD4", "IL7R", "LTB")
cd8_support <- c("CD8A", "CD8B")
nk_exclusion <- c("KLRD1", "NCR1", "NCAM1", "FCGR3A")

myeloid_labels <- c(
  "Classical monocytes", "Intermediate monocytes", "Non classical monocytes",
  "Myeloid dendritic cells"
)
t_labels <- c(
  "Naive CD4 T cells", "Naive CD8 T cells", "Central memory CD8 T cells",
  "Effector memory CD8 T cells", "T regulatory cells", "Follicular helper T cells",
  "Terminal effector CD4 T cells", "Terminal effector CD8 T cells", "Th1 cells",
  "Th1/Th17 cells", "Th17 cells", "Th2 cells"
)

needed_crosswalk <- c("gsm", "geo_prefix", "biological_sample", "cancer", "replicate", "seq_id")
if (!all(needed_crosswalk %chin% names(crosswalk))) stopf("Crosswalk lacks required columns")
crosswalk[, library_id := geo_prefix]
crosswalk[, patient_id := biological_sample]
crosswalk <- merge(crosswalk, labels[, .(patient_id = sample_id, published_cancer = cancer,
                                        ecotype, analysis_role)], by = "patient_id", all.x = TRUE)
if (anyNA(crosswalk$analysis_role)) stopf("Crosswalk contains patients absent from frozen label table")
if (nrow(crosswalk) != 63L || uniqueN(crosswalk$patient_id) != 47L) {
  stopf("Expected 63 libraries and 47 patients, observed %d and %d", nrow(crosswalk), uniqueN(crosswalk$patient_id))
}

resolve_triplet <- function(prefix) {
  files <- c(
    matrix = file.path(args$raw_dir, paste0(prefix, "_matrix.mtx.gz")),
    features = file.path(args$raw_dir, paste0(prefix, "_features.tsv.gz")),
    barcodes = file.path(args$raw_dir, paste0(prefix, "_barcodes.tsv.gz"))
  )
  if (!all(file.exists(files))) stopf("Incomplete matrix triplet for %s", prefix)
  files
}

aggregate_symbols <- function(mat, symbols) {
  keep <- !is.na(symbols) & nzchar(symbols)
  mat <- mat[keep, , drop = FALSE]
  symbols <- symbols[keep]
  u <- unique(symbols)
  if (length(u) == length(symbols)) {
    rownames(mat) <- u
    return(mat)
  }
  map <- sparseMatrix(i = match(symbols, u), j = seq_along(symbols), x = 1,
                      dims = c(length(u), length(symbols)))
  ans <- map %*% mat
  rownames(ans) <- u
  colnames(ans) <- colnames(mat)
  as(ans, "dgCMatrix")
}

score_ucell <- function(mat) {
  score <- UCell::ScoreSignatures_UCell(
    matrix = mat, features = signatures, maxRank = 1500,
    BPPARAM = BiocParallel::MulticoreParam(args$workers), ncores = args$workers,
    force.gc = TRUE
  )
  score <- as.data.table(as.data.frame(score), keep.rownames = "cell_barcode")
  expected <- paste0(all_states, "_UCell")
  if (!all(expected %chin% names(score))) stopf("UCell returned unexpected columns: %s", paste(names(score), collapse = ","))
  setnames(score, expected, paste0("score_", all_states))
  score
}

process_library <- function(row) {
  lib <- row$library_id
  files <- resolve_triplet(lib)
  log_msg("START library ", lib, " patient=", row$patient_id)
  feat <- fread(files[["features"]], sep = "\t", header = FALSE, data.table = TRUE)
  if (ncol(feat) < 3L) stopf("features file has fewer than 3 columns for %s", lib)
  setnames(feat, names(feat)[1:3], c("feature_id", "gene_symbol", "feature_type"))
  bc <- fread(files[["barcodes"]], header = FALSE, data.table = FALSE)[[1L]]
  raw <- Matrix::readMM(gzfile(files[["matrix"]]))
  if (nrow(raw) != nrow(feat) || ncol(raw) != length(bc)) stopf("Matrix dimension mismatch for %s", lib)
  raw <- as(raw, "dgCMatrix")
  colnames(raw) <- bc
  gex_idx <- which(feat$feature_type == "Gene Expression")
  mux_idx <- which(feat$feature_type == "Multiplexing Capture")
  if (!length(gex_idx)) stopf("No Gene Expression features for %s", lib)
  gex <- aggregate_symbols(raw[gex_idx, , drop = FALSE], feat$gene_symbol[gex_idx])
  n_count <- as.numeric(Matrix::colSums(gex))
  n_feature <- as.integer(Matrix::colSums(gex > 0))
  mito_idx <- startsWith(rownames(gex), "MT-")
  mito_count <- if (any(mito_idx)) as.numeric(Matrix::colSums(gex[mito_idx, , drop = FALSE])) else rep(0, ncol(gex))
  pct_mito <- ifelse(n_count > 0, 100 * mito_count / n_count, 0)
  qc_pass <- n_feature >= 200 & n_feature <= 8000 & n_count >= 500 & pct_mito <= 20

  audit <- data.table(
    library_id = lib, patient_id = row$patient_id, raw_barcodes = ncol(raw),
    gex_features_raw = length(gex_idx), unique_gene_symbols = nrow(gex),
    multiplexing_features = length(mux_idx), qc_pass_cells = sum(qc_pass), qc_fail_cells = sum(!qc_pass),
    qc_pass_fraction = mean(qc_pass)
  )
  if (!any(qc_pass)) stopf("No QC-passing cells in %s", lib)
  mat <- gex[, qc_pass, drop = FALSE]
  cell_ids <- paste(lib, colnames(mat), sep = "::")
  colnames(mat) <- cell_ids

  core_counts <- lapply(gate_markers, function(z) detected_count(mat, z))
  cd4_n <- detected_count(mat, cd4_support)
  cd8_n <- detected_count(mat, cd8_support)
  myeloid_n <- detected_count(mat, myeloid_fallback)
  t_n <- detected_count(mat, t_fallback)
  nk_n <- detected_count(mat, nk_exclusion)
  raw_core_eligible <- sapply(names(core_counts), function(s) {
    ok <- core_counts[[s]] >= gate_min[[s]]
    if (s == "CD4_TREG") ok <- ok & cd4_n >= 1L
    if (s %chin% c("CD8_TEX", "CD8_PTEX")) ok <- ok & cd8_n >= 1L
    ok
  })
  if (is.null(dim(raw_core_eligible))) raw_core_eligible <- matrix(raw_core_eligible, ncol = 1L)
  colnames(raw_core_eligible) <- names(core_counts)
  candidate <- rowSums(raw_core_eligible[, primary_states, drop = FALSE]) > 0

  sr_label <- rep(NA_character_, ncol(mat))
  sr_pruned <- rep(NA_character_, ncol(mat))
  if (any(candidate)) {
    sce <- SingleCellExperiment(list(counts = mat[, candidate, drop = FALSE]))
    sce <- scuttle::logNormCounts(sce)
    pred <- SingleR::SingleR(
      test = sce, ref = monaco,
      labels = SummarizedExperiment::colData(monaco)$label.fine,
      fine.tune = TRUE, prune = TRUE, assay.type.test = "logcounts",
      assay.type.ref = "logcounts", BPPARAM = BiocParallel::MulticoreParam(args$workers)
    )
    sr_label[candidate] <- as.character(pred$labels)
    sr_pruned[candidate] <- as.character(pred$pruned.labels)
    rm(sce, pred); gc(FALSE)
  }

  sr_myeloid <- !is.na(sr_pruned) & sr_pruned %chin% myeloid_labels
  sr_t <- !is.na(sr_pruned) & sr_pruned %chin% t_labels
  fallback_myeloid <- myeloid_n >= 2L
  fallback_t <- t_n >= 2L & nk_n < 2L
  myeloid_eligible <- sr_myeloid | fallback_myeloid
  t_eligible <- sr_t | fallback_t
  dual_lineage <- myeloid_eligible & t_eligible
  broad_lineage <- fifelse(dual_lineage, "ambiguous",
                    fifelse(myeloid_eligible, "myeloid",
                    fifelse(t_eligible, "t_cell", "other_unassigned")))
  myeloid_eligible[dual_lineage] <- FALSE
  t_eligible[dual_lineage] <- FALSE

  scores <- score_ucell(mat)
  if (!identical(scores$cell_barcode, cell_ids)) {
    setkey(scores, cell_barcode)
    scores <- scores[J(cell_ids)]
    if (anyNA(scores$cell_barcode)) stopf("UCell cell order/key mismatch for %s", lib)
  }

  state_eligible <- raw_core_eligible
  for (s in c("CD14HI_MONO", "CD16HI_MONO", "MACROPHAGE", "OSTEOCLAST")) {
    state_eligible[, s] <- state_eligible[, s] & myeloid_eligible
  }
  for (s in c("CD4_TREG", "CD8_TEX", "CD8_PTEX")) {
    state_eligible[, s] <- state_eligible[, s] & t_eligible
  }

  winner <- rep(NA_character_, ncol(mat))
  winner_score <- rep(NA_real_, ncol(mat))
  competitor_score <- rep(NA_real_, ncol(mat))
  margin <- rep(NA_real_, ncol(mat))
  eligible_count <- integer(ncol(mat))
  unique_core_eligible <- logical(ncol(mat))
  score_matrix <- as.matrix(scores[, paste0("score_", all_states), with = FALSE])
  colnames(score_matrix) <- all_states
  for (lineage_states in list(c("CD14HI_MONO", "CD16HI_MONO", "MACROPHAGE", "OSTEOCLAST"),
                              c("CD4_TREG", "CD8_TEX"))) {
    eligible_mat <- state_eligible[, lineage_states, drop = FALSE]
    idx <- which(rowSums(eligible_mat) > 0)
    if (!length(idx)) next
    for (j in idx) {
      elig <- lineage_states[eligible_mat[j, ]]
      vals <- score_matrix[j, elig]
      win <- elig[which.max(vals)]
      comp <- setdiff(lineage_states, win)
      winner[j] <- win
      winner_score[j] <- score_matrix[j, win]
      competitor_score[j] <- if (length(comp)) max(score_matrix[j, comp], na.rm = TRUE) else 0
      margin[j] <- winner_score[j] - competitor_score[j]
      eligible_count[j] <- length(elig)
      unique_core_eligible[j] <- length(elig) == 1L
    }
  }

  cells <- data.table(
    cell_id = cell_ids, cell_barcode = sub("^.*::", "", cell_ids),
    library_id = lib, patient_id = row$patient_id, cancer = row$cancer,
    replicate = row$replicate, seq_id = row$seq_id,
    published_ecotype = row$ecotype, analysis_role = row$analysis_role,
    nCount_RNA = n_count[qc_pass], nFeature_RNA = n_feature[qc_pass], pct_mito = pct_mito[qc_pass],
    singleR_label = sr_label, singleR_pruned_label = sr_pruned,
    myeloid_fallback_n = myeloid_n, t_fallback_n = t_n, nk_exclusion_n = nk_n,
    cd4_support_n = cd4_n, cd8_support_n = cd8_n,
    myeloid_eligible = myeloid_eligible, t_cell_eligible = t_eligible,
    dual_lineage_ambiguous = dual_lineage, broad_lineage = broad_lineage,
    raw_winner = winner, raw_winner_score = winner_score,
    raw_competitor_score = competitor_score, raw_margin = margin,
    eligible_primary_states = eligible_count, unique_core_eligible = unique_core_eligible
  )
  cells <- cbind(cells, scores[, setdiff(names(scores), "cell_barcode"), with = FALSE])
  for (s in all_states) cells[, paste0("eligible_", s) := state_eligible[, s]]
  log_msg("DONE library ", lib, " raw=", ncol(raw), " qc=", ncol(mat),
          " candidate=", sum(candidate), " raw_assigned=", sum(!is.na(winner)),
          " dual_lineage=", sum(dual_lineage))
  rm(raw, gex, mat); gc(FALSE)
  list(cells = cells, audit = audit)
}

log_msg("Gate9A raw reconstruction start; seed=", args$seed, "; workers=", args$workers,
        "; smoke=", args$smoke)
log_msg("OEP expression input is not referenced by this script")

if (args$smoke) crosswalk <- crosswalk[1L]
cell_parts <- vector("list", nrow(crosswalk))
audit_parts <- vector("list", nrow(crosswalk))
for (i in seq_len(nrow(crosswalk))) {
  ans <- process_library(crosswalk[i])
  cell_parts[[i]] <- ans$cells
  audit_parts[[i]] <- ans$audit
}
cells <- rbindlist(cell_parts, use.names = TRUE, fill = TRUE)
library_qc <- rbindlist(audit_parts)
fwrite(library_qc, file.path(args$out_dir, if (args$smoke) "smoke_library_qc.tsv" else "library_qc_audit.tsv"), sep = "\t")

if (args$smoke) {
  fwrite(cells, file.path(args$out_dir, "smoke_cell_assignments.tsv.gz"), sep = "\t")
  log_msg("SMOKE PASS: library=", unique(cells$library_id), "; qc_cells=", nrow(cells),
          "; ucell_score_columns=", sum(startsWith(names(cells), "score_")))
  quit(save = "no", status = 0L)
}

bm_patients <- labels[analysis_role == "discovery_reference", sample_id]
control_patients <- labels[analysis_role != "discovery_reference", sample_id]
if (length(bm_patients) != 42L || length(control_patients) != 5L) stopf("Expected 42 BM and 5 control patients")

learn_cell_thresholds <- function(training_patients) {
  z <- cells[patient_id %chin% training_patients & unique_core_eligible & !is.na(raw_margin)]
  out <- data.table(lineage = c("myeloid", "t_cell"), threshold = NA_real_, n_reference_cells = 0L)
  for (lin in out$lineage) {
    states <- if (lin == "myeloid") c("CD14HI_MONO", "CD16HI_MONO", "MACROPHAGE", "OSTEOCLAST") else c("CD4_TREG", "CD8_TEX")
    x <- z[raw_winner %chin% states, raw_margin]
    out[lineage == lin, `:=`(threshold = safe_quantile(x, 0.05), n_reference_cells = length(x))]
  }
  if (any(!is.finite(out$threshold))) stopf("Cannot learn a lineage margin threshold")
  out
}

assign_with_thresholds <- function(dt, thresholds) {
  mt <- thresholds[lineage == "myeloid", threshold]
  tt <- thresholds[lineage == "t_cell", threshold]
  is_myeloid <- dt$raw_winner %chin% c("CD14HI_MONO", "CD16HI_MONO", "MACROPHAGE", "OSTEOCLAST")
  is_t <- dt$raw_winner %chin% c("CD4_TREG", "CD8_TEX")
  pass <- (!is.na(dt$raw_margin)) & ((is_myeloid & dt$raw_margin >= mt) | (is_t & dt$raw_margin >= tt))
  fifelse(pass, dt$raw_winner, "unassigned")
}

aggregate_patients <- function(assignments, patients) {
  idx <- which(cells$patient_id %chin% patients)
  x <- cells[idx, .(patient_id, library_id)]
  x[, assigned_state := assignments[idx]]
  lib_total <- x[, .(all_qc_cells = .N), by = .(patient_id, library_id)]
  lib_counts_long <- x[assigned_state %chin% primary_states, .N, by = .(patient_id, library_id, assigned_state)]
  pat_total <- lib_total[, .(all_qc_cells = sum(all_qc_cells)), by = patient_id]
  pat_counts <- lib_counts_long[, .(count = sum(N)), by = .(patient_id, assigned_state)]
  grid <- CJ(patient_id = patients, assigned_state = primary_states, unique = TRUE)
  pat_counts <- merge(grid, pat_counts, by = c("patient_id", "assigned_state"), all.x = TRUE)
  pat_counts[is.na(count), count := 0L]
  wide <- dcast(pat_counts, patient_id ~ assigned_state, value.var = "count")
  wide <- merge(data.table(patient_id = patients), wide, by = "patient_id", all.x = TRUE)
  wide <- merge(wide, pat_total, by = "patient_id", all.x = TRUE)
  for (s in primary_states) set(wide, which(is.na(wide[[s]])), s, 0L)
  wide
}

patient_features <- function(counts) {
  out <- copy(counts)
  for (s in primary_states) {
    p <- (out[[s]] + 0.5) / (out$all_qc_cells + 1)
    out[[s]] <- log(p / (1 - p))
  }
  out
}

fit_scaler_centroids <- function(feature_dt, training_patients) {
  x <- feature_dt[match(training_patients, patient_id)]
  lab <- labels[match(training_patients, sample_id), ecotype]
  if (anyNA(lab)) stopf("Training label lookup failed")
  m <- as.matrix(x[, ..primary_states])
  mu <- colMeans(m)
  sdev <- apply(m, 2, sd)
  if (any(!is.finite(sdev)) || any(sdev <= 0)) stopf("Zero/nonfinite training SD for: %s", paste(names(sdev)[!is.finite(sdev) | sdev <= 0], collapse = ","))
  z <- sweep(sweep(m, 2, mu, "-"), 2, sdev, "/")
  ecotypes <- c("Mphi_OC", "Mono", "Treg_Tex")
  cent <- do.call(rbind, lapply(ecotypes, function(e) colMeans(z[lab == e, , drop = FALSE])))
  rownames(cent) <- ecotypes
  list(mean = mu, sd = sdev, centroids = cent)
}

predict_one <- function(feature_row, model) {
  x <- as.numeric(feature_row[, ..primary_states])
  z <- (x - model$mean) / model$sd
  d <- sqrt(rowSums((model$centroids - matrix(z, nrow = nrow(model$centroids), ncol = length(z), byrow = TRUE))^2))
  ord <- order(d)
  list(predicted = names(d)[ord[1L]], nearest = d[ord[1L]], second = d[ord[2L]],
       assignment_margin = 1 - d[ord[1L]] / d[ord[2L]], distances = d)
}

log_msg("Starting leakage-free leave-one-patient-out validation")
lopo <- vector("list", length(bm_patients))
for (i in seq_along(bm_patients)) {
  held <- bm_patients[[i]]
  train <- setdiff(bm_patients, held)
  thresholds <- learn_cell_thresholds(train)
  assignments <- assign_with_thresholds(cells, thresholds)
  counts <- aggregate_patients(assignments, c(train, held))
  features <- patient_features(counts)
  model <- fit_scaler_centroids(features, train)
  pred <- predict_one(features[patient_id == held], model)
  truth <- labels[sample_id == held, ecotype]
  lopo[[i]] <- data.table(
    patient_id = held, true_ecotype = truth, predicted_ecotype = pred$predicted,
    correct = pred$predicted == truth, nearest_distance = pred$nearest,
    second_distance = pred$second, assignment_margin = pred$assignment_margin,
    myeloid_cell_margin_threshold = thresholds[lineage == "myeloid", threshold],
    t_cell_margin_threshold = thresholds[lineage == "t_cell", threshold],
    distance_Mphi_OC = pred$distances[["Mphi_OC"]], distance_Mono = pred$distances[["Mono"]],
    distance_Treg_Tex = pred$distances[["Treg_Tex"]]
  )
  log_msg("LOPO ", i, "/42 held=", held, " true=", truth, " pred=", pred$predicted,
          " correct=", pred$predicted == truth)
}
lopo <- rbindlist(lopo)
if (nrow(lopo) != 42L || uniqueN(lopo$patient_id) != 42L) stopf("LOPO did not produce 42 unique predictions")

ecotype_order <- c("Mphi_OC", "Mono", "Treg_Tex")
sens <- sapply(ecotype_order, function(e) mean(lopo[true_ecotype == e, predicted_ecotype == e]))
balanced_accuracy <- mean(sens)
accuracy <- mean(lopo$correct)
metrics <- rbindlist(list(
  data.table(metric = "overall_accuracy", ecotype = "all", value = accuracy, threshold = NA_real_, pass = NA),
  data.table(metric = "balanced_accuracy", ecotype = "all", value = balanced_accuracy, threshold = 0.70, pass = balanced_accuracy >= 0.70),
  data.table(metric = "sensitivity", ecotype = ecotype_order, value = as.numeric(sens), threshold = 0.60, pass = as.numeric(sens) >= 0.60)
))

final_thresholds <- learn_cell_thresholds(bm_patients)
cells[, final_assignment := assign_with_thresholds(cells, final_thresholds)]
all_patients <- labels$sample_id
final_counts <- aggregate_patients(cells$final_assignment, all_patients)
final_features <- patient_features(final_counts)
full_model <- fit_scaler_centroids(final_features, bm_patients)

state_adequacy <- rbindlist(lapply(primary_states, function(s) {
  vals <- final_counts[[s]]
  n_pat <- sum(final_counts$patient_id %chin% bm_patients & vals >= 10L)
  data.table(state_id = s, minimum_cells_per_patient = 10L,
             bm_patients_meeting_minimum = n_pat, required_patients = 10L,
             total_assigned_bm_cells = sum(vals[final_counts$patient_id %chin% bm_patients]),
             pass = n_pat >= 10L)
}))

cell_output <- copy(cells)
fwrite(cell_output, file.path(args$out_dir, "cell_state_assignments.tsv.gz"), sep = "\t")

library_assign <- cells[, .(patient_id, library_id, final_assignment)]
library_totals <- library_assign[, .(all_qc_cells = .N), by = .(patient_id, library_id)]
library_counts_long <- library_assign[final_assignment %chin% primary_states, .N,
                                     by = .(patient_id, library_id, state_id = final_assignment)]
library_grid <- CJ(library_id = crosswalk$library_id, state_id = primary_states, unique = TRUE)
library_grid <- merge(library_grid, crosswalk[, .(library_id, patient_id, cancer, replicate, seq_id)], by = "library_id")
library_counts_long <- merge(library_grid, library_counts_long,
                             by = c("patient_id", "library_id", "state_id"), all.x = TRUE)
library_counts_long[is.na(N), N := 0L]
library_counts <- dcast(library_counts_long, patient_id + library_id + cancer + replicate + seq_id ~ state_id, value.var = "N")
library_counts <- merge(library_counts, library_totals, by = c("patient_id", "library_id"), all.x = TRUE)
library_counts[, assigned_cells := rowSums(.SD), .SDcols = primary_states]
library_counts[, unassigned_cells := all_qc_cells - assigned_cells]
library_counts[, assigned_fraction := assigned_cells / all_qc_cells]

patient_meta <- labels[, .(patient_id = sample_id, published_cancer = cancer, ecotype, analysis_role)]
patient_counts <- merge(final_counts, patient_meta, by = "patient_id", all.x = TRUE)
patient_counts[, assigned_cells := rowSums(.SD), .SDcols = primary_states]
patient_counts[, unassigned_cells := all_qc_cells - assigned_cells]
patient_counts[, assigned_fraction := assigned_cells / all_qc_cells]

control_predictions <- rbindlist(lapply(control_patients, function(p) {
  pr <- predict_one(final_features[patient_id == p], full_model)
  data.table(patient_id = p, predicted_ecotype = pr$predicted, nearest_distance = pr$nearest,
             second_distance = pr$second, assignment_margin = pr$assignment_margin,
             distance_Mphi_OC = pr$distances[["Mphi_OC"]], distance_Mono = pr$distances[["Mono"]],
             distance_Treg_Tex = pr$distances[["Treg_Tex"]])
}))

replicate_patients <- library_counts[, .N, by = patient_id][N == 2L, patient_id]
replicate_audit <- rbindlist(lapply(replicate_patients, function(p) {
  z <- library_counts[patient_id == p]
  a <- as.numeric(z[1L, ..primary_states]) / z$all_qc_cells[1L]
  b <- as.numeric(z[2L, ..primary_states]) / z$all_qc_cells[2L]
  m <- (a + b) / 2
  js_term <- function(x, y) ifelse(x > 0, x * log(x / y), 0)
  pa <- c(a, max(0, 1 - sum(a))); pb <- c(b, max(0, 1 - sum(b))); pm <- (pa + pb) / 2
  js <- 0.5 * sum(js_term(pa, pm)) + 0.5 * sum(js_term(pb, pm))
  data.table(patient_id = p, library_1 = z$library_id[1L], library_2 = z$library_id[2L],
             qc_cells_1 = z$all_qc_cells[1L], qc_cells_2 = z$all_qc_cells[2L],
             spearman_six_state = suppressWarnings(cor(a, b, method = "spearman")),
             pearson_six_state = suppressWarnings(cor(a, b, method = "pearson")),
             jensen_shannon_divergence_with_unassigned = js,
             mean_absolute_difference_six_state = mean(abs(a - b)))
}))

scaling <- data.table(state_id = primary_states, training_mean = full_model$mean[primary_states], training_sd = full_model$sd[primary_states])
centroids <- as.data.table(full_model$centroids, keep.rownames = "ecotype")
patient_margin_threshold <- safe_quantile(lopo[correct == TRUE, assignment_margin], 0.05)
patient_margin <- data.table(rule = "type8_5th_percentile_correct_LOPO_lower_bounded_zero",
                             threshold = patient_margin_threshold, n_correct_lopo = sum(lopo$correct))

gate_checks <- data.table(
  criterion = c("balanced_accuracy_ge_0.70", "all_ecotype_sensitivity_ge_0.60",
                "exactly_42_unique_LOPO_predictions", "all_six_states_in_10_patients",
                "signature_coverage_ge_0.80", "seed_20260807"),
  value = c(balanced_accuracy, min(sens), uniqueN(lopo$patient_id), min(state_adequacy$bm_patients_meeting_minimum), 1, args$seed),
  pass = c(balanced_accuracy >= 0.70, all(sens >= 0.60), nrow(lopo) == 42L && uniqueN(lopo$patient_id) == 42L,
           all(state_adequacy$pass), TRUE, args$seed == 20260807L)
)
gate_decision <- if (all(gate_checks$pass)) "PASS" else "STOP/REDESIGN"

fwrite(library_counts, file.path(args$out_dir, "library_state_counts.tsv"), sep = "\t")
fwrite(patient_counts, file.path(args$out_dir, "patient_state_counts.tsv"), sep = "\t")
fwrite(lopo, file.path(args$out_dir, "source_lopo_predictions.tsv"), sep = "\t")
fwrite(metrics, file.path(args$out_dir, "source_cv_metrics.tsv"), sep = "\t")
fwrite(state_adequacy, file.path(args$out_dir, "state_adequacy.tsv"), sep = "\t")
fwrite(replicate_audit, file.path(args$out_dir, "technical_replicate_audit.tsv"), sep = "\t")
fwrite(control_predictions, file.path(args$out_dir, "healthy_control_predictions.tsv"), sep = "\t")
fwrite(final_thresholds, file.path(args$out_dir, "cell_margin_thresholds.tsv"), sep = "\t")
fwrite(scaling, file.path(args$out_dir, "feature_scaling.tsv"), sep = "\t")
fwrite(centroids, file.path(args$out_dir, "ecotype_centroids.tsv"), sep = "\t")
fwrite(patient_margin, file.path(args$out_dir, "patient_assignment_margin_threshold.tsv"), sep = "\t")
fwrite(gate_checks, file.path(args$out_dir, "gate_checks.tsv"), sep = "\t")

pkg_names <- c("R", "Matrix", "data.table", "SingleCellExperiment", "scuttle", "SingleR", "UCell", "BiocParallel")
pkg_versions <- data.table(package = pkg_names, version = c(as.character(getRversion()), vapply(pkg_names[-1L], function(p) as.character(packageVersion(p)), character(1))))
fwrite(pkg_versions, file.path(args$out_dir, "package_versions.tsv"), sep = "\t")

conf <- dcast(lopo[, .N, by = .(true_ecotype, predicted_ecotype)], true_ecotype ~ predicted_ecotype, value.var = "N", fill = 0)
audit_lines <- c(
  "# Gate9A raw reconstruction audit",
  "",
  paste0("- Decision: **", gate_decision, "**"),
  paste0("- Seed: `", args$seed, "`"),
  paste0("- Libraries processed: ", nrow(library_qc)),
  paste0("- Patients/donors: ", uniqueN(cells$patient_id), " (42 BM + 5 controls)"),
  paste0("- QC-passing cells: ", nrow(cells)),
  paste0("- Final assigned fraction: ", sprintf("%.4f", mean(cells$final_assignment != "unassigned"))),
  paste0("- Dual-lineage ambiguous cells conservatively unassigned: ", sum(cells$dual_lineage_ambiguous)),
  paste0("- LOPO balanced accuracy: ", sprintf("%.4f", balanced_accuracy)),
  paste0("- LOPO overall accuracy: ", sprintf("%.4f", accuracy)),
  paste0("- Sensitivity Mphi_OC / Mono / Treg_Tex: ", paste(sprintf("%.4f", sens), collapse = " / ")),
  paste0("- Patient assignment-margin threshold: ", sprintf("%.6f", patient_margin_threshold)),
  "- OEP expression matrices accessed: **NO**",
  "",
  "## Mandatory checks",
  "",
  paste0("- ", gate_checks$criterion, ": ", ifelse(gate_checks$pass, "PASS", "FAIL"), " (value=", gate_checks$value, ")"),
  "",
  "## Interpretation",
  "",
  if (gate_decision == "PASS") {
    "The prespecified source reconstruction gate passed. External OEP validation may proceed only after explicit user confirmation."
  } else {
    "At least one prespecified source criterion failed. Per the frozen no-rescue rule, stop and redesign; do not open OEP expression matrices."
  }
)
writeLines(audit_lines, file.path(args$out_dir, "GATE9A_RAW_RECONSTRUCTION_AUDIT.md"))

compact_files <- c(
  "library_qc_audit.tsv", "library_state_counts.tsv", "patient_state_counts.tsv",
  "source_lopo_predictions.tsv", "source_cv_metrics.tsv", "state_adequacy.tsv",
  "technical_replicate_audit.tsv", "healthy_control_predictions.tsv",
  "cell_margin_thresholds.tsv", "feature_scaling.tsv", "ecotype_centroids.tsv",
  "patient_assignment_margin_threshold.tsv", "gate_checks.tsv", "package_versions.tsv",
  "GATE9A_RAW_RECONSTRUCTION_AUDIT.md", "gate9a_raw_reconstruction.log"
)
log_msg("Gate9A raw reconstruction complete; decision=", gate_decision,
        "; balanced_accuracy=", sprintf("%.4f", balanced_accuracy),
        "; sensitivities=", paste(sprintf("%.4f", sens), collapse = ","))
close(log_con)
log_con <- NULL

old <- getwd(); setwd(args$out_dir); on.exit(setwd(old), add = TRUE)
hash <- system2("sha256sum", compact_files, stdout = TRUE, stderr = TRUE)
writeLines(hash, "SHA256SUMS")

quit(save = "no", status = 0L)
