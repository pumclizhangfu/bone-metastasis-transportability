#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(SingleCellExperiment)
  library(scuttle)
  library(SingleR)
  library(UCell)
  library(BiocParallel)
  library(cluster)
})

options(stringsAsFactors = FALSE, warn = 1)

parse_args <- function(x) {
  out <- list(
    project = ".",
    raw_dir = "data/raw/oep005136/mbone_extracted",
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
  if (is.null(out$out_dir)) out$out_dir <- file.path(out$project, "results/gate9b_validation/external_validation_v1")
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
set.seed(args$seed)
dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(args$out_dir, if (args$smoke) "smoke_test.log" else "gate9b_external_validation.log")
log_con <- file(log_file, open = "wt")
on.exit(if (!is.null(log_con)) close(log_con), add = TRUE)
log_msg <- function(...) {
  z <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(z, "\n")
  writeLines(z, log_con); flush(log_con)
}
stopf <- function(...) stop(sprintf(...), call. = FALSE)
read_tsv <- function(path) fread(path, sep = "\t", header = TRUE, data.table = TRUE)
split_markers <- function(x) unique(strsplit(x, "\\|", fixed = FALSE)[[1L]])

required <- c(
  file.path(args$project, "results/gate9b_validation/input_audit/oep005136_mbone_crosswalk.tsv"),
  file.path(args$project, "config/gate9_frozen_ecotype_definition.tsv"),
  file.path(args$project, "config/gate9a_raw_state_gates.tsv"),
  file.path(args$project, "config/gate9b_validation_parameters.tsv"),
  file.path(args$project, "metadata/monaco_immune_data.rds"),
  file.path(args$project, "results/gate9a_reference/raw_reconstruction_v1/cell_margin_thresholds.tsv"),
  file.path(args$project, "results/gate9a_reference/raw_reconstruction_v1/feature_scaling.tsv"),
  file.path(args$project, "results/gate9a_reference/raw_reconstruction_v1/ecotype_centroids.tsv"),
  file.path(args$project, "results/gate9a_reference/raw_reconstruction_v1/patient_assignment_margin_threshold.tsv")
)
if (any(!file.exists(required))) stopf("Missing required inputs: %s", paste(required[!file.exists(required)], collapse = ", "))
if (!dir.exists(args$raw_dir)) stopf("Raw directory absent: %s", args$raw_dir)

crosswalk <- read_tsv(required[[1L]])
signature_table <- read_tsv(required[[2L]])
gate_table <- read_tsv(required[[3L]])
parameters <- read_tsv(required[[4L]])
monaco <- readRDS(required[[5L]])
source_cell_thresholds <- read_tsv(required[[6L]])
source_scaling <- read_tsv(required[[7L]])
source_centroids_dt <- read_tsv(required[[8L]])
source_patient_margin <- read_tsv(required[[9L]])$threshold[[1L]]

crosswalk[, primary_analysis := as.logical(primary_analysis)]
if (nrow(crosswalk) != 53L || sum(crosswalk$primary_analysis) != 49L || uniqueN(crosswalk[primary_analysis == TRUE, cancer_code]) != 11L) {
  stopf("Expected 53 archive samples, 49 primary patients, and 11 primary origins")
}

primary_states <- signature_table[primary_role == "primary", state_id]
expected_states <- c("CD14HI_MONO", "CD16HI_MONO", "MACROPHAGE", "OSTEOCLAST", "CD4_TREG", "CD8_TEX")
if (!identical(primary_states, expected_states)) stopf("Unexpected primary state order")
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
t_labels <- c(
  "Naive CD4 T cells", "Naive CD8 T cells", "Central memory CD8 T cells",
  "Effector memory CD8 T cells", "T regulatory cells", "Follicular helper T cells",
  "Terminal effector CD4 T cells", "Terminal effector CD8 T cells", "Th1 cells",
  "Th1/Th17 cells", "Th17 cells", "Th2 cells"
)

source_mean <- setNames(source_scaling$training_mean, source_scaling$state_id)[primary_states]
source_sd <- setNames(source_scaling$training_sd, source_scaling$state_id)[primary_states]
source_centroids <- as.matrix(source_centroids_dt[, ..primary_states])
rownames(source_centroids) <- source_centroids_dt$ecotype
ecotypes <- c("Mphi_OC", "Mono", "Treg_Tex")
source_centroids <- source_centroids[ecotypes, , drop = FALSE]
myeloid_margin <- source_cell_thresholds[lineage == "myeloid", threshold]
t_margin <- source_cell_thresholds[lineage == "t_cell", threshold]
if (length(myeloid_margin) != 1L || length(t_margin) != 1L || !is.finite(source_patient_margin)) stopf("Frozen margin inputs invalid")

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
  if (length(u) == length(symbols)) { rownames(mat) <- u; return(as(mat, "dgCMatrix")) }
  map <- sparseMatrix(i = match(symbols, u), j = seq_along(symbols), x = 1,
                      dims = c(length(u), length(symbols)))
  ans <- map %*% mat
  rownames(ans) <- u; colnames(ans) <- colnames(mat)
  as(ans, "dgCMatrix")
}

score_ucell <- function(mat) {
  z <- UCell::ScoreSignatures_UCell(
    matrix = mat, features = signatures, maxRank = 1500,
    BPPARAM = BiocParallel::MulticoreParam(args$workers), ncores = args$workers,
    force.gc = TRUE
  )
  z <- as.data.table(as.data.frame(z), keep.rownames = "cell_id")
  expected <- paste0(all_states, "_UCell")
  if (!all(expected %chin% names(z))) stopf("Unexpected UCell columns")
  setnames(z, expected, paste0("score_", all_states))
  z
}

process_library <- function(row) {
  d <- row$archive_dir
  path <- file.path(args$raw_dir, d)
  files <- c(
    matrix = file.path(path, "matrix.mtx.gz"),
    features = file.path(path, "features.tsv.gz"),
    barcodes = file.path(path, "barcodes.tsv.gz")
  )
  if (!all(file.exists(files))) stopf("Incomplete triplet for %s", d)
  log_msg("START library ", d, " patient=", row$patient_id, " primary=", row$primary_analysis)
  feat <- fread(files[["features"]], sep = "\t", header = FALSE, data.table = TRUE)
  if (ncol(feat) < 2L) stopf("Feature table invalid for %s", d)
  if (ncol(feat) == 2L) feat[, V3 := "Gene Expression"]
  setnames(feat, names(feat)[1:3], c("feature_id", "gene_symbol", "feature_type"))
  bc <- fread(files[["barcodes"]], header = FALSE, data.table = FALSE)[[1L]]
  raw <- as(Matrix::readMM(gzfile(files[["matrix"]])), "dgCMatrix")
  if (nrow(raw) != nrow(feat) || ncol(raw) != length(bc)) stopf("Dimension mismatch for %s", d)
  colnames(raw) <- bc
  gex_idx <- which(feat$feature_type == "Gene Expression")
  other_idx <- which(feat$feature_type != "Gene Expression")
  gex <- aggregate_symbols(raw[gex_idx, , drop = FALSE], feat$gene_symbol[gex_idx])
  n_count <- as.numeric(Matrix::colSums(gex))
  n_feature <- as.integer(Matrix::colSums(gex > 0))
  mito_idx <- startsWith(rownames(gex), "MT-")
  mito_count <- if (any(mito_idx)) as.numeric(Matrix::colSums(gex[mito_idx, , drop = FALSE])) else rep(0, ncol(gex))
  pct_mito <- ifelse(n_count > 0, 100 * mito_count / n_count, 0)
  qc_pass <- n_feature >= 200 & n_feature <= 8000 & n_count >= 500 & pct_mito <= 20
  if (!any(qc_pass)) stopf("No QC-passing cells for %s", d)

  mat <- gex[, qc_pass, drop = FALSE]
  cell_ids <- paste(d, colnames(mat), sep = "::")
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

  sr_label <- rep(NA_character_, ncol(mat)); sr_pruned <- rep(NA_character_, ncol(mat))
  if (any(candidate)) {
    sce <- SingleCellExperiment(list(counts = mat[, candidate, drop = FALSE]))
    sce <- scuttle::logNormCounts(sce)
    pred <- SingleR::SingleR(
      test = sce, ref = monaco,
      labels = SummarizedExperiment::colData(monaco)$label.fine,
      fine.tune = TRUE, prune = TRUE, assay.type.test = "logcounts", assay.type.ref = "logcounts",
      BPPARAM = BiocParallel::MulticoreParam(args$workers)
    )
    sr_label[candidate] <- as.character(pred$labels)
    sr_pruned[candidate] <- as.character(pred$pruned.labels)
    rm(sce, pred); gc(FALSE)
  }

  myeloid_eligible <- (!is.na(sr_pruned) & sr_pruned %chin% myeloid_labels) | myeloid_n >= 2L
  t_eligible <- (!is.na(sr_pruned) & sr_pruned %chin% t_labels) | (t_n >= 2L & nk_n < 2L)
  dual_lineage <- myeloid_eligible & t_eligible
  broad_lineage <- fifelse(dual_lineage, "ambiguous",
                    fifelse(myeloid_eligible, "myeloid", fifelse(t_eligible, "t_cell", "other_unassigned")))
  myeloid_eligible[dual_lineage] <- FALSE; t_eligible[dual_lineage] <- FALSE

  scores <- score_ucell(mat)
  if (!identical(scores$cell_id, cell_ids)) {
    setkey(scores, cell_id); scores <- scores[J(cell_ids)]
    if (anyNA(scores$cell_id)) stopf("UCell cell-key mismatch for %s", d)
  }
  state_eligible <- raw_core_eligible
  for (s in c("CD14HI_MONO", "CD16HI_MONO", "MACROPHAGE", "OSTEOCLAST")) state_eligible[, s] <- state_eligible[, s] & myeloid_eligible
  for (s in c("CD4_TREG", "CD8_TEX", "CD8_PTEX")) state_eligible[, s] <- state_eligible[, s] & t_eligible
  score_matrix <- as.matrix(scores[, paste0("score_", all_states), with = FALSE]); colnames(score_matrix) <- all_states
  winner <- rep(NA_character_, ncol(mat)); winner_score <- competitor_score <- margin <- rep(NA_real_, ncol(mat))
  eligible_count <- integer(ncol(mat)); unique_core_eligible <- logical(ncol(mat))
  for (lineage_states in list(c("CD14HI_MONO", "CD16HI_MONO", "MACROPHAGE", "OSTEOCLAST"), c("CD4_TREG", "CD8_TEX"))) {
    em <- state_eligible[, lineage_states, drop = FALSE]
    idx <- which(rowSums(em) > 0)
    for (j in idx) {
      elig <- lineage_states[em[j, ]]; win <- elig[which.max(score_matrix[j, elig])]
      comp <- setdiff(lineage_states, win)
      winner[j] <- win; winner_score[j] <- score_matrix[j, win]
      competitor_score[j] <- max(score_matrix[j, comp], na.rm = TRUE)
      margin[j] <- winner_score[j] - competitor_score[j]
      eligible_count[j] <- length(elig); unique_core_eligible[j] <- length(elig) == 1L
    }
  }
  is_myeloid <- winner %chin% c("CD14HI_MONO", "CD16HI_MONO", "MACROPHAGE", "OSTEOCLAST")
  is_t <- winner %chin% c("CD4_TREG", "CD8_TEX")
  assignment_pass <- !is.na(margin) & ((is_myeloid & margin >= myeloid_margin) | (is_t & margin >= t_margin))
  final_assignment <- fifelse(assignment_pass, winner, "unassigned")

  cells <- data.table(
    cell_id = cell_ids, cell_barcode = sub("^.*::", "", cell_ids), archive_dir = d,
    sample_id = row$sample_id, patient_id = row$patient_id, cancer_code = row$cancer_code,
    primary_analysis = row$primary_analysis, exclusion_reason = row$exclusion_reason,
    nCount_RNA = n_count[qc_pass], nFeature_RNA = n_feature[qc_pass], pct_mito = pct_mito[qc_pass],
    singleR_label = sr_label, singleR_pruned_label = sr_pruned,
    myeloid_fallback_n = myeloid_n, t_fallback_n = t_n, nk_exclusion_n = nk_n,
    cd4_support_n = cd4_n, cd8_support_n = cd8_n,
    myeloid_eligible = myeloid_eligible, t_cell_eligible = t_eligible,
    dual_lineage_ambiguous = dual_lineage, broad_lineage = broad_lineage,
    raw_winner = winner, raw_winner_score = winner_score, raw_competitor_score = competitor_score,
    raw_margin = margin, eligible_primary_states = eligible_count,
    unique_core_eligible = unique_core_eligible, final_assignment = final_assignment
  )
  cells <- cbind(cells, scores[, setdiff(names(scores), "cell_id"), with = FALSE])
  for (s in all_states) cells[, paste0("eligible_", s) := state_eligible[, s]]
  audit <- data.table(
    archive_dir = d, patient_id = row$patient_id, primary_analysis = row$primary_analysis,
    raw_barcodes = ncol(raw), gex_features_raw = length(gex_idx), other_features = length(other_idx),
    unique_gene_symbols = nrow(gex), qc_pass_cells = ncol(mat), qc_fail_cells = sum(!qc_pass),
    qc_pass_fraction = mean(qc_pass), candidate_cells = sum(candidate),
    assigned_cells = sum(final_assignment != "unassigned"), dual_lineage_cells = sum(dual_lineage)
  )
  log_msg("DONE library ", d, " raw=", ncol(raw), " qc=", ncol(mat), " candidate=", sum(candidate),
          " assigned=", sum(final_assignment != "unassigned"), " dual_lineage=", sum(dual_lineage))
  rm(raw, gex, mat); gc(FALSE)
  list(cells = cells, audit = audit)
}

log_msg("Gate9B external validation start; seed=", args$seed, "; workers=", args$workers, "; smoke=", args$smoke)
log_msg("Frozen source inputs loaded; OEP refitting prohibited")
if (args$smoke) crosswalk <- crosswalk[1L]
cell_parts <- vector("list", nrow(crosswalk)); audit_parts <- vector("list", nrow(crosswalk))
for (i in seq_len(nrow(crosswalk))) {
  ans <- process_library(crosswalk[i]); cell_parts[[i]] <- ans$cells; audit_parts[[i]] <- ans$audit
}
cells <- rbindlist(cell_parts, use.names = TRUE, fill = TRUE)
library_audit <- rbindlist(audit_parts)
fwrite(library_audit, file.path(args$out_dir, if (args$smoke) "smoke_library_audit.tsv" else "library_qc_audit.tsv"), sep = "\t")
if (args$smoke) {
  fwrite(cells, file.path(args$out_dir, "smoke_cell_assignments.tsv.gz"), sep = "\t")
  log_msg("SMOKE PASS library=", unique(cells$archive_dir), "; qc_cells=", nrow(cells), "; score_columns=", sum(startsWith(names(cells), "score_")))
  quit(save = "no", status = 0L)
}

aggregate_counts <- function(dt) {
  patients <- unique(dt$patient_id)
  totals <- dt[, .(all_qc_cells = .N), by = patient_id]
  counts_long <- dt[final_assignment %chin% primary_states, .N, by = .(patient_id, state_id = final_assignment)]
  grid <- CJ(patient_id = patients, state_id = primary_states, unique = TRUE)
  counts_long <- merge(grid, counts_long, by = c("patient_id", "state_id"), all.x = TRUE)
  counts_long[is.na(N), N := 0L]
  wide <- dcast(counts_long, patient_id ~ state_id, value.var = "N")
  wide <- merge(wide, totals, by = "patient_id", all.x = TRUE)
  wide
}

transform_counts <- function(count_dt) {
  z <- copy(count_dt)
  for (s in primary_states) {
    p <- (z[[s]] + 0.5) / (z$all_qc_cells + 1)
    logit <- log(p / (1 - p))
    z[[s]] <- (logit - source_mean[[s]]) / source_sd[[s]]
  }
  z
}

predict_frozen <- function(z_dt) {
  x <- as.matrix(z_dt[, ..primary_states])
  dmat <- sapply(ecotypes, function(e) sqrt(rowSums((x - matrix(source_centroids[e, ], nrow(x), length(primary_states), byrow = TRUE))^2)))
  colnames(dmat) <- ecotypes
  ord <- t(apply(dmat, 1, order))
  forced <- ecotypes[ord[, 1L]]
  nearest <- dmat[cbind(seq_len(nrow(dmat)), ord[, 1L])]
  second <- dmat[cbind(seq_len(nrow(dmat)), ord[, 2L])]
  margin <- ifelse(second > 0, 1 - nearest / second, 0)
  nonforced <- ifelse(is.finite(margin) & margin >= source_patient_margin, forced, "unassigned")
  data.table(patient_id = z_dt$patient_id, forced_label = forced, transferred_label = nonforced,
             nearest_distance = nearest, second_distance = second, assignment_margin = margin,
             distance_Mphi_OC = dmat[, "Mphi_OC"], distance_Mono = dmat[, "Mono"], distance_Treg_Tex = dmat[, "Treg_Tex"])
}

patient_counts <- aggregate_counts(cells)
patient_counts <- merge(patient_counts,
  crosswalk[, .(patient_id, sample_id, archive_dir, cancer_code, primary_analysis, exclusion_reason, sex, age, collection_date)],
  by = "patient_id", all.x = TRUE)
patient_counts[, assigned_cells := rowSums(.SD), .SDcols = primary_states]
patient_counts[, unassigned_cells := all_qc_cells - assigned_cells]
patient_counts[, assigned_fraction := assigned_cells / all_qc_cells]
patient_z <- transform_counts(patient_counts[, c("patient_id", primary_states, "all_qc_cells"), with = FALSE])
transfer <- predict_frozen(patient_z)
patient_assignments <- merge(patient_counts, patient_z[, c("patient_id", primary_states), with = FALSE], by = "patient_id", suffixes = c("_count", "_z"))
patient_assignments <- merge(patient_assignments, transfer, by = "patient_id")

# Durable checkpoint before consensus/deletion work. These files contain only
# frozen cell assignments and patient summaries and are also final deliverables.
fwrite(cells, file.path(args$out_dir, "oep005136_cell_state_assignments.tsv.gz"), sep = "\t")
fwrite(library_audit, file.path(args$out_dir, "library_qc_audit.tsv"), sep = "\t")
fwrite(patient_counts, file.path(args$out_dir, "patient_state_counts.tsv"), sep = "\t")
log_msg("CHECKPOINT cell reconstruction and patient counts written")

adjusted_rand <- function(a, b) {
  if (length(a) != length(b) || length(a) < 2L) return(NA_real_)
  tab <- table(a, b)
  choose2 <- function(x) x * (x - 1) / 2
  n2 <- choose2(sum(tab)); if (n2 == 0) return(NA_real_)
  index <- sum(choose2(tab)); aa <- sum(choose2(rowSums(tab))); bb <- sum(choose2(colSums(tab)))
  expected <- aa * bb / n2; upper <- 0.5 * (aa + bb)
  if (upper == expected) return(1)
  (index - expected) / (upper - expected)
}

consensus_cluster <- function(x, k, iterations, fraction, seed) {
  set.seed(seed)
  n <- nrow(x); take <- floor(fraction * n)
  co_sample <- matrix(0L, n, n); co_cluster <- matrix(0L, n, n)
  for (b in seq_len(iterations)) {
    idx <- sort(sample.int(n, take, replace = FALSE))
    labs <- cutree(hclust(dist(x[idx, , drop = FALSE], method = "euclidean"), method = "average"), k = k)
    co_sample[idx, idx] <- co_sample[idx, idx] + 1L
    for (g in split(idx, labs)) co_cluster[g, g] <- co_cluster[g, g] + 1L
  }
  cons <- co_cluster / co_sample; diag(cons) <- 1
  dimnames(cons) <- list(rownames(x), rownames(x))
  if (any(!is.finite(cons))) stopf("Consensus matrix contains unsampled pairs")
  final <- cutree(hclust(as.dist(1 - cons), method = "average"), k = k)
  names(final) <- rownames(x)
  upper <- cons[upper.tri(cons)]
  pac <- mean(upper > 0.10 & upper < 0.90)
  sil <- mean(cluster::silhouette(final, dist(x, method = "euclidean"))[, "sil_width"])
  list(consensus = cons, cluster = final, pac = pac, silhouette = sil,
       min_cluster_size = min(table(final)), sizes = paste(sort(as.integer(table(final))), collapse = "|"))
}

permutation_ari <- function(cluster_labels, transferred_labels, B, seed) {
  ok <- transferred_labels != "unassigned"
  a <- cluster_labels[ok]; b <- transferred_labels[ok]
  obs <- adjusted_rand(a, b)
  if (!is.finite(obs)) return(list(ari = NA_real_, p = NA_real_, n = sum(ok)))
  set.seed(seed)
  perm <- replicate(B, adjusted_rand(a, sample(b, replace = FALSE)))
  list(ari = obs, p = (1 + sum(perm >= obs)) / (B + 1), n = sum(ok))
}

axis_table <- function(z_dt, labels_vec) {
  q <- copy(z_dt)
  q[, transferred_label := labels_vec[match(patient_id, names(labels_vec))]]
  q <- q[transferred_label %chin% ecotypes]
  q[, axis_Mphi_OC := (MACROPHAGE + OSTEOCLAST) / 2]
  q[, axis_Mono := (CD14HI_MONO + CD16HI_MONO) / 2]
  q[, axis_Treg_Tex := (CD4_TREG + CD8_TEX) / 2]
  q[, .(axis_Mphi_OC = mean(axis_Mphi_OC), axis_Mono = mean(axis_Mono), axis_Treg_Tex = mean(axis_Treg_Tex), n = .N), by = transferred_label]
}

topology_audit <- function(z_dt, labels_vec) {
  axes <- axis_table(z_dt, labels_vec)
  out <- rbindlist(lapply(ecotypes, function(e) {
    row <- axes[transferred_label == e]
    if (!nrow(row)) return(data.table(ecotype = e, n = 0L, matched_axis = NA_real_, best_other_axis = NA_real_, pass = FALSE))
    vals <- as.numeric(row[, .(axis_Mphi_OC, axis_Mono, axis_Treg_Tex)])
    names(vals) <- ecotypes
    data.table(ecotype = e, n = row$n, matched_axis = vals[[e]], best_other_axis = max(vals[setdiff(ecotypes, e)]), pass = vals[[e]] > max(vals[setdiff(ecotypes, e)]))
  }))
  out
}

loco_audit <- function(z_dt, labels_vec, cancers_vec) {
  origins <- sort(unique(cancers_vec))
  rbindlist(lapply(origins, function(origin) {
    keep_ids <- names(cancers_vec)[cancers_vec != origin]
    labs <- labels_vec[keep_ids]
    retained <- all(ecotypes %chin% unique(labs[labs != "unassigned"]))
    topo <- topology_audit(z_dt[patient_id %chin% keep_ids], labs)
    data.table(omitted_cancer = origin, remaining_patients = length(keep_ids),
               assigned_Mphi_OC = sum(labs == "Mphi_OC"), assigned_Mono = sum(labs == "Mono"),
               assigned_Treg_Tex = sum(labs == "Treg_Tex"), all_three_labels_retained = retained,
               diagonal_topology_pass = all(topo$pass), pass = retained & all(topo$pass))
  }))
}

log_msg("Cell reconstruction complete; starting patient transfer and consensus validation")
primary_ids <- crosswalk[primary_analysis == TRUE, patient_id]
primary_assign <- patient_assignments[match(primary_ids, patient_id)]
if (anyNA(primary_assign$patient_id) || nrow(primary_assign) != 49L) stopf("Primary patient table mismatch")
X <- as.matrix(primary_assign[, paste0(primary_states, "_z"), with = FALSE]); colnames(X) <- primary_states; rownames(X) <- primary_assign$patient_id
transfer_labels <- setNames(primary_assign$transferred_label, primary_assign$patient_id)
cancers <- setNames(primary_assign$cancer_code, primary_assign$patient_id)

consensus_results <- lapply(2:5, function(k) consensus_cluster(X, k, 1000L, 0.80, args$seed + k))
names(consensus_results) <- as.character(2:5)
consensus_metrics <- rbindlist(lapply(2:5, function(k) {
  r <- consensus_results[[as.character(k)]]
  data.table(k = k, iterations = 1000L, patient_fraction = 0.80, pac = r$pac,
             mean_silhouette = r$silhouette, min_cluster_size = r$min_cluster_size, cluster_sizes = r$sizes)
}))
k3 <- consensus_results[["3"]]
k3_labels <- k3$cluster[primary_assign$patient_id]
agreement <- permutation_ari(k3_labels, primary_assign$transferred_label, 10000L, args$seed + 300L)
patient_assignments[, consensus_k3_cluster := NA_integer_]
patient_assignments[primary_analysis == TRUE, consensus_k3_cluster := as.integer(k3_labels[patient_id])]

log_msg("K=3 consensus complete; PAC=", sprintf("%.4f", k3$pac), "; silhouette=", sprintf("%.4f", k3$silhouette),
        "; ARI=", sprintf("%.4f", agreement$ari), "; permutation_p=", sprintf("%.6f", agreement$p))

subsample_counts_to_label <- function(assignments) {
  counts <- setNames(integer(length(primary_states)), primary_states)
  tab <- table(assignments)
  hit <- intersect(names(tab), primary_states); counts[hit] <- as.integer(tab[hit])
  N <- length(assignments)
  logit <- sapply(primary_states, function(s) {
    p <- (counts[[s]] + 0.5) / (N + 1); log(p / (1 - p))
  })
  z <- (logit - source_mean) / source_sd
  d <- sapply(ecotypes, function(e) sqrt(sum((z - source_centroids[e, ])^2)))
  ord <- order(d); forced <- names(d)[ord[1L]]; margin <- if (d[ord[2L]] > 0) 1 - d[ord[1L]] / d[ord[2L]] else 0
  if (is.finite(margin) && margin >= source_patient_margin) forced else "unassigned"
}

log_msg("Starting 100 x within-patient 80% cell subsampling")
set.seed(args$seed + 400L)
stability_rows <- vector("list", length(primary_ids))
for (i in seq_along(primary_ids)) {
  p <- primary_ids[[i]]; a <- cells[patient_id == p, final_assignment]; n <- length(a); take <- floor(0.80 * n)
  full_label <- transfer_labels[[p]]
  reps <- replicate(100L, subsample_counts_to_label(a[sample.int(n, take, replace = FALSE)]))
  stability_rows[[i]] <- data.table(patient_id = p, cancer_code = cancers[[p]], full_transferred_label = full_label,
                                    qc_cells = n, subsample_cells = take, iterations = 100L,
                                    matching_iterations = sum(reps == full_label), stability = mean(reps == full_label),
                                    subsample_Mphi_OC = sum(reps == "Mphi_OC"), subsample_Mono = sum(reps == "Mono"),
                                    subsample_Treg_Tex = sum(reps == "Treg_Tex"), subsample_unassigned = sum(reps == "unassigned"))
}
stability <- rbindlist(stability_rows)
assigned_stability <- stability[full_transferred_label %chin% ecotypes]
overall_stability <- median(assigned_stability$stability)
within_stability <- assigned_stability[, .(median_stability = median(stability), n_patients = .N), by = .(ecotype = full_transferred_label)]
within_stability <- merge(data.table(ecotype = ecotypes), within_stability, by = "ecotype", all.x = TRUE)
within_stability[is.na(n_patients), n_patients := 0L]

full_topology <- topology_audit(patient_z[patient_id %chin% primary_ids], transfer_labels)
loco <- loco_audit(patient_z[patient_id %chin% primary_ids], transfer_labels, cancers)

assoc_dt <- data.table(cancer_code = cancers, transferred_label = transfer_labels)
assoc_dt <- assoc_dt[transferred_label %chin% ecotypes]
assoc_tab <- table(assoc_dt$cancer_code, assoc_dt$transferred_label)
set.seed(args$seed + 500L)
fisher_result <- fisher.test(assoc_tab, simulate.p.value = TRUE, B = 10000L)
chi <- suppressWarnings(chisq.test(assoc_tab, correct = FALSE))
cramers_v <- sqrt(as.numeric(chi$statistic) / (sum(assoc_tab) * min(nrow(assoc_tab) - 1L, ncol(assoc_tab) - 1L)))
cancer_association <- data.table(test = c("Fisher_Freeman_Halton_simulated", "Cramers_V"),
                                 value = c(fisher_result$p.value, cramers_v),
                                 iterations = c(10000L, NA_integer_), n_assigned_patients = sum(assoc_tab))

endpoint1 <- k3$min_cluster_size >= 5L
endpoint2 <- k3$pac <= 0.20 && k3$silhouette >= 0.25
endpoint3 <- is.finite(agreement$ari) && agreement$ari >= 0.40 && is.finite(agreement$p) && agreement$p < 0.05
endpoint4 <- is.finite(overall_stability) && overall_stability >= 0.80 &&
             all(within_stability$n_patients > 0) && all(within_stability$median_stability >= 0.70)
endpoint5 <- all(loco$pass)
base_full_1to5 <- all(c(endpoint1, endpoint2, endpoint3, endpoint4, endpoint5))

evaluate_without_patient <- function(drop_id, index) {
  keep <- setdiff(primary_ids, drop_id)
  x <- X[keep, , drop = FALSE]
  labels_keep <- transfer_labels[keep]; cancers_keep <- cancers[keep]
  c3 <- consensus_cluster(x, 3L, 1000L, 0.80, args$seed + 1000L + index)
  agr <- permutation_ari(c3$cluster[keep], labels_keep, 10000L, args$seed + 2000L + index)
  stab <- stability[patient_id %chin% keep & full_transferred_label %chin% ecotypes]
  overall <- if (nrow(stab)) median(stab$stability) else NA_real_
  within <- stab[, .(median_stability = median(stability), n_patients = .N), by = .(ecotype = full_transferred_label)]
  within <- merge(data.table(ecotype = ecotypes), within, by = "ecotype", all.x = TRUE)
  within[is.na(n_patients), n_patients := 0L]
  loco_i <- loco_audit(patient_z[patient_id %chin% keep], labels_keep, cancers_keep)
  e1 <- c3$min_cluster_size >= 5L
  e2 <- c3$pac <= 0.20 && c3$silhouette >= 0.25
  e3 <- is.finite(agr$ari) && agr$ari >= 0.40 && is.finite(agr$p) && agr$p < 0.05
  e4 <- is.finite(overall) && overall >= 0.80 && all(within$n_patients > 0) && all(within$median_stability >= 0.70)
  e5 <- all(loco_i$pass)
  data.table(omitted_patient = drop_id, cancer_code = cancers[[drop_id]], remaining_patients = length(keep),
             k3_min_cluster_size = c3$min_cluster_size, k3_pac = c3$pac, k3_silhouette = c3$silhouette,
             ari = agr$ari, ari_permutation_p = agr$p, median_stability = overall,
             minimum_within_ecotype_stability = if (all(within$n_patients > 0L)) min(within$median_stability) else NA_real_,
             endpoint1_cluster_size_pass = e1, endpoint2_consensus_pass = e2,
             endpoint3_agreement_pass = e3, endpoint4_stability_pass = e4,
             endpoint5_loco_pass = e5, full_pass_endpoints_1_to_5 = all(c(e1, e2, e3, e4, e5)))
}

log_msg("Starting 49 leave-one-patient-out influence recomputations")
deletion <- vector("list", length(primary_ids))
for (i in seq_along(primary_ids)) {
  deletion[[i]] <- evaluate_without_patient(primary_ids[[i]], i)
  if (i %% 5L == 0L || i == length(primary_ids)) log_msg("Patient deletion ", i, "/", length(primary_ids), " complete")
}
deletion <- rbindlist(deletion)
fwrite(deletion, file.path(args$out_dir, "leave_one_patient_out_robustness.tsv"), sep = "\t")
log_msg("CHECKPOINT patient-deletion robustness written")
endpoint6 <- all(deletion$full_pass_endpoints_1_to_5 == base_full_1to5)

endpoint_table <- data.table(
  endpoint = c("K3_min_cluster_size", "K3_PAC_and_silhouette", "K3_ARI_and_permutation",
               "cell_subsample_stability", "leave_one_cancer_out", "single_patient_deletion_robustness"),
  value = c(k3$min_cluster_size, paste0(sprintf("%.6f", k3$pac), "|", sprintf("%.6f", k3$silhouette)),
            paste0(sprintf("%.6f", agreement$ari), "|", sprintf("%.6f", agreement$p)),
            paste0(sprintf("%.6f", overall_stability), "|", sprintf("%.6f", min(within_stability$median_stability))),
            sum(loco$pass), sum(deletion$full_pass_endpoints_1_to_5 == base_full_1to5)),
  threshold = c(">=5", "PAC<=0.20 & silhouette>=0.25", "ARI>=0.40 & P<0.05",
                "overall>=0.80 & each>=0.70", paste0("all ", nrow(loco)), paste0("all ", nrow(deletion))),
  pass = c(endpoint1, endpoint2, endpoint3, endpoint4, endpoint5, endpoint6)
)

all_three <- all(ecotypes %chin% unique(transfer_labels[transfer_labels != "unassigned"]))
partial_conditions <- all_three && endpoint4 && all(full_topology$pass)
gate_decision <- if (all(endpoint_table$pass)) "FULL PASS" else if (partial_conditions) "PARTIAL/EXPLORATORY" else "NO-GO"

k3_assign <- data.table(patient_id = names(k3_labels), consensus_k3_cluster = as.integer(k3_labels))
consensus_matrix_dt <- as.data.table(as.data.frame(k3$consensus, check.names = FALSE))
consensus_matrix_dt[, patient_id := rownames(k3$consensus)]
setcolorder(consensus_matrix_dt, c("patient_id", setdiff(names(consensus_matrix_dt), "patient_id")))

fwrite(patient_assignments, file.path(args$out_dir, "patient_assignments.tsv"), sep = "\t")
fwrite(consensus_metrics, file.path(args$out_dir, "consensus_metrics.tsv"), sep = "\t")
fwrite(consensus_matrix_dt, file.path(args$out_dir, "k3_consensus_matrix.tsv.gz"), sep = "\t")
fwrite(k3_assign, file.path(args$out_dir, "k3_cluster_assignments.tsv"), sep = "\t")
fwrite(data.table(ari = agreement$ari, permutation_p = agreement$p, n_nonforced_assigned = agreement$n, permutations = 10000L), file.path(args$out_dir, "k3_transfer_agreement.tsv"), sep = "\t")
fwrite(stability, file.path(args$out_dir, "cell_subsample_stability.tsv"), sep = "\t")
fwrite(within_stability, file.path(args$out_dir, "within_ecotype_stability.tsv"), sep = "\t")
fwrite(full_topology, file.path(args$out_dir, "full_data_topology.tsv"), sep = "\t")
fwrite(loco, file.path(args$out_dir, "leave_one_cancer_out.tsv"), sep = "\t")
fwrite(cancer_association, file.path(args$out_dir, "cancer_origin_association.tsv"), sep = "\t")
fwrite(endpoint_table, file.path(args$out_dir, "gate9b_endpoints.tsv"), sep = "\t")

pkg_names <- c("R", "Matrix", "data.table", "SingleCellExperiment", "scuttle", "SingleR", "UCell", "BiocParallel", "cluster")
pkg_versions <- data.table(package = pkg_names, version = c(as.character(getRversion()), vapply(pkg_names[-1L], function(p) as.character(packageVersion(p)), character(1))))
fwrite(pkg_versions, file.path(args$out_dir, "package_versions.tsv"), sep = "\t")

audit_lines <- c(
  "# Gate9B external validation final audit", "",
  paste0("- Decision: **", gate_decision, "**"),
  paste0("- Primary OEP patients: ", length(primary_ids), " from ", uniqueN(cancers), " cancer origins"),
  paste0("- QC-passing cells (all 53 archive samples): ", nrow(cells)),
  paste0("- Non-forced assigned primary patients: ", sum(transfer_labels != "unassigned"), "/", length(primary_ids)),
  paste0("- Transferred Mphi-OC / Mono / Treg-Tex / unassigned: ", paste(table(factor(transfer_labels, levels = c(ecotypes, "unassigned"))), collapse = " / ")),
  paste0("- K=3 cluster sizes: ", k3$sizes),
  paste0("- K=3 PAC: ", sprintf("%.6f", k3$pac)),
  paste0("- K=3 mean silhouette: ", sprintf("%.6f", k3$silhouette)),
  paste0("- K=3 transfer ARI: ", sprintf("%.6f", agreement$ari)),
  paste0("- K=3 transfer permutation P: ", sprintf("%.6f", agreement$p)),
  paste0("- Overall median cell-subsample stability: ", sprintf("%.6f", overall_stability)),
  paste0("- Minimum within-ecotype median stability: ", sprintf("%.6f", min(within_stability$median_stability))),
  paste0("- Leave-one-cancer-out passes: ", sum(loco$pass), "/", nrow(loco)),
  paste0("- Patient-deletion conclusion preserved: ", sum(deletion$full_pass_endpoints_1_to_5 == base_full_1to5), "/", nrow(deletion)),
  paste0("- Cancer-origin Fisher simulated P: ", sprintf("%.6f", fisher_result$p.value)),
  paste0("- Cancer-origin Cramer's V: ", sprintf("%.6f", cramers_v)),
  "- OEP-derived refitting or threshold tuning: **NO**", "", "## Endpoint checks", "",
  paste0("- ", endpoint_table$endpoint, ": ", ifelse(endpoint_table$pass, "PASS", "FAIL"), " (value=", endpoint_table$value, "; threshold=", endpoint_table$threshold, ")"),
  "", "## Interpretation", "",
  if (gate_decision == "FULL PASS") {
    "All frozen external-validation endpoints passed. Gate9C may be considered only after explicit user confirmation."
  } else if (gate_decision == "PARTIAL/EXPLORATORY") {
    "Stable biological gradients and all three transferred labels are present, but at least one frozen structural endpoint failed. The result is exploratory and does not support a claim of three conserved ecotypes."
  } else {
    "The frozen conditions for a stable external ecotype gradient were not met. Per the no-rescue rule, Gate9C must not proceed as confirmatory validation."
  }
)
writeLines(audit_lines, file.path(args$out_dir, "GATE9B_FINAL_AUDIT.md"))

compact_files <- c(
  "oep005136_cell_state_assignments.tsv.gz", "library_qc_audit.tsv", "patient_state_counts.tsv", "patient_assignments.tsv", "consensus_metrics.tsv",
  "k3_consensus_matrix.tsv.gz", "k3_cluster_assignments.tsv", "k3_transfer_agreement.tsv",
  "cell_subsample_stability.tsv", "within_ecotype_stability.tsv", "full_data_topology.tsv",
  "leave_one_cancer_out.tsv", "leave_one_patient_out_robustness.tsv", "cancer_origin_association.tsv",
  "gate9b_endpoints.tsv", "package_versions.tsv", "GATE9B_FINAL_AUDIT.md", "gate9b_external_validation.log"
)
log_msg("Gate9B complete; decision=", gate_decision)
close(log_con); log_con <- NULL
old <- getwd(); setwd(args$out_dir); on.exit(setwd(old), add = TRUE)
hash <- system2("sha256sum", compact_files, stdout = TRUE, stderr = TRUE)
writeLines(hash, "SHA256SUMS")
quit(save = "no", status = 0L)
