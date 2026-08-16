#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(SingleR)
  library(BiocParallel)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("Usage: run_gate12r_oep_paired_axis1.R <raw_dir> <crosswalk.tsv> <transfer_reference.rds> <gate12g_model.rds> <outdir>")
}
raw_dir <- args[[1L]]
crosswalk_file <- args[[2L]]
reference_file <- args[[3L]]
model_file <- args[[4L]]
outdir <- args[[5L]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

seed <- 1208L
n_subsample <- 499L
subsample_fraction <- 0.80
set.seed(seed)

crosswalk <- fread(crosswalk_file)
required_crosswalk <- c("archive_directory", "node_sample_id", "patient_id", "cancer", "tissue_class")
missing_crosswalk <- setdiff(required_crosswalk, names(crosswalk))
if (length(missing_crosswalk)) stop("Crosswalk missing: ", paste(missing_crosswalk, collapse = ", "))
if (nrow(crosswalk) != 9L || uniqueN(crosswalk$patient_id) != 4L) stop("Expected nine specimens from four paired patients")
if (anyDuplicated(crosswalk$archive_directory) || anyDuplicated(crosswalk$node_sample_id)) stop("Crosswalk identifiers must be unique")
crosswalk[, `:=`(archive_dir = archive_directory, sample_id = node_sample_id, cancer_code = cancer)]

transfer <- readRDS(reference_file)
model <- readRDS(model_file)
if (!identical(transfer$status, "PASS")) stop("Transfer reference calibration did not PASS")
if (model$selected_k != 4L) stop("Frozen model must contain four selected axes")

read_feature_table <- function(archive) {
  path <- file.path(raw_dir, archive, "features.tsv.gz")
  z <- fread(path, header = FALSE)
  if (ncol(z) == 2L) z[, V3 := "Gene Expression"]
  setnames(z, names(z)[1:3], c("feature_id", "gene_symbol", "feature_type"))
  z
}

feature_tables <- lapply(crosswalk$archive_dir, read_feature_table)
symbol_lists <- lapply(feature_tables, function(z) z[feature_type == "Gene Expression", gene_symbol])
common_genes <- Reduce(intersect, lapply(symbol_lists, unique))
if (length(common_genes) < 20000L) stop("Unexpectedly small common gene universe")

make_test_logcpm <- function(counts, genes) {
  lib <- Matrix::colSums(counts)
  if (any(lib <= 0)) stop("QC matrix contains zero-count cells")
  idx <- match(genes, rownames(counts), nomatch = 0L)
  out <- matrix(0, nrow = length(genes), ncol = ncol(counts), dimnames = list(genes, colnames(counts)))
  present <- which(idx > 0L)
  if (length(present)) out[present, ] <- as.matrix(counts[idx[present], , drop = FALSE])
  log1p(t(t(out) / lib) * 1e6)
}

classify_taxonomy <- function(counts, ref_obj) {
  test <- make_test_logcpm(counts, ref_obj$marker_genes)
  pred <- SingleR(test = test, ref = ref_obj$logcpm, labels = ref_obj$meta$label,
                  de.method = "classic", fine.tune = TRUE, prune = TRUE,
                  BPPARAM = SerialParam())
  data.table(cell_id = colnames(counts), label = as.character(pred$labels),
             pruned = as.character(pred$pruned.labels), delta_next = as.numeric(pred$delta.next))
}

cell_rows <- vector("list", nrow(crosswalk))
library_rows <- vector("list", nrow(crosswalk))

for (i in seq_len(nrow(crosswalk))) {
  row <- crosswalk[i]
  archive <- row$archive_dir
  path <- file.path(raw_dir, archive)
  feat <- feature_tables[[i]]
  bc <- fread(file.path(path, "barcodes.tsv.gz"), header = FALSE)[[1L]]
  raw <- as(readMM(gzfile(file.path(path, "matrix.mtx.gz"))), "dgCMatrix")
  if (nrow(raw) != nrow(feat) || ncol(raw) != length(bc)) stop("Raw matrix mismatch for ", archive)
  gex_rows <- which(feat$feature_type == "Gene Expression")
  symbols <- feat$gene_symbol[gex_rows]
  first_symbol <- !duplicated(symbols)
  symbol_index <- match(common_genes, symbols[first_symbol])
  if (anyNA(symbol_index)) stop("Common-gene alignment failed for ", archive)
  original_index <- gex_rows[which(first_symbol)[symbol_index]]
  mat <- raw[original_index, , drop = FALSE]
  rownames(mat) <- common_genes
  colnames(mat) <- bc

  n_count <- as.numeric(Matrix::colSums(mat))
  n_feature <- as.numeric(Matrix::colSums(mat > 0))
  mito_idx <- grep("^MT-", rownames(mat))
  pct_mito <- if (length(mito_idx)) 100 * as.numeric(Matrix::colSums(mat[mito_idx, , drop = FALSE])) / pmax(n_count, 1) else rep(0, ncol(mat))
  qc_keep <- n_feature >= 200 & n_feature <= 8000 & n_count >= 500 & pct_mito <= 20
  mat <- mat[, qc_keep, drop = FALSE]
  cell_ids <- paste(archive, colnames(mat), sep = "::")
  colnames(mat) <- cell_ids

  broad_pred <- classify_taxonomy(mat, transfer$references$Broad)
  broad_pred[, broad := fifelse(is.na(pruned) | !nzchar(pruned), "Unassigned", pruned)]
  state <- rep(NA_character_, ncol(mat))
  state_delta <- rep(NA_real_, ncol(mat))
  for (lineage in c("Myeloid", "T_NK")) {
    idx_lineage <- which(broad_pred$broad == lineage)
    if (!length(idx_lineage)) next
    state_pred <- classify_taxonomy(mat[, idx_lineage, drop = FALSE], transfer$references[[lineage]])
    unresolved <- if (lineage == "Myeloid") "Unresolved_myeloid" else "Unresolved_T_NK"
    state[idx_lineage] <- fifelse(is.na(state_pred$pruned) | !nzchar(state_pred$pruned), unresolved, state_pred$pruned)
    state_delta[idx_lineage] <- state_pred$delta_next
  }

  cell_rows[[i]] <- data.table(
    cell_id = cell_ids, cell_barcode = bc[qc_keep], archive_dir = archive,
    sample_id = row$sample_id, patient_id = row$patient_id, cancer_code = row$cancer_code,
    tissue_class = row$tissue_class, gate12g_broad = broad_pred$broad,
    gate12g_broad_delta = broad_pred$delta_next, gate12g_state = state,
    gate12g_state_delta = state_delta
  )
  library_rows[[i]] <- data.table(
    archive_dir = archive, sample_id = row$sample_id, patient_id = row$patient_id,
    cancer_code = row$cancer_code, tissue_class = row$tissue_class,
    n_raw_cells = ncol(raw), n_qc_cells = sum(qc_keep),
    median_counts_qc = median(n_count[qc_keep]), median_features_qc = median(n_feature[qc_keep]),
    median_pct_mito_qc = median(pct_mito[qc_keep]), common_gene_universe = length(common_genes)
  )
  message("COMPLETED ", archive, " qc=", sum(qc_keep), " broad_assigned=", sprintf("%.3f", mean(broad_pred$broad != "Unassigned")))
  rm(raw, mat); gc(verbose = FALSE)
}

cells <- rbindlist(cell_rows)
library_audit <- rbindlist(library_rows)
fwrite(cells, file.path(outdir, "paired_gate12g_cell_assignments.tsv.gz"), sep = "\t")
fwrite(library_audit, file.path(outdir, "paired_library_audit.tsv"), sep = "\t")

broad_levels <- sub("^Broad__", "", colnames(model$x_raw_clr)[model$blocks == "Broad"])
my_levels <- sub("^Myeloid__", "", colnames(model$x_raw_clr)[model$blocks == "Myeloid"])
tn_levels <- sub("^T_NK__", "", colnames(model$x_raw_clr)[model$blocks == "T_NK"])

make_count_wide <- function(dt, label_col, levels) {
  base <- unique(crosswalk[, .(sample_id, archive_dir, patient_id, cancer_code, tissue_class)])
  obs <- dt[!is.na(get(label_col)), .N, by = .(sample_id, label = get(label_col))]
  grid <- base[, .(label = levels), by = .(sample_id, archive_dir, patient_id, cancer_code, tissue_class)]
  z <- merge(grid, obs, by = c("sample_id", "label"), all.x = TRUE)
  z[is.na(N), N := 0L]
  wide <- dcast(z, sample_id + archive_dir + patient_id + cancer_code + tissue_class ~ label,
                value.var = "N", fill = 0)
  for (lev in levels) if (!lev %chin% names(wide)) wide[, (lev) := 0L]
  setcolorder(wide, c("sample_id", "archive_dir", "patient_id", "cancer_code", "tissue_class", levels))
  setorder(wide, sample_id)
  wide
}

broad_w <- make_count_wide(cells, "gate12g_broad", broad_levels)
my_w <- make_count_wide(cells[gate12g_broad == "Myeloid"], "gate12g_state", my_levels)
tn_w <- make_count_wide(cells[gate12g_broad == "T_NK"], "gate12g_state", tn_levels)
if (!identical(broad_w$sample_id, my_w$sample_id) || !identical(broad_w$sample_id, tn_w$sample_id)) stop("Composition blocks are misaligned")

close_clr <- function(m) {
  p <- as.matrix(m) + 0.5
  p <- p / rowSums(p)
  lp <- log(p)
  lp - rowMeans(lp)
}

project_counts <- function(bmat, mmat, tmat) {
  xb <- close_clr(bmat); xm <- close_clr(mmat); xt <- close_clr(tmat)
  colnames(xb) <- paste0("Broad__", colnames(bmat))
  colnames(xm) <- paste0("Myeloid__", colnames(mmat))
  colnames(xt) <- paste0("T_NK__", colnames(tmat))
  x <- cbind(xb, xm, xt)[, colnames(model$x_raw_clr), drop = FALSE]
  for (block in unique(model$blocks)) {
    j <- which(model$blocks == block)
    scale_value <- sqrt(sum(apply(model$x_raw_clr[, j, drop = FALSE], 2L, var)))
    x[, j] <- x[, j, drop = FALSE] / scale_value
  }
  centered <- sweep(x, 2L, model$pca$center, "-")
  score <- centered %*% model$pca$rotation[, seq_len(model$selected_k), drop = FALSE]
  colnames(score) <- paste0("Axis", seq_len(ncol(score)))
  score
}

bmat <- as.matrix(broad_w[, ..broad_levels])
mmat <- as.matrix(my_w[, ..my_levels])
tmat <- as.matrix(tn_w[, ..tn_levels])
scores <- project_counts(bmat, mmat, tmat)
score_dt <- cbind(broad_w[, .(sample_id, archive_dir, patient_id, cancer_code, tissue_class)],
                  data.table(qc_cells = rowSums(bmat), broad_assigned_fraction = 1 - bmat[, "Unassigned"] / rowSums(bmat),
                             myeloid_cells = rowSums(mmat), t_nk_cells = rowSums(tmat)),
                  as.data.table(scores))
score_dt[, projectable := qc_cells >= 500 & myeloid_cells >= 50 & t_nk_cells >= 50]
fwrite(score_dt, file.path(outdir, "paired_axis_scores.tsv"), sep = "\t")

paired_summary <- function(dt, axis, control) {
  w <- dcast(dt[projectable == TRUE & tissue_class %chin% c("bone_metastasis", control)],
             patient_id + cancer_code ~ tissue_class, value.var = axis)
  if (!all(c("bone_metastasis", control) %chin% names(w))) return(list(summary = data.table(), detail = data.table()))
  w <- w[is.finite(bone_metastasis) & is.finite(get(control))]
  w[, `:=`(control_score = get(control), difference = bone_metastasis - get(control))]
  can_exact <- nrow(w) > 0L && all(w$difference != 0) && !anyDuplicated(abs(w$difference))
  p <- if (nrow(w)) suppressWarnings(wilcox.test(w$difference, alternative = "greater", exact = can_exact)$p.value) else NA_real_
  contrast <- paste0("bm_vs_", control)
  list(
    summary = data.table(axis = axis, contrast = contrast, n_pairs = nrow(w),
                         n_positive = sum(w$difference > 0), positive_fraction = mean(w$difference > 0),
                         median_difference = median(w$difference), mean_difference = mean(w$difference),
                         wilcoxon_one_sided_p = p, exact_test = can_exact),
    detail = w[, .(axis = axis, contrast = contrast, patient_id, cancer_code,
                   bone_metastasis_score = bone_metastasis, control_score, difference)]
  )
}

pair_summaries <- list(); pair_details <- list()
for (axis in paste0("Axis", 1:4)) {
  for (control in c("normal_bone", "primary")) {
    z <- paired_summary(score_dt, axis, control)
    pair_summaries[[length(pair_summaries) + 1L]] <- z$summary
    pair_details[[length(pair_details) + 1L]] <- z$detail
  }
}
pair_summary <- rbindlist(pair_summaries, fill = TRUE)
pair_detail <- rbindlist(pair_details, fill = TRUE)
fwrite(pair_summary, file.path(outdir, "paired_axis_contrasts.tsv"), sep = "\t")
fwrite(pair_detail, file.path(outdir, "paired_axis_differences.tsv"), sep = "\t")

draw_multihyper <- function(counts, fraction = 0.8) {
  counts <- as.integer(counts); total <- sum(counts); draw <- floor(total * fraction)
  ans <- integer(length(counts)); remaining_total <- total; remaining_draw <- draw
  if (length(counts) > 1L) for (i in seq_len(length(counts) - 1L)) {
    ans[i] <- rhyper(1L, counts[i], remaining_total - counts[i], remaining_draw)
    remaining_total <- remaining_total - counts[i]; remaining_draw <- remaining_draw - ans[i]
  }
  ans[length(counts)] <- remaining_draw
  ans
}

sub_rows <- vector("list", n_subsample)
set.seed(seed + 1L)
for (iteration in seq_len(n_subsample)) {
  bs_b <- t(apply(bmat, 1L, draw_multihyper, fraction = subsample_fraction))
  bs_m <- t(apply(mmat, 1L, draw_multihyper, fraction = subsample_fraction))
  bs_t <- t(apply(tmat, 1L, draw_multihyper, fraction = subsample_fraction))
  colnames(bs_b) <- colnames(bmat); colnames(bs_m) <- colnames(mmat); colnames(bs_t) <- colnames(tmat)
  sc <- project_counts(bs_b, bs_m, bs_t)
  z <- copy(score_dt[, .(sample_id, patient_id, cancer_code, tissue_class, projectable)])
  z[, Axis1 := sc[, 1L]]
  normal <- paired_summary(z, "Axis1", "normal_bone")$summary
  primary <- paired_summary(z, "Axis1", "primary")$summary
  sub_rows[[iteration]] <- data.table(
    iteration = iteration,
    normal_n_positive = normal$n_positive, normal_n_pairs = normal$n_pairs,
    normal_median_difference = normal$median_difference,
    normal_rule_pass = normal$n_pairs == 2L && normal$n_positive == 2L && normal$median_difference > 0,
    primary_n_positive = primary$n_positive, primary_n_pairs = primary$n_pairs,
    primary_median_difference = primary$median_difference,
    primary_rule_pass = primary$n_pairs == 3L && primary$n_positive >= 2L && primary$median_difference > 0
  )
}
subsampling <- rbindlist(sub_rows)
fwrite(subsampling, file.path(outdir, "paired_axis1_subsampling.tsv.gz"), sep = "\t")

axis1_normal <- pair_summary[axis == "Axis1" & contrast == "bm_vs_normal_bone"]
axis1_primary <- pair_summary[axis == "Axis1" & contrast == "bm_vs_primary"]
normal_full_pass <- nrow(axis1_normal) == 1L && axis1_normal$n_pairs == 2L && axis1_normal$n_positive == 2L && axis1_normal$median_difference > 0
primary_full_pass <- nrow(axis1_primary) == 1L && axis1_primary$n_pairs == 3L && axis1_primary$n_positive >= 2L && axis1_primary$median_difference > 0
normal_sub_pass_fraction <- mean(subsampling$normal_rule_pass)
primary_sub_pass_fraction <- mean(subsampling$primary_rule_pass)
overall <- normal_full_pass && primary_full_pass && normal_sub_pass_fraction >= 0.80 && primary_sub_pass_fraction >= 0.80
decision <- data.table(
  all_specimens_projectable = all(score_dt$projectable),
  normal_full_rule_pass = normal_full_pass, primary_full_rule_pass = primary_full_pass,
  normal_subsample_pass_fraction = normal_sub_pass_fraction,
  primary_subsample_pass_fraction = primary_sub_pass_fraction,
  overall_paired_axis1_status = ifelse(overall, "PASS", "STOP")
)
fwrite(decision, file.path(outdir, "paired_axis1_decision.tsv"), sep = "\t")

plot_dt <- pair_detail[axis == "Axis1"]
plot_long <- melt(plot_dt, id.vars = c("contrast", "patient_id", "cancer_code"),
                  measure.vars = c("control_score", "bone_metastasis_score"),
                  variable.name = "specimen", value.name = "Axis1")
plot_long[, specimen := factor(specimen, levels = c("control_score", "bone_metastasis_score"),
                               labels = c("Paired control", "Bone metastasis"))]
p <- ggplot(plot_long, aes(specimen, Axis1, group = patient_id, colour = cancer_code)) +
  geom_line(linewidth = 0.7, alpha = 0.8) + geom_point(size = 2.5) +
  facet_wrap(~contrast, scales = "free_y") + theme_bw(base_size = 11) +
  labs(title = paste0("Frozen Axis1 in OEP005136 paired specimens: ", decision$overall_paired_axis1_status),
       x = NULL, y = "Frozen Axis1", colour = "Cancer")
ggsave(file.path(outdir, "FigureR1_OEP_paired_Axis1.pdf"), p, width = 8, height = 4.8)
ggsave(file.path(outdir, "FigureR1_OEP_paired_Axis1.png"), p, width = 8, height = 4.8, dpi = 300, bg = "white")

checkpoint <- c(
  "# Gate12R OEP005136 paired Axis1 checkpoint", "",
  paste0("- Projectable specimens: ", sum(score_dt$projectable), "/9"),
  paste0("- BM versus normal-bone full-data rule: ", normal_full_pass),
  paste0("- BM versus primary full-data rule: ", primary_full_pass),
  paste0("- Normal-bone rule preserved in 80% subsamples: ", sprintf("%.3f", normal_sub_pass_fraction)),
  paste0("- Primary-tumour rule preserved in 80% subsamples: ", sprintf("%.3f", primary_sub_pass_fraction)),
  paste0("- Frozen paired Axis1 decision: **", decision$overall_paired_axis1_status, "**"), "",
  "The paired cohorts contain only two BM/normal-bone and three BM/primary comparisons. Directional support is an independently defined biological endpoint but is not a powered clinical validation."
)
writeLines(checkpoint, file.path(outdir, "GATE12R_OEP_PAIRED_AXIS1_CHECKPOINT.md"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
cat("GATE12R_OEP_PAIRED_AXIS1_STATUS=", decision$overall_paired_axis1_status, "\n", sep = "")
print(decision)
