#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(SingleR)
  library(BiocParallel)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(paste("Usage: run_gate12g_oep_projection.R <raw_dir> <crosswalk.tsv>",
             "<transfer_reference.rds> <gate12g_full.rds> <gate9b_cell_assignments.tsv.gz> <outdir>"))
}

raw_dir <- args[[1L]]
crosswalk_file <- args[[2L]]
reference_file <- args[[3L]]
model_file <- args[[4L]]
gate9b_cell_file <- args[[5L]]
outdir <- args[[6L]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

seed <- 1208L
n_subsample <- 499L
subsample_fraction <- 0.80
set.seed(seed)

crosswalk <- fread(crosswalk_file)
crosswalk[, primary_analysis := as.logical(primary_analysis)]
crosswalk <- crosswalk[primary_analysis == TRUE]
if (nrow(crosswalk) != 49L || uniqueN(crosswalk$cancer_code) != 11L) {
  stop("Frozen external cohort must contain 49 patients from 11 cancer origins")
}
transfer <- readRDS(reference_file)
if (!identical(transfer$status, "PASS")) stop("Discovery reference calibration did not PASS")
model <- readRDS(model_file)
if (model$selected_k != 4L) stop("Frozen Gate12G model must have selected K=4")
old_cells <- fread(gate9b_cell_file, select = c("cell_id", "archive_dir", "cell_barcode",
                                                "patient_id", "cancer_code", "broad_lineage",
                                                "final_assignment"))
old_cells <- old_cells[archive_dir %chin% crosswalk$archive_dir]
setkey(old_cells, archive_dir, cell_barcode)

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
  if (any(lib <= 0)) stop("External QC matrix contains zero-count cells")
  idx <- match(genes, rownames(counts), nomatch = 0L)
  out <- matrix(0, nrow = length(genes), ncol = ncol(counts),
                dimnames = list(genes, colnames(counts)))
  present <- which(idx > 0L)
  if (length(present)) out[present, ] <- as.matrix(counts[idx[present], , drop = FALSE])
  log1p(t(t(out) / lib) * 1e6)
}

classify_taxonomy <- function(counts, ref_obj) {
  genes <- ref_obj$marker_genes
  test <- make_test_logcpm(counts, genes)
  pred <- SingleR(test = test, ref = ref_obj$logcpm, labels = ref_obj$meta$label,
                  de.method = "classic", fine.tune = TRUE, prune = TRUE,
                  BPPARAM = SerialParam())
  data.table(cell_id = colnames(counts),
             label = as.character(pred$labels),
             pruned = as.character(pred$pruned.labels),
             delta_next = as.numeric(pred$delta.next))
}

process_library <- function(row) {
  archive <- row$archive_dir
  path <- file.path(raw_dir, archive)
  files <- c(matrix = file.path(path, "matrix.mtx.gz"),
             features = file.path(path, "features.tsv.gz"),
             barcodes = file.path(path, "barcodes.tsv.gz"))
  if (!all(file.exists(files))) stop("Incomplete raw matrix for ", archive)
  feat <- fread(files[["features"]], header = FALSE)
  if (ncol(feat) == 2L) feat[, V3 := "Gene Expression"]
  setnames(feat, names(feat)[1:3], c("feature_id", "gene_symbol", "feature_type"))
  bc <- fread(files[["barcodes"]], header = FALSE)[[1L]]
  raw <- as(Matrix::readMM(gzfile(files[["matrix"]])), "dgCMatrix")
  colnames(raw) <- bc
  if (nrow(raw) != nrow(feat) || ncol(raw) != length(bc)) stop("Raw dimension mismatch for ", archive)
  gex_idx <- which(feat$feature_type == "Gene Expression")
  gex <- aggregate_symbols(raw[gex_idx, , drop = FALSE], feat$gene_symbol[gex_idx])

  old <- old_cells[J(archive)]
  idx <- match(old$cell_barcode, colnames(gex))
  if (anyNA(idx)) stop("Frozen Gate9B QC barcodes missing from raw matrix for ", archive)
  gex <- gex[, idx, drop = FALSE]
  cell_ids <- paste(archive, old$cell_barcode, sep = "::")
  colnames(gex) <- cell_ids
  if (!identical(cell_ids, old$cell_id)) stop("Frozen cell-id mismatch for ", archive)

  broad_pred <- classify_taxonomy(gex, transfer$references$Broad)
  broad_pred[, broad := fifelse(is.na(pruned) | !nzchar(pruned), "Unassigned", pruned)]
  state <- rep(NA_character_, ncol(gex))
  state_delta <- rep(NA_real_, ncol(gex))

  for (lineage in c("Myeloid", "T_NK")) {
    idx_lineage <- which(broad_pred$broad == lineage)
    if (!length(idx_lineage)) next
    state_pred <- classify_taxonomy(gex[, idx_lineage, drop = FALSE], transfer$references[[lineage]])
    unresolved <- if (lineage == "Myeloid") "Unresolved_myeloid" else "Unresolved_T_NK"
    state_value <- fifelse(is.na(state_pred$pruned) | !nzchar(state_pred$pruned), unresolved,
                           state_pred$pruned)
    state[idx_lineage] <- state_value
    state_delta[idx_lineage] <- state_pred$delta_next
  }

  ans <- data.table(cell_id = cell_ids, cell_barcode = old$cell_barcode,
                    archive_dir = archive, sample_id = row$sample_id,
                    patient_id = row$patient_id, cancer_code = row$cancer_code,
                    gate12g_broad = broad_pred$broad,
                    gate12g_broad_delta = broad_pred$delta_next,
                    gate12g_state = state, gate12g_state_delta = state_delta,
                    gate9b_broad_lineage = old$broad_lineage,
                    gate9b_final_assignment = old$final_assignment)
  message("COMPLETED ", archive, " cells=", nrow(ans),
          " broad_assigned=", sprintf("%.3f", mean(ans$gate12g_broad != "Unassigned")))
  ans
}

assignments <- vector("list", nrow(crosswalk))
for (i in seq_len(nrow(crosswalk))) assignments[[i]] <- process_library(crosswalk[i])
cells <- rbindlist(assignments)
if (nrow(cells) != nrow(old_cells)) stop("External cell count differs from frozen Gate9B QC universe")
fwrite(cells, file.path(outdir, "oep005136_gate12g_cell_assignments.tsv.gz"), sep = "\t")

broad_levels <- sub("^Broad__", "", colnames(model$x_raw_clr)[model$blocks == "Broad"])
my_levels <- sub("^Myeloid__", "", colnames(model$x_raw_clr)[model$blocks == "Myeloid"])
tn_levels <- sub("^T_NK__", "", colnames(model$x_raw_clr)[model$blocks == "T_NK"])

make_count_wide <- function(dt, label_col, levels, prefix = "") {
  base <- unique(crosswalk[, .(patient_id, sample_id, archive_dir, cancer_code)])
  obs <- dt[!is.na(get(label_col)), .N, by = .(patient_id, label = get(label_col))]
  grid <- base[, .(label = levels), by = .(patient_id, sample_id, archive_dir, cancer_code)]
  z <- merge(grid, obs, by = c("patient_id", "label"), all.x = TRUE)
  z[is.na(N), N := 0L]
  wide <- dcast(z, patient_id + sample_id + archive_dir + cancer_code ~ label,
                value.var = "N", fill = 0)
  for (lev in levels) if (!lev %chin% names(wide)) wide[, (lev) := 0L]
  setcolorder(wide, c("patient_id", "sample_id", "archive_dir", "cancer_code", levels))
  wide
}

broad_w <- make_count_wide(cells, "gate12g_broad", broad_levels)
my_w <- make_count_wide(cells[gate12g_broad == "Myeloid"], "gate12g_state", my_levels)
tn_w <- make_count_wide(cells[gate12g_broad == "T_NK"], "gate12g_state", tn_levels)
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
  xb <- close_clr(bmat)
  xm <- close_clr(mmat)
  xt <- close_clr(tmat)
  colnames(xb) <- paste0("Broad__", colnames(bmat))
  colnames(xm) <- paste0("Myeloid__", colnames(mmat))
  colnames(xt) <- paste0("T_NK__", colnames(tmat))
  x <- cbind(xb, xm, xt)
  x <- x[, colnames(model$x_raw_clr), drop = FALSE]
  block_scales <- setNames(numeric(length(unique(model$blocks))), unique(model$blocks))
  for (block in unique(model$blocks)) {
    j <- which(model$blocks == block)
    block_scales[[block]] <- sqrt(sum(apply(model$x_raw_clr[, j, drop = FALSE], 2L, var)))
    x[, j] <- x[, j, drop = FALSE] / block_scales[[block]]
  }
  centered <- sweep(x, 2L, model$pca$center, "-")
  score <- centered %*% model$pca$rotation[, seq_len(model$selected_k), drop = FALSE]
  colnames(score) <- paste0("Axis", seq_len(ncol(score)))
  score
}

bmat <- as.matrix(broad_w[, ..broad_levels])
mmat <- as.matrix(my_w[, ..my_levels])
tmat <- as.matrix(tn_w[, ..tn_levels])
full_scores <- project_counts(bmat, mmat, tmat)

patient_qc <- data.table(patient_id = broad_w$patient_id, sample_id = broad_w$sample_id,
                         archive_dir = broad_w$archive_dir, cancer_code = broad_w$cancer_code,
                         qc_cells = rowSums(bmat),
                         broad_assigned_fraction = 1 - bmat[, "Unassigned"] / rowSums(bmat),
                         myeloid_cells = rowSums(mmat), t_nk_cells = rowSums(tmat))
patient_qc[, projectable := qc_cells >= 500 & myeloid_cells >= 50 & t_nk_cells >= 50]
score_dt <- cbind(patient_qc, as.data.table(full_scores))
fwrite(score_dt, file.path(outdir, "oep005136_ecological_axis_scores.tsv"), sep = "\t")

# Transfer-reliable loading coverage.
class_metrics <- transfer$class_metrics
feature_reliability <- data.table(feature = colnames(model$x_raw_clr), block = model$blocks)
feature_reliability[, label := sub("^(Broad|Myeloid|T_NK)__", "", feature)]
feature_reliability[, taxonomy := block]
feature_reliability[, rejection_bin := label %chin% c("Unassigned", "Unresolved_myeloid", "Unresolved_T_NK")]
feature_reliability <- merge(feature_reliability,
  class_metrics[, .(taxonomy, label = truth, evaluable_patients, recall)],
  by = c("taxonomy", "label"), all.x = TRUE)
feature_reliability[, reliable := rejection_bin |
                      (!is.na(recall) & evaluable_patients >= 3L & recall >= 0.50)]
loading <- model$pca$rotation[, seq_len(model$selected_k), drop = FALSE]
coverage <- rbindlist(lapply(seq_len(model$selected_k), function(a) {
  data.table(axis = a,
             reliable_loading_coverage = sum(loading[feature_reliability$reliable, a]^2) /
               sum(loading[, a]^2))
}))

# Frozen sensitivity: collapse unreliable resolved states into lineage rejection bins.
bmat_s <- bmat; mmat_s <- mmat; tmat_s <- tmat
unreliable <- feature_reliability[reliable == FALSE]
for (i in seq_len(nrow(unreliable))) {
  lab <- unreliable$label[i]
  tax <- unreliable$taxonomy[i]
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
sensitivity_scores <- project_counts(bmat_s, mmat_s, tmat_s)
sensitivity <- data.table(axis = seq_len(model$selected_k),
  score_correlation = vapply(seq_len(model$selected_k), function(a) {
    cor(full_scores[, a], sensitivity_scores[, a], method = "spearman")
  }, numeric(1)))

draw_multihyper <- function(counts, fraction = 0.8) {
  counts <- as.integer(counts)
  draw_total <- floor(sum(counts) * fraction)
  ans <- integer(length(counts))
  remaining_total <- sum(counts)
  remaining_draw <- draw_total
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
boot_scores <- vector("list", n_subsample)
subsampling_rows <- vector("list", n_subsample)
set.seed(seed + 1L)
for (b in seq_len(n_subsample)) {
  bs_b <- t(apply(bmat, 1L, draw_multihyper, fraction = subsample_fraction))
  bs_m <- t(apply(mmat, 1L, draw_multihyper, fraction = subsample_fraction))
  bs_t <- t(apply(tmat, 1L, draw_multihyper, fraction = subsample_fraction))
  colnames(bs_b) <- colnames(bmat); colnames(bs_m) <- colnames(mmat); colnames(bs_t) <- colnames(tmat)
  sc <- project_counts(bs_b, bs_m, bs_t)
  boot_scores[[b]] <- sc
  subsampling_rows[[b]] <- data.table(iteration = b, axis = seq_len(model$selected_k),
    score_correlation = vapply(seq_len(model$selected_k), function(a) {
      cor(full_scores[projectable_idx, a], sc[projectable_idx, a], method = "spearman")
    }, numeric(1)))
}
subsampling <- rbindlist(subsampling_rows)
subsampling_summary <- subsampling[, .(
  median_score_correlation = median(score_correlation),
  q10_score_correlation = quantile(score_correlation, 0.10)
), by = axis]

eta_squared <- function(y, g) {
  ok <- is.finite(y) & !is.na(g)
  y <- y[ok]; g <- droplevels(factor(g[ok]))
  grand <- mean(y)
  ss_total <- sum((y - grand)^2)
  ss_between <- sum(vapply(split(y, g), function(z) length(z) * (mean(z) - grand)^2, numeric(1)))
  if (ss_total == 0) NA_real_ else ss_between / ss_total
}

axis_gate_for_subset <- function(idx, axis) {
  y <- full_scores[idx, axis]
  boot_cor <- vapply(boot_scores, function(z) cor(y, z[idx, axis], method = "spearman"), numeric(1))
  data.table(
    n_patients = length(idx),
    n_origins = uniqueN(patient_qc$cancer_code[idx]),
    median_subsample_correlation = median(boot_cor, na.rm = TRUE),
    rho_log10_qc_cells = cor(y, log10(patient_qc$qc_cells[idx]), method = "spearman"),
    rho_assignment_fraction = cor(y, patient_qc$broad_assigned_fraction[idx], method = "spearman"),
    cancer_eta2 = eta_squared(y, patient_qc$cancer_code[idx])
  )
}

gate_rows <- list()
for (a in seq_len(model$selected_k)) {
  base <- axis_gate_for_subset(projectable_idx, a)
  base[, `:=`(axis = a, deleted_origin = "NONE")]
  gate_rows[[length(gate_rows) + 1L]] <- base
  for (origin in unique(patient_qc$cancer_code[projectable_idx])) {
    idx <- projectable_idx[patient_qc$cancer_code[projectable_idx] != origin]
    z <- axis_gate_for_subset(idx, a)
    z[, `:=`(axis = a, deleted_origin = origin)]
    gate_rows[[length(gate_rows) + 1L]] <- z
  }
}
origin_deletion <- rbindlist(gate_rows)
origin_deletion[, pass := n_patients >= 30 & n_origins >= 6 &
                  median_subsample_correlation >= 0.80 &
                  abs(rho_log10_qc_cells) < 0.50 &
                  abs(rho_assignment_fraction) < 0.50 & cancer_eta2 < 0.50]

axis_gate <- merge(coverage, subsampling_summary, by = "axis")
axis_gate <- merge(axis_gate, sensitivity, by = "axis")
base_metrics <- origin_deletion[deleted_origin == "NONE",
  .(axis, n_projectable_patients = n_patients, n_projectable_origins = n_origins,
    rho_log10_qc_cells, rho_assignment_fraction, cancer_eta2)]
axis_gate <- merge(axis_gate, base_metrics, by = "axis")
deletion_summary <- origin_deletion[deleted_origin != "NONE",
                                    .(all_origin_deletions_pass = all(pass)), by = axis]
axis_gate <- merge(axis_gate, deletion_summary, by = "axis")
axis_gate[, discovery_eligible := axis %in% model$axis_summary[eligible == TRUE, axis]]
axis_gate[, externally_projectable := discovery_eligible &
             n_projectable_patients >= 30 & n_projectable_origins >= 6 &
             reliable_loading_coverage >= 0.80 &
             median_score_correlation >= 0.80 &
             abs(rho_log10_qc_cells) < 0.50 &
             abs(rho_assignment_fraction) < 0.50 & cancer_eta2 < 0.50 &
             score_correlation >= 0.80 & all_origin_deletions_pass]

# Secondary concordance with frozen Gate9B calls; not an independent endpoint.
concordance <- data.table(
  comparison = c("new_Myeloid_vs_gate9b_myeloid", "new_T_NK_vs_gate9b_t_cell",
                 "new_Classical_vs_gate9b_CD14HI", "new_Treg_vs_gate9b_CD4_TREG",
                 "new_CD8_exhausted_vs_gate9b_CD8_TEX"),
  denominator_cells = c(
    sum(cells$gate9b_broad_lineage == "myeloid"),
    sum(cells$gate9b_broad_lineage == "t_cell"),
    sum(cells$gate9b_final_assignment == "CD14HI_MONO"),
    sum(cells$gate9b_final_assignment == "CD4_TREG"),
    sum(cells$gate9b_final_assignment == "CD8_TEX")
  ),
  agreement_fraction = c(
    mean(cells$gate12g_broad[cells$gate9b_broad_lineage == "myeloid"] == "Myeloid"),
    mean(cells$gate12g_broad[cells$gate9b_broad_lineage == "t_cell"] == "T_NK"),
    mean(cells$gate12g_state[cells$gate9b_final_assignment == "CD14HI_MONO"] == "Classical_monocyte", na.rm = TRUE),
    mean(cells$gate12g_state[cells$gate9b_final_assignment == "CD4_TREG"] == "Treg", na.rm = TRUE),
    mean(cells$gate12g_state[cells$gate9b_final_assignment == "CD8_TEX"] == "CD8_exhausted", na.rm = TRUE)
  )
)

fwrite(patient_qc, file.path(outdir, "oep005136_projection_qc.tsv"), sep = "\t")
fwrite(feature_reliability, file.path(outdir, "transfer_feature_reliability.tsv"), sep = "\t")
fwrite(coverage, file.path(outdir, "axis_reliable_loading_coverage.tsv"), sep = "\t")
fwrite(sensitivity, file.path(outdir, "unreliable_state_collapse_sensitivity.tsv"), sep = "\t")
fwrite(subsampling, file.path(outdir, "cell_subsampling_stability.tsv.gz"), sep = "\t")
fwrite(origin_deletion, file.path(outdir, "leave_one_origin_out_gate.tsv"), sep = "\t")
fwrite(axis_gate, file.path(outdir, "external_axis_gate.tsv"), sep = "\t")
fwrite(concordance, file.path(outdir, "gate9b_secondary_concordance.tsv"), sep = "\t")
saveRDS(list(scores = score_dt, axis_gate = axis_gate, patient_qc = patient_qc,
             feature_reliability = feature_reliability, seed = seed,
             n_subsample = n_subsample),
        file.path(outdir, "gate12g_oep005136_projection.rds"), compress = "xz")

plot_scores <- melt(score_dt[projectable == TRUE],
                    id.vars = c("patient_id", "cancer_code", "qc_cells", "broad_assigned_fraction"),
                    measure.vars = paste0("Axis", model$axis_summary[eligible == TRUE, axis]),
                    variable.name = "axis_name", value.name = "score")
p1 <- ggplot(score_dt[projectable == TRUE], aes(Axis1, Axis2, colour = cancer_code)) +
  geom_point(size = 2, alpha = 0.85) + theme_bw(base_size = 10) +
  labs(title = "A  Frozen ecological-axis projection", colour = "Origin")
p2 <- ggplot(plot_scores, aes(cancer_code, score, colour = cancer_code)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.12) + geom_jitter(width = 0.12, size = 1) +
  facet_wrap(~axis_name, scales = "free_y") + theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(title = "B  Patient scores across cancer origins", x = NULL, y = "Frozen score")
p3 <- ggplot(subsampling[axis %in% model$axis_summary[eligible == TRUE, axis]],
             aes(factor(axis), score_correlation)) +
  geom_violin(fill = "#80B1D3", colour = "#377EB8") +
  geom_hline(yintercept = 0.80, linetype = 2, colour = "#B2182B") +
  coord_cartesian(ylim = c(0, 1)) + theme_bw(base_size = 10) +
  labs(title = "C  80% cell-subsampling stability", x = "Axis", y = "Patient-rank correlation")
p4 <- ggplot(plot_scores, aes(broad_assigned_fraction, score, colour = cancer_code)) +
  geom_point(size = 1.8, alpha = 0.8) + facet_wrap(~axis_name, scales = "free_y") +
  theme_bw(base_size = 10) + theme(legend.position = "none") +
  labs(title = "D  Assignment-fraction diagnostic", x = "Broad assigned fraction", y = "Frozen score")
fig <- (p1 | p2) / (p3 | p4) + plot_annotation(title = "Gate12G OEP005136 frozen external projection")
ggsave(file.path(outdir, "Figure7_gate12g_oep_projection.pdf"), fig, width = 14, height = 10)
ggsave(file.path(outdir, "Figure7_gate12g_oep_projection.png"), fig, width = 14, height = 10,
       dpi = 300, bg = "white")

passing <- axis_gate[externally_projectable == TRUE, paste0("Axis", axis, collapse = ", ")]
if (!nzchar(passing)) passing <- "none"
overall <- if (axis_gate[discovery_eligible == TRUE, all(externally_projectable)]) "PASS" else
  if (axis_gate[discovery_eligible == TRUE, any(externally_projectable)]) "PARTIAL" else "STOP"
checkpoint <- c(
  "# Gate12G OEP005136 external-projection checkpoint",
  "",
  paste0("- Frozen primary patients: ", nrow(crosswalk)),
  paste0("- Projectable patients: ", sum(patient_qc$projectable), "/", nrow(patient_qc)),
  paste0("- Projectable cancer origins: ", uniqueN(patient_qc[projectable == TRUE, cancer_code])),
  paste0("- Cell-subsampling iterations: ", n_subsample),
  paste0("- Externally projectable discovery axes: ", passing),
  paste0("- Overall decision: **", overall, "**"),
  "",
  "Projection used frozen discovery CLR transforms, block scales, PCA centre and loadings. No OEP-derived refitting, rotation, K selection or origin exclusion was performed.",
  "",
  "OEP005136 contains bone metastases only. PASS supports transportability of ecological coordinates, not anatomical progression, prognosis or a clinical subtype."
)
writeLines(checkpoint, file.path(outdir, "GATE12G_OEP_PROJECTION_CHECKPOINT.md"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))

cat("GATE12G_OEP_PROJECTION_STATUS=", overall, "\n", sep = "")
cat("EXTERNALLY_PROJECTABLE_AXES=", passing, "\n", sep = "")
print(axis_gate)
