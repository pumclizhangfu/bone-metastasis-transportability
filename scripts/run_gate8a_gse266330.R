#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(SingleCellExperiment)
  library(scuttle)
  library(SingleR)
  library(BiocParallel)
  library(ggplot2)
})

parse_args <- function(x) {
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    if (!startsWith(x[[i]], "--") || i == length(x)) stop("invalid arguments")
    out[[sub("^--", "", x[[i]])]] <- x[[i + 1L]]
    i <- i + 2L
  }
  out
}

`%||%` <- function(x, y) if (is.null(x)) y else x

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("raw-dir", "patient-info", "geo-crosswalk", "signatures", "markers", "output")
missing_args <- setdiff(required, names(cfg))
if (length(missing_args)) stop("missing arguments: ", paste(missing_args, collapse = ", "))
workers <- as.integer(cfg$workers %||% "8")
seed <- as.integer(cfg$seed %||% "20260806")
set.seed(seed)

raw_dir <- normalizePath(cfg$`raw-dir`, mustWork = TRUE)
output_dir <- cfg$output
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(output_dir, "gate8a_gse266330.log")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message", append = TRUE)
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

cat("Gate8A GSE266330 started", format(Sys.time()), "\n")
cat("raw_dir", raw_dir, "workers", workers, "seed", seed, "\n")

markers <- fread(cfg$markers)
signatures <- fread(cfg$signatures)
patient_info <- fread(cfg$`patient-info`, check.names = FALSE)
patient_info[, library_key := paste0(gsub("[^A-Za-z0-9]", "", cancer.id), as.integer(replicate))]
if (anyDuplicated(patient_info$library_key)) stop("patient metadata library keys are not unique")
geo_crosswalk <- fread(cfg$`geo-crosswalk`)
if (nrow(geo_crosswalk) != 63L || anyDuplicated(geo_crosswalk$geo_prefix) ||
    anyDuplicated(geo_crosswalk$metadata_library_key)) {
  stop("GEO crosswalk must contain 63 one-to-one library mappings")
}
if (!all(geo_crosswalk$title_key_matches_metadata) ||
    !all(geo_crosswalk$barcode_count == geo_crosswalk$expected_cells)) {
  stop("GEO crosswalk failed title or barcode-count validation")
}

matrix_files <- sort(list.files(raw_dir, pattern = "_matrix\\.mtx\\.gz$", full.names = TRUE))
if (length(matrix_files) != 63L) stop("expected 63 matrix files; found ", length(matrix_files))
prefix_paths <- sub("_matrix\\.mtx\\.gz$", "", matrix_files)
prefix_names <- basename(prefix_paths)
crosswalk_idx <- match(prefix_names, geo_crosswalk$geo_prefix)
if (anyNA(crosswalk_idx)) stop("archive prefixes are missing from GEO crosswalk")
library_keys <- geo_crosswalk$metadata_library_key[crosswalk_idx]
map_idx <- match(library_keys, patient_info$library_key)
if (anyNA(map_idx) || anyDuplicated(library_keys)) {
  stop("GEO library names do not map one-to-one to patient metadata")
}
library_meta <- copy(patient_info[map_idx])
library_meta[, `:=`(
  geo_prefix = prefix_names,
  library_key_from_file = library_keys,
  condition = ifelse(cancer == "ctrl", "healthy_bm", "bone_metastasis"),
  biological_sample = cancer.id
)]
library_meta[, `:=`(
  gsm = geo_crosswalk$gsm[crosswalk_idx],
  geo_file_key = geo_crosswalk$geo_file_key[crosswalk_idx],
  geo_sample_title = geo_crosswalk$sample_title[crosswalk_idx],
  barcode_count_expected = geo_crosswalk$barcode_count[crosswalk_idx]
)]
fwrite(
  geo_crosswalk[crosswalk_idx],
  file.path(output_dir, "geo_library_crosswalk_used.tsv"),
  sep = "\t"
)

feature_files <- paste0(prefix_paths, "_features.tsv.gz")
barcode_files <- paste0(prefix_paths, "_barcodes.tsv.gz")
if (!all(file.exists(feature_files)) || !all(file.exists(barcode_files))) {
  stop("missing feature or barcode files")
}

read_symbols <- function(path) {
  tab <- fread(path, header = FALSE, select = 2L, showProgress = FALSE)
  as.character(tab[[1L]])
}

cat("Reading feature universes", length(feature_files), "libraries\n")
symbol_lists <- lapply(feature_files, read_symbols)
common_genes <- Reduce(intersect, lapply(symbol_lists, unique))
feature_audit <- rbindlist(lapply(seq_along(symbol_lists), function(i) {
  data.table(
    geo_prefix = prefix_names[[i]],
    n_feature_rows = length(symbol_lists[[i]]),
    n_unique_symbols = uniqueN(symbol_lists[[i]]),
    n_duplicated_symbols = sum(duplicated(symbol_lists[[i]])),
    n_common_symbols = length(common_genes)
  )
}))
fwrite(feature_audit, file.path(output_dir, "feature_universe_audit.tsv"), sep = "\t")
if (length(common_genes) < 20000L) stop("unexpectedly small common gene universe: ", length(common_genes))
cat("Common gene universe", length(common_genes), "\n")

marker_sets <- split(markers$gene, markers$marker_group)
marker_sum <- function(mat, genes) {
  idx <- match(genes, rownames(mat), nomatch = 0L)
  idx <- idx[idx > 0L]
  if (!length(idx)) return(rep(0, ncol(mat)))
  as.numeric(Matrix::colSums(mat[idx, , drop = FALSE]))
}

cat("Loading Monaco immune reference\n")
if (!is.null(cfg$`reference-rds`)) {
  reference <- readRDS(cfg$`reference-rds`)
  cat("Loaded frozen Monaco reference RDS", cfg$`reference-rds`, "\n")
} else {
  if (!requireNamespace("celldex", quietly = TRUE)) {
    stop("celldex is required unless --reference-rds is supplied")
  }
  reference <- celldex::MonacoImmuneData()
}
param <- if (.Platform$OS.type == "unix" && workers > 1L) MulticoreParam(workers) else SerialParam()

map_lineage <- function(labels) {
  x <- tolower(labels)
  out <- rep(NA_character_, length(x))
  out[grepl("natural killer|(^|[^a-z])nk([^a-z]|$)|nkt", x)] <- "NK_NKT"
  out[is.na(out) & grepl("cd8|cytotoxic t", x)] <- "CD8_CTL"
  out[is.na(out) & grepl("cd4|treg|t regulatory|regulatory t|helper t|th1|th2|th17", x)] <- "CD4_T"
  out
}

candidate_obs <- vector("list", length(matrix_files))
hc_mats <- vector("list", length(matrix_files))
hc_obs_list <- vector("list", length(matrix_files))
library_audit <- vector("list", length(matrix_files))

for (i in seq_along(matrix_files)) {
  cat("READ_LIBRARY", i, library_meta$geo_prefix[[i]], "\n")
  symbols <- symbol_lists[[i]]
  barcodes <- as.character(fread(barcode_files[[i]], header = FALSE, showProgress = FALSE)[[1L]])
  if (length(barcodes) != library_meta$barcode_count_expected[[i]]) {
    stop("barcode count differs from frozen crosswalk for ", library_meta$geo_prefix[[i]])
  }
  mat <- as(readMM(gzfile(matrix_files[[i]])), "dgCMatrix")
  if (nrow(mat) != length(symbols) || ncol(mat) != length(barcodes)) stop("matrix dimension mismatch")
  first_symbol <- !duplicated(symbols)
  symbol_index <- match(common_genes, symbols[first_symbol])
  if (anyNA(symbol_index)) stop("common gene match failed")
  original_index <- which(first_symbol)[symbol_index]
  mat <- mat[original_index, , drop = FALSE]
  rownames(mat) <- common_genes
  colnames(mat) <- paste(library_meta$library_key[[i]], barcodes, sep = "::")

  n_count <- as.numeric(Matrix::colSums(mat))
  n_feature <- as.numeric(Matrix::colSums(mat > 0))
  mito_idx <- grep("^MT-", rownames(mat))
  pct_mito <- if (length(mito_idx)) {
    100 * as.numeric(Matrix::colSums(mat[mito_idx, , drop = FALSE])) / pmax(n_count, 1)
  } else rep(0, ncol(mat))
  qc_keep <- n_feature >= 200 & n_feature <= 8000 & n_count >= 500 & pct_mito <= 20
  mat_qc <- mat[, qc_keep, drop = FALSE]

  t_core <- marker_sum(mat_qc, marker_sets$T_core)
  cd4_support <- marker_sum(mat_qc, marker_sets$CD4_support)
  cd8_support <- marker_sum(mat_qc, marker_sets$CD8_support)
  nk_core <- marker_sum(mat_qc, marker_sets$NK_core)
  myeloid_audit <- marker_sum(mat_qc, marker_sets$Myeloid_audit)
  b_audit <- marker_sum(mat_qc, marker_sets$B_audit)
  epithelial_audit <- marker_sum(mat_qc, marker_sets$Epithelial_audit)
  candidate <- t_core >= 2 | nk_core >= 2

  obs <- data.table(
    cell_id = colnames(mat_qc),
    library_key = library_meta$library_key[[i]],
    biological_sample = library_meta$biological_sample[[i]],
    cancer = library_meta$cancer[[i]],
    condition = library_meta$condition[[i]],
    seq_id = library_meta$Seq.ID[[i]],
    n_count = n_count[qc_keep],
    n_feature = n_feature[qc_keep],
    pct_mito = pct_mito[qc_keep],
    t_core = t_core,
    cd4_support = cd4_support,
    cd8_support = cd8_support,
    nk_core = nk_core,
    myeloid_audit = myeloid_audit,
    b_audit = b_audit,
    epithelial_audit = epithelial_audit
  )
  cand_mat <- mat_qc[, candidate, drop = FALSE]
  cand_obs <- obs[candidate]
  if (ncol(cand_mat)) {
    cand_sce <- SingleCellExperiment(assays = list(counts = cand_mat), colData = as.data.frame(cand_obs))
    cand_sce <- logNormCounts(cand_sce)
    prediction <- SingleR(
      test = cand_sce,
      ref = reference,
      labels = reference$label.fine,
      assay.type.test = "logcounts",
      fine.tune = TRUE,
      prune = TRUE,
      BPPARAM = param
    )
    cand_obs[, `:=`(
      singler_label = as.character(prediction$labels),
      singler_pruned_label = as.character(prediction$pruned.labels),
      singler_delta_next = as.numeric(prediction$delta.next)
    )]
    cand_obs[, lineage_predicted := map_lineage(singler_pruned_label)]
    cand_obs[, canonical_support := fifelse(
      lineage_predicted == "CD4_T", t_core >= 2 & cd4_support >= 1,
      fifelse(
        lineage_predicted == "CD8_CTL", t_core >= 2 & cd8_support >= 1,
        fifelse(lineage_predicted == "NK_NKT", nk_core >= 2, FALSE)
      )
    )]
    cand_obs[, high_confidence := !is.na(lineage_predicted) & canonical_support]
    cand_obs[, lineage := fifelse(high_confidence, lineage_predicted, NA_character_)]
    hc_mats[[i]] <- cand_mat[, cand_obs$high_confidence, drop = FALSE]
    hc_obs_list[[i]] <- cand_obs[high_confidence == TRUE]
    rm(cand_sce, prediction)
  } else {
    cand_obs[, `:=`(
      singler_label = character(),
      singler_pruned_label = character(),
      singler_delta_next = numeric(),
      lineage_predicted = character(),
      canonical_support = logical(),
      high_confidence = logical(),
      lineage = character()
    )]
    hc_mats[[i]] <- cand_mat
    hc_obs_list[[i]] <- cand_obs
  }
  candidate_obs[[i]] <- cand_obs
  library_audit[[i]] <- data.table(
    geo_prefix = library_meta$geo_prefix[[i]],
    library_key = library_meta$library_key[[i]],
    biological_sample = library_meta$biological_sample[[i]],
    cancer = library_meta$cancer[[i]],
    condition = library_meta$condition[[i]],
    n_raw_cells = ncol(mat),
    n_qc_cells = sum(qc_keep),
    n_candidate_tnk = sum(candidate),
    n_high_confidence_tnk = sum(cand_obs$high_confidence),
    median_counts_qc = if (sum(qc_keep)) median(n_count[qc_keep]) else NA_real_,
    median_features_qc = if (sum(qc_keep)) median(n_feature[qc_keep]) else NA_real_,
    median_pct_mito_qc = if (sum(qc_keep)) median(pct_mito[qc_keep]) else NA_real_
  )
  cat(
    "ANNOTATED_LIBRARY", i,
    "candidates", ncol(cand_mat),
    "high_confidence", sum(cand_obs$high_confidence), "\n"
  )
  rm(mat, mat_qc, cand_mat, cand_obs)
  gc(verbose = FALSE)
}

library_audit <- rbindlist(library_audit)
fwrite(library_audit, file.path(output_dir, "input_library_audit.tsv"), sep = "\t")
candidate_obs <- rbindlist(candidate_obs)
cat(
  "Candidate/high-confidence T/NK cells",
  sum(library_audit$n_candidate_tnk),
  sum(library_audit$n_high_confidence_tnk), "\n"
)

label_audit <- candidate_obs[, .(
  n_cells = .N,
  n_canonical = sum(canonical_support),
  median_delta_next = median(singler_delta_next, na.rm = TRUE),
  median_t_core = median(t_core),
  median_cd4_support = median(cd4_support),
  median_cd8_support = median(cd8_support),
  median_nk_core = median(nk_core),
  median_myeloid_audit = median(myeloid_audit),
  median_b_audit = median(b_audit),
  median_epithelial_audit = median(epithelial_audit)
), by = .(singler_pruned_label, lineage_predicted)]
fwrite(label_audit, file.path(output_dir, "annotation_label_audit.tsv"), sep = "\t")
fwrite(candidate_obs, file.path(output_dir, "cell_annotations.tsv.gz"), sep = "\t", compress = "gzip")

hc_counts <- do.call(cbind, hc_mats)
hc_obs <- rbindlist(hc_obs_list)
stopifnot(identical(colnames(hc_counts), hc_obs$cell_id))
if (ncol(hc_counts) < 1000L) stop("too few high-confidence T/NK cells: ", ncol(hc_counts))
group_id <- paste(hc_obs$biological_sample, hc_obs$lineage, sep = "::")
groups <- sort(unique(group_id))
group_matrix <- sparseMatrix(
  i = seq_along(group_id),
  j = match(group_id, groups),
  x = 1,
  dims = c(length(group_id), length(groups)),
  dimnames = list(NULL, groups)
)
pseudobulk_counts <- as(hc_counts %*% group_matrix, "dgCMatrix")

group_split <- tstrsplit(groups, "::", fixed = TRUE)
pseudobulk_meta <- data.table(
  group_id = groups,
  biological_sample = group_split[[1L]],
  lineage = group_split[[2L]]
)
sample_meta <- unique(library_meta[, .(biological_sample, cancer, condition)])
sample_meta_index <- match(pseudobulk_meta$biological_sample, sample_meta$biological_sample)
if (anyNA(sample_meta_index)) stop("pseudobulk biological samples are missing from sample metadata")
pseudobulk_meta[, `:=`(
  cancer = sample_meta$cancer[sample_meta_index],
  condition = sample_meta$condition[sample_meta_index]
)]
cell_group_counts <- table(group_id)
pseudobulk_meta[, n_cells := as.integer(cell_group_counts[group_id])]
pseudobulk_meta[, total_umi := as.numeric(Matrix::colSums(pseudobulk_counts)[group_id])]
setcolorder(pseudobulk_meta, c("group_id", "biological_sample", "cancer", "condition", "lineage", "n_cells", "total_umi"))
if (!identical(colnames(pseudobulk_counts), pseudobulk_meta$group_id)) {
  stop("pseudobulk count columns and metadata rows are not identically ordered")
}
fwrite(pseudobulk_meta, file.path(output_dir, "pseudobulk_metadata.tsv"), sep = "\t")
saveRDS(pseudobulk_counts, file.path(output_dir, "pseudobulk_counts.rds"), compress = "xz")

signature_audit <- signatures[, {
  genes_present <- gene %in% rownames(pseudobulk_counts)
  .(
    n_frozen = .N,
    n_present = sum(genes_present),
    coverage = mean(genes_present),
    n_up_frozen = sum(direction == "up"),
    n_up_present = sum(direction == "up" & genes_present),
    up_coverage = sum(direction == "up" & genes_present) / sum(direction == "up"),
    n_down_frozen = sum(direction == "down"),
    n_down_present = sum(direction == "down" & genes_present),
    down_coverage = sum(direction == "down" & genes_present) / sum(direction == "down")
  )
}, by = .(signature, lineage)]
fwrite(signature_audit, file.path(output_dir, "signature_gene_audit.tsv"), sep = "\t")

library_sizes <- as.numeric(Matrix::colSums(pseudobulk_counts))
log_cpm <- log2(t(t(as.matrix(pseudobulk_counts)) / pmax(library_sizes, 1) * 1e6) + 0.5)
rank_matrix <- apply(log_cpm, 2L, function(x) rank(x, ties.method = "average") / length(x))
rownames(rank_matrix) <- rownames(pseudobulk_counts)

score_rows <- list()
score_i <- 0L
for (j in seq_len(nrow(pseudobulk_meta))) {
  lin <- pseudobulk_meta$lineage[[j]]
  wanted <- c("common_81", paste0("lineage_", lin))
  for (sig in wanted) {
    def <- signatures[signature == sig]
    up <- intersect(def[direction == "up", gene], rownames(rank_matrix))
    down <- intersect(def[direction == "down", gene], rownames(rank_matrix))
    audit <- signature_audit[signature == sig]
    minimum_coverage <- if (sig == "common_81") 0.80 else 0.70
    score_i <- score_i + 1L
    score_rows[[score_i]] <- cbind(
      pseudobulk_meta[j],
      data.table(
        signature = sig,
        score = mean(rank_matrix[up, j]) - mean(rank_matrix[down, j]),
        coverage = audit$coverage[[1L]],
        up_coverage = audit$up_coverage[[1L]],
        down_coverage = audit$down_coverage[[1L]],
        eligible_20 = pseudobulk_meta$n_cells[[j]] >= 20 &
          pseudobulk_meta$total_umi[[j]] >= 10000 &
          audit$coverage[[1L]] >= minimum_coverage &
          audit$up_coverage[[1L]] >= 0.70 &
          audit$down_coverage[[1L]] >= 0.70
      )
    )
  }
}
scores <- rbindlist(score_rows)
fwrite(scores, file.path(output_dir, "sample_lineage_scores.tsv"), sep = "\t")

two_sample_effect <- function(bm, ctrl) {
  if (length(bm) < 2L || length(ctrl) < 2L) {
    return(list(effect = NA_real_, p = NA_real_))
  }
  test <- suppressWarnings(wilcox.test(bm, ctrl, alternative = "greater", exact = FALSE))
  list(effect = median(as.vector(outer(bm, ctrl, "-"))), p = unname(test$p.value))
}

test_rows <- list()
for (sig in unique(scores$signature)) {
  for (lin in unique(scores[signature == sig, lineage])) {
    dat <- scores[signature == sig & lineage == lin & eligible_20]
    bm <- dat[condition == "bone_metastasis", score]
    ctrl <- dat[condition == "healthy_bm", score]
    effect <- two_sample_effect(bm, ctrl)
    test_rows[[length(test_rows) + 1L]] <- data.table(
      signature = sig,
      lineage = lin,
      n_bm = length(bm),
      n_control = length(ctrl),
      median_bm = if (length(bm)) median(bm) else NA_real_,
      median_control = if (length(ctrl)) median(ctrl) else NA_real_,
      hodges_lehmann_bm_minus_control = effect$effect,
      wilcoxon_one_sided_p = effect$p
    )
  }
}
tests <- rbindlist(test_rows)
tests[, endpoint_class := ifelse(signature == "common_81", "common", "lineage_specific")]
tests[, FDR := p.adjust(wilcoxon_one_sided_p, method = "BH"), by = endpoint_class]
fwrite(tests, file.path(output_dir, "cohort_effects.tsv"), sep = "\t")

common_scores <- scores[signature == "common_81" & eligible_20]
control_medians <- common_scores[condition == "healthy_bm", .(control_median = median(score)), by = lineage]
cancer_effects <- common_scores[condition == "bone_metastasis", .(
  n_bm = .N,
  median_bm = median(score)
), by = .(lineage, cancer)]
cancer_effects <- merge(cancer_effects, control_medians, by = "lineage", all.x = TRUE)
cancer_effects[, effect := median_bm - control_median]
cancer_effects[, eligible_cancer := n_bm >= 3L]
fwrite(cancer_effects, file.path(output_dir, "cancer_type_effects.tsv"), sep = "\t")

loco <- list()
bm_cancers <- sort(unique(common_scores[condition == "bone_metastasis", cancer]))
for (lin in unique(common_scores$lineage)) {
  ctrl <- common_scores[lineage == lin & condition == "healthy_bm", score]
  for (left_out in bm_cancers) {
    bm <- common_scores[lineage == lin & condition == "bone_metastasis" & cancer != left_out, score]
    effect <- two_sample_effect(bm, ctrl)
    loco[[length(loco) + 1L]] <- data.table(
      lineage = lin,
      left_out_cancer = left_out,
      n_bm = length(bm),
      n_control = length(ctrl),
      effect = effect$effect,
      p = effect$p,
      positive = effect$effect > 0
    )
  }
}
loco <- rbindlist(loco)
fwrite(loco, file.path(output_dir, "leave_one_cancer_out.tsv"), sep = "\t")

non_urologic <- list()
for (lin in unique(common_scores$lineage)) {
  ctrl <- common_scores[lineage == lin & condition == "healthy_bm", score]
  bm <- common_scores[lineage == lin & condition == "bone_metastasis" & !cancer %in% c("KC", "PC"), score]
  effect <- two_sample_effect(bm, ctrl)
  non_urologic[[length(non_urologic) + 1L]] <- data.table(
    lineage = lin,
    n_bm = length(bm),
    n_control = length(ctrl),
    effect = effect$effect,
    p = effect$p,
    positive = effect$effect > 0
  )
}
non_urologic <- rbindlist(non_urologic)
non_urologic[, FDR := p.adjust(p, method = "BH")]
fwrite(non_urologic, file.path(output_dir, "non_urologic_extension.tsv"), sep = "\t")

sensitivity <- list()
for (threshold in c(10L, 20L, 50L)) {
  for (lin in unique(scores$lineage)) {
    dat <- scores[signature == "common_81" & lineage == lin & n_cells >= threshold & total_umi >= 10000 &
      coverage >= 0.80 & up_coverage >= 0.70 & down_coverage >= 0.70]
    bm <- dat[condition == "bone_metastasis", score]
    ctrl <- dat[condition == "healthy_bm", score]
    effect <- two_sample_effect(bm, ctrl)
    sensitivity[[length(sensitivity) + 1L]] <- data.table(
      min_cells = threshold,
      lineage = lin,
      n_bm = length(bm),
      n_control = length(ctrl),
      effect = effect$effect,
      p = effect$p,
      positive = effect$effect > 0
    )
  }
}
sensitivity <- rbindlist(sensitivity)
fwrite(sensitivity, file.path(output_dir, "cell_threshold_sensitivity.tsv"), sep = "\t")

plot_data <- copy(common_scores)
plot_data[, condition_label := ifelse(condition == "healthy_bm", "Healthy BM", "Bone metastasis")]
p <- ggplot(plot_data, aes(x = condition_label, y = score, fill = condition_label)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.75) +
  facet_wrap(~ lineage, scales = "free_y") +
  scale_fill_manual(values = c("Healthy BM" = "#4C78A8", "Bone metastasis" = "#E45756")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none") +
  labs(x = NULL, y = "Frozen 81-gene directed rank score", title = "GSE266330 sample-level Gate8A validation")
ggsave(file.path(output_dir, "common_score_by_condition.pdf"), p, width = 9, height = 4.8)

common_test <- tests[signature == "common_81"]
lineage_specific_test <- tests[endpoint_class == "lineage_specific"]
positive_lineages <- sum(common_test$hodges_lehmann_bm_minus_control > 0, na.rm = TRUE)
significant_lineages <- sum(common_test$FDR < 0.10 & common_test$hodges_lehmann_bm_minus_control > 0, na.rm = TRUE)
lineage_specific_positive <- sum(lineage_specific_test$hodges_lehmann_bm_minus_control > 0, na.rm = TRUE)
lineage_specific_significant <- sum(
  lineage_specific_test$FDR < 0.10 & lineage_specific_test$hodges_lehmann_bm_minus_control > 0,
  na.rm = TRUE
)
eligible_cancer_summary <- cancer_effects[eligible_cancer == TRUE, .(
  positive_fraction = mean(effect > 0),
  n_cancers = .N
), by = lineage]
loco_summary <- loco[, .(positive_fraction = mean(positive), iterations = .N), by = lineage]
non_urologic_positive <- sum(non_urologic$positive, na.rm = TRUE)
cancer_consistent_lineages <- sum(eligible_cancer_summary$positive_fraction >= 0.75, na.rm = TRUE)
loco_stable_lineages <- sum(loco_summary$positive_fraction >= 0.80, na.rm = TRUE)
strong_provisional <- positive_lineages >= 2L && significant_lineages >= 2L &&
  non_urologic_positive >= 2L && cancer_consistent_lineages >= 2L &&
  loco_stable_lineages >= 2L && lineage_specific_positive >= 2L
directional_provisional <- positive_lineages >= 2L && non_urologic_positive >= 2L

decision <- c(
  "# Gate8A GSE266330 decision",
  "",
  paste0("- Common-score positive lineages: ", positive_lineages, "/3."),
  paste0("- Common-score FDR < 0.10 and positive: ", significant_lineages, "/3."),
  paste0("- Lineage-specific positive signatures: ", lineage_specific_positive, "/3; FDR < 0.10 and positive: ", lineage_specific_significant, "/3."),
  paste0("- Non-urologic positive lineages: ", non_urologic_positive, "/3."),
  paste0("- Eligible cancer-type positive fractions: ", paste(eligible_cancer_summary$lineage, sprintf("%.3f", eligible_cancer_summary$positive_fraction), collapse = "; "), "."),
  paste0("- Leave-one-cancer-out positive fractions: ", paste(loco_summary$lineage, sprintf("%.3f", loco_summary$positive_fraction), collapse = "; "), "."),
  paste0("- Cancer-consistent lineages (positive fraction >= 0.75): ", cancer_consistent_lineages, "/3."),
  paste0("- LOCO-stable lineages (positive fraction >= 0.80): ", loco_stable_lineages, "/3."),
  "",
  "This cohort alone cannot yield a full Gate8 GO because its healthy controls are run-confounded and Gate8B is required.",
  "",
  if (strong_provisional) {
    "**Decision: PROVISIONAL STRONG SUPPORT / PROCEED TO GATE8B.**"
  } else if (directional_provisional) {
    "**Decision: PROVISIONAL DIRECTIONAL SUPPORT / PROCEED TO GATE8B.**"
  } else {
    "**Decision: GATE8A DOES NOT SUPPORT THE FROZEN PROGRAM.**"
  }
)
writeLines(decision, file.path(output_dir, "gate8a_decision.md"))

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
cat("Gate8A completed", format(Sys.time()), "\n")
