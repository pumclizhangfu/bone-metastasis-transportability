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
required <- c("raw-dir", "crosswalk", "signatures", "markers", "output")
missing_args <- setdiff(required, names(cfg))
if (length(missing_args)) stop("missing arguments: ", paste(missing_args, collapse = ", "))
workers <- as.integer(cfg$workers %||% "8")
seed <- as.integer(cfg$seed %||% "20260806")
set.seed(seed)

raw_dir <- normalizePath(cfg$`raw-dir`, mustWork = TRUE)
output_dir <- cfg$output
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(output_dir, "gate8b_oep005136.log")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message", append = TRUE)
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

cat("Gate8B OEP005136 started", format(Sys.time()), "\n")
cat("raw_dir", raw_dir, "workers", workers, "seed", seed, "\n")

markers <- fread(cfg$markers)
signatures <- fread(cfg$signatures)
crosswalk <- fread(cfg$crosswalk)
expected_signatures <- c("common_81", "lineage_CD4_T", "lineage_CD8_CTL", "lineage_NK_NKT")
if (!setequal(unique(signatures$signature), expected_signatures)) {
  stop("frozen signature set is unexpected")
}
if (nrow(crosswalk) != 9L || anyDuplicated(crosswalk$archive_directory) ||
    anyDuplicated(crosswalk$node_sample_id)) {
  stop("Gate8B archive crosswalk must contain nine one-to-one specimens")
}
expected_tissues <- c(
  LUCA1 = "bone_metastasis,normal_bone,primary",
  LUCA2 = "bone_metastasis,normal_bone",
  PRAD1 = "bone_metastasis,primary",
  PRAD2 = "bone_metastasis,primary"
)
observed_tissues <- crosswalk[order(tissue_class), paste(tissue_class, collapse = ","), by = patient_id]
for (patient in names(expected_tissues)) {
  actual <- observed_tissues[patient_id == patient, V1]
  if (length(actual) != 1L || actual != expected_tissues[[patient]]) {
    stop("unexpected tissue set for ", patient, ": ", paste(actual, collapse = ","))
  }
}

sample_dirs <- sort(list.dirs(raw_dir, recursive = FALSE, full.names = TRUE))
dir_names <- basename(sample_dirs)
if (!setequal(dir_names, crosswalk$archive_directory) || length(sample_dirs) != 9L) {
  stop("extracted directories do not match the frozen nine-row crosswalk")
}
sample_dirs <- sample_dirs[match(crosswalk$archive_directory, dir_names)]
matrix_files <- file.path(sample_dirs, "matrix.mtx.gz")
feature_files <- file.path(sample_dirs, "features.tsv.gz")
barcode_files <- file.path(sample_dirs, "barcodes.tsv.gz")
if (!all(file.exists(matrix_files)) || !all(file.exists(feature_files)) || !all(file.exists(barcode_files))) {
  stop("one or more 10X triplets are incomplete")
}

read_symbols <- function(path) {
  tab <- fread(path, header = FALSE, select = 2L, showProgress = FALSE)
  as.character(tab[[1L]])
}

cat("Reading feature universes from", length(feature_files), "libraries\n")
symbol_lists <- lapply(feature_files, read_symbols)
common_genes <- Reduce(intersect, lapply(symbol_lists, unique))
feature_audit <- rbindlist(lapply(seq_along(symbol_lists), function(i) {
  data.table(
    archive_directory = crosswalk$archive_directory[[i]],
    node_sample_id = crosswalk$node_sample_id[[i]],
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

cat("Loading frozen Monaco immune reference\n")
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
  cat("READ_LIBRARY", i, crosswalk$archive_directory[[i]], "\n")
  symbols <- symbol_lists[[i]]
  barcodes <- as.character(fread(barcode_files[[i]], header = FALSE, showProgress = FALSE)[[1L]])
  if (anyDuplicated(barcodes)) stop("duplicated barcodes within ", crosswalk$archive_directory[[i]])
  mat <- as(readMM(gzfile(matrix_files[[i]])), "dgCMatrix")
  if (nrow(mat) != length(symbols) || ncol(mat) != length(barcodes)) stop("matrix dimension mismatch")
  first_symbol <- !duplicated(symbols)
  symbol_index <- match(common_genes, symbols[first_symbol])
  if (anyNA(symbol_index)) stop("common gene match failed")
  original_index <- which(first_symbol)[symbol_index]
  mat <- mat[original_index, , drop = FALSE]
  rownames(mat) <- common_genes
  colnames(mat) <- paste(crosswalk$node_sample_id[[i]], barcodes, sep = "::")

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
    archive_directory = crosswalk$archive_directory[[i]],
    specimen_id = crosswalk$node_sample_id[[i]],
    patient_id = crosswalk$patient_id[[i]],
    cancer = crosswalk$cancer[[i]],
    tissue_class = crosswalk$tissue_class[[i]],
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
    archive_directory = crosswalk$archive_directory[[i]],
    specimen_id = crosswalk$node_sample_id[[i]],
    patient_id = crosswalk$patient_id[[i]],
    cancer = crosswalk$cancer[[i]],
    tissue_class = crosswalk$tissue_class[[i]],
    n_raw_cells = ncol(mat),
    n_qc_cells = sum(qc_keep),
    n_candidate_tnk = sum(candidate),
    n_high_confidence_tnk = sum(cand_obs$high_confidence),
    median_counts_qc = if (sum(qc_keep)) median(n_count[qc_keep]) else NA_real_,
    median_features_qc = if (sum(qc_keep)) median(n_feature[qc_keep]) else NA_real_,
    median_pct_mito_qc = if (sum(qc_keep)) median(pct_mito[qc_keep]) else NA_real_
  )
  cat("ANNOTATED_LIBRARY", i, "candidates", ncol(cand_mat),
      "high_confidence", sum(cand_obs$high_confidence), "\n")
  rm(mat, mat_qc, cand_mat, cand_obs)
  gc(verbose = FALSE)
}

library_audit <- rbindlist(library_audit)
fwrite(library_audit, file.path(output_dir, "input_library_audit.tsv"), sep = "\t")
fwrite(crosswalk, file.path(output_dir, "archive_sample_crosswalk_used.tsv"), sep = "\t")
candidate_obs <- rbindlist(candidate_obs)
cat("Candidate/high-confidence T/NK cells", sum(library_audit$n_candidate_tnk),
    sum(library_audit$n_high_confidence_tnk), "\n")

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
group_id <- paste(hc_obs$specimen_id, hc_obs$lineage, sep = "::")
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
pseudobulk_meta <- data.table(group_id = groups, specimen_id = group_split[[1L]], lineage = group_split[[2L]])
sample_meta <- unique(crosswalk[, .(specimen_id = node_sample_id, patient_id, cancer, tissue_class)])
sample_idx <- match(pseudobulk_meta$specimen_id, sample_meta$specimen_id)
if (anyNA(sample_idx)) stop("pseudobulk specimens are missing from sample metadata")
pseudobulk_meta[, `:=`(
  patient_id = sample_meta$patient_id[sample_idx],
  cancer = sample_meta$cancer[sample_idx],
  tissue_class = sample_meta$tissue_class[sample_idx]
)]
cell_group_counts <- table(group_id)
pseudobulk_meta[, n_cells := as.integer(cell_group_counts[group_id])]
pseudobulk_meta[, total_umi := as.numeric(Matrix::colSums(pseudobulk_counts)[group_id])]
setcolorder(pseudobulk_meta, c("group_id", "specimen_id", "patient_id", "cancer", "tissue_class", "lineage", "n_cells", "total_umi"))
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
for (j in seq_len(nrow(pseudobulk_meta))) {
  lin <- pseudobulk_meta$lineage[[j]]
  wanted <- c("common_81", paste0("lineage_", lin))
  for (sig in wanted) {
    def <- signatures[signature == sig]
    up <- intersect(def[direction == "up", gene], rownames(rank_matrix))
    down <- intersect(def[direction == "down", gene], rownames(rank_matrix))
    audit <- signature_audit[signature == sig]
    minimum_coverage <- if (sig == "common_81") 0.80 else 0.70
    score_rows[[length(score_rows) + 1L]] <- cbind(
      pseudobulk_meta[j],
      data.table(
        signature = sig,
        score = mean(rank_matrix[up, j]) - mean(rank_matrix[down, j]),
        coverage = audit$coverage[[1L]],
        up_coverage = audit$up_coverage[[1L]],
        down_coverage = audit$down_coverage[[1L]],
        eligible_20 = pseudobulk_meta$n_cells[[j]] >= 20 & pseudobulk_meta$total_umi[[j]] >= 10000 &
          audit$coverage[[1L]] >= minimum_coverage & audit$up_coverage[[1L]] >= 0.70 &
          audit$down_coverage[[1L]] >= 0.70
      )
    )
  }
}
scores <- rbindlist(score_rows)
fwrite(scores, file.path(output_dir, "sample_lineage_scores.tsv"), sep = "\t")

paired_test <- function(dat, control_class, contrast_name) {
  wide <- dcast(dat[tissue_class %in% c("bone_metastasis", control_class)],
                patient_id + cancer ~ tissue_class, value.var = "score")
  if (!all(c("bone_metastasis", control_class) %in% names(wide))) {
    return(list(summary = data.table(contrast = contrast_name, n_pairs = 0L, n_positive = 0L,
                                    positive_fraction = NA_real_, median_difference = NA_real_,
                                    wilcoxon_one_sided_p = NA_real_, exact_test = NA),
                details = data.table()))
  }
  wide <- wide[!is.na(bone_metastasis) & !is.na(get(control_class))]
  wide[, `:=`(control_score = get(control_class), difference = bone_metastasis - get(control_class))]
  can_exact <- nrow(wide) > 0L && all(wide$difference != 0) && !anyDuplicated(abs(wide$difference))
  p <- if (nrow(wide)) suppressWarnings(wilcox.test(
    wide$difference, mu = 0, alternative = "greater", exact = can_exact
  )$p.value) else NA_real_
  list(
    summary = data.table(
      contrast = contrast_name,
      n_pairs = nrow(wide),
      n_positive = sum(wide$difference > 0),
      positive_fraction = if (nrow(wide)) mean(wide$difference > 0) else NA_real_,
      median_difference = if (nrow(wide)) median(wide$difference) else NA_real_,
      wilcoxon_one_sided_p = p,
      exact_test = can_exact
    ),
    details = wide[, .(contrast = contrast_name, patient_id, cancer,
                       bone_metastasis_score = bone_metastasis, control_score, difference)]
  )
}

test_rows <- list()
pair_rows <- list()
for (sig in unique(scores$signature)) {
  for (lin in unique(scores[signature == sig, lineage])) {
    dat <- scores[signature == sig & lineage == lin & eligible_20]
    for (contrast_spec in list(c("normal_bone", "bm_vs_normal_bone"), c("primary", "bm_vs_primary"))) {
      ans <- paired_test(dat, contrast_spec[[1L]], contrast_spec[[2L]])
      ans$summary[, `:=`(signature = sig, lineage = lin)]
      if (nrow(ans$details)) ans$details[, `:=`(signature = sig, lineage = lin)]
      test_rows[[length(test_rows) + 1L]] <- ans$summary
      pair_rows[[length(pair_rows) + 1L]] <- ans$details
    }
  }
}
tests <- rbindlist(test_rows, fill = TRUE)
pair_details <- rbindlist(pair_rows, fill = TRUE)
tests[, endpoint_class := ifelse(signature == "common_81", "common", "lineage_specific")]
tests[, FDR := p.adjust(wilcoxon_one_sided_p, method = "BH"), by = .(contrast, endpoint_class)]
setcolorder(tests, c("contrast", "signature", "lineage", "endpoint_class", "n_pairs", "n_positive",
                    "positive_fraction", "median_difference", "wilcoxon_one_sided_p", "FDR", "exact_test"))
fwrite(tests, file.path(output_dir, "paired_effects.tsv"), sep = "\t")
fwrite(pair_details, file.path(output_dir, "paired_differences.tsv"), sep = "\t")

sensitivity_rows <- list()
for (threshold in c(10L, 20L, 50L)) {
  for (lin in unique(scores$lineage)) {
    dat <- scores[
      signature == "common_81" & lineage == lin & n_cells >= threshold & total_umi >= 10000 &
        coverage >= 0.80 & up_coverage >= 0.70 & down_coverage >= 0.70
    ]
    for (contrast_spec in list(c("normal_bone", "bm_vs_normal_bone"), c("primary", "bm_vs_primary"))) {
      ans <- paired_test(dat, contrast_spec[[1L]], contrast_spec[[2L]])$summary
      ans[, `:=`(min_cells = threshold, lineage = lin)]
      sensitivity_rows[[length(sensitivity_rows) + 1L]] <- ans
    }
  }
}
sensitivity <- rbindlist(sensitivity_rows, fill = TRUE)
fwrite(sensitivity, file.path(output_dir, "cell_threshold_sensitivity.tsv"), sep = "\t")

plot_data <- pair_details[signature == "common_81"]
if (nrow(plot_data)) {
  plot_long <- melt(
    plot_data,
    id.vars = c("contrast", "patient_id", "cancer", "lineage"),
    measure.vars = c("control_score", "bone_metastasis_score"),
    variable.name = "tissue", value.name = "score"
  )
  plot_long[, tissue := factor(tissue, levels = c("control_score", "bone_metastasis_score"),
                               labels = c("Paired control", "Bone metastasis"))]
  p <- ggplot(plot_long, aes(x = tissue, y = score, group = patient_id, color = patient_id)) +
    geom_line(linewidth = 0.7, alpha = 0.8) +
    geom_point(size = 2.1) +
    facet_grid(contrast ~ lineage, scales = "free_y") +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom") +
    labs(x = NULL, y = "Frozen 81-gene directed rank score",
         title = "OEP005136 within-patient Gate8B validation")
  ggsave(file.path(output_dir, "common_score_paired_changes.pdf"), p, width = 9.5, height = 6.2)
}

primary_common <- tests[contrast == "bm_vs_normal_bone" & signature == "common_81"]
secondary_common <- tests[contrast == "bm_vs_primary" & signature == "common_81"]
primary_lineage <- tests[contrast == "bm_vs_normal_bone" & endpoint_class == "lineage_specific"]
primary_positive <- sum(primary_common$median_difference > 0, na.rm = TRUE)
primary_pair_support <- sum(primary_common$positive_fraction >= 0.75, na.rm = TRUE)
secondary_positive <- sum(secondary_common$median_difference > 0, na.rm = TRUE)
lineage_specific_positive <- sum(primary_lineage$median_difference > 0, na.rm = TRUE)

decision <- c(
  "# Gate8B OEP005136 decision",
  "",
  paste0("- Primary BM-versus-normal common-score positive lineages: ", primary_positive, "/3."),
  paste0("- Primary lineages with at least 75% positive evaluable pairs: ", primary_pair_support, "/3."),
  paste0("- Secondary BM-versus-primary common-score positive lineages: ", secondary_positive, "/3."),
  paste0("- Primary lineage-specific signatures with positive median paired effect: ", lineage_specific_positive, "/3."),
  "",
  "The primary normal-bone contrast contains only two patients; exact signed-rank P values are necessarily coarse. Direction and patient-level differences therefore carry more information than nominal significance.",
  "",
  if (primary_positive >= 2L && primary_pair_support >= 2L && lineage_specific_positive >= 2L) {
    "**Gate8B result: DIRECTIONALLY SUPPORTIVE.**"
  } else if (primary_positive >= 2L && primary_pair_support >= 2L) {
    "**Gate8B result: COMMON PROGRAM DIRECTIONALLY SUPPORTIVE, LINEAGE PROGRAM INCOMPLETE.**"
  } else {
    "**Gate8B result: DOES NOT SUPPORT THE FROZEN COMMON PROGRAM.**"
  },
  "",
  "Gate8B is interpreted jointly with the already frozen Gate8A result; no Gate8B expression result modifies the signatures or annotation markers."
)
writeLines(decision, file.path(output_dir, "gate8b_decision.md"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
cat("Gate8B completed", format(Sys.time()), "\n")
