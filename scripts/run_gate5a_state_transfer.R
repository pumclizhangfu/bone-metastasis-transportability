#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(SingleCellExperiment)
  library(scuttle)
  library(SingleR)
  library(BiocParallel)
})

parse_args <- function(x) {
  if (length(x) %% 2 != 0) stop("Arguments must be provided as --key value pairs")
  out <- list()
  for (i in seq(1, length(x), by = 2)) out[[sub("^--", "", x[[i]])]] <- x[[i + 1]]
  out
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("input", "gate4-genes", "work", "output", "workers", "seed", "min-cells")
missing_args <- setdiff(required, names(cfg))
if (length(missing_args)) stop("Missing arguments: ", paste(missing_args, collapse = ", "))

input_dir <- normalizePath(cfg$input, mustWork = TRUE)
gate4_gene_file <- normalizePath(cfg[["gate4-genes"]], mustWork = TRUE)
work_dir <- normalizePath(cfg$work, mustWork = FALSE)
output_dir <- normalizePath(cfg$output, mustWork = FALSE)
workers <- as.integer(cfg$workers)
seed <- as.integer(cfg$seed)
min_cells <- as.integer(cfg[["min-cells"]])
if (is.na(workers) || workers < 1 || workers > 20) stop("workers must be 1..20")
if (is.na(seed)) stop("seed must be an integer")
if (is.na(min_cells) || min_cells < 10) stop("min-cells must be at least 10")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(seed)
bp <- if (.Platform$OS.type == "unix" && workers > 1) {
  MulticoreParam(workers = workers, progressbar = FALSE)
} else {
  SerialParam(progressbar = FALSE)
}

message("CONFIG\tinput=", input_dir, "\tgate4_genes=", gate4_gene_file,
        "\twork=", work_dir, "\toutput=", output_dir,
        "\tworkers=", workers, "\tseed=", seed, "\tmin_cells=", min_cells)

lineages <- c("Myeloid", "T_NK")
accessions <- c("GSE143791", "GSE202813")
state_definition <- data.table(
  author_label = c(
    "Mono-1", "Mono-2", "Mono-3", "Monocyte pro", "Macrophage", "mDC", "PDC",
    "Thelper", "CD4 Naive", "Treg", "CTL-1", "CTL-2", "CTL-3", "CD8 Naive",
    "NKT", "NK1", "NK2"
  ),
  lineage = c(
    rep("Myeloid", 7), rep("T_NK", 10)
  ),
  meta_state = c(
    rep("Monocyte", 4), "Macrophage", rep("DC", 2),
    rep("CD4_T", 3), rep("CD8_CTL", 4), rep("NK_NKT", 3)
  ),
  primary_state = TRUE
)
state_levels <- list(
  Myeloid = c("Monocyte", "Macrophage", "DC"),
  T_NK = c("CD4_T", "CD8_CTL", "NK_NKT")
)
marker_modules <- list(
  Monocyte = c("FCN1", "S100A8", "S100A9", "CTSS", "LYZ", "VCAN"),
  Macrophage = c("C1QA", "C1QB", "C1QC", "APOE", "MSR1", "MRC1"),
  DC = c("FCER1A", "CD1C", "CLEC10A", "GZMB", "IRF7", "TCF4"),
  CD4_T = c("IL7R", "LTB", "MAL", "CCR7", "TCF7", "NOSIP"),
  CD8_CTL = c("CD8A", "CD8B", "CCL5", "GZMK", "CTSW", "NKG7"),
  NK_NKT = c("KLRD1", "FCGR3A", "GNLY", "PRF1", "TYROBP", "TRDC")
)
fwrite(state_definition, file.path(output_dir, "state_definition.tsv"), sep = "\t")
fwrite(rbindlist(lapply(names(marker_modules), function(x) {
  data.table(meta_state = x, gene = marker_modules[[x]])
})), file.path(output_dir, "state_marker_modules.tsv"), sep = "\t")

sce_files <- sort(list.files(input_dir, pattern = "\\.rds$", full.names = TRUE))
if (length(sce_files) != 42) stop("Expected 42 annotated SCE files; found ", length(sce_files))
sce_list <- lapply(sce_files, readRDS)
names(sce_list) <- sub("\\.rds$", "", basename(sce_files))
if (!all(vapply(sce_list, function(x) inherits(x, "SingleCellExperiment"), logical(1)))) {
  stop("Every input must be a SingleCellExperiment")
}
all_genes <- rownames(sce_list[[1]])
if (!all(vapply(sce_list, function(x) identical(rownames(x), all_genes), logical(1)))) {
  stop("Gene rows are not identical across annotated SCE files")
}

sample_design <- rbindlist(lapply(sce_list, function(sce) data.table(
  accession = unique(as.character(sce$accession)),
  cancer = unique(as.character(sce$cancer)),
  sample_id = unique(as.character(sce$sample_id)),
  patient_id = unique(as.character(sce$patient_id)),
  compartment = unique(as.character(sce$compartment))
)), use.names = TRUE)
if (anyDuplicated(sample_design$sample_id)) stop("Sample IDs are not unique")
complete_patients <- sample_design[, .(n_compartments = uniqueN(compartment)),
                                   by = .(accession, cancer, patient_id)][n_compartments == 3]

gate4_genes <- fread(gate4_gene_file)
if (!all(c("lineage", "gene") %in% names(gate4_genes))) stop("Gate4 robust-gene table lacks lineage/gene")
technical_genes <- all_genes[grepl("^MT-|^RPS|^RPL", all_genes)]

pool_cache <- new.env(parent = emptyenv())
make_pool <- function(accession_value, lineage_value) {
  key <- paste(accession_value, lineage_value, sep = "::")
  if (exists(key, envir = pool_cache, inherits = FALSE)) return(get(key, envir = pool_cache))
  mats <- list()
  metas <- list()
  counter <- 0L
  for (sample_id in names(sce_list)) {
    sce <- sce_list[[sample_id]]
    if (unique(as.character(sce$accession)) != accession_value) next
    idx <- which(as.character(sce$broad_class) == lineage_value)
    if (!length(idx)) next
    counter <- counter + 1L
    mats[[counter]] <- counts(sce)[, idx, drop = FALSE]
    metas[[counter]] <- data.table(
      accession = accession_value,
      cancer = unique(as.character(sce$cancer)),
      sample_id = unique(as.character(sce$sample_id)),
      patient_id = unique(as.character(sce$patient_id)),
      compartment = unique(as.character(sce$compartment)),
      barcode = colnames(sce)[idx],
      author_label = as.character(sce$author_label[idx]),
      broad_class = lineage_value
    )
  }
  if (!length(mats)) stop("No cells for ", accession_value, " ", lineage_value)
  mat <- as(do.call(cbind, mats), "dgCMatrix")
  meta <- rbindlist(metas, use.names = TRUE)
  if (!identical(colnames(mat), meta$barcode)) stop("Pool columns and metadata are not aligned")
  out <- list(counts = mat, meta = meta)
  assign(key, out, envir = pool_cache)
  out
}

collapse_truth <- function(author_label, lineage_value) {
  lookup <- state_definition[lineage == lineage_value]
  unname(setNames(lookup$meta_state, lookup$author_label)[author_label])
}

aggregate_reference <- function(count_matrix, metadata, labels, allowed_genes) {
  group <- paste(metadata$sample_id, labels, sep = "::")
  group_n <- table(group)
  keep_group <- names(group_n)[group_n >= min_cells]
  keep <- group %in% keep_group
  if (!any(keep)) stop("No reference groups meet min_cells")
  group <- factor(group[keep], levels = keep_group)
  design <- sparseMatrix(i = seq_along(group), j = as.integer(group), x = 1,
                         dims = c(length(group), nlevels(group)),
                         dimnames = list(NULL, levels(group)))
  agg <- count_matrix[allowed_genes, keep, drop = FALSE] %*% design
  info <- unique(data.table(
    reference_id = as.character(group),
    sample_id = metadata$sample_id[keep],
    patient_id = metadata$patient_id[keep],
    meta_state = labels[keep]
  ), by = "reference_id")
  info <- info[match(colnames(agg), reference_id)]
  if (anyNA(info$meta_state)) stop("Reference metadata alignment failed")
  state_reps <- info[, .N, by = meta_state]
  if (any(state_reps$N < 2)) stop("A reference state has fewer than two sample replicates")
  ref <- SingleCellExperiment(assays = list(counts = as(agg, "dgCMatrix")),
                              colData = S4Vectors::DataFrame(info))
  logNormCounts(ref)
}

classify_cells <- function(count_matrix, ref, run_seed) {
  test <- SingleCellExperiment(assays = list(counts = count_matrix))
  test <- logNormCounts(test)
  set.seed(run_seed)
  trained <- trainSingleR(
    ref = ref,
    labels = ref$meta_state,
    de.method = "classic",
    assay.type = "logcounts",
    BPPARAM = bp
  )
  pred <- classifySingleR(
    test = test,
    trained = trained,
    assay.type = "logcounts",
    fine.tune = TRUE,
    prune = TRUE,
    BPPARAM = bp
  )
  data.table(
    raw_state = as.character(pred$labels),
    pruned_state = as.character(pred$pruned.labels),
    delta_next = as.numeric(pred$delta.next)
  )
}

cv_rows <- list()
classifier_gene_rows <- list()
cv_counter <- 0L
for (lineage_value in lineages) {
  rcc <- make_pool("GSE202813", lineage_value)
  truth <- collapse_truth(rcc$meta$author_label, lineage_value)
  keep_truth <- !is.na(truth)
  rcc_counts <- rcc$counts[, keep_truth, drop = FALSE]
  rcc_meta <- rcc$meta[keep_truth]
  truth <- truth[keep_truth]
  excluded <- unique(c(technical_genes, gate4_genes[lineage == lineage_value, gene]))
  allowed <- setdiff(all_genes, excluded)
  if (length(allowed) < 5000) stop("Fewer than 5,000 classifier genes remain for ", lineage_value)
  classifier_gene_rows[[lineage_value]] <- data.table(
    lineage = lineage_value,
    gene = all_genes,
    excluded = all_genes %in% excluded,
    exclusion_reason = fifelse(all_genes %in% technical_genes, "technical",
                               fifelse(all_genes %in% gate4_genes[lineage == lineage_value, gene],
                                       "Gate4b_robust_gradient", "retained"))
  )
  patients <- sort(unique(rcc_meta$patient_id))
  for (left_out in patients) {
    cv_counter <- cv_counter + 1L
    train_idx <- rcc_meta$patient_id != left_out
    test_idx <- !train_idx
    ref <- aggregate_reference(rcc_counts[, train_idx, drop = FALSE], rcc_meta[train_idx],
                               truth[train_idx], allowed)
    pred <- classify_cells(rcc_counts[allowed, test_idx, drop = FALSE], ref,
                           seed + 100L + cv_counter)
    pred[, `:=`(
      lineage = lineage_value,
      left_out_patient = left_out,
      barcode = rcc_meta$barcode[test_idx],
      sample_id = rcc_meta$sample_id[test_idx],
      patient_id = rcc_meta$patient_id[test_idx],
      compartment = rcc_meta$compartment[test_idx],
      true_state = truth[test_idx]
    )]
    cv_rows[[paste(lineage_value, left_out, sep = "::")]] <- pred
    message("CV_FOLD_COMPLETE\t", lineage_value, "\tleft_out=", left_out,
            "\tcells=", nrow(pred), "\tcoverage=", sprintf("%.4f", mean(!is.na(pred$pruned_state))))
  }
}
classifier_genes <- rbindlist(classifier_gene_rows, use.names = TRUE)
fwrite(classifier_genes, file.path(output_dir, "classifier_gene_audit.tsv.gz"),
       sep = "\t", compress = "gzip")

cv_predictions <- rbindlist(cv_rows, use.names = TRUE)
cv_predictions[is.na(pruned_state) | pruned_state == "", pruned_state := "Unassigned_state"]
setcolorder(cv_predictions, c("lineage", "left_out_patient", "sample_id", "patient_id",
                              "compartment", "barcode", "true_state", "raw_state",
                              "pruned_state", "delta_next"))
fwrite(cv_predictions, file.path(output_dir, "cv_predictions.tsv.gz"), sep = "\t", compress = "gzip")

cv_confusion <- cv_predictions[, .N, by = .(lineage, true_state, pruned_state)]
setorder(cv_confusion, lineage, true_state, -N)
fwrite(cv_confusion, file.path(output_dir, "cv_confusion.tsv"), sep = "\t")
cv_state_metrics <- cv_predictions[, .(
  n_cells = .N,
  recall_all = mean(pruned_state == true_state),
  assigned_coverage = mean(pruned_state != "Unassigned_state")
), by = .(lineage, true_state)]
cv_state_metrics[, precision_assigned := vapply(seq_len(.N), function(i) {
  lin <- lineage[[i]]
  state <- true_state[[i]]
  tp <- cv_predictions[lineage == lin & true_state == state & pruned_state == state, .N]
  denom <- cv_predictions[lineage == lin & pruned_state == state, .N]
  if (denom > 0) tp / denom else NA_real_
}, numeric(1))]
fwrite(cv_state_metrics, file.path(output_dir, "cv_state_metrics.tsv"), sep = "\t")

cv_lineage_summary <- cv_predictions[, .(
  n_cells = .N,
  assigned_coverage = mean(pruned_state != "Unassigned_state"),
  accuracy_all = mean(pruned_state == true_state)
), by = lineage]
cv_lineage_summary[, balanced_accuracy := vapply(lineage, function(x) {
  mean(cv_state_metrics[lineage == x, recall_all])
}, numeric(1))]
cv_lineage_summary[, min_state_recall := vapply(lineage, function(x) {
  min(cv_state_metrics[lineage == x, recall_all])
}, numeric(1))]
cv_lineage_summary[, cv_pass := assigned_coverage >= 0.85 & balanced_accuracy >= 0.70 &
                       min_state_recall >= 0.60]
fwrite(cv_lineage_summary, file.path(output_dir, "cv_lineage_summary.tsv"), sep = "\t")

transfer_rows <- list()
pool_state_vectors <- list()
transfer_counter <- 0L
for (lineage_value in lineages) {
  rcc <- make_pool("GSE202813", lineage_value)
  rcc_truth <- collapse_truth(rcc$meta$author_label, lineage_value)
  keep_truth <- !is.na(rcc_truth)
  excluded <- unique(c(technical_genes, gate4_genes[lineage == lineage_value, gene]))
  allowed <- setdiff(all_genes, excluded)
  ref <- aggregate_reference(rcc$counts[, keep_truth, drop = FALSE], rcc$meta[keep_truth],
                             rcc_truth[keep_truth], allowed)

  prostate <- make_pool("GSE143791", lineage_value)
  transfer_counter <- transfer_counter + 1L
  pred <- classify_cells(prostate$counts[allowed, , drop = FALSE], ref,
                         seed + 1000L + transfer_counter)
  pred[is.na(pruned_state) | pruned_state == "", pruned_state := "Unassigned_state"]
  prostate_rows <- cbind(prostate$meta[, .(accession, cancer, sample_id, patient_id,
                                           compartment, barcode, broad_class)], pred)
  prostate_rows[, `:=`(lineage = lineage_value, meta_state = pruned_state,
                       state_source = "SingleR_RCC_meta_state")]
  transfer_rows[[paste(lineage_value, "prostate", sep = "::")]] <- prostate_rows
  pool_state_vectors[[paste("GSE143791", lineage_value, sep = "::")]] <- prostate_rows$meta_state

  full_rcc_state <- rep("Excluded_nonprimary_state", ncol(rcc$counts))
  full_rcc_state[keep_truth] <- rcc_truth[keep_truth]
  rcc_rows <- rcc$meta[, .(accession, cancer, sample_id, patient_id,
                            compartment, barcode, broad_class)]
  rcc_rows[, `:=`(
    raw_state = full_rcc_state,
    pruned_state = full_rcc_state,
    delta_next = NA_real_,
    lineage = lineage_value,
    meta_state = full_rcc_state,
    state_source = "RCC_author_collapsed"
  )]
  transfer_rows[[paste(lineage_value, "renal", sep = "::")]] <- rcc_rows
  pool_state_vectors[[paste("GSE202813", lineage_value, sep = "::")]] <- full_rcc_state
  message("TRANSFER_COMPLETE\t", lineage_value, "\tprostate_cells=", nrow(prostate_rows),
          "\tcoverage=", sprintf("%.4f", mean(prostate_rows$meta_state != "Unassigned_state")))
}

state_annotations <- rbindlist(transfer_rows, use.names = TRUE, fill = TRUE)
setcolorder(state_annotations, c("accession", "cancer", "sample_id", "patient_id",
                                 "compartment", "barcode", "broad_class", "lineage",
                                 "meta_state", "state_source", "raw_state", "pruned_state",
                                 "delta_next"))
if (anyDuplicated(state_annotations$barcode)) stop("State annotation barcodes are duplicated")
fwrite(state_annotations, file.path(output_dir, "state_annotations.tsv.gz"),
       sep = "\t", compress = "gzip", na = "NA")

transfer_coverage <- state_annotations[, .(
  lineage_cells = .N,
  assigned_cells = sum(!meta_state %in% c("Unassigned_state", "Excluded_nonprimary_state")),
  assigned_coverage = mean(!meta_state %in% c("Unassigned_state", "Excluded_nonprimary_state"))
), by = .(accession, cancer, lineage, state_source)]
fwrite(transfer_coverage, file.path(output_dir, "transfer_coverage.tsv"), sep = "\t")

valid_annotations <- state_annotations[!meta_state %in% c("Unassigned_state", "Excluded_nonprimary_state")]
state_counts <- valid_annotations[, .(n_cells = .N),
                                  by = .(accession, cancer, sample_id, patient_id,
                                         compartment, lineage, meta_state)]
lineage_totals <- state_annotations[, .(lineage_total_cells = .N),
                                    by = .(accession, sample_id, lineage)]
assigned_totals <- valid_annotations[, .(assigned_lineage_cells = .N),
                                     by = .(accession, sample_id, lineage)]
state_counts <- merge(state_counts, lineage_totals, by = c("accession", "sample_id", "lineage"), all.x = TRUE)
state_counts <- merge(state_counts, assigned_totals, by = c("accession", "sample_id", "lineage"), all.x = TRUE)
state_counts[, `:=`(
  fraction_of_all_lineage = n_cells / lineage_total_cells,
  fraction_of_assigned_lineage = n_cells / assigned_lineage_cells
)]
setorder(state_counts, accession, lineage, meta_state, patient_id, compartment)
fwrite(state_counts, file.path(output_dir, "state_counts.tsv"), sep = "\t")

eligibility_rows <- list()
for (lineage_value in lineages) {
  for (state_value in state_levels[[lineage_value]]) {
    for (accession_value in accessions) {
      design <- merge(
        complete_patients[accession == accession_value, .(accession, cancer, patient_id)],
        sample_design[accession == accession_value, .(accession, sample_id, patient_id, compartment)],
        by = c("accession", "patient_id"), all = FALSE
      )
      z <- state_counts[accession == accession_value & lineage == lineage_value & meta_state == state_value,
                        .(sample_id, n_cells)]
      design <- merge(design, z, by = "sample_id", all.x = TRUE)
      design[is.na(n_cells), n_cells := 0L]
      patient_ok <- design[, .(
        min_cells_across_triplet = min(n_cells),
        complete_eligible_triplet = .N == 3 & all(n_cells >= min_cells)
      ), by = .(accession, cancer, patient_id)]
      eligibility_rows[[paste(accession_value, lineage_value, state_value, sep = "::")]] <- data.table(
        accession = accession_value,
        cancer = unique(design$cancer),
        lineage = lineage_value,
        meta_state = state_value,
        complete_eligible_patients = sum(patient_ok$complete_eligible_triplet),
        total_complete_patients = nrow(patient_ok),
        min_cells_among_eligible = if (any(patient_ok$complete_eligible_triplet)) {
          min(patient_ok$min_cells_across_triplet[patient_ok$complete_eligible_triplet])
        } else NA_real_
      )
    }
  }
}
state_eligibility <- rbindlist(eligibility_rows, use.names = TRUE)
cross_ready <- state_eligibility[, .(
  cross_cancer_ready = all(complete_eligible_patients >= 3),
  min_complete_patients = min(complete_eligible_patients)
), by = .(lineage, meta_state)]
state_eligibility <- merge(state_eligibility, cross_ready, by = c("lineage", "meta_state"), all.x = TRUE)
setorder(state_eligibility, lineage, meta_state, accession)
fwrite(state_eligibility, file.path(output_dir, "state_eligibility.tsv"), sep = "\t")

# Create exact raw-count state pseudobulks for Gate5b.
pb_columns <- list()
pb_meta_rows <- list()
for (accession_value in accessions) {
  for (lineage_value in lineages) {
    pool <- make_pool(accession_value, lineage_value)
    states <- pool_state_vectors[[paste(accession_value, lineage_value, sep = "::")]]
    keep <- !states %in% c("Unassigned_state", "Excluded_nonprimary_state")
    groups <- unique(data.table(sample_id = pool$meta$sample_id[keep], meta_state = states[keep]))
    for (i in seq_len(nrow(groups))) {
      sample_value <- groups$sample_id[[i]]
      state_value <- groups$meta_state[[i]]
      idx <- which(keep & pool$meta$sample_id == sample_value & states == state_value)
      pb_id <- paste(accession_value, sample_value, lineage_value, state_value, sep = "::")
      pb_vec <- Matrix::rowSums(pool$counts[, idx, drop = FALSE])
      pb_columns[[pb_id]] <- pb_vec
      meta_row <- pool$meta[idx[[1]]]
      pb_meta_rows[[pb_id]] <- data.table(
        pseudobulk_id = pb_id,
        accession = accession_value,
        cancer = meta_row$cancer,
        sample_id = sample_value,
        patient_id = meta_row$patient_id,
        compartment = meta_row$compartment,
        lineage = lineage_value,
        meta_state = state_value,
        n_cells = length(idx),
        raw_umi_sum = as.numeric(sum(pool$counts[, idx, drop = FALSE])),
        aggregation_exact = isTRUE(all.equal(as.numeric(sum(pb_vec)),
                                             as.numeric(sum(pool$counts[, idx, drop = FALSE])),
                                             tolerance = 0))
      )
    }
  }
}
state_pb_counts <- as(do.call(cbind, pb_columns), "dgCMatrix")
state_pb_meta <- rbindlist(pb_meta_rows, use.names = TRUE)
state_pb_meta <- state_pb_meta[match(colnames(state_pb_counts), pseudobulk_id)]
if (!all(state_pb_meta$aggregation_exact)) stop("State pseudobulk aggregation failed")
saveRDS(list(
  counts = state_pb_counts,
  metadata = as.data.frame(state_pb_meta),
  state_definition = as.data.frame(state_definition),
  classifier_excluded_gate4_genes = TRUE,
  seed = seed,
  min_cells = min_cells
), file.path(output_dir, "state_pseudobulk_counts.rds"), compress = "xz")
fwrite(state_pb_meta, file.path(output_dir, "state_pseudobulk_metadata.tsv"), sep = "\t")

# Marker-module audit on aggregated state pseudobulks.
marker_rows <- list()
marker_counter <- 0L
for (accession_value in accessions) {
  for (lineage_value in lineages) {
    states_here <- state_levels[[lineage_value]]
    meta_idx <- which(state_pb_meta$accession == accession_value &
                        state_pb_meta$lineage == lineage_value &
                        state_pb_meta$meta_state %in% states_here)
    for (state_value in states_here) {
      cols <- meta_idx[state_pb_meta$meta_state[meta_idx] == state_value]
      if (!length(cols)) next
      summed <- Matrix::rowSums(state_pb_counts[, cols, drop = FALSE])
      log_cpm <- log2(1 + 1e6 * summed / sum(summed))
      names(log_cpm) <- rownames(state_pb_counts)
      for (module_value in states_here) {
        marker_counter <- marker_counter + 1L
        genes <- intersect(marker_modules[[module_value]], names(log_cpm))
        marker_rows[[marker_counter]] <- data.table(
          accession = accession_value,
          lineage = lineage_value,
          assigned_state = state_value,
          marker_module = module_value,
          n_marker_genes = length(genes),
          mean_log2_cpm = mean(log_cpm[genes])
        )
      }
    }
  }
}
marker_scores <- rbindlist(marker_rows, use.names = TRUE)
marker_scores[, expected_module := assigned_state == marker_module]
marker_scores[, module_rank := frank(-mean_log2_cpm, ties.method = "min"),
              by = .(accession, lineage, assigned_state)]
fwrite(marker_scores, file.path(output_dir, "state_marker_audit.tsv"), sep = "\t")

prostate_coverage <- transfer_coverage[accession == "GSE143791",
                                       .(lineage, prostate_coverage = assigned_coverage)]
gate_summary <- merge(cv_lineage_summary, prostate_coverage, by = "lineage", all = TRUE)
ready_counts <- cross_ready[cross_cancer_ready == TRUE, .(cross_ready_states = .N), by = lineage]
gate_summary <- merge(gate_summary, ready_counts, by = "lineage", all.x = TRUE)
gate_summary[is.na(cross_ready_states), cross_ready_states := 0L]
gate_summary[, transfer_pass := prostate_coverage >= 0.85]
gate_summary[, lineage_ready := cv_pass & transfer_pass & cross_ready_states >= 1]
fwrite(gate_summary, file.path(output_dir, "gate5a_summary.tsv"), sep = "\t")

myeloid_ready <- gate_summary[lineage == "Myeloid", lineage_ready]
tnk_ready <- gate_summary[lineage == "T_NK", lineage_ready]
myeloid_states <- gate_summary[lineage == "Myeloid", cross_ready_states]
tnk_states <- gate_summary[lineage == "T_NK", cross_ready_states]
decision <- if (isTRUE(myeloid_ready) && isTRUE(tnk_ready) &&
                myeloid_states >= 1 && tnk_states >= 2) {
  "GO"
} else if ((isTRUE(myeloid_ready) && myeloid_states >= 1) ||
           (isTRUE(tnk_ready) && tnk_states >= 1)) {
  "CONDITIONAL GO"
} else {
  "NO-GO"
}

decision_lines <- c(
  "# Gate5a state-transfer decision",
  "",
  paste0("**", decision, "**"),
  "",
  paste0("- Myeloid CV balanced accuracy: ", sprintf("%.3f", gate_summary[lineage == "Myeloid", balanced_accuracy])),
  paste0("- Myeloid minimum state recall: ", sprintf("%.3f", gate_summary[lineage == "Myeloid", min_state_recall])),
  paste0("- Myeloid prostate transfer coverage: ", sprintf("%.3f", gate_summary[lineage == "Myeloid", prostate_coverage])),
  paste0("- Myeloid cross-cancer-ready states: ", myeloid_states),
  paste0("- T/NK CV balanced accuracy: ", sprintf("%.3f", gate_summary[lineage == "T_NK", balanced_accuracy])),
  paste0("- T/NK minimum state recall: ", sprintf("%.3f", gate_summary[lineage == "T_NK", min_state_recall])),
  paste0("- T/NK prostate transfer coverage: ", sprintf("%.3f", gate_summary[lineage == "T_NK", prostate_coverage])),
  paste0("- T/NK cross-cancer-ready states: ", tnk_states),
  "",
  "Classifier features excluded Gate4b robust-gradient genes and technical mitochondrial/ribosomal genes.",
  "Gate5a validates state transfer and eligibility only; it does not test composition or within-state gradients."
)
writeLines(decision_lines, file.path(output_dir, "gate5a_decision.md"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("GATE5A_COMPLETE\tdecision=", decision,
        "\tmyeloid_bacc=", sprintf("%.4f", gate_summary[lineage == "Myeloid", balanced_accuracy]),
        "\ttnk_bacc=", sprintf("%.4f", gate_summary[lineage == "T_NK", balanced_accuracy]),
        "\tmyeloid_ready_states=", myeloid_states,
        "\ttnk_ready_states=", tnk_states)
