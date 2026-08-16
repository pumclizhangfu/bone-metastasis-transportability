#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(data.table)
  library(SingleR)
  library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: build_gate12g_transfer_reference.R <annotated_sce_dir> <gate12b_state_map.tsv.gz> <outdir>")
}

sce_dir <- args[[1L]]
state_file <- args[[2L]]
outdir <- args[[3L]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

seed <- 1208L
set.seed(seed)
min_cells_pb <- 20L
top_markers_per_label <- 60L

broad_train <- c("T_NK", "Myeloid", "B", "Progenitor", "Stromal", "Endothelial",
                 "Osteoclast", "Malignant", "Erythroid")
myeloid_train <- c("Classical_monocyte", "Inflammatory_monocyte", "C1QC_macrophage",
                   "Resident_macrophage", "cDC", "pDC", "Proliferating_myeloid")
tnk_train <- c("CD4_naive", "CD4_memory", "Treg", "CD8_effector", "CD8_exhausted",
               "NK_adaptive", "Proliferating_T_NK")

files <- sort(list.files(sce_dir, pattern = "\\.rds$", full.names = TRUE))
if (length(files) != 42L) stop("Expected 42 discovery SCE files; found ", length(files))
if (!file.exists(state_file)) stop("Missing Gate12B state map: ", state_file)

state_map <- fread(state_file, select = c("barcode", "sample_id", "patient_id", "lineage", "gate12b_state"))
if (anyDuplicated(state_map$barcode)) stop("Gate12B barcode is not unique")
setkey(state_map, barcode)

pb_store <- list(Broad = list(counts = list(), meta = list()),
                 Myeloid = list(counts = list(), meta = list()),
                 T_NK = list(counts = list(), meta = list()))
gene_order <- NULL

append_groups <- function(counts, labels, sample_id, patient_id, taxonomy, allowed) {
  labels <- as.character(labels)
  for (lab in allowed) {
    idx <- which(labels == lab)
    if (length(idx) < min_cells_pb) next
    key <- paste(sample_id, lab, sep = "::")
    pb_store[[taxonomy]]$counts[[key]] <<- Matrix::rowSums(counts[, idx, drop = FALSE])
    pb_store[[taxonomy]]$meta[[key]] <<- data.table(
      pseudobulk_id = key, sample_id = sample_id, patient_id = patient_id,
      taxonomy = taxonomy, label = lab, cells = length(idx)
    )
  }
}

for (i in seq_along(files)) {
  sce <- readRDS(files[[i]])
  if (is.null(gene_order)) gene_order <- rownames(sce)
  if (!identical(rownames(sce), gene_order)) stop("Gene order mismatch: ", files[[i]])
  counts <- as(assay(sce, "counts"), "dgCMatrix")
  md <- as.data.table(as.data.frame(colData(sce)))
  sample_id <- unique(md$sample_id)
  patient_id <- unique(md$patient_id)
  if (length(sample_id) != 1L || length(patient_id) != 1L) stop("Non-unique sample/patient in ", files[[i]])

  append_groups(counts, md$broad_class, sample_id, patient_id, "Broad", broad_train)

  sm <- state_map[J(colnames(sce))]
  if (!identical(sm$barcode, colnames(sce))) stop("State-map join failed for ", sample_id)
  append_groups(counts, fifelse(sm$lineage == "Myeloid", sm$gate12b_state, NA_character_),
                sample_id, patient_id, "Myeloid", myeloid_train)
  append_groups(counts, fifelse(sm$lineage == "T_NK", sm$gate12b_state, NA_character_),
                sample_id, patient_id, "T_NK", tnk_train)
  rm(sce, counts, md, sm); invisible(gc(FALSE))
}

finalize_pb <- function(store, taxonomy) {
  if (!length(store$counts)) stop("No pseudobulks for ", taxonomy)
  count_mat <- do.call(cbind, store$counts)
  rownames(count_mat) <- gene_order
  meta <- rbindlist(store$meta)
  meta <- meta[match(colnames(count_mat), pseudobulk_id)]
  if (!identical(meta$pseudobulk_id, colnames(count_mat))) stop("Pseudobulk metadata mismatch")
  lib <- colSums(count_mat)
  if (any(lib <= 0)) stop("Zero pseudobulk library")
  logcpm <- log1p(t(t(count_mat) / lib) * 1e6)
  list(logcpm = logcpm, meta = meta)
}

select_markers <- function(logcpm, labels, top_n = 60L) {
  labels <- factor(labels)
  labs <- levels(labels)
  means <- sapply(labs, function(lab) rowMeans(logcpm[, labels == lab, drop = FALSE]))
  if (length(labs) == 2L) means <- matrix(means, ncol = 2L, dimnames = list(rownames(logcpm), labs))
  marker_list <- lapply(seq_along(labs), function(j) {
    others <- apply(means[, -j, drop = FALSE], 1L, max)
    effect <- means[, j] - others
    ok <- is.finite(effect) & effect > 0 & means[, j] > log1p(1)
    genes <- rownames(means)[ok]
    genes[order(effect[ok], decreasing = TRUE)][seq_len(min(top_n, sum(ok)))]
  })
  names(marker_list) <- labs
  unique(unlist(marker_list, use.names = FALSE))
}

classify_matrix <- function(test, ref, ref_labels) {
  pred <- SingleR(test = test, ref = ref, labels = ref_labels,
                  de.method = "classic", fine.tune = TRUE, prune = TRUE,
                  BPPARAM = SerialParam())
  data.table(predicted = as.character(pred$labels),
             pruned = as.character(pred$pruned.labels),
             delta_next = as.numeric(pred$delta.next))
}

calibrate_taxonomy <- function(pb, taxonomy, allowed) {
  meta <- pb$meta
  patients <- unique(meta$patient_id)
  rows <- list()
  for (pid in patients) {
    train <- which(meta$patient_id != pid)
    test <- which(meta$patient_id == pid)
    if (!length(test) || uniqueN(meta$label[train]) < 2L) next
    fold_markers <- select_markers(pb$logcpm[, train, drop = FALSE], meta$label[train], top_markers_per_label)
    if (length(fold_markers) < 20L) stop("Too few fold markers for ", taxonomy, " patient ", pid)
    pr <- classify_matrix(pb$logcpm[fold_markers, test, drop = FALSE],
                          pb$logcpm[fold_markers, train, drop = FALSE], meta$label[train])
    pr <- cbind(meta[test, .(pseudobulk_id, sample_id, patient_id, truth = label, cells)], pr)
    pr[, taxonomy := taxonomy]
    rows[[length(rows) + 1L]] <- pr
  }
  predictions <- rbindlist(rows)
  predictions[, final := fifelse(is.na(pruned) | !nzchar(pruned), "REJECTED", pruned)]
  class_metrics <- predictions[, .(
    evaluable_patients = uniqueN(patient_id),
    pseudobulks = .N,
    recall = mean(final == truth),
    rejection_fraction = mean(final == "REJECTED")
  ), by = truth]
  class_metrics[, class_eligible := evaluable_patients >= 3L]
  eligible <- class_metrics[class_eligible == TRUE]
  summary <- data.table(
    taxonomy = taxonomy,
    classes_total = nrow(class_metrics),
    classes_evaluable = nrow(eligible),
    balanced_accuracy = if (nrow(eligible)) mean(eligible$recall) else NA_real_,
    median_class_recall = if (nrow(eligible)) median(eligible$recall) else NA_real_,
    threshold_balanced_accuracy = if (taxonomy == "Broad") 0.70 else 0.60,
    threshold_median_recall = 0.50
  )
  summary[, pass := balanced_accuracy >= threshold_balanced_accuracy &
                   median_class_recall >= threshold_median_recall]
  list(predictions = predictions, class_metrics = class_metrics, summary = summary)
}

references <- list()
calibration <- list()
allowed_map <- list(Broad = broad_train, Myeloid = myeloid_train, T_NK = tnk_train)
for (taxonomy in names(pb_store)) {
  message("Finalizing and calibrating ", taxonomy)
  pb <- finalize_pb(pb_store[[taxonomy]], taxonomy)
  cal <- calibrate_taxonomy(pb, taxonomy, allowed_map[[taxonomy]])
  final_markers <- select_markers(pb$logcpm, pb$meta$label, top_markers_per_label)
  references[[taxonomy]] <- list(
    logcpm = pb$logcpm[final_markers, , drop = FALSE],
    meta = pb$meta,
    marker_genes = final_markers,
    allowed_labels = allowed_map[[taxonomy]]
  )
  calibration[[taxonomy]] <- cal
  rm(pb); invisible(gc(FALSE))
}

predictions <- rbindlist(lapply(calibration, `[[`, "predictions"), fill = TRUE)
class_metrics <- rbindlist(lapply(names(calibration), function(nm) {
  z <- copy(calibration[[nm]]$class_metrics)
  z[, taxonomy := nm]
  z
}), fill = TRUE)
summary <- rbindlist(lapply(calibration, `[[`, "summary"), fill = TRUE)
overall_status <- if (all(summary$pass)) "PASS" else "STOP"

reference_manifest <- rbindlist(lapply(names(references), function(nm) {
  z <- references[[nm]]
  data.table(taxonomy = nm, marker_genes = length(z$marker_genes),
             reference_pseudobulks = ncol(z$logcpm),
             reference_patients = uniqueN(z$meta$patient_id),
             labels = uniqueN(z$meta$label))
}))

saveRDS(list(references = references, calibration_summary = summary,
             class_metrics = class_metrics, seed = seed,
             min_cells_pseudobulk = min_cells_pb,
             top_markers_per_label = top_markers_per_label,
             status = overall_status),
        file.path(outdir, "gate12g_transfer_reference.rds"), compress = "xz")
fwrite(predictions, file.path(outdir, "lopo_pseudobulk_predictions.tsv"), sep = "\t")
fwrite(class_metrics, file.path(outdir, "lopo_class_metrics.tsv"), sep = "\t")
fwrite(summary, file.path(outdir, "lopo_taxonomy_summary.tsv"), sep = "\t")
fwrite(reference_manifest, file.path(outdir, "reference_manifest.tsv"), sep = "\t")

checkpoint <- c(
  "# Gate12G discovery-reference transfer checkpoint",
  "",
  paste0("- Discovery SCE files: ", length(files)),
  paste0("- Discovery patients: ", uniqueN(rbindlist(lapply(pb_store, function(z) rbindlist(z$meta)), fill = TRUE)$patient_id)),
  paste0("- Minimum cells per reference pseudobulk: ", min_cells_pb),
  paste0("- Top positive markers per label: ", top_markers_per_label),
  paste0("- Overall transfer-calibration decision: **", overall_status, "**"),
  "",
  "External OEP005136 classification is permitted only when all three taxonomies pass their frozen patient-blocked calibration thresholds."
)
writeLines(checkpoint, file.path(outdir, "GATE12G_TRANSFER_REFERENCE_CHECKPOINT.md"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))

cat("GATE12G_TRANSFER_REFERENCE_STATUS=", overall_status, "\n", sep = "")
print(summary)
