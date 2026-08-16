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
  for (i in seq(1, length(x), by = 2)) {
    key <- sub("^--", "", x[[i]])
    out[[key]] <- x[[i + 1]]
  }
  out
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("input", "annotations", "ontology", "work", "output", "workers", "seed", "min-cells")
missing_args <- setdiff(required, names(cfg))
if (length(missing_args)) stop("Missing arguments: ", paste(missing_args, collapse = ", "))

input_dir <- normalizePath(cfg$input, mustWork = TRUE)
annotation_file <- normalizePath(cfg$annotations, mustWork = TRUE)
ontology_file <- normalizePath(cfg$ontology, mustWork = TRUE)
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
annotated_sce_dir <- file.path(work_dir, "annotated_sce")
dir.create(annotated_sce_dir, recursive = TRUE, showWarnings = FALSE)

message("CONFIG\tinput=", input_dir, "\tannotations=", annotation_file,
        "\tontology=", ontology_file, "\twork=", work_dir,
        "\toutput=", output_dir, "\tworkers=", workers,
        "\tseed=", seed, "\tmin_cells=", min_cells)

set.seed(seed)
bp <- if (.Platform$OS.type == "unix" && workers > 1) {
  MulticoreParam(workers = workers, progressbar = FALSE)
} else {
  SerialParam(progressbar = FALSE)
}

ontology <- fread(ontology_file)
required_ontology <- c("author_label", "broad_class", "harmonized_state", "primary_pseudobulk")
if (!all(required_ontology %in% names(ontology))) stop("Ontology is missing required columns")
if (anyDuplicated(ontology$author_label)) stop("Ontology author_label values must be unique")
label_to_broad <- setNames(ontology$broad_class, ontology$author_label)
label_to_state <- setNames(ontology$harmonized_state, ontology$author_label)

sce_files <- sort(list.files(input_dir, pattern = "\\.rds$", full.names = TRUE))
if (length(sce_files) != 42) stop("Expected 42 Gate3A SCE files; found ", length(sce_files))
sce_list <- lapply(sce_files, readRDS)
names(sce_list) <- sub("\\.rds$", "", basename(sce_files))
if (!all(vapply(sce_list, function(x) inherits(x, "SingleCellExperiment"), logical(1)))) {
  stop("Every Gate3A input must be a SingleCellExperiment")
}
if (!all(vapply(sce_list[-1], function(x) identical(rownames(x), rownames(sce_list[[1]])), logical(1)))) {
  stop("Gene rows are not identical across Gate3A SCE files")
}
all_barcodes <- unlist(lapply(sce_list, colnames), use.names = FALSE)
if (anyDuplicated(all_barcodes)) stop("Cell barcodes are not globally unique")
message("INPUT_AUDIT\tsamples=", length(sce_list), "\tcells=", length(all_barcodes),
        "\tgenes=", nrow(sce_list[[1]]))

author_labels <- readRDS(gzfile(annotation_file, open = "rb"))
if (!is.character(author_labels) || is.null(names(author_labels))) {
  stop("Author annotations must be a named character vector")
}
unknown_author_labels <- setdiff(unique(author_labels), ontology$author_label)
if (length(unknown_author_labels)) {
  stop("Ontology does not cover author labels: ", paste(unknown_author_labels, collapse = ", "))
}

marker_modules <- list(
  T_NK = c("CD3D", "CD3E", "TRAC", "CD247", "IL32", "LTB", "TRBC1", "TRBC2"),
  B = c("CD79A", "CD79B", "MS4A1", "CD19", "CD22", "CD37", "CD74", "CD83"),
  Myeloid = c("LST1", "TYROBP", "FCER1G", "CTSS", "CTSB", "LILRB1", "LYZ", "CTSD"),
  Progenitor = c("CD34", "SPINK2", "HLF", "GATA2", "MEIS1", "PROM1", "KIT", "SOX4"),
  Erythroid = c("HBB", "HBA1", "HBA2", "GYPA", "ALAS2", "AHSP", "CA1", "SLC4A1"),
  Stromal = c("COL1A1", "COL1A2", "DCN", "COL3A1", "CXCL12", "LEPR", "PDGFRA", "COL6A1"),
  Endothelial = c("PECAM1", "VWF", "EMCN", "CLDN5", "KDR", "ENG", "RAMP2", "ESAM"),
  Osteoclast = c("ACP5", "CTSK", "ATP6V0D2", "OCSTAMP", "TM7SF4", "MMP9", "CALCR", "NFATC1"),
  Osteoblast = c("BGLAP", "ALPL", "SP7", "IBSP", "RUNX2", "SPP1", "COL1A1", "BMP2"),
  Malignant = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "KLK2", "KLK3", "ACPP", "FOLH1", "CA9")
)
module_table <- rbindlist(lapply(names(marker_modules), function(module) {
  data.table(module = module, gene = marker_modules[[module]])
}))
fwrite(module_table, file.path(output_dir, "marker_modules.tsv"), sep = "\t")

empty_labels <- function(n) {
  data.table(
    author_label = rep(NA_character_, n),
    broad_class = rep("Unassigned", n),
    harmonized_state = rep("Unassigned", n),
    label_source = rep("Unassigned", n),
    singler_label = rep(NA_character_, n),
    singler_pruned_label = rep(NA_character_, n),
    singler_delta_next = rep(NA_real_, n),
    prostate_marker_genes = rep(NA_integer_, n),
    prostate_marker_umis = rep(NA_real_, n),
    epithelial_marker_genes = rep(NA_integer_, n)
  )
}

# First pass: preserve published RCC per-cell labels without imputing unmatched cells.
for (sample_id in names(sce_list)) {
  sce <- sce_list[[sample_id]]
  if (unique(sce$accession) != "GSE202813") next
  lab <- empty_labels(ncol(sce))
  matched_author <- unname(author_labels[colnames(sce)])
  matched <- !is.na(matched_author)
  lab$author_label[matched] <- matched_author[matched]
  lab$broad_class[matched] <- unname(label_to_broad[matched_author[matched]])
  lab$harmonized_state[matched] <- unname(label_to_state[matched_author[matched]])
  lab$label_source[matched] <- "RCC_author_per_cell"
  sce$author_label <- lab$author_label
  sce$broad_class <- lab$broad_class
  sce$harmonized_state <- lab$harmonized_state
  sce$label_source <- lab$label_source
  sce$singler_label <- lab$singler_label
  sce$singler_pruned_label <- lab$singler_pruned_label
  sce$singler_delta_next <- lab$singler_delta_next
  sce$prostate_marker_genes <- lab$prostate_marker_genes
  sce$prostate_marker_umis <- lab$prostate_marker_umis
  sce$epithelial_marker_genes <- lab$epithelial_marker_genes
  sce_list[[sample_id]] <- sce
}

# Build independent RCC author-labelled pseudobulk reference replicates.
ref_columns <- list()
ref_metadata <- list()
for (sample_id in names(sce_list)) {
  sce <- sce_list[[sample_id]]
  if (unique(sce$accession) != "GSE202813") next
  valid_classes <- setdiff(unique(sce$broad_class), "Unassigned")
  for (broad in valid_classes) {
    idx <- which(sce$broad_class == broad)
    if (length(idx) < min_cells) next
    ref_id <- paste(sample_id, broad, sep = "::")
    ref_columns[[ref_id]] <- Matrix::rowSums(counts(sce)[, idx, drop = FALSE])
    ref_metadata[[ref_id]] <- data.table(
      reference_id = ref_id,
      sample_id = sample_id,
      patient_id = unique(sce$patient_id),
      compartment = unique(sce$compartment),
      broad_class = broad,
      n_cells = length(idx)
    )
  }
}
if (!length(ref_columns)) stop("No RCC reference pseudobulks were generated")
ref_counts <- as(do.call(cbind, ref_columns), "dgCMatrix")
ref_meta <- rbindlist(ref_metadata, use.names = TRUE)
ref_meta <- ref_meta[match(colnames(ref_counts), reference_id)]
reference_classes <- sort(unique(ref_meta$broad_class))
if (length(reference_classes) < 7) stop("Fewer than 7 broad classes in RCC reference")
ref_sce <- SingleCellExperiment(
  assays = list(counts = ref_counts),
  colData = S4Vectors::DataFrame(ref_meta)
)
ref_sce <- logNormCounts(ref_sce)
fwrite(ref_meta, file.path(output_dir, "singler_reference_metadata.tsv"), sep = "\t")
saveRDS(ref_sce, file.path(work_dir, "rcc_author_pseudobulk_reference.rds"), compress = "xz")
trained_reference <- trainSingleR(
  ref = ref_sce,
  labels = ref_sce$broad_class,
  de.method = "classic",
  assay.type = "logcounts",
  BPPARAM = bp
)
saveRDS(trained_reference, file.path(work_dir, "rcc_author_trained_singler.rds"), compress = "xz")
message("REFERENCE_READY\tpseudobulks=", ncol(ref_sce),
        "\tclasses=", paste(reference_classes, collapse = ","))

prostate_specific <- c("KLK2", "KLK3", "KLK4", "ACPP", "AR", "FOLH1", "TMPRSS2", "NKX3-1")
epithelial <- c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "KRT17", "KRT4")

# Second pass: transfer only broad RCC classes to prostate cells.
for (sample_id in names(sce_list)) {
  sce <- sce_list[[sample_id]]
  if (unique(sce$accession) != "GSE143791") next
  message("ANNOTATE_START\t", sample_id, "\tcells=", ncol(sce))
  norm_sce <- logNormCounts(sce)
  set.seed(seed + match(sample_id, names(sce_list)))
  pred <- classifySingleR(
    test = norm_sce,
    trained = trained_reference,
    assay.type = "logcounts",
    fine.tune = TRUE,
    prune = TRUE,
    BPPARAM = bp
  )
  raw_label <- as.character(pred$labels)
  pruned_label <- as.character(pred$pruned.labels)
  assigned <- pruned_label
  assigned[is.na(assigned) | assigned == ""] <- "Unassigned"

  cts <- counts(sce)
  ps <- intersect(prostate_specific, rownames(sce))
  ep <- intersect(epithelial, rownames(sce))
  ps_detected <- if (length(ps)) Matrix::colSums(cts[ps, , drop = FALSE] > 0) else rep(0, ncol(sce))
  ps_umis <- if (length(ps)) Matrix::colSums(cts[ps, , drop = FALSE]) else rep(0, ncol(sce))
  ep_detected <- if (length(ep)) Matrix::colSums(cts[ep, , drop = FALSE] > 0) else rep(0, ncol(sce))
  malignant_supported <- ps_detected >= 2 | (ps_detected >= 1 & ps_umis >= 3 & ep_detected >= 1) | ep_detected >= 3
  unsupported_malignant <- assigned == "Malignant" & !malignant_supported
  assigned[unsupported_malignant] <- "Unassigned"
  marker_override <- malignant_supported & assigned != "Malignant"
  assigned[marker_override] <- "Malignant"

  source <- rep("SingleR_RCC_author_broad", ncol(sce))
  source[assigned == "Unassigned"] <- "Unassigned"
  source[marker_override] <- "prostate_marker_override"
  sce$author_label <- NA_character_
  sce$broad_class <- assigned
  sce$harmonized_state <- assigned
  sce$label_source <- source
  sce$singler_label <- raw_label
  sce$singler_pruned_label <- pruned_label
  sce$singler_delta_next <- as.numeric(pred$delta.next)
  sce$prostate_marker_genes <- as.integer(ps_detected)
  sce$prostate_marker_umis <- as.numeric(ps_umis)
  sce$epithelial_marker_genes <- as.integer(ep_detected)
  sce_list[[sample_id]] <- sce
  message("ANNOTATE_COMPLETE\t", sample_id,
          "\tlabelled=", sum(assigned != "Unassigned"),
          "\tunassigned=", sum(assigned == "Unassigned"),
          "\tmarker_overrides=", sum(marker_override))
  rm(norm_sce, pred, cts)
  gc(verbose = FALSE)
}

# Save annotated objects and a flat, auditable cell annotation table.
cell_annotation_rows <- list()
for (sample_id in names(sce_list)) {
  sce <- sce_list[[sample_id]]
  saveRDS(sce, file.path(annotated_sce_dir, paste0(sample_id, ".rds")), compress = "xz")
  cell_annotation_rows[[sample_id]] <- data.table(
    accession = as.character(sce$accession),
    cancer = as.character(sce$cancer),
    sample_id = as.character(sce$sample_id),
    patient_id = as.character(sce$patient_id),
    compartment = as.character(sce$compartment),
    barcode = colnames(sce),
    author_label = as.character(sce$author_label),
    broad_class = as.character(sce$broad_class),
    harmonized_state = as.character(sce$harmonized_state),
    label_source = as.character(sce$label_source),
    singler_label = as.character(sce$singler_label),
    singler_pruned_label = as.character(sce$singler_pruned_label),
    singler_delta_next = as.numeric(sce$singler_delta_next),
    prostate_marker_genes = as.integer(sce$prostate_marker_genes),
    prostate_marker_umis = as.numeric(sce$prostate_marker_umis),
    epithelial_marker_genes = as.integer(sce$epithelial_marker_genes)
  )
}
cell_annotations <- rbindlist(cell_annotation_rows, use.names = TRUE)
if (nrow(cell_annotations) != length(all_barcodes) || anyDuplicated(cell_annotations$barcode)) {
  stop("Cell annotation row conservation failed")
}
fwrite(cell_annotations, file.path(output_dir, "cell_annotations.tsv.gz"), sep = "\t", compress = "gzip", na = "NA")

coverage <- cell_annotations[, .(
  total_singlets = .N,
  labelled_singlets = sum(broad_class != "Unassigned"),
  unassigned_singlets = sum(broad_class == "Unassigned"),
  label_coverage = mean(broad_class != "Unassigned")
), by = .(accession, cancer)]
fwrite(coverage, file.path(output_dir, "annotation_coverage.tsv"), sep = "\t")

# Aggregate labelled groups and audit broad marker-module concordance.
marker_aggregate <- list()
marker_cells <- list()
for (sample_id in names(sce_list)) {
  sce <- sce_list[[sample_id]]
  for (broad in setdiff(unique(sce$broad_class), "Unassigned")) {
    key <- paste(unique(sce$accession), broad, sep = "::")
    idx <- which(sce$broad_class == broad)
    vec <- Matrix::rowSums(counts(sce)[, idx, drop = FALSE])
    if (is.null(marker_aggregate[[key]])) marker_aggregate[[key]] <- vec else marker_aggregate[[key]] <- marker_aggregate[[key]] + vec
    previous_n <- if (is.null(marker_cells[[key]])) 0L else marker_cells[[key]]
    marker_cells[[key]] <- previous_n + length(idx)
  }
}

marker_audit_rows <- list()
for (key in names(marker_aggregate)) {
  parts <- strsplit(key, "::", fixed = TRUE)[[1]]
  accession <- parts[[1]]
  expected <- parts[[2]]
  aggregate <- marker_aggregate[[key]]
  total <- sum(aggregate)
  log_cpm <- log2(1 + 1e6 * aggregate / max(total, 1))
  names(log_cpm) <- rownames(sce_list[[1]])
  scores <- vapply(marker_modules, function(genes) {
    present <- intersect(genes, names(log_cpm))
    if (!length(present)) return(NA_real_)
    mean(log_cpm[present])
  }, numeric(1))
  ranks <- rank(-scores, ties.method = "min", na.last = "keep")
  marker_audit_rows[[key]] <- data.table(
    accession = accession,
    cancer = ifelse(accession == "GSE143791", "prostate", "renal"),
    broad_class = expected,
    n_cells = as.integer(marker_cells[[key]]),
    expected_module_score = unname(scores[expected]),
    best_module = names(which.max(scores)),
    expected_module_rank = unname(ranks[expected]),
    pass_expected_top2 = unname(ranks[expected]) <= 2
  )
}
marker_audit <- rbindlist(marker_audit_rows, use.names = TRUE)
setorder(marker_audit, accession, broad_class)
fwrite(marker_audit, file.path(output_dir, "annotation_marker_audit.tsv"), sep = "\t")
marker_summary <- marker_audit[, .(
  labelled_cells_in_audit = sum(n_cells),
  cells_in_top2_groups = sum(n_cells[pass_expected_top2]),
  weighted_top2_fraction = sum(n_cells[pass_expected_top2]) / sum(n_cells)
), by = .(accession, cancer)]
fwrite(marker_summary, file.path(output_dir, "annotation_marker_summary.tsv"), sep = "\t")

# Aggregate raw counts to patient-compartment-broad-class pseudobulks.
pseudobulk_columns <- list()
pseudobulk_metadata_rows <- list()
for (sample_id in names(sce_list)) {
  sce <- sce_list[[sample_id]]
  for (broad in setdiff(unique(sce$broad_class), "Unassigned")) {
    idx <- which(sce$broad_class == broad)
    pb_id <- paste(sample_id, broad, sep = "::")
    cell_sum <- sum(counts(sce)[, idx, drop = FALSE])
    pb_vec <- Matrix::rowSums(counts(sce)[, idx, drop = FALSE])
    pseudobulk_columns[[pb_id]] <- pb_vec
    pseudobulk_metadata_rows[[pb_id]] <- data.table(
      pseudobulk_id = pb_id,
      accession = unique(sce$accession),
      cancer = unique(sce$cancer),
      sample_id = sample_id,
      patient_id = unique(sce$patient_id),
      compartment = unique(sce$compartment),
      broad_class = broad,
      n_cells = length(idx),
      raw_umi_sum = as.numeric(cell_sum),
      aggregation_exact = isTRUE(all.equal(as.numeric(sum(pb_vec)), as.numeric(cell_sum), tolerance = 0))
    )
  }
}
pseudobulk_counts <- as(do.call(cbind, pseudobulk_columns), "dgCMatrix")
pseudobulk_metadata <- rbindlist(pseudobulk_metadata_rows, use.names = TRUE)
pseudobulk_metadata <- pseudobulk_metadata[match(colnames(pseudobulk_counts), pseudobulk_id)]

patient_design <- unique(pseudobulk_metadata[, .(accession, cancer, patient_id, compartment)])
complete_patients <- patient_design[, .(
  complete_triplet = all(c("tumor", "involved", "distal") %in% compartment)
), by = .(accession, cancer, patient_id)]
pseudobulk_metadata <- merge(pseudobulk_metadata, complete_patients,
                             by = c("accession", "cancer", "patient_id"), all.x = TRUE, sort = FALSE)
pseudobulk_metadata[, eligible_min_cells := complete_triplet & n_cells >= min_cells]
pseudobulk_metadata <- pseudobulk_metadata[match(colnames(pseudobulk_counts), pseudobulk_id)]
if (!identical(pseudobulk_metadata$pseudobulk_id, colnames(pseudobulk_counts))) {
  stop("Pseudobulk metadata is not aligned to count columns")
}
if (!all(pseudobulk_metadata$aggregation_exact)) stop("Pseudobulk raw-count conservation failed")
fwrite(pseudobulk_metadata, file.path(output_dir, "pseudobulk_metadata.tsv"), sep = "\t")

patient_class_eligibility <- pseudobulk_metadata[, .(
  complete_eligible_triplet = all(c("tumor", "involved", "distal") %in% compartment[eligible_min_cells])
), by = .(accession, cancer, patient_id, broad_class)]
eligibility <- patient_class_eligibility[, .(
  n_complete_eligible_patients = sum(complete_eligible_triplet),
  paired_model_ready = sum(complete_eligible_triplet) >= 3
), by = .(accession, cancer, broad_class)]
cross_ready <- eligibility[, .(
  cancers_ready = uniqueN(accession[paired_model_ready]),
  cross_cancer_ready = uniqueN(accession[paired_model_ready]) == 2
), by = broad_class]
eligibility <- merge(eligibility, cross_ready, by = "broad_class", all.x = TRUE)
setorder(eligibility, broad_class, accession)
fwrite(eligibility, file.path(output_dir, "pseudobulk_eligibility.tsv"), sep = "\t")

saveRDS(
  list(
    counts = pseudobulk_counts,
    metadata = as.data.frame(pseudobulk_metadata),
    ontology = as.data.frame(ontology),
    marker_modules = marker_modules,
    seed = seed,
    min_cells = min_cells
  ),
  file.path(output_dir, "pseudobulk_counts.rds"),
  compress = "xz"
)

rcc_coverage <- coverage[accession == "GSE202813", label_coverage]
prostate_coverage <- coverage[accession == "GSE143791", label_coverage]
minimum_marker_concordance <- min(marker_summary$weighted_top2_fraction)
n_cross_cancer_classes <- uniqueN(eligibility[cross_cancer_ready == TRUE, broad_class])
criteria <- c(
  rcc_author_coverage_at_least_97pct = length(rcc_coverage) == 1 && rcc_coverage >= 0.97,
  prostate_label_coverage_at_least_85pct = length(prostate_coverage) == 1 && prostate_coverage >= 0.85,
  marker_concordance_at_least_90pct = is.finite(minimum_marker_concordance) && minimum_marker_concordance >= 0.90,
  cross_cancer_classes_at_least_5 = n_cross_cancer_classes >= 5,
  pseudobulk_aggregation_exact = all(pseudobulk_metadata$aggregation_exact)
)
decision <- if (all(criteria)) "GO" else "NO-GO"
decision_lines <- c(
  "# Gate 3B decision",
  "",
  paste0("**", decision, "**"),
  "",
  sprintf("- RCC author-label coverage: %.2f%% (threshold 97%%) - %s", 100 * rcc_coverage, ifelse(criteria[[1]], "PASS", "FAIL")),
  sprintf("- Prostate conservative broad-label coverage: %.2f%% (threshold 85%%) - %s", 100 * prostate_coverage, ifelse(criteria[[2]], "PASS", "FAIL")),
  sprintf("- Minimum cohort marker top-two concordance: %.2f%% (threshold 90%%) - %s", 100 * minimum_marker_concordance, ifelse(criteria[[3]], "PASS", "FAIL")),
  sprintf("- Broad classes ready in both cancers: %d (threshold 5) - %s", n_cross_cancer_classes, ifelse(criteria[[4]], "PASS", "FAIL")),
  paste0("- Exact raw-count aggregation: ", ifelse(criteria[[5]], "PASS", "FAIL")),
  "",
  "This gate validates annotation and statistical-unit construction only; it does not establish a biological gradient."
)
writeLines(decision_lines, file.path(output_dir, "gate3b_decision.md"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("GATE3B_COMPLETE\tdecision=", decision,
        "\trcc_coverage=", signif(rcc_coverage, 4),
        "\tprostate_coverage=", signif(prostate_coverage, 4),
        "\tmarker_concordance=", signif(minimum_marker_concordance, 4),
        "\tcross_cancer_classes=", n_cross_cancer_classes,
        "\tpseudobulks=", ncol(pseudobulk_counts))
