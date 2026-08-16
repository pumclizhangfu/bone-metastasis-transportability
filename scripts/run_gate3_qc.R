#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(SingleCellExperiment)
  library(scDblFinder)
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
required <- c("input", "work", "output", "workers", "seed")
missing_args <- setdiff(required, names(cfg))
if (length(missing_args)) stop("Missing arguments: ", paste(missing_args, collapse = ", "))

input_dir <- normalizePath(cfg$input, mustWork = TRUE)
work_dir <- normalizePath(cfg$work, mustWork = FALSE)
output_dir <- normalizePath(cfg$output, mustWork = FALSE)
workers <- as.integer(cfg$workers)
seed <- as.integer(cfg$seed)
if (is.na(workers) || workers < 1 || workers > 20) stop("workers must be 1..20")
if (is.na(seed)) stop("seed must be an integer")

dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
extract_dir <- file.path(work_dir, "extracted_counts")
sce_dir <- file.path(work_dir, "sce_singlets")
dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sce_dir, recursive = TRUE, showWarnings = FALSE)

message("CONFIG\tinput=", input_dir, "\twork=", work_dir, "\toutput=", output_dir,
        "\tworkers=", workers, "\tseed=", seed)

cohorts <- data.frame(
  accession = c("GSE143791", "GSE202813"),
  cancer = c("prostate", "renal"),
  stringsAsFactors = FALSE
)

manifest_rows <- list()
for (i in seq_len(nrow(cohorts))) {
  accession <- cohorts$accession[[i]]
  archive <- file.path(input_dir, accession, paste0(accession, "_RAW.tar"))
  if (!file.exists(archive)) stop("Missing archive: ", archive)
  members <- utils::untar(archive, list = TRUE)
  if (accession == "GSE143791") {
    selected <- members[grepl("^GSM[0-9]+_BMET[0-9]+-(Tumor|Involved|Distal)\\.count\\.csv\\.gz$", members)]
  } else {
    selected <- members[grepl("^GSM[0-9]+_RCC-BM[0-9]+-(Tumor|Involve|Noninvolved)\\.count\\.csv\\.gz$", members)]
  }
  if (!length(selected)) stop("No selected matrices in ", archive)
  cohort_extract <- file.path(extract_dir, accession)
  dir.create(cohort_extract, recursive = TRUE, showWarnings = FALSE)
  missing_files <- selected[!file.exists(file.path(cohort_extract, selected))]
  if (length(missing_files)) {
    message("EXTRACT\t", accession, "\tfiles=", length(missing_files))
    utils::untar(archive, files = missing_files, exdir = cohort_extract)
  }
  for (member in selected) {
    title <- sub("^GSM[0-9]+_", "", member)
    title <- sub("\\.count\\.csv\\.gz$", "", title)
    if (accession == "GSE143791") {
      patient <- sub("-.*$", "", title)
      compartment <- tolower(sub("^.*-", "", title))
    } else {
      patient <- sub("^RCC-", "", sub("-(Tumor|Involve|Noninvolved)$", "", title))
      raw_compartment <- sub("^.*-", "", title)
      compartment <- c(Tumor = "tumor", Involve = "involved", Noninvolved = "distal")[[raw_compartment]]
    }
    manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
      accession = accession,
      cancer = cohorts$cancer[[i]],
      sample_id = title,
      patient_id = patient,
      compartment = compartment,
      path = file.path(cohort_extract, member),
      stringsAsFactors = FALSE
    )
  }
}
manifest <- rbindlist(manifest_rows)
setorder(manifest, accession, patient_id, compartment)
fwrite(manifest, file.path(output_dir, "input_manifest.tsv"), sep = "\t")
message("MANIFEST\tsamples=", nrow(manifest), "\tpatients=", uniqueN(paste(manifest$accession, manifest$patient_id)))

read_matrix <- function(path) {
  command <- paste("gzip -dc", shQuote(path))
  tab <- fread(cmd = command, data.table = FALSE, check.names = FALSE, showProgress = FALSE)
  genes <- as.character(tab[[1]])
  tab[[1]] <- NULL
  dense <- as.matrix(tab)
  storage.mode(dense) <- "numeric"
  rownames(dense) <- make.unique(genes)
  sparse <- as(dense, "dgCMatrix")
  rm(tab, dense)
  sparse
}

process_sample <- function(meta, sample_index) {
  sample_id <- meta$sample_id[[1]]
  message("SAMPLE_START\t", sample_index, "/", nrow(manifest), "\t", sample_id)
  counts <- read_matrix(meta$path[[1]])
  raw_cells <- ncol(counts)
  n_count <- Matrix::colSums(counts)
  n_feature <- Matrix::colSums(counts > 0)
  mt <- grepl("^MT-", rownames(counts), ignore.case = FALSE)
  pct_mt <- if (any(mt)) 100 * Matrix::colSums(counts[mt, , drop = FALSE]) / pmax(n_count, 1) else rep(0, raw_cells)
  upper_feature <- median(n_feature) + 4 * mad(n_feature, constant = 1.4826)
  upper_feature <- min(10000, max(2500, upper_feature))
  pass_basic <- n_feature >= 200 & n_count >= 500 & n_feature <= upper_feature & pct_mt <= 25
  if (sum(pass_basic) < 100) stop("Fewer than 100 cells pass basic QC in ", sample_id)

  sce <- SingleCellExperiment(
    assays = list(counts = counts[, pass_basic, drop = FALSE]),
    colData = DataFrame(
      barcode = colnames(counts)[pass_basic],
      accession = meta$accession[[1]],
      cancer = meta$cancer[[1]],
      sample_id = sample_id,
      patient_id = meta$patient_id[[1]],
      compartment = meta$compartment[[1]],
      nCount = n_count[pass_basic],
      nFeature = n_feature[pass_basic],
      percent_mt = pct_mt[pass_basic]
    )
  )
  set.seed(seed + sample_index)
  bp <- if (.Platform$OS.type == "unix" && workers > 1) {
    MulticoreParam(workers = min(workers, 4), progressbar = FALSE)
  } else {
    SerialParam(progressbar = FALSE)
  }
  sce <- scDblFinder(sce, samples = rep(sample_id, ncol(sce)), BPPARAM = bp, verbose = FALSE)
  singlet <- sce$scDblFinder.class == "singlet"

  qc_all <- data.table(
    accession = meta$accession[[1]],
    cancer = meta$cancer[[1]],
    sample_id = sample_id,
    patient_id = meta$patient_id[[1]],
    compartment = meta$compartment[[1]],
    barcode = colnames(counts),
    nCount = as.numeric(n_count),
    nFeature = as.numeric(n_feature),
    percent_mt = as.numeric(pct_mt),
    upper_feature_threshold = upper_feature,
    pass_basic_qc = pass_basic,
    doublet_class = NA_character_,
    keep_singlet = FALSE
  )
  qc_all[pass_basic, doublet_class := as.character(sce$scDblFinder.class)]
  qc_all[pass_basic, keep_singlet := singlet]

  sce_singlet <- sce[, singlet, drop = FALSE]
  saveRDS(sce_singlet, file.path(sce_dir, paste0(sample_id, ".rds")), compress = "xz")
  summary <- data.table(
    accession = meta$accession[[1]],
    cancer = meta$cancer[[1]],
    sample_id = sample_id,
    patient_id = meta$patient_id[[1]],
    compartment = meta$compartment[[1]],
    raw_cells = raw_cells,
    basic_qc_cells = sum(pass_basic),
    singlet_cells = sum(singlet),
    basic_qc_fraction = sum(pass_basic) / raw_cells,
    singlet_fraction_of_basic = mean(singlet),
    median_nCount_singlet = median(sce_singlet$nCount),
    median_nFeature_singlet = median(sce_singlet$nFeature),
    median_percent_mt_singlet = median(sce_singlet$percent_mt),
    upper_feature_threshold = upper_feature,
    analyzable_100_cells = sum(singlet) >= 100
  )
  message("SAMPLE_COMPLETE\t", sample_id, "\traw=", raw_cells, "\tbasic=", sum(pass_basic), "\tsinglets=", sum(singlet))
  rm(counts, sce, sce_singlet)
  gc(verbose = FALSE)
  list(summary = summary, cells = qc_all)
}

sample_results <- vector("list", nrow(manifest))
for (i in seq_len(nrow(manifest))) {
  sample_results[[i]] <- process_sample(manifest[i], i)
}
sample_qc <- rbindlist(lapply(sample_results, `[[`, "summary"))
cell_qc <- rbindlist(lapply(sample_results, `[[`, "cells"))
fwrite(sample_qc, file.path(output_dir, "sample_qc.tsv"), sep = "\t")
fwrite(cell_qc, file.path(output_dir, "cell_qc.tsv.gz"), sep = "\t", compress = "gzip")

triplets <- sample_qc[, .(
  compartments = paste(sort(unique(compartment)), collapse = ","),
  n_compartments = uniqueN(compartment),
  min_singlets = min(singlet_cells),
  all_compartments_analyzable = all(analyzable_100_cells),
  complete_triplet = all(c("tumor", "involved", "distal") %in% compartment) && all(analyzable_100_cells)
), by = .(accession, cancer, patient_id)]
setorder(triplets, accession, patient_id)
fwrite(triplets, file.path(output_dir, "patient_triplet_audit.tsv"), sep = "\t")

criteria <- c(
  intended_samples_42 = nrow(sample_qc) == 42,
  complete_triplets_at_least_8 = sum(triplets$complete_triplet) >= 8,
  singlets_at_least_50000 = sum(sample_qc$singlet_cells) >= 50000,
  labels_complete = !anyNA(sample_qc[, .(accession, cancer, patient_id, compartment)]) &&
    all(sample_qc$accession != "" & sample_qc$cancer != "" & sample_qc$patient_id != "" & sample_qc$compartment != "")
)
decision <- if (all(criteria)) "GO" else "NO-GO"
decision_lines <- c(
  "# Gate 3A decision",
  "",
  paste0("**", decision, "**"),
  "",
  paste0("- Intended matrices readable: ", nrow(sample_qc), "/42 — ", ifelse(criteria[[1]], "PASS", "FAIL")),
  paste0("- Complete analyzable patient triplets: ", sum(triplets$complete_triplet), " (threshold 8) — ", ifelse(criteria[[2]], "PASS", "FAIL")),
  paste0("- Retained singlets: ", sum(sample_qc$singlet_cells), " (threshold 50,000) — ", ifelse(criteria[[3]], "PASS", "FAIL")),
  paste0("- Complete labels: ", ifelse(criteria[[4]], "PASS", "FAIL")),
  "",
  "This decision concerns data analyzability only; it is not evidence for a biological gradient."
)
writeLines(decision_lines, file.path(output_dir, "gate3a_decision.md"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("GATE3A_COMPLETE\tdecision=", decision, "\tsamples=", nrow(sample_qc),
        "\ttriplets=", sum(triplets$complete_triplet), "\tsinglets=", sum(sample_qc$singlet_cells))
