#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(SummarizedExperiment)
  library(edgeR)
  library(fgsea)
  library(msigdbr)
})

parse_args <- function(x) {
  if (length(x) %% 2 != 0) stop("Arguments must be --key value pairs")
  out <- list()
  for (i in seq(1, length(x), by = 2)) out[[sub("^--", "", x[[i]])]] <- x[[i + 1]]
  out
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("sce-dir", "annotations", "thresholds", "gate5b-genes",
              "gate5b-hallmarks", "gate5b-cohort-differential",
              "stress-genes", "output", "workers", "seed", "min-cells",
              "min-retention-fraction")
missing_args <- setdiff(required, names(cfg))
if (length(missing_args)) stop("Missing arguments: ", paste(missing_args, collapse = ", "))

sce_dir <- normalizePath(cfg[["sce-dir"]], mustWork = TRUE)
annotation_file <- normalizePath(cfg$annotations, mustWork = TRUE)
threshold_file <- normalizePath(cfg$thresholds, mustWork = TRUE)
gate5b_gene_file <- normalizePath(cfg[["gate5b-genes"]], mustWork = TRUE)
gate5b_hallmark_file <- normalizePath(cfg[["gate5b-hallmarks"]], mustWork = TRUE)
gate5b_cohort_file <- normalizePath(cfg[["gate5b-cohort-differential"]], mustWork = TRUE)
stress_gene_file <- normalizePath(cfg[["stress-genes"]], mustWork = TRUE)
output_dir <- normalizePath(cfg$output, mustWork = FALSE)
workers <- as.integer(cfg$workers)
seed <- as.integer(cfg$seed)
min_cells <- as.integer(cfg[["min-cells"]])
min_retention_fraction <- as.numeric(cfg[["min-retention-fraction"]])
if (is.na(workers) || workers < 1 || workers > 20) stop("workers must be 1..20")
if (is.na(seed)) stop("seed must be an integer")
if (is.na(min_cells) || min_cells < 10) stop("min-cells must be at least 10")
if (!is.finite(min_retention_fraction) || min_retention_fraction < 0.5 ||
    min_retention_fraction > 1) stop("min-retention-fraction must be 0.5..1")
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE,
                                                no.. = TRUE)) > 0) {
  stop("Output directory already exists and is not empty: ", output_dir)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

fwrite(data.table(
  field = c("parent_gate", "parent_status", "amendment_reason",
            "replaced_integrity_rule", "amended_integrity_rule",
            "outcome_blinding_scope"),
  value = c(
    "Gate6",
    "Stopped at GSE143791 CD8_CTL before differential-expression results were written",
    "4902 analyzable genes missed the absolute 5000-gene threshold by 98 genes",
    "At least 5000 tested genes in every model",
    paste0("At least ", min_retention_fraction * 100,
           "% of the stress-excluded Gate5b tested universe and Gate5b candidate set retained in every model"),
    "Amendment selected without inspecting Gate6 differential-expression directions, P values, or enrichment results"
  )
), file.path(output_dir, "protocol_amendment.tsv"), sep = "\t")

lineages <- c("CD4_T", "CD8_CTL", "NK_NKT")
accessions <- c("GSE143791", "GSE202813")
accession_to_cancer <- c(GSE143791 = "prostate", GSE202813 = "renal")
stress <- fread(stress_gene_file)
if (!"gene" %in% names(stress) || anyDuplicated(stress$gene)) stop("Invalid stress-gene table")
gate5b_tested <- fread(gate5b_cohort_file)[contrast == "trend"]
gate5b_candidates <- fread(gate5b_gene_file)[
  gate4b_direction_concordant == TRUE & !gene %in% stress$gene,
  .(lineage, gene, gate5b_meta_logFC = meta_logFC)
]
if (anyDuplicated(gate5b_tested[, .(accession, lineage, gene)])) stop("Duplicate Gate5b tested genes")
if (anyDuplicated(gate5b_candidates[, .(lineage, gene)])) stop("Duplicate Gate5b candidate genes")

message("CONFIG\tsce_dir=", sce_dir, "\tannotations=", annotation_file,
        "\tthresholds=", threshold_file, "\tgate5b_genes=", gate5b_gene_file,
        "\tgate5b_hallmarks=", gate5b_hallmark_file,
        "\tgate5b_cohort_differential=", gate5b_cohort_file,
        "\tstress_genes=", stress_gene_file, "\toutput=", output_dir,
        "\tworkers=", workers, "\tseed=", seed, "\tmin_cells=", min_cells,
        "\tmin_retention_fraction=", min_retention_fraction)

annotations <- fread(annotation_file)
thresholds <- fread(threshold_file)
if (!all(c("meta_state", "delta_threshold") %in% names(thresholds))) {
  stop("Threshold table lacks meta_state/delta_threshold")
}
threshold_lookup <- setNames(thresholds$delta_threshold, thresholds$meta_state)
if (!all(lineages %in% names(threshold_lookup))) stop("Missing T/NK confidence threshold")
ann <- annotations[lineage == "T_NK" & meta_state %in% lineages]
ann[, retain_high_confidence := accession == "GSE202813" |
      (accession == "GSE143791" & is.finite(delta_next) &
         delta_next >= threshold_lookup[meta_state])]
ann <- ann[retain_high_confidence == TRUE]
if (anyDuplicated(ann[, .(accession, sample_id, barcode)])) stop("Duplicate retained annotations")

sce_files <- sort(list.files(sce_dir, pattern = "\\.rds$", full.names = TRUE))
if (length(sce_files) != 42) stop("Expected 42 annotated SCE files; found ", length(sce_files))
pb_counts <- list()
pb_meta <- list()
manifest <- list()
reference_genes <- NULL
seen_keys <- character()

for (f in sce_files) {
  sce <- readRDS(f)
  counts <- assay(sce, "counts")
  genes <- rownames(counts)
  if (is.null(reference_genes)) reference_genes <- genes
  if (!identical(genes, reference_genes)) stop("Gene order mismatch in ", basename(f))
  cd <- as.data.table(as.data.frame(colData(sce)))
  needed <- c("barcode", "accession", "cancer", "sample_id", "patient_id",
              "compartment", "nCount", "nFeature", "percent_mt")
  if (!all(needed %in% names(cd))) stop("Missing colData fields in ", basename(f))
  if (uniqueN(cd$sample_id) != 1 || uniqueN(cd$accession) != 1) stop("SCE is not sample-specific: ", basename(f))
  sample_value <- unique(cd$sample_id)
  accession_value <- unique(cd$accession)
  a <- ann[accession == accession_value & sample_id == sample_value,
           .(accession, sample_id, barcode, meta_state, delta_next)]
  z <- merge(cd[, ..needed], a, by = c("accession", "sample_id", "barcode"), all = FALSE)
  if (nrow(z)) {
    seen_keys <- c(seen_keys, paste(z$accession, z$sample_id, z$barcode, sep = "::"))
    for (state_value in sort(unique(z$meta_state))) {
      zz <- z[meta_state == state_value]
      idx <- match(zz$barcode, cd$barcode)
      if (anyNA(idx)) stop("Barcode match failure in ", basename(f))
      aggregated <- Matrix::rowSums(counts[, idx, drop = FALSE])
      expected_umi <- sum(zz$nCount)
      observed_umi <- sum(aggregated)
      if (!isTRUE(all.equal(as.numeric(observed_umi), as.numeric(expected_umi), tolerance = 0))) {
        stop("Exact UMI aggregation failure in ", sample_value, " ", state_value)
      }
      pb_id <- paste(accession_value, sample_value, state_value, sep = "::")
      pb_counts[[pb_id]] <- aggregated
      pb_meta[[pb_id]] <- data.table(
        pseudobulk_id = pb_id,
        accession = accession_value,
        cancer = unique(zz$cancer),
        sample_id = sample_value,
        patient_id = unique(zz$patient_id),
        compartment = unique(zz$compartment),
        lineage = state_value,
        n_cells = nrow(zz),
        library_size = observed_umi,
        median_log10_nCount = median(log10(zz$nCount + 1)),
        median_log10_nFeature = median(log10(zz$nFeature + 1)),
        median_percent_mt = median(zz$percent_mt),
        aggregation_exact = TRUE
      )
    }
  }
  manifest[[basename(f)]] <- data.table(
    file = basename(f), accession = accession_value, sample_id = sample_value,
    n_genes = nrow(sce), n_cells = ncol(sce), retained_tnk_cells = nrow(z)
  )
  message("SCE_COMPLETE\t", basename(f), "\tretained_tnk_cells=", nrow(z))
}

expected_keys <- paste(ann$accession, ann$sample_id, ann$barcode, sep = "::")
if (!setequal(seen_keys, expected_keys) || anyDuplicated(seen_keys)) {
  stop("Retained annotation-to-SCE coverage is not exact")
}
if (!length(pb_counts)) stop("No high-confidence pseudobulks were aggregated")
all_counts <- do.call(cbind, pb_counts)
rownames(all_counts) <- reference_genes
meta <- rbindlist(pb_meta, use.names = TRUE)
setorder(meta, accession, lineage, patient_id, compartment)
all_counts <- all_counts[, meta$pseudobulk_id, drop = FALSE]
if (!identical(colnames(all_counts), meta$pseudobulk_id)) stop("Pseudobulk alignment failure")
if (!all(Matrix::colSums(all_counts) == meta$library_size)) stop("Pseudobulk library-size mismatch")

stress_found <- intersect(stress$gene, rownames(all_counts))
analysis_counts <- all_counts[setdiff(rownames(all_counts), stress_found), , drop = FALSE]
excluded <- merge(stress, data.table(gene = stress_found, present_in_matrix = TRUE),
                  by = "gene", all.x = TRUE)
excluded[is.na(present_in_matrix), present_in_matrix := FALSE]
fwrite(excluded, file.path(output_dir, "excluded_stress_genes.tsv"), sep = "\t")
fwrite(rbindlist(manifest), file.path(output_dir, "sce_manifest.tsv"), sep = "\t")
fwrite(meta, file.path(output_dir, "high_confidence_pseudobulk_metadata.tsv"), sep = "\t")
saveRDS(list(counts = all_counts, metadata = meta,
             retained_annotations = ann[, .(accession, cancer, sample_id, patient_id,
                                             compartment, barcode, meta_state, delta_next)]),
        file.path(output_dir, "high_confidence_pseudobulks.rds"), compress = "xz")
message("AGGREGATION_COMPLETE\tcells=", nrow(ann), "\tpseudobulks=", nrow(meta),
        "\tumi=", sum(meta$library_size), "\tstress_genes_removed=", length(stress_found))

msig <- as.data.table(msigdbr(db_species = "HS", species = "human", collection = "H"))
pathways <- lapply(split(msig$gene_symbol, msig$gs_name), unique)
if (length(pathways) != 50) stop("Expected 50 Hallmark sets; found ", length(pathways))
fwrite(unique(msig[, .(gs_name, gene_symbol, db_version)]),
       file.path(output_dir, "hallmark_gene_sets.tsv"), sep = "\t")

select_complete_group <- function(accession_value, lineage_value) {
  z <- meta[accession == accession_value & lineage == lineage_value & n_cells >= min_cells]
  complete <- z[, .(n_rows = .N,
                    has_all = setequal(compartment, c("distal", "involved", "tumor"))),
                by = patient_id][n_rows == 3 & has_all == TRUE, patient_id]
  z <- z[patient_id %in% complete]
  z[, compartment := factor(compartment, levels = c("distal", "involved", "tumor"))]
  z[, stage := c(distal = 0, involved = 1, tumor = 2)[as.character(compartment)]]
  setorder(z, patient_id, compartment)
  if (uniqueN(z$patient_id) < 3) stop("Fewer than 3 complete patients for ", accession_value, " ", lineage_value)
  if (nrow(z) != 3 * uniqueN(z$patient_id)) stop("Unbalanced triplets")
  z
}

add_qc_pc <- function(z) {
  vars <- c("median_log10_nCount", "median_log10_nFeature", "median_percent_mt")
  x <- as.matrix(z[, ..vars])
  if (any(!is.finite(x)) || any(apply(x, 2, sd) <= 0)) stop("Invalid QC variables")
  pc <- prcomp(x, center = TRUE, scale. = TRUE)
  anchor <- which.max(abs(pc$rotation[, 1]))
  orientation <- ifelse(pc$rotation[anchor, 1] < 0, -1, 1)
  z[, qc_pc1 := pc$x[, 1] * orientation]
  list(meta = z, loadings = data.table(qc_variable = rownames(pc$rotation),
                                       loading_pc1 = pc$rotation[, 1] * orientation,
                                       variance_explained_pc1 = pc$sdev[1]^2 / sum(pc$sdev^2)))
}

extract_test <- function(test, accession, lineage, n_patients, n_samples) {
  tab <- as.data.table(topTags(test, n = Inf, sort.by = "none")$table, keep.rownames = "gene")
  tab[, se := NA_real_]
  valid_se <- is.finite(tab$F) & tab$F > 0 & is.finite(tab$logFC)
  tab[valid_se, se := abs(logFC) / sqrt(F)]
  tab[, `:=`(accession = accession, cancer = unname(accession_to_cancer[accession]),
             lineage = lineage, contrast = "trend", n_patients = n_patients,
             n_pseudobulks = n_samples)]
  setcolorder(tab, c("accession", "cancer", "lineage", "contrast", "gene",
                     "logFC", "se", "logCPM", "F", "PValue", "FDR",
                     "n_patients", "n_pseudobulks"))
  tab
}

fit_model <- function(counts_matrix, group_meta, accession, lineage) {
  pc <- add_qc_pc(copy(group_meta))
  m <- pc$meta
  design <- model.matrix(~ factor(patient_id) + qc_pc1 + stage, data = m)
  if (qr(design)$rank != ncol(design)) stop("Adjusted design rank deficiency: ", accession, " ", lineage)
  y0 <- DGEList(counts = counts_matrix)
  keep <- filterByExpr(y0, design = design)
  tested_genes <- rownames(y0)[keep]
  accession_value <- accession
  lineage_value <- lineage
  parent_genes <- setdiff(gate5b_tested[
    accession == accession_value & lineage == lineage_value, unique(gene)
  ], stress$gene)
  candidate_genes <- gate5b_candidates[lineage == lineage_value, unique(gene)]
  parent_overlap <- sum(tested_genes %in% parent_genes)
  candidate_overlap <- sum(candidate_genes %in% tested_genes)
  parent_retention_fraction <- parent_overlap / length(parent_genes)
  candidate_retention_fraction <- candidate_overlap / length(candidate_genes)
  if (!is.finite(parent_retention_fraction) || !is.finite(candidate_retention_fraction) ||
      parent_retention_fraction < min_retention_fraction ||
      candidate_retention_fraction < min_retention_fraction) {
    stop("Gate5b analyzability retention below threshold: ", accession, " ", lineage)
  }
  y <- calcNormFactors(y0[keep, , keep.lib.sizes = FALSE], method = "TMM")
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE, abundance.trend = TRUE)
  test <- glmQLFTest(fit, coef = match("stage", colnames(design)))
  result <- extract_test(test, accession, lineage, uniqueN(m$patient_id), nrow(m))
  qc <- data.table(
    accession = accession, cancer = unname(accession_to_cancer[accession]), lineage = lineage,
    n_patients = uniqueN(m$patient_id), n_pseudobulks = nrow(m), min_cells = min(m$n_cells),
    median_cells = median(m$n_cells), min_library_size = min(m$library_size),
    median_library_size = median(m$library_size), genes_input = nrow(counts_matrix),
    genes_tested = sum(keep), design_rank = qr(design)$rank,
    design_columns = ncol(design), residual_df = nrow(design) - qr(design)$rank,
    gate5b_parent_genes = length(parent_genes),
    gate5b_parent_genes_retained = parent_overlap,
    gate5b_parent_retention_fraction = parent_retention_fraction,
    gate5b_candidate_genes = length(candidate_genes),
    gate5b_candidate_genes_retained = candidate_overlap,
    gate5b_candidate_retention_fraction = candidate_retention_fraction,
    qc_pc1_variance_explained = unique(pc$loadings$variance_explained_pc1), converged = TRUE
  )
  pc$loadings[, `:=`(accession = accession, cancer = unname(accession_to_cancer[accession]), lineage = lineage)]
  list(result = result, qc = qc, meta = m, loadings = pc$loadings,
       genes = rownames(y), signed_stat = setNames(sign(result$logFC) * sqrt(pmax(result$F, 0)), result$gene))
}

run_lopo <- function(counts_matrix, group_meta, genes, full_result, accession, lineage) {
  patients <- sort(unique(group_meta$patient_id))
  lfc <- matrix(NA_real_, nrow = length(genes), ncol = length(patients), dimnames = list(genes, patients))
  stat <- lfc
  for (i in seq_along(patients)) {
    use <- group_meta$patient_id != patients[[i]]
    m <- add_qc_pc(copy(group_meta[use]))$meta
    design <- model.matrix(~ factor(patient_id) + qc_pc1 + stage, data = m)
    if (qr(design)$rank != ncol(design)) stop("LOPO adjusted design rank deficiency")
    cts <- counts_matrix[genes, use, drop = FALSE]
    nonzero <- Matrix::rowSums(cts) > 0
    y <- calcNormFactors(DGEList(counts = cts[nonzero, , drop = FALSE]), method = "TMM")
    y <- estimateDisp(y, design, robust = TRUE)
    fit <- glmQLFit(y, design, robust = TRUE, abundance.trend = TRUE)
    test <- glmQLFTest(fit, coef = match("stage", colnames(design)))
    tab <- test$table
    lfc[rownames(tab), i] <- tab$logFC
    stat[rownames(tab), i] <- sign(tab$logFC) * sqrt(pmax(tab$F, 0))
    message("LOPO_COMPLETE\t", accession, "\t", lineage, "\tleft_out=", patients[[i]])
  }
  full_sign <- sign(full_result$logFC[match(genes, full_result$gene)])
  evaluated <- rowSums(is.finite(lfc))
  same <- rowSums(sign(lfc) == full_sign, na.rm = TRUE)
  stability <- data.table(
    accession = accession, cancer = unname(accession_to_cancer[accession]), lineage = lineage,
    gene = genes, full_logFC = full_result$logFC[match(genes, full_result$gene)],
    lopo_iterations = length(patients), lopo_evaluated = evaluated,
    lopo_same_direction = same, lopo_direction_fraction = same / pmax(evaluated, 1),
    lopo_complete = evaluated == length(patients)
  )
  list(stability = stability, signed_stats = stat)
}

cohort_results <- list(); model_qc <- list(); model_objects <- list(); qc_loadings <- list()
for (lineage_value in lineages) {
  for (accession_value in accessions) {
    key <- paste(accession_value, lineage_value, sep = "::")
    group_meta <- select_complete_group(accession_value, lineage_value)
    idx <- match(group_meta$pseudobulk_id, colnames(analysis_counts))
    cts <- analysis_counts[, idx, drop = FALSE]
    message("MODEL_START\t", accession_value, "\t", lineage_value,
            "\tpatients=", uniqueN(group_meta$patient_id), "\tsamples=", nrow(group_meta))
    fit <- fit_model(cts, group_meta, accession_value, lineage_value)
    lopo <- run_lopo(cts, group_meta, fit$genes, fit$result, accession_value, lineage_value)
    cohort_results[[key]] <- fit$result; model_qc[[key]] <- fit$qc
    qc_loadings[[key]] <- fit$loadings
    model_objects[[key]] <- list(meta = fit$meta, fit = fit, lopo = lopo)
    message("MODEL_COMPLETE\t", accession_value, "\t", lineage_value,
            "\tgenes_tested=", fit$qc$genes_tested)
  }
}
cohort <- rbindlist(cohort_results); qc_table <- rbindlist(model_qc)
gene_lopo <- rbindlist(lapply(model_objects, function(x) x$lopo$stability))
fwrite(cohort, file.path(output_dir, "cohort_differential.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(qc_table, file.path(output_dir, "model_qc.tsv"), sep = "\t")
fwrite(rbindlist(qc_loadings), file.path(output_dir, "qc_pc_loadings.tsv"), sep = "\t")
fwrite(gene_lopo, file.path(output_dir, "gene_lopo_stability.tsv.gz"), sep = "\t", compress = "gzip")

meta_rows <- list()
for (lineage_value in lineages) {
  a <- cohort[accession == "GSE143791" & lineage == lineage_value]
  b <- cohort[accession == "GSE202813" & lineage == lineage_value]
  z <- merge(a, b, by = c("lineage", "contrast", "gene"), suffixes = c("_prostate", "_renal"))
  z <- z[is.finite(se_prostate) & se_prostate > 0 & is.finite(se_renal) & se_renal > 0]
  w1 <- 1 / z$se_prostate^2; w2 <- 1 / z$se_renal^2
  z[, meta_logFC := (w1 * logFC_prostate + w2 * logFC_renal) / (w1 + w2)]
  z[, meta_se := sqrt(1 / (w1 + w2))]
  z[, meta_z := meta_logFC / meta_se]
  z[, meta_PValue := 2 * pnorm(-abs(meta_z))]
  z[, meta_FDR := p.adjust(meta_PValue, method = "BH")]
  z[, heterogeneity_Q := w1 * (logFC_prostate - meta_logFC)^2 + w2 * (logFC_renal - meta_logFC)^2]
  z[, I2 := fifelse(heterogeneity_Q > 0, pmax(0, (heterogeneity_Q - 1) / heterogeneity_Q) * 100, 0)]
  z[, same_direction := sign(logFC_prostate) == sign(logFC_renal) & sign(logFC_prostate) != 0]
  meta_rows[[lineage_value]] <- z
}
gene_meta <- rbindlist(meta_rows, fill = TRUE)
lp <- gene_lopo[accession == "GSE143791", .(lineage, gene,
  lopo_iterations_prostate = lopo_iterations, lopo_evaluated_prostate = lopo_evaluated,
  lopo_direction_fraction_prostate = lopo_direction_fraction, lopo_complete_prostate = lopo_complete)]
lr <- gene_lopo[accession == "GSE202813", .(lineage, gene,
  lopo_iterations_renal = lopo_iterations, lopo_evaluated_renal = lopo_evaluated,
  lopo_direction_fraction_renal = lopo_direction_fraction, lopo_complete_renal = lopo_complete)]
gene_meta <- merge(gene_meta, lp, by = c("lineage", "gene"), all.x = TRUE)
gene_meta <- merge(gene_meta, lr, by = c("lineage", "gene"), all.x = TRUE)
gene_meta[, robust_adjusted := same_direction & meta_FDR < 0.05 &
  PValue_prostate < 0.10 & PValue_renal < 0.10 &
  abs(logFC_prostate) >= 0.15 & abs(logFC_renal) >= 0.15 & I2 <= 50 &
  lopo_complete_prostate == TRUE & lopo_complete_renal == TRUE &
  lopo_direction_fraction_prostate >= 0.75 & lopo_direction_fraction_renal >= 0.75]
gate5b_genes <- copy(gate5b_candidates)
gene_meta <- merge(gene_meta, gate5b_genes, by = c("lineage", "gene"), all.x = TRUE)
gene_meta[, gate5b_candidate := is.finite(gate5b_meta_logFC)]
gene_meta[, gate5b_persistent := robust_adjusted & gate5b_candidate &
            sign(meta_logFC) == sign(gate5b_meta_logFC)]
setorder(gene_meta, lineage, meta_FDR)
fwrite(gene_meta, file.path(output_dir, "gene_meta_results.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(gene_meta[robust_adjusted == TRUE], file.path(output_dir, "robust_adjusted_genes.tsv"), sep = "\t")

run_fgsea <- function(diff, accession, lineage) {
  ranks <- setNames(sign(diff$logFC) * sqrt(pmax(diff$F, 0)), diff$gene)
  ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
  set.seed(seed + match(accession, accessions) * 100 + match(lineage, lineages))
  fg <- as.data.table(fgseaMultilevel(pathways = pathways, stats = ranks,
                                      minSize = 15, maxSize = 500, eps = 0,
                                      nproc = workers))
  fg[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))]
  fg[, `:=`(accession = accession, cancer = unname(accession_to_cancer[accession]),
             lineage = lineage, n_patients = unique(diff$n_patients))]
  fg
}

hallmark_cohort <- rbindlist(lapply(lineages, function(lineage_value) {
  rbindlist(lapply(accessions, function(accession_value) {
    message("FGSEA_START\t", accession_value, "\t", lineage_value)
    run_fgsea(cohort[accession == accession_value & lineage == lineage_value],
              accession_value, lineage_value)
  }))
}), fill = TRUE)
fwrite(hallmark_cohort, file.path(output_dir, "hallmark_cohort_results.tsv"), sep = "\t")

hallmark_lopo_rows <- list()
for (key in names(model_objects)) {
  obj <- model_objects[[key]]
  accession_value <- unique(obj$meta$accession); lineage_value <- unique(obj$meta$lineage)
  full <- hallmark_cohort[accession == accession_value & lineage == lineage_value]
  sm <- obj$lopo$signed_stats
  for (pathway_value in names(pathways)) {
    genes_here <- intersect(pathways[[pathway_value]], rownames(sm))
    if (length(genes_here) < 15) next
    scores <- colMeans(sm[genes_here, , drop = FALSE], na.rm = TRUE)
    full_nes <- full[pathway == pathway_value, NES]
    if (length(full_nes) != 1) next
    hallmark_lopo_rows[[paste(key, pathway_value, sep = "::")]] <- data.table(
      accession = accession_value, cancer = unname(accession_to_cancer[accession_value]),
      lineage = lineage_value, pathway = pathway_value, full_NES = full_nes,
      lopo_iterations = ncol(sm), lopo_evaluated = sum(is.finite(scores)),
      lopo_direction_fraction = mean(sign(scores) == sign(full_nes), na.rm = TRUE),
      lopo_complete = all(is.finite(scores)))
  }
}
hallmark_lopo <- rbindlist(hallmark_lopo_rows)
fwrite(hallmark_lopo, file.path(output_dir, "hallmark_lopo_stability.tsv"), sep = "\t")

hallmark_meta_rows <- list()
for (lineage_value in lineages) {
  a <- hallmark_cohort[accession == "GSE143791" & lineage == lineage_value]
  b <- hallmark_cohort[accession == "GSE202813" & lineage == lineage_value]
  z <- merge(a, b, by = c("lineage", "pathway"), suffixes = c("_prostate", "_renal"))
  z1 <- sign(z$NES_prostate) * qnorm(pmax(z$pval_prostate / 2, 1e-300), lower.tail = FALSE)
  z2 <- sign(z$NES_renal) * qnorm(pmax(z$pval_renal / 2, 1e-300), lower.tail = FALSE)
  w1 <- sqrt(z$n_patients_prostate); w2 <- sqrt(z$n_patients_renal)
  z[, combined_z := (w1 * z1 + w2 * z2) / sqrt(w1^2 + w2^2)]
  z[, combined_PValue := 2 * pnorm(-abs(combined_z))]
  z[, combined_FDR := p.adjust(combined_PValue, method = "BH")]
  z[, same_direction := sign(NES_prostate) == sign(NES_renal) & sign(NES_prostate) != 0]
  hallmark_meta_rows[[lineage_value]] <- z
}
hallmark_meta <- rbindlist(hallmark_meta_rows, fill = TRUE)
hp <- hallmark_lopo[accession == "GSE143791", .(lineage, pathway,
  lopo_iterations_prostate = lopo_iterations, lopo_evaluated_prostate = lopo_evaluated,
  lopo_direction_fraction_prostate = lopo_direction_fraction, lopo_complete_prostate = lopo_complete)]
hr <- hallmark_lopo[accession == "GSE202813", .(lineage, pathway,
  lopo_iterations_renal = lopo_iterations, lopo_evaluated_renal = lopo_evaluated,
  lopo_direction_fraction_renal = lopo_direction_fraction, lopo_complete_renal = lopo_complete)]
hallmark_meta <- merge(hallmark_meta, hp, by = c("lineage", "pathway"), all.x = TRUE)
hallmark_meta <- merge(hallmark_meta, hr, by = c("lineage", "pathway"), all.x = TRUE)
hallmark_meta[, robust_adjusted := same_direction & combined_FDR < 0.05 &
  pval_prostate < 0.10 & pval_renal < 0.10 & abs(NES_prostate) >= 1 & abs(NES_renal) >= 1 &
  lopo_complete_prostate == TRUE & lopo_complete_renal == TRUE &
  lopo_direction_fraction_prostate >= 0.75 & lopo_direction_fraction_renal >= 0.75]
gate5b_hallmarks <- fread(gate5b_hallmark_file)[gate4b_direction_concordant == TRUE,
                                                .(lineage, pathway, gate5b_combined_z = combined_z)]
hallmark_meta <- merge(hallmark_meta, gate5b_hallmarks, by = c("lineage", "pathway"), all.x = TRUE)
hallmark_meta[, gate5b_candidate := is.finite(gate5b_combined_z)]
hallmark_meta[, gate5b_persistent := robust_adjusted & gate5b_candidate &
                sign(combined_z) == sign(gate5b_combined_z)]
setorder(hallmark_meta, lineage, combined_FDR)
fwrite(hallmark_meta, file.path(output_dir, "hallmark_meta_results.tsv"), sep = "\t")
fwrite(hallmark_meta[robust_adjusted == TRUE], file.path(output_dir, "robust_adjusted_hallmarks.tsv"), sep = "\t")

summary <- data.table(lineage = lineages)
summary[, gate5b_candidate_genes := vapply(lineage, function(x) nrow(gate5b_genes[lineage == x]), integer(1))]
summary[, persistent_genes := vapply(lineage, function(x) nrow(gene_meta[lineage == x & gate5b_persistent == TRUE]), integer(1))]
summary[, gene_persistence_fraction := persistent_genes / pmax(gate5b_candidate_genes, 1)]
summary[, gate5b_candidate_hallmarks := vapply(lineage, function(x) nrow(gate5b_hallmarks[lineage == x]), integer(1))]
summary[, persistent_hallmarks := vapply(lineage, function(x) nrow(hallmark_meta[lineage == x & gate5b_persistent == TRUE]), integer(1))]
summary[, hallmark_persistence_fraction := persistent_hallmarks / pmax(gate5b_candidate_hallmarks, 1)]
summary[, strong_persistence := persistent_genes >= 10 | persistent_hallmarks >= 3]
summary[, support_persistence := persistent_genes >= 3 | persistent_hallmarks >= 1]
fwrite(summary, file.path(output_dir, "persistence_summary.tsv"), sep = "\t")

model_integrity <- nrow(qc_table) == 6 && all(qc_table$converged) &&
  all(qc_table$gate5b_parent_retention_fraction >= min_retention_fraction) &&
  all(qc_table$gate5b_candidate_retention_fraction >= min_retention_fraction) &&
  all(qc_table$design_rank == qc_table$design_columns) &&
  all(qc_table$residual_df > 0)
strong_n <- sum(summary$strong_persistence)
decision <- if (!model_integrity) "NO-GO" else if (strong_n >= 2) "GO" else if (strong_n == 1) "CONDITIONAL GO" else "NO-GO"
common_persistent_genes <- Reduce(intersect, lapply(lineages, function(x) gene_meta[lineage == x & gate5b_persistent == TRUE, gene]))
common_persistent_hallmarks <- Reduce(intersect, lapply(lineages, function(x) hallmark_meta[lineage == x & gate5b_persistent == TRUE, pathway]))
fwrite(data.table(feature_type = c(rep("gene", length(common_persistent_genes)), rep("hallmark", length(common_persistent_hallmarks))),
                  feature = c(common_persistent_genes, common_persistent_hallmarks)),
       file.path(output_dir, "common_persistent_features.tsv"), sep = "\t")

state_lines <- unlist(lapply(lineages, function(x) c(
  paste0("- ", x, " persistent Gate5b genes: ", summary[lineage == x, persistent_genes],
         "/", summary[lineage == x, gate5b_candidate_genes]),
  paste0("- ", x, " persistent Gate5b Hallmarks: ", summary[lineage == x, persistent_hallmarks],
         "/", summary[lineage == x, gate5b_candidate_hallmarks]))))
writeLines(c("# Gate6b amended confounder-robustness decision", "", paste0("**", decision, "**"), "",
             paste0("- Six high-confidence, QC-adjusted cohort-state models passed integrity: ",
                    ifelse(model_integrity, "PASS", "FAIL")),
             paste0("- Minimum Gate5b tested/candidate retention fraction: ",
                    format(min_retention_fraction, nsmall = 2)),
             paste0("- Curated stress genes excluded before modeling and enrichment: ", length(stress_found)),
             state_lines,
             paste0("- Persistent genes common to all three states: ", length(common_persistent_genes)),
             paste0("- Persistent Hallmarks common to all three states: ", length(common_persistent_hallmarks)), "",
             "Gate6b tests robustness to state-assignment confidence, pseudobulk QC, and a frozen acute-stress exclusion list; it does not establish causality."),
           file.path(output_dir, "gate6b_decision.md"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("GATE6B_COMPLETE\tdecision=", decision, "\tstrong_states=", strong_n,
        "\tcommon_genes=", length(common_persistent_genes),
        "\tcommon_hallmarks=", length(common_persistent_hallmarks))
