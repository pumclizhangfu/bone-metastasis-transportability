#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(edgeR)
  library(fgsea)
  library(msigdbr)
  library(BiocParallel)
})

parse_args <- function(x) {
  if (length(x) %% 2 != 0) stop("Arguments must be provided as --key value pairs")
  out <- list()
  for (i in seq(1, length(x), by = 2)) {
    out[[sub("^--", "", x[[i]])]] <- x[[i + 1]]
  }
  out
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("input", "annotations", "cv-predictions", "gate4-genes",
              "gate4-hallmarks", "output", "workers", "seed", "min-cells")
missing_args <- setdiff(required, names(cfg))
if (length(missing_args)) stop("Missing arguments: ", paste(missing_args, collapse = ", "))

input_file <- normalizePath(cfg$input, mustWork = TRUE)
annotation_file <- normalizePath(cfg$annotations, mustWork = TRUE)
cv_prediction_file <- normalizePath(cfg[["cv-predictions"]], mustWork = TRUE)
gate4_gene_file <- normalizePath(cfg[["gate4-genes"]], mustWork = TRUE)
gate4_hallmark_file <- normalizePath(cfg[["gate4-hallmarks"]], mustWork = TRUE)
output_dir <- normalizePath(cfg$output, mustWork = FALSE)
workers <- as.integer(cfg$workers)
seed <- as.integer(cfg$seed)
min_cells <- as.integer(cfg[["min-cells"]])
if (is.na(workers) || workers < 1 || workers > 20) stop("workers must be 1..20")
if (is.na(seed)) stop("seed must be an integer")
if (is.na(min_cells) || min_cells < 10) stop("min-cells must be at least 10")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(seed)
bp <- if (.Platform$OS.type == "unix" && workers > 1) {
  MulticoreParam(workers = workers, progressbar = FALSE)
} else {
  SerialParam(progressbar = FALSE)
}

message("CONFIG\tinput=", input_file, "\tannotations=", annotation_file,
        "\tcv_predictions=", cv_prediction_file,
        "\tgate4_genes=", gate4_gene_file, "\tgate4_hallmarks=", gate4_hallmark_file,
        "\toutput=", output_dir,
        "\tworkers=", workers, "\tseed=", seed, "\tmin_cells=", min_cells)

pb <- readRDS(input_file)
if (!all(c("counts", "metadata") %in% names(pb))) stop("Input RDS lacks counts or metadata")
all_counts <- pb$counts
meta <- as.data.table(pb$metadata)
if (!inherits(all_counts, "Matrix")) all_counts <- as(all_counts, "dgCMatrix")
if (!identical(colnames(all_counts), meta$pseudobulk_id)) stop("Count columns and metadata are not aligned")
if (anyDuplicated(rownames(all_counts))) stop("Gene identifiers are not unique")
if (!all(meta$aggregation_exact)) stop("Input contains failed aggregation checks")

gate4_robust_genes <- fread(gate4_gene_file)[lineage == "T_NK"]
gate4_robust_hallmarks <- fread(gate4_hallmark_file)[lineage == "T_NK"]
if (!all(c("gene", "meta_logFC") %in% names(gate4_robust_genes))) {
  stop("Gate4b robust-gene table lacks gene/meta_logFC")
}
if (!all(c("pathway", "combined_z") %in% names(gate4_robust_hallmarks))) {
  stop("Gate4b robust-Hallmark table lacks pathway/combined_z")
}
gate4_gene_direction <- setNames(sign(gate4_robust_genes$meta_logFC), gate4_robust_genes$gene)
gate4_hallmark_direction <- setNames(sign(gate4_robust_hallmarks$combined_z),
                                     gate4_robust_hallmarks$pathway)

# The internal `lineage` result column stores the validated T/NK meta-state.
lineages <- c("CD4_T", "CD8_CTL", "NK_NKT")
accessions <- c("GSE143791", "GSE202813")
accession_to_cancer <- c(GSE143791 = "prostate", GSE202813 = "renal")
contrasts <- c("trend", "involved_distal", "tumor_involved", "tumor_distal")

msig <- as.data.table(msigdbr(db_species = "HS", species = "human", collection = "H"))
pathways <- split(msig$gene_symbol, msig$gs_name)
pathways <- lapply(pathways, unique)
if (length(pathways) != 50) stop("Expected 50 MSigDB Hallmark sets; found ", length(pathways))
fwrite(unique(msig[, .(gs_name, gene_symbol, db_version)]),
       file.path(output_dir, "hallmark_gene_sets.tsv"), sep = "\t")
message("PATHWAYS_READY\tsets=", length(pathways), "\tdb_version=", unique(msig$db_version))

# Patient-paired state-composition analysis. The primary analysis uses all
# validated state assignments. A prostate-only sensitivity set removes cells
# below a state-specific delta threshold learned from correct RCC LOPO calls.
cv_predictions <- fread(cv_prediction_file)
state_annotations <- fread(annotation_file)
cv_tnk <- cv_predictions[lineage == "T_NK" & pruned_state %in% lineages]
confidence_thresholds <- cv_tnk[pruned_state == true_state, .(
  delta_threshold = as.numeric(quantile(delta_next, probs = 0.20, na.rm = TRUE, type = 8)),
  correct_cv_cells = .N
), by = .(meta_state = pruned_state)]
if (!setequal(confidence_thresholds$meta_state, lineages)) {
  stop("Could not derive a high-confidence threshold for every T/NK state")
}
fwrite(confidence_thresholds, file.path(output_dir, "composition_confidence_thresholds.tsv"), sep = "\t")

composition_primary_counts <- meta[lineage == "T_NK" & meta_state %in% lineages,
                                   .(accession, cancer, sample_id, patient_id,
                                     compartment, meta_state, n_cells)]
valid_compartments <- composition_primary_counts[, .(
  n_states = uniqueN(meta_state),
  all_states_ge_min = uniqueN(meta_state) == length(lineages) & all(n_cells >= min_cells)
), by = .(accession, cancer, patient_id, sample_id, compartment)]
composition_patients <- valid_compartments[all_states_ge_min == TRUE, .(
  valid_compartments = .N,
  has_all = all(c("distal", "involved", "tumor") %in% compartment)
), by = .(accession, cancer, patient_id)][valid_compartments == 3 & has_all == TRUE]
if (any(composition_patients[, .N, by = accession]$N < 3)) {
  stop("Fewer than three complete composition patients in a cohort")
}
composition_design <- merge(
  valid_compartments[, .(accession, cancer, patient_id, sample_id, compartment)],
  composition_patients[, .(accession, patient_id)],
  by = c("accession", "patient_id"), all = FALSE
)
composition_design[, stage := c(distal = 0, involved = 1, tumor = 2)[compartment]]
setorder(composition_design, accession, patient_id, stage)

threshold_lookup <- setNames(confidence_thresholds$delta_threshold, confidence_thresholds$meta_state)
high_confidence_annotations <- state_annotations[
  lineage == "T_NK" & meta_state %in% lineages &
    (accession == "GSE202813" |
       (accession == "GSE143791" & is.finite(delta_next) &
          delta_next >= threshold_lookup[meta_state]))
]
composition_sensitivity_counts <- high_confidence_annotations[, .(
  n_cells = .N
), by = .(accession, cancer, sample_id, patient_id, compartment, meta_state)]

fit_composition <- function(count_table, analysis_set) {
  rows <- list()
  counter <- 0L
  for (accession_value in accessions) {
    design <- composition_design[accession == accession_value]
    for (state_value in lineages) {
      counter <- counter + 1L
      counts_here <- count_table[accession == accession_value,
                                 .(sample_id, meta_state, n_cells)]
      design_grid <- copy(design)[, composition_join_key := 1L]
      state_grid <- data.table(meta_state = lineages, composition_join_key = 1L)
      grid <- merge(design_grid, state_grid, by = "composition_join_key",
                    allow.cartesian = TRUE)[, composition_join_key := NULL]
      if (!nrow(grid)) stop("Composition grid construction failed")
      grid <- merge(grid, counts_here, by = c("sample_id", "meta_state"), all.x = TRUE)
      grid[is.na(n_cells), n_cells := 0L]
      totals <- grid[, .(total_cells = sum(n_cells)),
                     by = .(accession, cancer, patient_id, sample_id, compartment, stage)]
      target <- grid[meta_state == state_value,
                     .(sample_id, state_cells = n_cells)]
      z <- merge(totals, target, by = "sample_id", all.x = TRUE)
      z[, response := log((state_cells + 0.5) / (total_cells - state_cells + 0.5))]
      patient <- factor(z$patient_id)
      stage <- as.numeric(z$stage)
      design_matrix <- model.matrix(~ patient + stage)
      if (qr(design_matrix)$rank != ncol(design_matrix)) stop("Composition design rank failure")
      fit <- lm.fit(design_matrix, z$response)
      residual_df <- nrow(design_matrix) - fit$rank
      sigma2 <- sum(fit$residuals^2) / residual_df
      covariance <- sigma2 * chol2inv(qr.R(fit$qr))
      coef_index <- match("stage", colnames(design_matrix))
      estimate <- unname(fit$coefficients[[coef_index]])
      se <- sqrt(covariance[coef_index, coef_index])
      statistic <- estimate / se
      p_value <- 2 * pt(-abs(statistic), df = residual_df)

      patients <- sort(unique(z$patient_id))
      lopo_estimates <- vapply(patients, function(left_out) {
        keep <- z$patient_id != left_out
        p2 <- factor(z$patient_id[keep])
        s2 <- as.numeric(z$stage[keep])
        d2 <- model.matrix(~ p2 + s2)
        if (qr(d2)$rank != ncol(d2)) return(NA_real_)
        f2 <- lm.fit(d2, z$response[keep])
        unname(f2$coefficients[[match("s2", colnames(d2))]])
      }, numeric(1))
      rows[[counter]] <- data.table(
        analysis_set = analysis_set,
        accession = accession_value,
        cancer = unname(accession_to_cancer[accession_value]),
        meta_state = state_value,
        n_patients = uniqueN(z$patient_id),
        n_samples = nrow(z),
        min_total_cells = min(z$total_cells),
        log_odds_per_stage = estimate,
        se = se,
        statistic = statistic,
        residual_df = residual_df,
        PValue = p_value,
        lopo_iterations = length(patients),
        lopo_evaluated = sum(is.finite(lopo_estimates)),
        lopo_direction_fraction = mean(sign(lopo_estimates) == sign(estimate), na.rm = TRUE),
        lopo_complete = all(is.finite(lopo_estimates))
      )
    }
  }
  rbindlist(rows, use.names = TRUE)
}

composition_primary <- fit_composition(composition_primary_counts, "primary")
composition_sensitivity <- fit_composition(composition_sensitivity_counts, "high_confidence")
composition_cohort <- rbindlist(list(composition_primary, composition_sensitivity), use.names = TRUE)
fwrite(composition_cohort, file.path(output_dir, "composition_cohort_results.tsv"), sep = "\t")

cp <- composition_primary[accession == "GSE143791"]
cr <- composition_primary[accession == "GSE202813"]
composition_meta <- merge(cp, cr, by = "meta_state", suffixes = c("_prostate", "_renal"))
w1 <- 1 / composition_meta$se_prostate^2
w2 <- 1 / composition_meta$se_renal^2
composition_meta[, meta_log_odds := (w1 * log_odds_per_stage_prostate +
                                       w2 * log_odds_per_stage_renal) / (w1 + w2)]
composition_meta[, meta_se := sqrt(1 / (w1 + w2))]
composition_meta[, meta_z := meta_log_odds / meta_se]
composition_meta[, meta_PValue := 2 * pnorm(-abs(meta_z))]
composition_meta[, meta_FDR := p.adjust(meta_PValue, method = "BH")]
composition_meta[, heterogeneity_Q := w1 * (log_odds_per_stage_prostate - meta_log_odds)^2 +
                                           w2 * (log_odds_per_stage_renal - meta_log_odds)^2]
composition_meta[, I2 := fifelse(heterogeneity_Q > 0,
                                 pmax(0, (heterogeneity_Q - 1) / heterogeneity_Q) * 100, 0)]
composition_meta[, same_direction := sign(log_odds_per_stage_prostate) ==
                                          sign(log_odds_per_stage_renal) &
                                          sign(log_odds_per_stage_prostate) != 0]
sens_p <- composition_sensitivity[accession == "GSE143791",
                                  .(meta_state,
                                    sensitivity_log_odds_prostate = log_odds_per_stage,
                                    sensitivity_PValue_prostate = PValue)]
composition_meta <- merge(composition_meta, sens_p, by = "meta_state", all.x = TRUE)
composition_meta[, sensitivity_direction_stable :=
                   sign(sensitivity_log_odds_prostate) == sign(log_odds_per_stage_prostate)]
composition_meta[, robust_composition := same_direction & meta_FDR < 0.05 &
                   PValue_prostate < 0.10 & PValue_renal < 0.10 &
                   abs(log_odds_per_stage_prostate) >= log(1.25) &
                   abs(log_odds_per_stage_renal) >= log(1.25) & I2 <= 50 &
                   lopo_complete_prostate == TRUE & lopo_complete_renal == TRUE &
                   lopo_direction_fraction_prostate >= 0.75 &
                   lopo_direction_fraction_renal >= 0.75 &
                   sensitivity_direction_stable == TRUE]
composition_meta[, abs_meta_log_odds_sort := abs(meta_log_odds)]
setorder(composition_meta, meta_FDR, -abs_meta_log_odds_sort)
composition_meta[, abs_meta_log_odds_sort := NULL]
fwrite(composition_meta, file.path(output_dir, "composition_meta_results.tsv"), sep = "\t")
fwrite(composition_meta[robust_composition == TRUE],
       file.path(output_dir, "robust_composition_shifts.tsv"), sep = "\t")

select_complete_group <- function(accession_value, lineage_value) {
  z <- meta[accession == accession_value & lineage == "T_NK" &
              meta_state == lineage_value & n_cells >= min_cells]
  complete <- z[, .(
    n_rows = .N,
    has_all = all(c("distal", "involved", "tumor") %in% compartment)
  ), by = patient_id][n_rows == 3 & has_all == TRUE, patient_id]
  z <- z[patient_id %in% complete]
  z[, compartment := factor(compartment, levels = c("distal", "involved", "tumor"))]
  z[, stage := c(distal = 0, involved = 1, tumor = 2)[as.character(compartment)]]
  setorder(z, patient_id, compartment)
  if (uniqueN(z$patient_id) < 3) stop("Fewer than 3 complete patients for ", accession_value, " ", lineage_value)
  if (nrow(z) != 3 * uniqueN(z$patient_id)) stop("Unbalanced triplets for ", accession_value, " ", lineage_value)
  z
}

extract_test <- function(test, accession, lineage, contrast, n_patients, n_samples) {
  tab <- as.data.table(topTags(test, n = Inf, sort.by = "none")$table, keep.rownames = "gene")
  tab[, se := fifelse(is.finite(F) & F > 0, abs(logFC) / sqrt(F), NA_real_)]
  tab[, `:=`(
    accession = accession,
    cancer = unname(accession_to_cancer[accession]),
    lineage = lineage,
    contrast = contrast,
    n_patients = n_patients,
    n_pseudobulks = n_samples
  )]
  setcolorder(tab, c("accession", "cancer", "lineage", "contrast", "gene",
                     "logFC", "se", "logCPM", "F", "PValue", "FDR",
                     "n_patients", "n_pseudobulks"))
  tab
}

fit_full_model <- function(counts_matrix, group_meta, accession, lineage) {
  patient <- factor(group_meta$patient_id)
  compartment <- factor(group_meta$compartment, levels = c("distal", "involved", "tumor"))
  stage <- as.numeric(group_meta$stage)
  design_cat <- model.matrix(~ patient + compartment)
  design_trend <- model.matrix(~ patient + stage)
  if (qr(design_cat)$rank != ncol(design_cat)) stop("Categorical design is rank deficient")
  if (qr(design_trend)$rank != ncol(design_trend)) stop("Trend design is rank deficient")

  y0 <- DGEList(counts = counts_matrix)
  keep <- filterByExpr(y0, design = design_cat)
  if (sum(keep) < 5000) stop("Fewer than 5,000 genes pass filtering for ", accession, " ", lineage)
  y <- y0[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")

  y_cat <- estimateDisp(y, design_cat, robust = TRUE)
  fit_cat <- glmQLFit(y_cat, design_cat, robust = TRUE, abundance.trend = TRUE)
  y_trend <- estimateDisp(y, design_trend, robust = TRUE)
  fit_trend <- glmQLFit(y_trend, design_trend, robust = TRUE, abundance.trend = TRUE)

  coef_involved <- match("compartmentinvolved", colnames(design_cat))
  coef_tumor <- match("compartmenttumor", colnames(design_cat))
  coef_stage <- match("stage", colnames(design_trend))
  if (anyNA(c(coef_involved, coef_tumor, coef_stage))) stop("Expected model coefficients were not found")
  tumor_involved <- numeric(ncol(design_cat))
  tumor_involved[coef_tumor] <- 1
  tumor_involved[coef_involved] <- -1

  n_patients <- uniqueN(group_meta$patient_id)
  n_samples <- nrow(group_meta)
  results <- rbindlist(list(
    extract_test(glmQLFTest(fit_trend, coef = coef_stage), accession, lineage, "trend", n_patients, n_samples),
    extract_test(glmQLFTest(fit_cat, coef = coef_involved), accession, lineage, "involved_distal", n_patients, n_samples),
    extract_test(glmQLFTest(fit_cat, contrast = tumor_involved), accession, lineage, "tumor_involved", n_patients, n_samples),
    extract_test(glmQLFTest(fit_cat, coef = coef_tumor), accession, lineage, "tumor_distal", n_patients, n_samples)
  ))

  qc <- data.table(
    accession = accession,
    cancer = unname(accession_to_cancer[accession]),
    lineage = lineage,
    n_patients = n_patients,
    n_pseudobulks = n_samples,
    min_cells = min(group_meta$n_cells),
    median_cells = median(group_meta$n_cells),
    min_library_size = min(colSums(counts_matrix)),
    median_library_size = median(colSums(counts_matrix)),
    genes_input = nrow(counts_matrix),
    genes_tested = sum(keep),
    trend_design_rank = qr(design_trend)$rank,
    trend_design_columns = ncol(design_trend),
    categorical_design_rank = qr(design_cat)$rank,
    categorical_design_columns = ncol(design_cat),
    converged = TRUE
  )
  list(results = results, qc = qc, genes = rownames(y), trend = results[contrast == "trend"])
}

run_lopo <- function(counts_matrix, group_meta, genes, full_trend, accession, lineage) {
  patients <- sort(unique(group_meta$patient_id))
  lfc <- matrix(NA_real_, nrow = length(genes), ncol = length(patients),
                dimnames = list(genes, patients))
  stat <- lfc
  for (i in seq_along(patients)) {
    left_out <- patients[[i]]
    use <- group_meta$patient_id != left_out
    m <- group_meta[use]
    cts <- counts_matrix[genes, use, drop = FALSE]
    patient <- factor(m$patient_id)
    stage <- as.numeric(m$stage)
    design <- model.matrix(~ patient + stage)
    if (qr(design)$rank != ncol(design)) stop("LOPO design rank deficiency: ", accession, " ", lineage, " ", left_out)
    nonzero <- Matrix::rowSums(cts) > 0
    y <- DGEList(counts = cts[nonzero, , drop = FALSE])
    y <- calcNormFactors(y, method = "TMM")
    y <- estimateDisp(y, design, robust = TRUE)
    fit <- glmQLFit(y, design, robust = TRUE, abundance.trend = TRUE)
    test <- glmQLFTest(fit, coef = match("stage", colnames(design)))
    tab <- test$table
    lfc[rownames(tab), i] <- tab$logFC
    stat[rownames(tab), i] <- sign(tab$logFC) * sqrt(pmax(tab$F, 0))
    message("LOPO_COMPLETE\t", accession, "\t", lineage, "\tleft_out=", left_out,
            "\tgenes=", nrow(tab))
  }
  full_sign <- sign(full_trend$logFC[match(genes, full_trend$gene)])
  evaluated <- rowSums(is.finite(lfc))
  same <- rowSums(sign(lfc) == full_sign, na.rm = TRUE)
  stability <- data.table(
    accession = accession,
    cancer = unname(accession_to_cancer[accession]),
    lineage = lineage,
    gene = genes,
    full_logFC = full_trend$logFC[match(genes, full_trend$gene)],
    lopo_iterations = length(patients),
    lopo_evaluated = evaluated,
    lopo_same_direction = same,
    lopo_direction_fraction = same / pmax(evaluated, 1),
    lopo_complete = evaluated == length(patients),
    lopo_min_logFC = apply(lfc, 1, function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)),
    lopo_max_logFC = apply(lfc, 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
  )
  list(stability = stability, signed_stats = stat)
}

run_fgsea <- function(diff, accession, lineage, contrast, run_seed) {
  ranks <- sign(diff$logFC) * sqrt(pmax(diff$F, 0))
  names(ranks) <- diff$gene
  ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
  set.seed(run_seed)
  fg <- as.data.table(fgseaMultilevel(
    pathways = pathways,
    stats = ranks,
    minSize = 15,
    maxSize = 500,
    eps = 0,
    nproc = workers
  ))
  fg[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))]
  fg[, `:=`(
    accession = accession,
    cancer = unname(accession_to_cancer[accession]),
    lineage = lineage,
    contrast = contrast,
    n_patients = unique(diff$n_patients)
  )]
  preferred <- c("accession", "cancer", "lineage", "contrast", "pathway",
                 "NES", "pval", "padj", "log2err", "size", "ES",
                 "n_patients", "leadingEdge")
  setcolorder(fg, intersect(preferred, names(fg)))
  fg
}

cohort_results <- list()
model_qc <- list()
group_objects <- list()
counter <- 0L
for (lineage in lineages) {
  for (accession in accessions) {
    counter <- counter + 1L
    key <- paste(accession, lineage, sep = "::")
    group_meta <- select_complete_group(accession, lineage)
    idx <- match(group_meta$pseudobulk_id, colnames(all_counts))
    counts_matrix <- all_counts[, idx, drop = FALSE]
    message("MODEL_START\t", accession, "\t", lineage,
            "\tpatients=", uniqueN(group_meta$patient_id), "\tsamples=", nrow(group_meta))
    fit <- fit_full_model(counts_matrix, group_meta, accession, lineage)
    cohort_results[[key]] <- fit$results
    model_qc[[key]] <- fit$qc
    lopo <- run_lopo(counts_matrix, group_meta, fit$genes, fit$trend, accession, lineage)
    group_objects[[key]] <- list(meta = group_meta, fit = fit, lopo = lopo)
    message("MODEL_COMPLETE\t", accession, "\t", lineage,
            "\tgenes_tested=", fit$qc$genes_tested)
  }
}
cohort_differential <- rbindlist(cohort_results, use.names = TRUE)
model_qc_table <- rbindlist(model_qc, use.names = TRUE)
gene_lopo <- rbindlist(lapply(group_objects, function(x) x$lopo$stability), use.names = TRUE)
fwrite(cohort_differential, file.path(output_dir, "cohort_differential.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(model_qc_table, file.path(output_dir, "model_qc.tsv"), sep = "\t")
fwrite(gene_lopo, file.path(output_dir, "gene_lopo_stability.tsv.gz"), sep = "\t", compress = "gzip")

meta_gene_rows <- list()
for (lineage in lineages) {
  for (contrast in contrasts) {
    lineage_value <- lineage
    contrast_value <- contrast
    a <- cohort_differential[accession == "GSE143791" & lineage == lineage_value & contrast == contrast_value]
    b <- cohort_differential[accession == "GSE202813" & lineage == lineage_value & contrast == contrast_value]
    z <- merge(a, b, by = c("lineage", "contrast", "gene"), suffixes = c("_prostate", "_renal"))
    valid <- is.finite(z$se_prostate) & z$se_prostate > 0 & is.finite(z$se_renal) & z$se_renal > 0
    z <- z[valid]
    w1 <- 1 / z$se_prostate^2
    w2 <- 1 / z$se_renal^2
    z[, meta_logFC := (w1 * logFC_prostate + w2 * logFC_renal) / (w1 + w2)]
    z[, meta_se := sqrt(1 / (w1 + w2))]
    z[, meta_z := meta_logFC / meta_se]
    z[, meta_PValue := 2 * pnorm(-abs(meta_z))]
    z[, heterogeneity_Q := w1 * (logFC_prostate - meta_logFC)^2 + w2 * (logFC_renal - meta_logFC)^2]
    z[, heterogeneity_PValue := pchisq(heterogeneity_Q, df = 1, lower.tail = FALSE)]
    z[, I2 := fifelse(heterogeneity_Q > 0, pmax(0, (heterogeneity_Q - 1) / heterogeneity_Q) * 100, 0)]
    z[, same_direction := sign(logFC_prostate) == sign(logFC_renal) & sign(logFC_prostate) != 0]
    z[, meta_FDR := p.adjust(meta_PValue, method = "BH")]
    meta_gene_rows[[paste(lineage, contrast, sep = "::")]] <- z[, .(
      lineage, contrast, gene,
      logFC_prostate, se_prostate, PValue_prostate, FDR_prostate,
      logFC_renal, se_renal, PValue_renal, FDR_renal,
      meta_logFC, meta_se, meta_z, meta_PValue, meta_FDR,
      same_direction, heterogeneity_Q, heterogeneity_PValue, I2,
      n_patients_prostate, n_patients_renal
    )]
  }
}
gene_meta <- rbindlist(meta_gene_rows, use.names = TRUE)

lopo_p <- gene_lopo[accession == "GSE143791"]
lopo_r <- gene_lopo[accession == "GSE202813"]
setnames(lopo_p, c("lopo_iterations", "lopo_evaluated", "lopo_direction_fraction", "lopo_complete"),
         c("lopo_iterations_prostate", "lopo_evaluated_prostate", "lopo_direction_fraction_prostate", "lopo_complete_prostate"))
setnames(lopo_r, c("lopo_iterations", "lopo_evaluated", "lopo_direction_fraction", "lopo_complete"),
         c("lopo_iterations_renal", "lopo_evaluated_renal", "lopo_direction_fraction_renal", "lopo_complete_renal"))
gene_meta <- merge(gene_meta, lopo_p[, .(lineage, gene, lopo_iterations_prostate,
                                         lopo_evaluated_prostate, lopo_direction_fraction_prostate,
                                         lopo_complete_prostate)],
                   by = c("lineage", "gene"), all.x = TRUE)
gene_meta <- merge(gene_meta, lopo_r[, .(lineage, gene, lopo_iterations_renal,
                                         lopo_evaluated_renal, lopo_direction_fraction_renal,
                                         lopo_complete_renal)],
                   by = c("lineage", "gene"), all.x = TRUE)
gene_meta[, robust_conserved := contrast == "trend" & same_direction & meta_FDR < 0.05 &
            PValue_prostate < 0.10 & PValue_renal < 0.10 &
            abs(logFC_prostate) >= 0.15 & abs(logFC_renal) >= 0.15 & I2 <= 50 &
            lopo_complete_prostate == TRUE & lopo_complete_renal == TRUE &
            lopo_direction_fraction_prostate >= 0.75 & lopo_direction_fraction_renal >= 0.75]
gene_meta[, abs_meta_logFC_sort := abs(meta_logFC)]
setorder(gene_meta, lineage, contrast, meta_FDR, -abs_meta_logFC_sort)
gene_meta[, abs_meta_logFC_sort := NULL]
fwrite(gene_meta, file.path(output_dir, "gene_meta_results.tsv.gz"), sep = "\t", compress = "gzip")
robust_genes <- gene_meta[robust_conserved == TRUE]
robust_genes[, gate4b_tnk_overlap := gene %in% names(gate4_gene_direction)]
robust_genes[, gate4b_direction_concordant := gate4b_tnk_overlap &
               sign(meta_logFC) == gate4_gene_direction[gene]]
fwrite(robust_genes, file.path(output_dir, "robust_conserved_genes.tsv"), sep = "\t")

hallmark_cohort_rows <- list()
path_counter <- 0L
for (lineage in lineages) {
  for (accession in accessions) {
    for (contrast in c("trend", "tumor_distal")) {
      path_counter <- path_counter + 1L
      accession_value <- accession
      lineage_value <- lineage
      contrast_value <- contrast
      diff <- cohort_differential[accession == accession_value & lineage == lineage_value & contrast == contrast_value]
      hallmark_cohort_rows[[paste(accession, lineage, contrast, sep = "::")]] <-
        run_fgsea(diff, accession, lineage, contrast, seed + 1000L + path_counter)
      message("FGSEA_COMPLETE\t", accession, "\t", lineage, "\t", contrast)
    }
  }
}
hallmark_cohort <- rbindlist(hallmark_cohort_rows, use.names = TRUE, fill = TRUE)
fwrite(hallmark_cohort, file.path(output_dir, "hallmark_cohort_results.tsv"), sep = "\t")

# Deterministic LOPO pathway-direction audit using the mean signed QL statistic
# among genes belonging to each Hallmark set.
hallmark_lopo_rows <- list()
for (key in names(group_objects)) {
  obj <- group_objects[[key]]
  accession <- unique(obj$meta$accession)
  lineage <- unique(obj$meta$meta_state)
  accession_value <- accession
  lineage_value <- lineage
  full_path <- hallmark_cohort[accession == accession_value & lineage == lineage_value & contrast == "trend"]
  stat_matrix <- obj$lopo$signed_stats
  patients <- colnames(stat_matrix)
  for (pathway in names(pathways)) {
    genes <- intersect(pathways[[pathway]], rownames(stat_matrix))
    if (length(genes) < 15) next
    scores <- colMeans(stat_matrix[genes, , drop = FALSE], na.rm = TRUE)
    pathway_value <- pathway
    full_nes <- full_path[pathway == pathway_value, NES]
    if (length(full_nes) != 1) next
    same <- sign(scores) == sign(full_nes)
    hallmark_lopo_rows[[paste(key, pathway, sep = "::")]] <- data.table(
      accession = accession,
      cancer = unname(accession_to_cancer[accession]),
      lineage = lineage,
      pathway = pathway,
      full_NES = full_nes,
      lopo_iterations = length(patients),
      lopo_evaluated = sum(is.finite(scores)),
      lopo_same_direction = sum(same, na.rm = TRUE),
      lopo_direction_fraction = mean(same, na.rm = TRUE),
      lopo_complete = all(is.finite(scores)),
      lopo_min_score = min(scores, na.rm = TRUE),
      lopo_max_score = max(scores, na.rm = TRUE)
    )
  }
}
hallmark_lopo <- rbindlist(hallmark_lopo_rows, use.names = TRUE)
fwrite(hallmark_lopo, file.path(output_dir, "hallmark_lopo_stability.tsv"), sep = "\t")

hallmark_meta_rows <- list()
for (lineage in lineages) {
  for (contrast in c("trend", "tumor_distal")) {
    lineage_value <- lineage
    contrast_value <- contrast
    a <- hallmark_cohort[accession == "GSE143791" & lineage == lineage_value & contrast == contrast_value]
    b <- hallmark_cohort[accession == "GSE202813" & lineage == lineage_value & contrast == contrast_value]
    z <- merge(a, b, by = c("lineage", "contrast", "pathway"), suffixes = c("_prostate", "_renal"))
    z1 <- sign(z$NES_prostate) * qnorm(pmax(z$pval_prostate / 2, 1e-300), lower.tail = FALSE)
    z2 <- sign(z$NES_renal) * qnorm(pmax(z$pval_renal / 2, 1e-300), lower.tail = FALSE)
    w1 <- sqrt(z$n_patients_prostate)
    w2 <- sqrt(z$n_patients_renal)
    z[, combined_z := (w1 * z1 + w2 * z2) / sqrt(w1^2 + w2^2)]
    z[, combined_PValue := 2 * pnorm(-abs(combined_z))]
    z[, combined_FDR := p.adjust(combined_PValue, method = "BH")]
    z[, same_direction := sign(NES_prostate) == sign(NES_renal) & sign(NES_prostate) != 0]
    hallmark_meta_rows[[paste(lineage, contrast, sep = "::")]] <- z[, .(
      lineage, contrast, pathway,
      NES_prostate, pval_prostate, padj_prostate, n_patients_prostate,
      NES_renal, pval_renal, padj_renal, n_patients_renal,
      combined_z, combined_PValue, combined_FDR, same_direction
    )]
  }
}
hallmark_meta <- rbindlist(hallmark_meta_rows, use.names = TRUE)
hp <- hallmark_lopo[accession == "GSE143791"]
hr <- hallmark_lopo[accession == "GSE202813"]
setnames(hp, c("lopo_iterations", "lopo_evaluated", "lopo_direction_fraction", "lopo_complete"),
         c("lopo_iterations_prostate", "lopo_evaluated_prostate", "lopo_direction_fraction_prostate", "lopo_complete_prostate"))
setnames(hr, c("lopo_iterations", "lopo_evaluated", "lopo_direction_fraction", "lopo_complete"),
         c("lopo_iterations_renal", "lopo_evaluated_renal", "lopo_direction_fraction_renal", "lopo_complete_renal"))
hallmark_meta <- merge(hallmark_meta, hp[, .(lineage, pathway, lopo_iterations_prostate,
                                             lopo_evaluated_prostate, lopo_direction_fraction_prostate,
                                             lopo_complete_prostate)],
                      by = c("lineage", "pathway"), all.x = TRUE)
hallmark_meta <- merge(hallmark_meta, hr[, .(lineage, pathway, lopo_iterations_renal,
                                             lopo_evaluated_renal, lopo_direction_fraction_renal,
                                             lopo_complete_renal)],
                      by = c("lineage", "pathway"), all.x = TRUE)
hallmark_meta[, robust_conserved := contrast == "trend" & same_direction & combined_FDR < 0.05 &
                 pval_prostate < 0.10 & pval_renal < 0.10 &
                 abs(NES_prostate) >= 1 & abs(NES_renal) >= 1 &
                 lopo_complete_prostate == TRUE & lopo_complete_renal == TRUE &
                 lopo_direction_fraction_prostate >= 0.75 & lopo_direction_fraction_renal >= 0.75]
hallmark_meta[, abs_combined_z_sort := abs(combined_z)]
setorder(hallmark_meta, lineage, contrast, combined_FDR, -abs_combined_z_sort)
hallmark_meta[, abs_combined_z_sort := NULL]
fwrite(hallmark_meta, file.path(output_dir, "hallmark_meta_results.tsv"), sep = "\t")
robust_hallmarks <- hallmark_meta[robust_conserved == TRUE]
robust_hallmarks[, gate4b_tnk_overlap := pathway %in% names(gate4_hallmark_direction)]
robust_hallmarks[, gate4b_direction_concordant := gate4b_tnk_overlap &
                   sign(combined_z) == gate4_hallmark_direction[pathway]]
fwrite(robust_hallmarks, file.path(output_dir, "robust_conserved_hallmarks.tsv"), sep = "\t")

robust_gene_table <- copy(robust_genes)
robust_hallmark_table <- copy(robust_hallmarks)
signal_summary <- data.table(lineage = lineages)
signal_summary[, robust_genes := vapply(
  lineage,
  function(x) nrow(robust_gene_table[robust_gene_table$lineage == x]),
  integer(1)
)]
signal_summary[, robust_hallmarks := vapply(
  lineage,
  function(x) nrow(robust_hallmark_table[robust_hallmark_table$lineage == x]),
  integer(1)
)]
signal_summary[, gate4b_overlap_genes := vapply(
  lineage,
  function(x) nrow(robust_gene_table[lineage == x & gate4b_direction_concordant == TRUE]),
  integer(1)
)]
signal_summary[, gate4b_overlap_hallmarks := vapply(
  lineage,
  function(x) nrow(robust_hallmark_table[lineage == x & gate4b_direction_concordant == TRUE]),
  integer(1)
)]
signal_summary[, strong_threshold := gate4b_overlap_genes >= 10 |
                 gate4b_overlap_hallmarks >= 3]
signal_summary[, support_threshold := gate4b_overlap_genes >= 3 |
                 gate4b_overlap_hallmarks >= 1]
fwrite(signal_summary, file.path(output_dir, "signal_summary.tsv"), sep = "\t")

model_integrity <- nrow(model_qc_table) == 6 && all(model_qc_table$converged) &&
  all(model_qc_table$genes_tested >= 5000) &&
  all(model_qc_table$trend_design_rank == model_qc_table$trend_design_columns) &&
  all(model_qc_table$categorical_design_rank == model_qc_table$categorical_design_columns)
strong_n <- sum(signal_summary$strong_threshold)
support_n <- sum(signal_summary$support_threshold)
robust_composition_n <- nrow(composition_meta[robust_composition == TRUE])
composition_detected <- robust_composition_n >= 1
within_state_detected <- strong_n >= 1
mechanism_class <- if (composition_detected && within_state_detected) {
  "MIXED_COMPOSITION_AND_WITHIN_STATE"
} else if (within_state_detected) {
  "WITHIN_STATE_DOMINANT"
} else if (composition_detected) {
  "COMPOSITION_DOMINANT"
} else {
  "UNRESOLVED"
}
decision <- if (!model_integrity) {
  "NO-GO"
} else if (strong_n >= 2) {
  "GO"
} else if (strong_n >= 1 || composition_detected) {
  "CONDITIONAL GO"
} else {
  "NO-GO"
}

state_lines <- unlist(lapply(lineages, function(x) c(
  paste0("- ", x, " robust within-state genes: ",
         signal_summary[lineage == x, robust_genes]),
  paste0("- ", x, " Gate4b-concordant genes: ",
         signal_summary[lineage == x, gate4b_overlap_genes]),
  paste0("- ", x, " robust within-state Hallmarks: ",
         signal_summary[lineage == x, robust_hallmarks]),
  paste0("- ", x, " Gate4b-concordant Hallmarks: ",
         signal_summary[lineage == x, gate4b_overlap_hallmarks])
)))
decision_lines <- c(
  "# Gate5b T/NK decomposition decision",
  "",
  paste0("**", decision, "**"),
  "",
  paste0("- Six cohort-state expression models passed integrity checks: ",
         ifelse(model_integrity, "PASS", "FAIL")),
  paste0("- Robust composition shifts: ", robust_composition_n),
  paste0("- Mechanism class: ", mechanism_class),
  state_lines,
  "",
  "Composition and expression models use the patient as the unit of inference.",
  "This gate distinguishes statistical sources of the broad T/NK gradient and does not establish causality."
)
writeLines(decision_lines, file.path(output_dir, "gate5b_decision.md"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("GATE5B_COMPLETE\tdecision=", decision,
        "\tmechanism=", mechanism_class,
        "\tcomposition_shifts=", robust_composition_n,
        "\tstrong_states=", strong_n,
        "\tsupport_states=", support_n)
