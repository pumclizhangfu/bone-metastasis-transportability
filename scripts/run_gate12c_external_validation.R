#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

options(stringsAsFactors = FALSE, warn = 1)
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) != 6L) {
  stop("Usage: run_gate12c_external_validation.R PROJECT GSE_RAW OEP_RAW LUNG_RAW BREAST_TARGET_DIR OUT_DIR")
}
project <- normalizePath(argv[[1]], mustWork = TRUE)
gse_raw <- normalizePath(argv[[2]], mustWork = TRUE)
oep_raw <- normalizePath(argv[[3]], mustWork = TRUE)
lung_raw <- normalizePath(argv[[4]], mustWork = TRUE)
breast_dir <- normalizePath(argv[[5]], mustWork = TRUE)
out <- normalizePath(argv[[6]], mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(20260808)

log_con <- file(file.path(out, "gate12c_external_validation.log"), open = "wt")
on.exit(try(close(log_con), silent = TRUE), add = TRUE)
log_msg <- function(...) {
  z <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(z, "\n"); writeLines(z, log_con); flush(log_con)
}

genes <- c("CCL4", "CXCL16", "SPP1", "CCR5", "CXCR6", "CD44")
axes <- data.table(
  ligand = c("CCL4", "CXCL16", "SPP1"),
  receptor = c("CCR5", "CXCR6", "CD44"),
  sender_state = "MACROPHAGE",
  receiver_state = c("CD8_TEX", "CD8_TEX", "CD4_TREG")
)
axes[, axis := paste(ligand, receptor, sep = " -> ")]
min_cells <- 20L
min_detect <- 0.05

aggregate_symbols <- function(mat, symbols) {
  keep <- !is.na(symbols) & nzchar(symbols)
  mat <- mat[keep, , drop = FALSE]; symbols <- symbols[keep]
  keep_gene <- symbols %chin% genes
  mat <- mat[keep_gene, , drop = FALSE]; symbols <- symbols[keep_gene]
  if (!length(symbols)) {
    ans <- Matrix(0, nrow = length(genes), ncol = ncol(mat), sparse = TRUE)
    rownames(ans) <- genes
    return(as(ans, "dgCMatrix"))
  }
  map <- sparseMatrix(i = match(symbols, genes), j = seq_along(symbols), x = 1,
                      dims = c(length(genes), length(symbols)))
  ans <- map %*% mat
  rownames(ans) <- genes
  as(ans, "dgCMatrix")
}

read_target_10x <- function(matrix_path, feature_path, barcode_path) {
  feat <- fread(feature_path, header = FALSE, showProgress = FALSE)
  symbols <- if (ncol(feat) >= 2L) as.character(feat[[2L]]) else as.character(feat[[1L]])
  gex <- if (ncol(feat) >= 3L) as.character(feat[[3L]]) == "Gene Expression" else rep(TRUE, nrow(feat))
  barcodes <- as.character(fread(barcode_path, header = FALSE, showProgress = FALSE)[[1L]])
  raw <- as(Matrix::readMM(gzfile(matrix_path)), "dgCMatrix")
  stopifnot(nrow(raw) == nrow(feat), ncol(raw) == length(barcodes))
  ans <- aggregate_symbols(raw[gex, , drop = FALSE], symbols[gex])
  colnames(ans) <- barcodes
  ans
}

summarise_library <- function(mat, ann, dataset, condition_value = NULL) {
  states <- c("MACROPHAGE", "CD8_TEX", "CD4_TREG")
  ann <- ann[final_assignment %chin% states]
  if (!nrow(ann)) return(NULL)
  idx <- match(ann$cell_barcode, colnames(mat))
  if (anyNA(idx)) stop("Assignment/barcode mismatch in ", dataset)
  mat <- mat[, idx, drop = FALSE]
  rows <- vector("list", length(states) * length(genes))
  k <- 0L
  for (state in states) {
    take <- which(ann$final_assignment == state)
    if (!length(take)) next
    denom <- pmax(as.numeric(ann$nCount_RNA[take]), 1)
    for (gene in genes) {
      k <- k + 1L
      x <- as.numeric(mat[gene, take, drop = TRUE])
      rows[[k]] <- data.table(
        dataset = dataset,
        patient_id = ann$patient_id[take],
        cancer = ann$cancer_code[take],
        condition = if (is.null(condition_value)) ann$condition[take] else condition_value,
        state = state, gene = gene,
        n_cells = 1L,
        detected_cells = as.integer(x > 0),
        sum_counts = x,
        sum_total = denom,
        sum_expression = log1p(1e4 * x / denom)
      )
    }
  }
  rbindlist(rows[seq_len(k)], use.names = TRUE)
}

collapse_patient <- function(x) {
  ans <- x[, .(
    n_cells = sum(n_cells),
    detected_cells = sum(detected_cells),
    sum_counts = sum(sum_counts),
    sum_total = sum(sum_total),
    mean_expression = sum(sum_expression) / sum(n_cells)
  ), by = .(dataset, patient_id, cancer, condition, state, gene)]
  ans[, `:=`(
    detected_fraction = detected_cells / n_cells,
    pooled_log1p_cpm = log1p(1e6 * sum_counts / pmax(sum_total, 1))
  )]
  ans
}

gse_cache <- file.path(out, "gse266330_state_gene_expression.tsv.gz")
if (file.exists(gse_cache)) {
  log_msg("RESUME GSE266330 cache")
  gse_state <- fread(gse_cache)
} else {
gse_ann <- fread(file.path(project, "results/gate9a_reference/raw_reconstruction_v1/cell_state_assignments.tsv.gz"))
gse_ann[, condition := fifelse(cancer == "ctrl", "healthy_bm", "bone_metastasis")]
gse_ann[, cancer_code := cancer]
gse_specs <- unique(gse_ann[, .(library_id)])
gse_rows <- vector("list", nrow(gse_specs))
for (i in seq_len(nrow(gse_specs))) {
  id <- gse_specs$library_id[[i]]
  prefix <- file.path(gse_raw, id)
  paths <- c(matrix = paste0(prefix, "_matrix.mtx.gz"),
             feature = paste0(prefix, "_features.tsv.gz"),
             barcode = paste0(prefix, "_barcodes.tsv.gz"))
  if (!all(file.exists(paths))) stop("Missing GSE266330 files for ", id)
  log_msg("GSE266330 ", i, "/", nrow(gse_specs), " ", id)
  mat <- read_target_10x(paths[["matrix"]], paths[["feature"]], paths[["barcode"]])
  gse_rows[[i]] <- summarise_library(mat, gse_ann[library_id == id], "GSE266330")
  rm(mat); gc(FALSE)
}
gse_state <- collapse_patient(rbindlist(gse_rows, use.names = TRUE))
fwrite(gse_state, file.path(out, "gse266330_state_gene_expression.tsv.gz"), sep = "\t")
}

oep_cache <- file.path(out, "oep005136_state_gene_expression.tsv.gz")
if (file.exists(oep_cache)) {
  log_msg("RESUME OEP005136 cache")
  oep_state <- fread(oep_cache)
} else {
oep_ann <- fread(file.path(project, "results/gate9b_validation/external_validation_v1/oep005136_cell_state_assignments.tsv.gz"))
oep_ann <- oep_ann[primary_analysis == TRUE]
oep_ann[, condition := "bone_metastasis"]
oep_specs <- unique(oep_ann[, .(archive_dir)])
oep_rows <- vector("list", nrow(oep_specs))
for (i in seq_len(nrow(oep_specs))) {
  id <- oep_specs$archive_dir[[i]]
  path <- file.path(oep_raw, id)
  paths <- c(matrix = file.path(path, "matrix.mtx.gz"),
             feature = file.path(path, "features.tsv.gz"),
             barcode = file.path(path, "barcodes.tsv.gz"))
  if (!all(file.exists(paths))) stop("Missing OEP files for ", id)
  log_msg("OEP005136 ", i, "/", nrow(oep_specs), " ", id)
  mat <- read_target_10x(paths[["matrix"]], paths[["feature"]], paths[["barcode"]])
  oep_rows[[i]] <- summarise_library(mat, oep_ann[archive_dir == id], "OEP005136", "bone_metastasis")
  rm(mat); gc(FALSE)
}
oep_state <- collapse_patient(rbindlist(oep_rows, use.names = TRUE))
fwrite(oep_state, file.path(out, "oep005136_state_gene_expression.tsv.gz"), sep = "\t")
}

support_ann <- fread(file.path(project, "results/gate10s_supportive_projection/full_v1/gate10s_cell_assignments.tsv.gz"))
support_ann[, condition := "bone_metastasis"]
lung_cache <- file.path(out, "gse225209_state_gene_expression.tsv.gz")
if (file.exists(lung_cache)) {
  log_msg("RESUME GSE225209 cache")
  lung_state <- fread(lung_cache)
} else {
lung_spec <- data.table(
  accession = c("GSM7041480", "GSM7041481", "GSM7041482"),
  stem = c("sg1", "sg2", "sg3"), sample_or_lesion = c("sz", "s13", "s14"),
  patient_id = c("LUCA_BM_01", "LUCA_BM_02", "LUCA_BM_03")
)
lung_rows <- vector("list", nrow(lung_spec))
for (i in seq_len(nrow(lung_spec))) {
  prefix <- file.path(lung_raw, paste(lung_spec$accession[[i]], lung_spec$stem[[i]], sep = "_"))
  paths <- c(matrix = paste0(prefix, "-matrix.mtx.gz"),
             feature = paste0(prefix, "-features.tsv.gz"),
             barcode = paste0(prefix, "-barcodes.tsv.gz"))
  log_msg("GSE225209 ", i, "/", nrow(lung_spec), " ", lung_spec$patient_id[[i]])
  mat <- read_target_10x(paths[["matrix"]], paths[["feature"]], paths[["barcode"]])
  ann <- support_ann[dataset == "GSE225209" & patient_id == lung_spec$patient_id[[i]]]
  ann[, cell_barcode := sub("^[^:]+::", "", cell_id)]
  ann[, cancer_code := "LUCA"]
  lung_rows[[i]] <- summarise_library(mat, ann, "GSE225209", "bone_metastasis")
  rm(mat); gc(FALSE)
}
lung_state <- collapse_patient(rbindlist(lung_rows, use.names = TRUE))
fwrite(lung_state, file.path(out, "gse225209_state_gene_expression.tsv.gz"), sep = "\t")
}

read_breast <- function(path) {
  z <- fread(path)
  mat <- as.matrix(z[, -1L, with = FALSE]); storage.mode(mat) <- "double"
  rownames(mat) <- z[[1L]]
  missing <- setdiff(genes, rownames(mat))
  if (length(missing)) {
    add <- matrix(0, nrow = length(missing), ncol = ncol(mat),
                  dimnames = list(missing, colnames(mat)))
    mat <- rbind(mat, add)
  }
  mat[genes, , drop = FALSE]
}
breast_files <- c(
  BoM12 = "GSE190772_BoM_logCounts.gate12c_genes.tsv.gz",
  BoM7 = "GSM6870693_BoM7_scRNA_LogCounts.gate12c_genes.tsv.gz",
  BoM8 = "GSM6870694_BoM8_scRNA_LogCounts.gate12c_genes.tsv.gz"
)
breast_mats <- setNames(lapply(file.path(breast_dir, breast_files), read_breast), names(breast_files))
breast_rows <- list()
for (lesion in c("BoM1", "BoM2", "BoM7", "BoM8")) {
  key <- if (lesion %chin% c("BoM1", "BoM2")) "BoM12" else lesion
  mat <- breast_mats[[key]]
  ann <- support_ann[dataset == "GSE190772" & sample_or_lesion == lesion &
                     final_assignment %chin% c("MACROPHAGE", "CD8_TEX", "CD4_TREG")]
  if (!nrow(ann)) next
  idx <- match(ann$cell_id, colnames(mat))
  if (anyNA(idx)) stop("Breast assignment/matrix mismatch: ", lesion)
  mat <- mat[, idx, drop = FALSE]
  for (state in c("MACROPHAGE", "CD8_TEX", "CD4_TREG")) {
    take <- which(ann$final_assignment == state)
    if (!length(take)) next
    for (gene in genes) {
      x <- as.numeric(mat[gene, take])
      breast_rows[[length(breast_rows) + 1L]] <- data.table(
        dataset = "GSE190772", lesion = lesion, patient_id = ann$patient_id[take][[1L]],
        cancer = "BRCA", condition = "bone_metastasis", state = state, gene = gene,
        n_cells = length(take), detected_cells = sum(x > 0),
        mean_expression = mean(x)
      )
    }
  }
}
breast_lesion <- rbindlist(breast_rows)
breast_state <- breast_lesion[, .(
  n_cells = sum(n_cells), detected_cells = sum(detected_cells),
  mean_expression = weighted.mean(mean_expression, n_cells)
), by = .(dataset, patient_id, cancer, condition, state, gene)]
breast_state[, `:=`(detected_fraction = detected_cells / n_cells,
                    sum_counts = NA_real_, sum_total = NA_real_, pooled_log1p_cpm = NA_real_)]
fwrite(breast_state, file.path(out, "gse190772_state_gene_expression.tsv.gz"), sep = "\t")

state_expression <- rbindlist(list(gse_state, oep_state, lung_state, breast_state),
                              use.names = TRUE, fill = TRUE)
setcolorder(state_expression, c("dataset", "patient_id", "cancer", "condition", "state", "gene",
                                "n_cells", "detected_cells", "detected_fraction", "mean_expression",
                                "sum_counts", "sum_total", "pooled_log1p_cpm"))
fwrite(state_expression, file.path(out, "external_state_gene_expression.tsv.gz"), sep = "\t")

axis_patients <- rbindlist(lapply(seq_len(nrow(axes)), function(i) {
  a <- axes[i]
  ligand <- state_expression[state == a$sender_state & gene == a$ligand,
    .(dataset, patient_id, cancer, condition,
      sender_cells = n_cells, ligand_detected_fraction = detected_fraction,
      ligand_expression = mean_expression)]
  receptor <- state_expression[state == a$receiver_state & gene == a$receptor,
    .(dataset, patient_id, cancer, condition,
      receiver_cells = n_cells, receptor_detected_fraction = detected_fraction,
      receptor_expression = mean_expression)]
  z <- merge(ligand, receptor, by = c("dataset", "patient_id", "cancer", "condition"), all = TRUE)
  z[, `:=`(axis = a$axis, ligand = a$ligand, receptor = a$receptor,
           sender_state = a$sender_state, receiver_state = a$receiver_state)]
  z[, eligible := !is.na(sender_cells) & !is.na(receiver_cells) &
                    sender_cells >= min_cells & receiver_cells >= min_cells]
  z[, patient_support := eligible & ligand_detected_fraction >= min_detect &
                           receptor_detected_fraction >= min_detect]
  z[, axis_score := ligand_expression + receptor_expression]
  z
}), use.names = TRUE, fill = TRUE)
setorder(axis_patients, axis, dataset, patient_id)
fwrite(axis_patients, file.path(out, "external_patient_axis_support.tsv"), sep = "\t")

dataset_grid <- CJ(axis = axes$axis, dataset = c("GSE266330", "OEP005136", "GSE225209", "GSE190772"))
dataset_summary <- axis_patients[condition == "bone_metastasis", .(
  eligible_patients = sum(eligible),
  supporting_patients = sum(patient_support),
  support_fraction = if (sum(eligible)) mean(patient_support[eligible]) else NA_real_,
  supporting_origins = uniqueN(cancer[patient_support])
), by = .(axis, dataset)]
dataset_summary <- merge(dataset_grid, dataset_summary, by = c("axis", "dataset"), all.x = TRUE)
dataset_summary[is.na(eligible_patients), `:=`(eligible_patients = 0L, supporting_patients = 0L,
                                               supporting_origins = 0L)]
dataset_summary[, dataset_pass := fifelse(
  dataset %chin% c("GSE266330", "OEP005136"),
  eligible_patients >= 10L & support_fraction >= 0.50 & supporting_origins >= 3L,
  eligible_patients >= 1L & supporting_patients >= 1L
)]
fwrite(dataset_summary, file.path(out, "external_dataset_axis_summary.tsv"), sep = "\t")

gse_contradiction <- rbindlist(lapply(axes$axis, function(ax) {
  z <- axis_patients[axis == ax & dataset == "GSE266330" & eligible == TRUE]
  bm <- z[condition == "bone_metastasis", axis_score]
  ctrl <- z[condition == "healthy_bm", axis_score]
  diff <- if (length(bm) && length(ctrl)) median(bm) - median(ctrl) else NA_real_
  evaluable <- length(bm) >= 3L && length(ctrl) >= 3L
  p_lower <- if (evaluable)
    suppressWarnings(wilcox.test(bm, ctrl, alternative = "less", exact = FALSE)$p.value) else NA_real_
  data.table(axis = ax, bm_n = length(bm), control_n = length(ctrl),
             median_difference = diff, p_bm_lower = p_lower,
             contradiction_evaluable = evaluable,
             contradiction = if (evaluable) diff < 0 && p_lower < 0.10 else NA)
}))
fwrite(gse_contradiction, file.path(out, "gse266330_control_contradiction.tsv"), sep = "\t")

wide <- dcast(dataset_summary, axis ~ dataset, value.var = c("dataset_pass", "support_fraction", "eligible_patients"))
decision <- merge(axes, wide, by = "axis", all.x = TRUE)
decision <- merge(decision, gse_contradiction[, .(axis,
                                                  gse_contradiction_evaluable = contradiction_evaluable,
                                                  gse_contradiction = contradiction,
                                                  gse_bm_control_difference = median_difference,
                                                  gse_bm_lower_p = p_bm_lower)], by = "axis", all.x = TRUE)
shortlist <- fread(file.path(project, "results/gate12c_sender_receiver/gate12c_provisional_shortlist.tsv"))
decision <- merge(decision, shortlist[, .(ligand, receptor, nichenet_rank, pearson)],
                  by = c("ligand", "receptor"), all.x = TRUE)
pass_cols <- grep("^dataset_pass_", names(decision), value = TRUE)
decision[, supporting_dataset_count := rowSums(.SD, na.rm = TRUE), .SDcols = pass_cols]
decision[, min_core_fraction := pmin(support_fraction_GSE266330, support_fraction_OEP005136, na.rm = FALSE)]
decision[, external_human_support := dataset_pass_GSE266330 & dataset_pass_OEP005136 &
             (is.na(gse_contradiction) | !gse_contradiction) &
             (dataset_pass_GSE225209 | dataset_pass_GSE190772)]
decision[, final_freeze_eligible := FALSE]
eligible <- decision[external_human_support == TRUE][order(-supporting_dataset_count,
                                                            -min_core_fraction,
                                                            -support_fraction_OEP005136,
                                                            nichenet_rank)]
if (nrow(eligible)) {
  top <- eligible[1L]
  tied <- eligible[supporting_dataset_count == top$supporting_dataset_count &
                   min_core_fraction == top$min_core_fraction &
                   support_fraction_OEP005136 == top$support_fraction_OEP005136 &
                   nichenet_rank == top$nichenet_rank]
  if (nrow(tied) == 1L) decision[axis == top$axis, final_freeze_eligible := TRUE]
}
setorder(decision, -final_freeze_eligible, -external_human_support,
         -supporting_dataset_count, -min_core_fraction, -support_fraction_OEP005136, nichenet_rank)
fwrite(decision, file.path(out, "gate12c_external_axis_decision.tsv"), sep = "\t")

freeze <- decision[final_freeze_eligible == TRUE,
  .(axis, ligand, receptor, sender_state, receiver_state,
    knockout_target = receptor, knockout_context = receiver_state,
    supporting_dataset_count, min_core_fraction, nichenet_rank,
    freeze_status = "PRIMARY_RULE_SELECTED_PENDING_ROBUSTNESS")]
fwrite(freeze, file.path(out, "gate12c_target_freeze.tsv"), sep = "\t")

theme_pub <- theme_classic(base_size = 9, base_family = "Arial") +
  theme(plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 8, colour = "#444444"),
        plot.tag = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 7),
        strip.background = element_rect(fill = "#F2F2F2", colour = NA),
        strip.text = element_text(face = "bold", size = 8))

plot_a <- dataset_summary[eligible_patients > 0]
plot_a[, dataset := factor(dataset, levels = c("GSE266330", "OEP005136", "GSE225209", "GSE190772"))]
p1 <- ggplot(plot_a, aes(support_fraction, axis, colour = dataset, size = eligible_patients)) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "#777777") +
  geom_point(alpha = 0.9) +
  scale_x_continuous(limits = c(0, 1), labels = scales::percent) +
  scale_colour_manual(values = c("#0072B2", "#D55E00", "#009E73", "#CC79A7")) +
  labs(title = "Independent patient-level co-detection", subtitle = "Dashed line: frozen 50% core-cohort threshold",
       x = "Supporting eligible patients", y = NULL, colour = "Dataset", size = "Eligible patients") + theme_pub

gse_origin <- axis_patients[dataset == "GSE266330" & condition == "bone_metastasis" & eligible == TRUE,
  .(eligible = .N, support_fraction = mean(patient_support)), by = .(axis, cancer)]
p2 <- ggplot(gse_origin, aes(cancer, axis, fill = support_fraction)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = eligible), size = 2.5) +
  scale_fill_gradient(low = "white", high = "#2166AC", limits = c(0, 1), labels = scales::percent) +
  labs(title = "GSE266330 cross-origin support",
       subtitle = "Numbers are eligible patients; healthy comparison was not evaluable",
       x = "Cancer origin", y = NULL, fill = "Support fraction") + theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

oep_origin <- axis_patients[dataset == "OEP005136" & eligible == TRUE,
  .(eligible = .N, support_fraction = mean(patient_support)), by = .(axis, cancer)]
p3 <- ggplot(oep_origin, aes(cancer, axis, fill = support_fraction)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = eligible), size = 2.5) +
  scale_fill_gradient(low = "white", high = "#B2182B", limits = c(0, 1), labels = scales::percent) +
  labs(title = "OEP005136 cross-origin support", subtitle = "Numbers are eligible patients",
       x = "Cancer origin", y = NULL, fill = "Support fraction") + theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_d <- copy(decision)
plot_d[, status := fifelse(final_freeze_eligible, "Primary-rule selected",
                          fifelse(external_human_support, "External pass", "Fail"))]
p4 <- ggplot(plot_d, aes(min_core_fraction, axis, colour = status, size = supporting_dataset_count)) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "#777777") +
  geom_point(alpha = 0.92) +
  scale_x_continuous(limits = c(0, 1), labels = scales::percent) +
  scale_colour_manual(values = c("Primary-rule selected" = "#B2182B", "External pass" = "#E69F00", "Fail" = "#777777")) +
  labs(title = "Prespecified primary-rule decision", x = "Minimum core-cohort support", y = NULL,
       colour = "Decision", size = "Supporting datasets") + theme_pub

fig <- (p1 | p2) / (p3 | p4) + plot_annotation(tag_levels = "A")
ggsave(file.path(out, "Figure5_external_axis_validation.pdf"), fig,
       width = 12.0, height = 8.8, device = cairo_pdf, bg = "white")
ggsave(file.path(out, "Figure5_external_axis_validation.png"), fig,
       width = 12.0, height = 8.8, dpi = 360, bg = "white")

status <- if (nrow(freeze) == 1L) "PRIMARY_RULE_SELECTION_COMPLETE_ROBUSTNESS_PENDING" else "NO_PRIMARY_TARGET_GATE12D_BLOCKED"
audit <- c(
  "# Gate12C-external validation audit", "",
  "- Plan version: gate12c_external_frozen_20260808",
  "- Statistical unit: patient/donor",
  "- Minimum sender and receiver cells per patient: 20",
  "- Cell-level detection threshold: 5% in each frozen state",
  "- Core cohort threshold: >=10 eligible patients, >=50% support, >=3 origins",
  "- GSE266330 technical replicates collapsed to patient/donor",
  "- OEP005136 restricted to 49 newly generated primary solid-tumour mBone patients",
  "- GSE225209 and GSE190772 remain separate supportive cohorts",
  "- GSE266330 healthy-control comparison: NOT EVALUABLE (0 controls met the paired state-size rule)",
  paste0("- Externally supported axes: ", paste(decision[external_human_support == TRUE, axis], collapse = "; ")),
  paste0("- Primary-rule selected target: ", if (nrow(freeze)) paste0(freeze$knockout_target, " in ", freeze$knockout_context) else "none"),
  paste0("- Status: ", status), "",
  "## Interpretation boundary", "",
  "Co-detection supports communication plausibility but does not demonstrate physical binding or causal signalling. Any subsequent knockout is an in-silico network perturbation only."
)
writeLines(audit, file.path(out, "GATE12C_EXTERNAL_AUDIT.md"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
log_msg("Complete status=", status)
cat("GATE12C_EXTERNAL_STATUS=", status, "\n", sep = "")
