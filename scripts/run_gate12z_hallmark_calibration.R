#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(edgeR)
  library(matrixStats)
  library(msigdbr)
  library(igraph)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(paste(
    "Usage: run_gate12z_hallmark_calibration.R <gate12g_model.rds>",
    "<state_pseudobulk_counts.rds> <gate12v2_activity.tsv.gz>",
    "<gate12v2_primary.tsv> <outdir> <workers>"
  ))
}

model_file <- normalizePath(args[[1L]], mustWork = TRUE)
pseudobulk_file <- normalizePath(args[[2L]], mustWork = TRUE)
activity_file <- normalizePath(args[[3L]], mustWork = TRUE)
primary_file <- normalizePath(args[[4L]], mustWork = TRUE)
outdir <- args[[5L]]
workers <- as.integer(args[[6L]])
if (!is.finite(workers) || workers < 1L) stop("workers must be >=1")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

seed <- 1208L
n_boot <- 1999L
n_null <- 999L
set.seed(seed)

zscore <- function(x) {
  s <- sd(x)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  (x - mean(x)) / s
}

model <- readRDS(model_file)
pb <- readRDS(pseudobulk_file)
frozen_activity <- fread(activity_file)
mixed_primary <- fread(primary_file)
if (ncol(model$x_raw_clr) != 26L || model$selected_k != 4L) stop("Unexpected frozen Axis1 model")
if (nrow(frozen_activity) != 228L || uniqueN(frozen_activity$sample_id) != 41L ||
    uniqueN(frozen_activity$patient_id) != 18L) stop("Unexpected frozen Hallmark activity scope")

pathways <- grep("^HALLMARK_", names(frozen_activity), value = TRUE)
if (length(pathways) != 50L) stop("Expected exactly 50 frozen Hallmark columns")
pathways <- sort(pathways)

# Rebuild the exact Gate12V2 standardized-gene universe and verify that the
# frozen Hallmark activity values are reproduced before any new inference.
pb_meta <- as.data.table(pb$metadata)
keep_pb <- pb_meta$n_cells >= 20L & pb_meta$sample_id %chin% model$meta$sample_id
pb_meta <- copy(pb_meta[keep_pb])
counts <- pb$counts[, keep_pb, drop = FALSE]
if (ncol(counts) != 228L) stop("Expected 228 retained pseudobulks")
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge, method = "TMM")
log_cpm <- cpm(dge, log = TRUE, prior.count = 2)
gene_mean <- rowMeans(log_cpm)
gene_sd <- matrixStats::rowSds(log_cpm)
valid_gene <- is.finite(gene_mean) & is.finite(gene_sd) & gene_sd > 0
log_cpm <- log_cpm[valid_gene, , drop = FALSE]
gene_z <- sweep(sweep(log_cpm, 1L, gene_mean[valid_gene], "-"),
                1L, gene_sd[valid_gene], "/")
upper_gene <- toupper(rownames(gene_z))

msig <- as.data.table(msigdbr(db_species = "HS", species = "human", collection = "H"))
hallmark_sets <- lapply(split(msig$gene_symbol, msig$gs_name), function(x) unique(toupper(x)))
hallmark_sets <- hallmark_sets[pathways]
if (length(hallmark_sets) != 50L || any(vapply(hallmark_sets, is.null, logical(1)))) {
  stop("Current Hallmark membership does not match frozen pathway names")
}
measured_idx <- lapply(hallmark_sets, function(g) unique(na.omit(match(g, upper_gene))))
if (any(lengths(measured_idx) < 15L)) stop("Every Hallmark must retain at least 15 measured genes")
reconstructed_activity <- vapply(measured_idx, function(idx) {
  colMeans(gene_z[idx, , drop = FALSE])
}, numeric(ncol(gene_z)))
colnames(reconstructed_activity) <- pathways
rownames(reconstructed_activity) <- colnames(gene_z)
frozen_order <- match(rownames(reconstructed_activity), frozen_activity$pseudobulk_id)
if (anyNA(frozen_order)) stop("Frozen pseudobulk IDs do not match reconstructed counts")
max_activity_error <- max(abs(
  reconstructed_activity - as.matrix(frozen_activity[frozen_order, ..pathways])
))
if (!is.finite(max_activity_error) || max_activity_error > 1e-10) {
  stop("Frozen Hallmark activity reconstruction error exceeds 1e-10: ", max_activity_error)
}

dat <- copy(frozen_activity[frozen_order])
dat[, `:=`(
  compartment = factor(compartment, levels = c("distal", "involved", "tumor")),
  cancer = factor(cancer),
  meta_state = factor(meta_state),
  patient_id_chr = as.character(patient_id),
  Axis1_z = zscore(Axis1_z)
)]
Y_raw <- as.matrix(dat[, ..pathways])
Y_z <- apply(Y_raw, 2L, zscore)
colnames(Y_z) <- pathways
X <- model.matrix(~ Axis1_z + compartment + cancer + meta_state, data = dat)
axis_col <- match("Axis1_z", colnames(X))
if (is.na(axis_col)) stop("Axis1 term absent from pseudobulk design matrix")

fit_multi <- function(Xb, Yb, coefficient_index) {
  fit <- lm.fit(Xb, Yb)
  cf <- fit$coefficients
  if (is.null(dim(cf))) cf <- matrix(cf, ncol = 1L)
  out <- cf[coefficient_index, ]
  as.numeric(out)
}

fixed_beta <- fit_multi(X, Y_z, axis_col)
names(fixed_beta) <- pathways

# Patient-cluster bootstrap: all pseudobulks from a patient move together.
patients <- unique(dat$patient_id_chr)
patient_rows <- split(seq_len(nrow(dat)), dat$patient_id_chr)
set.seed(seed + 10L)
pseudobulk_boot <- vector("list", n_boot)
for (b in seq_len(n_boot)) {
  draw <- sample(patients, length(patients), replace = TRUE)
  idx <- unlist(patient_rows[draw], use.names = FALSE)
  coef_b <- try(fit_multi(X[idx, , drop = FALSE], Y_z[idx, , drop = FALSE], axis_col), silent = TRUE)
  if (inherits(coef_b, "try-error") || length(coef_b) != length(pathways)) {
    coef_b <- rep(NA_real_, length(pathways))
  }
  pseudobulk_boot[[b]] <- data.table(iteration = b, pathway = pathways, beta = coef_b)
}
pseudobulk_boot <- rbindlist(pseudobulk_boot)
fwrite(pseudobulk_boot, file.path(outdir, "hallmark_patient_cluster_bootstrap.tsv.gz"), sep = "\t")

robust_summary <- pseudobulk_boot[, .(
  evaluable_bootstraps = sum(is.finite(beta)),
  ci_low = quantile(beta, 0.025, na.rm = TRUE),
  ci_high = quantile(beta, 0.975, na.rm = TRUE),
  median_boot_beta = median(beta, na.rm = TRUE)
), by = pathway]
robust_summary[, fixed_beta := fixed_beta[pathway]]
mixed_map <- mixed_primary[match(robust_summary$pathway, mixed_primary$pathway)]
robust_summary[, `:=`(
  mixed_beta = mixed_map$beta,
  mixed_q_value = mixed_map$q_value,
  original_mixed_singular = mixed_map$singular
)]
direction_fraction <- pseudobulk_boot[robust_summary[, .(pathway, primary_sign = sign(mixed_beta))],
                                     on = "pathway"][, .(
  bootstrap_same_mixed_direction_fraction = mean(sign(beta[is.finite(beta)]) ==
                                                    primary_sign[is.finite(beta)])
), by = pathway]
robust_summary <- merge(robust_summary, direction_fraction, by = "pathway", sort = FALSE)
robust_summary[, robust_direction_supported := evaluable_bootstraps >= 1900L &
                 ((ci_low > 0 & ci_high > 0) | (ci_low < 0 & ci_high < 0)) &
                 bootstrap_same_mixed_direction_fraction >= 0.90]

# Sample-level cell-count-weighted aggregation and patient-cluster bootstrap.
sample_weighted <- dat[, {
  w <- n_cells / sum(n_cells)
  vals <- colSums(sweep(as.matrix(.SD), 1L, w, "*"))
  as.list(vals)
}, by = .(sample_id, patient_id_chr, compartment, cancer, frozen_Axis1), .SDcols = pathways]
if (nrow(sample_weighted) != 41L || uniqueN(sample_weighted$patient_id_chr) != 18L) {
  stop("Sample aggregation must retain 41 samples and 18 patients")
}
sample_weighted[, Axis1_z := zscore(frozen_Axis1)]
SY_raw <- as.matrix(sample_weighted[, ..pathways])
SY_z <- apply(SY_raw, 2L, zscore)
colnames(SY_z) <- pathways
SX <- model.matrix(~ Axis1_z + compartment + cancer, data = sample_weighted)
s_axis_col <- match("Axis1_z", colnames(SX))
sample_beta <- fit_multi(SX, SY_z, s_axis_col)
names(sample_beta) <- pathways
sample_rows <- split(seq_len(nrow(sample_weighted)), sample_weighted$patient_id_chr)
set.seed(seed + 20L)
sample_boot <- vector("list", n_boot)
for (b in seq_len(n_boot)) {
  draw <- sample(names(sample_rows), length(sample_rows), replace = TRUE)
  idx <- unlist(sample_rows[draw], use.names = FALSE)
  coef_b <- try(fit_multi(SX[idx, , drop = FALSE], SY_z[idx, , drop = FALSE], s_axis_col), silent = TRUE)
  if (inherits(coef_b, "try-error") || length(coef_b) != length(pathways)) {
    coef_b <- rep(NA_real_, length(pathways))
  }
  sample_boot[[b]] <- data.table(iteration = b, pathway = pathways, beta = coef_b)
}
sample_boot <- rbindlist(sample_boot)
fwrite(sample_boot, file.path(outdir, "hallmark_sample_aggregated_bootstrap.tsv.gz"), sep = "\t")
sample_summary <- sample_boot[, .(
  evaluable_bootstraps = sum(is.finite(beta)),
  ci_low = quantile(beta, 0.025, na.rm = TRUE),
  ci_high = quantile(beta, 0.975, na.rm = TRUE),
  median_boot_beta = median(beta, na.rm = TRUE)
), by = pathway]
sample_summary[, fixed_beta := sample_beta[pathway]]
sample_summary[, sample_direction_supported := evaluable_bootstraps >= 1900L &
                 ((ci_low > 0 & ci_high > 0) | (ci_low < 0 & ci_high < 0))]

# Gene-level adjusted effects, frozen top-20 same-sign driver sets.
gene_fit <- lm.fit(X, t(gene_z))
gene_beta <- as.numeric(gene_fit$coefficients[axis_col, ])
names(gene_beta) <- upper_gene
driver_rows <- list()
driver_sets <- vector("list", length(pathways)); names(driver_sets) <- pathways
for (p in pathways) {
  idx <- measured_idx[[p]]
  genes <- upper_gene[idx]
  b <- gene_beta[genes]
  keep <- is.finite(b) & sign(b) == sign(fixed_beta[[p]]) & sign(b) != 0
  ord <- order(abs(b[keep]), decreasing = TRUE)
  chosen <- genes[keep][ord]
  if (length(chosen) > 20L) chosen <- chosen[seq_len(20L)]
  driver_sets[[p]] <- chosen
  if (length(chosen)) {
    driver_rows[[p]] <- data.table(
      pathway = p, rank = seq_along(chosen), gene = chosen,
      adjusted_gene_beta = gene_beta[chosen], pathway_fixed_beta = fixed_beta[[p]]
    )
  }
}
driver_table <- rbindlist(driver_rows, fill = TRUE)
fwrite(driver_table, file.path(outdir, "hallmark_top20_driver_genes.tsv"), sep = "\t")

jaccard <- function(a, b) {
  u <- union(a, b)
  if (!length(u)) return(0)
  length(intersect(a, b)) / length(u)
}
pair_idx <- t(combn(pathways, 2L))
overlap_pairs <- data.table(pathway_1 = pair_idx[, 1L], pathway_2 = pair_idx[, 2L])
overlap_pairs[, `:=`(
  gene_jaccard = mapply(function(a, b) jaccard(upper_gene[measured_idx[[a]]],
                                               upper_gene[measured_idx[[b]]]), pathway_1, pathway_2),
  driver_jaccard = mapply(function(a, b) jaccard(driver_sets[[a]], driver_sets[[b]]),
                           pathway_1, pathway_2),
  beta_1 = fixed_beta[pathway_1],
  beta_2 = fixed_beta[pathway_2]
)]
overlap_pairs[, `:=`(
  same_effect_sign = sign(beta_1) == sign(beta_2) & sign(beta_1) != 0,
  absolute_beta_difference = abs(beta_1 - beta_2)
)]
fwrite(overlap_pairs, file.path(outdir, "hallmark_pairwise_overlap.tsv.gz"), sep = "\t")

edges <- overlap_pairs[gene_jaccard >= 0.20 & same_effect_sign == TRUE,
                       .(from = pathway_1, to = pathway_2)]
g <- graph_from_data_frame(edges, directed = FALSE, vertices = data.frame(name = pathways))
membership <- components(g)$membership
component_table <- data.table(pathway = names(membership), redundancy_component = as.integer(membership))
component_summary <- component_table[, {
  members <- pathway
  member_beta <- fixed_beta[members]
  representative <- sort(members[abs(member_beta) == max(abs(member_beta))])[1L]
  all_genes <- lapply(members, function(p) upper_gene[measured_idx[[p]]])
  union_genes <- Reduce(union, all_genes)
  shared_two_plus <- if (length(members) > 1L) {
    tab <- table(unlist(all_genes, use.names = FALSE)); sum(tab >= 2L)
  } else 0L
  list(
    n_hallmarks = length(members),
    members = paste(sort(members), collapse = ";"),
    representative = representative,
    effect_sign = sign(fixed_beta[[representative]]),
    union_gene_count = length(union_genes),
    genes_shared_by_two_or_more = shared_two_plus,
    shared_gene_fraction = if (length(union_genes)) shared_two_plus / length(union_genes) else 0
  )
}, by = redundancy_component]
fwrite(component_table, file.path(outdir, "hallmark_redundancy_membership.tsv"), sep = "\t")
fwrite(component_summary, file.path(outdir, "hallmark_redundancy_components.tsv"), sep = "\t")

driver_recurrence <- driver_table[, .(
  n_hallmarks = uniqueN(pathway),
  hallmark_members = paste(sort(unique(pathway)), collapse = ";"),
  median_abs_adjusted_beta = median(abs(adjusted_gene_beta))
), by = gene][order(-n_hallmarks, -median_abs_adjusted_beta, gene)]
fwrite(driver_recurrence, file.path(outdir, "hallmark_driver_gene_recurrence.tsv"), sep = "\t")

# Matched-size random-gene-set null. Frisch-Waugh-Lovell residualization makes
# the coefficient identical to the full fixed-effects model while avoiding
# 49,950 repeated model-matrix fits.
X_cov <- model.matrix(~ compartment + cancer + meta_state, data = dat)
axis_resid <- lm.fit(X_cov, dat$Axis1_z)$residuals
axis_denom <- sum(axis_resid^2)
if (!is.finite(axis_denom) || axis_denom <= 0) stop("Axis1 residual denominator is invalid")

run_null_pathway <- function(i) {
  p <- pathways[[i]]
  set.seed(seed + i)
  m <- length(measured_idx[[p]])
  rows <- vector("list", n_null)
  for (j in seq_len(n_null)) {
    idx <- sample.int(nrow(gene_z), m, replace = FALSE)
    activity <- colMeans(gene_z[idx, , drop = FALSE])
    activity_z <- zscore(activity)
    beta <- if (all(is.finite(activity_z))) sum(axis_resid * activity_z) / axis_denom else NA_real_
    rows[[j]] <- data.table(pathway = p, iteration = j, set_size = m, null_beta = beta)
  }
  rbindlist(rows)
}
null_list <- if (.Platform$OS.type == "unix" && workers > 1L) {
  mclapply(seq_along(pathways), run_null_pathway, mc.cores = workers,
           mc.preschedule = FALSE, mc.set.seed = FALSE)
} else {
  lapply(seq_along(pathways), run_null_pathway)
}
null_dt <- rbindlist(null_list)
fwrite(null_dt, file.path(outdir, "hallmark_matched_gene_set_null.tsv.gz"), sep = "\t")
null_summary <- null_dt[, .(
  null_evaluable = sum(is.finite(null_beta)),
  null_median = median(null_beta, na.rm = TRUE),
  null_q025 = quantile(null_beta, 0.025, na.rm = TRUE),
  null_q975 = quantile(null_beta, 0.975, na.rm = TRUE)
), by = pathway]
null_summary[, observed_fixed_beta := fixed_beta[pathway]]
null_summary[, empirical_two_sided_p := vapply(pathway, function(p) {
  z <- null_dt[pathway == p & is.finite(null_beta), null_beta]
  (1 + sum(abs(z) >= abs(fixed_beta[[p]]))) / (n_null + 1)
}, numeric(1))]
null_summary[, empirical_q_value := p.adjust(empirical_two_sided_p, method = "BH")]
null_summary[, matched_null_fdr_pass := empirical_q_value < 0.05]
fwrite(null_summary, file.path(outdir, "hallmark_matched_null_summary.tsv"), sep = "\t")

final <- merge(robust_summary, sample_summary[, .(
  pathway, sample_fixed_beta = fixed_beta,
  sample_ci_low = ci_low, sample_ci_high = ci_high,
  sample_evaluable_bootstraps = evaluable_bootstraps,
  sample_direction_supported
)], by = "pathway", all = TRUE, sort = FALSE)
final <- merge(final, null_summary[, .(
  pathway, null_evaluable, null_q025, null_q975,
  empirical_two_sided_p, empirical_q_value, matched_null_fdr_pass
)], by = "pathway", all = TRUE, sort = FALSE)
final <- merge(final, component_table, by = "pathway", all.x = TRUE, sort = FALSE)
final[, concordant_all_new_layers := robust_direction_supported & sample_direction_supported &
        matched_null_fdr_pass & sign(fixed_beta) == sign(sample_fixed_beta)]
setorder(final, empirical_q_value, pathway)
fwrite(final, file.path(outdir, "hallmark_gate12z_calibration_summary.tsv"), sep = "\t")

summary <- data.table(
  metric = c(
    "pseudobulks", "samples", "patients", "hallmarks", "activity_reconstruction_max_error",
    "original_mixed_singular_models", "patient_cluster_supported",
    "sample_aggregated_supported", "matched_null_fdr", "concordant_all_new_layers",
    "redundancy_components", "multi_hallmark_components", "null_sets_per_hallmark"
  ),
  value = as.character(c(
    nrow(dat), nrow(sample_weighted), length(patients), length(pathways), max_activity_error,
    sum(mixed_primary$singular, na.rm = TRUE), sum(final$robust_direction_supported),
    sum(final$sample_direction_supported), sum(final$matched_null_fdr_pass),
    sum(final$concordant_all_new_layers), nrow(component_summary),
    sum(component_summary$n_hallmarks > 1L), n_null
  ))
)
fwrite(summary, file.path(outdir, "GATE12Z_HALLMARK_SUMMARY.tsv"), sep = "\t")

writeLines(c(
  "# Gate12Z Hallmark calibration checkpoint", "",
  paste0("- Frozen Hallmark activity reconstruction maximum error: ", format(max_activity_error, scientific = TRUE)),
  paste0("- Original mixed models flagged singular: ", sum(mixed_primary$singular, na.rm = TRUE), "/50"),
  paste0("- Patient-cluster bootstrap supported: ", sum(final$robust_direction_supported), "/50"),
  paste0("- Sample-aggregated patient-bootstrap supported: ", sum(final$sample_direction_supported), "/50"),
  paste0("- Matched-random-set empirical FDR supported: ", sum(final$matched_null_fdr_pass), "/50"),
  paste0("- Concordant across all three new layers: ", sum(final$concordant_all_new_layers), "/50"),
  paste0("- Gene/effect redundancy components: ", nrow(component_summary),
         " (multi-Hallmark components ", sum(component_summary$n_hallmarks > 1L), ")"),
  "",
  "These analyses recalibrate internal association only. They do not provide independent pathway validation, causal activation, a therapeutic target screen or evidence of 50 independent mechanisms."
), file.path(outdir, "GATE12Z_HALLMARK_CHECKPOINT.md"))

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
message("Gate12Z Hallmark calibration complete")
