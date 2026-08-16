#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(babelgene)
  library(UCell)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("Usage: run_gate12r_spatial_malignant_exclusion.R <transfer_ref.rds> <gate12g_model.rds> <raw_spatial_dir> <section_objects.rds> <outdir>")
}

reference_file <- args[[1L]]
model_file <- args[[2L]]
raw_dir <- args[[3L]]
section_file <- args[[4L]]
outdir <- args[[5L]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

seed <- 1208L
n_perm <- 999L
knn_k <- 6L
set.seed(seed)

transfer <- readRDS(reference_file)
model <- readRDS(model_file)
sections <- readRDS(section_file)
if (!identical(transfer$status, "PASS")) stop("Transfer reference did not PASS")
if (length(sections) != 4L) stop("Expected four spatial sections")

spatial_samples <- data.table(
  sample = c("GSM9564255", "GSM9564256", "GSM9564257", "GSM9564258"),
  library = c("24664", "24665", "24666", "24667")
)

axis_loading <- model$pca$rotation[, 1L]
feature_table <- data.table(
  feature = names(axis_loading), loading = as.numeric(axis_loading), taxonomy = unname(model$blocks)
)
feature_table[, label := sub("^(Broad|Myeloid|T_NK)__", "", feature)]
feature_table[, rejection_bin := label %chin% c("Unassigned", "Unresolved_myeloid", "Unresolved_T_NK")]
feature_table <- merge(
  feature_table,
  transfer$class_metrics[, .(taxonomy, label = truth, evaluable_patients, recall)],
  by = c("taxonomy", "label"), all.x = TRUE
)
feature_table[, reliable := rejection_bin | (!is.na(recall) & evaluable_patients >= 3L & recall >= 0.50)]
feature_table[, marker_eligible := reliable & !rejection_bin & !is.na(recall)]

build_human_signature <- function(features) {
  marker_rows <- list()
  for (i in which(features$marker_eligible)) {
    tax <- features$taxonomy[i]
    lab <- features$label[i]
    loading_value <- features$loading[i]
    ref <- transfer$references[[tax]]
    labels <- ref$meta$label
    if (!lab %chin% labels) next
    target_mean <- rowMeans(ref$logcpm[, labels == lab, drop = FALSE])
    other_labels <- setdiff(unique(labels), lab)
    other_means <- sapply(other_labels, function(z) rowMeans(ref$logcpm[, labels == z, drop = FALSE]))
    if (is.null(dim(other_means))) other_means <- matrix(other_means, ncol = 1L)
    effect <- target_mean - apply(other_means, 1L, max)
    ok <- is.finite(effect) & effect > 0 & target_mean > log1p(1)
    genes <- names(effect)[ok][order(effect[ok], decreasing = TRUE)]
    genes <- genes[seq_len(min(30L, length(genes)))]
    effects <- effect[genes]
    if (!length(genes) || sum(effects) <= 0) next
    marker_rows[[length(marker_rows) + 1L]] <- data.table(
      feature = features$feature[i], taxonomy = tax, label = lab,
      feature_loading = loading_value, human_gene = genes,
      specificity_effect = effects, coefficient = loading_value * effects / sum(effects)
    )
  }
  detail <- rbindlist(marker_rows)
  signature <- detail[, .(coefficient = sum(coefficient), contributing_features = uniqueN(feature)), by = human_gene]
  list(detail = detail, signature = signature[coefficient != 0])
}

map_mouse_signature <- function(human_signature) {
  orth <- as.data.table(babelgene::orthologs(human_signature$human_gene, species = "mouse"))
  orth <- orth[human_symbol %chin% human_signature$human_gene]
  orth[, human_mapping_n := uniqueN(symbol), by = human_symbol]
  orth[, mouse_mapping_n := uniqueN(human_symbol), by = symbol]
  orth <- orth[human_mapping_n == 1L & mouse_mapping_n == 1L]
  mouse <- merge(
    human_signature,
    orth[, .(human_gene = human_symbol, mouse_gene = symbol, mouse_ensembl = ensembl, support_n)],
    by = "human_gene", all = FALSE
  )
  mouse <- mouse[, .(
    coefficient = sum(coefficient), human_sources = uniqueN(human_gene), support_n = max(support_n)
  ), by = .(mouse_gene, mouse_ensembl)]
  mouse[coefficient != 0]
}

full_human <- build_human_signature(feature_table)
excluded_feature <- "Broad__Malignant"
if (!excluded_feature %chin% feature_table$feature) stop("Frozen malignant feature was not found")
excluded_features <- feature_table[feature != excluded_feature]
excluded_human <- build_human_signature(excluded_features)
full_mouse <- map_mouse_signature(full_human$signature)
excluded_mouse <- map_mouse_signature(excluded_human$signature)

mass <- feature_table[, .(
  n_features = .N,
  absolute_loading_mass = sum(abs(loading)),
  squared_loading_mass = sum(loading^2)
), by = taxonomy]
mass[, `:=`(
  absolute_fraction_of_total = absolute_loading_mass / sum(absolute_loading_mass),
  squared_fraction_of_total = squared_loading_mass / sum(squared_loading_mass)
)]
malignant_mass <- feature_table[feature == excluded_feature, .(
  feature, taxonomy, loading,
  absolute_loading_mass = abs(loading), squared_loading_mass = loading^2,
  absolute_fraction_of_total = abs(loading) / sum(abs(feature_table$loading)),
  squared_fraction_of_total = loading^2 / sum(feature_table$loading^2),
  absolute_fraction_of_broad = abs(loading) / sum(abs(feature_table[taxonomy == "Broad", loading])),
  squared_fraction_of_broad = loading^2 / sum(feature_table[taxonomy == "Broad", loading]^2)
)]

loading_coverage <- rbindlist(list(
  data.table(version = "full", feature_count = nrow(feature_table),
             represented_feature_count = sum(feature_table$marker_eligible),
             represented_absolute_loading_fraction = sum(abs(feature_table[marker_eligible == TRUE, loading])) / sum(abs(feature_table$loading)),
             represented_squared_loading_fraction = sum(feature_table[marker_eligible == TRUE, loading]^2) / sum(feature_table$loading^2)),
  data.table(version = "malignant_excluded", feature_count = nrow(excluded_features),
             represented_feature_count = sum(excluded_features$marker_eligible),
             represented_absolute_loading_fraction = sum(abs(excluded_features[marker_eligible == TRUE, loading])) / sum(abs(excluded_features$loading)),
             represented_squared_loading_fraction = sum(excluded_features[marker_eligible == TRUE, loading]^2) / sum(excluded_features$loading^2))
))

ortholog_audit <- rbindlist(list(
  data.table(version = "full", human_genes = nrow(full_human$signature), mouse_one_to_one_genes = nrow(full_mouse),
             positive_mouse_genes = sum(full_mouse$coefficient > 0), negative_mouse_genes = sum(full_mouse$coefficient < 0),
             ortholog_absolute_coefficient_coverage = sum(abs(full_mouse$coefficient)) / sum(abs(full_human$signature$coefficient))),
  data.table(version = "malignant_excluded", human_genes = nrow(excluded_human$signature), mouse_one_to_one_genes = nrow(excluded_mouse),
             positive_mouse_genes = sum(excluded_mouse$coefficient > 0), negative_mouse_genes = sum(excluded_mouse$coefficient < 0),
             ortholog_absolute_coefficient_coverage = sum(abs(excluded_mouse$coefficient)) / sum(abs(excluded_human$signature$coefficient)))
))

read_gz_lines <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  readLines(con, warn = FALSE)
}

read_spatial_counts <- function(sample, library, retained_barcodes) {
  prefix <- file.path(raw_dir, paste0(sample, "_Bone_ST_filtered_", library))
  con <- gzfile(paste0(prefix, "_matrix.mtx.gz"), "rb")
  on.exit(close(con))
  counts <- as(readMM(con), "dgCMatrix")
  ff <- strsplit(read_gz_lines(paste0(prefix, "_features.tsv.gz")), "\t", fixed = TRUE)
  feature_col <- if (all(lengths(ff) >= 2L)) 2L else 1L
  features <- vapply(ff, `[[`, character(1), feature_col)
  barcodes <- vapply(strsplit(read_gz_lines(paste0(prefix, "_barcodes.tsv.gz")), "\t", fixed = TRUE),
                     `[[`, character(1), 1L)
  if (nrow(counts) != length(features) || ncol(counts) != length(barcodes)) stop("10x mismatch for ", sample)
  rownames(counts) <- make.unique(features)
  colnames(counts) <- barcodes
  idx <- match(retained_barcodes, barcodes)
  if (anyNA(idx)) stop("Retained spots missing from raw counts for ", sample)
  counts[, idx, drop = FALSE]
}

score_signature <- function(counts, signature) {
  libsize <- Matrix::colSums(counts)
  gene_idx <- match(signature$mouse_gene, rownames(counts), nomatch = 0L)
  present <- which(gene_idx > 0L)
  sig <- signature[present]
  expr <- as.matrix(counts[gene_idx[present], , drop = FALSE])
  expr <- log1p(t(t(expr) / libsize) * 1e4)
  rownames(expr) <- sig$mouse_gene
  gene_sd <- apply(expr, 1L, sd)
  variable <- is.finite(gene_sd) & gene_sd > 0
  expr <- expr[variable, , drop = FALSE]
  sig <- sig[variable]
  z <- t(scale(t(expr)))
  weighted <- as.numeric(crossprod(sig$coefficient, z)) / sum(abs(sig$coefficient))
  positive <- sig[coefficient > 0, mouse_gene]
  negative <- sig[coefficient < 0, mouse_gene]
  ucell <- UCell::ScoreSignatures_UCell(
    matrix = counts,
    features = list(Axis1_positive = positive, Axis1_negative = negative),
    maxRank = 1500, ncores = 1, force.gc = TRUE
  )
  ucell <- as.data.frame(ucell)
  list(
    weighted = weighted,
    signed_ucell = ucell$Axis1_positive_UCell - ucell$Axis1_negative_UCell,
    genes_present = length(present), variable_genes = nrow(sig),
    absolute_coefficient_coverage = sum(abs(sig$coefficient)) / sum(abs(signature$coefficient)),
    positive_genes = sum(sig$coefficient > 0), negative_genes = sum(sig$coefficient < 0)
  )
}

sym_knn_graph <- function(coords, k = 6L) {
  d <- as.matrix(dist(coords))
  diag(d) <- Inf
  a <- matrix(0, nrow(d), nrow(d))
  for (i in seq_len(nrow(d))) a[i, order(d[i, ])[seq_len(k)]] <- 1
  a <- (a + t(a)) > 0
  diag(a) <- FALSE
  Matrix(a * 1, sparse = TRUE)
}

moran_stat <- function(x, adjacency) {
  xc <- x - mean(x)
  as.numeric((length(x) / sum(adjacency)) * (crossprod(xc, adjacency %*% xc) / crossprod(xc)))
}

rank_biserial <- function(x1, x0) {
  if (!length(x1) || !length(x0)) return(NA_real_)
  w <- suppressWarnings(wilcox.test(x1, x0, exact = FALSE)$statistic)
  as.numeric(2 * w / (length(x1) * length(x0)) - 1)
}

broad_mapping <- list(
  T_NK = "T cell", Myeloid = c("Myeloid", "Macrophage"),
  B = c("B cell", "plasma B cell", "pre-B cell"), Stromal = "MSCs",
  Endothelial = "Endothelial", Erythroid = "Erythroid"
)
broad_load <- feature_table[taxonomy == "Broad" & label != "Malignant", setNames(loading, label)]

spot_rows <- list()
metric_rows <- list()
for (ii in seq_len(nrow(spatial_samples))) {
  sample <- spatial_samples$sample[ii]
  library <- spatial_samples$library[ii]
  section <- sections[[sample]]
  barcodes <- rownames(section$weights)
  counts <- read_spatial_counts(sample, library, barcodes)
  full_score <- score_signature(counts, full_mouse)
  excl_score <- score_signature(counts, excluded_mouse)

  coords <- as.matrix(section$coords[barcodes, c("x", "y"), drop = FALSE])
  adjacency <- sym_knn_graph(coords, knn_k)
  observed_moran <- moran_stat(excl_score$weighted, adjacency)
  null_moran <- replicate(n_perm, moran_stat(sample(excl_score$weighted), adjacency))
  p_moran <- (1 + sum(null_moran >= observed_moran)) / (n_perm + 1)

  proxy <- rep(0, nrow(section$weights))
  represented_broad <- character()
  for (lab in names(broad_mapping)) {
    celltypes <- intersect(broad_mapping[[lab]], colnames(section$weights))
    if (!length(celltypes) || !lab %chin% names(broad_load)) next
    proxy <- proxy + broad_load[[lab]] * rowSums(section$weights[, celltypes, drop = FALSE])
    represented_broad <- c(represented_broad, lab)
  }
  observed_rho <- suppressWarnings(cor(excl_score$weighted, proxy, method = "spearman"))
  null_rho <- replicate(n_perm, suppressWarnings(cor(excl_score$weighted, sample(proxy), method = "spearman")))
  p_rho <- (1 + sum(null_rho >= observed_rho)) / (n_perm + 1)

  submitted <- as.character(section$submitted_label[barcodes])
  tumor <- submitted == "Tumor"
  adjacent <- !tumor & as.logical(adjacency %*% as.numeric(tumor) > 0)
  metric_rows[[ii]] <- data.table(
    sample = sample, spots = length(barcodes),
    score_rho_full_vs_malignant_excluded = cor(full_score$weighted, excl_score$weighted, method = "spearman"),
    malignant_excluded_ucell_rho = cor(excl_score$weighted, excl_score$signed_ucell, method = "spearman"),
    malignant_excluded_moran_i = observed_moran, moran_p_one_sided = p_moran,
    malignant_excluded_rctd_rho = observed_rho, rctd_p_one_sided = p_rho,
    tumor_spots = sum(tumor), tumor_vs_non_tumor_rb = rank_biserial(excl_score$weighted[tumor], excl_score$weighted[!tumor]),
    adjacent_non_tumor_spots = sum(adjacent),
    adjacent_vs_other_non_tumor_rb = rank_biserial(excl_score$weighted[adjacent], excl_score$weighted[!tumor & !adjacent]),
    signature_genes_present = excl_score$genes_present,
    variable_signature_genes = excl_score$variable_genes,
    absolute_coefficient_coverage = excl_score$absolute_coefficient_coverage,
    positive_genes_present = excl_score$positive_genes, negative_genes_present = excl_score$negative_genes,
    rctd_broad_absolute_loading_coverage = sum(abs(broad_load[represented_broad])) / sum(abs(broad_load)),
    rctd_broad_squared_loading_coverage = sum(broad_load[represented_broad]^2) / sum(broad_load^2)
  )
  spot_rows[[ii]] <- data.table(
    sample = sample, barcode = barcodes,
    x = section$coords[barcodes, "x"], y = section$coords[barcodes, "y"], submitted_label = submitted,
    axis1_full_recalculated = full_score$weighted,
    axis1_malignant_excluded = excl_score$weighted,
    malignant_excluded_signed_ucell = excl_score$signed_ucell,
    malignant_excluded_rctd_proxy = proxy
  )
}

metrics <- rbindlist(metric_rows)
metrics[, moran_q_bh := p.adjust(moran_p_one_sided, "BH")]
metrics[, rctd_q_bh := p.adjust(rctd_p_one_sided, "BH")]
spots <- rbindlist(spot_rows)

decision <- data.table(
  sections = nrow(metrics),
  all_full_vs_excluded_rho_ge_0_85 = all(metrics$score_rho_full_vs_malignant_excluded >= 0.85),
  moran_positive_q_lt_0_10_sections = sum(metrics$malignant_excluded_moran_i > 0 & metrics$moran_q_bh < 0.10),
  rctd_positive_q_lt_0_10_sections = sum(metrics$malignant_excluded_rctd_rho > 0 & metrics$rctd_q_bh < 0.10),
  all_sections_signature_coverage_ge_0_70 = all(metrics$absolute_coefficient_coverage >= 0.70),
  interpretation = "linked_descriptive_sensitivity_only",
  status = "COMPLETE"
)

fwrite(feature_table, file.path(outdir, "axis1_feature_reliability.tsv"), sep = "\t")
fwrite(mass, file.path(outdir, "axis1_block_loading_mass.tsv"), sep = "\t")
fwrite(malignant_mass, file.path(outdir, "axis1_malignant_loading_mass.tsv"), sep = "\t")
fwrite(loading_coverage, file.path(outdir, "axis1_represented_loading_coverage.tsv"), sep = "\t")
fwrite(ortholog_audit, file.path(outdir, "axis1_ortholog_audit_full_vs_excluded.tsv"), sep = "\t")
fwrite(excluded_human$detail, file.path(outdir, "axis1_malignant_excluded_human_marker_detail.tsv"), sep = "\t")
fwrite(excluded_human$signature, file.path(outdir, "axis1_malignant_excluded_human_gene_coefficients.tsv"), sep = "\t")
fwrite(excluded_mouse, file.path(outdir, "axis1_malignant_excluded_mouse_gene_coefficients.tsv"), sep = "\t")
fwrite(spots, file.path(outdir, "axis1_malignant_exclusion_spot_scores.tsv.gz"), sep = "\t")
fwrite(metrics, file.path(outdir, "axis1_malignant_exclusion_section_metrics.tsv"), sep = "\t")
fwrite(decision, file.path(outdir, "gate12r_spatial_malignant_exclusion_decision.tsv"), sep = "\t")

p1 <- ggplot(spots, aes(axis1_full_recalculated, axis1_malignant_excluded, colour = sample)) +
  geom_point(size = 0.35, alpha = 0.35) + theme_bw(base_size = 10) +
  labs(title = "A  Full versus malignant-excluded score", x = "Full recalculated Axis1", y = "Malignant-excluded Axis1")
p2 <- ggplot(spots, aes(x, y, colour = axis1_malignant_excluded)) +
  geom_point(size = 0.48) + scale_colour_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
  coord_equal() + facet_wrap(~sample) + theme_void(base_size = 10) +
  labs(title = "B  Malignant-excluded score at measured spots", colour = "Score")
p3 <- ggplot(metrics, aes(sample, malignant_excluded_moran_i)) +
  geom_col(fill = "#0072B2") + geom_hline(yintercept = 0, colour = "grey45") +
  theme_bw(base_size = 10) + labs(title = "C  Within-section Moran's I", x = NULL, y = "Moran's I")
p4 <- ggplot(metrics, aes(sample, malignant_excluded_rctd_rho)) +
  geom_col(fill = "#009E73") + geom_hline(yintercept = 0, colour = "grey45") +
  theme_bw(base_size = 10) + labs(title = "D  RCTD-proxy concordance", x = NULL, y = "Spearman rho")
fig <- p1 / p2 / (p3 | p4) + plot_annotation(
  title = "Gate12R spatial sensitivity after removing Broad__Malignant from frozen Axis1"
)
ggsave(file.path(outdir, "FigureR3_spatial_malignant_exclusion.pdf"), fig, width = 13, height = 12)
ggsave(file.path(outdir, "FigureR3_spatial_malignant_exclusion.png"), fig, width = 13, height = 12, dpi = 300, bg = "white")

writeLines(c(
  "# Gate12R spatial malignant-exclusion checkpoint", "",
  paste0("- Malignant absolute loading fraction of the full Axis1: ", sprintf("%.4f", malignant_mass$absolute_fraction_of_total)),
  paste0("- Malignant squared loading fraction of the full Axis1: ", sprintf("%.4f", malignant_mass$squared_fraction_of_total)),
  paste0("- Full-versus-excluded score rho range: ", sprintf("%.4f-%.4f", min(metrics$score_rho_full_vs_malignant_excluded), max(metrics$score_rho_full_vs_malignant_excluded))),
  paste0("- Positive Moran sections at BH q<0.10: ", decision$moran_positive_q_lt_0_10_sections, "/4"),
  paste0("- Positive RCTD-proxy sections at BH q<0.10: ", decision$rctd_positive_q_lt_0_10_sections, "/4"),
  "- All metrics are section-specific. Sections are not treated as independent animals.",
  "- This remains a linked, cross-species, descriptive sensitivity analysis and is not independent human spatial validation."
), file.path(outdir, "GATE12R_SPATIAL_MALIGNANT_EXCLUSION_CHECKPOINT.md"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
cat("GATE12R_SPATIAL_MALIGNANT_EXCLUSION_STATUS=COMPLETE\n")
print(decision)
