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
  stop("Usage: run_gate12g_spatial_axis1.R <transfer_ref.rds> <gate12g_model.rds> <raw_spatial_dir> <section_objects.rds> <outdir>")
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
if (!1L %in% model$axis_summary[eligible == TRUE, axis]) stop("Axis1 was not discovery-eligible")
if (length(sections) != 4L) stop("Expected four spatial sections")

spatial_samples <- data.table(
  sample = c("GSM9564255", "GSM9564256", "GSM9564257", "GSM9564258"),
  library = c("24664", "24665", "24666", "24667")
)

axis_loading <- model$pca$rotation[, 1L]
feature_table <- data.table(feature = names(axis_loading), loading = as.numeric(axis_loading),
                            taxonomy = model$blocks)
feature_table[, label := sub("^(Broad|Myeloid|T_NK)__", "", feature)]
feature_table[, rejection_bin := label %chin% c("Unassigned", "Unresolved_myeloid", "Unresolved_T_NK")]
class_metrics <- transfer$class_metrics
feature_table <- merge(feature_table,
  class_metrics[, .(taxonomy, label = truth, evaluable_patients, recall)],
  by = c("taxonomy", "label"), all.x = TRUE)
feature_table[, reliable := rejection_bin |
                (!is.na(recall) & evaluable_patients >= 3L & recall >= 0.50)]
feature_table[, marker_eligible := reliable & !rejection_bin & !is.na(recall)]

marker_rows <- list()
for (i in which(feature_table$marker_eligible)) {
  tax <- feature_table$taxonomy[i]
  lab <- feature_table$label[i]
  loading_value <- feature_table$loading[i]
  ref <- transfer$references[[tax]]
  labels <- ref$meta$label
  if (!lab %chin% labels) next
  target_mean <- rowMeans(ref$logcpm[, labels == lab, drop = FALSE])
  other_labels <- setdiff(unique(labels), lab)
  other_means <- sapply(other_labels, function(z) rowMeans(ref$logcpm[, labels == z, drop = FALSE]))
  if (is.null(dim(other_means))) other_means <- matrix(other_means, ncol = 1L)
  effect <- target_mean - apply(other_means, 1L, max)
  ok <- is.finite(effect) & effect > 0 & target_mean > log1p(1)
  ord <- order(effect[ok], decreasing = TRUE)
  genes <- names(effect)[ok][ord][seq_len(min(30L, sum(ok)))]
  effects <- effect[genes]
  if (!length(genes) || sum(effects) <= 0) next
  marker_rows[[length(marker_rows) + 1L]] <- data.table(
    feature = feature_table$feature[i], taxonomy = tax, label = lab,
    feature_loading = loading_value, human_gene = genes,
    specificity_effect = effects,
    coefficient = loading_value * effects / sum(effects)
  )
}
human_detail <- rbindlist(marker_rows)
human_signature <- human_detail[, .(coefficient = sum(coefficient),
                                    contributing_features = uniqueN(feature)), by = human_gene]
human_signature <- human_signature[coefficient != 0]

marker_feature_coverage <- feature_table[, sum(loading[marker_eligible]^2) / sum(loading^2)]

orth <- as.data.table(babelgene::orthologs(human_signature$human_gene, species = "mouse"))
orth <- orth[human_symbol %chin% human_signature$human_gene]
orth[, human_mapping_n := uniqueN(symbol), by = human_symbol]
orth[, mouse_mapping_n := uniqueN(human_symbol), by = symbol]
orth <- orth[human_mapping_n == 1L & mouse_mapping_n == 1L]
mouse_signature <- merge(human_signature, orth[, .(human_gene = human_symbol, mouse_gene = symbol,
                                                   mouse_ensembl = ensembl, support_n)],
                         by = "human_gene", all = FALSE)
mouse_signature <- mouse_signature[, .(
  coefficient = sum(coefficient),
  human_sources = uniqueN(human_gene),
  support_n = max(support_n)
), by = .(mouse_gene, mouse_ensembl)]
mouse_signature <- mouse_signature[coefficient != 0]

ortholog_abs_coverage <- sum(abs(mouse_signature$coefficient)) / sum(abs(human_signature$coefficient))
mapping_audit <- data.table(
  metric = c("marker_feature_loading_sq_coverage", "human_signature_genes",
             "mouse_one_to_one_genes", "positive_mouse_genes", "negative_mouse_genes",
             "ortholog_absolute_coefficient_coverage"),
  value = c(marker_feature_coverage, nrow(human_signature), nrow(mouse_signature),
            sum(mouse_signature$coefficient > 0), sum(mouse_signature$coefficient < 0),
            ortholog_abs_coverage)
)
mapping_pass <- nrow(mouse_signature) >= 50L &
  sum(mouse_signature$coefficient > 0) >= 15L &
  sum(mouse_signature$coefficient < 0) >= 15L & ortholog_abs_coverage >= 0.70
if (!mapping_pass) {
  fwrite(mapping_audit, file.path(outdir, "axis1_mapping_audit.tsv"), sep = "\t")
  stop("Frozen Axis1 orthologue mapping failed thresholds")
}

read_gz_lines <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con)); readLines(con, warn = FALSE)
}

read_spatial_counts <- function(sample, library, retained_barcodes) {
  prefix <- file.path(raw_dir, paste0(sample, "_Bone_ST_filtered_", library))
  con <- gzfile(paste0(prefix, "_matrix.mtx.gz"), "rb"); on.exit(close(con))
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

sym_knn_graph <- function(coords, k = 6L) {
  d <- as.matrix(dist(coords))
  diag(d) <- Inf
  n <- nrow(d)
  a <- matrix(0, n, n)
  for (i in seq_len(n)) a[i, order(d[i, ])[seq_len(k)]] <- 1
  a <- (a + t(a)) > 0
  diag(a) <- FALSE
  Matrix(a * 1, sparse = TRUE)
}

moran_stat <- function(x, adjacency) {
  n <- length(x); xc <- x - mean(x); s0 <- sum(adjacency)
  as.numeric((n / s0) * (crossprod(xc, adjacency %*% xc) / crossprod(xc)))
}

rank_biserial <- function(x1, x0) {
  if (!length(x1) || !length(x0)) return(NA_real_)
  w <- suppressWarnings(wilcox.test(x1, x0, exact = FALSE)$statistic)
  as.numeric(2 * w / (length(x1) * length(x0)) - 1)
}

broad_mapping <- list(
  T_NK = "T cell",
  Myeloid = c("Myeloid", "Macrophage"),
  B = c("B cell", "plasma B cell", "pre-B cell"),
  Stromal = "MSCs",
  Endothelial = "Endothelial",
  Erythroid = "Erythroid"
)

spot_rows <- list(); moran_rows <- list(); rctd_rows <- list(); tumor_rows <- list(); section_cov <- list()
positive_genes <- mouse_signature[coefficient > 0, mouse_gene]
negative_genes <- mouse_signature[coefficient < 0, mouse_gene]

for (ii in seq_len(nrow(spatial_samples))) {
  sample <- spatial_samples$sample[ii]
  library <- spatial_samples$library[ii]
  section <- sections[[sample]]
  barcodes <- rownames(section$weights)
  counts <- read_spatial_counts(sample, library, barcodes)
  libsize <- Matrix::colSums(counts)
  gene_idx <- match(mouse_signature$mouse_gene, rownames(counts), nomatch = 0L)
  present <- which(gene_idx > 0L)
  sig <- mouse_signature[present]
  expr <- as.matrix(counts[gene_idx[present], , drop = FALSE])
  expr <- log1p(t(t(expr) / libsize) * 1e4)
  rownames(expr) <- sig$mouse_gene
  gene_sd <- apply(expr, 1L, sd)
  variable <- is.finite(gene_sd) & gene_sd > 0
  expr <- expr[variable, , drop = FALSE]
  sig <- sig[variable]
  z <- t(scale(t(expr)))
  weighted_score <- as.numeric(crossprod(sig$coefficient, z)) / sum(abs(sig$coefficient))

  ucell <- UCell::ScoreSignatures_UCell(
    matrix = counts,
    features = list(Axis1_positive = intersect(positive_genes, rownames(counts)),
                    Axis1_negative = intersect(negative_genes, rownames(counts))),
    maxRank = 1500, ncores = 1, force.gc = TRUE
  )
  ucell <- as.data.frame(ucell)
  signed_ucell <- ucell$Axis1_positive_UCell - ucell$Axis1_negative_UCell

  coords <- as.matrix(section$coords[barcodes, c("x", "y"), drop = FALSE])
  adjacency <- sym_knn_graph(coords, knn_k)
  observed_moran <- moran_stat(weighted_score, adjacency)
  null_moran <- replicate(n_perm, moran_stat(sample(weighted_score), adjacency))
  p_moran <- (1 + sum(null_moran >= observed_moran)) / (n_perm + 1)

  broad_load <- feature_table[taxonomy == "Broad", setNames(loading, label)]
  proxy <- rep(0, nrow(section$weights))
  represented_broad <- character()
  for (lab in names(broad_mapping)) {
    celltypes <- intersect(broad_mapping[[lab]], colnames(section$weights))
    if (!length(celltypes) || !lab %chin% names(broad_load)) next
    proxy <- proxy + broad_load[[lab]] * rowSums(section$weights[, celltypes, drop = FALSE])
    represented_broad <- c(represented_broad, lab)
  }
  observed_rho <- suppressWarnings(cor(weighted_score, proxy, method = "spearman"))
  null_rho <- replicate(n_perm, suppressWarnings(cor(weighted_score, sample(proxy), method = "spearman")))
  p_rho <- (1 + sum(null_rho >= observed_rho)) / (n_perm + 1)

  submitted <- as.character(section$submitted_label[barcodes])
  tumor <- submitted == "Tumor"
  adjacent <- !tumor & as.logical(adjacency %*% as.numeric(tumor) > 0)

  coeff_coverage <- sum(abs(sig$coefficient)) / sum(abs(mouse_signature$coefficient))
  ucell_cor <- suppressWarnings(cor(weighted_score, signed_ucell, method = "spearman"))
  broad_sq_coverage <- sum(broad_load[represented_broad]^2) / sum(broad_load^2)

  section_cov[[ii]] <- data.table(sample = sample, spots = length(barcodes),
    signature_genes_present = length(present), variable_signature_genes = nrow(sig),
    absolute_coefficient_coverage = coeff_coverage,
    positive_genes_present = sum(sig$coefficient > 0), negative_genes_present = sum(sig$coefficient < 0),
    ucell_score_correlation = ucell_cor, rctd_broad_loading_sq_coverage = broad_sq_coverage)
  moran_rows[[ii]] <- data.table(sample = sample, moran_i = observed_moran,
                                 p_one_sided = p_moran, permutations = n_perm)
  rctd_rows[[ii]] <- data.table(sample = sample, spearman_rho = observed_rho,
                                p_one_sided = p_rho, permutations = n_perm)
  tumor_rows[[ii]] <- data.table(sample = sample, tumor_spots = sum(tumor),
    tumor_vs_non_tumor_rb = rank_biserial(weighted_score[tumor], weighted_score[!tumor]),
    adjacent_non_tumor_spots = sum(adjacent),
    adjacent_vs_other_non_tumor_rb = rank_biserial(weighted_score[adjacent], weighted_score[!tumor & !adjacent]))
  spot_rows[[ii]] <- data.table(sample = sample, barcode = barcodes,
    x = section$coords[barcodes, "x"], y = section$coords[barcodes, "y"],
    submitted_label = submitted, axis1_weighted_score = weighted_score,
    axis1_signed_ucell = signed_ucell, rctd_axis1_proxy = proxy,
    macrophage_weight = section$weights[barcodes, "Macrophage"],
    tcell_weight = section$weights[barcodes, "T cell"])
}

section_coverage <- rbindlist(section_cov)
moran <- rbindlist(moran_rows); moran[, q_bh := p.adjust(p_one_sided, "BH")]
moran[, pass := moran_i > 0 & q_bh < 0.10]
rctd <- rbindlist(rctd_rows); rctd[, q_bh := p.adjust(p_one_sided, "BH")]
rctd[, pass := spearman_rho > 0 & q_bh < 0.10]
tumor <- rbindlist(tumor_rows)
spots <- rbindlist(spot_rows)

loo <- rbindlist(lapply(unique(rctd$sample), function(omit) {
  z <- rctd[sample != omit]
  data.table(omitted_section = omit, sections_remaining = nrow(z),
             median_rho = median(z$spearman_rho), all_positive = all(z$spearman_rho > 0))
}))

mapping_gate <- mapping_pass & all(section_coverage$absolute_coefficient_coverage >= 0.70)
ucell_gate <- sum(section_coverage$ucell_score_correlation >= 0.60) >= 3L
moran_gate <- sum(moran$pass) >= 3L
rctd_gate <- sum(rctd$pass) >= 3L
loo_gate <- all(loo$median_rho > 0 & loo$all_positive)
overall <- if (mapping_gate && ucell_gate && moran_gate && rctd_gate && loo_gate) "PASS" else "STOP"

decision <- data.table(
  mapping_pass = mapping_gate,
  ucell_concordant_sections = sum(section_coverage$ucell_score_correlation >= 0.60),
  ucell_pass = ucell_gate,
  moran_pass_sections = sum(moran$pass), moran_pass = moran_gate,
  rctd_pass_sections = sum(rctd$pass), rctd_pass = rctd_gate,
  loo_positive = loo_gate, overall_status = overall
)

fwrite(feature_table, file.path(outdir, "axis1_feature_reliability.tsv"), sep = "\t")
fwrite(human_detail, file.path(outdir, "axis1_human_marker_detail.tsv"), sep = "\t")
fwrite(human_signature, file.path(outdir, "axis1_human_gene_coefficients.tsv"), sep = "\t")
fwrite(mouse_signature, file.path(outdir, "axis1_mouse_gene_coefficients.tsv"), sep = "\t")
fwrite(mapping_audit, file.path(outdir, "axis1_mapping_audit.tsv"), sep = "\t")
fwrite(section_coverage, file.path(outdir, "section_signature_coverage.tsv"), sep = "\t")
fwrite(spots, file.path(outdir, "axis1_spot_scores.tsv.gz"), sep = "\t")
fwrite(moran, file.path(outdir, "axis1_spatial_moran.tsv"), sep = "\t")
fwrite(rctd, file.path(outdir, "axis1_rctd_proxy_correlation.tsv"), sep = "\t")
fwrite(tumor, file.path(outdir, "axis1_tumor_descriptive.tsv"), sep = "\t")
fwrite(loo, file.path(outdir, "axis1_leave_one_section_out.tsv"), sep = "\t")
fwrite(decision, file.path(outdir, "axis1_spatial_decision.tsv"), sep = "\t")
saveRDS(list(mapping_audit = mapping_audit, section_coverage = section_coverage,
             moran = moran, rctd = rctd, tumor = tumor, loo = loo,
             decision = decision, seed = seed, permutations = n_perm),
        file.path(outdir, "gate12g_spatial_axis1_results.rds"), compress = "xz")

p1 <- ggplot(spots, aes(x, y, colour = axis1_weighted_score)) + geom_point(size = 0.5) +
  scale_colour_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
  coord_equal() + facet_wrap(~sample, scales = "free") + theme_void(base_size = 10) +
  labs(title = "A  Frozen Axis1 weighted spatial score", colour = "Axis1")
p2 <- ggplot(spots, aes(axis1_signed_ucell, axis1_weighted_score, colour = sample)) +
  geom_point(size = 0.35, alpha = 0.35) + theme_bw(base_size = 10) +
  labs(title = "B  Rank-based score sensitivity", x = "Signed UCell", y = "Weighted Axis1")
p3 <- ggplot(rctd, aes(sample, spearman_rho)) + geom_col(fill = "#009E73") +
  geom_hline(yintercept = 0, colour = "grey40") + theme_bw(base_size = 10) +
  labs(title = "C  Axis1 versus RCTD ecological proxy", x = NULL, y = "Spearman rho")
p4 <- ggplot(tumor, aes(sample, tumor_vs_non_tumor_rb)) + geom_col(fill = "#CC79A7") +
  geom_hline(yintercept = 0, colour = "grey40") + theme_bw(base_size = 10) +
  labs(title = "D  Tumour-label effect (descriptive)", x = NULL, y = "Rank-biserial effect")
fig <- (p1 / (p2 | p3 | p4)) + plot_annotation(
  title = paste0("Gate12G cross-species spatial projection of externally supported Axis1: ", overall))
ggsave(file.path(outdir, "FigureS_gate12g_spatial_axis1.pdf"), fig, width = 14, height = 10)
ggsave(file.path(outdir, "FigureS_gate12g_spatial_axis1.png"), fig, width = 14, height = 10,
       dpi = 300, bg = "white")

checkpoint <- c(
  "# Gate12G spatial Axis1 checkpoint", "",
  paste0("- Mouse one-to-one signature genes: ", nrow(mouse_signature)),
  paste0("- Orthologue absolute-coefficient coverage: ", sprintf("%.4f", ortholog_abs_coverage)),
  paste0("- UCell-concordant sections: ", decision$ucell_concordant_sections, "/4"),
  paste0("- Moran-pass sections: ", decision$moran_pass_sections, "/4"),
  paste0("- RCTD-proxy-pass sections: ", decision$rctd_pass_sections, "/4"),
  paste0("- Leave-one-section-out positive: ", decision$loo_positive),
  paste0("- Overall decision: **", overall, "**"), "",
  "This is cross-species, section-level spatial support. It is not human spatial validation, animal-level replication or causal evidence."
)
writeLines(checkpoint, file.path(outdir, "GATE12G_SPATIAL_AXIS1_CHECKPOINT.md"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
cat("GATE12G_SPATIAL_AXIS1_STATUS=", overall, "\n", sep = "")
print(decision)
