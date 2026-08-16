#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(Seurat)
  library(harmony)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(ggrastr)
})

args <- commandArgs(trailingOnly = TRUE)
project <- if (length(args) >= 1L) normalizePath(args[[1]], mustWork = TRUE) else
  "."
out <- if (length(args) >= 2L) args[[2]] else file.path(project, "results", "gate12b_cell_states")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

seed <- 20260808L
set.seed(seed)
sce_dir <- file.path(project, "data", "gate3b_work", "annotated_sce")
files <- sort(list.files(sce_dir, pattern = "\\.rds$", full.names = TRUE))
if (length(files) != 42L) stop("Expected 42 annotated SCE files; found ", length(files))

module_list <- list(
  Myeloid = list(
    Classical_monocyte = c("S100A8", "S100A9", "FCN1", "VCAN", "CCR2", "CTSS"),
    Inflammatory_monocyte = c("IL1B", "CXCL2", "CXCL3", "CCL20", "OSM", "NLRP3"),
    C1QC_macrophage = c("C1QA", "C1QB", "C1QC", "APOE", "LPL", "CTSD"),
    SPP1_macrophage = c("SPP1", "APOC1", "LGALS3", "TREM2", "CTSB", "GPNMB"),
    Resident_macrophage = c("SEPP1", "MRC1", "FOLR2", "MARCO", "CD163", "C1QC"),
    cDC = c("FCER1A", "CD1C", "CLEC10A", "HLA-DRA", "CST3"),
    pDC = c("GZMB", "JCHAIN", "IL3RA", "GZMB", "TCF4"),
    Proliferating_myeloid = c("MKI67", "TOP2A", "TYMS", "STMN1", "TUBA1B")
  ),
  T_NK = list(
    CD4_naive = c("CCR7", "MAL", "LTB", "TCF7", "LEF1", "IL7R"),
    CD4_memory = c("IL7R", "LTB", "MALAT1", "IL32", "CD40LG", "ICOS"),
    Treg = c("FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2", "TNFRSF18"),
    CD8_effector = c("CD8A", "CD8B", "CCL5", "NKG7", "GZMK", "GZMH"),
    CD8_exhausted = c("PDCD1", "LAG3", "TOX", "HAVCR2", "TIGIT", "CTLA4"),
    NK_cytotoxic = c("NKG7", "GNLY", "PRF1", "GZMB", "KLRD1", "FCGR3A"),
    NK_adaptive = c("KLRC2", "KLRD1", "TRDC", "TYROBP", "FCER1G", "XCL1"),
    Proliferating_T_NK = c("MKI67", "TOP2A", "TYMS", "STMN1", "TUBA1B")
  )
)

lineage_cols <- list(
  Myeloid = c(Classical_monocyte = "#E69F00", Inflammatory_monocyte = "#D55E00",
              C1QC_macrophage = "#0072B2", SPP1_macrophage = "#CC79A7",
              Resident_macrophage = "#009E73", cDC = "#56B4E9", pDC = "#7B3294",
              Proliferating_myeloid = "#111111", Unresolved_myeloid = "#999999"),
  T_NK = c(CD4_naive = "#56B4E9", CD4_memory = "#0072B2", Treg = "#CC79A7",
           CD8_effector = "#D55E00", CD8_exhausted = "#A50F15",
           NK_cytotoxic = "#009E73", NK_adaptive = "#E69F00",
           Proliferating_T_NK = "#111111", Unresolved_T_NK = "#999999")
)

message("Reading annotated SCE objects")
sces <- lapply(files, readRDS)
gene_order <- rownames(sces[[1]])
if (!all(vapply(sces, function(x) identical(rownames(x), gene_order), logical(1)))) {
  stop("Gene order differs among SCE objects")
}

safe_z <- function(x) {
  if (length(x) < 2L || is.na(sd(x)) || sd(x) == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

make_lineage_object <- function(lineage) {
  mats <- list()
  metas <- list()
  k <- 0L
  for (sce in sces) {
    keep <- as.character(colData(sce)$broad_class) == lineage
    if (!any(keep)) next
    k <- k + 1L
    mats[[k]] <- as(assay(sce, "counts")[, keep, drop = FALSE], "dgCMatrix")
    md <- as.data.table(as.data.frame(colData(sce)[keep, , drop = FALSE]))
    md[, barcode_gate12b := colnames(sce)[keep]]
    metas[[k]] <- md
  }
  counts <- do.call(cbind, mats)
  meta <- rbindlist(metas, fill = TRUE)
  stopifnot(ncol(counts) == nrow(meta), identical(colnames(counts), meta$barcode_gate12b))
  rownames(meta) <- meta$barcode_gate12b
  obj <- CreateSeuratObject(counts = counts, meta.data = as.data.frame(meta),
                            project = paste0("gate12b_", lineage), min.cells = 1,
                            min.features = 0)
  rm(counts, mats, metas, meta); invisible(gc())
  obj
}

annotate_clusters <- function(obj, lineage, resolution) {
  message("Processing ", lineage, " cells=", ncol(obj))
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,
                       verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2500,
                              verbose = FALSE)
  variable <- setdiff(VariableFeatures(obj), grep("^(MT-|RPL|RPS|HBA|HBB)",
                                                  VariableFeatures(obj), value = TRUE))
  VariableFeatures(obj) <- variable
  obj <- ScaleData(obj, features = variable, verbose = FALSE)
  obj <- RunPCA(obj, features = variable, npcs = 35, seed.use = seed, verbose = FALSE)
  obj <- RunHarmony(obj, group.by.vars = "accession", reduction.use = "pca",
                    dims.use = 1:30, assay.use = "RNA", max_iter = 30,
                    plot_convergence = FALSE, verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:25, k.param = 30,
                       verbose = FALSE)
  obj <- FindClusters(obj, resolution = resolution, random.seed = seed,
                      algorithm = 1, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = "harmony", dims = 1:25, n.neighbors = 35,
                 min.dist = 0.22, metric = "cosine", seed.use = seed, verbose = FALSE)

  modules <- lapply(module_list[[lineage]], intersect, y = rownames(obj))
  modules <- modules[lengths(modules) >= 3L]
  dat <- GetAssayData(obj, assay = "RNA", layer = "data")
  cl <- as.character(Idents(obj))
  clusters <- sort(unique(cl))
  scores <- rbindlist(lapply(names(modules), function(module_name) {
    genes <- modules[[module_name]]
    vals <- vapply(clusters, function(cluster_id) {
      cells <- which(cl == cluster_id)
      mean(Matrix::rowMeans(dat[genes, cells, drop = FALSE]))
    }, numeric(1))
    data.table(lineage = lineage, cluster = clusters, module = module_name,
               raw_module_score = vals, module_z = safe_z(vals),
               genes_present = paste(genes, collapse = ";"))
  }))
  setorder(scores, cluster, -module_z, -raw_module_score)
  labels <- scores[, .SD[1L], by = cluster][, .(cluster, state = module,
                                                winning_module_z = module_z,
                                                winning_raw_score = raw_module_score)]
  labels[winning_module_z < 0.5,
         state := if (lineage == "Myeloid") "Unresolved_myeloid" else "Unresolved_T_NK"]
  obj$gate12b_cluster <- as.character(Idents(obj))
  obj$gate12b_state <- labels$state[match(obj$gate12b_cluster, labels$cluster)]

  coords <- as.data.table(Embeddings(obj, "umap"), keep.rownames = "barcode")
  md <- as.data.table(obj@meta.data[coords$barcode, , drop = FALSE], keep.rownames = "barcode_meta")
  coords <- cbind(coords, md)
  coords[, lineage := lineage]

  cluster_audit <- coords[, .(
    cells = .N,
    patients = uniqueN(patient_id),
    prostate_fraction = mean(cancer == "prostate"),
    distal_fraction = mean(compartment == "distal"),
    involved_fraction = mean(compartment == "involved"),
    tumor_fraction = mean(compartment == "tumor"),
    top_existing_state = names(sort(table(harmonized_state), decreasing = TRUE))[1],
    top_existing_state_fraction = max(table(harmonized_state)) / .N
  ), by = .(lineage, cluster = gate12b_cluster, state = gate12b_state)]
  cluster_audit <- merge(cluster_audit, labels, by = c("cluster", "state"), all.x = TRUE)

  marker_dt <- tryCatch({
    mk <- FindAllMarkers(obj, assay = "RNA", only.pos = TRUE, min.pct = 0.10,
                         logfc.threshold = 0.25, test.use = "wilcox",
                         max.cells.per.ident = 1200, random.seed = seed, verbose = FALSE)
    as.data.table(mk)
  }, error = function(e) data.table(error = conditionMessage(e)))

  list(obj = obj, coords = coords, scores = scores, labels = labels,
       cluster_audit = cluster_audit, markers = marker_dt)
}

myeloid <- annotate_clusters(make_lineage_object("Myeloid"), "Myeloid", 0.35)
tnk <- annotate_clusters(make_lineage_object("T_NK"), "T_NK", 0.45)
rm(sces); invisible(gc())

coords <- rbindlist(list(myeloid$coords, tnk$coords), fill = TRUE)
fwrite(coords, file.path(out, "gate12b_cell_state_coordinates.tsv.gz"), sep = "\t")
fwrite(rbindlist(list(myeloid$scores, tnk$scores), fill = TRUE),
       file.path(out, "cluster_module_scores.tsv"), sep = "\t")
fwrite(rbindlist(list(myeloid$cluster_audit, tnk$cluster_audit), fill = TRUE),
       file.path(out, "cluster_annotation_audit.tsv"), sep = "\t")
fwrite(rbindlist(list(cbind(lineage = "Myeloid", myeloid$markers),
                           cbind(lineage = "T_NK", tnk$markers)), fill = TRUE),
       file.path(out, "descriptive_cluster_markers.tsv.gz"), sep = "\t")

sample_base <- unique(coords[, .(accession, cancer, patient_id, sample_id, compartment, lineage)])
state_levels <- unique(coords[, .(lineage, state = gate12b_state)])
composition <- rbindlist(lapply(unique(sample_base$lineage), function(lin) {
  sb <- sample_base[lineage == lin]
  states <- state_levels[lineage == lin, state]
  grid <- sb[, .(state = states),
             by = .(accession, cancer, patient_id, sample_id, compartment, lineage)]
  obs <- coords[lineage == lin, .(n_state = .N),
                by = .(accession, cancer, patient_id, sample_id, compartment,
                       lineage, state = gate12b_state)]
  z <- merge(grid, obs,
             by = c("accession", "cancer", "patient_id", "sample_id", "compartment",
                    "lineage", "state"), all.x = TRUE)
  z[is.na(n_state), n_state := 0L]
  totals <- coords[lineage == lin, .(lineage_total = .N),
                   by = .(accession, cancer, patient_id, sample_id, compartment, lineage)]
  z <- merge(z, totals,
             by = c("accession", "cancer", "patient_id", "sample_id", "compartment",
                    "lineage"), all.x = TRUE)
  z[, fraction := n_state / lineage_total]
  z[, empirical_logit := qlogis((n_state + 0.5) / (lineage_total + 1))]
  z
}))
composition[, compartment := factor(compartment, levels = c("distal", "involved", "tumor"))]
composition[, stage := as.integer(compartment) - 1L]
fwrite(composition, file.path(out, "patient_state_composition.tsv"), sep = "\t")

fit_state <- function(z, accession_value, lineage_value, state_value) {
  z <- copy(z[accession == accession_value & lineage == lineage_value & state == state_value])
  complete <- z[, .(n_compartments = uniqueN(compartment), rows = .N), by = patient_id][
    n_compartments == 3L & rows == 3L, patient_id]
  z <- z[patient_id %in% complete]
  if (uniqueN(z$patient_id) < 4L) return(NULL)
  fit <- lm(empirical_logit ~ factor(patient_id) + stage, data = z)
  co <- summary(fit)$coefficients
  if (!"stage" %in% rownames(co)) return(NULL)
  beta <- unname(co["stage", "Estimate"])
  se <- unname(co["stage", "Std. Error"])
  p <- unname(co["stage", "Pr(>|t|)"])
  patients <- sort(unique(z$patient_id))
  lopo <- vapply(patients, function(left) {
    zz <- z[patient_id != left]
    f <- lm(empirical_logit ~ factor(patient_id) + stage, data = zz)
    unname(coef(f)["stage"])
  }, numeric(1))
  data.table(accession = accession_value, cancer = unique(z$cancer),
             lineage = lineage_value, state = state_value,
             n_patients = uniqueN(z$patient_id), n_samples = nrow(z),
             beta_per_step = beta, se = se, p_value = p,
             odds_ratio_per_step = exp(beta),
             odds_ratio_distal_to_tumor = exp(2 * beta),
             ci_low_per_step = exp(beta - 1.96 * se), ci_high_per_step = exp(beta + 1.96 * se),
             lopo_same_direction = mean(sign(lopo) == sign(beta)),
             min_cells_state = min(z$n_state), total_cells_state = sum(z$n_state))
}

keys <- unique(composition[, .(accession, lineage, state)])
cohort_effects <- rbindlist(lapply(seq_len(nrow(keys)), function(i) {
  fit_state(composition, keys$accession[i], keys$lineage[i], keys$state[i])
}), fill = TRUE)
if (nrow(cohort_effects) == 0L) stop("No eligible patient-level state models")
cohort_effects[, q_value := p.adjust(p_value, method = "BH"), by = .(accession, lineage)]
fwrite(cohort_effects, file.path(out, "state_cohort_effects.tsv"), sep = "\t")

meta_effects <- cohort_effects[, {
  if (.N != 2L || uniqueN(accession) != 2L || any(!is.finite(se)) || any(se <= 0)) {
    list(k_cohorts = .N, beta_meta = NA_real_, se_meta = NA_real_, p_meta = NA_real_,
         odds_ratio_per_step = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
         same_direction = FALSE, min_lopo_stability = NA_real_, Q = NA_real_, I2 = NA_real_)
  } else {
    w <- 1 / se^2
    b <- sum(w * beta_per_step) / sum(w)
    s <- sqrt(1 / sum(w))
    q <- sum(w * (beta_per_step - b)^2)
    i2 <- max(0, (q - 1) / q) * 100
    list(k_cohorts = .N, beta_meta = b, se_meta = s,
         p_meta = 2 * pnorm(-abs(b / s)), odds_ratio_per_step = exp(b),
         ci_low = exp(b - 1.96 * s), ci_high = exp(b + 1.96 * s),
         same_direction = length(unique(sign(beta_per_step))) == 1L,
         min_lopo_stability = min(lopo_same_direction), Q = q, I2 = i2)
  }
}, by = .(lineage, state)]
meta_effects[, q_meta := p.adjust(p_meta, method = "BH"), by = lineage]
meta_effects[, cross_cancer_stable := k_cohorts == 2L & same_direction &
               min_lopo_stability >= 0.75 & is.finite(beta_meta)]
meta_effects[grepl("^Unresolved", state), cross_cancer_stable := FALSE]
meta_effects[, abs_beta_meta := abs(beta_meta)]
setorder(meta_effects, lineage, q_meta, -abs_beta_meta)
fwrite(meta_effects, file.path(out, "state_meta_effects.tsv"), sep = "\t")

theme_pub <- theme_classic(base_size = 9, base_family = "Arial") +
  theme(plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 8, colour = "#444444"),
        plot.tag = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 7),
        strip.background = element_rect(fill = "#F2F2F2", colour = NA),
        strip.text = element_text(face = "bold", size = 8))

plot_umap <- function(dt, lineage_value) {
  cols <- lineage_cols[[lineage_value]]
  centres <- dt[, .(umap_1 = median(umap_1), umap_2 = median(umap_2), cells = .N),
                by = gate12b_state]
  ggplot(dt, aes(umap_1, umap_2, colour = gate12b_state)) +
    ggrastr::geom_point_rast(size = 0.30, alpha = 0.62, raster.dpi = 450) +
    ggrepel::geom_label_repel(data = centres, aes(label = gsub("_", " ", gate12b_state)),
                              seed = seed, size = 2.2, fontface = "bold",
                              colour = "black", fill = scales::alpha("white", 0.88),
                              label.size = 0.12, min.segment.length = 0,
                              show.legend = FALSE) +
    scale_colour_manual(values = cols, drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    coord_equal() + labs(title = paste0(gsub("_", "/", lineage_value), " state map"),
                         x = "UMAP 1", y = "UMAP 2", colour = "State") + theme_pub
}

p_myeloid <- plot_umap(myeloid$coords, "Myeloid")
p_tnk <- plot_umap(tnk$coords, "T_NK")

forest_dt <- meta_effects[is.finite(beta_meta) & !grepl("^Unresolved", state)]
forest_dt[, state_label := factor(gsub("_", " ", state),
                                  levels = rev(gsub("_", " ", state)))]
p_forest <- ggplot(forest_dt, aes(odds_ratio_per_step, state_label, colour = lineage)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "#777777") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.18, linewidth = 0.45) +
  geom_point(aes(shape = cross_cancer_stable), size = 2.2) +
  scale_x_log10() +
  scale_colour_manual(values = c(Myeloid = "#D55E00", T_NK = "#0072B2")) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) +
  labs(title = "Cross-cancer compartment effects",
       subtitle = "Odds ratio per distal → involved → tumour step",
       x = "Meta-analytic odds ratio per step (95% CI)", y = NULL,
       colour = "Lineage", shape = "Directionally stable") + theme_pub

top_states <- meta_effects[cross_cancer_stable == TRUE & !grepl("^Unresolved", state)][
  order(-abs(beta_meta)), head(state, 4L)]
if (length(top_states) < 2L) top_states <- meta_effects[is.finite(beta_meta)][order(-abs(beta_meta)), head(state, 4L)]
paired_dt <- composition[state %in% top_states]
paired_dt[, state_label := factor(gsub("_", " ", state), levels = gsub("_", " ", top_states))]
p_paired <- ggplot(paired_dt, aes(compartment, fraction,
                                  group = interaction(accession, patient_id), colour = cancer)) +
  geom_line(alpha = 0.30, linewidth = 0.35) +
  geom_point(size = 1.1, alpha = 0.82, position = position_jitter(width = 0.035, height = 0)) +
  facet_wrap(~ state_label, scales = "free_y", nrow = 2) +
  scale_colour_manual(values = c(prostate = "#0072B2", renal = "#D55E00"),
                      labels = c(prostate = "Prostate", renal = "Renal")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Patient-level state remodelling", x = NULL, y = "Within-lineage fraction",
       colour = "Cohort") + theme_pub +
  theme(axis.text.x = element_text(angle = 28, hjust = 1))

fig <- (p_myeloid | p_tnk) / (p_paired | p_forest) +
  plot_layout(widths = c(1.05, 1), heights = c(1.05, 1)) +
  plot_annotation(tag_levels = "A")
ggsave(file.path(out, "Figure2_cell_state_remodelling_candidate.pdf"), fig,
       width = 12.2, height = 9.2, device = cairo_pdf, bg = "white")
ggsave(file.path(out, "Figure2_cell_state_remodelling_candidate.png"), fig,
       width = 12.2, height = 9.2, dpi = 360, bg = "white")

receipt <- list(seed = seed, myeloid_resolution = 0.35, tnk_resolution = 0.45,
                myeloid_cells = nrow(myeloid$coords), tnk_cells = nrow(tnk$coords),
                patients = uniqueN(coords$patient_id), samples = uniqueN(coords$sample_id),
                integration_batch = "accession_only",
                inferential_unit = "patient",
                state_label_method = "unsupervised_cluster_plus_frozen_marker_modules")
saveRDS(receipt, file.path(out, "gate12b_receipt.rds"))

audit <- c(
  "# Gate12B cell-state remodelling audit", "",
  paste0("- Myeloid cells: ", nrow(myeloid$coords)),
  paste0("- T/NK cells: ", nrow(tnk$coords)),
  paste0("- Myeloid clusters: ", uniqueN(myeloid$coords$gate12b_cluster)),
  paste0("- T/NK clusters: ", uniqueN(tnk$coords$gate12b_cluster)),
  paste0("- Eligible cohort state models: ", nrow(cohort_effects)),
  paste0("- Cross-cancer stable states: ", sum(meta_effects$cross_cancer_stable, na.rm = TRUE)),
  "- Batch correction variable: accession only",
  "- Inferential unit: patient; cells are not statistical replicates",
  "- Myeloid old transfer labels: audit layer only because Gate5A failed",
  "- T/NK old transfer labels: retained as a validated audit layer",
  "- Anatomical compartment order is not interpreted as lineage time",
  "- Status: COMPUTE_COMPLETE_VISUAL_AND_STATISTICAL_REVIEW_PENDING"
)
writeLines(audit, file.path(out, "GATE12B_AUDIT.md"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("GATE12B_COMPLETE=TRUE\n")
