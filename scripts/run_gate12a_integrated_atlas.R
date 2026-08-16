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
})

args <- commandArgs(trailingOnly = TRUE)
project <- if (length(args) >= 1L) normalizePath(args[[1]], mustWork = TRUE) else
  "."
out <- if (length(args) >= 2L) args[[2]] else file.path(project, "results", "gate12a_integrated_atlas")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

seed <- 20260808L
set.seed(seed)
sce_dir <- file.path(project, "data", "gate3b_work", "annotated_sce")
files <- sort(list.files(sce_dir, pattern = "\\.rds$", full.names = TRUE))
if (length(files) != 42L) stop("Expected 42 annotated discovery SCE files; found ", length(files))

message("Reading 42 annotated SCE objects")
sces <- lapply(files, readRDS)
gene_order <- rownames(sces[[1]])
if (!all(vapply(sces, function(x) identical(rownames(x), gene_order), logical(1)))) {
  stop("Gene order is not identical across SCE objects")
}

mats <- lapply(sces, function(x) as(assay(x, "counts"), "dgCMatrix"))
counts_mat <- do.call(cbind, mats)
meta <- rbindlist(lapply(sces, function(x) as.data.table(as.data.frame(colData(x)))), fill = TRUE)
if (ncol(counts_mat) != nrow(meta)) stop("Count/metadata cell mismatch")
if (anyDuplicated(colnames(counts_mat))) stop("Cell barcodes are not globally unique")
rownames(meta) <- colnames(counts_mat)
rm(mats, sces); invisible(gc())

meta[, compartment := factor(compartment, levels = c("distal", "involved", "tumor"))]
meta[, cancer := factor(cancer, levels = c("prostate", "renal"))]
meta[, broad_class := factor(broad_class,
  levels = c("T_NK", "Myeloid", "B", "Progenitor", "Stromal", "Endothelial",
             "Osteoclast", "Osteoblast", "Malignant", "Erythroid", "Unassigned"))]

obj <- CreateSeuratObject(counts = counts_mat, meta.data = as.data.frame(meta),
                          project = "paired_bone_metastasis", min.cells = 1, min.features = 0)
rm(counts_mat, meta); invisible(gc())

message("Normalizing and selecting variable features")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000,
                     verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 40, verbose = FALSE,
              seed.use = seed)

message("Running Harmony by source accession")
obj <- RunHarmony(obj, group.by.vars = "accession", reduction.use = "pca",
                  dims.use = 1:35, assay.use = "RNA", max_iter = 30,
                  plot_convergence = FALSE, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, n.neighbors = 40,
               min.dist = 0.25, metric = "cosine", seed.use = seed, verbose = FALSE)

coords <- as.data.table(Embeddings(obj, "umap"), keep.rownames = "barcode")
coords <- cbind(coords, as.data.table(obj@meta.data[coords$barcode, , drop = FALSE]))
fwrite(coords, file.path(out, "integrated_umap_coordinates.tsv.gz"), sep = "\t")

sample_audit <- coords[, .(
  cells = .N,
  median_nCount = median(nCount),
  median_nFeature = median(nFeature),
  median_percent_mt = median(percent_mt)
), by = .(accession, cancer, patient_id, sample_id, compartment)]
fwrite(sample_audit, file.path(out, "sample_atlas_audit.tsv"), sep = "\t")

class_counts <- coords[, .N, by = .(accession, cancer, patient_id, sample_id, compartment, broad_class)]
class_counts[, sample_total := sum(N), by = sample_id]
class_counts[, fraction := N / sample_total]
fwrite(class_counts, file.path(out, "broad_class_composition.tsv"), sep = "\t")

markers <- c("CD3D", "CD8A", "NKG7", "MS4A1", "LST1", "FCGR3A", "CD34", "COL1A1",
             "PECAM1", "ACP5", "BGLAP", "EPCAM", "KRT18", "HBB")
markers <- intersect(markers, rownames(obj))
Idents(obj) <- "broad_class"
classes <- levels(droplevels(obj$broad_class))
rna_data <- GetAssayData(obj, assay = "RNA", layer = "data")
rna_counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
avg <- sapply(classes, function(cl) {
  cells <- rownames(obj@meta.data)[obj$broad_class == cl]
  Matrix::rowMeans(rna_data[markers, cells, drop = FALSE])
})
pct <- sapply(classes, function(cl) {
  cells <- rownames(obj@meta.data)[obj$broad_class == cl]
  Matrix::rowMeans(rna_counts[markers, cells, drop = FALSE] > 0)
})
rownames(avg) <- rownames(pct) <- markers
colnames(avg) <- colnames(pct) <- classes
marker_summary <- rbindlist(lapply(seq_len(nrow(avg)), function(i) data.table(
  gene = rownames(avg)[i], broad_class = colnames(avg),
  average_log_normalized = as.numeric(avg[i, ]), detected_fraction = as.numeric(pct[i, ])
)))
fwrite(marker_summary, file.path(out, "marker_expression_summary.tsv"), sep = "\t")

class_cols <- c(
  T_NK = "#0077BB", Myeloid = "#EE7733", B = "#33BBEE", Progenitor = "#CC3311",
  Stromal = "#009988", Endothelial = "#EE3377", Osteoclast = "#AA4499",
  Osteoblast = "#BBBBBB", Malignant = "#000000", Erythroid = "#EE99AA",
  Unassigned = "#999999"
)
theme_atlas <- theme_classic(base_size = 8.5, base_family = "Arial") +
  theme(plot.title = element_text(face = "bold", size = 10),
        legend.title = element_text(size = 8), legend.text = element_text(size = 7),
        plot.tag = element_text(face = "bold", size = 11))

p1 <- DimPlot(obj, reduction = "umap", group.by = "broad_class", raster = TRUE,
              raster.dpi = c(300, 300), pt.size = 0.16, shuffle = TRUE,
              label = TRUE, repel = TRUE) +
  scale_color_manual(values = class_cols, drop = FALSE) +
  labs(title = "Cell atlas", color = "Cell class") + theme_atlas
p2 <- DimPlot(obj, reduction = "umap", group.by = "cancer", raster = TRUE,
              raster.dpi = c(300, 300), pt.size = 0.16, shuffle = TRUE) +
  scale_color_manual(values = c(prostate = "#0077BB", renal = "#EE7733")) +
  labs(title = "Cancer cohort", color = "Cancer") + theme_atlas
p3 <- DimPlot(obj, reduction = "umap", group.by = "compartment", raster = TRUE,
              raster.dpi = c(300, 300), pt.size = 0.16, shuffle = TRUE) +
  scale_color_manual(values = c(distal = "#33BBEE", involved = "#EEAA00", tumor = "#CC3311")) +
  labs(title = "Anatomical compartment", color = "Compartment") + theme_atlas
p4 <- DotPlot(obj, features = markers, group.by = "broad_class", assay = "RNA") +
  scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  scale_size(range = c(0, 5)) +
  coord_flip() + labs(title = "Canonical marker audit", x = NULL, y = NULL,
                      color = "Scaled mean", size = "% detected") + theme_atlas +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fig <- (p1 | p2) / (p3 | p4) + plot_annotation(tag_levels = "A")
ggsave(file.path(out, "Figure12A_integrated_single_cell_atlas.pdf"), fig,
       width = 10.5, height = 8.2, device = cairo_pdf, bg = "white")
ggsave(file.path(out, "Figure12A_integrated_single_cell_atlas.png"), fig,
       width = 10.5, height = 8.2, dpi = 300, type = "cairo", bg = "white")

reduction_out <- list(
  variable_features = VariableFeatures(obj),
  pca_stdev = Stdev(obj, reduction = "pca"),
  harmony_stdev = Stdev(obj, reduction = "harmony"),
  seed = seed,
  cells = ncol(obj), genes = nrow(obj), samples = uniqueN(obj$sample_id),
  patients = uniqueN(obj$patient_id)
)
saveRDS(reduction_out, file.path(out, "atlas_reduction_receipt.rds"))

audit <- c(
  "# Gate12A integrated single-cell atlas audit", "",
  paste0("- Cells: ", ncol(obj)), paste0("- Genes: ", nrow(obj)),
  paste0("- Samples: ", uniqueN(obj$sample_id)), paste0("- Patients: ", uniqueN(obj$patient_id)),
  "- Modality: scRNA-seq", "- Integration batch: accession only",
  "- Biological inference unit: patient; UMAP is descriptive only",
  "- Spatial transcriptomics claim: NO", "- Status: RENDER_COMPLETE_VISUAL_REVIEW_PENDING"
)
writeLines(audit, file.path(out, "GATE12A_AUDIT.md"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("GATE12A_RENDER_COMPLETE=TRUE\n")
