#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(data.table)
  library(edgeR)
  library(nichenetr)
  library(CellChat)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
project <- if (length(args) >= 1L) normalizePath(args[[1]], mustWork = TRUE) else
  "."
gate12b <- if (length(args) >= 2L) normalizePath(args[[2]], mustWork = TRUE) else
  file.path(project, "results", "gate12b_cell_states")
prior_dir <- if (length(args) >= 3L) normalizePath(args[[3]], mustWork = TRUE) else
  "data/external/nichenet_prior"
out <- if (length(args) >= 4L) args[[4]] else file.path(project, "results", "gate12c_sender_receiver")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

seed <- 20260808L
set.seed(seed)
lr_file <- file.path(prior_dir, "lr_network_human_21122021.rds")
lt_file <- file.path(prior_dir, "ligand_target_matrix_nsga2r_final.rds")
stopifnot(file.exists(lr_file), file.exists(lt_file))

message("Loading Gate12B state map and NicheNet priors")
state_map <- fread(file.path(gate12b, "gate12b_cell_state_coordinates.tsv.gz"),
                   select = c("barcode", "accession", "cancer", "patient_id", "sample_id",
                              "compartment", "lineage", "gate12b_state"))
if (anyDuplicated(state_map$barcode)) stop("Duplicated Gate12B barcodes")
lr_network <- unique(as.data.table(readRDS(lr_file))[, .(ligand = from, receptor = to)])
ligand_target_matrix <- readRDS(lt_file)

sce_dir <- file.path(project, "data", "gate3b_work", "annotated_sce")
files <- sort(list.files(sce_dir, pattern = "\\.rds$", full.names = TRUE))
if (length(files) != 42L) stop("Expected 42 SCE files")
sces <- lapply(files, readRDS)
gene_order <- rownames(sces[[1]])
if (!all(vapply(sces, function(x) identical(rownames(x), gene_order), logical(1)))) {
  stop("Gene order mismatch")
}
mats <- list(); metas <- list(); k <- 0L
for (sce in sces) {
  keep <- as.character(colData(sce)$broad_class) %in% c("Myeloid", "T_NK")
  if (!any(keep)) next
  k <- k + 1L
  mats[[k]] <- as(assay(sce, "counts")[, keep, drop = FALSE], "dgCMatrix")
  metas[[k]] <- data.table(barcode = colnames(sce)[keep])
}
counts <- do.call(cbind, mats)
cell_meta <- rbindlist(metas)
stopifnot(identical(colnames(counts), cell_meta$barcode))
cell_meta <- merge(cell_meta, state_map, by = "barcode", all.x = TRUE, sort = FALSE)
cell_meta <- cell_meta[match(colnames(counts), barcode)]
if (anyNA(cell_meta$gate12b_state)) stop("State map did not cover all cells")
unresolved_cells_excluded <- sum(grepl("^Unresolved", cell_meta$gate12b_state))
resolved <- !grepl("^Unresolved", cell_meta$gate12b_state)
counts <- counts[, resolved, drop = FALSE]
cell_meta <- cell_meta[resolved]
rm(sces, mats, metas, state_map); invisible(gc())

# Keep only simple single-gene ligand/receptor entries represented in the matrix.
lr_network <- lr_network[ligand %in% rownames(counts) & receptor %in% rownames(counts)]
lr_genes <- intersect(unique(c(lr_network$ligand, lr_network$receptor)), rownames(counts))
message("Simple LR pairs represented: ", nrow(lr_network), "; genes: ", length(lr_genes))

detect_by_group <- function(group_dt, genes) {
  group_dt[, group_id := .GRP,
           by = .(accession, cancer, compartment, lineage, gate12b_state)]
  design <- sparseMatrix(i = seq_len(nrow(group_dt)), j = group_dt$group_id, x = 1,
                         dims = c(nrow(group_dt), max(group_dt$group_id)))
  det <- (counts[genes, , drop = FALSE] > 0) %*% design
  sizes <- as.numeric(colSums(design))
  pct <- sweep(as.matrix(det), 2, sizes, "/")
  keys <- unique(group_dt[, .(group_id, accession, cancer, compartment, lineage,
                              state = gate12b_state)])
  setorder(keys, group_id)
  rbindlist(lapply(seq_len(ncol(pct)), function(j) data.table(
    keys[j], gene = rownames(pct), detected_fraction = pct[, j]
  )))
}

expression_support <- detect_by_group(copy(cell_meta), lr_genes)
fwrite(expression_support, file.path(out, "lr_expression_support_by_state.tsv.gz"), sep = "\t")

sender <- expression_support[lineage == "Myeloid" & compartment %in% c("involved", "tumor")]
sender <- sender[, .SD[which.max(detected_fraction)], by = .(accession, gene)][,
                 .(accession, ligand = gene, sender_state = state,
                   sender_detected_fraction = detected_fraction)]
receiver <- expression_support[lineage == "T_NK" & compartment %in% c("involved", "tumor")]
receiver <- receiver[, .SD[which.max(detected_fraction)], by = .(accession, gene)][,
                     .(accession, receptor = gene, receiver_state = state,
                       receiver_detected_fraction = detected_fraction)]

lr_support <- rbindlist(lapply(unique(cell_meta$accession), function(acc) {
  z <- merge(lr_network, sender[accession == acc], by = "ligand")
  z <- merge(z, receiver[accession == acc], by = "receptor")
  z[, accession := acc]
  z
}))
lr_support[, expressed_pair := sender_detected_fraction >= 0.10 &
             receiver_detected_fraction >= 0.10]
pair_cross <- dcast(lr_support, ligand + receptor ~ accession,
                    value.var = "expressed_pair", fill = FALSE)
accs <- sort(unique(cell_meta$accession))
if (!all(accs %in% names(pair_cross))) stop("Accession support columns missing")
pair_cross[, expression_both_cancers := Reduce(`&`, .SD), .SDcols = accs]
fwrite(lr_support, file.path(out, "lr_pair_expression_support.tsv"), sep = "\t")

# Broad-lineage pseudobulk counts and patient-fixed compartment effects.
cell_meta[, pseudobulk_id := paste(sample_id, lineage, sep = "::")]
pb_levels <- unique(cell_meta$pseudobulk_id)
pb_design <- sparseMatrix(i = seq_len(nrow(cell_meta)),
                          j = match(cell_meta$pseudobulk_id, pb_levels), x = 1,
                          dims = c(nrow(cell_meta), length(pb_levels)))
pb_counts <- counts %*% pb_design
colnames(pb_counts) <- pb_levels
pb_meta <- unique(cell_meta[, .(pseudobulk_id, accession, cancer, patient_id,
                                sample_id, compartment, lineage)])
pb_meta <- pb_meta[match(colnames(pb_counts), pseudobulk_id)]
pb_meta[, compartment := factor(compartment, levels = c("distal", "involved", "tumor"))]
pb_meta[, stage := as.integer(compartment) - 1L]

fit_pb <- function(accession_value, lineage_value) {
  md <- copy(pb_meta[accession == accession_value & lineage == lineage_value])
  complete <- md[, .(n_comp = uniqueN(compartment), rows = .N), by = patient_id][
    n_comp == 3L & rows == 3L, patient_id]
  md <- md[patient_id %in% complete]
  if (uniqueN(md$patient_id) < 4L) return(NULL)
  mat <- pb_counts[, md$pseudobulk_id, drop = FALSE]
  design <- model.matrix(~ factor(patient_id) + stage, data = md)
  y <- DGEList(mat)
  keep <- filterByExpr(y, design = design, min.count = 5, min.total.count = 15)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  stage_col <- which(colnames(design) == "stage")
  tab <- as.data.table(topTags(glmQLFTest(fit, coef = stage_col), n = Inf, sort.by = "none")$table,
                       keep.rownames = "gene")
  setnames(tab, c("logFC", "logCPM", "F", "PValue", "FDR"),
           c("log2fc_per_step", "logCPM", "F_stat", "p_value", "q_value"))

  patients <- sort(unique(md$patient_id))
  lopo <- rbindlist(lapply(patients, function(left) {
    md2 <- md[patient_id != left]
    mat2 <- pb_counts[rownames(y), md2$pseudobulk_id, drop = FALSE]
    d2 <- model.matrix(~ factor(patient_id) + stage, data = md2)
    yy <- DGEList(mat2)
    yy <- calcNormFactors(yy)
    yy <- estimateDisp(yy, d2, robust = TRUE)
    ff <- glmQLFit(yy, d2, robust = TRUE)
    tt <- as.data.table(topTags(glmQLFTest(ff, coef = which(colnames(d2) == "stage")),
                               n = Inf, sort.by = "none")$table, keep.rownames = "gene")
    tt[, .(gene, left_out = left, lopo_log2fc = logFC)]
  }))
  stab <- lopo[, .(lopo_same_direction = mean(sign(lopo_log2fc) ==
                                                sign(tab$log2fc_per_step[match(gene, tab$gene)]))),
               by = gene]
  tab <- merge(tab, stab, by = "gene", all.x = TRUE)
  tab[, `:=`(accession = accession_value, cancer = unique(md$cancer),
             lineage = lineage_value, n_patients = uniqueN(md$patient_id))]
  tab
}

effects <- rbindlist(lapply(accs, function(acc) {
  rbindlist(lapply(c("Myeloid", "T_NK"), function(lin) fit_pb(acc, lin)), fill = TRUE)
}), fill = TRUE)
fwrite(effects, file.path(out, "broad_lineage_pseudobulk_effects.tsv.gz"), sep = "\t")

# NicheNet target program is frozen from Gate6B before ligand ranking.
target_file <- file.path(project, "results", "gate6b_confounded_robustness",
                         "common_persistent_features.tsv")
target_dt <- fread(target_file)
geneset <- unique(target_dt[feature_type == "gene", feature])
geneset <- intersect(geneset, rownames(ligand_target_matrix))

tnk_detect <- rbindlist(lapply(accs, function(acc) {
  idx <- which(cell_meta$accession == acc & cell_meta$lineage == "T_NK" &
                 cell_meta$compartment %in% c("involved", "tumor"))
  pct <- Matrix::rowMeans(counts[, idx, drop = FALSE] > 0)
  data.table(accession = acc, gene = rownames(counts), detected_fraction = pct)
}))
tnk_background <- tnk_detect[, .(both = all(detected_fraction >= 0.05), cohorts = .N),
                             by = gene][cohorts == length(accs) & both == TRUE, gene]
background <- intersect(tnk_background, rownames(ligand_target_matrix))
geneset <- intersect(geneset, background)

potential_pairs <- pair_cross[expression_both_cancers == TRUE]
potential_ligands <- intersect(unique(potential_pairs$ligand), colnames(ligand_target_matrix))
if (length(geneset) < 20L || length(background) < 500L || length(potential_ligands) < 5L) {
  stop("Insufficient NicheNet inputs: geneset=", length(geneset),
       " background=", length(background), " ligands=", length(potential_ligands))
}
ligand_activity <- as.data.table(predict_ligand_activities(
  geneset = geneset, background_expressed_genes = background,
  ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands
))
setnames(ligand_activity, "test_ligand", "ligand")
setorder(ligand_activity, -pearson)
ligand_activity[, nichenet_rank := seq_len(.N)]
fwrite(ligand_activity, file.path(out, "nichenet_ligand_activity.tsv"), sep = "\t")

# CellChat is a descriptive screen only. Each condition is state-balanced by capped downsampling.
data(CellChatDB.human)
secreted_db <- subsetDB(CellChatDB.human, search = "Secreted Signaling", key = "annotation")
run_cellchat <- function(accession_value, compartment_value) {
  idx <- which(cell_meta$accession == accession_value &
                 cell_meta$compartment == compartment_value)
  if (length(idx) < 100L) return(NULL)
  md <- copy(cell_meta[idx])
  md[, local_idx := idx]
  chosen <- md[, {
    n_take <- min(.N, 450L)
    .SD[sample.int(.N, n_take)]
  }, by = gate12b_state]
  chosen <- chosen[, .SD[.N >= 20L], by = gate12b_state]
  if (uniqueN(chosen$gate12b_state) < 2L) return(NULL)
  sub_counts <- counts[, chosen$local_idx, drop = FALSE]
  lib <- Matrix::colSums(sub_counts)
  norm <- log1p(t(t(sub_counts) / pmax(lib, 1)) * 10000)
  cc_meta <- data.frame(state = chosen$gate12b_state,
                        row.names = colnames(sub_counts), stringsAsFactors = FALSE)
  cc <- createCellChat(object = norm, meta = cc_meta, group.by = "state")
  cc@DB <- secreted_db
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc, thresh.pc = 0.10, thresh.fc = 0.10,
                                   thresh.p = 0.05)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "truncatedMean", trim = 0.10,
                          population.size = FALSE, raw.use = TRUE,
                          distance.use = FALSE, nboot = 50, seed.use = seed)
  cc <- filterCommunication(cc, min.cells = 20)
  res <- as.data.table(subsetCommunication(cc))
  if (nrow(res) == 0L) return(NULL)
  res[, `:=`(accession = accession_value, compartment = compartment_value,
             cells_used = ncol(norm))]
  res
}

cellchat_results <- rbindlist(lapply(accs, function(acc) {
  rbindlist(lapply(c("distal", "involved", "tumor"), function(comp) {
    tryCatch(run_cellchat(acc, comp), error = function(e) {
      fwrite(data.table(accession = acc, compartment = comp,
                        error = conditionMessage(e)),
             file.path(out, paste0("cellchat_error_", acc, "_", comp, ".tsv")), sep = "\t")
      NULL
    })
  }), fill = TRUE)
}), fill = TRUE)
if (nrow(cellchat_results) > 0L) {
  fwrite(cellchat_results, file.path(out, "cellchat_condition_communications.tsv.gz"), sep = "\t")
}

myeloid_states <- unique(cell_meta[lineage == "Myeloid", gate12b_state])
tnk_states <- unique(cell_meta[lineage == "T_NK", gate12b_state])
if (nrow(cellchat_results) > 0L) {
  cc_axis <- cellchat_results[source %in% myeloid_states & target %in% tnk_states &
                                compartment %in% c("involved", "tumor")]
  cc_axis <- cc_axis[, .(cellchat_max_prob = max(prob, na.rm = TRUE),
                         cellchat_min_p = min(pval, na.rm = TRUE),
                         best_source = source[which.max(prob)],
                         best_target = target[which.max(prob)]),
                     by = .(accession, ligand, receptor)]
  cc_cross <- cc_axis[, .(cellchat_cohorts = uniqueN(accession),
                          cellchat_both_cancers = uniqueN(accession) == length(accs),
                          cellchat_max_prob = max(cellchat_max_prob),
                          cellchat_min_p = min(cellchat_min_p),
                          best_source = best_source[which.max(cellchat_max_prob)],
                          best_target = best_target[which.max(cellchat_max_prob)]),
                      by = .(ligand, receptor)]
} else {
  cc_cross <- data.table(ligand = character(), receptor = character(),
                         cellchat_cohorts = integer(), cellchat_both_cancers = logical(),
                         cellchat_max_prob = numeric(), cellchat_min_p = numeric(),
                         best_source = character(), best_target = character())
}

# Combine pre-specified evidence layers. External human validation remains a separate gate.
candidates <- merge(potential_pairs[, .(ligand, receptor, expression_both_cancers)],
                    ligand_activity, by = "ligand", all.x = TRUE)
candidates <- merge(candidates, cc_cross, by = c("ligand", "receptor"), all.x = TRUE)
lig_eff <- effects[lineage == "Myeloid" & gene %in% unique(candidates$ligand),
                   .(ligand = gene, accession, ligand_log2fc = log2fc_per_step,
                     ligand_q = q_value, ligand_lopo = lopo_same_direction)]
lig_wide <- dcast(lig_eff, ligand ~ accession,
                  value.var = c("ligand_log2fc", "ligand_q", "ligand_lopo"))
candidates <- merge(candidates, lig_wide, by = "ligand", all.x = TRUE)
fc_cols <- paste0("ligand_log2fc_", accs)
q_cols <- paste0("ligand_q_", accs)
lopo_cols <- paste0("ligand_lopo_", accs)
fc_mat <- as.matrix(candidates[, ..fc_cols])
q_mat <- as.matrix(candidates[, ..q_cols])
lopo_mat <- as.matrix(candidates[, ..lopo_cols])
candidates[, ligand_direction_both := rowSums(is.na(fc_mat)) == 0L &
             apply(fc_mat, 1, function(x) length(unique(sign(x))) == 1L)]
candidates[, ligand_lopo_both := rowSums(is.na(lopo_mat)) == 0L &
             apply(lopo_mat, 1, function(x) all(x >= 0.75))]
candidates[, ligand_effect_supported := rowSums(is.na(fc_mat)) == 0L &
             apply(fc_mat, 1, function(x) all(x > 0)) &
             rowSums(is.na(q_mat)) == 0L & apply(q_mat, 1, function(x) all(x <= 0.10))]
candidates[is.na(cellchat_both_cancers), cellchat_both_cancers := FALSE]
candidates[, provisional_gate12c_pass := expression_both_cancers &
             !is.na(nichenet_rank) & nichenet_rank <= 50L &
             ligand_direction_both & ligand_lopo_both & ligand_effect_supported &
             cellchat_both_cancers]
candidates[, external_human_support := NA]
candidates[, final_freeze_eligible := FALSE]
setorder(candidates, -provisional_gate12c_pass, nichenet_rank, -cellchat_max_prob)
fwrite(candidates, file.path(out, "gate12c_candidate_axes.tsv"), sep = "\t")

top_ligands <- unique(head(candidates[!is.na(nichenet_rank), ligand], 12L))
target_links <- rbindlist(lapply(top_ligands, function(lig) {
  as.data.table(get_weighted_ligand_target_links(ligand = lig, geneset = geneset,
                                                 ligand_target_matrix = ligand_target_matrix,
                                                 n = 100))
}), fill = TRUE)
fwrite(target_links, file.path(out, "nichenet_top_ligand_target_links.tsv"), sep = "\t")

theme_pub <- theme_classic(base_size = 9, base_family = "Arial") +
  theme(plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 8, colour = "#444444"),
        plot.tag = element_text(face = "bold", size = 12),
        strip.background = element_rect(fill = "#F2F2F2", colour = NA),
        strip.text = element_text(face = "bold", size = 8))

act_plot <- head(ligand_activity, 15L)
act_plot[, ligand := factor(ligand, levels = rev(ligand))]
p1 <- ggplot(act_plot, aes(pearson, ligand)) +
  geom_col(fill = "#D55E00", width = 0.75) +
  labs(title = "NicheNet ligand activity", subtitle = "Frozen Gate6B T/NK target program",
       x = "Ligand–target Pearson correlation", y = NULL) + theme_pub

axis_plot <- head(candidates[!is.na(nichenet_rank)], 18L)
axis_plot[, axis := factor(paste(ligand, receptor, sep = " → "),
                           levels = rev(paste(ligand, receptor, sep = " → ")))]
p2 <- ggplot(axis_plot, aes(nichenet_rank, axis)) +
  geom_point(aes(size = cellchat_max_prob, colour = provisional_gate12c_pass), alpha = 0.85) +
  scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC")) +
  labs(title = "Multi-layer ligand–receptor screen", x = "NicheNet rank", y = NULL,
       size = "CellChat probability", colour = "Provisional pass") + theme_pub

forest <- lig_eff[ligand %in% top_ligands]
forest[, label := factor(ligand, levels = rev(top_ligands))]
p3 <- ggplot(forest, aes(ligand_log2fc, label,
                         colour = cancer)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#777777") +
  geom_point(size = 2) +
  scale_colour_manual(values = c(prostate = "#0072B2", renal = "#D55E00")) +
  labs(title = "Patient-level myeloid ligand effects",
       subtitle = "Pseudobulk log2 fold change per anatomical step",
       x = "log2FC per step", y = NULL, colour = "Cohort") + theme_pub

link_plot <- target_links[target %in% head(names(sort(table(target_links$target), decreasing = TRUE)), 18L)]
if (nrow(link_plot) > 0L) {
  p4 <- ggplot(link_plot, aes(ligand, target, fill = weight)) +
    geom_tile(colour = "white", linewidth = 0.15) +
    scale_fill_gradient(low = "white", high = "#B2182B") +
    labs(title = "Predicted ligand–target links", x = NULL, y = NULL,
         fill = "Regulatory potential") + theme_pub +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
} else {
  p4 <- ggplot() + annotate("text", 0, 0, label = "No target links passed") + theme_void()
}

fig <- (p1 | p2) / (p3 | p4) + plot_annotation(tag_levels = "A")
ggsave(file.path(out, "Figure4_sender_receiver_candidate.pdf"), fig,
       width = 12.0, height = 8.8, device = cairo_pdf, bg = "white")
ggsave(file.path(out, "Figure4_sender_receiver_candidate.png"), fig,
       width = 12.0, height = 8.8, dpi = 360, bg = "white")

audit <- c(
  "# Gate12C sender–receiver audit", "",
  paste0("- Simple human LR pairs represented: ", nrow(lr_network)),
  paste0("- Cross-cancer expressed LR pairs: ", nrow(potential_pairs)),
  paste0("- Frozen receiver target genes available: ", length(geneset)),
  paste0("- Potential NicheNet ligands: ", length(potential_ligands)),
  paste0("- Unresolved Gate12B cells excluded: ", unresolved_cells_excluded),
  paste0("- CellChat communications retained: ", nrow(cellchat_results)),
  paste0("- Provisional Gate12C axes: ", sum(candidates$provisional_gate12c_pass, na.rm = TRUE)),
  "- Statistical unit for expression effects: patient pseudobulk",
  "- CellChat is used as descriptive interaction screening, not as inferential replication",
  "- NicheNet target program was frozen from Gate6B before ligand ranking",
  "- Final target freeze: BLOCKED_PENDING_EXTERNAL_HUMAN_EXPRESSION_SUPPORT",
  "- Status: COMPUTE_COMPLETE_REVIEW_PENDING"
)
writeLines(audit, file.path(out, "GATE12C_AUDIT.md"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("GATE12C_COMPLETE=TRUE\n")
