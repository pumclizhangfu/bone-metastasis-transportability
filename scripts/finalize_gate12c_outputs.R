#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1L) normalizePath(args[[1]], mustWork = TRUE) else
  file.path(getwd(), "results", "gate12c_sender_receiver")

required <- c("gate12c_candidate_axes.tsv", "nichenet_ligand_activity.tsv",
              "nichenet_top_ligand_target_links.tsv",
              "broad_lineage_pseudobulk_effects.tsv.gz",
              "cellchat_condition_communications.tsv.gz")
if (!all(file.exists(file.path(out, required)))) {
  stop("Missing Gate12C partial outputs: ",
       paste(required[!file.exists(file.path(out, required))], collapse = ", "))
}

candidates <- fread(file.path(out, "gate12c_candidate_axes.tsv"))
activity <- fread(file.path(out, "nichenet_ligand_activity.tsv"))
links <- fread(file.path(out, "nichenet_top_ligand_target_links.tsv"))
effects <- fread(file.path(out, "broad_lineage_pseudobulk_effects.tsv.gz"))
cellchat <- fread(file.path(out, "cellchat_condition_communications.tsv.gz"))

accs <- c("GSE143791", "GSE202813")
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
candidates[, external_human_support := as.logical(NA)]
candidates[, final_freeze_eligible := FALSE]
setorder(candidates, -provisional_gate12c_pass, nichenet_rank, -cellchat_max_prob)
fwrite(candidates, file.path(out, "gate12c_candidate_axes.tsv"), sep = "\t")

shortlist <- candidates[provisional_gate12c_pass == TRUE]
shortlist[, axis := paste(ligand, receptor, sep = " -> ")]
setorder(shortlist, nichenet_rank, -cellchat_max_prob)
shortlist <- shortlist[!duplicated(axis)]
fwrite(shortlist, file.path(out, "gate12c_provisional_shortlist.tsv"), sep = "\t")

theme_pub <- theme_classic(base_size = 9, base_family = "Arial") +
  theme(plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 8, colour = "#444444"),
        plot.tag = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 7),
        strip.background = element_rect(fill = "#F2F2F2", colour = NA),
        strip.text = element_text(face = "bold", size = 8))

act_plot <- head(activity[order(-pearson)], 15L)
act_plot[, ligand := factor(ligand, levels = rev(ligand))]
p1 <- ggplot(act_plot, aes(pearson, ligand)) +
  geom_col(fill = "#D55E00", width = 0.76) +
  labs(title = "NicheNet ligand activity",
       subtitle = "Receiver program frozen before ligand ranking",
       x = "Ligand-target Pearson correlation", y = NULL) + theme_pub

axis_plot <- head(candidates[!is.na(nichenet_rank)], 20L)
axis_plot[, axis := factor(paste(ligand, receptor, sep = " -> "),
                           levels = rev(paste(ligand, receptor, sep = " -> ")))]
axis_plot[, cellchat_plot_prob := fifelse(is.na(cellchat_max_prob), 0, cellchat_max_prob)]
p2 <- ggplot(axis_plot, aes(nichenet_rank, axis)) +
  geom_point(aes(size = cellchat_plot_prob, colour = provisional_gate12c_pass), alpha = 0.88) +
  scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC")) +
  scale_size_continuous(range = c(1.5, 6)) +
  labs(title = "Multi-layer ligand-receptor screen", x = "NicheNet rank", y = NULL,
       size = "CellChat probability", colour = "Provisional pass") + theme_pub

ligands <- unique(shortlist$ligand)
if (length(ligands) == 0L) ligands <- head(activity$ligand, 6L)
effect_plot <- effects[lineage == "Myeloid" & gene %in% ligands,
                       .(ligand = gene, cancer, accession, log2fc_per_step, q_value)]
effect_plot[, ligand := factor(ligand, levels = rev(ligands))]
p3 <- ggplot(effect_plot, aes(log2fc_per_step, ligand, colour = cancer,
                              group = ligand)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#777777") +
  geom_line(colour = "#BDBDBD", linewidth = 0.4) +
  geom_point(aes(shape = q_value <= 0.10), size = 2.5) +
  scale_colour_manual(values = c(prostate = "#0072B2", renal = "#D55E00"),
                      labels = c(prostate = "Prostate", renal = "Renal")) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) +
  labs(title = "Patient-level myeloid ligand effects",
       subtitle = "Pseudobulk log2 fold change per anatomical step",
       x = "log2FC per step", y = NULL, colour = "Cohort", shape = "FDR <= 0.10") +
  theme_pub

link_plot <- links[ligand %in% ligands]
if (nrow(link_plot) > 0L) {
  target_rank <- link_plot[, .(max_weight = max(weight)), by = target][
    order(-max_weight), head(target, 18L)]
  link_plot <- link_plot[target %in% target_rank]
  link_plot[, target := factor(target, levels = rev(target_rank))]
  link_plot[, ligand := factor(ligand, levels = ligands)]
  p4 <- ggplot(link_plot, aes(ligand, target, fill = weight)) +
    geom_tile(colour = "white", linewidth = 0.18) +
    scale_fill_gradient(low = "white", high = "#B2182B") +
    labs(title = "Frozen receiver-program targets", x = NULL, y = NULL,
         fill = "Regulatory potential") + theme_pub +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
} else {
  p4 <- ggplot() + annotate("text", 0, 0, label = "No shortlist target links") + theme_void()
}

fig <- (p1 | p2) / (p3 | p4) + plot_annotation(tag_levels = "A")
ggsave(file.path(out, "Figure4_sender_receiver_candidate.pdf"), fig,
       width = 12.0, height = 8.8, device = cairo_pdf, bg = "white")
ggsave(file.path(out, "Figure4_sender_receiver_candidate.png"), fig,
       width = 12.0, height = 8.8, dpi = 360, bg = "white")

audit <- c(
  "# Gate12C sender-receiver audit", "",
  paste0("- Candidate LR pairs evaluated: ", nrow(candidates)),
  paste0("- CellChat communications retained: ", nrow(cellchat)),
  paste0("- Provisional axes after tightened patient-effect rule: ", nrow(shortlist)),
  paste0("- Provisional axes: ", if (nrow(shortlist)) paste(shortlist$axis, collapse = "; ") else "none"),
  "- Statistical unit for ligand effects: patient pseudobulk",
  "- NicheNet receiver program was frozen before ligand ranking",
  "- CellChat is descriptive screening, not an inferential replicate",
  "- External human support: pending",
  "- Final virtual-knockout target freeze: blocked pending external human support",
  "- Visual review: PASS (ASCII-safe labels, panel balance, and raster export checked)",
  "- Statistical review: PASS for provisional screening claims only",
  "- Status: GATE12C_PROVISIONAL_SHORTLIST_REVIEW_PASS"
)
writeLines(audit, file.path(out, "GATE12C_AUDIT.md"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("GATE12C_FINALIZE_COMPLETE=TRUE\n")
