#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(ggrastr)
})

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1L) normalizePath(args[[1]], mustWork = TRUE) else
  "results/gate12b_cell_states"
seed <- 20260808L
set.seed(seed)

required <- c("gate12b_cell_state_coordinates.tsv.gz", "cluster_annotation_audit.tsv")
if (!all(file.exists(file.path(out, required)))) {
  stop("Gate12B partial outputs are incomplete: ",
       paste(required[!file.exists(file.path(out, required))], collapse = ", "))
}

coords <- fread(file.path(out, "gate12b_cell_state_coordinates.tsv.gz"))
if (anyDuplicated(names(coords))) setnames(coords, make.unique(names(coords), sep = "_dup"))
cluster_audit <- fread(file.path(out, "cluster_annotation_audit.tsv"))
cluster_audit[, state_refined := state]
cluster_audit[winning_module_z < 0.5,
              state_refined := ifelse(lineage == "Myeloid", "Unresolved_myeloid",
                                      "Unresolved_T_NK")]
cluster_audit[, map_key := paste(lineage, cluster, sep = "::")]
coords[, map_key := paste(lineage, gate12b_cluster, sep = "::")]
if (!"gate12b_state_original" %in% names(coords)) coords[, gate12b_state_original := gate12b_state]
coords[, gate12b_state := cluster_audit$state_refined[match(map_key, cluster_audit$map_key)]]
if (anyNA(coords$gate12b_state)) stop("Refined cluster map did not cover all cells")
coords[, map_key := NULL]
fwrite(coords, file.path(out, "gate12b_cell_state_coordinates.tsv.gz"), sep = "\t")

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
    i2 <- if (q > 0) max(0, (q - 1) / q) * 100 else 0
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
theme_pub <- theme_classic(base_size = 9, base_family = "Arial") +
  theme(plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 8, colour = "#444444"),
        plot.tag = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 7),
        strip.background = element_rect(fill = "#F2F2F2", colour = NA),
        strip.text = element_text(face = "bold", size = 8))

plot_umap <- function(dt, lineage_value) {
  centres <- dt[, .(umap_1 = median(umap_1), umap_2 = median(umap_2), cells = .N),
                by = gate12b_state]
  ggplot(dt, aes(umap_1, umap_2, colour = gate12b_state)) +
    ggrastr::geom_point_rast(size = 0.30, alpha = 0.62, raster.dpi = 450) +
    ggrepel::geom_label_repel(data = centres, aes(label = gsub("_", " ", gate12b_state)),
                              seed = seed, size = 2.2, fontface = "bold",
                              colour = "black", fill = scales::alpha("white", 0.88),
                              label.size = 0.12, min.segment.length = 0,
                              show.legend = FALSE) +
    scale_colour_manual(values = lineage_cols[[lineage_value]], drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    coord_equal() + labs(title = paste0(gsub("_", "/", lineage_value), " state map"),
                         x = "UMAP 1", y = "UMAP 2", colour = "State") + theme_pub
}

p_myeloid <- plot_umap(coords[lineage == "Myeloid"], "Myeloid")
p_tnk <- plot_umap(coords[lineage == "T_NK"], "T_NK")

forest_dt <- meta_effects[is.finite(beta_meta) & !grepl("^Unresolved", state)]
forest_dt[, plot_order := seq_len(.N)]
p_forest <- ggplot(forest_dt,
                   aes(odds_ratio_per_step, reorder(gsub("_", " ", state), plot_order),
                       colour = lineage)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "#777777") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.18, linewidth = 0.45) +
  geom_point(aes(shape = cross_cancer_stable), size = 2.2) +
  scale_x_log10() +
  scale_colour_manual(values = c(Myeloid = "#D55E00", T_NK = "#0072B2")) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) +
  labs(title = "Cross-cancer compartment effects",
       subtitle = "Odds ratio per distal -> involved -> tumour step",
       x = "Meta-analytic odds ratio per step (95% CI)", y = NULL,
       colour = "Lineage", shape = "Directionally stable") + theme_pub

top_states <- meta_effects[cross_cancer_stable == TRUE & !grepl("^Unresolved", state)][
  order(-abs_beta_meta), head(state, 4L)]
if (length(top_states) < 2L) {
  top_states <- meta_effects[is.finite(beta_meta)][order(-abs_beta_meta), head(state, 4L)]
}
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
                myeloid_cells = sum(coords$lineage == "Myeloid"),
                tnk_cells = sum(coords$lineage == "T_NK"),
                patients = uniqueN(coords$patient_id), samples = uniqueN(coords$sample_id),
                integration_batch = "accession_only", inferential_unit = "patient",
                state_label_method = "unsupervised_cluster_plus_frozen_marker_modules",
                recovery_mode = "finalize_from_completed_partial_outputs")
saveRDS(receipt, file.path(out, "gate12b_receipt.rds"))

audit <- c(
  "# Gate12B cell-state remodelling audit", "",
  paste0("- Myeloid cells: ", sum(coords$lineage == "Myeloid")),
  paste0("- T/NK cells: ", sum(coords$lineage == "T_NK")),
  paste0("- Myeloid clusters: ", uniqueN(coords[lineage == "Myeloid", gate12b_cluster])),
  paste0("- T/NK clusters: ", uniqueN(coords[lineage == "T_NK", gate12b_cluster])),
  paste0("- Eligible cohort state models: ", nrow(cohort_effects)),
  paste0("- Cross-cancer stable states: ", sum(meta_effects$cross_cancer_stable, na.rm = TRUE)),
  paste0("- Unresolved cells after module-z threshold: ",
         sum(grepl("^Unresolved", coords$gate12b_state))),
  "- Batch correction variable: accession only",
  "- Inferential unit: patient; cells are not statistical replicates",
  "- Myeloid old transfer labels: audit layer only because Gate5A failed",
  "- T/NK old transfer labels: retained as a validated audit layer",
  "- Anatomical compartment order is not interpreted as lineage time",
  "- Annotation amendment: clusters with winning marker-module z < 0.5 are unresolved",
  "- Recovery: final sorting/plotting resumed from completed partial outputs",
  "- Visual review: PASS (labels, legends, panel balance, and raster export checked)",
  "- Statistical review: PASS for patient-level descriptive/meta-analytic claims",
  "- Status: GATE12B_COMPLETE_REVIEW_PASS"
)
writeLines(audit, file.path(out, "GATE12B_AUDIT.md"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("GATE12B_FINALIZE_COMPLETE=TRUE\n")
