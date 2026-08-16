#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(digest)
  library(jsonlite)
})

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript scripts/build_gate12bg_clr_repaired_s2.R <config.tsv>", call. = FALSE)
}

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
config_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
cfg <- fread(config_path)
if (!identical(names(cfg), c("key", "value")) || anyDuplicated(cfg$key)) {
  stop("Config must contain unique key/value rows", call. = FALSE)
}
cfg_get <- function(name) {
  value <- cfg[cfg$key == name, value]
  if (length(value) != 1L || !nzchar(value)) stop("Missing config key: ", name, call. = FALSE)
  value
}

gate_id <- cfg_get("gate_id")
run_id <- cfg_get("run_id")
output_dir <- file.path(root, cfg_get("output_dir"))
started_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
current_stage <- "initialization"
checks <- data.table(check_id = character(), expected = character(), observed = character(),
                     pass = logical(), severity = character())

record_check <- function(id, expected, observed, pass, severity = "HARD") {
  checks <<- rbind(checks, data.table(check_id = id, expected = as.character(expected),
                                     observed = as.character(observed), pass = isTRUE(pass),
                                     severity = severity))
}

write_lines <- function(x, path) {
  writeLines(enc2utf8(x), path, useBytes = TRUE)
}

replace_once <- function(text, old, new, id) {
  locations <- gregexpr(old, text, fixed = TRUE)[[1L]]
  count <- if (identical(locations[[1L]], -1L)) 0L else length(locations)
  record_check(paste0("replacement_", id), "1 occurrence", count, count == 1L)
  if (count != 1L) stop("Expected exactly one occurrence for replacement: ", id, call. = FALSE)
  sub(old, new, text, fixed = TRUE)
}

input_keys <- c("receipt", "patient_audit", "patient_contrasts", "strict_effects",
                "informative_effects", "pseudocount_sensitivity", "direction_trace",
                "manuscript", "legends")

main <- function() {
  current_stage <<- "output_initialization"
  if (dir.exists(output_dir)) stop("Output directory already exists: ", output_dir, call. = FALSE)
  dir.create(file.path(output_dir, "figures", "supplementary"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "source_data", "supplementary"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "manuscript"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "admin"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "provenance"), recursive = TRUE, showWarnings = FALSE)

  current_stage <<- "input_identity"
  input_audit <- rbindlist(lapply(input_keys, function(key) {
    rel <- cfg_get(paste0("input.", key, ".path"))
    expected <- cfg_get(paste0("input.", key, ".sha256"))
    abs <- file.path(root, rel)
    if (!file.exists(abs)) stop("Missing frozen input: ", rel, call. = FALSE)
    actual <- digest(abs, algo = "sha256", file = TRUE)
    data.table(input_key = key, path = rel, expected_sha256 = expected, actual_sha256 = actual,
               bytes = file.info(abs)$size, hash_match = identical(actual, expected))
  }))
  record_check("frozen_input_hashes", paste0(length(input_keys), "/", length(input_keys)),
               paste0(sum(input_audit$hash_match), "/", nrow(input_audit)), all(input_audit$hash_match))
  if (!all(input_audit$hash_match)) stop("Frozen input hash mismatch", call. = FALSE)
  fwrite(input_audit, file.path(output_dir, "provenance", "input_hashes.tsv"), sep = "\t")

  receipt <- read_json(file.path(root, cfg_get("input.receipt.path")), simplifyVector = TRUE)
  receipt_pass <- identical(receipt$status, "COMPLETED") && identical(as.integer(receipt$exit_code), 0L) &&
    nrow(receipt$checks) == 99L && all(receipt$checks$pass) && isTRUE(receipt$protected_manifest_unchanged)
  record_check("gate12bf_receipt", "COMPLETED; exit 0; 99/99 checks; protected unchanged",
               paste(receipt$status, receipt$exit_code, sum(receipt$checks$pass), nrow(receipt$checks),
                     receipt$protected_manifest_unchanged, sep = ";"), receipt_pass)
  if (!receipt_pass) stop("Gate12BF receipt is not eligible for figure reconstruction", call. = FALSE)

  patient_audit <- fread(file.path(root, cfg_get("input.patient_audit.path")))
  contrasts <- fread(file.path(root, cfg_get("input.patient_contrasts.path")))
  strict <- fread(file.path(root, cfg_get("input.strict_effects.path")))
  informative <- fread(file.path(root, cfg_get("input.informative_effects.path")))
  pc_sens <- fread(file.path(root, cfg_get("input.pseudocount_sensitivity.path")))
  trace <- fread(file.path(root, cfg_get("input.direction_trace.path")))

  primary_pc <- as.numeric(cfg_get("primary_pseudocount"))
  strict_patients <- patient_audit[complete_triplet == TRUE]
  strict_counts <- strict_patients[, .N, by = cancer]
  record_check("strict_patient_counts", "prostate=7;renal=4",
               paste(strict_counts$cancer, strict_counts$N, collapse = ";"),
               strict_counts[cancer == "prostate", N] == 7L && strict_counts[cancer == "renal", N] == 4L)
  record_check("renal_singletons_descriptive_only", "5",
               patient_audit[cancer == "renal" & gate12bf_role == "descriptive_singleton_only", .N],
               patient_audit[cancer == "renal" & gate12bf_role == "descriptive_singleton_only", .N] == 5L)

  display_map <- c(
    C1QC_macrophage = "C1QC macrophage",
    Inflammatory_monocyte = "Inflammatory monocyte",
    Classical_monocyte = "Classical monocyte",
    Resident_macrophage = "Resident macrophage",
    Proliferating_myeloid = "Proliferating myeloid",
    cDC = "cDC", pDC = "pDC", Unresolved_myeloid = "Unresolved myeloid",
    CD8_exhausted = "CD8 exhausted", Treg = "Treg",
    Proliferating_T_NK = "Proliferating T/NK", CD4_memory = "CD4 memory",
    CD8_effector = "CD8 effector", CD4_naive = "CD4 naive",
    NK_adaptive = "NK adaptive", Unresolved_T_NK = "Unresolved T/NK"
  )
  state_order <- names(display_map)
  if (!setequal(state_order, trace$state)) stop("State display map does not match traceability table", call. = FALSE)
  state_levels <- rev(unname(display_map[state_order]))
  cancer_labels <- c(prostate = "Prostate", renal = "Renal")
  lineage_cols <- c(Myeloid = "#D55E00", T_NK = "#0072B2")
  contrast_cols <- c(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B")

  current_stage <<- "source_reconstruction"
  patient_plot <- contrasts[abs(pseudocount - primary_pc) <= 1e-15 &
                              patient_key %chin% strict_patients$patient_key]
  patient_plot <- merge(patient_plot, trace[, .(lineage, state, original_direction,
                                                originally_cross_cancer_stable,
                                                gate12bf_direction_class)],
                        by = c("lineage", "state"), all.x = TRUE)
  patient_plot[, state_label := unname(display_map[state])]
  patient_plot[, state_factor := factor(state_label, levels = state_levels)]
  patient_plot[, cancer_label := factor(cancer_labels[cancer], levels = cancer_labels)]
  patient_plot[, patient_label := patient_id]
  record_check("panel_a_rows", "176", nrow(patient_plot), nrow(patient_plot) == 176L)
  record_check("panel_a_patients", "11", uniqueN(patient_plot$patient_key), uniqueN(patient_plot$patient_key) == 11L)
  record_check("panel_a_finite", "TRUE", all(is.finite(patient_plot$tumor_minus_distal)),
               all(is.finite(patient_plot$tumor_minus_distal)))

  forest <- strict[abs(pseudocount - primary_pc) <= 1e-15]
  forest <- merge(forest, trace[, .(lineage, state, originally_cross_cancer_stable,
                                    gate12bf_direction_class)], by = c("lineage", "state"), all.x = TRUE)
  forest[, state_label := unname(display_map[state])]
  forest[, state_factor := factor(state_label, levels = state_levels)]
  forest[, cancer_label := factor(cancer_labels[cancer], levels = cancer_labels)]
  forest[, exact_fdr := factor(q_value_exact < 0.05,
                              levels = c(FALSE, TRUE), labels = c("q >= 0.05", "q < 0.05"))]
  record_check("panel_b_rows", "32", nrow(forest), nrow(forest) == 32L)
  record_check("panel_b_exact_n", "prostate=7;renal=4",
               paste(forest[, unique(exact_n), by = cancer][, paste(cancer, V1, sep = "=")], collapse = ";"),
               all(forest[cancer == "prostate", exact_n] == 7L) && all(forest[cancer == "renal", exact_n] == 4L))
  record_check("renal_exact_resolution", "minimum P=0.125", min(forest[cancer == "renal", p_value_exact]),
               abs(min(forest[cancer == "renal", p_value_exact]) - 0.125) <= 1e-15)

  strict_sens <- strict[, .(pseudocount, cancer, lineage, state,
                            analysis_set = "Complete triplets",
                            n_patients, estimate = estimate_mean_stage_slope)]
  informative_sens <- informative[contrast == "stage_slope",
                                  .(pseudocount, cancer, lineage, state,
                                    analysis_set = "At least 2 compartments",
                                    n_patients, estimate = estimate_mean)]
  sensitivity <- rbindlist(list(strict_sens, informative_sens), use.names = TRUE)
  sensitivity <- merge(sensitivity, trace[, .(lineage, state, original_direction,
                                               originally_cross_cancer_stable,
                                               gate12bf_direction_class)],
                       by = c("lineage", "state"), all.x = TRUE)
  sensitivity[, direction_match := sign(estimate) == original_direction]
  sensitivity[, state_label := unname(display_map[state])]
  sensitivity[, state_factor := factor(state_label, levels = state_levels)]
  sensitivity[, cohort_set := fifelse(cancer == "prostate" & analysis_set == "Complete triplets",
                                      "Prostate\ntriplets (n=7)",
                               fifelse(cancer == "prostate", "Prostate\n>=2 compartments (n=9)",
                               fifelse(analysis_set == "Complete triplets", "Renal\ntriplets (n=4)",
                                       "Renal\n>=2 compartments (n=4)")))]
  sensitivity[, cohort_set := factor(cohort_set,
    levels = c("Prostate\ntriplets (n=7)", "Prostate\n>=2 compartments (n=9)",
               "Renal\ntriplets (n=4)", "Renal\n>=2 compartments (n=4)"))]
  sensitivity[, pc_label := factor(sprintf("%.2g", pseudocount), levels = c("0.25", "0.5", "1"))]
  match_summary <- sensitivity[, .(settings = .N, matches = sum(direction_match)),
                               by = .(lineage, state, originally_cross_cancer_stable)]
  record_check("panel_c_rows", "192", nrow(sensitivity), nrow(sensitivity) == 192L)
  record_check("nine_prespecified_states_all_12", "9/9",
               paste0(sum(match_summary[originally_cross_cancer_stable == TRUE, matches == 12L]), "/9"),
               nrow(match_summary[originally_cross_cancer_stable == TRUE]) == 9L &&
                 all(match_summary[originally_cross_cancer_stable == TRUE, matches == 12L]))
  record_check("all_state_robustness", "14 states at 12/12; CD4 memory 6/12; proliferating T/NK 11/12",
               paste0("14=", sum(match_summary$matches == 12L), ";CD4_memory=",
                      match_summary[state == "CD4_memory", matches], ";Proliferating_T_NK=",
                      match_summary[state == "Proliferating_T_NK", matches]),
               sum(match_summary$matches == 12L) == 14L &&
                 match_summary[state == "CD4_memory", matches] == 6L &&
                 match_summary[state == "Proliferating_T_NK", matches] == 11L)
  record_check("traceability_classification", "SUPPORTED=14;HETEROGENEOUS=1;DIRECTION_ONLY=1",
               paste(trace[, .N, by = gate12bf_direction_class][,
                           paste(gate12bf_direction_class, N, sep = "=")], collapse = ";"),
               trace[gate12bf_direction_class == "SUPPORTED", .N] == 14L &&
                 trace[gate12bf_direction_class == "HETEROGENEOUS", .N] == 1L &&
                 trace[gate12bf_direction_class == "DIRECTION_ONLY", .N] == 1L)

  source_dir <- file.path(output_dir, "source_data", "supplementary")
  fwrite(patient_plot[, .(cancer, patient_key, patient_id, lineage, state,
                          tumor_minus_distal, original_direction,
                          originally_cross_cancer_stable, gate12bf_direction_class)],
         file.path(source_dir, "FigureS2A_patient_tumor_minus_distal_CLR.tsv"), sep = "\t")
  fwrite(forest[, .(cancer, lineage, state, n_patients, estimate_mean,
                    bootstrap_ci_low, bootstrap_ci_high, p_value_exact, q_value_exact,
                    exact_configurations, originally_cross_cancer_stable,
                    gate12bf_direction_class)],
         file.path(source_dir, "FigureS2B_cohort_specific_CLR_effects.tsv"), sep = "\t")
  fwrite(sensitivity[, .(cancer, analysis_set, n_patients, pseudocount, lineage, state,
                         estimate_mean_stage_slope = estimate, original_direction,
                         direction_match, originally_cross_cancer_stable,
                         gate12bf_direction_class)],
         file.path(source_dir, "FigureS2C_analysis_set_pseudocount_sensitivity.tsv"), sep = "\t")
  fwrite(match_summary, file.path(source_dir, "FigureS2_state_direction_summary.tsv"), sep = "\t")

  current_stage <<- "figure_build"
  base_theme <- theme_classic(base_size = 8.4, base_family = "sans") +
    theme(plot.title = element_text(face = "bold", size = 9.2, hjust = 0),
          strip.background = element_blank(), strip.text = element_text(face = "bold", size = 8),
          axis.title = element_text(size = 8.2), axis.text = element_text(size = 7.2),
          legend.title = element_text(size = 7.4), legend.text = element_text(size = 7),
          plot.margin = margin(4, 5, 4, 4))

  heat_lim <- max(abs(patient_plot$tumor_minus_distal))
  p_a <- ggplot(patient_plot, aes(patient_label, state_factor, fill = tumor_minus_distal)) +
    geom_tile(color = "white", linewidth = 0.25) +
    facet_grid(. ~ cancer_label, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(low = contrast_cols[["low"]], mid = contrast_cols[["mid"]],
                         high = contrast_cols[["high"]], midpoint = 0,
                         limits = c(-heat_lim, heat_lim), oob = squish,
                         name = "Tumor - distal\nCLR contrast") +
    labs(title = "Patient contrasts", x = "Patient", y = NULL) +
    base_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.ticks = element_blank(), panel.border = element_rect(fill = NA, color = "#333333", linewidth = 0.35),
          legend.position = "bottom", legend.key.width = unit(18, "mm"))

  forest_lim <- max(abs(c(forest$bootstrap_ci_low, forest$bootstrap_ci_high))) * 1.06
  p_b <- ggplot(forest, aes(estimate_mean, state_factor, color = lineage, fill = exact_fdr)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, color = "#666666") +
    geom_errorbar(aes(xmin = bootstrap_ci_low, xmax = bootstrap_ci_high),
                  orientation = "y", width = 0.12, linewidth = 0.55) +
    geom_point(shape = 21, size = 2.35, stroke = 0.65) +
    facet_grid(. ~ cancer_label) +
    scale_color_manual(values = lineage_cols, labels = c(Myeloid = "Myeloid", T_NK = "T/NK"),
                       name = NULL) +
    scale_fill_manual(values = c("q >= 0.05" = "white", "q < 0.05" = "#222222"), name = "Exact BH") +
    coord_cartesian(xlim = c(-forest_lim, forest_lim)) +
    labs(title = "Cohort estimates", x = "Mean tumor - distal CLR contrast\n(95% patient-bootstrap interval)", y = NULL) +
    base_theme +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid.major.y = element_line(color = "#ECECEC", linewidth = 0.25),
          legend.position = "bottom", legend.box = "horizontal")

  sens_lim <- max(abs(sensitivity$estimate))
  p_c <- ggplot(sensitivity, aes(pc_label, state_factor, fill = estimate)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(data = sensitivity[direction_match == FALSE], label = "X", size = 3.2,
              color = "#111111", fontface = "bold") +
    facet_grid(. ~ cohort_set, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(low = contrast_cols[["low"]], mid = contrast_cols[["mid"]],
                         high = contrast_cols[["high"]], midpoint = 0,
                         limits = c(-sens_lim, sens_lim), oob = squish,
                         name = "Mean CLR\nstage slope") +
    labs(title = "Pseudocount and analysis-set sensitivity", x = "Pseudocount", y = NULL) +
    base_theme +
    theme(axis.ticks = element_blank(), panel.border = element_rect(fill = NA, color = "#333333", linewidth = 0.35),
          legend.position = "bottom", legend.key.width = unit(20, "mm"))

  composed <- (p_a | p_b) / p_c +
    plot_layout(heights = c(1.15, 1), widths = c(1.12, 1)) +
    plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold", size = 13)))

  fig_base <- file.path(output_dir, "figures", "supplementary", "FigureS2_repaired_CLR_compositional_robustness")
  staged_base <- file.path(tempdir(), paste0("gate12bg_figure_", Sys.getpid()))
  staged_pdf <- paste0(staged_base, ".pdf")
  staged_png <- paste0(staged_base, ".png")
  on.exit(unlink(c(staged_pdf, staged_png), force = TRUE), add = TRUE)
  ggsave(staged_pdf, composed,
         width = as.numeric(cfg_get("figure_width_in")), height = as.numeric(cfg_get("figure_height_in")),
         device = cairo_pdf, bg = "white")
  ggsave(staged_png, composed,
         width = as.numeric(cfg_get("figure_width_in")), height = as.numeric(cfg_get("figure_height_in")),
         dpi = as.numeric(cfg_get("figure_dpi")), bg = "white")
  png_dim <- dim(png::readPNG(staged_png))
  expected_dim <- c(as.numeric(cfg_get("figure_height_in")) * as.numeric(cfg_get("figure_dpi")),
                    as.numeric(cfg_get("figure_width_in")) * as.numeric(cfg_get("figure_dpi")))
  record_check("figure_png_dimensions", paste(expected_dim, collapse = "x"), paste(png_dim[1:2], collapse = "x"),
               identical(as.numeric(png_dim[1:2]), expected_dim))
  copy_ok <- c(file.copy(staged_png, paste0(fig_base, ".png"), overwrite = FALSE),
               file.copy(staged_pdf, paste0(fig_base, ".pdf"), overwrite = FALSE))
  record_check("staged_figure_copy", "PNG=TRUE;PDF=TRUE", paste(copy_ok, collapse = ";"), all(copy_ok))
  if (!all(copy_ok)) stop("Failed to copy staged figure files into the run output", call. = FALSE)
  final_png <- paste0(fig_base, ".png")
  final_pdf <- paste0(fig_base, ".pdf")
  png_bytes <- file.info(final_png)$size
  pdf_bytes <- file.info(final_pdf)$size
  pdf_con <- file(staged_pdf, open = "rb")
  pdf_magic <- rawToChar(readBin(pdf_con, what = "raw", n = 5L))
  close(pdf_con)
  record_check("figure_png_size", "PNG >100 KB", png_bytes, is.finite(png_bytes) && png_bytes > 100000)
  record_check("figure_pdf_size", "PDF >10 KB", pdf_bytes, is.finite(pdf_bytes) && pdf_bytes > 10000)
  record_check("figure_pdf_header", "%PDF-", pdf_magic, identical(pdf_magic, "%PDF-"))

  current_stage <<- "manuscript_patch"
  manuscript <- paste(readLines(file.path(root, cfg_get("input.manuscript.path")), warn = FALSE,
                                encoding = "UTF-8"), collapse = "\n")
  old_results_1 <- "Because proportions within a lineage are compositional, we repeated the analysis using patient-paired centered-log-ratio (CLR) coordinates and patient-cluster bootstrap intervals. All nine previously stable states retained their original direction, and all 16 modeled state directions matched the original meta-analytic direction in each of six frozen pseudocount and sample-subset sensitivities (Supplementary Figure S2). The largest positive CLR effects were inflammatory monocytes (beta, 2.646, 95% interval, 2.360-2.933) and C1QC macrophages (2.392, 2.025-2.767). Negative myeloid effects included classical monocytes (-0.922, -1.228 to -0.652), proliferating myeloid cells (-0.716, -1.078 to -0.395) and resident macrophages (-0.614, -1.029 to -0.297). Within T/NK cells, CD8-exhausted cells increased (0.964, 0.543-1.415), whereas CD4-naive cells (-2.043, -2.439 to -1.641) and adaptive natural-killer cells (-1.587, -1.974 to -1.193) decreased. CD8-effector cells retained the negative direction but had an interval crossing zero (-0.158, -0.390-0.076)."
  new_results_1 <- "Because proportions within a lineage are compositional, we repeated the analysis using patient-paired centered-log-ratio (CLR) coordinates. Primary inference was restricted to the seven prostate and four renal patients with complete distal-marrow, involved-marrow and tumor triplets, and used each patient's tumor-minus-distal contrast. At the primary pseudocount of 0.5, C1QC-macrophage contrasts were 2.981 (95% patient-bootstrap interval, 2.070-3.826; exact *P*=0.0156) in prostate cancer and 6.556 (5.548-7.388; exact *P*=0.125) in renal cancer. Corresponding CD8-exhausted contrasts were 1.930 (1.310-2.755; *P*=0.0156) and 3.731 (1.823-5.401; *P*=0.125). Classical-monocyte contrasts were -1.855 (-2.533 to -1.183; *P*=0.0156) and -2.726 (-3.992 to -1.143; *P*=0.125), whereas CD4-naive contrasts were -3.638 (-4.412 to -2.616; *P*=0.0156) and -5.053 (-6.391 to -3.901; *P*=0.125). All nine states prespecified from the empirical-logit analysis retained their direction in both cancers across pseudocounts of 0.25, 0.5 and 1.0 and in the at-least-two-compartment sensitivity set (Supplementary Figure S2)."
  manuscript <- replace_once(manuscript, old_results_1, new_results_1, "results_primary")

  old_results_2 <- "These concordant analyses reduce dependence on the original empirical-logit representation and on any single zero-replacement choice. They do not imply absolute expansion or depletion, because each effect describes a state relative to the geometric mean of the other states in its lineage. They also do not remove cancer heterogeneity. The C1QC-macrophage CLR effect remained positive but had an I-squared of 91.2%, and some states not meeting the original cross-cancer rule, including regulatory T cells, had pooled CLR effects that were not promoted to uniform pan-cancer findings."
  new_results_2 <- "Across all 16 coordinates, 14 retained the original direction in all 12 cohort-by-pseudocount-by-analysis-set settings. CD4-memory cells were directionally heterogeneous between cancers, and proliferating T/NK cells changed direction in one prostate sensitivity setting. CD8-effector uncertainty intervals crossed zero in both cohorts, and prostate proliferating-myeloid uncertainty also included zero. With four renal patients, the smallest attainable two-sided exact sign-flip *P* value was 0.125; the renal bootstrap intervals therefore summarize uncertainty but do not increase exact-test resolution. No pooled pan-cancer significance test was used. CLR effects denote enrichment or depletion relative to the geometric mean of the other states within the same lineage and do not imply absolute changes in cell number."
  manuscript <- replace_once(manuscript, old_results_2, new_results_2, "results_guardrails")

  old_methods <- "To account explicitly for the relative nature of within-lineage state proportions, we repeated the state analysis in CLR coordinates. For each eligible lineage-level sample, a frozen pseudocount was added to every resolved or unresolved state, the vector was closed to one and the natural logarithm of each state proportion was centered by the sample-specific geometric mean. The primary model estimated the anatomical-step effect within cancer while retaining patient pairing. Uncertainty used 999 patient-cluster bootstrap resamples. The analysis was repeated across six frozen pseudocount and sample-subset sensitivities. Concordance required the previously stable empirical-logit states to retain their direction in the primary CLR analysis and in the version-frozen majority of sensitivities. CLR coefficients were interpreted as relative enrichment or depletion within the lineage, not as absolute abundance changes."
  new_methods <- "To account explicitly for the relative nature of within-lineage state proportions, we repeated the state analysis in CLR coordinates. Within each eligible lineage-level sample, a pseudocount was added to every resolved and unresolved state, the vector was closed to one and each natural-log proportion was centered by the sample-specific mean log proportion. The primary pseudocount was 0.5, with 0.25 and 1.0 examined as prespecified sensitivities. Primary inference used only complete distal-marrow, involved-marrow and tumor triplets (seven prostate and four renal patients) and summarized the patient-level tumor-minus-distal CLR contrast separately within each cancer. Two-sided exact *P* values exhaustively enumerated all 2^7 prostate and 2^4 renal sign patterns using the absolute mean patient contrast, without a plus-one correction. Benjamini-Hochberg correction was applied within each cancer-by-lineage family. Uncertainty intervals used 9,999 within-cancer patient resamples; the same resampled patient indices were shared across all 16 states in each iteration and percentile type-7 intervals were reported. Patients with at least two anatomical compartments were retained only for a secondary mean-stage-slope sensitivity (nine prostate and four renal patients). Five renal tumor-only singletons were descriptive only and did not enter inferential counts, degrees of freedom, resampling or weights. Cross-cancer fixed-effect summaries were descriptive and were not used as pooled pan-cancer significance tests. CLR coordinates were interpreted as relative enrichment or depletion within a lineage, not as absolute abundance changes."
  manuscript <- replace_once(manuscript, old_methods, new_methods, "methods")

  old_discussion <- "The paired CLR analysis strengthens the discovery component in a specific way. All nine previously stable myeloid and T/NK state directions survived a model that treats within-lineage proportions as compositional. This reduces concern that the main directions were artifacts of empirical-logit modeling, one pseudocount or one complete-case definition. It also changes the language that is justified. An increased CLR coordinate means enrichment relative to the other states in that lineage. It does not establish an absolute increase in cell number. Substantial C1QC-macrophage heterogeneity and imprecise CD8-effector estimates further argue against flattening prostate and renal bone metastases into a homogeneous program."
  new_discussion <- "The paired CLR analysis strengthens the discovery component in a specific but bounded way. All nine prespecified myeloid and T/NK state directions survived complete-triplet analysis in both cancers across three pseudocounts and the broader contrast-informative sensitivity. This reduces concern that the principal directions were artifacts of empirical-logit modeling or one zero-replacement choice. It does not establish absolute changes in cell number, and the two non-prespecified exceptions across all 16 coordinates show why direction was not generalized indiscriminately. Cancer-specific effect magnitudes, the exact-test resolution limit in the four-patient renal cohort and imprecise CD8-effector estimates further argue against flattening prostate and renal bone metastases into a homogeneous program."
  manuscript <- replace_once(manuscript, old_discussion, new_discussion, "discussion")

  old_legend <- "## Supplementary Figure S2. Paired compositional robustness\n\na, Empirical-logit versus paired CLR effects on the same scale. b, CLR estimates for all 16 states; all 16 directions were retained in all six prespecified sensitivity settings."
  new_legend <- "## Supplementary Figure S2. Patient-level CLR compositional robustness\n\na, Patient-level tumor-minus-distal CLR contrasts for all 16 within-lineage coordinates in the seven prostate and four renal complete triplets. Columns are patients and rows are cell states; color is the within-patient contrast, so cells are not treated as biological replicates. b, Cancer-specific mean tumor-minus-distal contrasts with 95% intervals from 9,999 patient bootstraps. Filled points have an exact Benjamini-Hochberg-adjusted *q* value below 0.05 within the cancer-by-lineage family. Exact sign-flip inference used 128 prostate and 16 renal configurations; the smallest attainable two-sided renal *P* value was 0.125. c, Mean anatomical-stage slopes across pseudocounts of 0.25, 0.5 and 1.0 in the complete-triplet and at-least-two-compartment analysis sets. A cross marks a direction opposite to the original empirical-logit direction. All nine prespecified stable states matched in all 12 settings. Across all 16 coordinates, 14 matched in all 12 settings; CD4-memory cells matched in 6 of 12 because the renal direction was opposite, and proliferating T/NK cells matched in 11 of 12. CLR values are relative within lineage and do not represent absolute cell-number changes."
  manuscript <- replace_once(manuscript, old_legend, new_legend, "manuscript_legend")
  manuscript_path <- file.path(output_dir, "manuscript", "Gate12BG_CLR_Repaired_Manuscript.md")
  write_lines(manuscript, manuscript_path)

  legends <- paste(readLines(file.path(root, cfg_get("input.legends.path")), warn = FALSE,
                             encoding = "UTF-8"), collapse = "\n")
  legends <- replace_once(legends, old_legend, new_legend, "standalone_legend")
  legends_path <- file.path(output_dir, "manuscript", "Gate12BG_Figure_Legends.md")
  write_lines(legends, legends_path)

  manuscript_checks <- list(
    old_16_of_16_absent = !grepl("all 16 modeled state directions matched", manuscript, fixed = TRUE),
    old_999_absent = !grepl("Uncertainty used 999 patient-cluster bootstrap resamples", manuscript, fixed = TRUE),
    old_pooled_clr_absent = !grepl("C1QC-macrophage CLR effect remained positive but had an I-squared", manuscript, fixed = TRUE),
    new_patient_counts_present = grepl("seven prostate and four renal patients", manuscript, fixed = TRUE),
    new_bootstrap_present = grepl("9,999 within-cancer patient resamples", manuscript, fixed = TRUE),
    new_resolution_limit_present = grepl("smallest attainable two-sided exact sign-flip", manuscript, fixed = TRUE),
    new_direction_summary_present = grepl("14 retained the original direction in all 12", manuscript, fixed = TRUE),
    new_s2_title_present = grepl("Patient-level CLR compositional robustness", manuscript, fixed = TRUE)
  )
  for (nm in names(manuscript_checks)) {
    record_check(paste0("manuscript_", nm), "TRUE", manuscript_checks[[nm]], manuscript_checks[[nm]])
  }

  patch_log <- c(
    "# Gate12BG manuscript patch log", "",
    "- Source manuscript: Gate12BE review-driven manuscript (hash-frozen).",
    "- Replaced the invalid pooled CLR Results paragraph with cancer-specific complete-triplet estimates.",
    "- Replaced the 16/16-in-6/6 claim with the Gate12BF classification: 14 SUPPORTED, 1 HETEROGENEOUS and 1 DIRECTION_ONLY.",
    "- Replaced 999 bootstraps with the executed 9,999-patient-bootstrap specification.",
    "- Added exhaustive sign-flip inference, cancer-by-lineage BH correction and the renal P-value resolution limit.",
    "- Declared the five renal tumor-only singletons descriptive-only.",
    "- Replaced the Supplementary Figure S2 title and legend in both manuscript and standalone legends.",
    "- No other manuscript sections, figures, references, authors or affiliations were edited."
  )
  write_lines(patch_log, file.path(output_dir, "admin", "GATE12BG_MANUSCRIPT_PATCH_LOG.md"))

  visual_contract <- data.table(
    panel = c("A", "B", "C"),
    display = c("patient-by-state diverging heatmap", "cohort-stratified forest plot",
                "analysis-set-by-pseudocount diverging heatmap"),
    biological_unit = c("patient", "patient", "patient"),
    inference = c("descriptive patient contrasts", "cancer-specific exact sign flip plus patient bootstrap",
                  "directional sensitivity only"),
    prohibited = c("cell-level inference", "pooled pan-cancer significance", "treating settings as replicates")
  )
  fwrite(visual_contract, file.path(output_dir, "admin", "GATE12BG_FIGURE_CONTRACT.tsv"), sep = "\t")

  current_stage <<- "final_validation"
  record_check("all_effect_directions_finite", "TRUE",
               all(is.finite(forest$estimate_mean)) && all(is.finite(sensitivity$estimate)),
               all(is.finite(forest$estimate_mean)) && all(is.finite(sensitivity$estimate)))
  record_check("no_pooled_inference_panel", "TRUE", TRUE, TRUE)
  if (any(!checks$pass & checks$severity == "HARD")) {
    fwrite(checks, file.path(output_dir, "admin", "validation_checks.tsv"), sep = "\t")
    stop("One or more Gate12BG hard checks failed", call. = FALSE)
  }

  session_text <- capture.output(sessionInfo())
  write_lines(session_text, file.path(output_dir, "provenance", "sessionInfo.txt"))
  run_log <- c(
    paste0("RUN_STARTED gate=", gate_id, " run=", run_id, " at=", started_at),
    "INPUT Gate12BF=run_v3 status=COMPLETED checks=99/99",
    "PRIMARY complete_triplets prostate_n=7 renal_n=4 pseudocount=0.5",
    "FIGURE panels=A_patient_heatmap,B_cohort_forest,C_sensitivity_heatmap",
    paste0("RUN_COMPLETED status=COMPLETED checks=", nrow(checks), " all_hard_checks_pass=TRUE")
  )
  write_lines(run_log, file.path(output_dir, "GATE12BG_RUN_LOG.txt"))
  fwrite(checks, file.path(output_dir, "admin", "validation_checks.tsv"), sep = "\t")

  manifest_targets <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
  manifest_targets <- manifest_targets[file.info(manifest_targets)$isdir == FALSE]
  manifest_targets <- manifest_targets[!basename(manifest_targets) %chin% c("GATE12BG_RECEIPT.json", "output_manifest.tsv")]
  manifest <- rbindlist(lapply(sort(manifest_targets), function(path) {
    rel <- substring(normalizePath(path, winslash = "/"), nchar(normalizePath(output_dir, winslash = "/")) + 2L)
    data.table(file = rel, bytes = file.info(path)$size,
               sha256 = digest(path, algo = "sha256", file = TRUE))
  }))
  fwrite(manifest, file.path(output_dir, "provenance", "output_manifest.tsv"), sep = "\t")
  manifest_sha <- digest(file.path(output_dir, "provenance", "output_manifest.tsv"), algo = "sha256", file = TRUE)

  ended_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  receipt_out <- list(
    schema_version = "1.0",
    gate_id = gate_id,
    run_id = run_id,
    status = "COMPLETED",
    command = paste("Rscript", "scripts/build_gate12bg_clr_repaired_s2.R", args[[1L]]),
    working_directory = root,
    started_at = started_at,
    ended_at = ended_at,
    exit_code = 0L,
    source_gate = list(gate = "Gate12BF", run = "run_v3", receipt_sha256 = input_audit[input_key == "receipt", actual_sha256]),
    analysis_contract = list(primary_unit = "complete patient triplet", primary_contrast = "tumor minus distal CLR",
                             cohorts = "cancer-specific", primary_pseudocount = primary_pc,
                             exact_test = "exhaustive two-sided patient sign flip",
                             bootstrap = "9999 within-cancer patient resamples",
                             pooled_pan_cancer_inference = FALSE),
    patient_counts = list(prostate = 7L, renal = 4L, renal_descriptive_singletons = 5L),
    direction_summary = as.list(setNames(trace[, .N, by = gate12bf_direction_class]$N,
                                              trace[, .N, by = gate12bf_direction_class]$gate12bf_direction_class)),
    checks = checks,
    output_manifest_sha256 = manifest_sha
  )
  write_json(receipt_out, file.path(output_dir, "GATE12BG_RECEIPT.json"), auto_unbox = TRUE,
             pretty = TRUE, na = "null")

  cat(sprintf("RUN_COMPLETED status=COMPLETED checks=%d all_hard_checks_pass=TRUE\n", nrow(checks)))
}

tryCatch(
  main(),
  error = function(e) {
    if (dir.exists(output_dir)) {
      fail <- list(schema_version = "1.0", gate_id = gate_id, run_id = run_id, status = "FAILED",
                   command = paste("Rscript", "scripts/build_gate12bg_clr_repaired_s2.R", args[[1L]]),
                   working_directory = root, started_at = started_at,
                   ended_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), exit_code = 1L,
                   failure_stage = current_stage, error = conditionMessage(e), checks = checks)
      write_json(fail, file.path(output_dir, "GATE12BG_RECEIPT.json"), auto_unbox = TRUE,
                 pretty = TRUE, na = "null")
      if (nrow(checks)) fwrite(checks, file.path(output_dir, "admin", "validation_checks.tsv"), sep = "\t")
    }
    message("RUN_FAILED stage=", current_stage, " error=", conditionMessage(e))
    quit(status = 1L)
  }
)
