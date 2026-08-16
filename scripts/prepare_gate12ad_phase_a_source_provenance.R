#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(jsonlite)
  library(png)
})

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) normalizePath(args[[1L]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
out_rel <- if (length(args) >= 2L) args[[2L]] else "results/gate12ad_figure_restructure/phase_a_source_provenance"
out_dir <- file.path(root, out_rel)
source_dir <- file.path(out_dir, "source_data")
provenance_dir <- file.path(out_dir, "provenance")
admin_dir <- file.path(out_dir, "admin")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(provenance_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(admin_dir, recursive = TRUE, showWarnings = FALSE)

path <- function(x) file.path(root, x)
rel_path <- function(x) {
  z <- normalizePath(x, mustWork = TRUE)
  prefix <- paste0(root, .Platform$file.sep)
  if (!startsWith(z, prefix)) stop("Path is outside project root: ", z)
  substring(z, nchar(prefix) + 1L)
}
sha256_file <- function(x) digest(x, file = TRUE, algo = "sha256")
write_tsv <- function(x, filename, dir = source_dir, compress = FALSE) {
  target <- file.path(dir, filename)
  fwrite(x, target, sep = "\t", quote = FALSE, na = "NA", compress = if (compress) "gzip" else "none")
  invisible(target)
}
assert <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)

input_rel <- c(
  "results/gate12a_integrated_atlas/integrated_umap_coordinates.tsv.gz",
  "results/gate12a_integrated_atlas/broad_class_composition.tsv",
  "results/gate12a_integrated_atlas/marker_expression_summary.tsv",
  "results/gate12b_cell_states/patient_state_composition.tsv",
  "results/gate9b_validation/input_audit/oep005136_mbone_crosswalk.tsv",
  "results/gate9b_validation/external_validation_v1/patient_assignments.tsv",
  "results/gate9b_validation/external_validation_v1/k3_transfer_agreement.tsv",
  "scripts/run_gate9b_external_validation.R",
  "scripts/prepare_gate12ad_phase_a_source_provenance.R",
  "scripts/audit_gate12ad_phase_a.py",
  "results/gate12ab_minor_revision_closure/analysis/block_context/axis1_block_lopo_predictions.tsv",
  "results/gate12ab_minor_revision_closure/analysis/block_context/axis1_block_lopo_summary.tsv",
  "results/gate12ab_minor_revision_closure/analysis/block_context/axis1_block_lopo_paired_direction.tsv",
  "results/gate12ab_minor_revision_closure/source_data/Figure5_GSE266330_all_patient_scores.tsv",
  "results/gate12ab_minor_revision_closure/source_data/Figure5_hallmark_calibration_representatives.tsv",
  "results/gate12ab_minor_revision_closure/source_data/Figure7_all_axis1_paired_differences.tsv",
  "results/gate12v3_spatial_geometry/spot_spatial_geometry.tsv.gz",
  "results/gate12ab_minor_revision_closure/source_data/Figure6_distance_curves.tsv",
  "results/gate12ab_minor_revision_closure/source_data/Figure6_section_effects.tsv",
  "results/gate12ab_minor_revision_closure/source_data/Figure6_neighborhood_class_effects.tsv",
  "results/gate12e_histology_assets/GSM9564255_Bone_ST_24664_tissue_lowres_image.png",
  "results/gate12e_histology_assets/GSM9564255_Bone_ST_24664_scalefactors_json.json",
  "results/gate12e_histology_assets/GSM9564256_Bone_ST_24665_tissue_lowres_image.png",
  "results/gate12e_histology_assets/GSM9564256_Bone_ST_24665_scalefactors_json.json",
  "results/gate12e_histology_assets/GSM9564257_Bone_ST_24666_tissue_lowres_image.png",
  "results/gate12e_histology_assets/GSM9564257_Bone_ST_24666_scalefactors_json.json",
  "results/gate12e_histology_assets/GSM9564258_Bone_ST_24667_tissue_lowres_image.png",
  "results/gate12e_histology_assets/GSM9564258_Bone_ST_24667_scalefactors_json.json",
  "results/gate12ab_minor_revision_closure/GATE12AB_BASE_SOURCE_MANIFEST.tsv",
  "results/gate12ab_minor_revision_closure/reproducibility/GATE12AB_RELEASE_MANIFEST.tsv",
  "results/gate12ab_minor_revision_closure/admin/GATE12AB_RELEASE_BUILD_RECEIPT.json"
)
input_abs <- path(input_rel)
assert(all(file.exists(input_abs)), paste("Missing inputs:", paste(input_rel[!file.exists(input_abs)], collapse = ", ")))

input_manifest <- data.table(
  relative_path = input_rel,
  bytes = as.numeric(file.info(input_abs)$size),
  sha256 = vapply(input_abs, sha256_file, character(1L))
)
write_tsv(input_manifest, "GATE12AD_A_INPUT_MANIFEST.tsv", provenance_dir)

## Figure 1: direct panel-level data from frozen discovery assets.
class_order <- c("T_NK", "Myeloid", "B", "Progenitor", "Stromal", "Endothelial",
                 "Osteoclast", "Osteoblast", "Malignant", "Erythroid", "Unassigned")
class_label <- c(T_NK = "T/NK", Myeloid = "Myeloid", B = "B cell", Progenitor = "Progenitor",
                 Stromal = "Stromal", Endothelial = "Endothelial", Osteoclast = "Osteoclast",
                 Osteoblast = "Osteoblast", Malignant = "Malignant", Erythroid = "Erythroid",
                 Unassigned = "Unassigned")

coords_raw <- fread(path(input_rel[[1L]]), check.names = TRUE)
assert(all(c("barcode", "umap_1", "umap_2", "accession", "cancer", "sample_id", "patient_id",
             "compartment", "broad_class", "harmonized_state", "label_source") %in% names(coords_raw)),
       "Unexpected integrated UMAP schema")
coords <- coords_raw[, .(
  cell_id = barcode,
  umap_1,
  umap_2,
  accession,
  cancer,
  sample_id,
  patient_id,
  compartment,
  broad_class,
  broad_class_label = unname(class_label[broad_class]),
  harmonized_state,
  label_source
)]
coords[, broad_class_order := match(broad_class, class_order)]
assert(nrow(coords) == 107886L && uniqueN(coords$cell_id) == 107886L, "Figure 1 UMAP cell scope mismatch")
setorder(coords, broad_class_order, cell_id)
write_tsv(coords, "Figure1B_umap_coordinates.tsv.gz", compress = TRUE)

umap_labels <- coords[, .(
  label_x = median(umap_1),
  label_y = median(umap_2),
  n_cells = .N,
  n_samples = uniqueN(sample_id),
  n_patients = uniqueN(patient_id)
), by = .(broad_class, broad_class_label, broad_class_order)]
setorder(umap_labels, broad_class_order)
write_tsv(umap_labels, "Figure1B_umap_label_positions.tsv")

composition_raw <- fread(path(input_rel[[2L]]))
sample_meta <- unique(composition_raw[, .(accession, cancer, patient_id, sample_id, compartment, sample_total)])
assert(nrow(sample_meta) == 42L && uniqueN(sample_meta$patient_id) == 18L, "Discovery sample scope mismatch")
sample_meta[, cancer_order := match(cancer, c("prostate", "renal"))]
sample_meta[, compartment_order := match(compartment, c("distal", "involved", "tumor"))]
setorder(sample_meta, cancer_order, compartment_order, patient_id, sample_id)
sample_meta[, sample_panel_order := seq_len(.N)]

composition <- merge(
  CJ(sample_id = sample_meta$sample_id, broad_class = class_order, unique = TRUE),
  composition_raw[, .(sample_id, broad_class, N)],
  by = c("sample_id", "broad_class"), all.x = TRUE
)
composition[is.na(N), N := 0L]
composition <- merge(composition, sample_meta, by = "sample_id", all.x = TRUE)
composition[, fraction := N / sample_total]
composition[, broad_class_label := unname(class_label[broad_class])]
composition[, broad_class_order := match(broad_class, class_order)]
setorder(composition, sample_panel_order, broad_class_order)
assert(nrow(composition) == 42L * length(class_order), "Complete composition grid is incomplete")
assert(max(abs(composition[, sum(fraction), by = sample_id]$V1 - 1)) < 1e-12, "Composition fractions do not sum to one")
write_tsv(composition, "Figure1D_all_sample_broad_composition.tsv")

selected_broad <- c("T_NK", "Myeloid", "Stromal", "Malignant")
paired_fractions <- composition[broad_class %chin% selected_broad]
paired_fractions[, compartments_observed_for_patient := uniqueN(compartment), by = .(patient_id, broad_class)]
paired_fractions[, connected_trajectory := compartments_observed_for_patient > 1L]
setorder(paired_fractions, broad_class_order, cancer_order, patient_id, compartment_order)
write_tsv(paired_fractions, "Figure1E_patient_connected_broad_fractions.tsv")

marker_raw <- fread(path(input_rel[[3L]]))
gene_order <- unique(marker_raw$gene)
coverage <- composition[, .(
  n_cells = sum(N),
  n_samples_present = uniqueN(sample_id[N > 0]),
  n_patients_present = uniqueN(patient_id[N > 0])
), by = broad_class]
markers <- merge(marker_raw, coverage, by = "broad_class", all.x = TRUE)
markers[, mean_scaled_expression := {
  s <- sd(average_log_normalized)
  if (!is.finite(s) || s == 0) rep(0, .N) else (average_log_normalized - mean(average_log_normalized)) / s
}, by = gene]
markers[, detected_percent := 100 * detected_fraction]
markers[, marker_gene_order := match(gene, gene_order)]
markers[, broad_class_order := match(broad_class, class_order)]
markers[, broad_class_label := unname(class_label[broad_class])]
setorder(markers, marker_gene_order, broad_class_order)
assert(uniqueN(markers$gene) == 13L && uniqueN(markers$broad_class) == 11L && nrow(markers) == 143L,
       "Marker summary scope mismatch")
write_tsv(markers, "Figure1C_canonical_marker_dotplot.tsv")

## Figure 4: exact ordered patient assignments and permutation replay materials.
crosswalk <- fread(path("results/gate9b_validation/input_audit/oep005136_mbone_crosswalk.tsv"))
assignments <- fread(path("results/gate9b_validation/external_validation_v1/patient_assignments.tsv"))
primary_order <- crosswalk[primary_analysis == TRUE, .(
  analysis_order = seq_len(.N), archive_dir, sample_id, patient_id, sample_name,
  disease_name, cancer_code, sex, age, collection_date
)]
assignment_drop <- c("sex", "age", "collection_date", "primary_analysis", "exclusion_reason")
assignment_keep <- assignments[, setdiff(names(assignments), assignment_drop), with = FALSE]
ordered <- merge(primary_order, assignment_keep, by = c("patient_id", "sample_id", "archive_dir", "cancer_code"),
                 all.x = TRUE, sort = FALSE)
setorder(ordered, analysis_order)
assert(nrow(ordered) == 49L && !anyNA(ordered$transferred_label) && !anyNA(ordered$consensus_k3_cluster),
       "Ordered Figure 4 patient assignments are incomplete")
ordered[, included_in_ari := transferred_label != "unassigned"]
ordered[, transferred_label_display := fifelse(transferred_label == "Mphi_OC", "Macrophage/osteoclast",
                                        fifelse(transferred_label == "Mono", "Monocyte",
                                        fifelse(transferred_label == "Treg_Tex", "Treg/Tex", "Unassigned")))]
ordered[, consensus_cluster_display := paste("Consensus", consensus_k3_cluster)]
ordered_export <- ordered[, .(
  analysis_order, patient_id, sample_id, archive_dir, cancer_code, disease_name, sex, age,
  assigned_cells, unassigned_cells, assigned_fraction, transferred_label,
  transferred_label_display, consensus_k3_cluster, consensus_cluster_display,
  nearest_distance, second_distance, assignment_margin, included_in_ari
)]
write_tsv(ordered_export, "Figure4C_ordered_patient_assignments.tsv")

agreement_all <- ordered_export[, .(N = .N),
  by = .(transferred_label, transferred_label_display, consensus_k3_cluster, consensus_cluster_display)]
agreement_all[, row_total := sum(N), by = transferred_label]
agreement_all[, column_total := sum(N), by = consensus_k3_cluster]
agreement_all[, `:=`(row_fraction = N / row_total, column_fraction = N / column_total)]
setorder(agreement_all, transferred_label_display, consensus_k3_cluster)
write_tsv(agreement_all, "Figure4C_transfer_consensus_matrix_all_patients.tsv")

adjusted_rand <- function(a, b) {
  assert(length(a) == length(b) && length(a) >= 2L, "ARI vectors are invalid")
  tab <- table(a, b)
  choose2 <- function(x) x * (x - 1) / 2
  n2 <- choose2(sum(tab))
  index <- sum(choose2(tab))
  aa <- sum(choose2(rowSums(tab)))
  bb <- sum(choose2(colSums(tab)))
  expected <- aa * bb / n2
  upper <- 0.5 * (aa + bb)
  if (upper == expected) return(1)
  (index - expected) / (upper - expected)
}

ari_input <- ordered_export[included_in_ari == TRUE]
observed_ari <- adjusted_rand(ari_input$consensus_k3_cluster, ari_input$transferred_label)
permutation_seed <- 20260807L + 300L
permutations <- 10000L
set.seed(permutation_seed)
ari_null <- replicate(permutations, adjusted_rand(
  ari_input$consensus_k3_cluster,
  sample(ari_input$transferred_label, replace = FALSE)
))
permutation_p <- (1 + sum(ari_null >= observed_ari)) / (permutations + 1)
frozen_agreement <- fread(path("results/gate9b_validation/external_validation_v1/k3_transfer_agreement.tsv"))
assert(abs(observed_ari - frozen_agreement$ari) < 1e-15, "ARI does not replay")
assert(abs(permutation_p - frozen_agreement$permutation_p) < 1e-15, "Permutation P does not replay")
null_table <- data.table(
  iteration = seq_len(permutations),
  null_ari = ari_null,
  null_ge_observed = ari_null >= observed_ari
)
write_tsv(null_table, "Figure4C_ari_permutation_null.tsv.gz", compress = TRUE)
permutation_receipt <- data.table(
  test = "upper-tail label-permutation test of adjusted Rand index",
  observed_ari = observed_ari,
  n_included_patients = nrow(ari_input),
  excluded_unassigned_patients = nrow(ordered_export) - nrow(ari_input),
  permutations = permutations,
  seed = permutation_seed,
  exceedances_ge_observed = sum(ari_null >= observed_ari),
  empirical_p = permutation_p,
  p_formula = "(1 + sum(null_ari >= observed_ari)) / (B + 1)",
  order_source = "OEP005136 crosswalk primary_analysis rows in submitted archive order"
)
write_tsv(permutation_receipt, "Figure4C_ari_permutation_receipt.tsv", provenance_dir)

## Figure 5: all point-level data for the restructured main figure.
pred <- fread(path("results/gate12ab_minor_revision_closure/analysis/block_context/axis1_block_lopo_predictions.tsv"))
model_order <- c("broad_only", "myeloid_only", "t_nk_only", "broad_plus_depth")
model_label <- c(broad_only = "Broad block", myeloid_only = "Myeloid block",
                 t_nk_only = "T/NK block", broad_plus_depth = "Broad block + depth")
pred[, model_order := match(model, model_order)]
pred[, model_label := unname(model_label[model])]
setorder(pred, model_order, held_out_patient, sample_id)
assert(nrow(pred) == 164L && uniqueN(pred$model) == 4L && all(pred[, .N, by = model]$N == 41L),
       "LOPO prediction scope mismatch")
write_tsv(pred, "Figure5B_block_lopo_predictions.tsv")

metrics <- pred[, .(
  n_samples = .N,
  n_patients = uniqueN(patient_id),
  cross_validated_R2_training_mean = 1 - sum(residual^2) / sum(baseline_residual^2),
  spearman_rho = cor(observed_Axis1, predicted_Axis1, method = "spearman"),
  pearson_r = cor(observed_Axis1, predicted_Axis1),
  mean_absolute_error = mean(abs(residual)),
  root_mean_squared_error = sqrt(mean(residual^2))
), by = .(model, model_order, model_label)]
paired_direction <- fread(path("results/gate12ab_minor_revision_closure/analysis/block_context/axis1_block_lopo_paired_direction.tsv"))
metrics <- merge(metrics, paired_direction[, .(model, n_positive, n_pairs, positive_fraction, median_difference)],
                 by = "model", all.x = TRUE)
setorder(metrics, model_order)
frozen_metrics <- fread(path("results/gate12ab_minor_revision_closure/analysis/block_context/axis1_block_lopo_summary.tsv"))
metric_check <- merge(metrics, frozen_metrics, by = "model", suffixes = c("_new", "_frozen"))
assert(max(abs(metric_check$cross_validated_R2_training_mean_new - metric_check$cross_validated_R2_training_mean_frozen)) < 1e-15,
       "LOPO metric mismatch")
write_tsv(metrics, "Figure5B_block_lopo_metrics.tsv")

external_scores <- fread(path("results/gate12ab_minor_revision_closure/source_data/Figure5_GSE266330_all_patient_scores.tsv"))
write_tsv(external_scores, "Figure5A_GSE266330_patient_scores.tsv")
hallmark <- fread(path("results/gate12ab_minor_revision_closure/source_data/Figure5_hallmark_calibration_representatives.tsv"))
write_tsv(hallmark, "Figure5C_hallmark_calibration_representatives.tsv")
paired_axis <- fread(path("results/gate12ab_minor_revision_closure/source_data/Figure7_all_axis1_paired_differences.tsv"))
paired_axis[, comparison_label := fifelse(contrast == "bm_vs_normal_bone", "BM - matched normal bone", "BM - matched primary tumor")]
paired_axis[, positive_difference := difference > 0]
paired_axis[, `:=`(
  required_positive = fifelse(contrast == "bm_vs_normal_bone", 2L, 2L),
  required_evaluable = fifelse(contrast == "bm_vs_normal_bone", 2L, 3L)
)]
paired_axis[, observed_positive := sum(positive_difference), by = contrast]
paired_axis[, endpoint_pass := observed_positive >= required_positive]
write_tsv(paired_axis, "Figure5D_paired_axis1_differences.tsv")

## Figure 6: exact spatial overlay coordinates, common scale and image receipts.
spots <- fread(path("results/gate12v3_spatial_geometry/spot_spatial_geometry.tsv.gz"))
sample_manifest <- data.table(
  sample = paste0("GSM956425", 5:8),
  library = as.character(24664:24667)
)
full_limit <- unname(quantile(abs(spots$full_axis1), 0.98, na.rm = TRUE))
excluded_limit <- unname(quantile(abs(spots$malignant_excluded_axis1), 0.98, na.rm = TRUE))
shared_limit <- unname(quantile(c(abs(spots$full_axis1), abs(spots$malignant_excluded_axis1)), 0.98, na.rm = TRUE))

overlay_rows <- list()
asset_rows <- list()
scale_rows <- list()
for (i in seq_len(nrow(sample_manifest))) {
  sample_value <- sample_manifest$sample[[i]]
  prefix <- paste0(sample_value, "_Bone_ST_", sample_manifest$library[[i]])
  image_rel <- paste0("results/gate12e_histology_assets/", prefix, "_tissue_lowres_image.png")
  json_rel <- paste0("results/gate12e_histology_assets/", prefix, "_scalefactors_json.json")
  image_path <- path(image_rel)
  json_path <- path(json_rel)
  image <- readPNG(image_path)
  scale_json <- fromJSON(json_path)
  width <- dim(image)[2L]
  height <- dim(image)[1L]
  scale_factor <- scale_json$tissue_lowres_scalef
  section <- copy(spots[sample == sample_value])
  section[, `:=`(
    plot_x_lowres = pixel_col * scale_factor,
    plot_y_lowres = height - pixel_row * scale_factor,
    image_width_lowres = width,
    image_height_lowres = height,
    shared_score_limit = shared_limit,
    full_axis1_display = pmax(pmin(full_axis1, shared_limit), -shared_limit),
    malignant_excluded_axis1_display = pmax(pmin(malignant_excluded_axis1, shared_limit), -shared_limit),
    full_axis1_clipped = abs(full_axis1) > shared_limit,
    malignant_excluded_axis1_clipped = abs(malignant_excluded_axis1) > shared_limit,
    source_annotated_tumor = submitted_label == "Tumor"
  )]
  assert(min(section$plot_x_lowres) >= 0 && max(section$plot_x_lowres) <= width &&
           min(section$plot_y_lowres) >= 0 && max(section$plot_y_lowres) <= height,
         paste(sample_value, "contains coordinates outside the H&E image"))
  overlay_rows[[i]] <- section

  spot_55_px <- scale_json$spot_diameter_fullres * scale_factor
  bar_500_px <- spot_55_px * (500 / 55)
  scale_rows[[i]] <- data.table(
    sample = sample_value,
    library = sample_manifest$library[[i]],
    image_width_lowres_px = width,
    image_height_lowres_px = height,
    tissue_lowres_scalef = scale_factor,
    spot_diameter_fullres_px = scale_json$spot_diameter_fullres,
    pixels_per_55um_spot_lowres = spot_55_px,
    pixels_per_500um_lowres = bar_500_px,
    bar_fraction_of_image_width = bar_500_px / width,
    coordinate_bounds_pass = TRUE
  )
  asset_rows[[i]] <- data.table(
    sample = sample_value,
    library = sample_manifest$library[[i]],
    image_relative_path = image_rel,
    image_bytes = as.numeric(file.info(image_path)$size),
    image_sha256 = sha256_file(image_path),
    scalefactors_relative_path = json_rel,
    scalefactors_bytes = as.numeric(file.info(json_path)$size),
    scalefactors_sha256 = sha256_file(json_path)
  )
}
overlay <- rbindlist(overlay_rows)
setorder(overlay, sample, barcode)
assert(nrow(overlay) == 8190L && uniqueN(overlay$sample) == 4L, "Spatial spot scope mismatch")
write_tsv(overlay, "Figure6A_C_spot_overlay_coordinates_and_scores.tsv.gz", compress = TRUE)
write_tsv(rbindlist(asset_rows), "Figure6_histology_asset_index.tsv", provenance_dir)
write_tsv(rbindlist(scale_rows), "Figure6_scale_bar_audit.tsv", provenance_dir)

scale_audit <- rbindlist(list(
  data.table(layer = "Full Axis1", original_layer_abs_q98 = full_limit,
             shared_abs_q98_limit = shared_limit, n_finite = sum(is.finite(spots$full_axis1)),
             n_clipped_at_shared_limit = sum(abs(spots$full_axis1) > shared_limit, na.rm = TRUE),
             raw_min = min(spots$full_axis1, na.rm = TRUE), raw_max = max(spots$full_axis1, na.rm = TRUE)),
  data.table(layer = "Malignant-excluded Axis1", original_layer_abs_q98 = excluded_limit,
             shared_abs_q98_limit = shared_limit, n_finite = sum(is.finite(spots$malignant_excluded_axis1)),
             n_clipped_at_shared_limit = sum(abs(spots$malignant_excluded_axis1) > shared_limit, na.rm = TRUE),
             raw_min = min(spots$malignant_excluded_axis1, na.rm = TRUE),
             raw_max = max(spots$malignant_excluded_axis1, na.rm = TRUE))
))
scale_audit[, clipped_fraction := n_clipped_at_shared_limit / n_finite]
scale_audit[, statistics_use_unclipped_values := TRUE]
write_tsv(scale_audit, "Figure6_shared_color_scale_audit.tsv", provenance_dir)

image_processing <- data.table(
  operation = "uniform global histology lightening for score-overlay visibility",
  lightening_fraction = 0.12,
  transform = "display_rgb = source_rgb * (1 - 0.12) + 0.12",
  applied_to = "Full and malignant-excluded Axis1 overlay backgrounds",
  local_adjustment_applied = FALSE,
  local_masking_applied = FALSE,
  source_annotation_term = "source-annotated tumor spots/regions"
)
write_tsv(image_processing, "Figure6_image_processing_receipt.tsv", provenance_dir)

write_tsv(fread(path("results/gate12ab_minor_revision_closure/source_data/Figure6_distance_curves.tsv")),
          "Figure6D_distance_curves.tsv")
write_tsv(fread(path("results/gate12ab_minor_revision_closure/source_data/Figure6_section_effects.tsv")),
          "Figure6E_section_effects.tsv")
write_tsv(fread(path("results/gate12ab_minor_revision_closure/source_data/Figure6_neighborhood_class_effects.tsv")),
          "Figure6F_neighborhood_class_effects.tsv")

readme <- c(
  "# Gate12AD-A panel-level source and provenance package",
  "",
  "This package prepares the frozen inputs required to rebuild the six-main-figure architecture.",
  "Gate12AB inputs are read-only; no threshold, label, endpoint or numerical result is refit.",
  "",
  "## Figure 1",
  "",
  "- Panel B: all 107,886 UMAP coordinates plus deterministic median label positions.",
  "- Panel C: 13 canonical markers across 11 broad classes; color is the gene-wise z-score of the frozen average log-normalized expression and size is detected percent.",
  "- Panel D: complete 42-sample by 11-class grid with explicit zeros.",
  "- Panel E: complete sample-level fractions for four prespecified broad classes, with patient trajectory eligibility.",
  "",
  "## Figure 4",
  "",
  "The exact 49-person primary-analysis order is published. ARI excludes only seven frozen unassigned labels, leaving 42 people. The saved 10,000-value null distribution uses seed 20261107 and exactly reproduces the frozen empirical P value.",
  "",
  "## Figure 5",
  "",
  "The four patient-held-out models contain 41 observations each. All observed and predicted Axis1 values are provided together with R2, Spearman, Pearson, MAE and RMSE summaries. Favorable external scores and the failed paired endpoint remain in the same figure source package.",
  "",
  "## Figure 6",
  "",
  "All 8,190 spots are retained at measured coordinates. Full and malignant-excluded display values use one pooled absolute 98th-percentile limit. Clipping is for visualization only; inference uses unclipped scores. Histology is uniformly lightened by 12% for overlay visibility, with no local adjustment or masking.",
  "",
  "## Reproduction",
  "",
  "Run from the project root:",
  "",
  "```bash",
  "Rscript scripts/prepare_gate12ad_phase_a_source_provenance.R .",
  "python3 scripts/audit_gate12ad_phase_a.py .",
  "```"
)
writeLines(readme, file.path(out_dir, "README.md"))
writeLines(capture.output(sessionInfo()), file.path(provenance_dir, "sessionInfo.txt"))

receipt <- list(
  gate = "Gate12AD-A",
  status = "GENERATED_PENDING_INDEPENDENT_AUDIT",
  parent = "Gate12AB frozen minor-revision closure",
  figure1 = list(cells = nrow(coords), samples = nrow(sample_meta), patients = uniqueN(sample_meta$patient_id),
                 marker_genes = uniqueN(markers$gene), broad_classes = uniqueN(markers$broad_class)),
  figure4 = list(primary_patients = nrow(ordered_export), ari_patients = nrow(ari_input),
                 ari = observed_ari, permutation_p = permutation_p, permutations = permutations,
                 seed = permutation_seed),
  figure5 = list(prediction_rows = nrow(pred), models = uniqueN(pred$model), samples_per_model = 41L),
  figure6 = list(spots = nrow(overlay), sections = uniqueN(overlay$sample), shared_score_limit = shared_limit,
                 full_clipped = scale_audit[layer == "Full Axis1", n_clipped_at_shared_limit],
                 excluded_clipped = scale_audit[layer == "Malignant-excluded Axis1", n_clipped_at_shared_limit],
                 global_histology_lightening = 0.12),
  frozen_inputs_modified = FALSE
)
write_json(receipt, file.path(admin_dir, "GATE12AD_A_SOURCE_RECEIPT.json"),
           auto_unbox = TRUE, pretty = TRUE, digits = 16)

payload <- sort(list.files(out_dir, recursive = TRUE, full.names = TRUE))
payload <- payload[file.info(payload)$isdir == FALSE]
payload <- payload[!basename(payload) %chin% c(
  "GATE12AD_A_SOURCE_MANIFEST.tsv",
  "GATE12AD_A_AUDIT.tsv",
  "GATE12AD_A_AUDIT.md"
)]
output_manifest <- data.table(
  relative_path = vapply(payload, rel_path, character(1L)),
  bytes = as.numeric(file.info(payload)$size),
  sha256 = vapply(payload, sha256_file, character(1L))
)
write_tsv(output_manifest, "GATE12AD_A_SOURCE_MANIFEST.tsv", provenance_dir)

cat("GATE12AD_A_SOURCE_PREPARATION_STATUS=COMPLETE\n")
cat("OUTPUT_DIR=", rel_path(out_dir), "\n", sep = "")
cat("SOURCE_FILES=", length(list.files(source_dir, full.names = TRUE)), "\n", sep = "")
cat("OBSERVED_ARI=", format(observed_ari, digits = 16), "\n", sep = "")
cat("PERMUTATION_P=", format(permutation_p, digits = 16), "\n", sep = "")
cat("SHARED_SPATIAL_LIMIT=", format(shared_limit, digits = 16), "\n", sep = "")
