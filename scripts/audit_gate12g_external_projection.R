#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: audit_gate12g_external_projection.R <full_loadings.tsv> <gate9b_patient_counts.tsv> <outdir>")
}

loading_file <- args[[1L]]
count_file <- args[[2L]]
outdir <- args[[3L]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

loadings <- fread(loading_file)
counts <- fread(count_file)
required_loading <- c("feature", "axis", "loading", "block")
required_counts <- c("patient_id", "cancer_code", "primary_analysis", "all_qc_cells",
                     "assigned_cells", "assigned_fraction", "CD14HI_MONO", "CD16HI_MONO",
                     "CD4_TREG", "CD8_TEX", "MACROPHAGE", "OSTEOCLAST")
if (!all(required_loading %chin% names(loadings))) stop("Loading input lacks required columns")
if (!all(required_counts %chin% names(counts))) stop("Gate9B count input lacks required columns")

counts[, primary_analysis := as.logical(primary_analysis)]
primary <- counts[primary_analysis == TRUE]

# This table is frozen before coverage results are calculated. `proxy` means that either
# the biological definition or the denominator differs from Gate12G and cannot support
# a full frozen-coordinate projection.
crosswalk <- data.table(feature = unique(loadings$feature))
crosswalk[, `:=`(gate9b_source = NA_character_, mapping_status = "absent",
                 mapping_note = "not represented by the six-state Gate9B classifier")]

set_map <- function(feature, source, status, note) {
  feature_value <- feature
  crosswalk[feature == feature_value,
            `:=`(gate9b_source = source, mapping_status = status, mapping_note = note)]
}

set_map("Broad__T_NK", "broad_lineage=t_cell", "proxy",
        "Gate9B broad lineage excludes NK and is not the frozen Gate12A broad classifier")
set_map("Broad__Myeloid", "broad_lineage=myeloid", "proxy",
        "Gate9B eligibility gate differs from the frozen Gate12A broad classifier")
set_map("Broad__Osteoclast", "OSTEOCLAST/all_qc_cells", "exact",
        "direct osteoclast count over the all-cell denominator")
set_map("Myeloid__Classical_monocyte", "CD14HI_MONO", "proxy",
        "marker-compatible state but the complete within-myeloid denominator is unavailable")
set_map("Myeloid__Inflammatory_monocyte", "CD16HI_MONO", "proxy",
        "CD16-high and inflammatory-monocyte definitions are not equivalent")
set_map("Myeloid__C1QC_macrophage", "MACROPHAGE", "proxy",
        "Gate9B macrophages do not separate C1QC and resident states")
set_map("Myeloid__Resident_macrophage", "MACROPHAGE", "proxy",
        "Gate9B macrophages do not separate C1QC and resident states")
set_map("T_NK__Treg", "CD4_TREG", "proxy",
        "state is compatible but the complete within-T/NK denominator is unavailable")
set_map("T_NK__CD8_exhausted", "CD8_TEX", "proxy",
        "state is compatible but the complete within-T/NK denominator is unavailable")

coverage <- merge(loadings, crosswalk, by = "feature", all.x = TRUE)
coverage[, loading_sq := loading^2]
coverage_summary <- coverage[, .(
  total_loading_sq = sum(loading_sq),
  exact_loading_sq = sum(loading_sq[mapping_status == "exact"]),
  exact_or_proxy_loading_sq = sum(loading_sq[mapping_status %chin% c("exact", "proxy")])
), by = axis]
coverage_summary[, `:=`(
  exact_coverage = exact_loading_sq / total_loading_sq,
  exact_or_proxy_coverage = exact_or_proxy_loading_sq / total_loading_sq,
  required_coverage = 0.80
)]
coverage_summary[, existing_gate9b_full_projection_eligible := exact_coverage >= required_coverage]

cohort_audit <- data.table(
  metric = c("archive_rows", "primary_patients", "primary_cancer_origins",
             "primary_min_qc_cells", "primary_median_qc_cells",
             "primary_min_assigned_fraction", "primary_median_assigned_fraction",
             "primary_patients_ge_500_qc"),
  value = c(nrow(counts), nrow(primary), uniqueN(primary$cancer_code), min(primary$all_qc_cells),
            median(primary$all_qc_cells), min(primary$assigned_fraction),
            median(primary$assigned_fraction), sum(primary$all_qc_cells >= 500))
)

decision <- if (all(!coverage_summary$existing_gate9b_full_projection_eligible)) {
  "FULL_REANNOTATION_REQUIRED"
} else {
  "EXISTING_GATE9B_MAY_BE_SUFFICIENT"
}

fwrite(crosswalk, file.path(outdir, "gate9b_to_gate12g_feature_crosswalk.tsv"), sep = "\t")
fwrite(coverage, file.path(outdir, "gate9b_loading_coverage_detail.tsv"), sep = "\t")
fwrite(coverage_summary, file.path(outdir, "gate9b_loading_coverage_summary.tsv"), sep = "\t")
fwrite(cohort_audit, file.path(outdir, "oep005136_cohort_audit.tsv"), sep = "\t")

checkpoint <- c(
  "# Gate12G OEP005136 input and coverage audit",
  "",
  paste0("- Primary patients: ", nrow(primary)),
  paste0("- Cancer origins: ", uniqueN(primary$cancer_code)),
  paste0("- Patients with >=500 QC cells: ", sum(primary$all_qc_cells >= 500), "/", nrow(primary)),
  paste0("- Existing Gate9B full-projection decision: **", decision, "**"),
  "",
  "The six-state Gate9B output is retained for secondary concordance only. Proxy mappings do not count toward the frozen 80% full-projection coverage threshold.",
  "",
  "A new discovery-reference transfer is required before any OEP005136 Gate12G axis score is calculated."
)
writeLines(checkpoint, file.path(outdir, "GATE12G_EXTERNAL_INPUT_AUDIT.md"))

cat("GATE12G_EXTERNAL_INPUT_DECISION=", decision, "\n", sep = "")
