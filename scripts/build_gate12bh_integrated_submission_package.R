#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(jsonlite)
  library(magick)
  library(png)
})

options(stringsAsFactors = FALSE, warn = 1)
started_at <- Sys.time()
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript scripts/build_gate12bh_integrated_submission_package.R config/gate12bh_integrated_submission_package_v1.tsv")
}

root <- normalizePath(".", mustWork = TRUE)
config_path <- normalizePath(args[[1L]], mustWork = TRUE)
cfg_dt <- fread(config_path, sep = "\t", header = TRUE)
if (!identical(names(cfg_dt), c("key", "value")) || anyDuplicated(cfg_dt$key)) {
  stop("Configuration must contain unique key/value rows")
}
cfg <- setNames(as.list(cfg_dt$value), cfg_dt$key)
required_cfg <- c(
  "source_package", "source_repair", "output_dir", "s2_target_stem",
  "expected_be_files", "expected_main_plates", "expected_supplementary_plates",
  "expected_main_panels", "expected_supplementary_panels",
  "expected_be_checkpoint_sha256", "expected_be_figure_manifest_sha256",
  "expected_bg_receipt_sha256", "expected_bg_output_manifest_sha256",
  "expected_bg_manuscript_sha256", "expected_bg_legends_sha256"
)
if (!all(required_cfg %in% names(cfg))) stop("Configuration is missing required keys")

src_be <- normalizePath(file.path(root, cfg$source_package), mustWork = TRUE)
src_bg <- normalizePath(file.path(root, cfg$source_repair), mustWork = TRUE)
out <- file.path(root, cfg$output_dir)
out_parent <- dirname(out)
staging <- file.path(out_parent, paste0(".", basename(out), ".staging"))
if (dir.exists(out) || file.exists(out)) stop("Refusing to overwrite existing Gate12BH output: ", out)
if (dir.exists(staging) || file.exists(staging)) stop("Refusing to overwrite prior Gate12BH staging directory: ", staging)
dir.create(out_parent, recursive = TRUE, showWarnings = FALSE)

sha256_file <- function(path) digest(path, file = TRUE, algo = "sha256")
assert_hash <- function(path, expected, label) {
  observed <- sha256_file(path)
  if (!identical(observed, expected)) {
    stop(label, " input hash changed: expected ", expected, "; observed ", observed)
  }
}

be_checkpoint <- file.path(src_be, "admin", "GATE12BE_CHECKPOINT.md")
be_figure_manifest <- file.path(src_be, "provenance", "GATE12BE_FIGURE_MANIFEST.tsv")
bg_receipt <- file.path(src_bg, "GATE12BG_RECEIPT.json")
bg_output_manifest <- file.path(src_bg, "provenance", "output_manifest.tsv")
bg_manuscript <- file.path(src_bg, "manuscript", "Gate12BG_CLR_Repaired_Manuscript.md")
bg_legends <- file.path(src_bg, "manuscript", "Gate12BG_Figure_Legends.md")
assert_hash(be_checkpoint, cfg$expected_be_checkpoint_sha256, "Gate12BE checkpoint")
assert_hash(be_figure_manifest, cfg$expected_be_figure_manifest_sha256, "Gate12BE figure manifest")
assert_hash(bg_receipt, cfg$expected_bg_receipt_sha256, "Gate12BG receipt")
assert_hash(bg_output_manifest, cfg$expected_bg_output_manifest_sha256, "Gate12BG output manifest")
assert_hash(bg_manuscript, cfg$expected_bg_manuscript_sha256, "Gate12BG manuscript")
assert_hash(bg_legends, cfg$expected_bg_legends_sha256, "Gate12BG legends")

list_tree <- function(path) {
  sort(list.files(path, recursive = TRUE, full.names = TRUE, all.files = TRUE,
                  include.dirs = TRUE, no.. = TRUE))
}
source_files <- function(path) {
  x <- list_tree(path)
  x[file.exists(x) & !file.info(x)$isdir]
}
relative_to <- function(paths, base) substring(paths, nchar(base) + 2L)
copy_file <- function(src, dst, overwrite = FALSE) {
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(src, dst, overwrite = overwrite, copy.mode = TRUE, copy.date = TRUE)
  if (!isTRUE(ok)) stop("Failed to copy ", src, " -> ", dst)
}
copy_tree <- function(src, dst) {
  dir.create(dst, recursive = TRUE, showWarnings = FALSE)
  entries <- list_tree(src)
  for (p in entries) {
    rel <- relative_to(p, src)
    target <- file.path(dst, rel)
    if (isTRUE(file.info(p)$isdir)) {
      dir.create(target, recursive = TRUE, showWarnings = FALSE)
    } else {
      copy_file(p, target, overwrite = FALSE)
    }
  }
}
move_within_staging <- function(rel, parent_subdir = "gate12be") {
  src <- file.path(staging, rel)
  if (!file.exists(src)) stop("Expected inherited file is missing: ", rel)
  dst <- file.path(staging, "provenance", "parents", parent_subdir, rel)
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  if (!file.rename(src, dst)) stop("Failed to relocate inherited metadata: ", rel)
}

be_files <- source_files(src_be)
if (length(be_files) != as.integer(cfg$expected_be_files)) {
  stop("Gate12BE file count changed: expected ", cfg$expected_be_files, "; observed ", length(be_files))
}
copy_tree(src_be, staging)

## Relocate superseded Gate12BE package-level metadata so only Gate12BH files are active.
superseded_be <- c(
  "GATE12BE_Main_Figure_Contact_Sheet.png",
  "GATE12BE_Supplementary_Figure_Contact_Sheet.png",
  "GATE12BE_REVISION_MATRIX.md",
  "admin/GATE12BE_CHECKPOINT.md",
  "admin/GATE12BE_CROSS_REFERENCE_AUDIT.tsv",
  "admin/GATE12BE_FIGURE_CONTRACT.tsv",
  "admin/GATE12BE_MAIN_FIGURE_FINAL_VISUAL_AUDIT.md",
  "manuscript/Gate12BE_Review_Driven_Manuscript.md",
  "manuscript/Gate12BE_Figure_Legends.md",
  "provenance/GATE12BE_FIGURE_MANIFEST.tsv"
)
invisible(lapply(superseded_be, move_within_staging))

## Retain Gate12BG validation metadata as parent provenance.
bg_parent_files <- c(
  "GATE12BG_RECEIPT.json", "GATE12BG_RUN_LOG.txt",
  "admin/GATE12BG_FIGURE_CONTRACT.tsv", "admin/GATE12BG_MANUSCRIPT_PATCH_LOG.md",
  "admin/validation_checks.tsv", "provenance/input_hashes.tsv",
  "provenance/output_manifest.tsv", "provenance/sessionInfo.txt"
)
for (rel in bg_parent_files) {
  copy_file(file.path(src_bg, rel),
            file.path(staging, "provenance", "parents", "gate12bg", rel), overwrite = FALSE)
}

## Controlled active-content replacement: S2 plate, four source tables, manuscript and legends.
s2_stem <- cfg$s2_target_stem
copy_file(file.path(src_bg, "figures", "supplementary", "FigureS2_repaired_CLR_compositional_robustness.png"),
          file.path(staging, "figures", "supplementary", paste0(s2_stem, ".png")), overwrite = TRUE)
copy_file(file.path(src_bg, "figures", "supplementary", "FigureS2_repaired_CLR_compositional_robustness.pdf"),
          file.path(staging, "figures", "supplementary", paste0(s2_stem, ".pdf")), overwrite = TRUE)
old_s2_source <- file.path(staging, "source_data", "supplementary", "FigureS2_composition_robustness.tsv")
if (!file.exists(old_s2_source)) stop("Expected superseded S2 source table is missing")
if (!file.remove(old_s2_source)) stop("Failed to remove superseded S2 source table from staging")
new_s2_sources <- c(
  "FigureS2A_patient_tumor_minus_distal_CLR.tsv",
  "FigureS2B_cohort_specific_CLR_effects.tsv",
  "FigureS2C_analysis_set_pseudocount_sensitivity.tsv",
  "FigureS2_state_direction_summary.tsv"
)
for (f in new_s2_sources) {
  copy_file(file.path(src_bg, "source_data", "supplementary", f),
            file.path(staging, "source_data", "supplementary", f), overwrite = FALSE)
}
active_manuscript <- file.path(staging, "manuscript", "Gate12BH_Integrated_Manuscript.md")
active_legends <- file.path(staging, "manuscript", "Gate12BH_Figure_Legends.md")
copy_file(bg_manuscript, active_manuscript, overwrite = FALSE)
legend_lines <- readLines(bg_legends, warn = FALSE)
legend_lines[1L] <- "# Gate12BH figure legends"
writeLines(legend_lines, active_legends, useBytes = TRUE)

## Parent-input inventory freezes every source file consumed by the integration.
input_rows <- rbindlist(list(
  data.table(parent = "Gate12BE", relative_path = relative_to(be_files, src_be), absolute_path = be_files),
  data.table(parent = "Gate12BG", relative_path = relative_to(source_files(src_bg), src_bg),
             absolute_path = source_files(src_bg))
))
input_rows[, `:=`(bytes = file.info(absolute_path)$size,
                  sha256 = vapply(absolute_path, sha256_file, character(1)))]
input_rows[, absolute_path := NULL]
setorder(input_rows, parent, relative_path)
fwrite(input_rows, file.path(staging, "provenance", "GATE12BH_PARENT_INPUT_HASHES.tsv"), sep = "\t")

main_stems <- c(
  "Figure1_literature_standard_atlas", "Figure2_literature_standard_myeloid",
  "Figure3_literature_standard_tnk", "Figure4_review_driven_communication_hypotheses",
  "Figure5_review_driven_representation_stress_test", "Figure6_review_driven_spatial_organization"
)
supp_stems <- c(
  "FigureS1_cohort_qc_annotation", "FigureS2_composition_robustness",
  "FigureS3_tnk_transcription_robustness", "FigureS4_fixed_representation_diagnostics",
  "FigureS5_axis1_context_diagnostics", "FigureS6_paired_endpoint_robustness",
  "FigureS7_internal_functional_annotation", "FigureS7_internal_functional_annotation_continued",
  "FigureS8_external_heterogeneity_sensitivity", "FigureS9_spatial_map_archive",
  "FigureS10_spatial_sensitivity", "FigureS11_axis1_construction_redundancy_endpoint",
  "FigureS11_axis1_construction_redundancy_endpoint_continued",
  "FigureS12_sender_receiver_external_sensitivity",
  "FigureS12_sender_receiver_external_sensitivity_continued"
)
expected_panels <- c(5L, 5L, 5L, 5L, 6L, 5L,
                     4L, 3L, 3L, 2L, 3L, 2L, 2L, 1L, 3L, 3L, 2L, 2L, 4L, 4L, 1L)
names(expected_panels) <- c(main_stems, supp_stems)

figure_manifest <- rbindlist(lapply(c(main_stems, supp_stems), function(stem) {
  type <- if (stem %chin% main_stems) "main" else "supplementary"
  d <- file.path(staging, "figures", type)
  png_path <- file.path(d, paste0(stem, ".png"))
  pdf_path <- file.path(d, paste0(stem, ".pdf"))
  if (!file.exists(png_path) || !file.exists(pdf_path)) stop("Missing figure pair for ", stem)
  x <- readPNG(png_path, info = TRUE)
  dpi <- attr(x, "info")$dpi
  data.table(
    figure = stem, type = type, panels = expected_panels[[stem]],
    relative_png = file.path("figures", type, paste0(stem, ".png")),
    relative_pdf = file.path("figures", type, paste0(stem, ".pdf")),
    png_exists = TRUE, pdf_exists = TRUE,
    width_px = dim(x)[2L], height_px = dim(x)[1L], dpi_x = dpi[1L], dpi_y = dpi[2L],
    png_bytes = file.info(png_path)$size, pdf_bytes = file.info(pdf_path)$size,
    png_sha256 = sha256_file(png_path), pdf_sha256 = sha256_file(pdf_path)
  )
}))
fwrite(figure_manifest, file.path(staging, "provenance", "GATE12BH_FIGURE_MANIFEST.tsv"), sep = "\t")

## Compact contact sheets are rebuilt so the S2 thumbnail cannot remain stale.
make_plate <- function(path, label, geometry = "1050x780") {
  im <- image_read(path)
  im <- image_scale(im, geometry)
  im <- image_extent(im, "1050x780", gravity = "center", color = "white")
  lab <- image_blank(width = 1050, height = 58, color = "white")
  lab <- image_annotate(lab, label, gravity = "center", size = 27, font = "Arial", color = "#222222")
  image_append(c(lab, im), stack = TRUE)
}
make_sheet <- function(stems, dir, labels, ncol, out_path) {
  plates <- Map(function(stem, label) make_plate(file.path(dir, paste0(stem, ".png")), label), stems, labels)
  nrow_sheet <- ceiling(length(plates) / ncol)
  blank <- image_blank(width = 1050, height = 838, color = "white")
  while (length(plates) < nrow_sheet * ncol) plates[[length(plates) + 1L]] <- blank
  rows <- lapply(seq_len(nrow_sheet), function(i) {
    idx <- ((i - 1L) * ncol + 1L):(i * ncol)
    image_append(do.call(c, plates[idx]), stack = FALSE)
  })
  image_write(image_append(do.call(c, rows), stack = TRUE), out_path, format = "png")
}
make_sheet(main_stems, file.path(staging, "figures", "main"), paste("Figure", 1:6), 2L,
           file.path(staging, "GATE12BH_Main_Figure_Contact_Sheet.png"))
make_sheet(supp_stems, file.path(staging, "figures", "supplementary"),
           c(paste0("Figure S", 1:7), "Figure S7 continued", paste0("Figure S", 8:11),
             "Figure S11 continued", "Figure S12", "Figure S12 continued"), 3L,
           file.path(staging, "GATE12BH_Supplementary_Figure_Contact_Sheet.png"))

## Hard audit: inventory, inherited hashes, text claims, source tables and active-file hygiene.
checks <- list()
add_check <- function(domain, check, observed, expected, pass, note = "") {
  checks[[length(checks) + 1L]] <<- data.table(
    domain = domain, check = check, observed = as.character(observed),
    expected = as.character(expected), status = if (isTRUE(pass)) "PASS" else "FAIL", note = note
  )
}
add_check("parents", "Gate12BE source-file count", length(be_files), cfg$expected_be_files,
          length(be_files) == as.integer(cfg$expected_be_files))
add_check("parents", "Gate12BG receipt hash", sha256_file(bg_receipt), cfg$expected_bg_receipt_sha256,
          sha256_file(bg_receipt) == cfg$expected_bg_receipt_sha256)
add_check("figures", "all PNG/PDF pairs", nrow(figure_manifest), 21L, nrow(figure_manifest) == 21L)
add_check("figures", "main plates", sum(figure_manifest$type == "main"), cfg$expected_main_plates,
          sum(figure_manifest$type == "main") == as.integer(cfg$expected_main_plates))
add_check("figures", "supplementary plates", sum(figure_manifest$type == "supplementary"),
          cfg$expected_supplementary_plates,
          sum(figure_manifest$type == "supplementary") == as.integer(cfg$expected_supplementary_plates))
add_check("figures", "main panels", sum(figure_manifest[type == "main", panels]), cfg$expected_main_panels,
          sum(figure_manifest[type == "main", panels]) == as.integer(cfg$expected_main_panels))
add_check("figures", "supplementary panels", sum(figure_manifest[type == "supplementary", panels]),
          cfg$expected_supplementary_panels,
          sum(figure_manifest[type == "supplementary", panels]) == as.integer(cfg$expected_supplementary_panels))
add_check("figures", "main PNG resolution", min(figure_manifest[type == "main", dpi_x]), ">=450 dpi",
          all(figure_manifest[type == "main", dpi_x] >= 449))
add_check("figures", "supplementary PNG resolution", min(figure_manifest[type == "supplementary", dpi_x]), ">=400 dpi",
          all(figure_manifest[type == "supplementary", dpi_x] >= 399))
s2_row <- figure_manifest[figure == s2_stem]
add_check("S2", "repaired S2 panel count", s2_row$panels, 3L, s2_row$panels == 3L)
add_check("S2", "repaired S2 dimensions", paste(s2_row$width_px, s2_row$height_px, sep = "x"), "5400x4950",
          s2_row$width_px == 5400L && s2_row$height_px == 4950L)
add_check("S2", "repaired S2 PNG inherited exactly", s2_row$png_sha256,
          sha256_file(file.path(src_bg, "figures", "supplementary", "FigureS2_repaired_CLR_compositional_robustness.png")),
          s2_row$png_sha256 == sha256_file(file.path(src_bg, "figures", "supplementary", "FigureS2_repaired_CLR_compositional_robustness.png")))
add_check("S2", "repaired S2 PDF inherited exactly", s2_row$pdf_sha256,
          sha256_file(file.path(src_bg, "figures", "supplementary", "FigureS2_repaired_CLR_compositional_robustness.pdf")),
          s2_row$pdf_sha256 == sha256_file(file.path(src_bg, "figures", "supplementary", "FigureS2_repaired_CLR_compositional_robustness.pdf")))

be_manifest <- fread(be_figure_manifest)
inherited_new <- figure_manifest[figure != s2_stem]
inherited_old <- be_manifest[figure != s2_stem]
setkey(inherited_new, figure)
setkey(inherited_old, figure)
hash_compare <- inherited_new[inherited_old]
add_check("inheritance", "all 20 non-S2 figure PNG hashes unchanged", sum(hash_compare$png_sha256 == hash_compare$i.png_sha256),
          20L, nrow(hash_compare) == 20L && all(hash_compare$png_sha256 == hash_compare$i.png_sha256))
add_check("inheritance", "all 20 non-S2 figure PDF hashes unchanged", sum(hash_compare$pdf_sha256 == hash_compare$i.pdf_sha256),
          20L, nrow(hash_compare) == 20L && all(hash_compare$pdf_sha256 == hash_compare$i.pdf_sha256))

compare_source_tree <- function(src, dst, exclude_pattern = NULL) {
  sf <- source_files(src)
  if (!is.null(exclude_pattern)) sf <- sf[!grepl(exclude_pattern, basename(sf))]
  rel <- relative_to(sf, src)
  df <- file.path(dst, rel)
  data.table(relative_path = rel, source_exists = file.exists(sf), destination_exists = file.exists(df),
             hash_match = file.exists(df) & vapply(seq_along(sf), function(i) sha256_file(sf[i]) == sha256_file(df[i]), logical(1)))
}
main_cmp <- compare_source_tree(file.path(src_be, "source_data", "main"),
                                file.path(staging, "source_data", "main"))
supp_cmp <- compare_source_tree(file.path(src_be, "source_data", "supplementary"),
                                file.path(staging, "source_data", "supplementary"),
                                "^FigureS2_composition_robustness\\.tsv$")
add_check("inheritance", "all main source-data hashes unchanged", sum(main_cmp$hash_match), nrow(main_cmp), all(main_cmp$hash_match))
add_check("inheritance", "all non-S2 supplementary source-data hashes unchanged", sum(supp_cmp$hash_match), nrow(supp_cmp),
          all(supp_cmp$hash_match))
new_s2_match <- vapply(new_s2_sources, function(f) {
  sha256_file(file.path(src_bg, "source_data", "supplementary", f)) ==
    sha256_file(file.path(staging, "source_data", "supplementary", f))
}, logical(1))
add_check("S2", "four repaired source tables inherited exactly", sum(new_s2_match), 4L, all(new_s2_match))
add_check("S2", "obsolete two-panel source table absent", file.exists(old_s2_source), FALSE, !file.exists(old_s2_source))

s2a <- fread(file.path(staging, "source_data", "supplementary", new_s2_sources[1L]))
s2b <- fread(file.path(staging, "source_data", "supplementary", new_s2_sources[2L]))
s2c <- fread(file.path(staging, "source_data", "supplementary", new_s2_sources[3L]))
s2s <- fread(file.path(staging, "source_data", "supplementary", new_s2_sources[4L]))
add_check("S2", "panel A biological units", paste0(nrow(s2a), " rows; ", uniqueN(s2a$patient_key), " patients"),
          "176 rows; 11 patients", nrow(s2a) == 176L && uniqueN(s2a$patient_key) == 11L)
add_check("S2", "panel B cancer-specific effects", paste0(nrow(s2b), " rows; ", uniqueN(s2b$cancer), " cancers"),
          "32 rows; 2 cancers", nrow(s2b) == 32L && uniqueN(s2b$cancer) == 2L)
add_check("S2", "panel B exact renal resolution", min(s2b[cancer == "renal", p_value_exact]), 0.125,
          isTRUE(all.equal(min(s2b[cancer == "renal", p_value_exact]), 0.125)))
add_check("S2", "panel C frozen settings", paste0(nrow(s2c), " rows; ", uniqueN(paste(s2c$cancer, s2c$analysis_set, s2c$pseudocount))),
          "192 rows; 12 settings", nrow(s2c) == 192L && uniqueN(paste(s2c$cancer, s2c$analysis_set, s2c$pseudocount)) == 12L)
add_check("S2", "direction summary", paste0(sum(s2s$matches == 12L), ";", s2s[state == "CD4_memory", matches], ";",
                                                       s2s[state == "Proliferating_T_NK", matches]),
          "14;6;11", sum(s2s$matches == 12L) == 14L && s2s[state == "CD4_memory", matches] == 6L &&
            s2s[state == "Proliferating_T_NK", matches] == 11L)

main_src_n <- length(source_files(file.path(staging, "source_data", "main")))
supp_src_n <- length(source_files(file.path(staging, "source_data", "supplementary")))
add_check("source-data", "main source-data files", main_src_n, 42L, main_src_n == 42L)
add_check("source-data", "supplementary source-data files", supp_src_n, 43L, supp_src_n == 43L)

man_lines <- readLines(active_manuscript, warn = FALSE)
leg_lines <- readLines(active_legends, warn = FALSE)
man_text <- paste(man_lines, collapse = "\n")
leg_text <- paste(leg_lines, collapse = "\n")
add_check("manuscript", "active manuscript inherited from Gate12BG", sha256_file(active_manuscript),
          cfg$expected_bg_manuscript_sha256, sha256_file(active_manuscript) == cfg$expected_bg_manuscript_sha256)
add_check("manuscript", "standalone legend relabelled Gate12BH", leg_lines[1L], "# Gate12BH figure legends",
          identical(leg_lines[1L], "# Gate12BH figure legends"))
required_text <- c(
  "seven prostate and four renal patients", "9,999 within-cancer patient resamples",
  "Across all 16 coordinates, 14 retained the original direction in all 12",
  "smallest attainable two-sided exact sign-flip *P* value was 0.125",
  "No pooled pan-cancer significance test was used",
  "## Supplementary Figure S2. Patient-level CLR compositional robustness"
)
for (s in required_text) {
  add_check("claims", paste0("required repaired text: ", s), grepl(s, man_text, fixed = TRUE), TRUE,
            grepl(s, man_text, fixed = TRUE))
}
stale_text <- c(
  "all 16 modeled state directions matched", "all 16 directions were retained in all six",
  "Uncertainty used 999 patient-cluster bootstrap resamples", "I-squared of 91.2%",
  "## Supplementary Figure S2. Paired compositional robustness"
)
for (s in stale_text) {
  add_check("claims", paste0("obsolete text absent: ", s), grepl(s, man_text, fixed = TRUE), FALSE,
            !grepl(s, man_text, fixed = TRUE))
}
for (i in 1:6) {
  heading <- paste0("## Figure ", i, ".")
  add_check("legends", paste0("Figure ", i, " heading once in manuscript"), sum(startsWith(man_lines, heading)), 1L,
            sum(startsWith(man_lines, heading)) == 1L)
  add_check("legends", paste0("Figure ", i, " heading once in standalone legends"), sum(startsWith(leg_lines, heading)), 1L,
            sum(startsWith(leg_lines, heading)) == 1L)
}
for (i in 1:12) {
  heading <- paste0("## Supplementary Figure S", i, ".")
  add_check("legends", paste0("Figure S", i, " heading once in manuscript"), sum(startsWith(man_lines, heading)), 1L,
            sum(startsWith(man_lines, heading)) == 1L)
  add_check("legends", paste0("Figure S", i, " heading once in standalone legends"), sum(startsWith(leg_lines, heading)), 1L,
            sum(startsWith(leg_lines, heading)) == 1L)
}
ref_numbers <- as.integer(sub("^([0-9]+)\\..*", "\\1", grep("^[0-9]+\\.", man_lines, value = TRUE)))
add_check("citations", "numbered references are consecutive", paste0(length(ref_numbers), " references"), "42 consecutive",
          identical(ref_numbers, seq_len(42L)))
placeholder_n <- sum(grepl("TO BE COMPLETED", man_lines))
add_check("metadata", "author-side placeholders intentionally retained", placeholder_n, 5L, placeholder_n == 5L,
          "Author, affiliation, correspondence, funding and CRediT fields remain intentionally unresolved.")
active_manuscript_files <- list.files(file.path(staging, "manuscript"), full.names = FALSE)
add_check("hygiene", "only Gate12BH manuscript files active", paste(sort(active_manuscript_files), collapse = ";"),
          "Gate12BH_Figure_Legends.md;Gate12BH_Integrated_Manuscript.md",
          identical(sort(active_manuscript_files), sort(c("Gate12BH_Figure_Legends.md", "Gate12BH_Integrated_Manuscript.md"))))
add_check("visual-audit", "main contact sheet rebuilt", file.exists(file.path(staging, "GATE12BH_Main_Figure_Contact_Sheet.png")),
          TRUE, file.exists(file.path(staging, "GATE12BH_Main_Figure_Contact_Sheet.png")))
add_check("visual-audit", "supplementary contact sheet rebuilt", file.exists(file.path(staging, "GATE12BH_Supplementary_Figure_Contact_Sheet.png")),
          TRUE, file.exists(file.path(staging, "GATE12BH_Supplementary_Figure_Contact_Sheet.png")))

audit <- rbindlist(checks, fill = TRUE)
dir.create(file.path(staging, "admin"), recursive = TRUE, showWarnings = FALSE)
fwrite(audit, file.path(staging, "admin", "GATE12BH_CROSS_REFERENCE_AUDIT.tsv"), sep = "\t")
failed <- audit[status == "FAIL"]
checkpoint <- c(
  "# Gate12BH integrated submission-package checkpoint", "",
  paste0("- Audit status: **", if (nrow(failed)) "FAIL" else "PASS", "**"),
  paste0("- Machine checks: ", nrow(audit), " total; ", nrow(audit) - nrow(failed), " passed; ", nrow(failed), " failed."),
  "- Integration scope: Gate12BE package plus the Gate12BG v4 patient-level CLR repair; no other analysis or figure was changed.",
  "- Figure inventory: 6 main figures (31 panels) and 12 supplementary figures across 15 plates (39 panels).",
  "- S2 replacement: three panels, 11 complete-triplet patients, cancer-specific exact sign-flip inference, 9,999 patient bootstraps and 12 frozen sensitivity settings.",
  "- Inheritance: all 20 non-S2 figure pairs and every non-S2 source-data file retain the Gate12BE hashes.",
  "- Active text: the Gate12BG repaired manuscript and S2 legend are the only active manuscript files, relabelled as Gate12BH.",
  "- Claim boundary: no pooled pan-cancer CLR significance test; CLR effects remain within-lineage relative quantities.",
  "- Human-review aids: both contact sheets were rebuilt, including the repaired S2 thumbnail.",
  "- Outstanding by author instruction: author list, affiliations, correspondence, funding and CRediT remain placeholders.",
  "- Scope limitation: references were preserved from Gate12BE/Gate12BG; this gate checks numbering and cross-references but does not re-search the literature.",
  "", "## Failed checks", "",
  if (!nrow(failed)) "None." else paste0("- ", failed$domain, " / ", failed$check, ": observed ", failed$observed,
                                          "; expected ", failed$expected)
)
writeLines(checkpoint, file.path(staging, "admin", "GATE12BH_CHECKPOINT.md"), useBytes = TRUE)

writeLines(capture.output(sessionInfo()), file.path(staging, "provenance", "GATE12BH_SESSION_INFO.txt"), useBytes = TRUE)
run_log <- c(
  "Gate12BH integrated submission package",
  paste0("Started: ", format(started_at, "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Ended: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "Command: Rscript scripts/build_gate12bh_integrated_submission_package.R config/gate12bh_integrated_submission_package_v1.tsv",
  paste0("Checks: ", nrow(audit), " total; ", nrow(failed), " failed"),
  paste0("Status: ", if (nrow(failed)) "FAILED" else "COMPLETED")
)
writeLines(run_log, file.path(staging, "GATE12BH_RUN_LOG.txt"), useBytes = TRUE)

if (nrow(failed)) {
  cat("GATE12BH_AUDIT_STATUS=FAIL\n")
  cat("STAGING=", staging, "\n", sep = "")
  cat("CHECKS=", nrow(audit), " PASSED=", nrow(audit) - nrow(failed), " FAILED=", nrow(failed), "\n", sep = "")
  quit(status = 1L)
}

## Payload manifest intentionally excludes itself and the receipt; the receipt records its hash.
payload_files <- source_files(staging)
payload_files <- payload_files[!basename(payload_files) %chin% c("GATE12BH_OUTPUT_MANIFEST.tsv", "GATE12BH_RECEIPT.json")]
payload_manifest <- data.table(
  file = relative_to(payload_files, staging),
  bytes = file.info(payload_files)$size,
  sha256 = vapply(payload_files, sha256_file, character(1))
)
setorder(payload_manifest, file)
manifest_path <- file.path(staging, "provenance", "GATE12BH_OUTPUT_MANIFEST.tsv")
fwrite(payload_manifest, manifest_path, sep = "\t")

receipt <- list(
  schema_version = "1.0",
  gate_id = "Gate12BH",
  run_id = "gate12bh_run_v1",
  status = "COMPLETED",
  command = "Rscript scripts/build_gate12bh_integrated_submission_package.R config/gate12bh_integrated_submission_package_v1.tsv",
  working_directory = root,
  started_at = format(started_at, "%Y-%m-%d %H:%M:%S %Z"),
  ended_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  exit_code = 0L,
  parents = list(
    Gate12BE = list(checkpoint_sha256 = cfg$expected_be_checkpoint_sha256,
                    figure_manifest_sha256 = cfg$expected_be_figure_manifest_sha256,
                    files = length(be_files)),
    Gate12BG = list(receipt_sha256 = cfg$expected_bg_receipt_sha256,
                    output_manifest_sha256 = cfg$expected_bg_output_manifest_sha256)
  ),
  integration_scope = list(
    replaced_figure = "Supplementary Figure S2 only",
    replaced_source_tables = new_s2_sources,
    active_manuscript = "manuscript/Gate12BH_Integrated_Manuscript.md",
    active_legends = "manuscript/Gate12BH_Figure_Legends.md",
    other_figures_or_analyses_changed = FALSE
  ),
  inventory = list(main_figures = 6L, main_panels = 31L, supplementary_figures = 12L,
                   supplementary_plates = 15L, supplementary_panels = 39L,
                   payload_files = nrow(payload_manifest)),
  audit = list(total = nrow(audit), passed = nrow(audit), failed = 0L),
  payload_manifest = list(path = "provenance/GATE12BH_OUTPUT_MANIFEST.tsv",
                          sha256 = sha256_file(manifest_path))
)
write_json(receipt, file.path(staging, "GATE12BH_RECEIPT.json"), auto_unbox = TRUE, pretty = TRUE)

if (!file.rename(staging, out)) stop("Failed to finalize Gate12BH staging directory")
cat("GATE12BH_AUDIT_STATUS=PASS\n")
cat("OUTPUT=", out, "\n", sep = "")
cat("CHECKS=", nrow(audit), " PASSED=", nrow(audit), " FAILED=0\n", sep = "")
cat("PAYLOAD_FILES=", nrow(payload_manifest), "\n", sep = "")
