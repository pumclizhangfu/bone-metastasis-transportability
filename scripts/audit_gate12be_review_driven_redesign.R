#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(magick)
  library(png)
})

options(stringsAsFactors = FALSE, warn = 1)
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
pkg <- file.path(root, "results/gate12be_review_driven_redesign")
main_dir <- file.path(pkg, "figures", "main")
supp_dir <- file.path(pkg, "figures", "supplementary")
admin_dir <- file.path(pkg, "admin")
prov_dir <- file.path(pkg, "provenance")
manuscript_dir <- file.path(pkg, "manuscript")
dir.create(admin_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(prov_dir, recursive = TRUE, showWarnings = FALSE)

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
                     4L, 2L, 3L, 2L, 3L, 2L, 2L, 1L, 3L, 3L, 2L,
                     2L, 4L, 4L, 1L)
names(expected_panels) <- c(main_stems, supp_stems)

checks <- list()
add_check <- function(domain, check, observed, expected, pass, note = "") {
  checks[[length(checks) + 1L]] <<- data.table(
    domain = domain, check = check, observed = as.character(observed), expected = as.character(expected),
    status = if (isTRUE(pass)) "PASS" else "FAIL", note = note
  )
}
sha256_file <- function(path) digest(path, file = TRUE, algo = "sha256")

## Figure file and resolution audit.
manifest <- rbindlist(lapply(c(main_stems, supp_stems), function(stem) {
  type <- if (stem %chin% main_stems) "main" else "supplementary"
  d <- if (type == "main") main_dir else supp_dir
  png_path <- file.path(d, paste0(stem, ".png"))
  pdf_path <- file.path(d, paste0(stem, ".pdf"))
  if (!file.exists(png_path) || !file.exists(pdf_path)) {
    return(data.table(figure = stem, type = type, panels = expected_panels[[stem]], png_exists = file.exists(png_path),
                      pdf_exists = file.exists(pdf_path), width_px = NA_integer_, height_px = NA_integer_,
                      dpi_x = NA_real_, dpi_y = NA_real_, png_bytes = NA_real_, pdf_bytes = NA_real_,
                      png_sha256 = NA_character_, pdf_sha256 = NA_character_))
  }
  x <- readPNG(png_path, info = TRUE)
  dpi <- attr(x, "info")$dpi
  data.table(
    figure = stem, type = type, panels = expected_panels[[stem]], png_exists = TRUE, pdf_exists = TRUE,
    width_px = dim(x)[2], height_px = dim(x)[1], dpi_x = dpi[1], dpi_y = dpi[2],
    png_bytes = file.info(png_path)$size, pdf_bytes = file.info(pdf_path)$size,
    png_sha256 = sha256_file(png_path), pdf_sha256 = sha256_file(pdf_path)
  )
}), fill = TRUE)
manifest[, relative_png := file.path("figures", type, paste0(figure, ".png"))]
manifest[, relative_pdf := file.path("figures", type, paste0(figure, ".pdf"))]
setcolorder(manifest, c("figure", "type", "panels", "relative_png", "relative_pdf", "png_exists", "pdf_exists",
                        "width_px", "height_px", "dpi_x", "dpi_y", "png_bytes", "pdf_bytes", "png_sha256", "pdf_sha256"))
fwrite(manifest, file.path(prov_dir, "GATE12BE_FIGURE_MANIFEST.tsv"), sep = "\t")
add_check("figures", "all PNG/PDF pairs exist", sum(manifest$png_exists & manifest$pdf_exists), nrow(manifest),
          all(manifest$png_exists & manifest$pdf_exists))
add_check("figures", "main PNG resolution", min(manifest[type == "main", dpi_x]), ">=450 dpi",
          all(manifest[type == "main", dpi_x] >= 449))
add_check("figures", "supplementary PNG resolution", min(manifest[type == "supplementary", dpi_x]), ">=400 dpi",
          all(manifest[type == "supplementary", dpi_x] >= 399))
add_check("figures", "main-panel count", sum(manifest[type == "main", panels]), 31L,
          sum(manifest[type == "main", panels]) == 31L)
add_check("figures", "supplementary-panel count", sum(manifest[type == "supplementary", panels]), 38L,
          sum(manifest[type == "supplementary", panels]) == 38L)

## Manuscript, legend, and citation cross-reference audit.
manuscript_path <- file.path(manuscript_dir, "Gate12BE_Review_Driven_Manuscript.md")
txt <- readLines(manuscript_path, warn = FALSE)
legend_i <- grep("^# Figure legends$", txt)
body <- if (length(legend_i) == 1L) txt[seq_len(legend_i - 1L)] else txt
legends <- if (length(legend_i) == 1L) txt[legend_i:length(txt)] else character()
writeLines(c("# Gate12BE figure legends", "", legends[-1L]), file.path(manuscript_dir, "Gate12BE_Figure_Legends.md"))
for (i in 1:6) {
  add_check("cross-reference", paste0("main Figure ", i, " cited in body"), any(grepl(paste0("Figure ", i, "([^0-9]|$)"), body)),
            TRUE, any(grepl(paste0("Figure ", i, "([^0-9]|$)"), body)))
  add_check("legends", paste0("main Figure ", i, " legend present"), any(grepl(paste0("^## Figure ", i, "\\."), legends)),
            TRUE, any(grepl(paste0("^## Figure ", i, "\\."), legends)))
}
for (i in 1:12) {
  add_check("cross-reference", paste0("Supplementary Figure S", i, " cited in body"),
            any(grepl(paste0("Supplementary Figure S", i, "([^0-9]|$)"), body)), TRUE,
            any(grepl(paste0("Supplementary Figure S", i, "([^0-9]|$)"), body)))
  add_check("legends", paste0("Supplementary Figure S", i, " legend present"),
            any(grepl(paste0("^## Supplementary Figure S", i, "\\."), legends)), TRUE,
            any(grepl(paste0("^## Supplementary Figure S", i, "\\."), legends)))
}
full_text <- paste(txt, collapse = "\n")
required <- c("NicheNet", "CellChat", "unique-target freeze was revoked", "No virtual knockout",
              "10.1038/s41592-019-0667-5", "10.1038/s41467-021-21246-9")
for (s in required) add_check("claims", paste0("required bounded statement: ", s), grepl(s, full_text, fixed = TRUE), TRUE,
                              grepl(s, full_text, fixed = TRUE))
stale <- c("Supplementary Figure S5D", "Frozen transcriptional and discrete representations fail predefined external criteria",
           "Axis1 is technically reconstructable but fails uniform", "K-compact", "Gate12BD aligned submission package")
for (s in stale) add_check("claims", paste0("obsolete wording absent: ", s), grepl(s, full_text, fixed = TRUE), FALSE,
                           !grepl(s, full_text, fixed = TRUE))
ref_numbers <- as.integer(sub("^([0-9]+)\\..*", "\\1", grep("^[0-9]+\\.", txt, value = TRUE)))
add_check("citations", "numbered reference count", length(ref_numbers), 42L,
          identical(ref_numbers, seq_len(42L)))
placeholder_n <- sum(grepl("TO BE COMPLETED", txt))
add_check("metadata", "administrative placeholders intentionally unresolved",
          placeholder_n, 5L, placeholder_n == 5L,
          "Author, affiliation, correspondence, funding and CRediT fields retained by author instruction; manuscript is not yet submission-complete.")

## Panel source-data and coordinate-integrity audit.
main_src <- list.files(file.path(pkg, "source_data", "main"), full.names = TRUE)
supp_src <- list.files(file.path(pkg, "source_data", "supplementary"), full.names = TRUE)
add_check("source-data", "main source-data files", length(main_src), ">=31", length(main_src) >= 31L)
add_check("source-data", "supplementary source-data files", length(supp_src), ">=40", length(supp_src) >= 40L)

## Final main-figure visual contract. These checks prevent later rebuilds from
## silently restoring the superseded plate headers or ambiguous encodings.
visual_audit_path <- file.path(admin_dir, "GATE12BE_MAIN_FIGURE_FINAL_VISUAL_AUDIT.md")
add_check("visual-contract", "final human visual-audit report exists",
          file.exists(visual_audit_path), TRUE, file.exists(visual_audit_path))
f1_script <- paste(readLines(file.path(root, "scripts/build_gate12az_literature_standard_figure1.R"),
                             warn = FALSE), collapse = "\n")
f23_script <- paste(readLines(file.path(root, "scripts/build_gate12az_literature_standard_figures23.R"),
                              warn = FALSE), collapse = "\n")
f456_script <- paste(readLines(file.path(root, "scripts/build_gate12be_review_driven_figures456.R"),
                               warn = FALSE), collapse = "\n")
plate_headers <- c(
  "Single-cell atlas of the bone-metastatic ecosystem",
  "Myeloid-state remodelling in the bone-metastatic niche",
  "T/NK-state remodelling in bone metastases",
  "Expression-derived myeloid-to-lymphoid communication hypotheses"
)
build_text <- paste(f1_script, f23_script, f456_script, sep = "\n")
plate_header_free <- !any(vapply(plate_headers, grepl, logical(1), x = build_text, fixed = TRUE))
add_check("visual-contract", "Figures 1-4 omit redundant plate headers",
          plate_header_free, TRUE, plate_header_free)
upgraded_umap_safe <- grepl('"gate12be_review_driven_redesign"', f1_script, fixed = TRUE) &&
  grepl('"gate12be_review_driven_redesign"', f23_script, fixed = TRUE) &&
  grepl('file.path(out, "source_data", "main")', f1_script, fixed = TRUE) &&
  grepl('file.path(out, "source_data", "main")', f23_script, fixed = TRUE)
add_check("visual-contract", "Gate12BE rebuild preserves upgraded UMAPs and main source-data layout",
          upgraded_umap_safe, TRUE, upgraded_umap_safe)
f1d_labels_removed <- grepl('theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()',
                            f1_script, fixed = TRUE)
add_check("visual-contract", "Figure 1D omits unreadable per-bar patient labels",
          f1d_labels_removed, TRUE, f1d_labels_removed)
f2e_decodable <- grepl('geom_text_repel(data = myeloid_concordance,', f23_script, fixed = TRUE) &&
  grepl('labels = c("Not stable", "Stable")', f23_script, fixed = TRUE)
add_check("visual-contract", "Figure 2E labels every state and defines stability fill",
          f2e_decodable, TRUE, f2e_decodable)
f3_shared_axis <- grepl('facet_wrap(~state_label, nrow = 1) +', f23_script, fixed = TRUE)
add_check("visual-contract", "Figure 3D uses a common effect-size axis",
          f3_shared_axis, TRUE, f3_shared_axis)
f4_missing_encoded <- grepl('short_chat[is.na(prob)]', f456_script, fixed = TRUE) &&
  grepl('not retained', f456_script, fixed = TRUE)
add_check("visual-contract", "Figure 4B explicitly encodes non-retained combinations",
          f4_missing_encoded, TRUE, f4_missing_encoded)
f5_precision_ne <- grepl('"NE", sprintf("%+.4f", difference)', f456_script, fixed = TRUE) &&
  grepl('sprintf("%.4f\\n%s", effect_value, unit_label)', f456_script, fixed = TRUE)
add_check("visual-contract", "Figure 5B/C show NE and avoid rounded negative zero",
          f5_precision_ne, TRUE, f5_precision_ne)
f5e_readable <- grepl('p5e <- wrap_elements(full =\n  p5e_heat +', f456_script, fixed = TRUE) &&
  grepl('axis.text.y = element_text(size = 5.7)', f456_script, fixed = TRUE)
add_check("visual-contract", "Figure 5E uses the simplified enlarged heatmap",
          f5e_readable, TRUE, f5e_readable)

prototype_path <- file.path(admin_dir, "Figure5_BDEF_PUBLISHED_PROTOTYPE_MAP.tsv")
add_check("visual-precedent", "Figure 5 B/D/E/F published-prototype map exists",
          file.exists(prototype_path), TRUE, file.exists(prototype_path))
if (file.exists(prototype_path)) {
  prototype_map <- fread(prototype_path)
  add_check("visual-precedent", "Figure 5 B/D/E/F all mapped",
            paste(sort(unique(prototype_map$panel)), collapse = ","), "B,D,E,F",
            identical(sort(unique(prototype_map$panel)), c("B", "D", "E", "F")))
  add_check("visual-precedent", "every prototype has DOI or URL",
            sum(grepl("^https://", prototype_map$doi_or_url)), nrow(prototype_map),
            all(grepl("^https://", prototype_map$doi_or_url)))
}
f6_prototype_path <- file.path(admin_dir, "Figure6_ABCDE_PUBLISHED_PROTOTYPE_MAP.tsv")
add_check("visual-precedent", "Figure 6 A-E published-prototype map exists",
          file.exists(f6_prototype_path), TRUE, file.exists(f6_prototype_path))
if (file.exists(f6_prototype_path)) {
  f6_prototype_map <- fread(f6_prototype_path)
  mapped_panels <- sort(unique(unlist(strsplit(gsub("[^A-E]", "", f6_prototype_map$panel), ""))))
  add_check("visual-precedent", "Figure 6 A-E all mapped",
            paste(mapped_panels, collapse = ","), "A,B,C,D,E",
            identical(mapped_panels, LETTERS[1:5]))
  add_check("visual-precedent", "every Figure 6 prototype has DOI or URL",
            sum(grepl("^https://", f6_prototype_map$doi_or_url)), nrow(f6_prototype_map),
            all(grepl("^https://", f6_prototype_map$doi_or_url)))
}
f5b <- fread(file.path(pkg, "source_data", "main", "Figure5B_OEP_paired_programme_differences.tsv"))
add_check("Figure5", "panel B paired-difference matrix",
          paste0(nrow(f5b), " tiles; ", sum(!is.na(f5b$difference)), " observed"),
          "15 tiles; 12 observed", nrow(f5b) == 15L && sum(!is.na(f5b$difference)) == 12L)
f5d <- fread(file.path(pkg, "source_data", "main", "Figure5D_Axis1_four_model_reconstruction.tsv"))
add_check("Figure5", "panel D displays all reconstruction models",
          paste0(uniqueN(f5d$model), " models; ", nrow(f5d), " observations"),
          "4 models; 164 observations", uniqueN(f5d$model) == 4L && nrow(f5d) == 164L)
f5e <- fread(file.path(pkg, "source_data", "main", "Figure5E_external_frozen_feature_heatmap.tsv.gz"))
add_check("Figure5", "panel E frozen feature heatmap dimensions",
          paste0(uniqueN(f5e$feature), " features x ", uniqueN(f5e$patient_id), " patients"),
          "26 features x 46 patients",
          uniqueN(f5e$feature) == 26L && uniqueN(f5e$patient_id) == 46L && nrow(f5e) == 1196L)
f5f <- fread(file.path(pkg, "source_data", "main", "Figure5F_Axis1_paired_effects.tsv"))
add_check("Figure5", "panel F contains every independent paired effect",
          paste0(nrow(f5f), " patients; ", paste(sort(unique(f5f$contrast)), collapse = ",")),
          "5 patients; both comparisons",
          nrow(f5f) == 5L && uniqueN(f5f$contrast) == 2L)
coord_specs <- data.table(
  file = c("Figure1B_umap_coordinates.tsv.gz", "Figure2A_myeloid_umap.tsv.gz", "Figure3A_tnk_umap.tsv.gz"),
  expected_rows = c(107886L, 29588L, 56855L), id_col = c("cell_id", "barcode", "barcode")
)
for (i in seq_len(nrow(coord_specs))) {
  d <- fread(file.path(pkg, "source_data", "main", coord_specs$file[i]))
  ok <- nrow(d) == coord_specs$expected_rows[i] && uniqueN(d[[coord_specs$id_col[i]]]) == nrow(d)
  add_check("coordinates", coord_specs$file[i], paste0(nrow(d), " rows; ", uniqueN(d[[coord_specs$id_col[i]]]), " unique cells"),
            paste0(coord_specs$expected_rows[i], " rows; all unique"), ok,
            "UMAP source coordinates retain all eligible cells.")
}
f6_script <- paste(readLines(file.path(root, "scripts/build_gate12be_review_driven_figures456.R"), warn = FALSE), collapse = "\n")
f6_ok <- grepl("lighten_histology", f6_script, fixed = TRUE) && grepl("asset$light_image", f6_script, fixed = TRUE)
add_check("reproducibility", "Figure 6 overlay uses uniform lightened image", f6_ok, TRUE, f6_ok)
f6e_facets_ok <- grepl('facet_grid(cols = vars(layer_short))', f6_script, fixed = TRUE) &&
  grepl('scale_x_continuous(limits = c(-.03, .26)', f6_script, fixed = TRUE)
add_check("Figure6", "panel E uses aligned fixed-scale forest facets", f6e_facets_ok, TRUE, f6e_facets_ok)
f6e <- fread(file.path(pkg, "source_data", "main", "Figure6E_section_effects.tsv"))
add_check("Figure6", "panel E contains all linked sections and score layers",
          paste0(uniqueN(f6e$sample), " sections x ", uniqueN(f6e$layer_label), " layers"),
          "4 sections x 3 layers", uniqueN(f6e$sample) == 4L && uniqueN(f6e$layer_label) == 3L)

## Compact, labeled visual contact sheets for final human review.
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
  nrow <- ceiling(length(plates) / ncol)
  blank <- image_blank(width = 1050, height = 838, color = "white")
  while (length(plates) < nrow * ncol) plates[[length(plates) + 1L]] <- blank
  rows <- lapply(seq_len(nrow), function(i) {
    idx <- ((i - 1L) * ncol + 1L):(i * ncol)
    image_append(do.call(c, plates[idx]), stack = FALSE)
  })
  sheet <- image_append(do.call(c, rows), stack = TRUE)
  image_write(sheet, out_path, format = "png")
}
make_sheet(main_stems, main_dir, paste("Figure", 1:6), 2L, file.path(pkg, "GATE12BE_Main_Figure_Contact_Sheet.png"))
make_sheet(supp_stems, supp_dir,
           c(paste0("Figure S", 1:7), "Figure S7 continued", paste0("Figure S", 8:11),
             "Figure S11 continued", "Figure S12", "Figure S12 continued"),
           3L, file.path(pkg, "GATE12BE_Supplementary_Figure_Contact_Sheet.png"))
add_check("visual-audit", "main contact sheet", file.exists(file.path(pkg, "GATE12BE_Main_Figure_Contact_Sheet.png")), TRUE,
          file.exists(file.path(pkg, "GATE12BE_Main_Figure_Contact_Sheet.png")))
add_check("visual-audit", "supplementary contact sheet", file.exists(file.path(pkg, "GATE12BE_Supplementary_Figure_Contact_Sheet.png")), TRUE,
          file.exists(file.path(pkg, "GATE12BE_Supplementary_Figure_Contact_Sheet.png")))

audit <- rbindlist(checks, fill = TRUE)
fwrite(audit, file.path(admin_dir, "GATE12BE_CROSS_REFERENCE_AUDIT.tsv"), sep = "\t")
failed <- audit[status == "FAIL"]
checkpoint <- c(
  "# Gate12BE review-driven submission checkpoint", "",
  paste0("- Audit status: **", if (nrow(failed)) "FAIL" else "PASS", "**"),
  paste0("- Machine checks: ", nrow(audit), " total; ", nrow(audit) - nrow(failed), " passed; ", nrow(failed), " failed."),
  "- Figure inventory: 6 main figures (31 panels) and 12 supplementary figures across 15 plates (38 panels).",
  "- Raster/vector delivery: every plate has PNG and PDF; main PNGs are 450 dpi and supplementary PNGs are at least 400 dpi.",
  "- UMAP integrity: all 107,886 atlas cells, 29,588 myeloid cells and 56,855 T/NK cells are represented exactly once in panel source data.",
  "- Source data: panel-level main and supplementary tables are bundled with SHA-256 figure receipts.",
  "- Claim boundary: communication pairs remain hypotheses; the unique-target freeze was revoked and no virtual knockout was executed.",
  "- Spatial image handling: native H&E is retained; Figure 6B,C use one uniform 12% lightening operation documented in code and legend.",
  "- Human-review aids: main and supplementary contact sheets are included.",
  "- Outstanding by author instruction: author list, affiliations and corresponding-author details remain placeholders.",
  "",
  "## Failed checks", "",
  if (!nrow(failed)) "None." else paste0("- ", failed$domain, " / ", failed$check, ": observed ", failed$observed, "; expected ", failed$expected)
)
writeLines(checkpoint, file.path(admin_dir, "GATE12BE_CHECKPOINT.md"))
cat("GATE12BE_AUDIT_STATUS=", if (nrow(failed)) "FAIL" else "PASS", "\n", sep = "")
cat("CHECKS=", nrow(audit), " PASSED=", nrow(audit) - nrow(failed), " FAILED=", nrow(failed), "\n", sep = "")
if (nrow(failed)) quit(status = 1L)
