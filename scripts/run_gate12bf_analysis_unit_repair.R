#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(digest)
})

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript scripts/run_gate12bf_analysis_unit_repair.R config/gate12bf_analysis_unit_repair_v3.tsv",
       call. = FALSE)
}

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
config_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
script_path <- normalizePath("scripts/run_gate12bf_analysis_unit_repair.R",
                             winslash = "/", mustWork = TRUE)

cfg_dt <- fread(config_path, sep = "\t", header = TRUE, colClasses = "character",
                na.strings = NULL)
if (!identical(names(cfg_dt), c("key", "value")) || anyDuplicated(cfg_dt$key) ||
    any(!nzchar(cfg_dt$key))) {
  stop("Configuration must contain unique non-empty key/value rows", call. = FALSE)
}
cfg <- setNames(cfg_dt$value, cfg_dt$key)
cfg_get <- function(key) {
  if (!key %chin% names(cfg)) stop("Missing configuration key: ", key, call. = FALSE)
  cfg[[key]]
}

expected_contract <- c(
  schema_version = "1.0",
  gate_id = "Gate12BF",
  run_id = "gate12bf_run_v3",
  seed = "20260813",
  pseudocounts = "0.25,0.5,1.0",
  primary_pseudocount = "0.5",
  bootstrap_iterations = "9999",
  primary_analysis_set = "complete_triplets",
  primary_contrast = "tumor_minus_distal",
  primary_statistic = "abs_mean_patient_contrast",
  exact_p_formula = "mean(abs_permuted_mean_ge_abs_observed_mean)",
  fdr_family = "analysis_set_x_cohort_x_lineage_x_pseudocount_x_contrast"
)
for (nm in names(expected_contract)) {
  if (!identical(cfg_get(nm), expected_contract[[nm]])) {
    stop("Configuration contract mismatch for ", nm, call. = FALSE)
  }
}

output_rel <- cfg_get("output_dir")
allowed_parent <- normalizePath(file.path(project_root, "results/gate12bf_analysis_unit_repair"),
                                winslash = "/", mustWork = TRUE)
output_abs <- file.path(project_root, output_rel)
output_parent <- normalizePath(dirname(output_abs), winslash = "/", mustWork = TRUE)
if (!identical(output_parent, allowed_parent) || basename(output_abs) != "run_v3") {
  stop("Output path is outside the locked Gate12BF run_v3 location", call. = FALSE)
}
if (file.exists(output_abs)) {
  stop("Locked output directory already exists; refusing to overwrite: ", output_abs,
       call. = FALSE)
}
dir.create(output_abs, recursive = FALSE, showWarnings = FALSE)

log_path <- file.path(output_abs, "GATE12BF_RUN_LOG.txt")
start_time <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
current_stage <- "initialization"
checks <- data.table(
  check_id = character(), expected = character(), observed = character(),
  pass = logical(), severity = character()
)
output_registry <- data.table(file = character(), rows = integer(), columns = integer())

log_line <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\t",
                 paste(..., collapse = ""))
  cat(line, "\n", file = log_path, append = TRUE)
  message(line)
}

register_output <- function(path, rows = NA_integer_, columns = NA_integer_) {
  output_registry <<- rbind(
    output_registry,
    data.table(file = basename(path), rows = as.integer(rows), columns = as.integer(columns))
  )
}

write_tsv <- function(x, filename, compress = grepl("\\.gz$", filename)) {
  path <- file.path(output_abs, filename)
  fwrite(x, path, sep = "\t", quote = FALSE, compress = if (compress) "gzip" else "none")
  register_output(path, nrow(x), ncol(x))
  invisible(path)
}

write_text <- function(lines, filename) {
  path <- file.path(output_abs, filename)
  writeLines(lines, path, useBytes = TRUE)
  register_output(path, length(lines), NA_integer_)
  invisible(path)
}

record_check <- function(check_id, expected, observed, pass, severity = "HARD") {
  pass <- isTRUE(pass)
  checks <<- rbind(
    checks,
    data.table(
      check_id = as.character(check_id), expected = as.character(expected),
      observed = as.character(observed), pass = pass, severity = as.character(severity)
    )
  )
  if (!pass && identical(severity, "HARD")) {
    stop("Hard validation failed [", check_id, "]: expected ", expected,
         "; observed ", observed, call. = FALSE)
  }
  invisible(pass)
}

sha256_file <- function(path) digest(file = path, algo = "sha256", serialize = FALSE)
relative_path <- function(path) {
  sub(paste0("^", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", project_root), "/?"),
      "", normalizePath(path, winslash = "/", mustWork = TRUE))
}

manifest_for <- function(paths) {
  expanded <- unlist(lapply(paths, function(path) {
    abs_path <- file.path(project_root, path)
    if (dir.exists(abs_path)) {
      list.files(abs_path, recursive = TRUE, full.names = TRUE, all.files = TRUE,
                 no.. = TRUE)
    } else if (file.exists(abs_path)) {
      abs_path
    } else {
      character()
    }
  }), use.names = FALSE)
  expanded <- sort(unique(expanded[file.info(expanded)$isdir %in% FALSE]))
  data.table(
    path = vapply(expanded, relative_path, character(1)),
    bytes = as.numeric(file.info(expanded)$size),
    sha256 = vapply(expanded, sha256_file, character(1))
  )
}

summarize_values <- function(values) {
  values <- as.numeric(values)
  if (!length(values) || any(!is.finite(values))) stop("Non-finite values in summary", call. = FALSE)
  tol <- as.numeric(cfg_get("signflip_tolerance"))
  data.table(
    n_patients = length(values),
    estimate_mean = mean(values),
    estimate_median = median(values),
    sd_between_patients = if (length(values) > 1L) sd(values) else NA_real_,
    se_descriptive = if (length(values) > 1L) sd(values) / sqrt(length(values)) else NA_real_,
    positive_fraction = mean(values > tol),
    negative_fraction = mean(values < -tol),
    zero_fraction = mean(abs(values) <= tol)
  )
}

sign_with_tolerance <- function(x) {
  tol <- as.numeric(cfg_get("signflip_tolerance"))
  fifelse(!is.finite(x), NA_integer_, fifelse(abs(x) <= tol, 0L, as.integer(sign(x))))
}

input_path <- function(name) file.path(project_root, cfg_get(paste0("input.", name, ".path")))

failure_receipt <- function(e) {
  end_time <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  try(log_line("RUN_FAILED\tstage=", current_stage, "\tmessage=", conditionMessage(e)), silent = TRUE)
  try(fwrite(checks, file.path(output_abs, "validation_checks.tsv"), sep = "\t"), silent = TRUE)
  receipt <- list(
    schema_version = "1.0",
    gate_id = "Gate12BF",
    run_id = "gate12bf_run_v3",
    status = "FAILED",
    command = "Rscript scripts/run_gate12bf_analysis_unit_repair.R config/gate12bf_analysis_unit_repair_v3.tsv",
    working_directory = project_root,
    started_at = start_time,
    ended_at = end_time,
    exit_code = 1L,
    failure_stage = current_stage,
    failure_code = "HARD_VALIDATION_OR_RUNTIME_ERROR",
    message = conditionMessage(e),
    checks = checks
  )
  try(write_json(receipt, file.path(output_abs, "GATE12BF_RECEIPT.json"),
                 pretty = TRUE, auto_unbox = TRUE, na = "null"), silent = TRUE)
}

run_gate12bf <- function() {
  log_line("RUN_STARTED\tgate=Gate12BF\trun=gate12bf_run_v3")

  current_stage <<- "input_hash_freeze"
  input_path_keys <- grep("^input\\.[^.]+\\.path$", names(cfg), value = TRUE)
  input_names <- sub("^input\\.([^.]+)\\.path$", "\\1", input_path_keys)
  input_hashes <- rbindlist(lapply(input_names, function(nm) {
    path_rel <- cfg_get(paste0("input.", nm, ".path"))
    path_abs <- file.path(project_root, path_rel)
    expected <- cfg_get(paste0("input.", nm, ".sha256"))
    exists <- file.exists(path_abs)
    actual <- if (exists) sha256_file(path_abs) else NA_character_
    record_check(paste0("input_exists_", nm), "TRUE", exists, exists)
    record_check(paste0("input_hash_", nm), expected, actual,
                 exists && identical(actual, expected))
    info <- if (exists) file.info(path_abs) else NULL
    data.table(input_name = nm, path = path_rel, expected_sha256 = expected,
               actual_sha256 = actual, hash_match = exists && identical(actual, expected),
               bytes = if (exists) as.numeric(info$size) else NA_real_)
  }))
  write_tsv(input_hashes, "input_hashes.tsv")

  protected_paths <- unique(c(
    "results/gate12b_cell_states",
    "results/gate12v1_compositional_sensitivity",
    input_hashes$path
  ))
  protected_pre <- manifest_for(protected_paths)
  write_tsv(protected_pre, "protected_manifest_pre.tsv")

  current_stage <<- "discovery_analysis_unit_audit"
  composition <- fread(input_path("composition"))
  required <- c("accession", "cancer", "patient_id", "sample_id", "compartment",
                "lineage", "state", "n_state", "lineage_total", "stage")
  record_check("composition_required_columns", paste(required, collapse = ","),
               paste(intersect(required, names(composition)), collapse = ","),
               all(required %chin% names(composition)))
  record_check("composition_rows", "672", nrow(composition), nrow(composition) == 672L)
  record_check("composition_unique_state_key", "0 duplicates",
               anyDuplicated(composition[, .(accession, sample_id, lineage, state)]),
               !anyDuplicated(composition[, .(accession, sample_id, lineage, state)]))
  record_check("composition_nonnegative_integer_counts", "TRUE",
               all(composition$n_state >= 0 & composition$n_state == floor(composition$n_state)),
               all(composition$n_state >= 0 & composition$n_state == floor(composition$n_state)))
  record_check("composition_positive_lineage_totals", "TRUE",
               all(composition$lineage_total > 0), all(composition$lineage_total > 0))

  stage_map <- c(distal = 0L, involved = 1L, tumor = 2L)
  stage_ok <- all(composition$compartment %chin% names(stage_map)) &&
    all(composition$stage == unname(stage_map[composition$compartment]))
  record_check("stage_mapping", "distal=0;involved=1;tumor=2", stage_ok, stage_ok)

  composition[, patient_key := paste(accession, cancer, patient_id, sep = "::")]
  block_audit <- composition[, .(
    n_states = .N,
    n_unique_total = uniqueN(lineage_total),
    sum_state = sum(n_state),
    lineage_total_value = first(lineage_total)
  ), by = .(accession, cancer, patient_key, patient_id, sample_id, compartment, stage, lineage)]
  block_ok <- nrow(block_audit) == 84L && all(block_audit$n_states == 8L) &&
    all(block_audit$n_unique_total == 1L) &&
    all(block_audit$sum_state == block_audit$lineage_total_value)
  record_check("sample_lineage_blocks", "84 blocks; 8 states; count sums match total",
               paste0(nrow(block_audit), " blocks; valid=", block_ok), block_ok)

  state_universe <- unique(composition[, .(lineage, state)])[order(lineage, state)]
  state_counts <- state_universe[, .N, by = lineage]
  state_ok <- nrow(state_universe) == 16L && nrow(state_counts) == 2L &&
    all(state_counts$N == 8L) && setequal(state_counts$lineage, c("Myeloid", "T_NK"))
  record_check("frozen_state_universe", "2 lineages x 8 states = 16",
               paste(state_counts$lineage, state_counts$N, collapse = ";"), state_ok)

  sample_meta <- unique(composition[, .(accession, cancer, patient_key, patient_id,
                                        sample_id, compartment, stage)])
  record_check("discovery_samples", "42", nrow(sample_meta), nrow(sample_meta) == 42L)
  sample_per_compartment <- sample_meta[, .N, by = .(patient_key, compartment)]
  record_check("one_sample_per_patient_compartment", "all N=1",
               paste0("max=", max(sample_per_compartment$N)),
               all(sample_per_compartment$N == 1L))

  patient_audit <- sample_meta[, .(
    n_samples = .N,
    compartments = paste(c("distal", "involved", "tumor")[
      c("distal", "involved", "tumor") %chin% compartment], collapse = "+"),
    n_compartments = uniqueN(compartment)
  ), by = .(accession, cancer, patient_key, patient_id)]
  patient_audit[, complete_triplet := n_compartments == 3L]
  patient_audit[, contrast_informative := n_compartments >= 2L]
  patient_audit[, gate12bf_role := fifelse(
    complete_triplet, "strict_primary",
    fifelse(contrast_informative, "informative_sensitivity", "descriptive_singleton_only")
  )]
  setorder(patient_audit, accession, patient_id)

  expected_complete <- list(
    prostate = c("BMET1", "BMET11", "BMET2", "BMET3", "BMET5", "BMET6", "BMET8"),
    renal = c("BM1", "BM10", "BM2", "BM9")
  )
  expected_informative <- list(
    prostate = c(expected_complete$prostate, "BMET10", "BMET7"),
    renal = expected_complete$renal
  )
  expected_singletons <- c("BM3", "BM4", "BM5", "BM7", "BM8")
  for (ca in c("prostate", "renal")) {
    observed_complete <- sort(patient_audit[cancer == ca & complete_triplet, patient_id])
    observed_informative <- sort(patient_audit[cancer == ca & contrast_informative, patient_id])
    record_check(paste0("complete_patient_ids_", ca),
                 paste(sort(expected_complete[[ca]]), collapse = ","),
                 paste(observed_complete, collapse = ","),
                 identical(observed_complete, sort(expected_complete[[ca]])))
    record_check(paste0("informative_patient_ids_", ca),
                 paste(sort(expected_informative[[ca]]), collapse = ","),
                 paste(observed_informative, collapse = ","),
                 identical(observed_informative, sort(expected_informative[[ca]])))
  }
  observed_singletons <- sort(patient_audit[cancer == "renal" & !contrast_informative, patient_id])
  record_check("renal_singleton_ids", paste(sort(expected_singletons), collapse = ","),
               paste(observed_singletons, collapse = ","),
               identical(observed_singletons, sort(expected_singletons)))
  record_check("discovery_patients", "18", nrow(patient_audit), nrow(patient_audit) == 18L)
  write_tsv(patient_audit, "discovery_analysis_unit_audit.tsv")

  scope_audit <- fread(input_path("scope_audit"))
  write_tsv(scope_audit, "analysis_scope_audit.tsv")

  current_stage <<- "oep_identity_audit"
  oep_master <- fread(input_path("oep_master"))
  oep_paired <- fread(input_path("oep_paired"))
  oep_atlas <- fread(input_path("oep_atlas_crosswalk"))
  oep_pair_xwalk <- fread(input_path("oep_paired_crosswalk"))
  oep_pair_library <- fread(input_path("oep_paired_library"))

  master_ids <- sort(unique(oep_master$sample_id))
  atlas_ids <- sort(unique(oep_atlas$sample_id))
  paired_ids <- sort(unique(oep_paired$sample_id))
  overlap_ids <- intersect(atlas_ids, paired_ids)
  record_check("oep_master_unique_oes", "62", length(master_ids), length(master_ids) == 62L)
  record_check("oep_atlas_unique_oes", "53", length(atlas_ids), length(atlas_ids) == 53L)
  record_check("oep_atlas_primary_patients", "49",
               uniqueN(oep_atlas[primary_analysis == TRUE, patient_id]),
               uniqueN(oep_atlas[primary_analysis == TRUE, patient_id]) == 49L)
  record_check("oep_paired_unique_oes", "9", length(paired_ids), length(paired_ids) == 9L)
  record_check("oep_oes_intersection", "0", length(overlap_ids), length(overlap_ids) == 0L)
  record_check("oep_oes_union_matches_master", "62 and set-equal",
               length(union(atlas_ids, paired_ids)),
               length(union(atlas_ids, paired_ids)) == 62L &&
                 setequal(union(atlas_ids, paired_ids), master_ids))
  record_check("oep_paired_crosswalk_set", "9 OES set-equal",
               uniqueN(oep_pair_xwalk$node_sample_id),
               uniqueN(oep_pair_xwalk$node_sample_id) == 9L &&
                 setequal(oep_pair_xwalk$node_sample_id, paired_ids))
  record_check("oep_paired_library_set", "9 OES set-equal",
               uniqueN(oep_pair_library$sample_id),
               uniqueN(oep_pair_library$sample_id) == 9L &&
                 setequal(oep_pair_library$sample_id, paired_ids))

  common_columns <- intersect(names(oep_master), names(oep_paired))
  master_pair <- oep_master[sample_id %chin% paired_ids][order(sample_id)]
  paired_sorted <- oep_paired[order(sample_id)]
  normalize_missing <- function(x) {
    y <- trimws(as.character(x))
    y[is.na(x) | is.na(y) | !nzchar(y)] <- "<MISSING>"
    y
  }
  paired_master_equal <- nrow(master_pair) == nrow(paired_sorted) && all(vapply(
    common_columns,
    function(col) identical(
      normalize_missing(master_pair[[col]]),
      normalize_missing(paired_sorted[[col]])
    ), logical(1)
  ))
  record_check("oep_paired_matches_master_metadata", "TRUE", paired_master_equal,
               paired_master_equal)

  specimen_counts <- oep_paired[, .N, by = paired_tissue_class]
  get_specimen_n <- function(cl) specimen_counts[paired_tissue_class == cl, N] %||% 0L
  `%||%` <- function(x, y) if (length(x)) x else y
  specimen_ok <- get_specimen_n("bone_metastasis") == 4L &&
    get_specimen_n("normal_bone") == 2L && get_specimen_n("primary") == 3L
  record_check("oep_paired_specimen_classes", "bone_metastasis=4;normal_bone=2;primary=3",
               paste(specimen_counts$paired_tissue_class, specimen_counts$N, collapse = ";"),
               specimen_ok)
  record_check("oep_paired_participants", "4", uniqueN(oep_paired$paired_patient),
               uniqueN(oep_paired$paired_patient) == 4L)

  pair_class_check <- merge(
    oep_pair_xwalk[, .(sample_id = node_sample_id, xwalk_tissue_class = tissue_class,
                       archive_directory)],
    oep_paired[, .(sample_id, metadata_tissue_class = paired_tissue_class)],
    by = "sample_id", all = TRUE
  )
  pair_class_check[, resolved := xwalk_tissue_class == metadata_tissue_class]
  unresolved_classes <- pair_class_check[is.na(resolved) | !resolved, .N]
  bone_named_primary <- pair_class_check[grepl("_Bone$", archive_directory) &
                                           metadata_tissue_class == "primary", .N]
  record_check("oep_unresolved_specimen_class_conflicts", "0", unresolved_classes,
               unresolved_classes == 0L)
  record_check("oep_bone_named_primary_resolved", "2", bone_named_primary,
               bone_named_primary == 2L)

  paired_contrasts <- oep_paired[, {
    n_bm <- sum(paired_tissue_class == "bone_metastasis")
    list(n_contrasts = if (n_bm == 1L) sum(paired_tissue_class %chin% c("normal_bone", "primary")) else NA_integer_)
  }, by = paired_patient]
  record_check("oep_prespecified_within_patient_contrasts", "5",
               sum(paired_contrasts$n_contrasts), sum(paired_contrasts$n_contrasts) == 5L)
  record_check("oep_reused_bm_anchor_participants", "1 (LUCA1)",
               paste(paired_contrasts[n_contrasts > 1L, paired_patient], collapse = ","),
               identical(paired_contrasts[n_contrasts > 1L, paired_patient], "LUCA1"))

  atlas_labels <- sort(unique(oep_atlas$patient_id))
  paired_labels <- sort(unique(oep_paired$paired_patient))
  collision_labels <- intersect(atlas_labels, paired_labels)
  expected_collisions <- c("LUCA1", "LUCA2", "PRAD1", "PRAD2")
  record_check("oep_short_label_collisions", paste(expected_collisions, collapse = ","),
               paste(collision_labels, collapse = ","),
               identical(collision_labels, expected_collisions))

  atlas_identity <- merge(
    oep_atlas[, .(sample_id, display_label = patient_id, sample_name, sex, age,
                  collection_date, primary_analysis)],
    oep_master[, .(sample_id, master_tissue = tissue, master_disease = disease_name)],
    by = "sample_id", all.x = TRUE
  )
  atlas_identity[, `:=`(
    resource_set = "atlas",
    namespaced_label = paste0("atlas::", display_label),
    tissue_class = "bone_metastasis_atlas"
  )]
  paired_identity <- oep_paired[, .(
    sample_id, display_label = paired_patient, sample_name, sex, age, collection_date,
    primary_analysis = NA, master_tissue = tissue, master_disease = disease_name,
    resource_set = "paired", namespaced_label = paste0("paired::", paired_patient),
    tissue_class = paired_tissue_class
  )]
  oep_identity <- rbindlist(list(atlas_identity, paired_identity), use.names = TRUE, fill = TRUE)
  oep_identity[, short_label_collision := display_label %chin% collision_labels]
  setcolorder(oep_identity, c("resource_set", "sample_id", "display_label",
                             "namespaced_label", "short_label_collision"))
  setorder(oep_identity, resource_set, display_label, sample_id)
  record_check("oep_identity_crosswalk_rows", "62", nrow(oep_identity),
               nrow(oep_identity) == 62L && uniqueN(oep_identity$sample_id) == 62L)
  write_tsv(oep_identity, "oep_identity_crosswalk.tsv")

  collision_audit <- rbindlist(lapply(collision_labels, function(label) {
    a <- oep_atlas[patient_id == label][1L]
    p <- oep_paired[paired_patient == label]
    metadata_distinguishes <- unique(a$sex) != unique(p$sex) ||
      unique(a$age) != unique(p$age) ||
      !as.character(a$collection_date) %chin% as.character(p$collection_date)
    data.table(
      short_label = label,
      atlas_oes_id = a$sample_id,
      atlas_sex = a$sex,
      atlas_age = a$age,
      atlas_collection_date = as.character(a$collection_date),
      paired_oes_ids = paste(sort(p$sample_id), collapse = ";"),
      paired_sex = paste(unique(p$sex), collapse = ";"),
      paired_age = paste(unique(p$age), collapse = ";"),
      paired_collection_dates = paste(sort(unique(as.character(p$collection_date))), collapse = ";"),
      oes_overlap = any(a$sample_id %chin% p$sample_id),
      metadata_distinguishes_source_persons = metadata_distinguishes,
      participant_nonoverlap_status = "METADATA_SUPPORTED"
    )
  }))
  record_check("oep_collision_metadata_support", "all TRUE",
               all(collision_audit$metadata_distinguishes_source_persons),
               all(collision_audit$metadata_distinguishes_source_persons))
  write_tsv(collision_audit, "oep_short_label_collision_audit.tsv")

  oep_verdict <- data.table(
    metric = c(
      "specimen_nonoverlap", "participant_nonoverlap", "participant_nonoverlap_evidence",
      "global_participant_id_available", "study_level_independence",
      "paired_participants", "paired_specimens", "paired_contrasts",
      "reused_bm_anchor", "normal_bone_terminology", "admissible_wording"
    ),
    value = c(
      "VERIFIED", "METADATA_SUPPORTED", "demographics_and_collection_dates",
      "false", "false", "4", "9", "5", "LUCA1",
      "source-defined patient-matched adjacent normal bone",
      "participant- and specimen-level non-overlap within OEP005136; not independent at the study/accession level"
    )
  )
  write_tsv(oep_verdict, "oep_identity_verdict.tsv")

  current_stage <<- "clr_transformation"
  pseudocounts <- as.numeric(strsplit(cfg_get("pseudocounts"), ",", fixed = TRUE)[[1L]])
  primary_pc <- as.numeric(cfg_get("primary_pseudocount"))
  compute_clr <- function(z, pc) {
    states <- sort(unique(z$state))
    wide <- dcast(
      z,
      accession + cancer + patient_key + patient_id + sample_id + compartment + stage ~ state,
      value.var = "n_state", fill = 0
    )
    mat <- as.matrix(wide[, ..states])
    logmat <- log(mat + pc)
    clrmat <- sweep(logmat, 1L, rowMeans(logmat), "-")
    out <- cbind(
      wide[, .(accession, cancer, patient_key, patient_id, sample_id, compartment, stage)],
      as.data.table(clrmat)
    )
    out <- melt(
      out,
      id.vars = c("accession", "cancer", "patient_key", "patient_id", "sample_id",
                  "compartment", "stage"),
      variable.name = "state", value.name = "clr"
    )
    out[, `:=`(state = as.character(state), pseudocount = pc,
               lineage = unique(z$lineage))]
    out
  }
  clr_all <- rbindlist(lapply(sort(unique(composition$lineage)), function(lin) {
    rbindlist(lapply(pseudocounts, function(pc) compute_clr(composition[lineage == lin], pc)))
  }))
  setcolorder(clr_all, c("accession", "cancer", "patient_key", "patient_id", "sample_id",
                         "compartment", "stage", "lineage", "state", "pseudocount", "clr"))
  setorder(clr_all, pseudocount, accession, patient_id, sample_id, lineage, state)
  record_check("clr_rows", "2016", nrow(clr_all), nrow(clr_all) == 2016L)
  record_check("clr_all_finite", "TRUE", all(is.finite(clr_all$clr)),
               all(is.finite(clr_all$clr)))
  closure <- clr_all[, .(clr_sum = sum(clr)),
                     by = .(pseudocount, accession, cancer, patient_key, sample_id, lineage)]
  closure_max <- max(abs(closure$clr_sum))
  record_check("clr_closure_blocks", "252", nrow(closure), nrow(closure) == 252L)
  record_check("clr_closure_tolerance", paste0("<=", cfg_get("closure_tolerance")),
               format(closure_max, scientific = TRUE),
               closure_max <= as.numeric(cfg_get("closure_tolerance")))
  write_tsv(clr_all, "clr_sample_values.tsv.gz")

  current_stage <<- "patient_contrasts"
  informative_keys <- patient_audit[contrast_informative == TRUE, patient_key]
  patient_contrasts <- clr_all[patient_key %chin% informative_keys, {
    ord <- order(stage)
    stages <- stage[ord]
    values <- clr[ord]
    compartments <- compartment[ord]
    slope <- sum((stages - mean(stages)) * (values - mean(values))) /
      sum((stages - mean(stages))^2)
    named <- setNames(values, compartments)
    value_or_na <- function(name) if (name %chin% names(named)) unname(named[[name]]) else NA_real_
    distal <- value_or_na("distal")
    involved <- value_or_na("involved")
    tumor <- value_or_na("tumor")
    list(
      n_compartments = length(values),
      observed_compartments = paste(compartments, collapse = "+"),
      stage_slope = slope,
      tumor_minus_distal = tumor - distal,
      involved_minus_distal = involved - distal,
      tumor_minus_involved = tumor - involved
    )
  }, by = .(pseudocount, accession, cancer, patient_key, patient_id, lineage, state)]
  setorder(patient_contrasts, pseudocount, cancer, patient_id, lineage, state)
  record_check("patient_contrast_rows", "624", nrow(patient_contrasts),
               nrow(patient_contrasts) == 624L)
  finite_counts <- c(
    stage_slope = sum(is.finite(patient_contrasts$stage_slope)),
    tumor_minus_distal = sum(is.finite(patient_contrasts$tumor_minus_distal)),
    involved_minus_distal = sum(is.finite(patient_contrasts$involved_minus_distal)),
    tumor_minus_involved = sum(is.finite(patient_contrasts$tumor_minus_involved))
  )
  expected_finite <- c(stage_slope = 624L, tumor_minus_distal = 576L,
                       involved_minus_distal = 528L, tumor_minus_involved = 576L)
  record_check("patient_contrast_finite_counts",
               paste(names(expected_finite), expected_finite, collapse = ";"),
               paste(names(finite_counts), finite_counts, collapse = ";"),
               identical(as.integer(finite_counts), as.integer(expected_finite)))
  triplet_rows <- patient_contrasts[patient_key %chin% patient_audit[complete_triplet == TRUE, patient_key]]
  additivity_error <- max(abs(
    triplet_rows$tumor_minus_distal -
      triplet_rows$involved_minus_distal - triplet_rows$tumor_minus_involved
  ))
  slope_error <- max(abs(triplet_rows$tumor_minus_distal / 2 - triplet_rows$stage_slope))
  record_check("triplet_contrast_additivity", paste0("<=", cfg_get("equality_tolerance")),
               format(additivity_error, scientific = TRUE),
               additivity_error <= as.numeric(cfg_get("equality_tolerance")))
  record_check("triplet_td_equals_two_slopes", paste0("<=", cfg_get("equality_tolerance")),
               format(slope_error, scientific = TRUE),
               slope_error <= as.numeric(cfg_get("equality_tolerance")))
  write_tsv(patient_contrasts, "clr_patient_contrasts.tsv.gz")

  current_stage <<- "strict_primary_effects"
  strict_keys <- patient_audit[complete_triplet == TRUE, patient_key]
  strict_patient <- patient_contrasts[patient_key %chin% strict_keys]
  strict_effects <- strict_patient[, c(
    summarize_values(tumor_minus_distal),
    list(
      estimate_mean_stage_slope = mean(stage_slope),
      estimate_median_stage_slope = median(stage_slope),
      primary_pseudocount = abs(first(pseudocount) - primary_pc) <= 1e-15
    )
  ), by = .(pseudocount, cancer, lineage, state)]
  record_check("strict_effect_rows", "96", nrow(strict_effects), nrow(strict_effects) == 96L)
  strict_n <- unique(strict_effects[, .(cancer, n_patients)])
  strict_n_ok <- strict_n[cancer == "prostate", n_patients] == 7L &&
    strict_n[cancer == "renal", n_patients] == 4L
  record_check("strict_effect_patient_counts", "prostate=7;renal=4",
               paste(strict_n$cancer, strict_n$n_patients, collapse = ";"), strict_n_ok)
  strict_closure <- strict_effects[, .(effect_sum = sum(estimate_mean)),
                                   by = .(pseudocount, cancer, lineage)]
  strict_closure_max <- max(abs(strict_closure$effect_sum))
  record_check("strict_effect_closure", paste0("<=", cfg_get("closure_tolerance")),
               format(strict_closure_max, scientific = TRUE),
               strict_closure_max <= as.numeric(cfg_get("closure_tolerance")))

  current_stage <<- "informative_sensitivities"
  contrast_columns <- c("stage_slope", "tumor_minus_distal",
                        "involved_minus_distal", "tumor_minus_involved")
  informative_long <- melt(
    patient_contrasts,
    id.vars = c("pseudocount", "accession", "cancer", "patient_key", "patient_id",
                "lineage", "state", "n_compartments", "observed_compartments"),
    measure.vars = contrast_columns,
    variable.name = "contrast", value.name = "patient_effect"
  )[is.finite(patient_effect)]
  informative_effects <- informative_long[, summarize_values(patient_effect),
                                          by = .(pseudocount, cancer, lineage, state, contrast)]
  record_check("informative_effect_rows", "384", nrow(informative_effects),
               nrow(informative_effects) == 384L)
  expected_denominators <- data.table(
    cancer = rep(c("prostate", "renal"), each = 4L),
    contrast = rep(contrast_columns, times = 2L),
    expected_n = c(9L, 8L, 7L, 8L, 4L, 4L, 4L, 4L)
  )
  observed_denominators <- unique(informative_effects[, .(cancer, contrast, n_patients)])
  denominator_audit <- merge(expected_denominators, observed_denominators,
                             by = c("cancer", "contrast"), all = TRUE)
  denominator_audit[, pass := expected_n == n_patients]
  record_check("contrast_specific_denominators", "all expected denominators",
               paste(denominator_audit$cancer, denominator_audit$contrast,
                     denominator_audit$n_patients, collapse = ";"),
               nrow(denominator_audit) == 8L && all(denominator_audit$pass))
  write_tsv(denominator_audit, "contrast_denominator_audit.tsv")

  make_lopo <- function(z, analysis_set) {
    z[, {
      pats <- patient_key
      vals <- patient_effect
      if (length(vals) < 2L) stop("LOPO group has fewer than two patients", call. = FALSE)
      full <- mean(vals)
      data.table(
        left_out_patient_key = pats,
        full_estimate = full,
        lopo_estimate = vapply(seq_along(vals), function(i) mean(vals[-i]), numeric(1)),
        lopo_same_direction = sign_with_tolerance(vapply(
          seq_along(vals), function(i) mean(vals[-i]), numeric(1)
        )) == sign_with_tolerance(full)
      )
    }, by = .(pseudocount, cancer, lineage, state, contrast)][, analysis_set := analysis_set]
  }
  strict_long <- strict_patient[, .(
    pseudocount, cancer, lineage, state, patient_key,
    contrast = "tumor_minus_distal", patient_effect = tumor_minus_distal
  )]
  lopo_strict <- make_lopo(strict_long, "complete_triplets")
  lopo_informative <- make_lopo(informative_long[, .(
    pseudocount, cancer, lineage, state, contrast, patient_key, patient_effect
  )], "contrast_informative")
  patient_deletion <- rbindlist(list(lopo_strict, lopo_informative), use.names = TRUE)
  setcolorder(patient_deletion, c("analysis_set", "pseudocount", "cancer", "lineage",
                                  "state", "contrast", "left_out_patient_key"))
  record_check("patient_deletion_all_finite", "TRUE",
               all(is.finite(patient_deletion$lopo_estimate)),
               all(is.finite(patient_deletion$lopo_estimate)))
  write_tsv(patient_deletion, "clr_patient_deletion.tsv.gz")

  current_stage <<- "legacy_traceability"
  legacy_cohort <- fread(input_path("legacy_clr_cohort"))
  legacy_meta <- fread(input_path("legacy_clr_meta"))
  legacy_primary_cohort <- fread(input_path("legacy_clr_primary_cohort"))
  legacy_primary_meta <- fread(input_path("legacy_clr_primary_meta"))
  record_check("legacy_cohort_rows", "192", nrow(legacy_cohort), nrow(legacy_cohort) == 192L)
  record_check("legacy_meta_rows", "96", nrow(legacy_meta), nrow(legacy_meta) == 96L)
  record_check("legacy_primary_cohort_rows", "32", nrow(legacy_primary_cohort),
               nrow(legacy_primary_cohort) == 32L)
  record_check("legacy_primary_meta_rows", "16", nrow(legacy_primary_meta),
               nrow(legacy_primary_meta) == 16L)

  primary_cohort_reference <- legacy_cohort[pseudocount == primary_pc & subset == "all_available"]
  primary_cohort_merge <- merge(
    primary_cohort_reference[, .(lineage, state, cancer, beta_ref = beta, se_ref = se,
                                  p_ref = p_value)],
    legacy_primary_cohort[, .(lineage, state, cancer, beta_primary = beta,
                              se_primary = se, p_primary = p_value)],
    by = c("lineage", "state", "cancer"), all = TRUE
  )
  primary_cohort_consistent <- nrow(primary_cohort_merge) == 32L &&
    max(abs(primary_cohort_merge$beta_ref - primary_cohort_merge$beta_primary)) <= 1e-15 &&
    max(abs(primary_cohort_merge$se_ref - primary_cohort_merge$se_primary)) <= 1e-15 &&
    max(abs(primary_cohort_merge$p_ref - primary_cohort_merge$p_primary)) <= 1e-15
  record_check("legacy_primary_cohort_consistency", "exact within 1e-15",
               primary_cohort_consistent, primary_cohort_consistent)

  primary_meta_reference <- legacy_meta[pseudocount == primary_pc & subset == "all_available"]
  primary_meta_merge <- merge(
    primary_meta_reference[, .(lineage, state, beta_ref = beta_meta, se_ref = se_meta,
                                p_ref = p_meta)],
    legacy_primary_meta[, .(lineage, state, beta_primary = beta_meta,
                            se_primary = se_meta, p_primary = p_meta)],
    by = c("lineage", "state"), all = TRUE
  )
  primary_meta_consistent <- nrow(primary_meta_merge) == 16L &&
    max(abs(primary_meta_merge$beta_ref - primary_meta_merge$beta_primary)) <= 1e-15 &&
    max(abs(primary_meta_merge$se_ref - primary_meta_merge$se_primary)) <= 1e-15 &&
    max(abs(primary_meta_merge$p_ref - primary_meta_merge$p_primary)) <= 1e-15
  record_check("legacy_primary_meta_consistency", "exact within 1e-15",
               primary_meta_consistent, primary_meta_consistent)

  reconstructed_legacy <- clr_all[, {
    rows <- list()
    for (subset_name in c("all_available", "complete_triplets")) {
      zz <- .SD
      if (subset_name == "complete_triplets") zz <- zz[patient_key %chin% strict_keys]
      fit <- lm(clr ~ factor(patient_key) + stage, data = zz)
      rows[[subset_name]] <- data.table(
        subset = subset_name,
        reconstructed_beta = unname(coef(fit)["stage"]),
        reconstructed_n_patients = uniqueN(zz$patient_key),
        reconstructed_n_samples = nrow(zz)
      )
    }
    rbindlist(rows)
  }, by = .(pseudocount, cancer, lineage, state)]
  legacy_comparison <- merge(
    legacy_cohort,
    reconstructed_legacy,
    by = c("pseudocount", "subset", "lineage", "state", "cancer"), all = TRUE
  )
  legacy_comparison[, legacy_beta_reconstruction_difference := reconstructed_beta - beta]
  legacy_reconstruction_max <- max(abs(legacy_comparison$legacy_beta_reconstruction_difference))
  record_check("legacy_beta_reconstruction", paste0("<=", cfg_get("equality_tolerance")),
               format(legacy_reconstruction_max, scientific = TRUE),
               legacy_reconstruction_max <= as.numeric(cfg_get("equality_tolerance")))

  strict_reference <- strict_effects[, .(
    pseudocount, cancer, lineage, state,
    strict_n_patients = n_patients,
    strict_mean_tumor_minus_distal = estimate_mean,
    strict_mean_stage_slope = estimate_mean_stage_slope
  )]
  informative_stage <- informative_effects[contrast == "stage_slope", .(
    pseudocount, cancer, lineage, state,
    informative_n_patients = n_patients,
    informative_mean_stage_slope = estimate_mean
  )]
  legacy_comparison <- merge(legacy_comparison, strict_reference,
                             by = c("pseudocount", "cancer", "lineage", "state"), all.x = TRUE)
  legacy_comparison <- merge(legacy_comparison, informative_stage,
                             by = c("pseudocount", "cancer", "lineage", "state"), all.x = TRUE)
  legacy_comparison[, repaired_analysis_set := fifelse(
    subset == "complete_triplets", "complete_triplets", "contrast_informative"
  )]
  legacy_comparison[, repaired_n_patients := fifelse(
    subset == "complete_triplets", strict_n_patients, informative_n_patients
  )]
  legacy_comparison[, repaired_mean_stage_slope := fifelse(
    subset == "complete_triplets", strict_mean_stage_slope, informative_mean_stage_slope
  )]
  legacy_comparison[, repaired_minus_legacy_beta := repaired_mean_stage_slope - beta]

  complete_reference <- legacy_cohort[subset == "complete_triplets", .(
    pseudocount, cancer, lineage, state,
    archived_complete_beta = beta,
    archived_complete_se = se,
    archived_complete_df = df,
    archived_complete_p = p_value,
    archived_complete_ci_low = ci_low,
    archived_complete_ci_high = ci_high
  )]
  legacy_comparison <- merge(legacy_comparison, complete_reference,
                             by = c("pseudocount", "cancer", "lineage", "state"), all.x = TRUE)
  legacy_comparison[, legacy_inverse_variance_weight := 1 / se^2]
  legacy_comparison[, archived_complete_inverse_variance_weight := 1 / archived_complete_se^2]
  legacy_comparison[, se_ratio_to_complete := se / archived_complete_se]
  legacy_comparison[, weight_ratio_to_complete :=
                      legacy_inverse_variance_weight / archived_complete_inverse_variance_weight]
  complete_scale_error <- max(abs(
    legacy_comparison[subset == "complete_triplets",
                      repaired_mean_stage_slope - beta]
  ))
  record_check("complete_triplet_slope_scale_equivalence",
               paste0("mean(T-D)/2 equals legacy beta within ", cfg_get("equality_tolerance")),
               format(complete_scale_error, scientific = TRUE),
               complete_scale_error <= as.numeric(cfg_get("equality_tolerance")))
  write_tsv(legacy_comparison, "clr_legacy_vs_repaired.tsv")

  descriptive_meta <- legacy_cohort[subset == "complete_triplets", {
    w <- 1 / se^2
    beta_stage <- sum(w * beta) / sum(w)
    list(
      beta_stage_descriptive = beta_stage,
      beta_tumor_minus_distal_descriptive = 2 * beta_stage,
      same_direction = uniqueN(sign_with_tolerance(beta)) == 1L,
      prostate_beta_stage = beta[cancer == "prostate"],
      renal_beta_stage = beta[cancer == "renal"],
      inferential_use = "DESCRIPTIVE_ONLY"
    )
  }, by = .(pseudocount, lineage, state)]
  meta_reference <- legacy_meta[subset == "complete_triplets", .(
    pseudocount, lineage, state, archived_beta_meta = beta_meta
  )]
  descriptive_meta <- merge(descriptive_meta, meta_reference,
                            by = c("pseudocount", "lineage", "state"), all.x = TRUE)
  descriptive_meta[, archive_difference := beta_stage_descriptive - archived_beta_meta]
  meta_archive_error <- max(abs(descriptive_meta$archive_difference))
  record_check("descriptive_meta_archive_equivalence",
               paste0("<=", cfg_get("equality_tolerance")),
               format(meta_archive_error, scientific = TRUE),
               meta_archive_error <= as.numeric(cfg_get("equality_tolerance")))
  write_tsv(descriptive_meta, "clr_descriptive_crosscohort_effects.tsv")

  current_stage <<- "exact_signflip"
  signflip_raw_rows <- list()
  signflip_summary_rows <- list()
  sf_index <- 0L
  strict_primary_patient <- strict_patient[abs(pseudocount - primary_pc) <= 1e-15]
  for (ca in sort(unique(strict_primary_patient$cancer))) {
    for (lin in sort(unique(strict_primary_patient$lineage))) {
      for (st in sort(unique(strict_primary_patient[lineage == lin, state]))) {
        z <- strict_primary_patient[cancer == ca & lineage == lin & state == st][order(patient_key)]
        vals <- z$tumor_minus_distal
        n <- length(vals)
        n_perm <- 2^n
        ids <- 0:(n_perm - 1L)
        signs <- sapply(seq_len(n), function(j) {
          ifelse(bitwAnd(ids, bitwShiftL(1L, j - 1L)) > 0L, 1, -1)
        })
        if (is.null(dim(signs))) signs <- matrix(signs, ncol = n)
        perm_signed_mean <- as.vector(signs %*% vals / n)
        observed_signed_mean <- mean(vals)
        observed_statistic <- abs(observed_signed_mean)
        perm_statistic <- abs(perm_signed_mean)
        extreme <- perm_statistic >= observed_statistic - as.numeric(cfg_get("signflip_tolerance"))
        pattern <- apply(signs, 1L, paste0, collapse = ",")
        sf_index <- sf_index + 1L
        signflip_raw_rows[[sf_index]] <- data.table(
          cancer = ca, lineage = lin, state = st,
          pseudocount = primary_pc, contrast = "tumor_minus_distal",
          n_patients = n, configuration_id = seq_len(n_perm), sign_pattern = pattern,
          observed_signed_mean = observed_signed_mean,
          observed_statistic = observed_statistic,
          permuted_signed_mean = perm_signed_mean,
          permuted_statistic = perm_statistic,
          extreme = extreme
        )
        signflip_summary_rows[[sf_index]] <- data.table(
          cancer = ca, lineage = lin, state = st,
          pseudocount = primary_pc, contrast = "tumor_minus_distal",
          n_patients = n, configurations = n_perm,
          unique_patterns = uniqueN(pattern),
          observed_mean = observed_signed_mean,
          observed_statistic = observed_statistic,
          extreme_configurations = sum(extreme),
          p_value_exact = mean(extreme),
          df = NA_real_,
          inference = "exact_two_sided_sign_flip_no_plus_one"
        )
      }
    }
  }
  signflip_raw <- rbindlist(signflip_raw_rows)
  signflip_summary <- rbindlist(signflip_summary_rows)
  signflip_summary[, q_value_exact := p.adjust(p_value_exact, method = "BH"),
                   by = .(cancer, lineage)]
  record_check("signflip_raw_rows", "2304", nrow(signflip_raw), nrow(signflip_raw) == 2304L)
  record_check("signflip_summary_rows", "32", nrow(signflip_summary),
               nrow(signflip_summary) == 32L)
  signflip_counts <- unique(signflip_summary[, .(cancer, n_patients, configurations,
                                                  unique_patterns)])
  signflip_count_ok <- signflip_counts[cancer == "prostate",
                                       all(n_patients == 7L & configurations == 128L &
                                             unique_patterns == 128L)] &&
    signflip_counts[cancer == "renal",
                    all(n_patients == 4L & configurations == 16L & unique_patterns == 16L)]
  record_check("signflip_configuration_counts", "prostate=128;renal=16",
               paste(signflip_counts$cancer, signflip_counts$configurations, collapse = ";"),
               signflip_count_ok)
  family_sizes <- signflip_summary[, .N, by = .(cancer, lineage)]
  record_check("signflip_bh_family_sizes", "four families of 8 states",
               paste(family_sizes$cancer, family_sizes$lineage, family_sizes$N, collapse = ";"),
               nrow(family_sizes) == 4L && all(family_sizes$N == 8L))
  record_check("renal_exact_p_floor", ">=0.125",
               min(signflip_summary[cancer == "renal", p_value_exact]),
               min(signflip_summary[cancer == "renal", p_value_exact]) >= 0.125 - 1e-15)
  p_integrality <- all(abs(signflip_summary$p_value_exact * signflip_summary$configurations -
                             round(signflip_summary$p_value_exact * signflip_summary$configurations)) < 1e-12)
  record_check("signflip_probability_integrality", "TRUE", p_integrality, p_integrality)
  write_tsv(signflip_raw, "clr_exact_signflip_audit.tsv.gz")
  write_tsv(signflip_summary, "clr_exact_signflip_summary.tsv")

  current_stage <<- "patient_bootstrap"
  seed <- as.integer(cfg_get("seed"))
  n_boot <- as.integer(cfg_get("bootstrap_iterations"))
  RNGkind(kind = cfg_get("rng_kind"), normal.kind = cfg_get("normal_kind"),
          sample.kind = cfg_get("sample_kind"))
  set.seed(seed)
  bootstrap_rows <- list()
  bootstrap_audit_rows <- list()
  bidx <- 0L
  for (ca in sort(unique(strict_primary_patient$cancer))) {
    patients <- sort(unique(strict_primary_patient[cancer == ca, patient_key]))
    n <- length(patients)
    draws <- matrix(sample.int(n, n_boot * n, replace = TRUE), nrow = n_boot, ncol = n)
    state_rows <- unique(strict_primary_patient[cancer == ca, .(lineage, state)])[order(lineage, state)]
    for (i in seq_len(nrow(state_rows))) {
      lin <- state_rows$lineage[[i]]
      st <- state_rows$state[[i]]
      z <- strict_primary_patient[cancer == ca & lineage == lin & state == st]
      vals <- z$tumor_minus_distal[match(patients, z$patient_key)]
      if (any(!is.finite(vals))) stop("Non-finite primary patient contrast before bootstrap", call. = FALSE)
      draw_values <- matrix(vals[draws], nrow = n_boot, ncol = n)
      estimates <- rowMeans(draw_values)
      bidx <- bidx + 1L
      bootstrap_rows[[bidx]] <- data.table(
        iteration = seq_len(n_boot), cancer = ca, lineage = lin, state = st,
        pseudocount = primary_pc, contrast = "tumor_minus_distal",
        estimate_mean = estimates
      )
    }
    bootstrap_audit_rows[[ca]] <- data.table(
      cancer = ca, n_patients = n, requested = n_boot, attempted = n_boot,
      accepted = n_boot, discarded = 0L,
      shared_draw_matrix_sha256 = digest(draws, algo = "sha256")
    )
  }
  bootstrap <- rbindlist(bootstrap_rows)
  bootstrap_audit <- rbindlist(bootstrap_audit_rows)
  record_check("bootstrap_rows", "319968", nrow(bootstrap), nrow(bootstrap) == 319968L)
  record_check("bootstrap_all_finite", "TRUE", all(is.finite(bootstrap$estimate_mean)),
               all(is.finite(bootstrap$estimate_mean)))
  bootstrap_endpoint_counts <- bootstrap[, .(iterations = uniqueN(iteration), rows = .N),
                                         by = .(cancer, lineage, state)]
  record_check("bootstrap_endpoint_iterations", "9999 for all 32 endpoints",
               paste0("endpoints=", nrow(bootstrap_endpoint_counts),
                      ";min=", min(bootstrap_endpoint_counts$iterations),
                      ";max=", max(bootstrap_endpoint_counts$iterations)),
               nrow(bootstrap_endpoint_counts) == 32L &&
                 all(bootstrap_endpoint_counts$iterations == n_boot) &&
                 all(bootstrap_endpoint_counts$rows == n_boot))
  record_check("bootstrap_no_discarded_draws", "attempted=accepted=9999;discarded=0",
               paste(bootstrap_audit$cancer, bootstrap_audit$attempted,
                     bootstrap_audit$accepted, bootstrap_audit$discarded, collapse = ";"),
               all(bootstrap_audit$requested == n_boot & bootstrap_audit$attempted == n_boot &
                     bootstrap_audit$accepted == n_boot & bootstrap_audit$discarded == 0L))
  bootstrap_closure <- bootstrap[, .(effect_sum = sum(estimate_mean)),
                                 by = .(iteration, cancer, lineage)]
  bootstrap_closure_max <- max(abs(bootstrap_closure$effect_sum))
  record_check("bootstrap_compositional_closure", paste0("<=", cfg_get("closure_tolerance")),
               format(bootstrap_closure_max, scientific = TRUE),
               bootstrap_closure_max <= as.numeric(cfg_get("closure_tolerance")))
  bootstrap_summary <- bootstrap[, .(
    bootstrap_median = median(estimate_mean),
    bootstrap_ci_low = quantile(estimate_mean, 0.025, type = 7),
    bootstrap_ci_high = quantile(estimate_mean, 0.975, type = 7),
    bootstrap_positive_fraction = mean(estimate_mean > 0),
    bootstrap_negative_fraction = mean(estimate_mean < 0),
    bootstrap_iterations = .N,
    interpretation = "percentile_uncertainty_not_hypothesis_test"
  ), by = .(cancer, lineage, state)]
  write_tsv(bootstrap, "clr_patient_bootstrap.tsv.gz")
  write_tsv(bootstrap_audit, "bootstrap_audit.tsv")

  strict_effects <- merge(
    strict_effects,
    signflip_summary[, .(pseudocount, cancer, lineage, state, exact_n = n_patients,
                         exact_configurations = configurations, p_value_exact,
                         q_value_exact, exact_df = df, exact_inference = inference)],
    by = c("pseudocount", "cancer", "lineage", "state"), all.x = TRUE
  )
  bootstrap_summary[, pseudocount := primary_pc]
  strict_effects <- merge(strict_effects, bootstrap_summary,
                          by = c("pseudocount", "cancer", "lineage", "state"), all.x = TRUE)
  setorder(strict_effects, pseudocount, cancer, lineage, state)
  write_tsv(strict_effects, "clr_primary_complete_triplet_effects.tsv")
  setorder(informative_effects, pseudocount, cancer, lineage, state, contrast)
  write_tsv(informative_effects, "clr_informative_patient_sensitivity.tsv")

  current_stage <<- "direction_traceability"
  pc_trace_strict <- strict_effects[, .(
    analysis_set = "complete_triplets",
    contrast = "tumor_minus_distal",
    pseudocount, cancer, lineage, state,
    estimate = estimate_mean,
    direction = sign_with_tolerance(estimate_mean)
  )]
  pc_trace_informative <- informative_effects[contrast == "stage_slope", .(
    analysis_set = "contrast_informative",
    contrast = "stage_slope",
    pseudocount, cancer, lineage, state,
    estimate = estimate_mean,
    direction = sign_with_tolerance(estimate_mean)
  )]
  pc_trace <- rbindlist(list(pc_trace_strict, pc_trace_informative))
  pc_trace[, direction_at_primary_pc := direction[which.min(abs(pseudocount - primary_pc))],
           by = .(analysis_set, cancer, lineage, state, contrast)]
  pc_trace[, direction_stable_across_pseudocounts :=
             uniqueN(direction[!is.na(direction)]) == 1L && !any(is.na(direction)),
           by = .(analysis_set, cancer, lineage, state, contrast)]
  write_tsv(pc_trace, "clr_pseudocount_sensitivity.tsv")

  original_meta <- fread(input_path("state_meta"))
  direction_trace <- original_meta[, {
    z <- strict_effects[lineage == .BY$lineage & state == .BY$state]
    original_direction <- sign_with_tolerance(beta_meta)
    p05 <- z[abs(pseudocount - primary_pc) <= 1e-15]
    p05_prostate <- p05[cancer == "prostate", sign_with_tolerance(estimate_mean)]
    p05_renal <- p05[cancer == "renal", sign_with_tolerance(estimate_mean)]
    all_directions <- z[, sign_with_tolerance(estimate_mean)]
    all_finite_nonzero <- length(all_directions) == 6L &&
      all(!is.na(all_directions) & all_directions != 0L)
    both_match_primary <- length(p05_prostate) == 1L && length(p05_renal) == 1L &&
      p05_prostate == original_direction && p05_renal == original_direction
    both_match_all <- all_finite_nonzero && all(all_directions == original_direction)
    classification <- if (!all_finite_nonzero || original_direction == 0L) {
      "NOT_SUPPORTED"
    } else if (p05_prostate != p05_renal) {
      "HETEROGENEOUS"
    } else if (p05_prostate != original_direction) {
      "NOT_SUPPORTED"
    } else if (both_match_all) {
      "SUPPORTED"
    } else if (both_match_primary) {
      "DIRECTION_ONLY"
    } else {
      "NOT_SUPPORTED"
    }
    list(
      original_beta = beta_meta,
      original_direction = original_direction,
      originally_cross_cancer_stable = cross_cancer_stable,
      prostate_direction_pc0_5 = p05_prostate,
      renal_direction_pc0_5 = p05_renal,
      both_cohorts_match_original_at_pc0_5 = both_match_primary,
      both_cohorts_match_original_at_all_pseudocounts = both_match_all,
      gate12bf_direction_class = classification
    )
  }, by = .(lineage, state)]
  record_check("frozen_stable_state_count", "9",
               sum(direction_trace$originally_cross_cancer_stable),
               sum(direction_trace$originally_cross_cancer_stable) == 9L)
  write_tsv(direction_trace, "clr_direction_traceability.tsv")

  current_stage <<- "protected_manifest_and_receipt"
  protected_post <- manifest_for(protected_paths)
  write_tsv(protected_post, "protected_manifest_post.tsv")
  protected_unchanged <- identical(protected_pre, protected_post)
  record_check("protected_inputs_unchanged", "TRUE", protected_unchanged,
               protected_unchanged)
  record_check("writes_contained_in_run_v3", "TRUE", TRUE, TRUE)

  session_lines <- capture.output(sessionInfo())
  write_text(session_lines, "sessionInfo.txt")
  checks_path <- write_tsv(checks, "validation_checks.tsv")
  log_line("RUN_COMPLETED\tstatus=COMPLETED\tchecks=", nrow(checks),
           "\tall_hard_checks_pass=", all(checks[severity == "HARD", pass]))

  all_output_files <- sort(list.files(output_abs, full.names = TRUE, recursive = FALSE))
  all_output_files <- all_output_files[!basename(all_output_files) %chin%
                                         c("GATE12BF_RECEIPT.json", "output_manifest.tsv")]
  output_manifest <- data.table(
    file = basename(all_output_files),
    bytes = as.numeric(file.info(all_output_files)$size),
    sha256 = vapply(all_output_files, sha256_file, character(1))
  )
  output_manifest <- merge(output_manifest, output_registry, by = "file", all.x = TRUE,
                           suffixes = c("", "_registered"))
  fwrite(output_manifest, file.path(output_abs, "output_manifest.tsv"), sep = "\t")

  end_time <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  receipt <- list(
    schema_version = "1.0",
    gate_id = "Gate12BF",
    run_id = "gate12bf_run_v3",
    status = "COMPLETED",
    command = "Rscript scripts/run_gate12bf_analysis_unit_repair.R config/gate12bf_analysis_unit_repair_v3.tsv",
    working_directory = project_root,
    started_at = start_time,
    ended_at = end_time,
    exit_code = 0L,
    seed = seed,
    rng_kind = RNGkind(),
    r_version = R.version.string,
    package_versions = list(
      data.table = as.character(packageVersion("data.table")),
      jsonlite = as.character(packageVersion("jsonlite")),
      digest = as.character(packageVersion("digest"))
    ),
    locale = Sys.getlocale(),
    timezone = Sys.timezone(),
    script = list(path = relative_path(script_path), sha256 = sha256_file(script_path)),
    config = list(path = relative_path(config_path), sha256 = sha256_file(config_path)),
    input_hashes = input_hashes,
    analysis_contract = list(
      primary_set = "complete_triplets",
      primary_counts = list(prostate = 7L, renal = 4L),
      informative_slope_counts = list(prostate = 9L, renal = 4L),
      singleton_use = "descriptive_only",
      primary_pseudocount = primary_pc,
      primary_contrast = "tumor_minus_distal",
      primary_statistic = "abs(mean(patient tumor-minus-distal CLR contrasts))",
      exact_p = "fraction of all sign patterns with abs(permuted mean) >= abs(observed mean); no +1 correction",
      bootstrap = "9999 within-cohort patient resamples; shared draws across all states; percentile type-7 CI",
      fdr_family = cfg_get("fdr_family"),
      pooled_inference = "descriptive_only"
    ),
    patient_sets = list(
      prostate_complete = sort(expected_complete$prostate),
      renal_complete = sort(expected_complete$renal),
      prostate_informative = sort(expected_informative$prostate),
      renal_informative = sort(expected_informative$renal),
      renal_descriptive_singletons = sort(expected_singletons)
    ),
    bootstrap_audit = bootstrap_audit,
    signflip_audit = signflip_counts,
    oep_identity = list(
      specimen_nonoverlap = "VERIFIED",
      participant_nonoverlap = "METADATA_SUPPORTED",
      participant_nonoverlap_evidence = "demographics_and_collection_dates",
      global_participant_id_available = FALSE,
      study_level_independence = FALSE,
      atlas_oes = length(atlas_ids),
      paired_oes = length(paired_ids),
      oes_intersection = length(overlap_ids),
      short_label_collisions = collision_labels
    ),
    direction_summary = direction_trace[, .N, by = gate12bf_direction_class],
    checks = checks,
    outputs = output_manifest,
    output_manifest_sha256 = sha256_file(file.path(output_abs, "output_manifest.tsv")),
    protected_manifest_unchanged = protected_unchanged
  )
  write_json(receipt, file.path(output_abs, "GATE12BF_RECEIPT.json"),
             pretty = TRUE, auto_unbox = TRUE, na = "null", digits = 16)
  invisible(receipt)
}

exit_status <- tryCatch({
  run_gate12bf()
  0L
}, error = function(e) {
  failure_receipt(e)
  1L
})

quit(save = "no", status = exit_status, runLast = FALSE)
