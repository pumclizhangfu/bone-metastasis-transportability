#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else
  normalizePath("results/gate12c_external_validation", mustWork = TRUE)

primary <- fread(file.path(out, "gate12c_target_freeze.tsv"))
sensitivity <- fread(file.path(out, "external_threshold_sensitivity.tsv"))
winner_grid <- unique(sensitivity[, .(min_cells, min_detect, scenario_winner)])
primary_axis <- if (nrow(primary)) primary$axis[[1L]] else NA_character_
retention <- if (nrow(primary)) mean(winner_grid$scenario_winner == primary_axis, na.rm = TRUE) else 0

final <- data.table(
  primary_rule_axis = primary_axis,
  primary_rule_target = if (nrow(primary)) primary$knockout_target[[1L]] else NA_character_,
  primary_rule_context = if (nrow(primary)) primary$knockout_context[[1L]] else NA_character_,
  sensitivity_scenarios = nrow(winner_grid),
  primary_winner_scenarios = sum(winner_grid$scenario_winner == primary_axis, na.rm = TRUE),
  primary_winner_retention = retention,
  robustness_threshold = 0.75,
  final_gate12d_authorized = retention >= 0.75,
  final_status = if (retention >= 0.75) "ROBUST_TARGET_FREEZE" else "PRIMARY_SELECTION_REVOKED_BY_SENSITIVITY"
)
fwrite(final, file.path(out, "gate12c_target_freeze_final.tsv"), sep = "\t")

main_audit <- readLines(file.path(out, "GATE12C_EXTERNAL_AUDIT.md"), warn = FALSE)
audit <- c(
  "# Gate12C-external final review", "",
  main_audit, "",
  "## Post-primary robustness review", "",
  paste0("- Primary-rule selected axis: ", primary_axis),
  paste0("- Primary-winner retention: ", final$primary_winner_scenarios, "/", final$sensitivity_scenarios,
         " (", sprintf("%.1f%%", 100 * retention), ")"),
  "- Required retention for final Gate12D authorisation: >=75%",
  paste0("- Final Gate12D authorised: ", toupper(as.character(final$final_gate12d_authorized))),
  paste0("- Final status: ", final$final_status), "",
  "The sensitivity audit may revoke the primary target but may not substitute a different gene."
)
writeLines(audit, file.path(out, "GATE12C_EXTERNAL_FINAL_AUDIT.md"))
cat("GATE12C_EXTERNAL_FINAL_STATUS=", final$final_status, "\n", sep = "")

