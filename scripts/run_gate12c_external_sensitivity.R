#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else
  normalizePath("results/gate12c_external_validation", mustWork = TRUE)
patients <- fread(file.path(out, "external_patient_axis_support.tsv"))
primary <- fread(file.path(out, "gate12c_external_axis_decision.tsv"))
grid <- CJ(min_cells = c(10L, 20L, 30L), min_detect = c(0.02, 0.05, 0.10))

one_scenario <- function(min_cells, min_detect) {
  z <- copy(patients)
  z[, eligible_s := !is.na(sender_cells) & !is.na(receiver_cells) &
                      sender_cells >= min_cells & receiver_cells >= min_cells]
  z[, support_s := eligible_s & ligand_detected_fraction >= min_detect &
                     receptor_detected_fraction >= min_detect]
  ds <- z[condition == "bone_metastasis", .(
    eligible = sum(eligible_s), supporting = sum(support_s),
    fraction = if (sum(eligible_s)) mean(support_s[eligible_s]) else NA_real_,
    origins = uniqueN(cancer[support_s])
  ), by = .(axis, dataset)]
  full <- merge(CJ(axis = unique(z$axis),
                   dataset = c("GSE266330", "OEP005136", "GSE225209", "GSE190772")),
                ds, by = c("axis", "dataset"), all.x = TRUE)
  full[is.na(eligible), `:=`(eligible = 0L, supporting = 0L, origins = 0L)]
  full[, pass := fifelse(dataset %chin% c("GSE266330", "OEP005136"),
                         eligible >= 10L & fraction >= 0.50 & origins >= 3L,
                         eligible >= 1L & supporting >= 1L)]
  w <- dcast(full, axis ~ dataset, value.var = c("pass", "fraction"))
  w <- merge(w, primary[, .(axis, nichenet_rank)], by = "axis")
  pass_cols <- grep("^pass_", names(w), value = TRUE)
  w[, supporting_datasets := rowSums(.SD, na.rm = TRUE), .SDcols = pass_cols]
  w[, min_core := pmin(fraction_GSE266330, fraction_OEP005136)]
  w[, external_pass := pass_GSE266330 & pass_OEP005136 & (pass_GSE225209 | pass_GSE190772)]
  e <- w[external_pass == TRUE][order(-supporting_datasets, -min_core,
                                      -fraction_OEP005136, nichenet_rank)]
  winner <- if (nrow(e)) e$axis[[1L]] else NA_character_
  w[, `:=`(min_cells = min_cells, min_detect = min_detect,
           scenario_winner = winner)]
  w
}

res <- rbindlist(lapply(seq_len(nrow(grid)), function(i)
  one_scenario(grid$min_cells[[i]], grid$min_detect[[i]])))
setorder(res, min_cells, min_detect, axis)
fwrite(res, file.path(out, "external_threshold_sensitivity.tsv"), sep = "\t")

winner <- unique(res[, .(min_cells, min_detect, scenario_winner)])
primary_axis <- primary[final_freeze_eligible == TRUE, axis]
retention <- mean(winner$scenario_winner == primary_axis, na.rm = TRUE)

winner[, min_detect_label := paste0(round(100 * min_detect), "%")]
p <- ggplot(winner, aes(factor(min_cells), factor(min_detect_label,
                   levels = c("10%", "5%", "2%")), fill = scenario_winner)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sub(" -> ", "\n", scenario_winner, fixed = TRUE)), size = 3) +
  scale_fill_manual(values = c("CXCL16 -> CXCR6" = "#B2182B",
                               "SPP1 -> CD44" = "#E69F00",
                               "CCL4 -> CCR5" = "#2166AC"), na.value = "#D9D9D9") +
  labs(title = "External target-freeze sensitivity",
       subtitle = paste0("Primary winner retention: ", scales::percent(retention, accuracy = 1)),
       x = "Minimum cells in sender and receiver state", y = "Detection threshold",
       fill = "Scenario winner") +
  theme_classic(base_size = 10, base_family = "Arial") +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
ggsave(file.path(out, "FigureS5_external_threshold_sensitivity.pdf"), p,
       width = 6.8, height = 4.8, device = cairo_pdf, bg = "white")
ggsave(file.path(out, "FigureS5_external_threshold_sensitivity.png"), p,
       width = 6.8, height = 4.8, dpi = 360, bg = "white")

audit <- c(
  "# Gate12C-external threshold-sensitivity audit", "",
  "- Primary rule remains 20 cells and 5% detection.",
  "- Sensitivity grid: 10/20/30 cells by 2%/5%/10% detection.",
  paste0("- Primary frozen axis: ", primary_axis),
  paste0("- Primary-winner retention: ", sum(winner$scenario_winner == primary_axis, na.rm = TRUE),
         "/", nrow(winner), " (", sprintf("%.1f%%", 100 * retention), ")"),
  paste0("- Scenario winners: ", paste(unique(na.omit(winner$scenario_winner)), collapse = "; ")),
  paste0("- Status: ", if (retention >= 0.75) "ROBUSTNESS_SUPPORTS_FREEZE" else "FREEZE_SENSITIVITY_FLAG")
)
writeLines(audit, file.path(out, "GATE12C_EXTERNAL_SENSITIVITY_AUDIT.md"))
cat("GATE12C_EXTERNAL_SENSITIVITY=", if (retention >= 0.75) "PASS" else "FLAG", "\n", sep = "")

