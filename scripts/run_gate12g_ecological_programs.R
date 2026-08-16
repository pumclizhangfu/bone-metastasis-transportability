#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: run_gate12g_ecological_programs.R <broad.tsv> <states.tsv> <outdir> [n_boot=199] [seed=1208]")
}

broad_file <- args[[1L]]
state_file <- args[[2L]]
outdir <- args[[3L]]
n_boot <- if (length(args) >= 4L) as.integer(args[[4L]]) else 199L
seed <- if (length(args) >= 5L) as.integer(args[[5L]]) else 1208L
run_label <- if (n_boot >= 1999L) "full" else "smoke"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

required <- function(x, cols, label) {
  miss <- setdiff(cols, names(x))
  if (length(miss)) stop(label, " missing columns: ", paste(miss, collapse = ", "))
}

close_clr <- function(m, pseudocount = 0.5) {
  m <- as.matrix(m)
  storage.mode(m) <- "double"
  p <- m + pseudocount
  p <- p / rowSums(p)
  lp <- log(p)
  lp - rowMeans(lp)
}

block_scale <- function(x, blocks, train_rows = seq_len(nrow(x))) {
  out <- x
  scales <- setNames(numeric(length(unique(blocks))), unique(blocks))
  for (b in unique(blocks)) {
    j <- which(blocks == b)
    total_var <- sum(apply(x[train_rows, j, drop = FALSE], 2L, var))
    if (!is.finite(total_var) || total_var <= 0) stop("Non-positive total variance in block ", b)
    scales[[b]] <- sqrt(total_var)
    out[, j] <- out[, j, drop = FALSE] / scales[[b]]
  }
  list(x = out, scales = scales)
}

eta_squared <- function(y, g) {
  ok <- is.finite(y) & !is.na(g)
  y <- y[ok]
  g <- droplevels(factor(g[ok]))
  if (length(y) < 3L || nlevels(g) < 2L || var(y) == 0) return(NA_real_)
  grand <- mean(y)
  ss_between <- sum(vapply(split(y, g), function(z) length(z) * (mean(z) - grand)^2, numeric(1)))
  ss_total <- sum((y - grand)^2)
  if (ss_total == 0) NA_real_ else ss_between / ss_total
}

all_permutations <- function(v) {
  if (length(v) == 1L) return(matrix(v, nrow = 1L))
  do.call(rbind, lapply(seq_along(v), function(i) {
    rest <- all_permutations(v[-i])
    cbind(v[i], rest)
  }))
}

align_loading_correlations <- function(reference, candidate) {
  k <- min(ncol(reference), ncol(candidate))
  ref <- reference[, seq_len(k), drop = FALSE]
  can <- candidate[, seq_len(k), drop = FALSE]
  cm <- abs(cor(ref, can, method = "pearson"))
  perms <- all_permutations(seq_len(k))
  totals <- apply(perms, 1L, function(p) sum(cm[cbind(seq_len(k), p)]))
  best <- perms[which.max(totals), ]
  data.table(axis = seq_len(k), matched_axis = best,
             loading_cor = cm[cbind(seq_len(k), best)])
}

theme_gate <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

broad <- fread(broad_file)
states <- fread(state_file)
required(broad, c("accession", "cancer", "patient_id", "sample_id", "compartment",
                  "broad_class", "N", "sample_total"), "broad composition")
required(states, c("accession", "cancer", "patient_id", "sample_id", "compartment",
                   "lineage", "state", "n_state", "lineage_total"), "state composition")

sample_meta <- unique(broad[, .(accession, cancer, patient_id, sample_id, compartment,
                                retained_cells = sample_total)])
if (sample_meta[, anyDuplicated(sample_id)]) stop("sample_id is not unique in broad-class metadata")

broad_feature_audit <- broad[, .(
  total_cells = sum(N),
  samples_present = sum(N > 0),
  sample_fraction_present = mean(N > 0)
), by = broad_class]
broad_feature_audit[, eligible := total_cells >= 50 & sample_fraction_present >= 0.20]

state_feature_audit <- states[, .(
  lineage = unique(lineage),
  total_cells = sum(n_state),
  samples_present = sum(n_state > 0),
  sample_fraction_present = mean(n_state > 0)
), by = state]
state_feature_audit[, eligible := total_cells >= 100 & sample_fraction_present >= 0.20]

eligible_broad <- broad_feature_audit[eligible == TRUE, broad_class]
eligible_states <- state_feature_audit[eligible == TRUE, state]

lineage_depth <- unique(states[, .(sample_id, lineage, lineage_total)])
bad_samples <- lineage_depth[lineage %chin% c("Myeloid", "T_NK") & lineage_total < 50, unique(sample_id)]
analysis_samples <- setdiff(sample_meta$sample_id, bad_samples)
meta <- sample_meta[sample_id %chin% analysis_samples]
meta[, compartment_order := match(compartment, c("distal", "involved", "tumor"))]
setorder(meta, accession, patient_id, compartment_order)
meta[, compartment_order := NULL]

bw <- dcast(broad[broad_class %chin% eligible_broad & sample_id %chin% analysis_samples],
            sample_id ~ broad_class, value.var = "N", fill = 0)
my_states <- state_feature_audit[lineage == "Myeloid" & eligible == TRUE, state]
tn_states <- state_feature_audit[lineage == "T_NK" & eligible == TRUE, state]
mw <- dcast(states[lineage == "Myeloid" & state %chin% my_states & sample_id %chin% analysis_samples],
            sample_id ~ state, value.var = "n_state", fill = 0)
tw <- dcast(states[lineage == "T_NK" & state %chin% tn_states & sample_id %chin% analysis_samples],
            sample_id ~ state, value.var = "n_state", fill = 0)

setkey(meta, sample_id)
setkey(bw, sample_id)
setkey(mw, sample_id)
setkey(tw, sample_id)
if (!identical(meta$sample_id, bw$sample_id) || !identical(meta$sample_id, mw$sample_id) ||
    !identical(meta$sample_id, tw$sample_id)) stop("Wide matrices do not align to sample metadata")

broad_names <- setdiff(names(bw), "sample_id")
my_names <- setdiff(names(mw), "sample_id")
tn_names <- setdiff(names(tw), "sample_id")

xb <- close_clr(bw[, ..broad_names])
xm <- close_clr(mw[, ..my_names])
xt <- close_clr(tw[, ..tn_names])
colnames(xb) <- paste0("Broad__", broad_names)
colnames(xm) <- paste0("Myeloid__", my_names)
colnames(xt) <- paste0("T_NK__", tn_names)
x_raw <- cbind(xb, xm, xt)
rownames(x_raw) <- meta$sample_id
blocks <- c(rep("Broad", ncol(xb)), rep("Myeloid", ncol(xm)), rep("T_NK", ncol(xt)))
names(blocks) <- colnames(x_raw)

full_scaled <- block_scale(x_raw, blocks)$x
pca <- prcomp(full_scaled, center = TRUE, scale. = FALSE, rank. = 4)
scores <- as.data.table(pca$x[, seq_len(min(4L, ncol(pca$x))), drop = FALSE])
setnames(scores, paste0("Axis", seq_len(ncol(scores))))
scores <- cbind(copy(meta), scores)

variance <- pca$sdev^2 / sum(pca$sdev^2)
axis_summary <- data.table(axis = seq_len(min(4L, length(variance))),
                           explained_variance = variance[seq_len(min(4L, length(variance)))])

# Leave-one-patient-out reconstruction error. All samples from a patient are held out together.
cv_rows <- list()
patients <- unique(meta$patient_id)
for (pid in patients) {
  train <- which(meta$patient_id != pid)
  test <- which(meta$patient_id == pid)
  scaled <- block_scale(x_raw, blocks, train_rows = train)$x
  fit <- prcomp(scaled[train, , drop = FALSE], center = TRUE, scale. = FALSE, rank. = 4)
  centered_test <- sweep(scaled[test, , drop = FALSE], 2L, fit$center, "-")
  baseline <- mean(centered_test^2)
  for (k in 2:4) {
    rot <- fit$rotation[, seq_len(k), drop = FALSE]
    sc <- centered_test %*% rot
    recon_centered <- sc %*% t(rot)
    err <- mean((centered_test - recon_centered)^2)
    cv_rows[[length(cv_rows) + 1L]] <- data.table(patient_id = pid, K = k,
                                                  normalized_mse = err / baseline,
                                                  n_test_samples = length(test))
  }
}
cv_patient <- rbindlist(cv_rows)
cv_summary <- cv_patient[, .(
  mean_normalized_mse = mean(normalized_mse),
  se_normalized_mse = sd(normalized_mse) / sqrt(.N),
  n_patients = .N
), by = K]
best_row <- cv_summary[which.min(mean_normalized_mse)]
one_se_limit <- best_row$mean_normalized_mse + best_row$se_normalized_mse
selected_k <- cv_summary[mean_normalized_mse <= one_se_limit, min(K)]
cv_summary[, `:=`(one_se_limit = one_se_limit, selected = K == selected_k)]

# Patient-cluster bootstrap loading stability.
set.seed(seed)
boot_rows <- vector("list", n_boot)
reference_loading <- pca$rotation[, 1:4, drop = FALSE]
for (b in seq_len(n_boot)) {
  sampled_patients <- sample(patients, length(patients), replace = TRUE)
  idx <- unlist(lapply(sampled_patients, function(pid) which(meta$patient_id == pid)), use.names = FALSE)
  scaled_b <- block_scale(x_raw[idx, , drop = FALSE], blocks)$x
  fit_b <- try(prcomp(scaled_b, center = TRUE, scale. = FALSE, rank. = 4), silent = TRUE)
  if (inherits(fit_b, "try-error") || ncol(fit_b$rotation) < 4L) {
    boot_rows[[b]] <- data.table(iteration = b, axis = 1:4, matched_axis = NA_integer_, loading_cor = NA_real_)
  } else {
    z <- align_loading_correlations(reference_loading, fit_b$rotation[, 1:4, drop = FALSE])
    z[, iteration := b]
    boot_rows[[b]] <- z
  }
}
boot <- rbindlist(boot_rows)
stability <- boot[, .(
  successful_bootstraps = sum(is.finite(loading_cor)),
  median_loading_cor = median(loading_cor, na.rm = TRUE),
  q10_loading_cor = quantile(loading_cor, 0.10, na.rm = TRUE),
  q90_loading_cor = quantile(loading_cor, 0.90, na.rm = TRUE)
), by = axis]

# Confounding and paired anatomical-direction audit.
confound_rows <- lapply(seq_len(4L), function(a) {
  y <- scores[[paste0("Axis", a)]]
  data.table(axis = a,
             rho_log10_retained_cells = cor(y, log10(scores$retained_cells), method = "spearman"),
             cancer_accession_eta2 = eta_squared(y, scores$cancer),
             compartment_eta2 = eta_squared(y, scores$compartment))
})
confounding <- rbindlist(confound_rows)

contrast_rows <- list()
for (a in seq_len(4L)) {
  axis_col <- paste0("Axis", a)
  for (ca in unique(scores$cancer)) {
    q <- scores[cancer == ca, .(patient_id, compartment, value = get(axis_col))]
    qw <- dcast(q, patient_id ~ compartment, value.var = "value")
    for (ref in c("distal", "involved")) {
      if (all(c("tumor", ref) %in% names(qw))) {
        d <- qw[["tumor"]] - qw[[ref]]
        d <- d[is.finite(d)]
        contrast_rows[[length(contrast_rows) + 1L]] <- data.table(
          axis = a, cancer = ca, contrast = paste0("tumor_minus_", ref),
          n_pairs = length(d), median_difference = if (length(d)) median(d) else NA_real_,
          mean_difference = if (length(d)) mean(d) else NA_real_
        )
      }
    }
  }
}
contrasts <- rbindlist(contrast_rows)
direction <- contrasts[n_pairs >= 3 & is.finite(median_difference), .(
  cancers_evaluable = uniqueN(cancer),
  direction_concordant = uniqueN(sign(median_difference)) == 1L,
  direction = if (uniqueN(sign(median_difference)) == 1L) unique(sign(median_difference)) else 0
), by = .(axis, contrast)]
cross_direction <- direction[cancers_evaluable >= 2, .(cross_cancer_direction = any(direction_concordant)), by = axis]

axis_summary <- merge(axis_summary, stability, by = "axis", all.x = TRUE)
axis_summary <- merge(axis_summary, confounding, by = "axis", all.x = TRUE)
axis_summary <- merge(axis_summary, cross_direction, by = "axis", all.x = TRUE)
axis_summary[is.na(cross_cancer_direction), cross_cancer_direction := FALSE]
axis_summary[, within_selected_K := axis <= selected_k]
axis_summary[, eligible := within_selected_K &
               explained_variance >= 0.075 &
               median_loading_cor >= 0.75 &
               q10_loading_cor >= 0.50 &
               abs(rho_log10_retained_cells) < 0.50 &
               cancer_accession_eta2 < 0.50 &
               cross_cancer_direction]

stable_any <- axis_summary[within_selected_K == TRUE,
                           any(median_loading_cor >= 0.75 & q10_loading_cor >= 0.50)]
if (axis_summary[, any(eligible)]) {
  gate_status <- "PASS"
} else if (stable_any) {
  gate_status <- "CONDITIONAL"
} else {
  gate_status <- "STOP"
}

# Feature manifest and loadings.
feature_manifest <- data.table(feature = colnames(x_raw), block = blocks,
                               full_block_scale = block_scale(x_raw, blocks)$scales[blocks])
loadings <- as.data.table(pca$rotation[, 1:4, drop = FALSE], keep.rownames = "feature")
setnames(loadings, paste0("PC", 1:4), paste0("Axis", 1:4))
loadings <- melt(loadings, id.vars = "feature", variable.name = "axis_name", value.name = "loading")
loadings[, axis := as.integer(sub("Axis", "", axis_name))]
loadings <- merge(loadings, feature_manifest[, .(feature, block)], by = "feature", all.x = TRUE)

input_audit <- data.table(
  metric = c("input_samples", "analysis_samples", "excluded_low_lineage_samples",
             "patients", "cancers", "broad_features", "myeloid_features", "t_nk_features",
             "bootstrap_iterations", "selected_K", "eligible_axes", "gate_status"),
  value = c(nrow(sample_meta), nrow(meta), length(bad_samples), uniqueN(meta$patient_id),
            uniqueN(meta$cancer), length(broad_names), length(my_names), length(tn_names),
            n_boot, selected_k, sum(axis_summary$eligible), gate_status)
)

fwrite(input_audit, file.path(outdir, "input_audit.tsv"), sep = "\t")
fwrite(broad_feature_audit, file.path(outdir, "broad_feature_audit.tsv"), sep = "\t")
fwrite(state_feature_audit, file.path(outdir, "state_feature_audit.tsv"), sep = "\t")
fwrite(lineage_depth, file.path(outdir, "sample_lineage_depth.tsv"), sep = "\t")
fwrite(feature_manifest, file.path(outdir, "feature_manifest.tsv"), sep = "\t")
fwrite(scores, file.path(outdir, "ecological_axis_scores.tsv"), sep = "\t")
fwrite(loadings, file.path(outdir, "ecological_axis_loadings.tsv"), sep = "\t")
fwrite(cv_patient, file.path(outdir, "lopo_reconstruction_patient.tsv"), sep = "\t")
fwrite(cv_summary, file.path(outdir, "lopo_reconstruction_summary.tsv"), sep = "\t")
fwrite(boot, file.path(outdir, "bootstrap_loading_correlations.tsv.gz"), sep = "\t")
fwrite(axis_summary, file.path(outdir, "axis_eligibility.tsv"), sep = "\t")
fwrite(contrasts, file.path(outdir, "paired_anatomical_contrasts.tsv"), sep = "\t")
fwrite(direction, file.path(outdir, "cross_cancer_direction.tsv"), sep = "\t")
saveRDS(list(meta = meta, x_raw_clr = x_raw, blocks = blocks, pca = pca,
             selected_k = selected_k, axis_summary = axis_summary, plan_seed = seed,
             n_boot = n_boot, run_label = run_label),
        file.path(outdir, paste0("gate12g_ecological_program_", run_label, ".rds")))

# Smoke figure: visual summaries only; not used for axis selection.
plot_x <- full_scaled
plot_z <- scale(plot_x)
ord <- order(scores$Axis1)
heat <- as.data.table(plot_z[ord, , drop = FALSE], keep.rownames = "sample_id")
heat <- melt(heat, id.vars = "sample_id", variable.name = "feature", value.name = "z")
heat[, sample_id := factor(sample_id, levels = scores$sample_id[ord])]
heat[, feature := factor(feature, levels = colnames(plot_z))]
p1 <- ggplot(heat, aes(feature, sample_id, fill = pmax(-2.5, pmin(2.5, z)))) +
  geom_tile() + scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-2.5, 2.5)) +
  labs(title = "A  Compositionally transformed discovery matrix", x = NULL, y = NULL, fill = "row z") +
  theme_gate + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                     axis.text.y = element_text(size = 5))

top_load <- loadings[axis <= selected_k][order(axis, -abs(loading)), head(.SD, 8), by = axis]
top_load[, feature_label := sub("^[^_]+__", "", feature)]
top_load[, feature_label := reorder(feature_label, loading)]
p2 <- ggplot(top_load, aes(loading, feature_label, fill = block)) +
  geom_col() + facet_wrap(~axis, scales = "free_y") +
  labs(title = "B  Leading ecological-axis loadings", x = "Loading", y = NULL) + theme_gate +
  theme(legend.position = "bottom")

score_long <- melt(scores, id.vars = c("accession", "cancer", "patient_id", "sample_id",
                                       "compartment", "retained_cells"),
                   measure.vars = paste0("Axis", seq_len(selected_k)),
                   variable.name = "axis_name", value.name = "score")
p3 <- ggplot(score_long, aes(compartment, score, colour = cancer, group = interaction(cancer, patient_id))) +
  geom_line(alpha = 0.35) + geom_point(size = 1.6) + facet_wrap(~axis_name, scales = "free_y") +
  scale_x_discrete(limits = c("distal", "involved", "tumor")) +
  labs(title = "C  Paired anatomical observations", x = NULL, y = "Axis score") + theme_gate +
  theme(legend.position = "bottom")

p4 <- ggplot(boot[axis <= selected_k], aes(factor(axis), loading_cor)) +
  geom_violin(fill = "#80B1D3", colour = "#377EB8", scale = "width") +
  geom_hline(yintercept = c(0.50, 0.75), linetype = c(3, 2), colour = c("#666666", "#B2182B")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "D  Patient-cluster bootstrap stability", x = "Axis", y = "|loading correlation|") + theme_gate

fig <- (p1 | p2) / (p3 | p4) + plot_annotation(
  title = paste0("Gate12G ecological-program ", run_label, ": ", gate_status, "; selected K=", selected_k),
  subtitle = paste0(nrow(meta), " samples from ", uniqueN(meta$patient_id),
                    " patients; ", n_boot, " patient-cluster bootstraps"))
figure_stem <- paste0("FigureS_gate12g_ecological_", run_label)
ggsave(file.path(outdir, paste0(figure_stem, ".pdf")), fig, width = 14, height = 10)
png_target <- file.path(outdir, paste0(figure_stem, ".png"))
png_tmp <- tempfile(pattern = "gate12g_", fileext = ".png")
ggsave(png_tmp, fig, width = 14, height = 10, dpi = 300, bg = "white")
if (!file.copy(png_tmp, png_target, overwrite = TRUE)) stop("Failed to copy rendered PNG to output directory")
unlink(png_tmp)

eligible_text <- axis_summary[eligible == TRUE, paste0("Axis", axis, collapse = ", ")]
if (!nzchar(eligible_text)) eligible_text <- "none"
checkpoint <- c(
  paste0("# Gate12G ecological-program ", run_label, " checkpoint"),
  "",
  paste0("- Date: 2026-08-08"),
  paste0("- Frozen plan: `GATE12G_ECOLOGICAL_PROGRAM_PLAN.md`"),
  paste0("- Input: ", nrow(meta), " eligible samples from ", uniqueN(meta$patient_id), " patients"),
  paste0("- Excluded samples: ", length(bad_samples), " (lineage depth <50)"),
  paste0("- Bootstrap iterations: ", n_boot, " (", run_label, ")"),
  paste0("- LOPO one-SE selected K: ", selected_k),
  paste0("- Eligible axes: ", eligible_text),
  paste0("- Gate decision: **", gate_status, "**"),
  "",
  if (run_label == "smoke") {
    "Smoke results establish feasibility only. Axis naming, external projection and manuscript claims remain prohibited until the 1,999-bootstrap full run and frozen OEP005136 projection are complete."
  } else {
    "Full-run loadings are frozen for external projection. Final biological naming and pan-cancer claims remain prohibited until the frozen OEP005136 projection is complete."
  },
  "",
  "Cancer and accession are perfectly linked in discovery data; `cancer_accession_eta2` is therefore treated as a combined dominance diagnostic."
)
writeLines(checkpoint, file.path(outdir, paste0("GATE12G_", toupper(run_label), "_CHECKPOINT.md")))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))

cat("GATE12G_STATUS=", gate_status, "\n", sep = "")
cat("SELECTED_K=", selected_k, "\n", sep = "")
cat("ELIGIBLE_AXES=", eligible_text, "\n", sep = "")
