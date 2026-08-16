#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: run_gate12e_spatial_axis_analysis.R SECTION_RDS OUT_DIR [N_PERM]")
}
section_rds <- normalizePath(args[[1]], mustWork = TRUE)
out_dir <- normalizePath(args[[2]], mustWork = FALSE)
n_perm <- if (length(args) >= 3L) as.integer(args[[3]]) else 999L
if (!is.finite(n_perm) || n_perm < 19L) stop("N_PERM must be an integer >= 19")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

base_seed <- 20260808L
knn_k <- 6L
alpha_q <- 0.10

axes <- data.table(
  axis_id = c("A1", "A2", "A3"),
  axis = c("CCL4 -> CCR5", "CXCL16 -> CXCR6", "SPP1 -> CD44"),
  human_ligand = c("CCL4", "CXCL16", "SPP1"),
  human_receptor = c("CCR5", "CXCR6", "CD44"),
  ligand = c("Ccl4", "Cxcl16", "Spp1"),
  receptor = c("Ccr5", "Cxcr6", "Cd44")
)

write_tsv <- function(x, path) {
  fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
}

build_symmetric_knn <- function(coords, k = 6L) {
  coords <- as.matrix(coords[, c("x", "y"), drop = FALSE])
  if (nrow(coords) <= k) stop("Not enough spots for ", k, "-NN graph")
  distance <- as.matrix(dist(coords))
  diag(distance) <- Inf
  neighbours <- t(vapply(
    seq_len(nrow(distance)),
    function(i) order(distance[i, ], method = "radix")[seq_len(k)],
    integer(k)
  ))
  directed <- sparseMatrix(
    i = rep(seq_len(nrow(coords)), each = k),
    j = as.vector(t(neighbours)),
    x = 1,
    dims = c(nrow(coords), nrow(coords))
  )
  adjacency <- as((directed + t(directed)) > 0, "dMatrix")
  diag(adjacency) <- 0
  adjacency <- drop0(adjacency)
  degree <- as.numeric(rowSums(adjacency))
  if (any(degree == 0)) stop("Symmetric k-NN graph contains isolated spots")
  list(adjacency = adjacency, degree = degree, undirected_edges = sum(adjacency) / 2)
}

moran_permutation <- function(x, adjacency, permutations, seed) {
  x <- as.numeric(x)
  n <- length(x)
  centered <- x - mean(x)
  denominator <- sum(centered^2)
  s0 <- sum(adjacency)
  if (!is.finite(denominator) || denominator <= 0 || s0 <= 0) {
    return(c(estimate = NA_real_, p_one_sided = NA_real_))
  }
  observed <- (n / s0) * sum(centered * as.numeric(adjacency %*% centered)) / denominator
  set.seed(seed)
  permuted <- replicate(permutations, sample(centered, replace = FALSE))
  spatial_lag <- adjacency %*% permuted
  null <- (n / s0) * colSums(permuted * spatial_lag) / denominator
  p_value <- (1 + sum(null >= observed, na.rm = TRUE)) / (permutations + 1)
  c(estimate = observed, p_one_sided = p_value)
}

neighbour_colocalization <- function(sender, receiver, adjacency, degree, permutations, seed) {
  sender <- as.numeric(sender)
  receiver <- as.numeric(receiver)
  neighbour_receiver <- as.numeric(adjacency %*% receiver) / degree
  observed <- suppressWarnings(cor(sender, neighbour_receiver, method = "spearman"))
  if (!is.finite(observed)) {
    return(list(rho = NA_real_, p_one_sided = NA_real_, neighbour_receiver = neighbour_receiver))
  }
  set.seed(seed)
  receiver_permuted <- replicate(permutations, sample(receiver, replace = FALSE))
  neighbour_permuted <- adjacency %*% receiver_permuted
  neighbour_permuted <- sweep(neighbour_permuted, 1L, degree, "/")
  sender_rank <- rank(sender, ties.method = "average")
  null <- apply(neighbour_permuted, 2L, function(z) {
    suppressWarnings(cor(sender_rank, rank(z, ties.method = "average"), method = "pearson"))
  })
  p_value <- (1 + sum(null >= observed, na.rm = TRUE)) / (permutations + 1)
  list(rho = observed, p_one_sided = p_value, neighbour_receiver = neighbour_receiver)
}

rank_biserial <- function(values, group) {
  values <- as.numeric(values)
  group <- as.logical(group)
  keep <- is.finite(values) & !is.na(group)
  values <- values[keep]
  group <- group[keep]
  n1 <- sum(group)
  n0 <- sum(!group)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  ranks <- rank(values, ties.method = "average")
  u1 <- sum(ranks[group]) - n1 * (n1 + 1) / 2
  2 * u1 / (n1 * n0) - 1
}

lexicographic_order <- function(dt) {
  order(
    -dt$coloc_sig_sections,
    -dt$min_moran_sig_sections,
    -dt$median_coloc_rho,
    -dt$median_tumor_adjacency_rb,
    -dt$gate12c_robustness_wins,
    dt$axis_id,
    na.last = TRUE
  )
}

sections <- readRDS(section_rds)
if (length(sections) != 4L) stop("Expected four spatial sections")
required_section_fields <- c(
  "coords", "submitted_label", "weights", "candidate_expression",
  "candidate_detection"
)
for (sample in names(sections)) {
  missing_fields <- setdiff(required_section_fields, names(sections[[sample]]))
  if (length(missing_fields)) stop(sample, ": missing fields: ", paste(missing_fields, collapse = ", "))
}

reconstruction_dir <- dirname(section_rds)
reference_detection_path <- file.path(reconstruction_dir, "scrna_candidate_gene_detection.tsv")
if (!file.exists(reference_detection_path)) stop("Missing reference detection receipt")
reference_detection <- fread(reference_detection_path)

reference_support <- axes[, {
  ligand_fraction <- reference_detection[
    gene == .BY$ligand & cell_type == "Macrophage", detection_fraction
  ]
  receptor_fraction <- reference_detection[
    gene == .BY$receptor & cell_type == "T cell", detection_fraction
  ]
  .(
    ligand_macrophage_detection = if (length(ligand_fraction)) ligand_fraction[[1L]] else NA_real_,
    receptor_tcell_detection = if (length(receptor_fraction)) receptor_fraction[[1L]] else NA_real_
  )
}, by = .(axis_id, axis, ligand, receptor)]
reference_support[, `:=`(
  ligand_reference_pass = ligand_macrophage_detection >= 0.05,
  receptor_reference_pass = receptor_tcell_detection >= 0.05
)]
reference_support[, reference_support_pass := ligand_reference_pass & receptor_reference_pass]

graph_receipt <- list()
detection_rows <- list()
moran_rows <- list()
coloc_rows <- list()
adjacency_rows <- list()
plot_rows <- list()

for (section_index in seq_along(sections)) {
  sample <- names(sections)[[section_index]]
  section <- sections[[sample]]
  spot_ids <- rownames(section$weights)
  if (!identical(spot_ids, rownames(section$coords))) stop(sample, ": weight/coordinate order mismatch")
  if (!identical(spot_ids, colnames(section$candidate_expression))) stop(sample, ": weight/expression order mismatch")
  if (length(section$submitted_label) != length(spot_ids)) stop(sample, ": label length mismatch")

  graph <- build_symmetric_knn(section$coords, knn_k)
  adjacency <- graph$adjacency
  degree <- graph$degree
  graph_receipt[[sample]] <- data.table(
    sample = sample,
    spots = length(spot_ids),
    k = knn_k,
    symmetric_edges = graph$undirected_edges,
    min_degree = min(degree),
    median_degree = median(degree),
    max_degree = max(degree)
  )

  detection_matrix <- section$candidate_detection
  for (gene_index in seq_along(unique(c(axes$ligand, axes$receptor)))) {
    gene <- unique(c(axes$ligand, axes$receptor))[[gene_index]]
    detected <- if (gene %in% rownames(detection_matrix)) {
      as.logical(detection_matrix[gene, , drop = TRUE])
    } else {
      rep(FALSE, length(spot_ids))
    }
    expression <- if (gene %in% rownames(section$candidate_expression)) {
      as.numeric(section$candidate_expression[gene, ])
    } else {
      rep(0, length(spot_ids))
    }
    moran <- moran_permutation(
      expression,
      adjacency,
      n_perm,
      base_seed + section_index * 1000L + gene_index
    )
    detection_rows[[length(detection_rows) + 1L]] <- data.table(
      sample = sample,
      gene = gene,
      spots = length(spot_ids),
      detected_spots = sum(detected),
      detection_fraction = mean(detected),
      availability_pass = mean(detected) >= 0.02
    )
    gene_axes <- axes[ligand == gene | receptor == gene]
    moran_rows[[length(moran_rows) + 1L]] <- data.table(
      sample = sample,
      axis_id = gene_axes$axis_id[[1L]],
      axis = gene_axes$axis[[1L]],
      role = if (gene_axes$ligand[[1L]] == gene) "ligand" else "receptor",
      gene = gene,
      moran_i = unname(moran[["estimate"]]),
      p_one_sided = unname(moran[["p_one_sided"]]),
      permutations = n_perm
    )
  }

  labels <- unname(section$submitted_label[spot_ids])
  if (anyNA(labels)) labels <- unname(section$submitted_label)
  tumor <- labels == "Tumor"
  tumor_adjacent <- !tumor & as.numeric(adjacency %*% as.numeric(tumor)) > 0
  other_non_tumor <- !tumor & !tumor_adjacent

  for (axis_index in seq_len(nrow(axes))) {
    axis_row <- axes[axis_index]
    ligand_expression <- as.numeric(section$candidate_expression[axis_row$ligand, ])
    receptor_expression <- as.numeric(section$candidate_expression[axis_row$receptor, ])
    sender <- ligand_expression * section$weights[, "Macrophage"]
    receiver <- receptor_expression * section$weights[, "T cell"]
    coloc <- neighbour_colocalization(
      sender,
      receiver,
      adjacency,
      degree,
      n_perm,
      base_seed + 100000L + section_index * 1000L + axis_index
    )
    colocalization_signal <- sender * coloc$neighbour_receiver
    rb <- rank_biserial(colocalization_signal[!tumor], tumor_adjacent[!tumor])
    coloc_rows[[length(coloc_rows) + 1L]] <- data.table(
      sample = sample,
      axis_id = axis_row$axis_id,
      axis = axis_row$axis,
      ligand = axis_row$ligand,
      receptor = axis_row$receptor,
      spearman_rho = coloc$rho,
      p_one_sided = coloc$p_one_sided,
      permutations = n_perm
    )
    adjacency_rows[[length(adjacency_rows) + 1L]] <- data.table(
      sample = sample,
      axis_id = axis_row$axis_id,
      axis = axis_row$axis,
      tumor_spots = sum(tumor),
      tumor_adjacent_non_tumor_spots = sum(tumor_adjacent),
      other_non_tumor_spots = sum(other_non_tumor),
      rank_biserial_tumor_adjacent_vs_other = rb,
      adjacent_median_colocalization = median(colocalization_signal[tumor_adjacent], na.rm = TRUE),
      other_median_colocalization = median(colocalization_signal[other_non_tumor], na.rm = TRUE)
    )
  }

  plot_rows[[sample]] <- data.table(
    sample = sample,
    x = section$coords$x,
    y = section$coords$y,
    submitted_label = labels,
    macrophage_weight = section$weights[, "Macrophage"],
    tcell_weight = section$weights[, "T cell"]
  )
}

graph_receipt <- rbindlist(graph_receipt)
spatial_detection <- rbindlist(detection_rows)
moran <- rbindlist(moran_rows)
coloc <- rbindlist(coloc_rows)
tumor_adjacency <- rbindlist(adjacency_rows)
plot_data <- rbindlist(plot_rows)

# Frozen multiplicity rule: adjust the three candidate-axis tests within each
# section and statistical family. Ligand and receptor Moran tests are separate
# families; section-level P values are never pooled across sections.
moran[, q_bh := p.adjust(p_one_sided, method = "BH"), by = .(sample, role)]
moran[, positive_q_pass := moran_i > 0 & q_bh < alpha_q]
coloc[, q_bh := p.adjust(p_one_sided, method = "BH"), by = sample]
coloc[, positive_q_pass := spearman_rho > 0 & q_bh < alpha_q]

availability_summary <- axes[, {
  ligand_rows <- spatial_detection[gene == .BY$ligand]
  receptor_rows <- spatial_detection[gene == .BY$receptor]
  .(
    ligand_available_sections = sum(ligand_rows$availability_pass),
    receptor_available_sections = sum(receptor_rows$availability_pass)
  )
}, by = .(axis_id, axis, ligand, receptor)]
availability_summary[, spatial_availability_pass :=
  ligand_available_sections >= 3L & receptor_available_sections >= 3L]

moran_summary <- axes[, {
  ligand_rows <- moran[axis_id == .BY$axis_id & role == "ligand"]
  receptor_rows <- moran[axis_id == .BY$axis_id & role == "receptor"]
  .(
    ligand_moran_sig_sections = sum(ligand_rows$positive_q_pass),
    receptor_moran_sig_sections = sum(receptor_rows$positive_q_pass)
  )
}, by = .(axis_id, axis)]
moran_summary[, `:=`(
  min_moran_sig_sections = pmin(ligand_moran_sig_sections, receptor_moran_sig_sections),
  spatial_autocorrelation_pass = ligand_moran_sig_sections >= 3L & receptor_moran_sig_sections >= 3L
)]

coloc_summary <- coloc[, .(
  coloc_sig_sections = sum(positive_q_pass),
  median_coloc_rho = median(spearman_rho, na.rm = TRUE)
), by = .(axis_id, axis)]
coloc_summary[, neighborhood_colocalization_pass := coloc_sig_sections >= 3L]

adjacency_summary <- tumor_adjacency[, .(
  median_tumor_adjacency_rb = median(rank_biserial_tumor_adjacent_vs_other, na.rm = TRUE)
), by = .(axis_id, axis)]

gate12c_path <- file.path(
  dirname(dirname(reconstruction_dir)),
  "gate12c_external_validation",
  "external_threshold_sensitivity.tsv"
)
if (!file.exists(gate12c_path)) {
  gate12c_path <- file.path(
    dirname(reconstruction_dir), "gate12c_external_validation", "external_threshold_sensitivity.tsv"
  )
}
gate12c_wins <- data.table(axis = axes$axis, gate12c_robustness_wins = 0L)
if (file.exists(gate12c_path)) {
  sensitivity <- unique(fread(gate12c_path)[, .(min_cells, min_detect, scenario_winner)])
  wins <- sensitivity[, .(gate12c_robustness_wins = .N), by = .(axis = scenario_winner)]
  gate12c_wins <- merge(gate12c_wins[, .(axis)], wins, by = "axis", all.x = TRUE)
  gate12c_wins[is.na(gate12c_robustness_wins), gate12c_robustness_wins := 0L]
}

axis_decision <- Reduce(
  function(x, y) merge(x, y, by = intersect(names(x), names(y)), all = TRUE),
  list(
    reference_support,
    availability_summary,
    moran_summary,
    coloc_summary,
    adjacency_summary,
    gate12c_wins
  )
)
axis_decision[, spatially_eligible :=
  reference_support_pass & spatial_availability_pass &
  spatial_autocorrelation_pass & neighborhood_colocalization_pass]
axis_decision[, spatial_rank := NA_integer_]
if (any(axis_decision$spatially_eligible)) {
  eligible_index <- which(axis_decision$spatially_eligible)
  ranked <- eligible_index[lexicographic_order(axis_decision[eligible_index])]
  axis_decision[ranked, spatial_rank := seq_along(ranked)]
}
setorder(axis_decision, spatial_rank, axis_id, na.last = TRUE)

eligible_axes <- axis_decision[spatially_eligible == TRUE]
loo_rows <- list()
for (omitted in names(sections)) {
  remaining <- setdiff(names(sections), omitted)
  for (axis_index in seq_len(nrow(eligible_axes))) {
    axis_row <- eligible_axes[axis_index]
    moran_sub <- moran[sample %in% remaining & axis_id == axis_row$axis_id]
    coloc_sub <- coloc[sample %in% remaining & axis_id == axis_row$axis_id]
    adj_sub <- tumor_adjacency[sample %in% remaining & axis_id == axis_row$axis_id]
    loo_rows[[length(loo_rows) + 1L]] <- data.table(
      omitted_section = omitted,
      axis_id = axis_row$axis_id,
      axis = axis_row$axis,
      coloc_sig_sections = sum(coloc_sub$positive_q_pass),
      min_moran_sig_sections = min(
        sum(moran_sub[role == "ligand"]$positive_q_pass),
        sum(moran_sub[role == "receptor"]$positive_q_pass)
      ),
      median_coloc_rho = median(coloc_sub$spearman_rho, na.rm = TRUE),
      median_tumor_adjacency_rb = median(
        adj_sub$rank_biserial_tumor_adjacent_vs_other,
        na.rm = TRUE
      ),
      gate12c_robustness_wins = axis_row$gate12c_robustness_wins
    )
  }
}
leave_one_out <- rbindlist(loo_rows, fill = TRUE)
if (nrow(leave_one_out)) {
  leave_one_out[, loo_rank := {
    local_order <- lexicographic_order(.SD)
    result <- integer(.N)
    result[local_order] <- seq_along(local_order)
    result
  }, by = omitted_section]
  leave_one_out[, loo_winner := loo_rank == 1L]
}

global_winner <- axis_decision[spatial_rank == 1L]
freeze_status <- "BLOCKED_NO_SPATIALLY_ELIGIBLE_AXIS"
frozen_axis <- NA_character_
frozen_mouse_receptor <- NA_character_
frozen_human_receptor <- NA_character_
loo_stable <- FALSE
if (nrow(global_winner) == 1L) {
  frozen_axis <- global_winner$axis[[1L]]
  frozen_mouse_receptor <- axes[axis == frozen_axis, receptor][[1L]]
  frozen_human_receptor <- axes[axis == frozen_axis, human_receptor][[1L]]
  loo_stable <- nrow(leave_one_out) > 0L && all(
    leave_one_out[loo_rank == 1L, axis] == frozen_axis
  ) && length(unique(leave_one_out$omitted_section)) == length(sections)
  freeze_status <- if (loo_stable) {
    "PASS_TARGET_FROZEN"
  } else {
    "BLOCKED_LEAVE_ONE_SECTION_OUT_INSTABILITY"
  }
}

target_freeze_audit <- data.table(
  permutations = n_perm,
  sections = length(sections),
  independent_animal_inference = FALSE,
  eligible_axes = nrow(eligible_axes),
  global_top_axis = frozen_axis,
  mouse_receptor_evidence_gene = frozen_mouse_receptor,
  human_virtual_knockout_target = frozen_human_receptor,
  leave_one_section_out_stable = loo_stable,
  final_status = freeze_status
)

write_tsv(graph_receipt, file.path(out_dir, "spatial_graph_receipt.tsv"))
write_tsv(reference_support, file.path(out_dir, "reference_axis_support.tsv"))
write_tsv(spatial_detection, file.path(out_dir, "spatial_gene_detection.tsv"))
write_tsv(moran, file.path(out_dir, "spatial_moran.tsv"))
write_tsv(coloc, file.path(out_dir, "neighborhood_colocalization.tsv"))
write_tsv(tumor_adjacency, file.path(out_dir, "tumor_adjacency.tsv"))
write_tsv(axis_decision, file.path(out_dir, "axis_decision.tsv"))
write_tsv(leave_one_out, file.path(out_dir, "leave_one_section_out.tsv"))
write_tsv(target_freeze_audit, file.path(out_dir, "target_freeze_audit.tsv"))

saveRDS(
  list(
    seed = base_seed,
    permutations = n_perm,
    knn_k = knn_k,
    graph_receipt = graph_receipt,
    reference_support = reference_support,
    spatial_detection = spatial_detection,
    moran = moran,
    colocalization = coloc,
    tumor_adjacency = tumor_adjacency,
    axis_decision = axis_decision,
    leave_one_section_out = leave_one_out,
    target_freeze_audit = target_freeze_audit
  ),
  file.path(out_dir, "gate12e_spatial_axis_results.rds"),
  compress = "xz"
)

theme_gate12 <- theme_bw(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "#ECEFF4", colour = NA),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 11),
    legend.title = element_text(face = "bold")
  )

p_tissue <- ggplot(plot_data, aes(x, -y, colour = submitted_label)) +
  geom_point(size = 0.38, alpha = 0.9) +
  facet_wrap(~sample, nrow = 1) +
  scale_colour_manual(values = c(
    "Tumor" = "#D73027", "MSCs" = "#FC8D59", "Others" = "#91BFDB",
    "Distantial Tissue" = "#BDBDBD"
  )) +
  coord_equal() +
  labs(title = "A  Submitted tissue labels", x = NULL, y = NULL, colour = "Label") +
  theme_gate12 + theme(axis.text = element_blank(), axis.ticks = element_blank())

p_macrophage <- ggplot(plot_data, aes(x, -y, colour = macrophage_weight)) +
  geom_point(size = 0.42) + facet_wrap(~sample, nrow = 1) + coord_equal() +
  scale_colour_viridis_c(option = "magma", limits = c(0, 1)) +
  labs(title = "B  RCTD macrophage weight", x = NULL, y = NULL, colour = "Weight") +
  theme_gate12 + theme(axis.text = element_blank(), axis.ticks = element_blank())

p_tcell <- ggplot(plot_data, aes(x, -y, colour = tcell_weight)) +
  geom_point(size = 0.42) + facet_wrap(~sample, nrow = 1) + coord_equal() +
  scale_colour_viridis_c(option = "viridis", limits = c(0, 1)) +
  labs(title = "C  RCTD T-cell weight", x = NULL, y = NULL, colour = "Weight") +
  theme_gate12 + theme(axis.text = element_blank(), axis.ticks = element_blank())

detection_plot <- merge(
  spatial_detection,
  melt(axes, id.vars = c("axis_id", "axis"), measure.vars = c("ligand", "receptor"),
       variable.name = "role", value.name = "gene"),
  by = "gene"
)
p_detection <- ggplot(detection_plot, aes(sample, gene, size = detection_fraction,
                                          colour = availability_pass)) +
  geom_point() +
  scale_size_continuous(range = c(1, 5), labels = scales::percent_format(accuracy = 1)) +
  scale_colour_manual(values = c(`TRUE` = "#1B9E77", `FALSE` = "#D95F02")) +
  labs(title = "D  Candidate-gene spatial availability", x = NULL, y = NULL,
       size = "Detected", colour = ">=2%") + theme_gate12 +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

p_moran <- ggplot(moran, aes(sample, paste0(gene, " (", role, ")"), fill = moran_i)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = ifelse(positive_q_pass, "*", "")), size = 4) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(title = "E  Spatial autocorrelation", subtitle = "* positive Moran's I, BH q<0.10",
       x = NULL, y = NULL, fill = "Moran's I") + theme_gate12 +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

p_coloc <- ggplot(coloc, aes(sample, axis, fill = spearman_rho)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = ifelse(positive_q_pass, "*", "")), size = 4) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(title = "F  Macrophage-to-T-cell neighbourhood colocalization",
       subtitle = "* positive Spearman rho, BH q<0.10", x = NULL, y = NULL,
       fill = "Spearman rho") + theme_gate12 +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

decision_plot <- melt(
  axis_decision,
  id.vars = c("axis_id", "axis", "spatially_eligible"),
  measure.vars = c("coloc_sig_sections", "ligand_moran_sig_sections", "receptor_moran_sig_sections"),
  variable.name = "evidence", value.name = "supporting_sections"
)
decision_plot[, evidence := factor(evidence,
  levels = c("coloc_sig_sections", "ligand_moran_sig_sections", "receptor_moran_sig_sections"),
  labels = c("Neighbour colocalization", "Ligand Moran", "Receptor Moran")
)]
p_decision <- ggplot(decision_plot, aes(axis, supporting_sections, colour = evidence, group = evidence)) +
  geom_hline(yintercept = 3, linetype = 2, colour = "grey55") +
  geom_point(size = 2.5) + geom_line(linewidth = 0.5) +
  scale_y_continuous(breaks = 0:4, limits = c(0, 4.15)) +
  labs(title = "G  Frozen decision criteria", x = NULL, y = "Supporting sections",
       colour = "Evidence family") + theme_gate12 +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "bottom")

figure <- p_tissue / p_macrophage / p_tcell / ((p_detection | p_moran) / (p_coloc | p_decision)) +
  plot_annotation(
    title = "Gate12E cross-species spatial support in early mouse bone metastasis",
    subtitle = paste0("Four tissue sections; symmetric 6-NN graph; ", n_perm,
                      " within-section permutations; sections are not independent animals"),
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(size = 10))
  )
png_path <- file.path(out_dir, "Figure6_gate12e_spatial_axis.png")
png_tmp <- tempfile(pattern = "gate12e_figure6_", fileext = ".png")
ggsave(png_tmp, figure, width = 15, height = 13, dpi = 300,
       device = ragg::agg_png, bg = "white")
if (!file.copy(png_tmp, png_path, overwrite = TRUE)) {
  stop("Failed to copy rendered Figure 6 PNG into the output directory")
}
unlink(png_tmp)
ggsave(file.path(out_dir, "Figure6_gate12e_spatial_axis.pdf"), figure,
       width = 15, height = 13, device = cairo_pdf, bg = "white")

capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))

cat("GATE12E_SPATIAL_AXIS_ANALYSIS=PASS\n")
cat("PERMUTATIONS=", n_perm, "\n", sep = "")
cat("ELIGIBLE_AXES=", nrow(eligible_axes), "\n", sep = "")
cat("GLOBAL_TOP_AXIS=", ifelse(is.na(frozen_axis), "NONE", frozen_axis), "\n", sep = "")
cat("FINAL_STATUS=", freeze_status, "\n", sep = "")
