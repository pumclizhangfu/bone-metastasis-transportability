#!/usr/bin/env python3
"""Independent structural and numeric audit for Gate12AD-A source/provenance outputs."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import sys
from pathlib import Path


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.suffix == ".gz" else path.open("r", encoding="utf-8", newline="")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def as_bool(value: str) -> bool:
    return value.upper() == "TRUE"


def main() -> int:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    out_rel = Path(sys.argv[2] if len(sys.argv) > 2 else "results/gate12ad_figure_restructure/phase_a_source_provenance")
    out = root / out_rel
    source = out / "source_data"
    provenance = out / "provenance"
    admin = out / "admin"
    admin.mkdir(parents=True, exist_ok=True)

    checks: list[tuple[str, bool, str, str]] = []

    def record(name: str, condition: bool, observed, expected) -> None:
        checks.append((name, bool(condition), str(observed), str(expected)))

    required = [
        source / "Figure1B_umap_coordinates.tsv.gz",
        source / "Figure1B_umap_label_positions.tsv",
        source / "Figure1C_canonical_marker_dotplot.tsv",
        source / "Figure1D_all_sample_broad_composition.tsv",
        source / "Figure1E_patient_connected_broad_fractions.tsv",
        source / "Figure4C_ordered_patient_assignments.tsv",
        source / "Figure4C_transfer_consensus_matrix_all_patients.tsv",
        source / "Figure4C_ari_permutation_null.tsv.gz",
        source / "Figure5A_GSE266330_patient_scores.tsv",
        source / "Figure5B_block_lopo_predictions.tsv",
        source / "Figure5B_block_lopo_metrics.tsv",
        source / "Figure5C_hallmark_calibration_representatives.tsv",
        source / "Figure5D_paired_axis1_differences.tsv",
        source / "Figure6A_C_spot_overlay_coordinates_and_scores.tsv.gz",
        source / "Figure6D_distance_curves.tsv",
        source / "Figure6E_section_effects.tsv",
        source / "Figure6F_neighborhood_class_effects.tsv",
        provenance / "Figure4C_ari_permutation_receipt.tsv",
        provenance / "Figure6_histology_asset_index.tsv",
        provenance / "Figure6_scale_bar_audit.tsv",
        provenance / "Figure6_shared_color_scale_audit.tsv",
        provenance / "Figure6_image_processing_receipt.tsv",
        provenance / "GATE12AD_A_INPUT_MANIFEST.tsv",
        provenance / "GATE12AD_A_SOURCE_MANIFEST.tsv",
        admin / "GATE12AD_A_SOURCE_RECEIPT.json",
        out / "README.md",
    ]
    missing = [str(p.relative_to(root)) for p in required if not p.is_file()]
    record("required_files", not missing, missing, "[]")
    if missing:
        write_audit(admin, checks)
        return 1

    umap = read_tsv(source / "Figure1B_umap_coordinates.tsv.gz")
    record("figure1_umap_rows", len(umap) == 107886, len(umap), 107886)
    record("figure1_umap_unique_cells", len({r["cell_id"] for r in umap}) == 107886,
           len({r["cell_id"] for r in umap}), 107886)
    record("figure1_umap_classes", len({r["broad_class"] for r in umap}) == 11,
           len({r["broad_class"] for r in umap}), 11)

    labels = read_tsv(source / "Figure1B_umap_label_positions.tsv")
    record("figure1_umap_labels", len(labels) == 11, len(labels), 11)

    markers = read_tsv(source / "Figure1C_canonical_marker_dotplot.tsv")
    record("figure1_marker_rows", len(markers) == 143, len(markers), 143)
    record("figure1_marker_genes", len({r["gene"] for r in markers}) == 13,
           len({r["gene"] for r in markers}), 13)
    record("figure1_marker_classes", len({r["broad_class"] for r in markers}) == 11,
           len({r["broad_class"] for r in markers}), 11)
    record("figure1_marker_required_fields",
           all(k in markers[0] for k in ["mean_scaled_expression", "detected_percent", "n_samples_present", "n_patients_present"]),
           sorted(markers[0].keys()), "scaled expression, detected percent, sample/patient coverage")

    composition = read_tsv(source / "Figure1D_all_sample_broad_composition.tsv")
    samples = {r["sample_id"] for r in composition}
    record("figure1_composition_rows", len(composition) == 462, len(composition), 462)
    record("figure1_composition_samples", len(samples) == 42, len(samples), 42)
    totals = {sample: 0.0 for sample in samples}
    for row in composition:
        totals[row["sample_id"]] += float(row["fraction"])
    max_sum_error = max(abs(value - 1.0) for value in totals.values())
    record("figure1_composition_closure", max_sum_error < 1e-12, max_sum_error, "<1e-12")

    patient_fractions = read_tsv(source / "Figure1E_patient_connected_broad_fractions.tsv")
    record("figure1_patient_fraction_rows", len(patient_fractions) == 168, len(patient_fractions), 168)
    record("figure1_patient_fraction_classes", len({r["broad_class"] for r in patient_fractions}) == 4,
           len({r["broad_class"] for r in patient_fractions}), 4)

    assignments = read_tsv(source / "Figure4C_ordered_patient_assignments.tsv")
    ari_assignments = [row for row in assignments if as_bool(row["included_in_ari"])]
    record("figure4_primary_patients", len(assignments) == 49, len(assignments), 49)
    record("figure4_analysis_order", [int(r["analysis_order"]) for r in assignments] == list(range(1, 50)),
           f"{assignments[0]['analysis_order']}..{assignments[-1]['analysis_order']}", "1..49")
    record("figure4_ari_patients", len(ari_assignments) == 42, len(ari_assignments), 42)

    matrix_all = read_tsv(source / "Figure4C_transfer_consensus_matrix_all_patients.tsv")
    matrix_total = sum(int(row["N"]) for row in matrix_all)
    record("figure4_matrix_patient_total", matrix_total == 49, matrix_total, 49)

    null = read_tsv(source / "Figure4C_ari_permutation_null.tsv.gz")
    record("figure4_null_iterations", len(null) == 10000, len(null), 10000)
    receipt_rows = read_tsv(provenance / "Figure4C_ari_permutation_receipt.tsv")
    receipt = receipt_rows[0]
    observed_ari = float(receipt["observed_ari"])
    permutation_p = float(receipt["empirical_p"])
    exceedances = sum(as_bool(row["null_ge_observed"]) for row in null)
    replay_p = (1 + exceedances) / 10001
    record("figure4_observed_ari", math.isclose(observed_ari, -0.00412667085314435, abs_tol=1e-15),
           observed_ari, -0.00412667085314435)
    record("figure4_permutation_seed", int(receipt["seed"]) == 20261107, receipt["seed"], 20261107)
    record("figure4_permutation_p", math.isclose(permutation_p, 0.491650834916508, abs_tol=1e-15),
           permutation_p, 0.491650834916508)
    record("figure4_null_replays_p", math.isclose(replay_p, permutation_p, abs_tol=1e-15), replay_p, permutation_p)

    predictions = read_tsv(source / "Figure5B_block_lopo_predictions.tsv")
    models = sorted({row["model"] for row in predictions})
    model_counts = {model: sum(row["model"] == model for row in predictions) for model in models}
    record("figure5_prediction_rows", len(predictions) == 164, len(predictions), 164)
    record("figure5_prediction_models", len(models) == 4, models, "4 frozen models")
    record("figure5_samples_per_model", all(v == 41 for v in model_counts.values()), model_counts, "41 each")
    metrics = read_tsv(source / "Figure5B_block_lopo_metrics.tsv")
    record("figure5_metric_rows", len(metrics) == 4, len(metrics), 4)
    metric_r2 = {row["model"]: float(row["cross_validated_R2_training_mean"]) for row in metrics}
    expected_r2 = {
        "broad_only": 0.889906207121955,
        "myeloid_only": 0.930683064644437,
        "t_nk_only": 0.907753838000451,
        "broad_plus_depth": 0.885629594534671,
    }
    record("figure5_r2_anchors", all(math.isclose(metric_r2[k], v, abs_tol=1e-14) for k, v in expected_r2.items()),
           metric_r2, expected_r2)
    external = read_tsv(source / "Figure5A_GSE266330_patient_scores.tsv")
    record("figure5_external_people", len(external) == 47, len(external), 47)
    record("figure5_nonprojectable_retained", sum(not as_bool(r["projectable"]) for r in external) == 1,
           sum(not as_bool(r["projectable"]) for r in external), 1)
    paired = read_tsv(source / "Figure5D_paired_axis1_differences.tsv")
    record("figure5_paired_differences", len(paired) == 5, len(paired), 5)
    paired_status = {r["contrast"]: as_bool(r["endpoint_pass"]) for r in paired}
    record("figure5_paired_endpoint_asymmetry",
           paired_status.get("bm_vs_primary") is True and paired_status.get("bm_vs_normal_bone") is False,
           paired_status, "primary PASS; normal bone FAIL")

    overlays = read_tsv(source / "Figure6A_C_spot_overlay_coordinates_and_scores.tsv.gz")
    record("figure6_spot_rows", len(overlays) == 8190, len(overlays), 8190)
    record("figure6_sections", len({r["sample"] for r in overlays}) == 4,
           len({r["sample"] for r in overlays}), 4)
    record("figure6_unreachable_retained", any(r["distance_ring"] == "unreachable" for r in overlays),
           sum(r["distance_ring"] == "unreachable" for r in overlays), ">0")

    assets = read_tsv(provenance / "Figure6_histology_asset_index.tsv")
    record("figure6_asset_rows", len(assets) == 4, len(assets), 4)
    asset_hash_pass = True
    for row in assets:
        image = root / row["image_relative_path"]
        scale_json = root / row["scalefactors_relative_path"]
        asset_hash_pass &= image.is_file() and scale_json.is_file()
        asset_hash_pass &= sha256(image) == row["image_sha256"] and sha256(scale_json) == row["scalefactors_sha256"]
    record("figure6_asset_hashes", asset_hash_pass, asset_hash_pass, True)

    scale_bars = read_tsv(provenance / "Figure6_scale_bar_audit.tsv")
    record("figure6_scale_bar_rows", len(scale_bars) == 4, len(scale_bars), 4)
    record("figure6_scale_bar_bounds", all(as_bool(r["coordinate_bounds_pass"]) for r in scale_bars),
           [r["coordinate_bounds_pass"] for r in scale_bars], "all TRUE")
    bar_values = [float(r["pixels_per_500um_lowres"]) for r in scale_bars]
    record("figure6_scale_bar_numeric_range", min(bar_values) > 38.0 and max(bar_values) < 41.0,
           [min(bar_values), max(bar_values)], "38-41 low-resolution pixels")

    scales = read_tsv(provenance / "Figure6_shared_color_scale_audit.tsv")
    shared_limits = {float(r["shared_abs_q98_limit"]) for r in scales}
    clipped = {r["layer"]: int(r["n_clipped_at_shared_limit"]) for r in scales}
    record("figure6_one_shared_limit", len(shared_limits) == 1, shared_limits, "one exact value")
    shared_limit = next(iter(shared_limits))
    record("figure6_shared_limit_anchor", math.isclose(shared_limit, 0.329220075, abs_tol=5e-10),
           shared_limit, 0.329220075)
    record("figure6_clip_counts", clipped == {"Full Axis1": 156, "Malignant-excluded Axis1": 172},
           clipped, {"Full Axis1": 156, "Malignant-excluded Axis1": 172})
    image_processing = read_tsv(provenance / "Figure6_image_processing_receipt.tsv")[0]
    record("figure6_lightening_disclosed", math.isclose(float(image_processing["lightening_fraction"]), 0.12),
           image_processing["lightening_fraction"], 0.12)
    record("figure6_no_local_adjustment",
           not as_bool(image_processing["local_adjustment_applied"]) and not as_bool(image_processing["local_masking_applied"]),
           (image_processing["local_adjustment_applied"], image_processing["local_masking_applied"]), "FALSE, FALSE")

    for manifest_name in ["GATE12AD_A_INPUT_MANIFEST.tsv", "GATE12AD_A_SOURCE_MANIFEST.tsv"]:
        manifest_path = provenance / manifest_name
        manifest = read_tsv(manifest_path)
        failures = []
        for row in manifest:
            target = root / row["relative_path"]
            if not target.is_file() or target.stat().st_size != int(float(row["bytes"])) or sha256(target) != row["sha256"]:
                failures.append(row["relative_path"])
        record(f"manifest_verify:{manifest_name}", not failures, failures, "[]")

    with (admin / "GATE12AD_A_SOURCE_RECEIPT.json").open("r", encoding="utf-8") as handle:
        source_receipt = json.load(handle)
    record("receipt_parent_immutability", source_receipt.get("frozen_inputs_modified") is False,
           source_receipt.get("frozen_inputs_modified"), False)

    failed = write_audit(admin, checks)
    print(f"GATE12AD_A_AUDIT_STATUS={'PASS' if not failed else 'FAIL'}")
    print(f"CHECKS_PASSED={len(checks) - len(failed)}/{len(checks)}")
    if failed:
        print("FAILED_CHECKS=" + ",".join(failed))
        return 1
    return 0


def write_audit(admin: Path, checks: list[tuple[str, bool, str, str]]) -> list[str]:
    failed = [name for name, ok, _, _ in checks if not ok]
    with (admin / "GATE12AD_A_AUDIT.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["check", "status", "observed", "expected"])
        for name, ok, observed, expected in checks:
            writer.writerow([name, "PASS" if ok else "FAIL", observed, expected])
    lines = [
        "# Gate12AD-A source/provenance audit",
        "",
        f"- Status: **{'PASS' if not failed else 'FAIL'}**",
        f"- Checks passed: {len(checks) - len(failed)}/{len(checks)}",
        f"- Failed checks: {len(failed)}",
        "- Scope: Figure 1 panel-level source data, Figure 4 exact permutation replay, Figure 5 patient-held-out predictions, Figure 6 image/scale/clipping provenance, and input/output hash manifests.",
        "",
    ]
    if failed:
        lines.extend(["## Failed checks", "", *[f"- {name}" for name in failed]])
    else:
        lines.extend([
            "## Decision",
            "",
            "Gate12AD-A satisfies the frozen source/provenance contract and may proceed to the six-main-figure rebuild. No frozen Gate12AB input was changed.",
        ])
    (admin / "GATE12AD_A_AUDIT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return failed


if __name__ == "__main__":
    raise SystemExit(main())
