"""Post-hoc diagnostics for the frozen Kuramoto audit winner.

The diagnostics in this module are deliberately downstream of program
selection.  They interpret the already-frozen Kuramoto result; they do not
participate in synthesis and do not upgrade it to prospective evidence.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from itertools import combinations
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


from .canonical import normalize_trajectory, permutation_entropy
from .datasets import Trajectory, kuramoto_domain
from .programs import Program
from .statistics import (
    fixed_squash,
    independent_cross_label_pairs,
    orientation_free_auc,
    stable_seed,
)
from .synthesis import representation


EPS = 1.0e-9
AUDIT_SEED = 20260822 + 193
FIRST_ORDER_BOOTSTRAP_SEED = 9922
SECOND_ORDER_BOOTSTRAP_SEED = 9923
COMPONENT_DIFFERENCE_BOOTSTRAP_SEED = 9924
PAIRING_SEED = 555
FROZEN_KURAMOTO_PROGRAM = "log_ratio(increment_cv,radial_permutation_entropy)"


@dataclass(frozen=True)
class KuramotoHiddenRecord:
    trajectory: Trajectory
    coupling: float
    mean_order: float
    mean_second_order: float
    std_real_second_order: float


def selected_scalar(values: np.ndarray) -> tuple[float, float, float]:
    """Evaluate the frozen log-ratio after subset-specific canonicalization."""

    x = normalize_trajectory(values)
    radial = np.linalg.norm(x, axis=1)
    increments = np.linalg.norm(np.diff(x, axis=0), axis=1)
    increment_cv = float(np.std(increments) / (np.mean(increments) + 1.0e-12))
    radial_pe = float(permutation_entropy(radial))
    value = math.log((abs(increment_cv) + EPS) / (abs(radial_pe) + EPS))
    return value, increment_cv, radial_pe


def sine_norm_from_second_order(theta: np.ndarray) -> tuple[float, float]:
    """Return both sides of the noiseless sine-sensor / Z_2 identity."""

    phase = np.asarray(theta, dtype=np.float64)
    left = float(np.sum(np.sin(phase) ** 2))
    z2 = np.mean(np.exp(2j * phase))
    right = float(0.5 * len(phase) * (1.0 - np.real(z2)))
    return left, right


def centered_sine_norm_from_second_order(
    theta: np.ndarray,
    channel_means: np.ndarray,
) -> tuple[float, float]:
    """Return both sides of the exact channel-centering correction."""

    phase = np.asarray(theta, dtype=np.float64)
    means = np.asarray(channel_means, dtype=np.float64)
    if phase.shape != means.shape:
        raise ValueError("theta and channel_means must have the same shape")
    sine = np.sin(phase)
    raw_left, raw_right = sine_norm_from_second_order(phase)
    if not np.isclose(raw_left, raw_right, rtol=0.0, atol=1.0e-12):
        raise AssertionError("raw second-harmonic identity failed")
    left = float(np.sum((sine - means) ** 2))
    right = float(raw_right - 2.0 * np.sum(means * sine) + np.sum(means**2))
    return left, right


def kuramoto_with_hidden_order(
    n_per_class: int = 120,
    seed: int = AUDIT_SEED,
) -> list[KuramotoHiddenRecord]:
    """Mirror the released generator while retaining its unobserved phases."""

    rng = np.random.default_rng(seed)
    result: list[KuramotoHiddenRecord] = []
    oscillators = 5
    for label, coupling_interval in ((0, (2.8, 4.0)), (1, (0.05, 0.45))):
        for index in range(n_per_class):
            coupling = float(rng.uniform(*coupling_interval))
            omega = rng.normal(1.0, 0.25, size=oscillators)
            theta = rng.uniform(-np.pi, np.pi, size=oscillators)
            kept_values: list[np.ndarray] = []
            kept_order: list[float] = []
            kept_second_order: list[float] = []
            kept_real_second_order: list[float] = []
            for step in range(2300):
                differences = theta[None, :] - theta[:, None]
                theta += 0.025 * (
                    omega + coupling * np.mean(np.sin(differences), axis=1)
                )
                if step >= 1020 and step % 5 == 0:
                    kept_values.append(np.sin(theta.copy()))
                    kept_order.append(float(abs(np.mean(np.exp(1j * theta)))))
                    second_order = np.mean(np.exp(2j * theta))
                    kept_second_order.append(float(abs(second_order)))
                    kept_real_second_order.append(float(np.real(second_order)))
            clean = np.asarray(kept_values[:256])
            values = clean + rng.normal(0.0, 0.015, size=clean.shape)
            trajectory = Trajectory(
                f"audit-kuramoto-{label}-{index:03d}",
                "kuramoto",
                label,
                values,
            )
            result.append(
                KuramotoHiddenRecord(
                    trajectory=trajectory,
                    coupling=coupling,
                    mean_order=float(np.mean(kept_order[:256])),
                    mean_second_order=float(np.mean(kept_second_order[:256])),
                    std_real_second_order=float(np.std(kept_real_second_order[:256])),
                )
            )
    return result


def exact_regeneration_error(records: list[KuramotoHiddenRecord]) -> float:
    released = [
        replace(item, identifier=f"audit-{item.identifier}")
        for item in kuramoto_domain(n_per_class=len(records) // 2, seed=AUDIT_SEED)
    ]
    if [item.trajectory.identifier for item in records] != [item.identifier for item in released]:
        raise AssertionError("post-hoc identifiers do not match the frozen audit")
    return float(
        max(
            np.max(np.abs(hidden.trajectory.values - observed.values))
            for hidden, observed in zip(records, released, strict=True)
        )
    )


def verify_frozen_program(
    selected: tuple[Program, ...],
    audit_records: list[Trajectory],
    table: dict[str, dict[str, float]],
    manual_values: np.ndarray,
) -> float:
    """Fail closed unless the interpreted scalar is the evaluated winner."""

    selected_names = tuple(program.name for program in selected)
    if selected_names != (FROZEN_KURAMOTO_PROGRAM,):
        raise AssertionError(
            "Kuramoto post-hoc diagnostic requires exactly the frozen program "
            f"{FROZEN_KURAMOTO_PROGRAM!r}; received {selected_names!r}"
        )
    evaluated = np.asarray(
        [selected[0].evaluate(table[item.identifier]) for item in audit_records],
        dtype=np.float64,
    )
    manual = np.asarray(manual_values, dtype=np.float64)
    if evaluated.shape != manual.shape:
        raise AssertionError(
            f"manual/evaluated scalar shape mismatch: {manual.shape} != {evaluated.shape}"
        )
    error = float(np.max(np.abs(evaluated - manual)))
    if error > 1.0e-12:
        raise AssertionError(
            "manual Kuramoto scalar disagrees with frozen Program.evaluate: "
            f"max abs difference {error}"
        )
    return error


def stratified_spearman_interval(
    values: np.ndarray,
    order: np.ndarray,
    labels: np.ndarray,
    *,
    seed: int,
    repetitions: int = 10_000,
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    strata = [np.flatnonzero(labels == label) for label in sorted(set(labels.tolist()))]
    samples = np.empty(repetitions, dtype=np.float64)
    for repetition in range(repetitions):
        indices = np.concatenate(
            [rng.choice(stratum, size=len(stratum), replace=True) for stratum in strata]
        )
        samples[repetition] = float(spearmanr(values[indices], order[indices]).statistic)
    return float(np.quantile(samples, 0.025)), float(np.quantile(samples, 0.975))


def _weak_component_difference_bootstrap(
    selected: np.ndarray,
    increment_cv: np.ndarray,
    radial_pe: np.ndarray,
    second_order: np.ndarray,
    labels: np.ndarray,
    *,
    repetitions: int = 10_000,
) -> dict[str, object]:
    """Compare dependent rank correlations inside the weak-coupling class."""

    weak = np.flatnonzero(labels == 1)
    rng = np.random.default_rng(COMPONENT_DIFFERENCE_BOOTSTRAP_SEED)
    selected_minus_increment = np.empty(repetitions, dtype=np.float64)
    selected_minus_oriented_entropy = np.empty(repetitions, dtype=np.float64)
    for repetition in range(repetitions):
        indices = rng.choice(weak, size=len(weak), replace=True)
        selected_rho = float(spearmanr(selected[indices], second_order[indices]).statistic)
        increment_rho = float(
            spearmanr(increment_cv[indices], second_order[indices]).statistic
        )
        oriented_entropy_rho = float(
            spearmanr(-radial_pe[indices], second_order[indices]).statistic
        )
        selected_minus_increment[repetition] = selected_rho - increment_rho
        selected_minus_oriented_entropy[repetition] = (
            selected_rho - oriented_entropy_rho
        )

    selected_rho = float(spearmanr(selected[weak], second_order[weak]).statistic)
    increment_rho = float(spearmanr(increment_cv[weak], second_order[weak]).statistic)
    oriented_entropy_rho = float(spearmanr(-radial_pe[weak], second_order[weak]).statistic)
    return {
        "seed": COMPONENT_DIFFERENCE_BOOTSTRAP_SEED,
        "repetitions": repetitions,
        "sampling_unit": "weak-coupling audit record",
        "resampling": "sample 120 weak-coupling records with replacement; recompute dependent Spearman correlations on the same resample",
        "selected_minus_increment_cv": {
            "point_difference": selected_rho - increment_rho,
            "95_interval": [
                float(np.quantile(selected_minus_increment, 0.025)),
                float(np.quantile(selected_minus_increment, 0.975)),
            ],
        },
        "selected_minus_oriented_radial_permutation_entropy": {
            "point_difference": selected_rho - oriented_entropy_rho,
            "95_interval": [
                float(np.quantile(selected_minus_oriented_entropy, 0.025)),
                float(np.quantile(selected_minus_oriented_entropy, 0.975)),
            ],
            "orientation": "negative radial permutation entropy, chosen to align its observed association with mean |Z_2|",
        },
    }


def _channel_subset_rows(
    records: list[KuramotoHiddenRecord],
) -> tuple[list[dict[str, object]], dict[tuple[int, ...], np.ndarray]]:
    labels = np.asarray([item.trajectory.label for item in records], dtype=int)
    rows: list[dict[str, object]] = []
    values_by_subset: dict[tuple[int, ...], np.ndarray] = {}
    for size in range(1, 6):
        for subset in combinations(range(5), size):
            values = np.asarray(
                [selected_scalar(item.trajectory.values[:, subset])[0] for item in records]
            )
            values_by_subset[subset] = values
            rows.append(
                {
                    "channel_count": size,
                    "channels": "-".join(str(index) for index in subset),
                    "orientation_free_auc": orientation_free_auc(values, labels),
                }
            )
    return rows, values_by_subset


def _pairing_sensitivity(
    audit_records: list[Trajectory],
    table: dict[str, dict[str, float]],
    base_names: tuple[str, ...],
    selected: tuple[Program, ...],
    threshold: float,
) -> tuple[dict[str, object], np.ndarray]:
    labels = np.asarray([item.label for item in audit_records], dtype=int)
    before = fixed_squash(representation(audit_records, table, base_names, ()))
    after = fixed_squash(representation(audit_records, table, base_names, selected))
    label_values = sorted(set(labels.tolist()))
    left = np.flatnonzero(labels == label_values[0])
    right = np.flatnonzero(labels == label_values[1])
    before_all = np.max(np.abs(before[left, None, :] - before[right][None, :, :]), axis=2) <= threshold
    after_all = np.max(np.abs(after[left, None, :] - after[right][None, :, :]), axis=2) <= threshold

    rng = np.random.default_rng(PAIRING_SEED)
    risks = np.empty(20_000, dtype=np.float64)
    for repetition in range(len(risks)):
        risks[repetition] = float(np.mean(after_all[np.arange(len(left)), rng.permutation(len(right))]))

    fixed_left, fixed_right = independent_cross_label_pairs(labels, stable_seed("kuramoto"))
    fixed_after = (
        np.max(np.abs(after[fixed_left] - after[fixed_right]), axis=1) <= threshold
    )
    residual_right = fixed_right[fixed_after]
    summary = {
        "all_cross_label_pairs": {
            "pairs": int(before_all.size),
            "risk_before": float(np.mean(before_all)),
            "risk_after": float(np.mean(after_all)),
            "absolute_reduction": float(np.mean(before_all) - np.mean(after_all)),
        },
        "random_perfect_matchings": {
            "seed": PAIRING_SEED,
            "repetitions": int(len(risks)),
            "post_risk_95_range": [
                float(np.quantile(risks, 0.025)),
                float(np.quantile(risks, 0.975)),
            ],
            "reduction_95_range": [
                float(np.quantile(1.0 - risks, 0.025)),
                float(np.quantile(1.0 - risks, 0.975)),
            ],
            "interpretation": "pairing-sensitivity range, not a confidence interval",
        },
        "released_matching": {
            "seed": stable_seed("kuramoto"),
            "post_collisions": int(np.sum(fixed_after)),
            "trials": int(len(fixed_after)),
        },
    }
    return summary, residual_right


def _plot(
    records: list[dict[str, object]],
    subset_rows: list[dict[str, object]],
    residual_identifiers: set[str],
    output: Path,
) -> None:
    frame = pd.DataFrame(records)
    subset_frame = pd.DataFrame(subset_rows)
    colors = {0: "#1D5F8A", 1: "#B95F2A"}
    labels = {0: "strong coupling", 1: "weak coupling"}
    fig, axes = plt.subplots(1, 2, figsize=(8.3, 3.55), gridspec_kw={"width_ratios": [1.15, 0.85]})

    for label in (0, 1):
        part = frame[frame["label"] == label]
        axes[0].scatter(
            part["mean_second_order"],
            part["selected_scalar"],
            s=17,
            alpha=0.68,
            color=colors[label],
            edgecolors="none",
            label=labels[label],
        )
    residual = frame[frame["identifier"].isin(residual_identifiers)]
    axes[0].scatter(
        residual["mean_second_order"],
        residual["selected_scalar"],
        s=36,
        facecolors="none",
        edgecolors="#9F3D45",
        linewidths=1.0,
        label="weak records in 10 residual pairs",
    )
    axes[0].set_xlabel(r"latent second-harmonic order $\overline{R}_2$")
    axes[0].set_ylabel(r"frozen scalar $f_\star$")
    axes[0].grid(color="#D5DDE1", linewidth=0.6, alpha=0.8)
    axes[0].spines[["top", "right"]].set_visible(False)
    axes[0].legend(frameon=False, fontsize=7, loc="lower right")
    axes[0].text(
        0.03,
        0.96,
        "pooled " + r"$\rho_S=0.907$" + "\nweak class " + r"$\rho_S=0.944$",
        transform=axes[0].transAxes,
        ha="left",
        va="top",
        fontsize=8,
        color="#17232D",
    )

    for size in range(1, 6):
        values = subset_frame.loc[
            subset_frame["channel_count"] == size, "orientation_free_auc"
        ].to_numpy()
        jitter = np.linspace(-0.09, 0.09, len(values)) if len(values) > 1 else np.zeros(1)
        axes[1].scatter(
            size + jitter,
            values,
            s=24,
            color="#69747D",
            alpha=0.72,
            edgecolors="none",
        )
        axes[1].scatter(size, np.mean(values), s=43, color="#1D5F8A", zorder=3)
    means = subset_frame.groupby("channel_count")["orientation_free_auc"].mean()
    axes[1].plot(means.index, means.values, color="#1D5F8A", linewidth=1.2, zorder=2)
    axes[1].axhline(0.5, color="#B95F2A", linestyle="--", linewidth=0.9)
    axes[1].set_xlim(0.65, 5.35)
    axes[1].set_ylim(0.47, 1.02)
    axes[1].set_xticks(range(1, 6))
    axes[1].set_xlabel("channels retained")
    axes[1].set_ylabel("orientation-free label AUC")
    axes[1].grid(axis="y", color="#D5DDE1", linewidth=0.6, alpha=0.8)
    axes[1].spines[["top", "right"]].set_visible(False)

    fig.tight_layout(w_pad=2.2)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    fig.savefig(output.with_suffix(".png"), dpi=220)
    plt.close(fig)


def run_posthoc_diagnostic(
    audit_records: list[Trajectory],
    table: dict[str, dict[str, float]],
    base_names: tuple[str, ...],
    selected: tuple[Program, ...],
    *,
    threshold: float,
    output_dir: Path,
    figure_dir: Path,
) -> dict[str, object]:
    hidden = kuramoto_with_hidden_order()
    regeneration_error = exact_regeneration_error(hidden)
    if regeneration_error != 0.0:
        raise AssertionError(f"Kuramoto post-hoc regeneration mismatch: {regeneration_error}")
    if [item.identifier for item in audit_records] != [item.trajectory.identifier for item in hidden]:
        raise AssertionError("post-hoc audit order differs from the evaluated split")

    labels = np.asarray([item.trajectory.label for item in hidden], dtype=int)
    scalar_parts = [selected_scalar(item.trajectory.values) for item in hidden]
    selected_values = np.asarray([parts[0] for parts in scalar_parts])
    increment_cv = np.asarray([parts[1] for parts in scalar_parts])
    radial_pe = np.asarray([parts[2] for parts in scalar_parts])
    program_agreement_error = verify_frozen_program(
        selected,
        audit_records,
        table,
        selected_values,
    )
    mean_order = np.asarray([item.mean_order for item in hidden])
    mean_second_order = np.asarray([item.mean_second_order for item in hidden])
    std_real_second_order = np.asarray(
        [item.std_real_second_order for item in hidden]
    )
    coupling = np.asarray([item.coupling for item in hidden])
    first_order_rho = float(spearmanr(selected_values, mean_order).statistic)
    first_order_interval = stratified_spearman_interval(
        selected_values,
        mean_order,
        labels,
        seed=FIRST_ORDER_BOOTSTRAP_SEED,
    )
    second_order_rho = float(
        spearmanr(selected_values, mean_second_order).statistic
    )
    second_order_interval = stratified_spearman_interval(
        selected_values,
        mean_second_order,
        labels,
        seed=SECOND_ORDER_BOOTSTRAP_SEED,
    )

    subset_rows, _ = _channel_subset_rows(hidden)
    subset_frame = pd.DataFrame(subset_rows)
    subset_summary = {
        str(size): {
            "subset_count": int(np.sum(subset_frame["channel_count"] == size)),
            "mean_auc": float(
                subset_frame.loc[
                    subset_frame["channel_count"] == size, "orientation_free_auc"
                ].mean()
            ),
            "minimum_auc": float(
                subset_frame.loc[
                    subset_frame["channel_count"] == size, "orientation_free_auc"
                ].min()
            ),
            "maximum_auc": float(
                subset_frame.loc[
                    subset_frame["channel_count"] == size, "orientation_free_auc"
                ].max()
            ),
        }
        for size in range(1, 6)
    }

    pairing, residual_indices = _pairing_sensitivity(
        audit_records,
        table,
        base_names,
        selected,
        threshold,
    )
    residual_identifiers = {audit_records[index].identifier for index in residual_indices}
    weak_mask = labels == 1
    residual_mask = np.asarray(
        [item.trajectory.identifier in residual_identifiers for item in hidden], dtype=bool
    )
    residual_summary = {
        "weak_records_in_released_residual_pairs": int(np.sum(residual_mask)),
        "mean_first_order_residual": float(np.mean(mean_order[residual_mask])),
        "mean_first_order_other_weak": float(
            np.mean(mean_order[weak_mask & ~residual_mask])
        ),
        "mean_second_order_residual": float(
            np.mean(mean_second_order[residual_mask])
        ),
        "mean_second_order_other_weak": float(
            np.mean(mean_second_order[weak_mask & ~residual_mask])
        ),
        "mean_coupling_residual": float(np.mean(coupling[residual_mask])),
        "mean_coupling_other_weak": float(np.mean(coupling[weak_mask & ~residual_mask])),
    }

    weak_second_order_interval = stratified_spearman_interval(
        selected_values[weak_mask],
        mean_second_order[weak_mask],
        labels[weak_mask],
        seed=SECOND_ORDER_BOOTSTRAP_SEED,
    )
    component_differences = _weak_component_difference_bootstrap(
        selected_values,
        increment_cv,
        radial_pe,
        mean_second_order,
        labels,
    )

    record_rows = [
        {
            "identifier": item.trajectory.identifier,
            "label": item.trajectory.label,
            "coupling": item.coupling,
            "mean_first_order": item.mean_order,
            "mean_second_order": item.mean_second_order,
            "std_real_second_order": item.std_real_second_order,
            "selected_scalar": float(selected_values[index]),
            "increment_cv": float(increment_cv[index]),
            "radial_permutation_entropy": float(radial_pe[index]),
            "weak_record_in_released_residual_pair": bool(residual_mask[index]),
        }
        for index, item in enumerate(hidden)
    ]

    result: dict[str, object] = {
        "status": "COMPUTED_POST_HOC_NO_CLAIM_UPGRADE",
        "interpretation": "post-hoc order-sensitive multichannel surrogate within the named Kuramoto generator",
        "prohibited_interpretations": [
            "a conserved quantity",
            "the Kuramoto order parameter itself",
            "a critical-coupling estimate",
            "a universal order coordinate",
            "evidence beyond the named controlled generator",
        ],
        "audit_seed": AUDIT_SEED,
        "exact_regeneration_max_abs_difference": regeneration_error,
        "selected_program_verification": {
            "expected_program": FROZEN_KURAMOTO_PROGRAM,
            "selected_programs": [program.name for program in selected],
            "manual_vs_program_evaluate_max_abs_difference": program_agreement_error,
            "failure_tolerance": 1.0e-12,
        },
        "records": len(hidden),
        "selected_scalar_orientation_free_auc": orientation_free_auc(selected_values, labels),
        "latent_order_diagnostics": {
            "first_harmonic_mean_abs_z1": {
                "orientation_free_auc": orientation_free_auc(mean_order, labels),
                "selected_scalar_spearman_pooled": first_order_rho,
                "selected_scalar_spearman_strong_class": float(
                    spearmanr(selected_values[~weak_mask], mean_order[~weak_mask]).statistic
                ),
                "selected_scalar_spearman_weak_class": float(
                    spearmanr(selected_values[weak_mask], mean_order[weak_mask]).statistic
                ),
                "stratified_record_bootstrap": {
                    "seed": FIRST_ORDER_BOOTSTRAP_SEED,
                    "repetitions": 10_000,
                    "sampling": "resample 120 records with replacement within each label stratum",
                    "pooled_spearman_95_interval": list(first_order_interval),
                },
            },
            "second_harmonic_mean_abs_z2": {
                "orientation_free_auc": orientation_free_auc(mean_second_order, labels),
                "selected_scalar_spearman_pooled": second_order_rho,
                "selected_scalar_spearman_strong_class": float(
                    spearmanr(
                        selected_values[~weak_mask], mean_second_order[~weak_mask]
                    ).statistic
                ),
                "selected_scalar_spearman_weak_class": float(
                    spearmanr(
                        selected_values[weak_mask], mean_second_order[weak_mask]
                    ).statistic
                ),
                "record_bootstrap": {
                    "seed": SECOND_ORDER_BOOTSTRAP_SEED,
                    "repetitions": 10_000,
                    "pooled_sampling": "resample 120 records with replacement within each label stratum",
                    "pooled_spearman_95_interval": list(second_order_interval),
                    "weak_class_sampling": "resample the 120 weak-coupling records with replacement",
                    "weak_class_spearman_95_interval": list(
                        weak_second_order_interval
                    ),
                },
            },
            "std_real_z2": {
                "selected_scalar_spearman_pooled": float(
                    spearmanr(selected_values, std_real_second_order).statistic
                ),
                "selected_scalar_spearman_weak_class": float(
                    spearmanr(
                        selected_values[weak_mask], std_real_second_order[weak_mask]
                    ).statistic
                ),
            },
        },
        "component_diagnostics": {
            "increment_cv": {
                "orientation_free_auc": orientation_free_auc(increment_cv, labels),
                "spearman_vs_mean_abs_z1_pooled": float(
                    spearmanr(increment_cv, mean_order).statistic
                ),
                "spearman_vs_mean_abs_z2_weak_class": float(
                    spearmanr(
                        increment_cv[weak_mask], mean_second_order[weak_mask]
                    ).statistic
                ),
            },
            "radial_permutation_entropy": {
                "orientation_free_auc": orientation_free_auc(radial_pe, labels),
                "spearman_vs_mean_abs_z1_pooled": float(
                    spearmanr(radial_pe, mean_order).statistic
                ),
                "spearman_vs_mean_abs_z2_weak_class": float(
                    spearmanr(
                        radial_pe[weak_mask], mean_second_order[weak_mask]
                    ).statistic
                ),
            },
            "weak_class_dependent_correlation_bootstrap": component_differences,
        },
        "channel_subset_protocol": {
            "subsets": "all 31 nonempty subsets of the five observed sine channels",
            "subset_counts_by_size": {str(size): math.comb(5, size) for size in range(1, 6)},
            "recomputation": "for every subset, rerun centering and global RMS canonicalization before recomputing both primitives and the frozen ratio",
        },
        "channel_subset_auc": subset_summary,
        "pairing_sensitivity": pairing,
        "released_residual_pairs": residual_summary,
        "sine_channel_identity": {
            "noiseless_raw_sensor": "||sin(theta)||^2 = N/2 * (1 - Re Z_2)",
            "noiseless_centered_sensor": "||sin(theta)-mean_t sin(theta)||^2 = N/2*(1-Re Z_2) - 2*sum_j(mean_t sin(theta_j))*sin(theta_j) + sum_j(mean_t sin(theta_j))^2",
            "boundary": "the identity precedes added Gaussian measurement noise; the released noisy array then undergoes channel centering and global RMS normalization, and the frozen ratio is not Z_1, Z_2, or a conserved quantity",
        },
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "kuramoto_posthoc.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    pd.DataFrame(record_rows).to_csv(output_dir / "kuramoto_posthoc_records.csv", index=False)
    subset_frame.to_csv(output_dir / "kuramoto_channel_subsets.csv", index=False)
    _plot(
        record_rows,
        subset_rows,
        residual_identifiers,
        figure_dir / "kuramoto_posthoc.pdf",
    )
    return result
