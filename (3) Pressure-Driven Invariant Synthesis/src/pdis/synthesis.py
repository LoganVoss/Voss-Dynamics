"""Pressure-driven observable selection and frozen-domain evaluation."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .datasets import Trajectory
from .programs import Program
from .statistics import (
    all_cross_label_collision_indicators,
    collision_indicators,
    collision_summary,
    orientation_free_auc,
    paired_bootstrap_reduction,
    stable_seed,
)


@dataclass(frozen=True)
class SynthesisResult:
    base_names: tuple[str, ...]
    selected: tuple[Program, ...]
    source_before: dict[str, dict[str, float | int]]
    source_after: dict[str, dict[str, float | int]]
    rounds: tuple[dict[str, object], ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "base_names": list(self.base_names),
            "selected": [program.to_dict() for program in self.selected],
            "source_before": self.source_before,
            "source_after": self.source_after,
            "rounds": list(self.rounds),
        }


def _program_values(program: Program, records: list[Trajectory], table: dict[str, dict[str, float]]) -> np.ndarray:
    return np.asarray([program.evaluate(table[item.identifier]) for item in records], dtype=np.float64)


def representation(
    records: list[Trajectory],
    table: dict[str, dict[str, float]],
    base_names: tuple[str, ...],
    selected: tuple[Program, ...] | list[Program],
) -> np.ndarray:
    columns = [np.asarray([table[item.identifier][name] for item in records]) for name in base_names]
    columns.extend(_program_values(program, records, table) for program in selected)
    return np.column_stack(columns)


def _similarity_error(
    program: Program,
    records: list[Trajectory],
    table: dict[str, dict[str, float]],
    similarity_table: dict[str, dict[str, float]],
) -> float:
    original = _program_values(program, records, table)
    transformed = _program_values(program, records, similarity_table)
    error = np.abs(original - transformed) / (1.0 + np.abs(original))
    return float(np.max(error))


def _program_error(
    program: Program,
    records: list[Trajectory],
    table: dict[str, dict[str, float]],
    perturbed_table: dict[str, dict[str, float]],
    *,
    quantile: float,
) -> float:
    original = _program_values(program, records, table)
    perturbed = _program_values(program, records, perturbed_table)
    error = np.abs(original - perturbed) / (1.0 + np.abs(original))
    return float(np.quantile(error, quantile))


def _stratified_split(records: list[Trajectory], seed: int) -> tuple[list[Trajectory], list[Trajectory]]:
    labels = np.asarray([item.label for item in records], dtype=int)
    rng = np.random.default_rng(seed)
    selection: list[Trajectory] = []
    audit: list[Trajectory] = []
    for label in sorted(set(labels.tolist())):
        indices = np.flatnonzero(labels == label)
        rng.shuffle(indices)
        cut = max(2, int(round(0.6 * len(indices))))
        selection.extend(records[index] for index in indices[:cut])
        audit.extend(records[index] for index in indices[cut:])
    return selection, audit


def evaluate_domain(
    records: list[Trajectory],
    table: dict[str, dict[str, float]],
    base_names: tuple[str, ...],
    selected: tuple[Program, ...] | list[Program],
    *,
    threshold: float,
) -> dict[str, object]:
    domain = records[0].domain
    labels = np.asarray([item.label for item in records], dtype=int)
    values = representation(records, table, base_names, selected)
    summary = collision_summary(values, labels, threshold=threshold, seed=stable_seed(domain))
    program_auc = {
        program.name: orientation_free_auc(_program_values(program, records, table), labels)
        for program in selected
    }
    return {**summary, "program_auc": program_auc}


def synthesize(
    domains: dict[str, list[Trajectory]],
    table: dict[str, dict[str, float]],
    similarity_table: dict[str, dict[str, float]],
    noise_table: dict[str, dict[str, float]],
    grammar: list[Program],
    *,
    base_names: tuple[str, ...] = ("radial_lag1", "path_tortuosity"),
    threshold: float = 0.10,
    maximum_rounds: int = 3,
    invariance_tolerance: float = 5.0e-6,
    noise_tolerance: float = 0.20,
    minimum_audit_reduction: float = 0.15,
    minimum_audit_auc: float = 0.75,
    minimum_selection_auc: float = 0.70,
) -> SynthesisResult:
    selected: list[Program] = []
    source_before = {
        name: evaluate_domain(records, table, base_names, selected, threshold=threshold)
        for name, records in domains.items()
    }
    rounds: list[dict[str, object]] = []
    inner_splits = {
        name: _stratified_split(records, stable_seed(name + "-inner"))
        for name, records in domains.items()
    }

    all_records = [item for records in domains.values() for item in records]
    remaining = list(grammar)
    for round_index in range(maximum_rounds):
        candidate_rows: list[tuple[tuple[float, float, float, int, str], Program, dict[str, object]]] = []
        for program in remaining:
            invariance_error = _similarity_error(program, all_records, table, similarity_table)
            if invariance_error > invariance_tolerance:
                continue
            noise_error = _program_error(program, all_records, table, noise_table, quantile=0.95)
            if noise_error > noise_tolerance:
                continue
            selection_reductions: dict[str, float] = {}
            audit_reductions: dict[str, float] = {}
            selection_aucs: dict[str, float] = {}
            audit_aucs: dict[str, float] = {}
            for name, (selection_records, audit_records) in inner_splits.items():
                for split_name, split_records, destination in (
                    ("selection", selection_records, selection_reductions),
                    ("audit", audit_records, audit_reductions),
                ):
                    labels = np.asarray([item.label for item in split_records], dtype=int)
                    before_rep = representation(split_records, table, base_names, selected)
                    after_rep = representation(split_records, table, base_names, [*selected, program])
                    before = all_cross_label_collision_indicators(before_rep, labels, threshold=threshold)
                    after = all_cross_label_collision_indicators(after_rep, labels, threshold=threshold)
                    destination[name] = float(np.mean(before) - np.mean(after))
                selection_labels = np.asarray([item.label for item in selection_records], dtype=int)
                audit_labels = np.asarray([item.label for item in audit_records], dtype=int)
                selection_aucs[name] = orientation_free_auc(
                    _program_values(program, selection_records, table), selection_labels
                )
                audit_aucs[name] = orientation_free_auc(
                    _program_values(program, audit_records, table), audit_labels
                )
            minimum = min(audit_reductions.values())
            mean = float(np.mean(list(audit_reductions.values())))
            selection_minimum = min(selection_reductions.values())
            selection_auc_minimum = min(selection_aucs.values())
            audit_auc_minimum = min(audit_aucs.values())
            mean_auc = float(np.mean(list(audit_aucs.values())))
            if (
                minimum < minimum_audit_reduction
                or audit_auc_minimum < minimum_audit_auc
                or selection_auc_minimum < minimum_selection_auc
            ):
                continue
            utility = (
                minimum
                + 0.50 * mean
                + 0.15 * selection_minimum
                + 0.05 * (mean_auc - 0.5)
                - 0.003 * program.complexity
                - 0.10 * noise_error
            )
            key = (utility, minimum, mean, mean_auc, -program.complexity, program.name)
            candidate_rows.append(
                (
                    key,
                    program,
                    {
                        "program": program.to_dict(),
                        "invariance_error": invariance_error,
                        "noise_error_p95": noise_error,
                        "selection_risk_reduction": selection_reductions,
                        "audit_risk_reduction": audit_reductions,
                        "selection_orientation_free_auc": selection_aucs,
                        "audit_orientation_free_auc": audit_aucs,
                        "minimum_audit_auc": audit_auc_minimum,
                        "minimum_selection_auc": selection_auc_minimum,
                        "worst_domain_reduction": minimum,
                        "mean_reduction": mean,
                        "utility": utility,
                    },
                )
            )
        if not candidate_rows:
            break
        candidate_rows.sort(key=lambda row: row[0], reverse=True)
        _, winner, audit = candidate_rows[0]
        if float(audit["mean_reduction"]) <= 0.0 or float(audit["worst_domain_reduction"]) < 0.0:
            break
        selected.append(winner)
        remaining = [program for program in remaining if program != winner]
        rounds.append(
            {
                "round": round_index + 1,
                "winner": audit,
                "top_candidates": [row[2] for row in candidate_rows[:8]],
            }
        )

    source_after = {
        name: evaluate_domain(records, table, base_names, selected, threshold=threshold)
        for name, records in domains.items()
    }
    return SynthesisResult(base_names, tuple(selected), source_before, source_after, tuple(rounds))


def paired_domain_reduction(
    records: list[Trajectory],
    table: dict[str, dict[str, float]],
    base_names: tuple[str, ...],
    selected: tuple[Program, ...],
    *,
    threshold: float,
) -> dict[str, object]:
    labels = np.asarray([item.label for item in records], dtype=int)
    seed = stable_seed(records[0].domain)
    before = collision_indicators(
        representation(records, table, base_names, ()), labels, threshold=threshold, seed=seed
    )
    after = collision_indicators(
        representation(records, table, base_names, selected), labels, threshold=threshold, seed=seed
    )
    mean, low, high = paired_bootstrap_reduction(before, after, seed=seed + 1)
    return {
        "before_collisions": int(np.sum(before)),
        "after_collisions": int(np.sum(after)),
        "trials": int(len(before)),
        "risk_before": float(np.mean(before)),
        "risk_after": float(np.mean(after)),
        "absolute_reduction": mean,
        "bootstrap_95": [low, high],
    }
