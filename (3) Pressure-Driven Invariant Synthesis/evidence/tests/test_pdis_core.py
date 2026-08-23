from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np

from pdis.canonical import (
    canonical_features,
    normalized_laplacian_eigenvalues,
    recurrence_adjacency,
)
from pdis.datasets import Trajectory, sealed_ucr_domain, similarity_transform
from pdis.programs import Program
from pdis.statistics import collision_indicators, one_sided_binomial_upper
from pdis.synthesis import representation, synthesize


def toy_records() -> list[Trajectory]:
    values = np.sin(np.linspace(0, 8 * np.pi, 96))[:, None]
    return [
        Trajectory(f"toy-{index}", "toy", index % 2, values + 0.01 * index)
        for index in range(8)
    ]


def test_quantile_recurrence_similarity_invariance() -> None:
    t = np.linspace(0, 12, 180)
    values = np.column_stack([np.sin(t), np.cos(1.7 * t), np.sin(0.4 * t + 0.3)])
    transformed = similarity_transform(values, 123)
    a, _, _ = recurrence_adjacency(values)
    b, _, _ = recurrence_adjacency(transformed)
    assert np.array_equal(a, b)


def test_lambda2_uses_second_smallest_including_zero_multiplicity() -> None:
    path4 = np.zeros((4, 4))
    for left, right in ((0, 1), (1, 2), (2, 3)):
        path4[left, right] = path4[right, left] = 1
    disconnected = np.zeros((4, 4))
    disconnected[0, 1] = disconnected[1, 0] = 1
    disconnected[2, 3] = disconnected[3, 2] = 1
    assert np.isclose(normalized_laplacian_eigenvalues(path4)[1], 0.5)
    assert np.isclose(normalized_laplacian_eigenvalues(disconnected)[1], 0.0)


def test_synthesized_program_reenters_unary_representation() -> None:
    records = toy_records()
    table = {
        item.identifier: {"base_a": 0.0, "base_b": 0.0, "separator": float(item.label)}
        for item in records
    }
    before = representation(records, table, ("base_a", "base_b"), ())
    program = Program("primitive", "separator")
    after = representation(records, table, ("base_a", "base_b"), (program,))
    assert before.shape == (8, 2)
    assert after.shape == (8, 3)
    assert np.array_equal(after[:, -1], np.asarray([item.label for item in records]))


def test_threshold_collisions_are_nested_after_append() -> None:
    records = toy_records()
    labels = np.asarray([item.label for item in records])
    table = {
        item.identifier: {"base_a": 0.0, "base_b": 0.0, "separator": float(item.label)}
        for item in records
    }
    before = representation(records, table, ("base_a", "base_b"), ())
    after = representation(records, table, ("base_a", "base_b"), (Program("primitive", "separator"),))
    c0 = collision_indicators(before, labels, threshold=0.2, seed=7)
    c1 = collision_indicators(after, labels, threshold=0.2, seed=7)
    assert np.all(c1 <= c0)
    assert np.sum(c1) == 0


def test_finite_sample_synthesis_splits_a_collision_class() -> None:
    records = toy_records()
    table = {
        item.identifier: {"base_a": 0.0, "base_b": 0.0, "separator": float(item.label)}
        for item in records
    }
    result = synthesize(
        {"toy": records},
        table,
        table,
        table,
        [Program("primitive", "separator")],
        base_names=("base_a", "base_b"),
        threshold=0.2,
        maximum_rounds=1,
    )
    assert [program.name for program in result.selected] == ["separator"]
    assert result.source_after["toy"]["collisions"] == 0


def test_zero_collision_bound_is_exact_rule_of_three_generalization() -> None:
    assert np.isclose(one_sided_binomial_upper(0, 299), 1 - 0.05 ** (1 / 299))
    assert one_sided_binomial_upper(0, 299) < 0.0101


def test_feature_extraction_is_deterministic() -> None:
    values = np.sin(np.linspace(0, 20, 180))[:, None]
    assert canonical_features(values) == canonical_features(values)


def test_real_loader_reads_only_declared_test_member() -> None:
    root = Path(__file__).resolve().parents[2]
    archive = root / "evidence" / "data" / "raw" / "ECGFiveDays.zip"
    records = sealed_ucr_domain("ECGFiveDays", archive, per_class=3, seed=9)
    assert len(records) == 6
    assert all(item.domain == "ecgfivedays" for item in records)
    assert {item.label for item in records} == {0, 1}
