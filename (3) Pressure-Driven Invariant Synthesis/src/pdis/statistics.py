"""Leakage-resistant collision and separation statistics."""

from __future__ import annotations

import hashlib

import numpy as np
from scipy.stats import beta


def stable_seed(text: str, base: int = 20260822) -> int:
    digest = hashlib.sha256(text.encode("utf-8")).digest()
    return (int.from_bytes(digest[:8], "big") + base) % (2**32 - 1)


def fixed_squash(values: np.ndarray) -> np.ndarray:
    """Apply a fixed, source-independent monotone map to every coordinate."""

    x = np.asarray(values, dtype=np.float64)
    if x.ndim == 1:
        x = x[:, None]
    return (2.0 / np.pi) * np.arctan(x)


def independent_cross_label_pairs(labels: np.ndarray, seed: int) -> tuple[np.ndarray, np.ndarray]:
    unique = sorted(set(np.asarray(labels, dtype=int).tolist()))
    if len(unique) != 2:
        raise ValueError("the benchmark collision witness is binary")
    left = np.flatnonzero(labels == unique[0]).copy()
    right = np.flatnonzero(labels == unique[1]).copy()
    rng = np.random.default_rng(seed)
    rng.shuffle(left)
    rng.shuffle(right)
    count = min(len(left), len(right))
    return left[:count], right[:count]


def collision_indicators(
    representation: np.ndarray,
    labels: np.ndarray,
    *,
    threshold: float,
    seed: int,
) -> np.ndarray:
    ranked = fixed_squash(representation)
    left, right = independent_cross_label_pairs(labels, seed)
    distances = np.max(np.abs(ranked[left] - ranked[right]), axis=1)
    return distances <= threshold


def all_cross_label_collision_indicators(
    representation: np.ndarray,
    labels: np.ndarray,
    *,
    threshold: float,
) -> np.ndarray:
    """Dependent all-pairs U-statistic indicators for source selection only."""

    transformed = fixed_squash(representation)
    unique = sorted(set(np.asarray(labels, dtype=int).tolist()))
    if len(unique) != 2:
        raise ValueError("the benchmark collision witness is binary")
    left = transformed[labels == unique[0]]
    right = transformed[labels == unique[1]]
    distances = np.max(np.abs(left[:, None, :] - right[None, :, :]), axis=2)
    return (distances <= threshold).reshape(-1)


def collision_summary(
    representation: np.ndarray,
    labels: np.ndarray,
    *,
    threshold: float,
    seed: int,
    alpha: float = 0.05,
) -> dict[str, float | int]:
    indicators = collision_indicators(representation, labels, threshold=threshold, seed=seed)
    collisions = int(np.sum(indicators))
    trials = int(len(indicators))
    return {
        "collisions": collisions,
        "trials": trials,
        "risk": float(collisions / trials) if trials else 0.0,
        "upper_95": one_sided_binomial_upper(collisions, trials, alpha=alpha),
    }


def one_sided_binomial_upper(successes: int, trials: int, alpha: float = 0.05) -> float:
    if trials <= 0:
        return 1.0
    if successes == trials:
        return 1.0
    if successes == 0:
        return float(1.0 - alpha ** (1.0 / trials))
    return float(beta.ppf(1.0 - alpha, successes + 1, trials - successes))


def orientation_free_auc(values: np.ndarray, labels: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float64)
    labels = np.asarray(labels, dtype=int)
    unique = sorted(set(labels.tolist()))
    if len(unique) != 2:
        return 0.5
    a = values[labels == unique[0]]
    b = values[labels == unique[1]]
    wins = 0.0
    for left in a:
        wins += float(np.sum(left < b)) + 0.5 * float(np.sum(left == b))
    auc = wins / max(len(a) * len(b), 1)
    return float(max(auc, 1.0 - auc))


def paired_bootstrap_reduction(
    before: np.ndarray,
    after: np.ndarray,
    *,
    seed: int,
    repetitions: int = 4000,
) -> tuple[float, float, float]:
    before = np.asarray(before, dtype=float)
    after = np.asarray(after, dtype=float)
    if len(before) != len(after):
        raise ValueError("paired arrays must have the same length")
    difference = before - after
    if len(difference) == 0:
        return 0.0, 0.0, 0.0
    rng = np.random.default_rng(seed)
    samples = rng.choice(difference, size=(repetitions, len(difference)), replace=True).mean(axis=1)
    return float(np.mean(difference)), float(np.quantile(samples, 0.025)), float(np.quantile(samples, 0.975))
