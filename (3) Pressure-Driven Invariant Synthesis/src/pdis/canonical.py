"""Canonical unary trajectory descriptors.

The recurrence construction is exactly invariant to translation, orthogonal
channel transformations, and positive global scaling (up to floating-point
tie handling). Every function in this module maps one trajectory to scalars;
pairwise comparison happens only after these coordinates are assembled.
"""

from __future__ import annotations

import math
from typing import Iterable

import numpy as np


EPS = 1.0e-12


def as_2d(values: np.ndarray | Iterable[float]) -> np.ndarray:
    x = np.asarray(values, dtype=np.float64)
    if x.ndim == 1:
        x = x[:, None]
    if x.ndim != 2 or x.shape[0] < 12:
        raise ValueError("a trajectory must have shape (time, channels) with at least 12 samples")
    if not np.all(np.isfinite(x)):
        x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
    return x


def normalize_trajectory(values: np.ndarray | Iterable[float]) -> np.ndarray:
    """Center each channel and divide by one global RMS scale."""

    x = as_2d(values)
    centered = x - np.mean(x, axis=0, keepdims=True)
    scale = float(np.sqrt(np.mean(centered * centered)))
    if scale <= EPS:
        return np.zeros_like(centered)
    return centered / scale


def delay_cloud(values: np.ndarray, dimension: int = 3, lag: int = 2) -> np.ndarray:
    x = normalize_trajectory(values)
    start = (dimension - 1) * lag
    cloud = np.concatenate([x[start - j * lag : len(x) - j * lag] for j in range(dimension)], axis=1)
    return cloud


def deterministic_subsample(cloud: np.ndarray, maximum: int = 72) -> np.ndarray:
    if len(cloud) <= maximum:
        return cloud
    index = np.linspace(0, len(cloud) - 1, maximum).round().astype(int)
    return cloud[index]


def pairwise_distances(cloud: np.ndarray) -> np.ndarray:
    gram = cloud @ cloud.T
    squared = np.maximum(np.diag(gram)[:, None] + np.diag(gram)[None, :] - 2.0 * gram, 0.0)
    return np.sqrt(squared)


def recurrence_adjacency(
    values: np.ndarray,
    *,
    quantile: float = 0.12,
    maximum_points: int = 72,
) -> tuple[np.ndarray, float, np.ndarray]:
    """Return quantile-recurrence adjacency, threshold, and distance matrix."""

    if not 0.0 < quantile < 1.0:
        raise ValueError("quantile must lie strictly between zero and one")
    cloud = deterministic_subsample(delay_cloud(values), maximum_points)
    distances = pairwise_distances(cloud)
    upper = distances[np.triu_indices(len(distances), 1)]
    threshold = float(np.quantile(upper, quantile))
    adjacency = (distances <= threshold).astype(np.float64)
    np.fill_diagonal(adjacency, 0.0)
    return adjacency, threshold, distances


def normalized_laplacian_eigenvalues(adjacency: np.ndarray) -> np.ndarray:
    degree = adjacency.sum(axis=1)
    inv_sqrt = np.zeros_like(degree)
    positive = degree > 0
    inv_sqrt[positive] = 1.0 / np.sqrt(degree[positive])
    laplacian = np.eye(len(adjacency)) - inv_sqrt[:, None] * adjacency * inv_sqrt[None, :]
    if np.any(~positive):
        laplacian[~positive, ~positive] = 0.0
    return np.sort(np.linalg.eigvalsh(laplacian))


def graph_transitivity(adjacency: np.ndarray) -> float:
    degrees = adjacency.sum(axis=1)
    triples = float(np.sum(degrees * np.maximum(degrees - 1.0, 0.0)))
    if triples <= EPS:
        return 0.0
    closed_oriented = float(np.trace(adjacency @ adjacency @ adjacency))
    return closed_oriented / triples


def prim_mst_edges(distances: np.ndarray) -> np.ndarray:
    """MST edge lengths, which are the finite H0 persistence lifetimes."""

    n = len(distances)
    selected = np.zeros(n, dtype=bool)
    selected[0] = True
    best = distances[0].copy()
    best[0] = np.inf
    edges: list[float] = []
    for _ in range(n - 1):
        masked = np.where(selected, np.inf, best)
        index = int(np.argmin(masked))
        length = float(masked[index])
        if not np.isfinite(length):
            break
        edges.append(length)
        selected[index] = True
        best = np.minimum(best, distances[index])
    return np.asarray(edges, dtype=np.float64)


def entropy_from_positive(values: np.ndarray) -> float:
    v = np.asarray(values, dtype=np.float64)
    v = v[v > EPS]
    if len(v) <= 1:
        return 0.0
    p = v / np.sum(v)
    return float(-np.sum(p * np.log(p + EPS)) / np.log(len(p)))


def h1_persistence(cloud: np.ndarray) -> tuple[float, float, float]:
    """Return H1 count, maximum lifetime, and total lifetime."""

    try:
        from ripser import ripser
    except ImportError as exc:  # pragma: no cover - dependency is pinned by the project
        raise RuntimeError("ripser is required for H1 persistence features") from exc
    diagram = ripser(cloud, maxdim=1, thresh=np.inf)["dgms"][1]
    if len(diagram) == 0:
        return 0.0, 0.0, 0.0
    finite = diagram[np.isfinite(diagram[:, 1])]
    if len(finite) == 0:
        return 0.0, 0.0, 0.0
    lifetimes = np.maximum(finite[:, 1] - finite[:, 0], 0.0)
    scale = float(np.median(pairwise_distances(cloud)[np.triu_indices(len(cloud), 1)]))
    scale = max(scale, EPS)
    normalized = lifetimes / scale
    strong = normalized[normalized > 0.03]
    if len(strong) == 0:
        return 0.0, 0.0, 0.0
    return float(len(strong)), float(np.max(strong)), float(np.sum(strong))


def spectral_entropy(signal: np.ndarray) -> float:
    centered = signal - np.mean(signal)
    power = np.abs(np.fft.rfft(centered)) ** 2
    if len(power) > 0:
        power[0] = 0.0
    return entropy_from_positive(power)


def permutation_entropy(signal: np.ndarray) -> float:
    if len(signal) < 5:
        return 0.0
    patterns = np.argsort(np.lib.stride_tricks.sliding_window_view(signal, 3), axis=1)
    codes = patterns[:, 0] * 9 + patterns[:, 1] * 3 + patterns[:, 2]
    _, counts = np.unique(codes, return_counts=True)
    p = counts / counts.sum()
    return float(-np.sum(p * np.log(p + EPS)) / math.log(6.0))


def canonical_features(values: np.ndarray) -> dict[str, float]:
    """Compute the declared unary feature library for one trajectory."""

    x = normalize_trajectory(values)
    cloud = deterministic_subsample(delay_cloud(x), 72)
    adjacency, _, distances = recurrence_adjacency(x, maximum_points=72)
    eig = normalized_laplacian_eigenvalues(adjacency)
    degrees = adjacency.sum(axis=1)
    positive_eig = eig[eig > 1.0e-10]

    mst = prim_mst_edges(distances)
    distance_scale = float(np.median(distances[np.triu_indices(len(distances), 1)]))
    distance_scale = max(distance_scale, EPS)
    mst_normalized = mst / distance_scale
    h1_count, h1_max, h1_total = h1_persistence(cloud)

    radial = np.linalg.norm(x, axis=1)
    increments = np.linalg.norm(np.diff(x, axis=0), axis=1)
    path_length = float(np.sum(increments))
    endpoint = float(np.linalg.norm(x[-1] - x[0]))
    covariance = np.cov(x, rowvar=False)
    covariance = np.atleast_2d(covariance)
    covariance_eig = np.maximum(np.linalg.eigvalsh(covariance), 0.0)

    radial_std = float(np.std(radial))
    lag1 = 0.0
    if radial_std > EPS:
        lag1 = float(np.corrcoef(radial[:-1], radial[1:])[0, 1])
    turning = np.diff(radial)
    turning_rate = float(np.mean(turning[:-1] * turning[1:] < 0.0)) if len(turning) > 1 else 0.0

    features = {
        # Baseline coordinates.
        "radial_lag1": lag1,
        "path_tortuosity": path_length / (endpoint + 0.05 * path_length + EPS),
        # Recurrence geometry. lambda_2 is literally eig[1], including zero multiplicity.
        "recurrence_lambda2": float(eig[1]) if len(eig) > 1 else 0.0,
        "recurrence_lambda3": float(eig[2]) if len(eig) > 2 else 0.0,
        "recurrence_spectral_entropy": entropy_from_positive(positive_eig),
        "recurrence_degree_cv": float(np.std(degrees) / (np.mean(degrees) + EPS)),
        "recurrence_transitivity": graph_transitivity(adjacency),
        # Persistent topology.
        "h0_total_persistence": float(np.sum(mst_normalized)),
        "h0_max_persistence": float(np.max(mst_normalized)) if len(mst_normalized) else 0.0,
        "h0_persistence_entropy": entropy_from_positive(mst_normalized),
        "h1_feature_count": h1_count,
        "h1_max_persistence": h1_max,
        "h1_total_persistence": h1_total,
        # Similarity-invariant temporal and covariance coordinates.
        "radial_spectral_entropy": spectral_entropy(radial),
        "radial_permutation_entropy": permutation_entropy(radial),
        "radial_turning_rate": turning_rate,
        "increment_cv": float(np.std(increments) / (np.mean(increments) + EPS)),
        "covariance_entropy": entropy_from_positive(covariance_eig),
    }
    return {name: float(np.nan_to_num(value, nan=0.0, posinf=1.0e6, neginf=-1.0e6)) for name, value in features.items()}
