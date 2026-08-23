"""
Geometric phrase readout and recursive promotion.

The historical implementation called its auxiliary quantity an MDL score:

    L_total(P) = L(P) + L(field | P, clusters)

However, the code deterministically constructs one nearest-anchor clustering
and scores that one phrase. It does not search a candidate model class or
implement a decodable MDL code. The legacy function and result-field names are
retained for compatibility. The demonstrated result is that the phrase readout
is many-to-one on microstates.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .psi import Field, default_anchors


ANCHOR_NAMES = ("NORTH", "EAST", "SOUTH", "WEST", "NORTHEAST", "ZENITH")


@dataclass
class Cluster:
    indices: np.ndarray
    centroid: np.ndarray
    nearest_anchor: int
    lock_radius: float


def assign_clusters(field: Field, lock: float = 1.15) -> list[Cluster]:
    if field.n == 0:
        return []
    # nearest-anchor basins, then split by gap if a basin is empty of lock
    d = np.linalg.norm(field.r[:, None, :] - field.p[None, :, :], axis=-1)
    nearest = np.argmin(d, axis=1)
    mind = d[np.arange(field.n), nearest]
    clusters: list[Cluster] = []
    for a in range(len(field.Q)):
        idx = np.where((nearest == a) & (mind < lock))[0]
        if len(idx) == 0:
            continue
        centroid = field.r[idx].mean(axis=0)
        clusters.append(
            Cluster(
                indices=idx,
                centroid=centroid,
                nearest_anchor=a,
                lock_radius=lock,
            )
        )
    free = np.where(mind >= lock)[0]
    if len(free) >= 2:
        clusters.append(
            Cluster(
                indices=free,
                centroid=field.r[free].mean(axis=0),
                nearest_anchor=-1,
                lock_radius=lock,
            )
        )
    elif len(free) == 1:
        clusters.append(
            Cluster(
                indices=free,
                centroid=field.r[free][0],
                nearest_anchor=-1,
                lock_radius=lock,
            )
        )
    return clusters


def _anchor_name(idx: int) -> str:
    if 0 <= idx < len(ANCHOR_NAMES):
        return ANCHOR_NAMES[idx]
    return f"A{idx}"


def phrase_from_clusters(clusters: list[Cluster]) -> str:
    parts = []
    for c in sorted(clusters, key=lambda x: (-len(x.indices), x.nearest_anchor)):
        if c.nearest_anchor < 0:
            parts.append(f"FREE:{len(c.indices)}")
        else:
            parts.append(f"{_anchor_name(c.nearest_anchor)}:{len(c.indices)}")
    return " ".join(parts) if parts else "VOID"


def mdl_score(field: Field, clusters: list[Cluster]) -> dict:
    """
    Legacy geometric score (historically named ``mdl_score``).

    ``L_phrase``: heuristic name cost for used anchors / free basins.
    ``L_data``: Gaussian-style position residual about phrase centers.

    This is not a candidate search or a complete, decodable prefix code.
    """
    if not clusters:
        return {
            "L_phrase": 0.0,
            "L_data": 0.0,
            "L_total": 0.0,
            "phrase": "VOID",
            "residual_var": 0.0,
        }
    alphabet = len(field.Q) + 1
    L_phrase = float(len(clusters) * np.log2(alphabet))
    resid = []
    for c in clusters:
        center = field.p[c.nearest_anchor] if c.nearest_anchor >= 0 else c.centroid
        delta = field.r[c.indices] - center
        resid.append(np.sum(delta**2))
    sse = float(np.sum(resid))
    var = sse / max(field.n, 1)
    # bits for N 2D residuals ~ (N) * 0.5 log2(2π e σ²) * 2, floor σ²
    sigma2 = max(var, 1e-8)
    L_data = float(field.n * np.log2(2.0 * np.pi * np.e * sigma2))
    phrase = phrase_from_clusters(clusters)
    return {
        "L_phrase": L_phrase,
        "L_data": L_data,
        "L_total": L_phrase + L_data,
        "phrase": phrase,
        "residual_var": var,
        "n_clusters": len(clusters),
    }


def observe(field: Field, lock: float = 1.15) -> dict:
    clusters = assign_clusters(field, lock=lock)
    score = mdl_score(field, clusters)
    score["occupancy"] = occupancy(field)
    score["clusters"] = clusters
    return score


def occupancy(field: Field) -> list[int]:
    d = np.linalg.norm(field.r[:, None, :] - field.p[None, :, :], axis=-1)
    nearest = np.argmin(d, axis=1)
    return [int(np.sum(nearest == a)) for a in range(len(field.Q))]


def promote(field: Field, clusters: list[Cluster], charge_scale: float = 1.3) -> None:
    """Original recursion: selected phrases become new high-Q anchors."""
    if not clusters:
        return
    new_p = []
    new_Q = []
    new_th = []
    new_ph = []
    for c in clusters:
        if c.nearest_anchor >= 0:
            new_p.append(c.centroid)
            new_Q.append(field.Q[c.nearest_anchor] * charge_scale)
            new_th.append(field.theta_a[c.nearest_anchor])
            new_ph.append(field.phi_a[c.nearest_anchor])
        else:
            new_p.append(c.centroid)
            new_Q.append(1.2 * charge_scale)
            new_th.append(float(np.mean(field.theta[c.indices])))
            new_ph.append(float(np.mean(field.phi[c.indices])))
    field.p = np.vstack([field.p, np.asarray(new_p)])
    field.Q = np.concatenate([field.Q, np.asarray(new_Q)])
    field.theta_a = np.concatenate([field.theta_a, np.asarray(new_th)])
    field.phi_a = np.concatenate([field.phi_a, np.asarray(new_ph)])


def phrase_centers(field: Field, clusters: list[Cluster]) -> np.ndarray:
    if not clusters:
        return field.p.copy()
    return np.stack(
        [
            field.p[c.nearest_anchor] if c.nearest_anchor >= 0 else c.centroid
            for c in clusters
        ]
    )


def reset_tokens(field: Field, seed: int, spread: float = 3.0) -> None:
    rng = np.random.default_rng(seed)
    n = field.n
    field.r = rng.uniform(-spread, spread, size=(n, 2))
    field.v = rng.normal(0.0, 0.05, size=(n, 2))
    field.theta = rng.uniform(0.0, 2 * np.pi, size=n)
    field.phi = rng.uniform(0.0, 2 * np.pi, size=n)


# silence unused import if default_anchors is unused in some flows
_ = default_anchors
