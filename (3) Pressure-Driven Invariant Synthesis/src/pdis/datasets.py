"""Controlled and sealed-domain trajectory datasets."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import zipfile

import numpy as np


@dataclass(frozen=True)
class Trajectory:
    identifier: str
    domain: str
    label: int
    values: np.ndarray


def _rk4_step(state: np.ndarray, dt: float, derivative) -> np.ndarray:
    k1 = derivative(state)
    k2 = derivative(state + 0.5 * dt * k1)
    k3 = derivative(state + 0.5 * dt * k2)
    k4 = derivative(state + dt * k3)
    return state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def logistic_domain(n_per_class: int = 28, seed: int = 2026082201) -> list[Trajectory]:
    rng = np.random.default_rng(seed)
    result: list[Trajectory] = []
    for label, interval in ((0, (3.20, 3.42)), (1, (3.90, 3.99))):
        for index in range(n_per_class):
            r = rng.uniform(*interval)
            value = rng.uniform(0.1, 0.9)
            series = []
            for step in range(640):
                value = r * value * (1.0 - value)
                if step >= 384:
                    series.append(value)
            values = np.asarray(series)[:, None]
            values += rng.normal(0.0, 0.003, size=values.shape)
            result.append(Trajectory(f"logistic-{label}-{index:03d}", "logistic", label, values))
    return result


def lorenz_domain(n_per_class: int = 28, seed: int = 2026082202) -> list[Trajectory]:
    rng = np.random.default_rng(seed)
    result: list[Trajectory] = []
    for label, rho_interval in ((0, (17.0, 21.0)), (1, (27.0, 31.0))):
        for index in range(n_per_class):
            rho = rng.uniform(*rho_interval)
            state = rng.normal(size=3)

            def derivative(z: np.ndarray) -> np.ndarray:
                x, y, zz = z
                return np.asarray([10.0 * (y - x), x * (rho - zz) - y, x * y - (8.0 / 3.0) * zz])

            kept = []
            for step in range(4400):
                state = _rk4_step(state, 0.008, derivative)
                if step >= 2864 and step % 6 == 0:
                    kept.append(state.copy())
            values = np.asarray(kept[:256])
            values += rng.normal(0.0, 0.01, size=values.shape)
            result.append(Trajectory(f"lorenz-{label}-{index:03d}", "lorenz", label, values))
    return result


def kuramoto_domain(n_per_class: int = 28, seed: int = 2026082203) -> list[Trajectory]:
    rng = np.random.default_rng(seed)
    result: list[Trajectory] = []
    oscillators = 5
    for label, coupling_interval in ((0, (2.8, 4.0)), (1, (0.05, 0.45))):
        for index in range(n_per_class):
            coupling = rng.uniform(*coupling_interval)
            omega = rng.normal(1.0, 0.25, size=oscillators)
            theta = rng.uniform(-np.pi, np.pi, size=oscillators)
            kept = []
            for step in range(2300):
                differences = theta[None, :] - theta[:, None]
                theta += 0.025 * (omega + coupling * np.mean(np.sin(differences), axis=1))
                if step >= 1020 and step % 5 == 0:
                    kept.append(np.sin(theta.copy()))
            values = np.asarray(kept[:256])
            values += rng.normal(0.0, 0.015, size=values.shape)
            result.append(Trajectory(f"kuramoto-{label}-{index:03d}", "kuramoto", label, values))
    return result


def synthetic_domains(n_per_class: int = 28) -> dict[str, list[Trajectory]]:
    return {
        "logistic": logistic_domain(n_per_class=n_per_class),
        "lorenz": lorenz_domain(n_per_class=n_per_class),
        "kuramoto": kuramoto_domain(n_per_class=n_per_class),
    }


def _read_ts_from_zip(path: Path, member: str) -> tuple[np.ndarray, np.ndarray]:
    with zipfile.ZipFile(path) as archive:
        text = archive.read(member).decode("utf-8")
    data_started = False
    rows: list[np.ndarray] = []
    labels: list[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if not data_started:
            if line.lower() == "@data":
                data_started = True
            continue
        signal, label = line.rsplit(":", 1)
        rows.append(np.asarray([float(value) for value in signal.split(",")], dtype=np.float64))
        labels.append(label.strip())
    unique = {label: index for index, label in enumerate(sorted(set(labels)))}
    return np.asarray(rows), np.asarray([unique[label] for label in labels], dtype=int)


def sealed_ucr_domain(
    name: str,
    archive_path: Path,
    *,
    per_class: int = 32,
    seed: int = 2026082210,
) -> list[Trajectory]:
    """Load only the archive TEST split and select a reproducible balanced subset."""

    member = f"{name}_TEST.ts"
    values, labels = _read_ts_from_zip(archive_path, member)
    rng = np.random.default_rng(seed)
    selected: list[int] = []
    for label in sorted(set(labels.tolist())):
        candidates = np.flatnonzero(labels == label)
        take = min(per_class, len(candidates))
        selected.extend(rng.choice(candidates, size=take, replace=False).tolist())
    selected.sort()
    return [
        Trajectory(f"{name.lower()}-{index:04d}", name.lower(), int(labels[index]), values[index, :, None])
        for index in selected
    ]


def similarity_transform(values: np.ndarray, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    x = np.asarray(values, dtype=np.float64)
    channels = x.shape[1] if x.ndim == 2 else 1
    matrix = rng.normal(size=(channels, channels))
    q, _ = np.linalg.qr(matrix)
    scale = float(np.exp(rng.uniform(np.log(0.2), np.log(5.0))))
    translation = rng.normal(0.0, 10.0, size=(1, channels))
    return scale * np.atleast_2d(x).reshape(len(x), channels) @ q + translation


def noisy_transform(values: np.ndarray, seed: int, relative_sigma: float = 0.02) -> np.ndarray:
    rng = np.random.default_rng(seed)
    x = np.asarray(values, dtype=np.float64)
    scale = float(np.std(x))
    return x + rng.normal(0.0, relative_sigma * max(scale, 1.0e-9), size=x.shape)

