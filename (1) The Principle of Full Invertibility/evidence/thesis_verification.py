#!/usr/bin/env python3
"""
Supplemental self-consistency suite for The Principle of Full Invertibility.

This script supplements the baseline batteries with tests chosen specifically
to separate exact statements from numerical evidence:

1. determinant factorization of the complete (r, v, phi) Euler map;
2. the configuration-dependent contraction certificate for the phase inverse;
3. the exact two-clock phase-difference law;
4. replicated finite-precision recovery horizons;
5. discrete Shannon entropy under a permutation versus a many-to-one map;
6. an explicit warning that the visualization wall clamp is non-injective.

The output is a single JSON ledger suitable for archival with the manuscript.
"""

from __future__ import annotations

import json
import math
import platform
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from engine.discrete import force_1d, forward  # noqa: E402
from engine.psi import (  # noqa: E402
    FieldParams,
    euler_step,
    forces,
    make_field,
    pairwise,
    run,
    run_inverse,
    state_distance,
    wrap,
)


def angular_difference(a: float, b: float) -> float:
    return float(np.arctan2(np.sin(a - b), np.cos(a - b)))


def phase_jacobian(field, r: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """Jacobian D_phi G_r for one Kuramoto phase step."""
    prm = field.params
    _, dist = pairwise(r, prm.eps)
    n = field.n
    J = np.eye(n)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            c = math.cos(float(phi[j] - phi[i])) / float(dist[i, j])
            J[i, j] += prm.beta * prm.dt * c
            J[i, i] -= prm.beta * prm.dt * c
    return J


def phase_contraction_certificate(field, r: np.ndarray) -> float:
    """
    Sufficient infinity-norm Lipschitz bound for the fixed-point inverse.

        q = 2 |beta| dt max_i sum_{j != i} 1 / d_ij.

    q < 1 certifies a contraction. q >= 1 is inconclusive, not a proof of
    non-invertibility.
    """
    _, dist = pairwise(r, field.params.eps)
    reciprocal = 1.0 / dist
    np.fill_diagonal(reciprocal, 0.0)
    return float(
        2.0
        * abs(field.params.beta)
        * field.params.dt
        * np.max(np.sum(reciprocal, axis=1))
    )


def rvphi_vector(field) -> np.ndarray:
    return np.concatenate([field.r.ravel(), field.v.ravel(), field.phi])


def set_rvphi(field, x: np.ndarray) -> None:
    n = field.n
    field.r = x[: 2 * n].reshape(n, 2).copy()
    field.v = x[2 * n : 4 * n].reshape(n, 2).copy()
    field.phi = x[4 * n : 5 * n].copy()


def determinant_factorization() -> dict:
    """
    Check

        det D Psi = gamma^(2N) det D_phi G_r

    for the complete (r, v, phi) step with theta held fixed. The factorization
    is exact away from force cutoffs and phase wrapping boundaries.
    """
    rows = []
    epsilon = 2e-7
    for n in (2, 3, 4):
        for seed in range(5):
            prm = FieldParams(
                gamma=0.982,
                beta=0.08,
                dt=0.21,
                theta_drift=0.0,
                wall_bounce=0.0,
            )
            field = make_field(n=n, seed=100 + seed + 10 * n, params=prm, spread=1.2)
            # Keep the finite-difference chart away from the 2pi branch cut.
            field.phi = np.linspace(0.5, 2.2, n) + 0.01 * seed
            x0 = rvphi_vector(field)
            base = field.copy()
            euler_step(base)
            y0 = rvphi_vector(base)
            Jfd = np.zeros((x0.size, x0.size))
            for k in range(x0.size):
                plus = field.copy()
                minus = field.copy()
                xp = x0.copy()
                xm = x0.copy()
                xp[k] += epsilon
                xm[k] -= epsilon
                set_rvphi(plus, xp)
                set_rvphi(minus, xm)
                euler_step(plus)
                euler_step(minus)
                Jfd[:, k] = (rvphi_vector(plus) - rvphi_vector(minus)) / (
                    2.0 * epsilon
                )

            r_prime = field.r + field.v * prm.dt
            Jphase = phase_jacobian(field, r_prime, field.phi)
            det_fd = float(np.linalg.det(Jfd))
            det_phase = float(np.linalg.det(Jphase))
            prediction = float((prm.gamma ** (2 * n)) * det_phase)
            rel = abs(det_fd - prediction) / max(abs(prediction), 1e-15)
            rows.append(
                {
                    "n": n,
                    "seed": seed,
                    "det_full_finite_difference": det_fd,
                    "det_phase": det_phase,
                    "gamma_pow_2n": float(prm.gamma ** (2 * n)),
                    "factorized_prediction": prediction,
                    "relative_error": rel,
                    "phase_contraction_q_bound": phase_contraction_certificate(
                        field, r_prime
                    ),
                    "unused_base_norm": float(np.linalg.norm(y0)),
                }
            )
    return {
        "identity": "det D(Psi) = gamma^(2N) det D_phi(G_r)",
        "rows": rows,
        "max_relative_error": max(row["relative_error"] for row in rows),
        "interpretation": (
            "The positional force Jacobian and phase-dependence of the force "
            "cancel from the determinant factorization. The phase block and "
            "gamma are the independent local singularity channels."
        ),
    }


def two_clock_law() -> dict:
    """
    Verify the exact reduced equation

      delta' = delta + (omega_1-omega_2)dt
               - 2 beta dt sin(delta)/D

    where delta = phi_1 - phi_2 and D is the softened separation.
    """
    rng = np.random.default_rng(20260822)
    residuals = []
    cases = []
    for k in range(200):
        prm = FieldParams(
            gamma=1.0,
            beta=float(rng.uniform(0.01, 0.6)),
            dt=float(rng.uniform(0.05, 0.4)),
            alpha=0.0,
            wave_amp=0.0,
            mutual_rep=0.0,
            mutual_att=0.0,
            theta_drift=0.0,
            wall_bounce=0.0,
        )
        field = make_field(n=2, seed=k, params=prm)
        distance = float(rng.uniform(0.2, 40.0))
        field.r = np.array([[-distance / 2.0, 0.0], [distance / 2.0, 0.0]])
        field.v[:] = 0.0
        field.phi = rng.uniform(0.2, 5.8, size=2)
        field.omega = rng.uniform(0.01, 0.1, size=2)

        delta = angular_difference(float(field.phi[0]), float(field.phi[1]))
        _, dist = pairwise(field.r, prm.eps)
        D = float(dist[0, 1])
        predicted = (
            delta
            + float(field.omega[0] - field.omega[1]) * prm.dt
            - 2.0 * prm.beta * prm.dt * math.sin(delta) / D
        )
        euler_step(field)
        actual = angular_difference(float(field.phi[0]), float(field.phi[1]))
        residual = abs(angular_difference(actual, predicted))
        residuals.append(residual)
        if k < 5:
            cases.append(
                {
                    "distance_softened": D,
                    "beta": prm.beta,
                    "dt": prm.dt,
                    "delta_before": delta,
                    "delta_predicted": predicted,
                    "delta_actual": actual,
                    "angular_residual": residual,
                }
            )

    # Exact synchronized manifold: equal omega and delta = 0 remain equal.
    sync_errors = []
    for distance in (0.24, 1.0, 30.0, 1e6):
        prm = FieldParams(
            gamma=1.0,
            beta=0.45,
            dt=0.25,
            alpha=0.0,
            wave_amp=0.0,
            mutual_rep=0.0,
            mutual_att=0.0,
            theta_drift=0.0,
            wall_bounce=0.0,
        )
        field = make_field(n=2, seed=9, params=prm)
        frozen_r = np.array([[-distance / 2.0, 0.0], [distance / 2.0, 0.0]])
        field.r = frozen_r.copy()
        field.v[:] = 0.0
        field.phi[:] = 1.234
        field.omega[:] = 0.05
        max_delta = 0.0
        for _ in range(1000):
            euler_step(field)
            field.r = frozen_r.copy()
            field.v[:] = 0.0
            max_delta = max(
                max_delta,
                abs(
                    angular_difference(
                        float(field.phi[0]), float(field.phi[1])
                    )
                ),
            )
        sync_errors.append({"distance": distance, "max_abs_delta": max_delta})

    return {
        "reduced_law": (
            "delta' = delta + Delta_omega*dt "
            "- (2*beta*dt/D)*sin(delta) (mod 2pi)"
        ),
        "max_law_residual": max(residuals),
        "sample_cases": cases,
        "equal_omega_synchronized_manifold": sync_errors,
        "two_clock_global_condition": (
            "For signed a = 2*beta*dt/D with |a| < 1, the lifted delta map "
            "is an orientation-preserving circle diffeomorphism. For |a| > 1 "
            "it is non-monotone and non-injective. Synchronization is "
            "hyperbolically stable for 0 < a < 2, nonhyperbolically stable at "
            "a = 2, and unstable for a > 2."
        ),
    }


def replicated_information_horizon() -> dict:
    rows = []
    seeds = list(range(21, 33))
    for gamma in (1.0, 0.982):
        for steps in (40, 80, 100, 120, 160, 200):
            errors = []
            samples = []
            for seed in seeds:
                prm = FieldParams(
                    gamma=gamma,
                    beta=0.08,
                    dt=0.35,
                    theta_drift=0.0,
                    wall_bounce=0.0,
                )
                field = make_field(n=8, seed=seed, params=prm)
                origin = field.copy()
                max_q = 0.0
                first_uncertified_step = None
                for step in range(1, steps + 1):
                    euler_step(field)
                    q = phase_contraction_certificate(field, field.r)
                    max_q = max(max_q, q)
                    if q >= 1.0 and first_uncertified_step is None:
                        first_uncertified_step = step
                inverse_residuals = run_inverse(field, steps, kind="euler")
                max_inverse_residual = max(inverse_residuals)
                error = state_distance(origin, field)["l2"]
                errors.append(error)
                samples.append(
                    {
                        "seed": seed,
                        "l2": error,
                        "max_phase_contraction_q": max_q,
                        "contraction_certified": max_q < 1.0,
                        "max_inverse_phase_residual": max_inverse_residual,
                        "inverse_solver_converged": max_inverse_residual < 1e-12,
                        "roundtrip_certified": (
                            max_q < 1.0 and max_inverse_residual < 1e-12
                        ),
                        "first_uncertified_step": first_uncertified_step,
                    }
                )
            e = np.asarray(errors, dtype=float)
            logs = np.log10(np.maximum(e, np.finfo(float).tiny))
            certified_samples = [
                sample for sample in samples if sample["roundtrip_certified"]
            ]
            certified_e = np.asarray(
                [sample["l2"] for sample in certified_samples], dtype=float
            )
            certified_logs = np.log10(
                np.maximum(certified_e, np.finfo(float).tiny)
            )
            rows.append(
                {
                    "gamma": gamma,
                    "steps": steps,
                    "seeds": len(seeds),
                    "median_l2": float(np.median(e)),
                    "p10_l2": float(np.quantile(e, 0.10)),
                    "p90_l2": float(np.quantile(e, 0.90)),
                    "max_l2": float(np.max(e)),
                    "median_log10_l2": float(np.median(logs)),
                    "p10_log10_l2": float(np.quantile(logs, 0.10)),
                    "p90_log10_l2": float(np.quantile(logs, 0.90)),
                    "fraction_below_1e6": float(np.mean(e < 1e-6)),
                    "fraction_below_1e4": float(np.mean(e < 1e-4)),
                    "certified_seed_count": len(certified_samples),
                    "certified_seeds": [
                        sample["seed"] for sample in certified_samples
                    ],
                    "certified_median_l2": float(np.median(certified_e)),
                    "certified_p10_l2": float(np.quantile(certified_e, 0.10)),
                    "certified_p90_l2": float(np.quantile(certified_e, 0.90)),
                    "certified_max_l2": float(np.max(certified_e)),
                    "certified_median_log10_l2": float(
                        np.median(certified_logs)
                    ),
                    "certified_p10_log10_l2": float(
                        np.quantile(certified_logs, 0.10)
                    ),
                    "certified_p90_log10_l2": float(
                        np.quantile(certified_logs, 0.90)
                    ),
                    "certified_fraction_below_1e6": float(
                        np.mean(certified_e < 1e-6)
                    ),
                    "certified_fraction_below_1e4": float(
                        np.mean(certified_e < 1e-4)
                    ),
                    "samples": samples,
                }
            )
    return {
        "n_tokens": 8,
        "seeds": seeds,
        "rows": rows,
        "interpretation": (
            "Round-trip summaries include only trajectories with max q < 1 "
            "at every step and wrapped forward phase residual below 1e-12. "
            "They measure float64 round-trip recovery, not an isolated "
            "condition number. Uncertified or unconverged trajectories are "
            "reported separately."
        ),
    }


def euler_inverse_newton_step(field, max_iterations: int = 30) -> dict:
    """
    Invert one step with Newton's method for the phase block.

    This is a diagnostic implementation, not a replacement for the production
    inverse. Convergence recovers one phase branch; it does not certify that the
    branch is unique.
    """
    prm = field.params
    r_new = field.r.copy()
    v_new = field.v.copy()
    phi_target = field.phi.copy()
    phi_old = wrap(phi_target - field.omega * prm.dt)
    _, dist = pairwise(r_new, prm.eps)
    residual_norm = math.inf
    condition = math.inf
    for iteration in range(max_iterations):
        s = np.sin(phi_old[None, :] - phi_old[:, None])
        couple = np.sum(s / dist, axis=1)
        predicted = phi_old + field.omega * prm.dt + prm.beta * prm.dt * couple
        residual = np.arctan2(
            np.sin(predicted - phi_target), np.cos(predicted - phi_target)
        )
        J = phase_jacobian(field, r_new, phi_old)
        condition = float(np.linalg.cond(J))
        step = np.linalg.solve(J, residual)
        phi_old = wrap(phi_old - step)
        residual_norm = float(np.max(np.abs(residual)))
        if residual_norm < 1e-13:
            break

    F = forces(field, r_new, field.theta, phi_old)
    if abs(prm.gamma) < 1e-15:
        raise ZeroDivisionError("gamma = 0")
    v_old = (v_new - (prm.dt / field.m[:, None]) * F) / prm.gamma
    r_old = r_new - v_old * prm.dt
    field.r = r_old
    field.v = v_old
    field.phi = phi_old
    return {
        "phase_residual": residual_norm,
        "phase_jacobian_condition": condition,
        "iterations": iteration + 1,
    }


def newton_phase_rescue() -> dict:
    rows = []
    for beta, steps, n, seed in (
        (0.08, 40, 8, 22),
        (0.08, 200, 8, 22),
        (0.40, 30, 10, 6),
        (0.50, 30, 10, 6),
    ):
        prm = FieldParams(
            gamma=0.982,
            beta=beta,
            dt=0.35,
            theta_drift=0.0,
            wall_bounce=0.0,
        )
        field = make_field(n=n, seed=seed, params=prm)
        origin = field.copy()
        run(field, steps, kind="euler")

        fixed = field.copy()
        fixed_residuals = run_inverse(fixed, steps, kind="euler")

        newton = field.copy()
        newton_diagnostics = []
        for _ in range(steps):
            newton_diagnostics.append(euler_inverse_newton_step(newton))

        rows.append(
            {
                "beta": beta,
                "steps": steps,
                "n": n,
                "seed": seed,
                "fixed_point_l2": state_distance(origin, fixed)["l2"],
                "newton_l2": state_distance(origin, newton)["l2"],
                "fixed_point_max_reported_residual": max(fixed_residuals),
                "newton_max_phase_residual": max(
                    d["phase_residual"] for d in newton_diagnostics
                ),
                "newton_max_phase_jacobian_condition": max(
                    d["phase_jacobian_condition"] for d in newton_diagnostics
                ),
            }
        )
    return {
        "rows": rows,
        "interpretation": (
            "Newton can recover a locally valid branch when the fixed-point "
            "iteration fails. This does not establish global injectivity. "
            "For seed 22, Newton recovers the generating branch even though "
            "two additional complete preimages exist."
        ),
    }


def phase_multiroot_collision() -> dict:
    """
    Construct three complete preimages of the seed-22, step-9 target.

    This is a direct counterexample to global invertibility of the unrestricted
    default phase map.
    """
    prm = FieldParams(
        gamma=0.982,
        beta=0.08,
        dt=0.35,
        theta_drift=0.0,
        wall_bounce=0.0,
    )
    field = make_field(n=8, seed=22, params=prm)
    for _ in range(8):
        euler_step(field)
    generating_prestate = field.copy()
    euler_step(field)
    target = field.copy()
    _, dist = pairwise(target.r, prm.eps)

    starts = [
        generating_prestate.phi.copy(),
        np.array(
            [
                3.803563,
                6.181529,
                5.802384,
                4.779631,
                4.879962,
                6.065629,
                4.973750,
                3.681332,
            ]
        ),
        np.array(
            [
                4.950865,
                6.181770,
                5.802543,
                4.779626,
                4.880205,
                6.065831,
                3.825819,
                3.681122,
            ]
        ),
    ]

    def refine(phi):
        phi = phi.copy()
        residual_norm = math.inf
        for _ in range(40):
            s = np.sin(phi[None, :] - phi[:, None])
            couple = np.sum(s / dist, axis=1)
            predicted = phi + target.omega * prm.dt + prm.beta * prm.dt * couple
            residual = np.arctan2(
                np.sin(predicted - target.phi),
                np.cos(predicted - target.phi),
            )
            residual_norm = float(np.max(np.abs(residual)))
            if residual_norm < 1e-13:
                break
            J = phase_jacobian(target, target.r, phi)
            phi = wrap(phi - np.linalg.solve(J, residual))
        return phi, residual_norm

    def complete_prestate(phi):
        candidate = target.copy()
        candidate.theta = target.theta.copy()
        candidate.phi = phi.copy()
        F = forces(candidate, target.r, candidate.theta, phi)
        candidate.v = (
            target.v - (prm.dt / candidate.m[:, None]) * F
        ) / prm.gamma
        candidate.r = target.r - candidate.v * prm.dt
        return candidate

    roots = []
    for index, start in enumerate(starts):
        phi, residual = refine(start)
        candidate = complete_prestate(phi)
        output = candidate.copy()
        euler_step(output)
        roots.append(
            {
                "root": index,
                "phi": phi.tolist(),
                "phase_equation_residual": residual,
                "det_phase_jacobian": float(
                    np.linalg.det(phase_jacobian(target, target.r, phi))
                ),
                "prestate_distance_from_generating": state_distance(
                    candidate, generating_prestate
                )["l2"],
                "forward_collision_residual": state_distance(output, target)[
                    "l2"
                ],
            }
        )

    pairwise_root_distances = []
    refined_phases = [np.asarray(root["phi"]) for root in roots]
    for i in range(len(refined_phases)):
        for j in range(i + 1, len(refined_phases)):
            delta = np.arctan2(
                np.sin(refined_phases[i] - refined_phases[j]),
                np.cos(refined_phases[i] - refined_phases[j]),
            )
            pairwise_root_distances.append(
                {"root_i": i, "root_j": j, "phase_l2": float(np.linalg.norm(delta))}
            )

    return {
        "seed": 22,
        "target_step": 9,
        "q_at_target": phase_contraction_certificate(target, target.r),
        "root_count_constructed": len(roots),
        "roots": roots,
        "pairwise_root_phase_distances": pairwise_root_distances,
        "interpretation": (
            "The unrestricted default map is not globally injective. "
            "Each phase root, paired with its algebraically recovered "
            "velocity and position, is a distinct complete prestate with the "
            "same target."
        ),
    }


def shannon_entropy(probabilities: np.ndarray) -> float:
    p = probabilities[probabilities > 0]
    return float(-np.sum(p * np.log2(p)))


def finite_entropy_contrast() -> dict:
    m = 31
    rng = np.random.default_rng(31)
    weights = rng.random(m * m)
    p = weights / weights.sum()

    image = np.empty(m * m, dtype=int)
    for r in range(m):
        for v in range(m):
            rp, vp = forward(r, v, m, 3)
            image[r * m + v] = rp * m + vp
    pushed = np.zeros_like(p)
    for source, target in enumerate(image):
        pushed[target] += p[source]

    uniform = np.full(m * m, 1.0 / (m * m))
    collapsed = np.zeros(m * m)
    for r in range(m):
        for v in range(m):
            rp = (r + v) % m
            vp = force_1d(rp, m)
            collapsed[rp * m + vp] += uniform[r * m + v]

    return {
        "M": m,
        "bijective_g": 3,
        "random_distribution_entropy_before_bits": shannon_entropy(p),
        "random_distribution_entropy_after_bijection_bits": shannon_entropy(pushed),
        "bijection_entropy_abs_difference_bits": abs(
            shannon_entropy(p) - shannon_entropy(pushed)
        ),
        "uniform_entropy_before_g0_bits": shannon_entropy(uniform),
        "uniform_entropy_after_g0_bits": shannon_entropy(collapsed),
        "g0_unique_images": int(np.count_nonzero(collapsed)),
        "interpretation": (
            "Discrete Shannon entropy is invariant under relabeling by a "
            "bijection. The g=0 map aggregates probability masses and lowers "
            "the uniform entropy by log2(M) bits."
        ),
    }


def fixed_partition_marginal_entropy() -> dict:
    """
    Sum one-dimensional marginal entropies under one fixed quantizer.

    This is an observer-dependent diagnostic, not a joint-state entropy.
    Unlike the legacy battery, bin edges are declared once and held fixed
    across coordinates, ensembles, and times.
    """
    position_edges = np.array(
        [
            -np.inf,
            -20,
            -10,
            -5,
            -3,
            -2,
            -1,
            0,
            1,
            2,
            3,
            5,
            10,
            20,
            np.inf,
        ],
        dtype=float,
    )
    velocity_edges = np.array(
        [
            -np.inf,
            -20,
            -10,
            -5,
            -2,
            -1,
            -0.5,
            0,
            0.5,
            1,
            2,
            5,
            10,
            20,
            np.inf,
        ],
        dtype=float,
    )
    angle_edges = np.linspace(0.0, 2.0 * np.pi, 13)

    def marginal_sum(samples, edges):
        total = 0.0
        for coordinate in range(samples.shape[1]):
            counts, _ = np.histogram(samples[:, coordinate], bins=edges)
            probabilities = counts[counts > 0].astype(float) / samples.shape[0]
            total += shannon_entropy(probabilities)
        return total

    rows = []
    ensemble = 120
    n = 6
    times = (0, 10, 30, 80, 150)
    for gamma in (1.0, 0.982):
        prm = FieldParams(gamma=gamma, theta_drift=0.0, wall_bounce=0.0)
        fields = [make_field(n=n, seed=200 + i, params=prm) for i in range(ensemble)]
        previous = 0
        for target in times:
            if target > previous:
                for field in fields:
                    run(field, target - previous, kind="euler")
            r = np.stack([field.r.ravel() for field in fields])
            v = np.stack([field.v.ravel() for field in fields])
            theta = np.stack([field.theta for field in fields])
            phi = np.stack([field.phi for field in fields])
            h_r = marginal_sum(r, position_edges)
            h_v = marginal_sum(v, velocity_edges)
            h_theta = marginal_sum(theta, angle_edges)
            h_phi = marginal_sum(phi, angle_edges)
            rows.append(
                {
                    "gamma": gamma,
                    "t": target,
                    "position_marginal_sum_bits": h_r,
                    "velocity_marginal_sum_bits": h_v,
                    "theta_marginal_sum_bits": h_theta,
                    "phase_marginal_sum_bits": h_phi,
                    "full_marginal_sum_bits": h_r + h_v + h_theta + h_phi,
                }
            )
            previous = target
    return {
        "ensemble": ensemble,
        "n_tokens": n,
        "position_edges": [
            "-inf" if np.isneginf(value) else "inf" if np.isposinf(value) else value
            for value in position_edges.tolist()
        ],
        "velocity_edges": [
            "-inf" if np.isneginf(value) else "inf" if np.isposinf(value) else value
            for value in velocity_edges.tolist()
        ],
        "angle_edges": angle_edges.tolist(),
        "rows": rows,
        "interpretation": (
            "These are sums of one-dimensional marginal entropies under a "
            "fixed observer partition. They are not joint entropies and may "
            "rise or fall without complete-state collisions."
        ),
    }


def decoupled_clock_memory_and_readout() -> dict:
    """Test exact post-decoupling memory and controlled local reset effects."""
    lock_params = FieldParams(
        gamma=1.0,
        beta=0.45,
        dt=0.25,
        alpha=0.0,
        wave_amp=0.0,
        mutual_rep=0.0,
        mutual_att=0.0,
        theta_drift=0.0,
        wall_bounce=0.0,
    )
    locked = make_field(n=2, seed=2, params=lock_params)
    hold_near = np.array([[-0.12, 0.0], [0.12, 0.0]])
    locked.r = hold_near.copy()
    locked.v[:] = 0.0
    locked.phi = np.array([0.4, 2.9])
    locked.omega[:] = 0.05
    # Stop before machine-exact synchrony to preserve a nonzero remembered
    # phase difference.
    for _ in range(4):
        euler_step(locked)
        locked.r = hold_near.copy()
        locked.v[:] = 0.0
    remembered_delta = angular_difference(
        float(locked.phi[0]), float(locked.phi[1])
    )

    def post_decoupling(equal_frequency):
        field = locked.copy()
        field.params = replace(lock_params, beta=0.0)
        field.r = np.array([[-15.0, 0.0], [15.0, 0.0]])
        field.v[:] = 0.0
        if equal_frequency:
            field.omega[:] = 0.05
        else:
            field.omega = np.array([0.03, 0.08])
        initial = angular_difference(float(field.phi[0]), float(field.phi[1]))
        trace = []
        for step in range(1, 301):
            euler_step(field)
            field.r = np.array([[-15.0, 0.0], [15.0, 0.0]])
            field.v[:] = 0.0
            if step in (1, 25, 50, 100, 200, 300):
                actual = angular_difference(
                    float(field.phi[0]), float(field.phi[1])
                )
                expected = angular_difference(
                    initial + step * float(field.omega[0] - field.omega[1])
                    * field.params.dt,
                    0.0,
                )
                trace.append(
                    {
                        "t": step,
                        "delta_actual": actual,
                        "delta_expected": expected,
                        "angular_residual": abs(
                            angular_difference(actual, expected)
                        ),
                    }
                )
        return trace

    # Controlled reset: clone one phase state and alter only separation.
    reset_base = locked.copy()
    reset_base.params = replace(lock_params, beta=0.4)
    reset_rows = []
    for distance in (0.3, 1.0, 3.0, 24.0):
        state = reset_base.copy()
        state.r = np.array(
            [[-distance / 2.0, 0.0], [distance / 2.0, 0.0]]
        )
        state.v[:] = 0.0
        control = state.copy()
        reset_state = state.copy()
        reset_state.phi[0] = 0.0
        euler_step(control)
        euler_step(reset_state)
        effect = abs(
            angular_difference(
                float(reset_state.phi[1]), float(control.phi[1])
            )
        )
        softened_distance = distance + state.params.eps
        predicted = abs(
            state.params.beta
            * state.params.dt
            / softened_distance
            * (
                math.sin(-float(state.phi[1]))
                - math.sin(float(state.phi[0] - state.phi[1]))
            )
        )
        reset_rows.append(
            {
                "distance": distance,
                "softened_distance": softened_distance,
                "reset_vs_control_effect": effect,
                "analytic_prediction": predicted,
                "absolute_residual": abs(effect - predicted),
            }
        )

    equal_trace = post_decoupling(True)
    unequal_trace = post_decoupling(False)
    return {
        "remembered_delta_after_4_lock_steps": remembered_delta,
        "coupling_after_lock": 0.0,
        "equal_omega_trace": equal_trace,
        "unequal_omega_trace": unequal_trace,
        "max_equal_omega_residual": max(
            row["angular_residual"] for row in equal_trace
        ),
        "max_unequal_omega_residual": max(
            row["angular_residual"] for row in unequal_trace
        ),
        "controlled_reset_rows": reset_rows,
        "interpretation": (
            "With beta set exactly to zero, equal-frequency phase difference "
            "is conserved and unequal-frequency drift follows the exact "
            "linear law. Reset effects are measured against an identical "
            "no-reset control with only distance changed."
        ),
    }


def wall_collision() -> dict:
    """
    Construct two distinct pre-states collapsed by the optional wall clamp.

    This operation belongs to the visualization/bounded experiments, not to
    the invertible core theorem.
    """
    prm = FieldParams(
        gamma=0.982,
        beta=0.08,
        dt=0.35,
        alpha=0.0,
        wave_amp=0.0,
        mutual_rep=0.0,
        mutual_att=0.0,
        theta_drift=0.0,
        wall=3.6,
        wall_bounce=0.4,
    )
    a = make_field(n=1, seed=1, params=prm)
    b = a.copy()
    # Both prestates begin inside the declared wall and overshoot it during the
    # drift. The clamp discards their different overshoot distances.
    a.r[:] = np.array([[3.0, 0.0]])
    b.r[:] = np.array([[3.2, 0.0]])
    a.v[:] = np.array([[3.0, 0.0]])
    b.v[:] = np.array([[3.0, 0.0]])
    before = state_distance(a, b)
    euler_step(a)
    euler_step(b)
    after = state_distance(a, b)
    return {
        "before_l2": before["l2"],
        "after_l2": after["l2"],
        "collision": after["l2"] < 1e-12,
        "wall": prm.wall,
        "wall_bounce": prm.wall_bounce,
        "interpretation": (
            "Clamping all overshoots to the same wall coordinate is "
            "many-to-one. Invertibility claims must exclude wall collisions "
            "or replace the clamp with an invertible boundary rule."
        ),
    }


def main() -> None:
    determinant = determinant_factorization()
    two_clock = two_clock_law()
    horizon = replicated_information_horizon()
    newton = newton_phase_rescue()
    multiroot = phase_multiroot_collision()
    finite_entropy = finite_entropy_contrast()
    fixed_entropy = fixed_partition_marginal_entropy()
    corrected_clock = decoupled_clock_memory_and_readout()
    wall = wall_collision()

    assert determinant["max_relative_error"] < 1e-6
    assert two_clock["max_law_residual"] < 1e-12
    assert multiroot["root_count_constructed"] >= 3
    assert max(
        root["forward_collision_residual"] for root in multiroot["roots"]
    ) < 1e-10
    assert min(
        root["prestate_distance_from_generating"]
        for root in multiroot["roots"][1:]
    ) > 1.0
    assert corrected_clock["max_equal_omega_residual"] < 1e-12
    assert corrected_clock["max_unequal_omega_residual"] < 1e-12
    assert max(
        row["absolute_residual"]
        for row in corrected_clock["controlled_reset_rows"]
    ) < 1e-12
    assert wall["collision"]

    ledger = {
        "metadata": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "purpose": "verification supplement for The Principle of Full Invertibility",
        },
        "determinant_factorization": determinant,
        "two_clock_law": two_clock,
        "replicated_information_horizon": horizon,
        "newton_phase_rescue": newton,
        "phase_multiroot_collision": multiroot,
        "finite_entropy_contrast": finite_entropy,
        "fixed_partition_marginal_entropy": fixed_entropy,
        "decoupled_clock_memory_and_readout": corrected_clock,
        "wall_collision": wall,
    }
    output = ROOT / "evidence" / "thesis_verification.json"
    output.write_text(json.dumps(ledger, indent=2))
    print(json.dumps(ledger, indent=2))
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
