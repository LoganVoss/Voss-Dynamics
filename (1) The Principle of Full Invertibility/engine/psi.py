"""
Engine for the drift--kick--phase model analyzed in the thesis.

The implemented advance is

    Z(t + Δt) = Ψ(Z(t); A, λ)

with token state

    z_i = (r_i, v_i, θ_i, φ_i)

and total force

    F = F_compass + F_clock/interference + F_mutual + F_observation (optional).

Two integrators are provided:

* damped semi-implicit Euler — the practical map (γ ≈ 0.982)
* velocity Verlet — the ideal mechanical limit (γ = 1)

The translational inverse of the Euler map is algebraic. Phase inversion uses
a fixed-point iteration on the Kuramoto residual and is unique only in a
certified contraction chamber. The unrestricted default phase map can have
multiple preimages; the routine returns one attempted branch.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

TWOPI = 2.0 * np.pi


def wrap(a: np.ndarray | float) -> np.ndarray | float:
    return np.mod(a, TWOPI)


@dataclass
class FieldParams:
    alpha: float = 0.055
    gamma: float = 0.982
    beta: float = 0.08
    dt: float = 0.35
    eps: float = 1e-3
    r_min: float = 0.12
    r_max: float = 4.2
    k_wave: float = 2.4
    wave_amp: float = 0.12
    mutual_rep: float = 0.18
    mutual_att: float = 0.035
    att_soft: float = 0.4
    theta_drift: float = 0.0
    wall: float = 3.6
    wall_bounce: float = 0.0
    obs_strength: float = 0.04


@dataclass
class Field:
    r: np.ndarray
    v: np.ndarray
    theta: np.ndarray
    phi: np.ndarray
    omega: np.ndarray
    m: np.ndarray
    p: np.ndarray
    Q: np.ndarray
    theta_a: np.ndarray
    phi_a: np.ndarray
    params: FieldParams = field(default_factory=FieldParams)
    phrase_centers: Optional[np.ndarray] = None

    @property
    def n(self) -> int:
        return int(self.r.shape[0])

    def copy(self) -> "Field":
        return Field(
            r=self.r.copy(),
            v=self.v.copy(),
            theta=self.theta.copy(),
            phi=self.phi.copy(),
            omega=self.omega.copy(),
            m=self.m.copy(),
            p=self.p.copy(),
            Q=self.Q.copy(),
            theta_a=self.theta_a.copy(),
            phi_a=self.phi_a.copy(),
            params=self.params,
            phrase_centers=None
            if self.phrase_centers is None
            else self.phrase_centers.copy(),
        )

    def flatten(self) -> np.ndarray:
        return np.concatenate(
            [self.r.ravel(), self.v.ravel(), self.theta, self.phi]
        )

    def unflatten_into(self, x: np.ndarray) -> None:
        n = self.n
        self.r = x[: 2 * n].reshape(n, 2).copy()
        self.v = x[2 * n : 4 * n].reshape(n, 2).copy()
        self.theta = wrap(x[4 * n : 5 * n].copy())
        self.phi = wrap(x[5 * n : 6 * n].copy())


def default_anchors() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    p = np.array(
        [
            [0.0, 2.6],
            [2.6, 0.0],
            [0.0, -2.6],
            [-2.6, 0.0],
            [1.85, 1.85],
            [0.0, 2.0],
        ],
        dtype=float,
    )
    Q = np.array([2.2, 2.0, 1.6, 1.6, 2.4, 1.4], dtype=float)
    theta_a = np.array(
        [np.pi / 2, 0.0, -np.pi / 2, np.pi, np.pi / 4, np.pi / 2]
    )
    phi_a = np.array([0.0, np.pi / 2, np.pi, 3 * np.pi / 2, np.pi / 6, 0.0])
    return p, Q, theta_a, phi_a


def make_field(
    n: int = 16,
    seed: int = 0,
    params: Optional[FieldParams] = None,
    spread: float = 3.0,
) -> Field:
    rng = np.random.default_rng(seed)
    p, Q, theta_a, phi_a = default_anchors()
    return Field(
        r=rng.uniform(-spread, spread, size=(n, 2)),
        v=rng.normal(0.0, 0.05, size=(n, 2)),
        theta=rng.uniform(0.0, TWOPI, size=n),
        phi=rng.uniform(0.0, TWOPI, size=n),
        omega=rng.uniform(0.02, 0.08, size=n),
        m=np.ones(n),
        p=p,
        Q=Q,
        theta_a=theta_a,
        phi_a=phi_a,
        params=params or FieldParams(),
    )


def pairwise(r: np.ndarray, eps: float) -> tuple[np.ndarray, np.ndarray]:
    dvec = r[:, None, :] - r[None, :, :]
    dist = np.linalg.norm(dvec, axis=-1)
    dist = dist + np.eye(len(r)) * 1e9
    return dvec, dist + eps


def forces(field: Field, r: np.ndarray, theta: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """Total force at given (r, θ, φ). Force depends on position and phase, not velocity."""
    prm = field.params
    n = r.shape[0]
    F = np.zeros((n, 2))

    for a in range(len(field.Q)):
        delta = field.p[a] - r
        dist = np.linalg.norm(delta, axis=1) + prm.eps
        mask = (dist > prm.r_min) & (dist < prm.r_max)
        f_mod = 0.55 + 0.45 * np.cos(theta - field.theta_a[a])
        phase_mod = 0.6 + 0.4 * np.cos(phi - field.phi_a[a])
        strength = prm.alpha * field.Q[a] * field.m * f_mod * phase_mod / dist
        F[mask] += (delta[mask].T * strength[mask]).T

    dvec, dist = pairwise(r, prm.eps)
    # explicit interference scalar field U = Σ A_j/(d+ε) sin(k d + φ_j)
    # force contribution -∇U on each particle from every other source
    kd = prm.k_wave * dist
    s = np.sin(kd + phi[None, :])
    c = np.cos(kd + phi[None, :])
    u_amp = prm.wave_amp / dist
    # U_i = Σ_j A/d sin(k d + φ_j)
    # dU/dd = A[ -sin(kd+φ)/d² + k cos(kd+φ)/d ]
    # F_i = -∇_{r_i} U = -(dU/dd) (r_i - r_j)/d
    dU_dd = u_amp * (-s / dist + prm.k_wave * c)
    F += np.sum((-dU_dd)[:, :, None] * (dvec / dist[:, :, None]), axis=1)

    rep = prm.mutual_rep / (dist**2)
    att = prm.mutual_att / (dist + prm.att_soft)
    coeff = att - rep
    F += np.sum((-dvec) * coeff[:, :, None], axis=1)

    if field.phrase_centers is not None and prm.obs_strength != 0.0:
        for c in field.phrase_centers:
            delta = c - r
            distc = np.linalg.norm(delta, axis=1) + prm.eps
            F += prm.obs_strength * (delta.T / distc).T

    return F


def apply_walls(field: Field) -> None:
    prm = field.params
    if prm.wall <= 0 or prm.wall_bounce == 0.0:
        return
    for d in range(2):
        over = field.r[:, d] > prm.wall
        under = field.r[:, d] < -prm.wall
        field.r[over, d] = prm.wall
        field.v[over, d] *= -prm.wall_bounce
        field.r[under, d] = -prm.wall
        field.v[under, d] *= -prm.wall_bounce


def euler_step(field: Field) -> None:
    """Original practical map: semi-implicit Euler with damping γ."""
    prm = field.params
    field.r = field.r + field.v * prm.dt
    apply_walls(field)
    F = forces(field, field.r, field.theta, field.phi)
    field.v = prm.gamma * field.v + (prm.dt / field.m[:, None]) * F
    _advance_phases(field)


def _advance_phases(field: Field) -> None:
    prm = field.params
    _, dist = pairwise(field.r, prm.eps)
    # φ_i ← φ_i + ω_i Δt + β Σ_j sin(φ_j - φ_i) / (d_ij + ε)
    s = np.sin(field.phi[None, :] - field.phi[:, None])
    couple = np.sum(s / dist, axis=1)
    field.phi = wrap(field.phi + field.omega * prm.dt + prm.beta * couple * prm.dt)
    if prm.theta_drift != 0.0:
        field.theta = wrap(
            field.theta + prm.theta_drift * np.sin(field.phi - field.theta)
        )


def euler_inverse_step(
    field: Field, phase_iters: int = 200, phase_tolerance: float = 1e-15
) -> dict:
    """
    Attempt one inverse branch of an Euler step, assuming no wall collisions.

    Order of the forward map:
        r' = r + v Δt
        F  = F(r', θ, φ)          # old phases
        v' = γ v + (Δt/m) F
        φ' = φ + ω Δt + β K(φ; r')
        θ' = θ + η sin(φ' - θ)    # if drift enabled

    The translational recovery is exact once a valid old phase is selected.
    Fixed-point phase iteration is guaranteed only when its contraction bound
    is below one; outside that chamber the target may have multiple preimages.
    """
    prm = field.params
    r_new = field.r.copy()
    v_new = field.v.copy()
    theta_new = field.theta.copy()
    phi_new = field.phi.copy()

    # invert θ drift if present: θ' = θ + η sin(φ' - θ)
    theta_old = theta_new.copy()
    if prm.theta_drift != 0.0:
        for _ in range(20):
            residual = theta_new - theta_old - prm.theta_drift * np.sin(
                phi_new - theta_old
            )
            # d/dθ_old of residual = -1 - η * cos(φ'-θ_old) * (-1) = -1 + η cos(...)
            deriv = -1.0 + prm.theta_drift * np.cos(phi_new - theta_old)
            theta_old = theta_old - residual / deriv
        theta_old = wrap(theta_old)

    # invert Kuramoto: φ' = φ + ω Δt + β Δt K(φ; r')
    # fixed point: φ = φ' - ω Δt - β Δt K(φ; r')
    phi_old = wrap(phi_new - field.omega * prm.dt)
    _, dist = pairwise(r_new, prm.eps)
    last_err = math.inf
    phase_iterations = 0
    for phase_iterations in range(1, phase_iters + 1):
        s = np.sin(phi_old[None, :] - phi_old[:, None])
        couple = np.sum(s / dist, axis=1)
        phi_next = wrap(phi_new - field.omega * prm.dt - prm.beta * prm.dt * couple)
        phi_old = phi_next
        # Convergence is judged by the wrapped residual of the forward phase
        # equation, not merely by the difference between fixed-point iterates.
        s_check = np.sin(phi_old[None, :] - phi_old[:, None])
        couple_check = np.sum(s_check / dist, axis=1)
        predicted = wrap(
            phi_old
            + field.omega * prm.dt
            + prm.beta * prm.dt * couple_check
        )
        phase_residual = np.arctan2(
            np.sin(predicted - phi_new), np.cos(predicted - phi_new)
        )
        last_err = float(np.max(np.abs(phase_residual)))
        if last_err < phase_tolerance:
            break

    F = forces(field, r_new, theta_old, phi_old)
    if abs(prm.gamma) < 1e-15:
        raise ZeroDivisionError("γ = 0 destroys invertibility of the translational block")
    v_old = (v_new - (prm.dt / field.m[:, None]) * F) / prm.gamma
    r_old = r_new - v_old * prm.dt

    field.r = r_old
    field.v = v_old
    field.theta = theta_old
    field.phi = phi_old
    return {
        "phase_residual": last_err,
        "phase_iterations": phase_iterations,
        "phase_converged": last_err < phase_tolerance,
    }


def verlet_step(field: Field) -> None:
    """Ideal reversible limit: γ = 1, velocity Verlet on (r, v), then phases."""
    prm = field.params
    a = forces(field, field.r, field.theta, field.phi) / field.m[:, None]
    field.r = field.r + field.v * prm.dt + 0.5 * a * prm.dt**2
    a_new = forces(field, field.r, field.theta, field.phi) / field.m[:, None]
    field.v = field.v + 0.5 * (a + a_new) * prm.dt
    _advance_phases(field)


def verlet_inverse_step(
    field: Field, phase_iters: int = 200, phase_tolerance: float = 1e-15
) -> dict:
    """
    Reverse a Verlet step by inverting phases, negating time via
    velocity flip + forward Verlet + velocity flip — after phases
    have been inverted to the pre-step values.

    More direct: store is expensive, so we invert phases first, then
    use the standard Verlet time-reversal on (r, v) with those phases
    held to the values the forward step used. Because phases are
    updated AFTER (r, v), Verlet on (r, v) used (θ, φ)_old for both
    accelerations only approximately (a_new uses same old phases in
    our split). Our split Verlet uses old phases for a and a_new —
    wait: a_new = F(r_new, θ_old, φ_old)/m. Then phases advance.
    So (r, v) Verlet is closed given old phases. Inverse:
      invert phases to old
      then invert Verlet with those phases held fixed.
    """
    prm = field.params
    r_new = field.r.copy()
    v_new = field.v.copy()
    theta_new = field.theta.copy()
    phi_new = field.phi.copy()

    theta_old = theta_new.copy()
    if prm.theta_drift != 0.0:
        for _ in range(20):
            residual = theta_new - theta_old - prm.theta_drift * np.sin(
                phi_new - theta_old
            )
            deriv = -1.0 + prm.theta_drift * np.cos(phi_new - theta_old)
            theta_old = theta_old - residual / deriv
        theta_old = wrap(theta_old)

    phi_old = wrap(phi_new - field.omega * prm.dt)
    _, dist = pairwise(r_new, prm.eps)
    last_err = math.inf
    phase_iterations = 0
    for phase_iterations in range(1, phase_iters + 1):
        s = np.sin(phi_old[None, :] - phi_old[:, None])
        couple = np.sum(s / dist, axis=1)
        phi_next = wrap(phi_new - field.omega * prm.dt - prm.beta * prm.dt * couple)
        phi_old = phi_next
        s_check = np.sin(phi_old[None, :] - phi_old[:, None])
        couple_check = np.sum(s_check / dist, axis=1)
        predicted = wrap(
            phi_old
            + field.omega * prm.dt
            + prm.beta * prm.dt * couple_check
        )
        phase_residual = np.arctan2(
            np.sin(predicted - phi_new), np.cos(predicted - phi_new)
        )
        last_err = float(np.max(np.abs(phase_residual)))
        if last_err < phase_tolerance:
            break

    # Verlet inverse with phases held at old values:
    #   r' = r + v dt + 1/2 a(r,θ,φ) dt²
    #   v' = v + 1/2 (a(r)+a(r')) dt
    # Given r', v', θ, φ:
    #   a' = a(r',θ,φ)
    #   v_half = v' - 1/2 a' dt
    #   r = r' - v_half dt   ... then a = a(r), and v = v_half - 1/2 a dt
    # The position update used a(r), not a(r'), so r = r' - v dt - 1/2 a(r) dt²
    # is implicit. Iterate.
    a_new = forces(field, r_new, theta_old, phi_old) / field.m[:, None]
    v_half = v_new - 0.5 * a_new * prm.dt
    r_old = r_new - v_half * prm.dt
    for _ in range(12):
        a_old = forces(field, r_old, theta_old, phi_old) / field.m[:, None]
        r_old = r_new - v_half * prm.dt - 0.0 * a_old
        # exact: r' = r + v dt + 1/2 a(r) dt², v' = v + 1/2 (a+a') dt
        # v = v' - 1/2 (a+a') dt
        # r = r' - v dt - 1/2 a dt²
        v_old = v_new - 0.5 * (a_old + a_new) * prm.dt
        r_old = r_new - v_old * prm.dt - 0.5 * a_old * prm.dt**2
    a_old = forces(field, r_old, theta_old, phi_old) / field.m[:, None]
    v_old = v_new - 0.5 * (a_old + a_new) * prm.dt

    field.r = r_old
    field.v = v_old
    field.theta = theta_old
    field.phi = phi_old
    return {
        "phase_residual": last_err,
        "phase_iterations": phase_iterations,
        "phase_converged": last_err < phase_tolerance,
    }


def run(field: Field, steps: int, kind: str = "euler") -> None:
    step = euler_step if kind == "euler" else verlet_step
    for _ in range(steps):
        step(field)


def run_inverse(field: Field, steps: int, kind: str = "euler") -> list[float]:
    inv = euler_inverse_step if kind == "euler" else verlet_inverse_step
    residuals = []
    for _ in range(steps):
        residuals.append(inv(field)["phase_residual"])
    return residuals


def kinetic(field: Field) -> float:
    return float(0.5 * np.sum(field.m[:, None] * field.v**2))


def kuramoto_R(field: Field) -> float:
    return float(np.abs(np.mean(np.exp(1j * field.phi))))


def state_distance(a: Field, b: Field) -> dict:
    dr = float(np.max(np.linalg.norm(a.r - b.r, axis=1)))
    dv = float(np.max(np.linalg.norm(a.v - b.v, axis=1)))
    theta_delta = np.arctan2(
        np.sin(a.theta - b.theta), np.cos(a.theta - b.theta)
    )
    phi_delta = np.arctan2(np.sin(a.phi - b.phi), np.cos(a.phi - b.phi))
    dth = float(np.max(np.abs(theta_delta)))
    dph = float(np.max(np.abs(phi_delta)))
    return {
        "max_dr": dr,
        "max_dv": dv,
        "max_dtheta": dth,
        "max_dphi": dph,
        "l2": float(
            np.sqrt(
                np.sum((a.r - b.r) ** 2)
                + np.sum((a.v - b.v) ** 2)
                + np.sum(theta_delta**2)
                + np.sum(phi_delta**2)
            )
        ),
    }
