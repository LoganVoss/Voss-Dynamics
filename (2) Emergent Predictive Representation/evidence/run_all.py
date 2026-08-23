#!/usr/bin/env python3
"""Run every analytic, numerical, and figure-producing thesis check.

All random generators use fixed seeds.  Synthetic target data are generated
only to test identifiability, freezing, scoring, and power; they are never
presented as laboratory evidence.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from scipy.linalg import expm
from scipy.stats import chi2, norm

from evidence.epr_core import (
    I2,
    PAULIS,
    SIGMA_X,
    SIGMA_Z,
    TargetConfig,
    binary_gaussian_mutual_information,
    bloch_density,
    born_probability,
    candidate_fibre_variation,
    candidate_global_fibre_variation,
    chsh_value,
    coherent_memory_homodyne_tv,
    compatible_peirce_degrees,
    density_from_spinor,
    epr_record_drift,
    equal_variance_gaussian_tv,
    fit_frozen_epsilon,
    geometric_fibre_shift,
    hopf_map,
    qnd_coordinates,
    qnd_stabilizer_dimension,
    relational_density,
    shots_per_orientation,
    simulate_target_records,
    singlet_joint_probability,
    spinor,
    target_rate_difference,
    zeno_excited_probability,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "evidence" / "outputs"
FIGURE_DIR = ROOT / "thesis" / "figures"
SEED = 20260821

COLORS = {
    "ink": "#18222c",
    "blue": "#1d5f8a",
    "orange": "#b95f2a",
    "green": "#3f7356",
    "red": "#9f3d45",
    "gray": "#70777d",
    "light": "#d8dde1",
}


def configure_plots() -> None:
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["STIX Two Text", "Times New Roman", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "font.size": 9.0,
            "axes.titlesize": 10.0,
            "axes.labelsize": 9.0,
            "legend.fontsize": 8.0,
            "xtick.labelsize": 8.0,
            "ytick.labelsize": 8.0,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "axes.prop_cycle": mpl.cycler(
                color=[
                    COLORS["blue"],
                    COLORS["orange"],
                    COLORS["green"],
                    COLORS["red"],
                    COLORS["gray"],
                ]
            ),
        }
    )


def save_figure(fig: mpl.figure.Figure, stem: str) -> None:
    for suffix in ("pdf", "png"):
        fig.savefig(FIGURE_DIR / f"{stem}.{suffix}", facecolor="white")
    plt.close(fig)


def symbolic_qnd_checks() -> dict[str, str | float]:
    z = sp.symbols("z", real=True)
    a, phi = sp.symbols("a phi", positive=True, finite=True)
    bz = a * sp.cos(phi) * (1 - z**2)
    btheta = a * sp.sin(phi)
    mutheta = a**2 * sp.sin(phi) * sp.cos(phi) * z

    gprime = sp.simplify(-btheta / bz)
    gpp = sp.diff(gprime, z)
    noise_residual = sp.simplify(btheta + gprime * bz)
    drift_residual = sp.simplify(
        mutheta + sp.Rational(1, 2) * gpp * bz**2
    )

    fprime = 1 / (1 - z**2)
    fpp = sp.diff(fprime, z)
    quotient_diffusion = sp.simplify(fprime * bz - a * sp.cos(phi))
    quotient_drift = sp.simplify(
        sp.Rational(1, 2) * fpp * bz**2
        - a**2 * sp.cos(phi) ** 2 * z
    )

    return {
        "invariant_noise_residual": str(noise_residual),
        "invariant_drift_residual": str(drift_residual),
        "quotient_diffusion_residual": str(quotient_diffusion),
        "quotient_drift_residual": str(quotient_drift),
        "invariant_gradient": str(gprime),
        "quotient_gradient": str(sp.simplify(fprime)),
    }


def verify_qm_limits() -> tuple[dict[str, float], dict[str, np.ndarray]]:
    rng = np.random.default_rng(SEED)

    born_trace_errors = []
    born_formula_values = []
    born_trace_values = []
    for _ in range(2_000):
        r = rng.normal(size=3)
        r *= rng.random() ** (1.0 / 3.0) / np.linalg.norm(r)
        axis = rng.normal(size=3)
        axis /= np.linalg.norm(axis)
        rho = bloch_density(r)
        effect = 0.5 * (I2 + np.tensordot(axis, PAULIS, axes=1))
        p_formula = born_probability(r, axis)
        p_trace = float(np.real(np.trace(rho @ effect)))
        born_formula_values.append(p_formula)
        born_trace_values.append(p_trace)
        born_trace_errors.append(abs(p_formula - p_trace))

    unitary_errors = []
    for _ in range(300):
        r = rng.normal(size=3)
        r *= rng.random() ** (1.0 / 3.0) / np.linalg.norm(r)
        rho = bloch_density(r)
        axis = rng.normal(size=3)
        axis /= np.linalg.norm(axis)
        angle = rng.uniform(-np.pi, np.pi)
        hamiltonian = 0.5 * np.tensordot(axis, PAULIS, axes=1)
        unitary = expm(-1.0j * angle * hamiltonian)
        evolved = unitary @ rho @ unitary.conj().T
        unitary_errors.append(
            max(
                abs(np.trace(evolved) - 1.0),
                abs(np.linalg.eigvalsh(evolved).min() - np.linalg.eigvalsh(rho).min()),
            )
        )

    times = np.linspace(0.0, 5.0, 301)
    gamma_phi = 0.7
    coherence_numeric = []
    for t in times:
        generator = np.diag([-gamma_phi, -gamma_phi, 0.0])
        coherence_numeric.append((expm(generator * t) @ np.array([1.0, 0.0, 0.0]))[0])
    coherence_numeric = np.asarray(coherence_numeric)
    coherence_exact = np.exp(-gamma_phi * times)

    omega = 1.3
    rabi_numeric = []
    for t in times:
        unitary = expm(-1.0j * omega * t * SIGMA_X / 2.0)
        rho = unitary @ bloch_density(np.array([0.0, 0.0, 1.0])) @ unitary.conj().T
        rabi_numeric.append(float(np.real(np.trace(rho @ SIGMA_Z))))
    rabi_numeric = np.asarray(rabi_numeric)
    rabi_exact = np.cos(omega * times)

    gamma_grid = np.logspace(-2.0, 2.5, 120)
    zeno_survival = np.array(
        [
            zeno_excited_probability(gamma, drive_rate=1.0, duration=np.pi)
            for gamma in gamma_grid
        ]
    )

    angles = np.linspace(0.0, np.pi, 301)
    singlet_curve = -np.cos(angles)
    marginal_errors = []
    for _ in range(200):
        axis_a = rng.normal(size=3)
        axis_b = rng.normal(size=3)
        for outcome_a in (-1, 1):
            marginal = sum(
                singlet_joint_probability(outcome_a, b, axis_a, axis_b)
                for b in (-1, 1)
            )
            marginal_errors.append(abs(marginal - 0.5))

    metrics = {
        "born_max_abs_error": float(max(born_trace_errors)),
        "unitary_max_invariant_error": float(max(unitary_errors)),
        "dephasing_max_abs_error": float(
            np.max(np.abs(coherence_numeric - coherence_exact))
        ),
        "rabi_max_abs_error": float(np.max(np.abs(rabi_numeric - rabi_exact))),
        "chsh_abs": float(abs(chsh_value())),
        "chsh_tsirelson_error": float(abs(abs(chsh_value()) - 2.0 * np.sqrt(2.0))),
        "no_signalling_max_abs_error": float(max(marginal_errors)),
        "zeno_survival_low_gamma": float(zeno_survival[0]),
        "zeno_survival_high_gamma": float(zeno_survival[-1]),
    }
    arrays = {
        "born_formula": np.asarray(born_formula_values),
        "born_trace": np.asarray(born_trace_values),
        "times": times,
        "coherence": coherence_numeric,
        "coherence_exact": coherence_exact,
        "rabi": rabi_numeric,
        "rabi_exact": rabi_exact,
        "gamma_grid": gamma_grid,
        "zeno_survival": zeno_survival,
        "angles": angles,
        "singlet_curve": singlet_curve,
    }
    return metrics, arrays


def run_tiny_target() -> tuple[dict[str, object], dict[str, np.ndarray]]:
    rng = np.random.default_rng(SEED + 1)
    measurement_scale = 1.0
    epsilon_true = 0.02
    probe_time = 0.05
    z0 = 0.0

    calibration_angles = np.array([np.pi / 3.0, 2.0 * np.pi / 3.0, 4.0 * np.pi / 3.0])
    holdout_angles = np.array([np.pi / 2.0, np.pi, 3.0 * np.pi / 2.0])
    shots_each = 1_200_000
    rate_standard_error = np.sqrt(2.0 / (shots_each * probe_time))

    calibration_truth = target_rate_difference(
        calibration_angles, measurement_scale, epsilon_true, z0
    )
    calibration_observed = calibration_truth + rng.normal(
        scale=rate_standard_error, size=calibration_angles.size
    )
    epsilon_hat, epsilon_standard_error = fit_frozen_epsilon(
        calibration_angles,
        calibration_observed,
        measurement_scale=measurement_scale,
        z0=z0,
        weights=np.full(calibration_angles.size, 1.0 / rate_standard_error**2),
    )
    calibration_frozen = target_rate_difference(
        calibration_angles, measurement_scale, epsilon_hat, z0
    )

    holdout_frozen = target_rate_difference(
        holdout_angles, measurement_scale, epsilon_hat, z0
    )
    holdout_truth = target_rate_difference(
        holdout_angles, measurement_scale, epsilon_true, z0
    )
    holdout_observed = holdout_truth + rng.normal(
        scale=rate_standard_error, size=holdout_angles.size
    )
    residuals = holdout_observed - holdout_frozen
    holdout_design = target_rate_difference(
        holdout_angles, measurement_scale, 1.0, z0
    )
    holdout_covariance = (
        rate_standard_error**2 * np.eye(holdout_angles.size)
        + epsilon_standard_error**2
        * np.outer(holdout_design, holdout_design)
    )
    holdout_precision = np.linalg.inv(holdout_covariance)
    holdout_chi_square = float(residuals @ holdout_precision @ residuals)
    holdout_p_value = float(
        1.0 - chi2.cdf(holdout_chi_square, holdout_angles.size)
    )
    null_chi_square = float(
        np.sum((holdout_observed / rate_standard_error) ** 2)
    )
    _, epr_logdet = np.linalg.slogdet(holdout_covariance)
    null_logdet = holdout_angles.size * np.log(rate_standard_error**2)
    heldout_log_likelihood_gain = 0.5 * (
        null_chi_square
        + null_logdet
        - holdout_chi_square
        - epr_logdet
    )

    config = TargetConfig(
        measurement_scale=measurement_scale,
        epsilon=epsilon_true,
        probe_time=probe_time,
        time_step=0.001,
        z0=z0,
        solid_angle=np.pi,
    )
    y_plus, y_minus = simulate_target_records(
        config, paths_per_orientation=80_000, seed=SEED + 2
    )
    simulated_rate_difference = float(
        (np.mean(y_plus) - np.mean(y_minus)) / probe_time
    )
    simulated_rate_standard_error = float(
        np.sqrt(np.var(y_plus, ddof=1) / y_plus.size + np.var(y_minus, ddof=1) / y_minus.size)
        / probe_time
    )
    analytic_rate_difference = float(
        target_rate_difference(np.pi, measurement_scale, epsilon_true, z0)
    )

    epsilon_grid = np.logspace(-4.0, -1.0, 100)
    shots_grid = shots_per_orientation(
        epsilon_grid,
        measurement_scale=measurement_scale,
        probe_time=probe_time,
        solid_angle=np.pi,
        z0=z0,
    )
    mean_integrated = (
        measurement_scale
        * epsilon_grid
        * probe_time
        * np.sin(np.pi / 2.0)
    )
    information_bits = binary_gaussian_mutual_information(
        mean_integrated, variance=probe_time
    )
    fibre_variation_grid = np.array(
        [
            candidate_global_fibre_variation(
                measurement_scale, epsilon, probe_time
            )
            for epsilon in epsilon_grid
        ]
    )
    z_grid = np.linspace(-1.0, 1.0, 301)
    fibre_variation_profile = candidate_fibre_variation(
        z_grid,
        measurement_scale=measurement_scale,
        epsilon=epsilon_true,
        probe_time=probe_time,
    )
    gamma_grid = np.linspace(0.0, 2.0 * np.pi, 721, endpoint=False)
    sine_span = float(np.max(np.sin(gamma_grid)) - np.min(np.sin(gamma_grid)))
    brute_mean_separations = (
        measurement_scale
        * abs(epsilon_true)
        * (1.0 - z_grid**2)
        * probe_time
        * sine_span
    )
    brute_fibre_variation = equal_variance_gaussian_tv(
        brute_mean_separations, variance=probe_time
    )
    brute_global_variation = float(np.max(brute_fibre_variation))

    prediction_grid = np.linspace(0.0, 2.0 * np.pi, 401)
    frozen_curve = target_rate_difference(
        prediction_grid, measurement_scale, epsilon_hat, z0
    )

    metrics: dict[str, object] = {
        "synthetic_only": True,
        "measurement_scale": measurement_scale,
        "epsilon_true": epsilon_true,
        "epsilon_frozen": epsilon_hat,
        "epsilon_frozen_standard_error": epsilon_standard_error,
        "probe_time": probe_time,
        "shots_per_orientation_per_angle": shots_each,
        "rate_standard_error": rate_standard_error,
        "holdout_chi_square": holdout_chi_square,
        "holdout_degrees_of_freedom": int(holdout_angles.size),
        "holdout_goodness_p_value": holdout_p_value,
        "heldout_log_likelihood_gain_vs_memoryless_reduced_qubit_null": (
            heldout_log_likelihood_gain
        ),
        "analytic_rate_difference_at_pi": analytic_rate_difference,
        "sde_rate_difference_at_pi": simulated_rate_difference,
        "sde_rate_standard_error_at_pi": simulated_rate_standard_error,
        "sde_minus_short_probe_in_standard_errors": (
            simulated_rate_difference - analytic_rate_difference
        )
        / simulated_rate_standard_error,
        "candidate_global_fibre_variation": (
            candidate_global_fibre_variation(
                measurement_scale, epsilon_true, probe_time
            )
        ),
        "quantum_null_global_fibre_variation": (
            candidate_global_fibre_variation(
                measurement_scale, 0.0, probe_time
            )
        ),
        "brute_force_global_fibre_variation": brute_global_variation,
        "brute_force_vs_analytic_variation_error": abs(
            brute_global_variation
            - candidate_global_fibre_variation(
                measurement_scale, epsilon_true, probe_time
            )
        ),
    }
    arrays = {
        "calibration_angles": calibration_angles,
        "calibration_truth": calibration_truth,
        "calibration_observed": calibration_observed,
        "calibration_frozen": calibration_frozen,
        "holdout_angles": holdout_angles,
        "holdout_truth": holdout_truth,
        "holdout_frozen": holdout_frozen,
        "holdout_observed": holdout_observed,
        "prediction_grid": prediction_grid,
        "frozen_curve": frozen_curve,
        "epsilon_grid": epsilon_grid,
        "shots_grid": shots_grid,
        "information_bits": information_bits,
        "fibre_variation_grid": fibre_variation_grid,
        "z_grid": z_grid,
        "fibre_variation_profile": fibre_variation_profile,
        "y_plus": y_plus,
        "y_minus": y_minus,
    }
    return metrics, arrays


def figure_scalar_selection() -> None:
    degrees = np.arange(0, 9)
    dimensions = np.array([qnd_stabilizer_dimension(int(d)) for d in degrees])
    x = degrees

    fig, ax = plt.subplots(figsize=(6.7, 3.25))
    bar_colors = [COLORS["blue"] if d == 2 else "#aab3b9" for d in degrees]
    bars = ax.bar(x, dimensions, width=0.68, color=bar_colors, alpha=0.95)
    ax.axhline(
        2.0,
        color=COLORS["orange"],
        linewidth=2.0,
        label=r"minimal detector reference: $\dim W_P=2$",
    )
    family_labels = {
        1: (r"$\mathbb{R}$", 3.2),
        2: (r"$\mathbb{C}$", 6.2),
        4: (r"$\mathbb{H}$", 8.2),
        8: (r"$\mathbb{O}$", 30.2),
    }
    for degree, (label, y_position) in family_labels.items():
        ax.text(
            degree,
            y_position,
            label,
            ha="center",
            va="bottom",
            fontsize=9,
        )
    ax.set_xticks(x)
    ax.set_xlabel(r"Peirce degree $d$ (including arbitrary spin factors)")
    ax.set_ylabel(r"filter-stabilizer dimension $1+d(d-1)/2$")
    ax.set_title("Conditional match to unnormalized QND-filter generators")
    ax.set_ylim(0.0, 33.0)
    ax.legend(frameon=False, loc="upper left")
    ax.text(
        2,
        3.65,
        "unique dimensional match",
        ha="center",
        color=COLORS["orange"],
        fontsize=8,
    )
    fig.tight_layout()
    save_figure(fig, "scalar_selection")


def figure_relational_geometry() -> None:
    fig = plt.figure(figsize=(7.05, 3.35))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122, projection="3d")

    phi = 0.42
    direction = np.array([np.cos(phi), np.sin(phi)])
    transverse = np.array([-np.sin(phi), np.cos(phi)])
    x_values = np.linspace(-3.0, 3.0, 200)
    for c_value, color in zip((-1.2, 0.0, 1.2), ("#8c969d", COLORS["blue"], "#8c969d")):
        points = (
            x_values[:, None] * direction[None, :]
            + c_value * transverse[None, :]
        )
        ax1.plot(points[:, 0], points[:, 1], color=color, linewidth=1.8)
        ax1.text(points[-1, 0], points[-1, 1], rf"$c_\varphi={c_value:g}$", fontsize=7)
    ax1.arrow(
        0.0,
        0.0,
        direction[0],
        direction[1],
        width=0.025,
        color=COLORS["orange"],
        length_includes_head=True,
    )
    ax1.arrow(
        0.0,
        0.0,
        transverse[0],
        transverse[1],
        width=0.025,
        color=COLORS["green"],
        length_includes_head=True,
    )
    ax1.text(*(1.2 * direction), r"$x_\varphi$ (record-driven)", color=COLORS["orange"])
    ax1.text(*(1.25 * transverse), r"$c_\varphi$ (invariant)", color=COLORS["green"])
    ax1.set_xlabel(r"$s=\operatorname{artanh}z$")
    ax1.set_ylabel(r"relative phase $\vartheta$")
    ax1.set_title("QND trajectories are translations")
    ax1.set_aspect("equal", adjustable="box")
    ax1.set_xlim(-3.6, 3.6)
    ax1.set_ylim(-3.2, 3.2)

    u = np.linspace(0.0, 2.0 * np.pi, 80)
    v = np.linspace(0.0, np.pi, 40)
    sphere_x = np.outer(np.cos(u), np.sin(v))
    sphere_y = np.outer(np.sin(u), np.sin(v))
    sphere_z = np.outer(np.ones_like(u), np.cos(v))
    ax2.plot_wireframe(
        sphere_x,
        sphere_y,
        sphere_z,
        rstride=8,
        cstride=5,
        linewidth=0.35,
        color="#b9c0c5",
        alpha=0.55,
    )
    c_value = 0.25
    q = x_values[:, None] * direction[None, :] + c_value * transverse[None, :]
    s_values = q[:, 0]
    theta_values = q[:, 1]
    curve = np.column_stack(
        [
            1.0 / np.cosh(s_values) * np.cos(theta_values),
            1.0 / np.cosh(s_values) * np.sin(theta_values),
            np.tanh(s_values),
        ]
    )
    ax2.plot(
        curve[:, 0],
        curve[:, 1],
        curve[:, 2],
        color=COLORS["blue"],
        linewidth=2.4,
    )
    ax2.scatter(*curve[::35].T, color=COLORS["orange"], s=9)
    ax2.set_xlabel("$r_x$")
    ax2.set_ylabel("$r_y$")
    ax2.set_zlabel("$r_z$")
    ax2.set_title("The same path in the Bloch sphere")
    ax2.set_box_aspect((1.0, 1.0, 1.0))
    ax2.view_init(elev=23.0, azim=38.0)
    fig.tight_layout()
    save_figure(fig, "relational_geometry")


def figure_hopf_projection() -> None:
    fig, axes = plt.subplots(1, 2, figsize=(7.05, 3.1))
    ax = axes[0]
    chi = np.linspace(0.0, 2.0 * np.pi, 400)
    ax.plot(np.cos(chi), np.sin(chi), color=COLORS["blue"], linewidth=2.2)
    selected = np.array([0.15 * np.pi, 0.9 * np.pi, 1.55 * np.pi])
    ax.scatter(np.cos(selected), np.sin(selected), color=COLORS["orange"], zorder=3)
    for index, phase in enumerate(selected, 1):
        ax.text(
            1.13 * np.cos(phase),
            1.13 * np.sin(phase),
            rf"$\Lambda_{index}$",
            ha="center",
            va="center",
        )
    ax.annotate(
        "",
        xy=(1.95, 0.0),
        xytext=(1.25, 0.0),
        arrowprops={"arrowstyle": "->", "color": COLORS["ink"], "lw": 1.4},
    )
    ax.text(1.60, 0.18, r"$\mathcal{Q}$", ha="center", va="center")
    ax.text(2.30, 0.0, r"one $\rho$", va="center", ha="center", fontsize=11)
    ax.text(
        1.60,
        -0.34,
        r"$\rho=\psi\psi^\dagger$",
        ha="center",
        va="center",
        fontsize=9,
    )
    ax.set_aspect("equal")
    ax.set_xlim(-1.35, 2.75)
    ax.set_ylim(-1.35, 1.35)
    ax.axis("off")
    ax.set_title("A redundant $U(1)$ fibre projects to one state")

    ax = axes[1]
    omega = np.linspace(-2.0 * np.pi, 2.0 * np.pi, 401)
    shift = geometric_fibre_shift(omega)
    ax.plot(omega / np.pi, shift / np.pi, color=COLORS["green"], linewidth=2.2)
    ax.scatter(
        [-1.0, 1.0],
        [0.5, -0.5],
        color=[COLORS["orange"], COLORS["blue"]],
        zorder=3,
    )
    ax.axhline(0.0, color=COLORS["light"], linewidth=0.8)
    ax.axvline(0.0, color=COLORS["light"], linewidth=0.8)
    ax.set_xlabel(r"oriented solid angle $\Omega/\pi$")
    ax.set_ylabel(r"relational holonomy $\Gamma_{\mathrm{rel}}/\pi$")
    ax.set_title(r"Closed loops: $\Gamma_{\mathrm{rel}}=-\Omega/2$")
    fig.tight_layout()
    save_figure(fig, "hopf_projection")


def figure_qm_limits(arrays: dict[str, np.ndarray]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(7.05, 5.25))

    ax = axes[0, 0]
    ax.scatter(
        arrays["born_formula"][::5],
        arrays["born_trace"][::5],
        s=6,
        alpha=0.45,
        color=COLORS["blue"],
    )
    ax.plot([0.0, 1.0], [0.0, 1.0], color=COLORS["orange"], linewidth=1.2)
    ax.set_xlabel("Bloch-form probability")
    ax.set_ylabel(r"$\operatorname{Tr}(\rho E)$")
    ax.set_title("Born probabilities")
    ax.set_aspect("equal", adjustable="box")

    ax = axes[0, 1]
    ax.plot(arrays["times"], arrays["rabi"], label=r"$\langle Z\rangle$ under Rabi drive")
    ax.plot(
        arrays["times"],
        arrays["coherence"],
        label=r"$|\rho_{01}|/|\rho_{01}(0)|$ under dephasing",
    )
    ax.set_xlabel("time (inverse-rate units)")
    ax.set_ylabel("normalized observable")
    ax.set_title("Unitary rotation and dephasing")
    ax.legend(frameon=False)

    ax = axes[1, 0]
    ax.semilogx(
        arrays["gamma_grid"], arrays["zeno_survival"], color=COLORS["green"], linewidth=2.0
    )
    ax.set_xlabel(r"measurement dephasing $\Gamma/\Omega_R$")
    ax.set_ylabel("excited-state survival")
    ax.set_title("Quantum-Zeno suppression")
    ax.set_ylim(0.0, 1.02)

    ax = axes[1, 1]
    ax.plot(
        arrays["angles"] / np.pi,
        arrays["singlet_curve"],
        color=COLORS["red"],
        linewidth=2.0,
        label=r"$E(\mathbf{a},\mathbf{b})=-\cos\theta$",
    )
    ax.axhline(0.0, color=COLORS["light"], linewidth=0.8)
    ax.set_xlabel(r"analyzer separation $\theta/\pi$")
    ax.set_ylabel("correlation")
    ax.set_title(r"Singlet law; $|S_{\mathrm{CHSH}}|=2\sqrt{2}$")
    ax.legend(frameon=False)

    fig.tight_layout()
    save_figure(fig, "quantum_limits")


def figure_tiny_target(arrays: dict[str, np.ndarray], metrics: dict[str, object]) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.75))

    ax = axes[0]
    angle = np.linspace(0.0, 2.0 * np.pi, 300)
    radius = 0.48
    centers = (-0.62, 0.62)
    ax.plot(
        centers[0] + radius * np.cos(angle),
        radius * np.sin(angle),
        color=COLORS["blue"],
        linewidth=2.2,
        label=r"$C_+$",
    )
    ax.plot(
        centers[1] + radius * np.cos(angle),
        radius * np.sin(angle),
        color=COLORS["orange"],
        linewidth=2.2,
        label=r"$C_-$",
    )
    ax.scatter(
        [centers[0] + radius, centers[1] + radius],
        [0.0, 0.0],
        color=COLORS["ink"],
        s=10,
        zorder=3,
    )
    ax.annotate(
        "",
        xy=(
            centers[0] + radius * np.cos(0.9),
            radius * np.sin(0.9),
        ),
        xytext=(
            centers[0] + radius * np.cos(0.55),
            radius * np.sin(0.55),
        ),
        arrowprops={"arrowstyle": "->", "color": COLORS["blue"]},
    )
    ax.annotate(
        "",
        xy=(
            centers[1] + radius * np.cos(-0.9),
            radius * np.sin(-0.9),
        ),
        xytext=(
            centers[1] + radius * np.cos(-0.55),
            radius * np.sin(-0.55),
        ),
        arrowprops={"arrowstyle": "->", "color": COLORS["orange"]},
    )
    ax.set_aspect("equal")
    ax.set_xlim(-1.25, 1.45)
    ax.set_ylim(-0.72, 0.72)
    ax.axis("off")
    ax.set_title("Opposite closed histories\nsame projected endpoint")
    ax.legend(frameon=False, loc="lower center", ncol=2)

    ax = axes[1]
    ax.plot(
        arrays["prediction_grid"] / np.pi,
        arrays["frozen_curve"],
        color=COLORS["blue"],
        linewidth=2.0,
        label="frozen VD-Hopf-1 curve",
    )
    ax.errorbar(
        arrays["calibration_angles"] / np.pi,
        arrays["calibration_observed"],
        yerr=float(metrics["rate_standard_error"]),
        fmt="o",
        color=COLORS["gray"],
        ms=4,
        label="calibration",
    )
    ax.errorbar(
        arrays["holdout_angles"] / np.pi,
        arrays["holdout_observed"],
        yerr=float(metrics["rate_standard_error"]),
        fmt="s",
        color=COLORS["orange"],
        ms=4,
        label="synthetic holdout",
    )
    ax.axhline(
        0.0,
        color=COLORS["red"],
        linestyle="--",
        linewidth=1.0,
        label="memoryless reduced-qubit null",
    )
    ax.set_xlabel(r"loop solid angle $\Omega/\pi$")
    ax.set_ylabel(r"orientation-odd rate $\Delta_\star$")
    ax.set_title("One-parameter frozen prediction")
    ax.legend(frameon=False, fontsize=6.8)

    ax = axes[2]
    bins = np.linspace(
        min(np.quantile(arrays["y_plus"], 0.005), np.quantile(arrays["y_minus"], 0.005)),
        max(np.quantile(arrays["y_plus"], 0.995), np.quantile(arrays["y_minus"], 0.995)),
        75,
    )
    ax.hist(
        arrays["y_plus"],
        bins=bins,
        density=True,
        histtype="step",
        color=COLORS["blue"],
        linewidth=1.4,
        label=r"$C_+$",
    )
    ax.hist(
        arrays["y_minus"],
        bins=bins,
        density=True,
        histtype="step",
        color=COLORS["orange"],
        linewidth=1.4,
        label=r"$C_-$",
    )
    ax.set_xlabel(r"integrated record $Y_T$")
    ax.set_ylabel("density")
    ax.set_title("Synthetic short-probe records")
    ax.legend(frameon=False)
    fig.tight_layout()
    save_figure(fig, "tiny_target")


def figure_power(arrays: dict[str, np.ndarray]) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.65))

    ax = axes[0]
    ax.loglog(
        arrays["epsilon_grid"],
        2.0 * arrays["shots_grid"],
        color=COLORS["blue"],
        linewidth=2.2,
    )
    ax.set_xlabel(r"fibre coupling $|\varepsilon|$")
    ax.set_ylabel("total records, both orientations")
    ax.set_title(r"$5\sigma$, 90% power at $\Omega=\pi$")
    ax.grid(which="both", color="#e3e6e8", linewidth=0.45)

    ax = axes[1]
    ax.loglog(
        arrays["epsilon_grid"],
        arrays["information_bits"],
        color=COLORS["green"],
        linewidth=2.2,
    )
    ax.set_xlabel(r"fibre coupling $|\varepsilon|$")
    ax.set_ylabel(r"$I(H;Y_T\mid\rho)$ (bits/record)")
    ax.set_title("History-information witness")
    ax.grid(which="both", color="#e3e6e8", linewidth=0.45)

    ax = axes[2]
    ax.plot(
        arrays["z_grid"],
        arrays["fibre_variation_profile"],
        color=COLORS["orange"],
        linewidth=2.2,
    )
    ax.axhline(0.0, color=COLORS["red"], linestyle="--", linewidth=1.0)
    ax.set_xlabel(r"projected latitude $z$")
    ax.set_ylabel(r"$V_{\mathrm{p}}(\rho;\eta_0)$ (total variation)")
    ax.set_title(r"Matched protocol slice; max at $z=0$")
    ax.set_ylim(bottom=-0.001)
    fig.tight_layout()
    save_figure(fig, "power_and_information")


def write_csv_outputs(
    target_arrays: dict[str, np.ndarray], target_metrics: dict[str, object]
) -> None:
    with (OUTPUT_DIR / "frozen_predictions.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "split",
                "solid_angle_rad",
                "observed_rate_difference",
                "frozen_rate_prediction",
                "standard_error",
                "synthetic",
            ]
        )
        for angle, observed, frozen in zip(
            target_arrays["calibration_angles"],
            target_arrays["calibration_observed"],
            target_arrays["calibration_frozen"],
        ):
            writer.writerow(
                [
                    "calibration",
                    f"{angle:.12g}",
                    f"{observed:.12g}",
                    f"{frozen:.12g}",
                    f"{target_metrics['rate_standard_error']:.12g}",
                    "true",
                ]
            )
        for angle, observed, frozen in zip(
            target_arrays["holdout_angles"],
            target_arrays["holdout_observed"],
            target_arrays["holdout_frozen"],
        ):
            writer.writerow(
                [
                    "synthetic_holdout",
                    f"{angle:.12g}",
                    f"{observed:.12g}",
                    f"{frozen:.12g}",
                    f"{target_metrics['rate_standard_error']:.12g}",
                    "true",
                ]
            )

    selected_epsilons = np.array([0.05, 0.02, 0.01, 0.005, 0.001])
    selected_shots = shots_per_orientation(
        selected_epsilons, measurement_scale=1.0, probe_time=0.05
    )
    selected_information = binary_gaussian_mutual_information(
        selected_epsilons * 0.05, variance=0.05
    )
    with (OUTPUT_DIR / "power_table.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "epsilon",
                "records_per_orientation",
                "total_records",
                "mutual_information_bits_per_record",
            ]
        )
        for epsilon, shots, information in zip(
            selected_epsilons, selected_shots, selected_information
        ):
            writer.writerow(
                [
                    f"{epsilon:.8g}",
                    str(int(np.ceil(shots))),
                    str(int(2 * np.ceil(shots))),
                    f"{information:.12g}",
                ]
            )


def write_manifest(paths: list[Path]) -> None:
    rows = []
    for path in sorted(paths):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        rows.append(f"{digest}  {path.relative_to(ROOT)}")
    (OUTPUT_DIR / "manifest.sha256").write_text("\n".join(rows) + "\n")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    configure_plots()

    symbolic = symbolic_qnd_checks()
    qm_metrics, qm_arrays = verify_qm_limits()
    target_metrics, target_arrays = run_tiny_target()

    scalar_metrics = {
        "detector_reference_dimension": 2,
        "integer_solutions": compatible_peirce_degrees(),
        "qnd_stabilizer_dimensions": {
            str(d): qnd_stabilizer_dimension(d) for d in (1, 2, 4, 8)
        },
        "theorem_status": (
            "conditional on rank>=2 irreducible EJA geometry and an accessible "
            "injective-surjective response onto an unnormalized filter algebra"
        ),
    }
    hopf_metrics = {
        "global_phase_projection_error": float(
            np.max(
                [
                    np.linalg.norm(
                        density_from_spinor(spinor(0.37, -0.81, chi))
                        - density_from_spinor(spinor(0.37, -0.81, 0.0))
                    )
                    for chi in np.linspace(0.0, 2.0 * np.pi, 101)
                ]
            )
        ),
        "opposite_pi_loop_density_error": float(
            np.linalg.norm(
                density_from_spinor(
                    np.exp(1j * float(geometric_fibre_shift(np.pi)))
                    * spinor(0.0, 0.0)
                )
                - density_from_spinor(
                    np.exp(1j * float(geometric_fibre_shift(-np.pi)))
                    * spinor(0.0, 0.0)
                )
            )
        ),
        "court5_pair": (
            "algebraically inserted opposite spinor phases have identical "
            "projected rho; no control-loop or memory carrier is simulated"
        ),
    }

    metrics = {
        "evidence_status": {
            "analytic_and_computational": True,
            "laboratory_data": False,
            "physical_deviation_confirmed": False,
        },
        "random_seed": SEED,
        "symbolic_qnd": symbolic,
        "scalar_selection": scalar_metrics,
        "hopf_projection": hopf_metrics,
        "known_quantum_limits": qm_metrics,
        "ordinary_quantum_memory_countermodel": {
            "coherent_amplitude_magnitude": 0.3,
            "both_arms_photon_number": 0.09,
            "optimal_homodyne_total_variation": float(
                coherent_memory_homodyne_tv(0.3)
            ),
            "interpretation": (
                "same reduced qubit state and cavity photon number do not "
                "imply the same enlarged quantum state"
            ),
        },
        "tiny_target_synthetic_freeze_test": target_metrics,
    }
    metrics_path = OUTPUT_DIR / "metrics.json"
    metrics_path.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")

    freeze_specification = {
        "version": "VD-Hopf-1",
        "status": (
            "pure-state phenomenological record ansatz; no complete instrument "
            "and no laboratory confirmation"
        ),
        "memoryless_reduced_qubit_null": "epsilon = 0",
        "candidate_record_law": (
            "dY = a[z + epsilon(1-z^2)sin(Gamma_rel)]dt + dW"
        ),
        "closed_loop_holonomy": (
            "Gamma_rel_plus_minus = minus_plus Omega/2"
        ),
        "frozen_rate_difference": (
            "Delta(Omega) = -2 a epsilon (1-z0^2) "
            "cos(Gamma0) sin(Omega/2)"
        ),
        "measurement_scale": target_metrics["measurement_scale"],
        "probe_time": target_metrics["probe_time"],
        "z0": 0.0,
        "gamma_rel_offset": 0.0,
        "epsilon_frozen_from_calibration": target_metrics["epsilon_frozen"],
        "epsilon_standard_error": target_metrics[
            "epsilon_frozen_standard_error"
        ],
        "shots_per_orientation_per_angle": target_metrics[
            "shots_per_orientation_per_angle"
        ],
        "calibration_angles_rad": target_arrays["calibration_angles"].tolist(),
        "synthetic_holdout_angles_rad": target_arrays["holdout_angles"].tolist(),
        "implemented_score": (
            "cell-level Gaussian summary log likelihood with propagated "
            "calibration covariance; no target refit"
        ),
        "not_implemented": [
            "path likelihood",
            "process tensor",
            "qutrit-cavity artifact model",
            "external timestamped commitment",
        ],
        "synthetic_protocol_validation_only": True,
    }
    freeze_path = OUTPUT_DIR / "frozen_prediction.json"
    freeze_path.write_text(
        json.dumps(freeze_specification, indent=2, sort_keys=True) + "\n"
    )

    write_csv_outputs(target_arrays, target_metrics)
    figure_scalar_selection()
    figure_relational_geometry()
    figure_hopf_projection()
    figure_qm_limits(qm_arrays)
    figure_tiny_target(target_arrays, target_metrics)
    figure_power(target_arrays)

    manifest_inputs = [
        Path(__file__),
        ROOT / "evidence" / "epr_core.py",
        ROOT / "evidence" / "tests" / "test_epr_core.py",
        metrics_path,
        freeze_path,
        OUTPUT_DIR / "frozen_predictions.csv",
        OUTPUT_DIR / "power_table.csv",
        ROOT / "README.md",
        ROOT / "CLAIM_LEDGER.md",
        ROOT / "pyproject.toml",
        ROOT / "uv.lock",
        ROOT / "blind_qnd_research.py",
        ROOT / "experiment" / "PREREGISTRATION.md",
        ROOT / "thesis" / "main.tex",
        ROOT / "thesis" / "references.bib",
        ROOT / "Emergent-Predictive-Representation.pdf",
        *sorted(FIGURE_DIR.glob("*.pdf")),
    ]
    write_manifest(manifest_inputs)

    print(json.dumps(metrics, indent=2, sort_keys=True))
    print(f"Wrote evidence to {OUTPUT_DIR}")
    print(f"Wrote figures to {FIGURE_DIR}")


if __name__ == "__main__":
    main()
