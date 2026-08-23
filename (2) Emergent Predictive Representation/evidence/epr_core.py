"""Core mathematics for the Emergent Predictive Representation evidence suite.

The module keeps three logically different layers separate:

1. exact complex-qubit identities and the QND coordinate chart;
2. a conditional Jordan-algebra scalar-selection calculation;
3. a one-parameter phenomenological history-memory record ansatz.

The deformation is a candidate model, not an observed law of nature.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.linalg import expm
from scipy.special import expit, ndtr


ComplexArray = NDArray[np.complex128]
RealArray = NDArray[np.float64]

I2 = np.eye(2, dtype=complex)
SIGMA_X = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
SIGMA_Y = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex)
SIGMA_Z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
PAULIS = np.stack([SIGMA_X, SIGMA_Y, SIGMA_Z])


def spinor(s: float, theta: float, chi: float = 0.0) -> ComplexArray:
    """Return the normalized spinor in the relational chart.

    beta / alpha = exp(-s + i theta), while chi is the global Hopf-fibre
    coordinate forgotten by the density-matrix projection.
    """

    ratio = np.exp(-s + 1.0j * theta)
    psi = np.array([1.0, ratio], dtype=complex)
    psi /= np.linalg.norm(psi)
    return np.exp(1.0j * chi) * psi


def density_from_spinor(psi: ArrayLike) -> ComplexArray:
    """Hopf/ray projection Q(psi) = |psi><psi|."""

    vector = np.asarray(psi, dtype=complex)
    vector = vector / np.linalg.norm(vector)
    return np.outer(vector, vector.conj())


def relational_density(
    s: float, theta: float, purity_radius: float = 1.0
) -> ComplexArray:
    """Map real relational coordinates to a valid qubit density matrix."""

    if not 0.0 <= purity_radius <= 1.0:
        raise ValueError("purity_radius must lie in [0, 1]")
    direction = np.array(
        [
            1.0 / np.cosh(s) * np.cos(theta),
            1.0 / np.cosh(s) * np.sin(theta),
            np.tanh(s),
        ]
    )
    return bloch_density(purity_radius * direction)


def bloch_density(r: ArrayLike) -> ComplexArray:
    """Return rho = (I + r.sigma)/2 for a physical Bloch vector."""

    vector = np.asarray(r, dtype=float)
    if vector.shape != (3,):
        raise ValueError("Bloch vector must have shape (3,)")
    if np.linalg.norm(vector) > 1.0 + 1e-12:
        raise ValueError("Bloch vector lies outside the Bloch ball")
    return 0.5 * (I2 + np.tensordot(vector, PAULIS, axes=1))


def bloch_vector(rho: ArrayLike) -> RealArray:
    """Recover the three real Bloch coordinates from rho."""

    matrix = np.asarray(rho, dtype=complex)
    return np.real(np.array([np.trace(matrix @ p) for p in PAULIS]))


def hopf_map(psi: ArrayLike) -> RealArray:
    """Map a normalized spinor in S^3 to its point on S^2."""

    return bloch_vector(density_from_spinor(psi))


def qnd_coordinates(s: float, theta: float, phi: float) -> tuple[float, float]:
    """Longitudinal and transverse coordinates for homodyne angle phi."""

    x_phi = s * np.cos(phi) + theta * np.sin(phi)
    c_phi = -s * np.sin(phi) + theta * np.cos(phi)
    return float(x_phi), float(c_phi)


def qnd_coordinate_inverse(
    x_phi: float, c_phi: float, phi: float
) -> tuple[float, float]:
    """Invert the orthogonal QND chart."""

    s = x_phi * np.cos(phi) - c_phi * np.sin(phi)
    theta = x_phi * np.sin(phi) + c_phi * np.cos(phi)
    return float(s), float(theta)


def qnd_stabilizer_dimension(peirce_degree: int) -> int:
    """Dimension of R K plus so(d), modulo global cone dilation."""

    if peirce_degree < 0:
        raise ValueError("Peirce degree must be nonnegative")
    return 1 + peirce_degree * (peirce_degree - 1) // 2


def compatible_peirce_degrees(
    detector_reference_dimension: int = 2, search: Iterable[int] = range(0, 33)
) -> list[int]:
    """Solve dim(W_P) = 1 + d(d-1)/2 over candidate integer degrees."""

    return [
        d
        for d in search
        if qnd_stabilizer_dimension(d) == detector_reference_dimension
    ]


def born_probability(r: ArrayLike, axis: ArrayLike, outcome: int = 1) -> float:
    """Projective qubit probability in Bloch form."""

    vector = np.asarray(r, dtype=float)
    direction = np.asarray(axis, dtype=float)
    direction /= np.linalg.norm(direction)
    sign = 1.0 if outcome > 0 else -1.0
    return float(0.5 * (1.0 + sign * np.dot(vector, direction)))


def singlet_joint_probability(
    outcome_a: int, outcome_b: int, axis_a: ArrayLike, axis_b: ArrayLike
) -> float:
    """P(a,b|x,y) for the two-qubit singlet."""

    a = 1 if outcome_a > 0 else -1
    b = 1 if outcome_b > 0 else -1
    x = np.asarray(axis_a, dtype=float)
    y = np.asarray(axis_b, dtype=float)
    x /= np.linalg.norm(x)
    y /= np.linalg.norm(y)
    return float((1.0 - a * b * np.dot(x, y)) / 4.0)


def singlet_correlation(axis_a: ArrayLike, axis_b: ArrayLike) -> float:
    """E(a,b) = -a.b for the singlet."""

    x = np.asarray(axis_a, dtype=float)
    y = np.asarray(axis_b, dtype=float)
    return float(-np.dot(x / np.linalg.norm(x), y / np.linalg.norm(y)))


def chsh_value() -> float:
    """Tsirelson-saturating CHSH value for standard coplanar settings."""

    a0 = np.array([1.0, 0.0, 0.0])
    a1 = np.array([0.0, 1.0, 0.0])
    b0 = np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0)
    b1 = np.array([1.0, -1.0, 0.0]) / np.sqrt(2.0)
    return (
        singlet_correlation(a0, b0)
        + singlet_correlation(a0, b1)
        + singlet_correlation(a1, b0)
        - singlet_correlation(a1, b1)
    )


def zeno_excited_probability(
    measurement_dephasing: float,
    drive_rate: float,
    duration: float,
    initial_excited: bool = True,
) -> float:
    """Solve a driven qubit master equation with Z-dephasing.

    The Bloch generator is dx=-Gamma*x, dy=-Gamma*y-Omega*z,
    dz=Omega*y.  Increasing Gamma suppresses transitions.
    """

    gamma = float(measurement_dephasing)
    omega = float(drive_rate)
    generator = np.array(
        [[-gamma, 0.0, 0.0], [0.0, -gamma, -omega], [0.0, omega, 0.0]]
    )
    z0 = 1.0 if initial_excited else -1.0
    r0 = np.array([0.0, 0.0, z0])
    zt = float((expm(generator * duration) @ r0)[2])
    return 0.5 * (1.0 + zt)


def geometric_fibre_shift(solid_angle: ArrayLike) -> RealArray:
    """Geometric holonomy label under the convention Gamma_rel=-Omega/2.

    This function does not simulate a control loop or a physical memory carrier.
    """

    return -0.5 * np.asarray(solid_angle, dtype=float)


def deformation_term(z: ArrayLike, gamma_rel: ArrayLike) -> RealArray:
    """Selected finite-basis target term under the stated truncation.

    It is smooth, reference-covariant, odd under orientation reversal, and
    vanishes on the two sharp QND eigenstates.
    """

    z_array = np.asarray(z, dtype=float)
    gamma_array = np.asarray(gamma_rel, dtype=float)
    return (1.0 - z_array**2) * np.sin(gamma_array)


def epr_record_drift(
    z: ArrayLike, gamma_rel: ArrayLike, measurement_scale: float, epsilon: float
) -> RealArray:
    """Candidate record drift; epsilon=0 is the memoryless reduced-state null."""

    return measurement_scale * (
        np.asarray(z, dtype=float)
        + epsilon * deformation_term(z, gamma_rel)
    )


def target_rate_difference(
    solid_angle: ArrayLike,
    measurement_scale: float,
    epsilon: float,
    z0: float = 0.0,
    gamma_rel_offset: float = 0.0,
) -> RealArray:
    """Short-probe mean-rate difference, positive minus negative loop.

    gamma_rel_offset is the relational holonomy before the mirror loop pair.
    """

    omega = np.asarray(solid_angle, dtype=float)
    return (
        -2.0
        * measurement_scale
        * epsilon
        * (1.0 - z0**2)
        * np.cos(gamma_rel_offset)
        * np.sin(omega / 2.0)
    )


def equal_variance_gaussian_tv(
    mean_difference: ArrayLike, variance: float
) -> RealArray:
    """Total variation between two Gaussians with common variance.

    For N(mu_1, variance) and N(mu_2, variance),
    TV = 2 Phi(|mu_1-mu_2| / (2 sqrt(variance))) - 1.
    """

    if variance <= 0.0:
        raise ValueError("variance must be positive")
    separation = np.abs(np.asarray(mean_difference, dtype=float))
    return 2.0 * ndtr(separation / (2.0 * np.sqrt(variance))) - 1.0


def coherent_memory_homodyne_tv(
    alpha: complex | ArrayLike, quadrature_phase: float = 0.0
) -> RealArray:
    """TV distance for homodyne readout of coherent memories |+alpha>,|-alpha>.

    Both memories have the same photon number |alpha|^2.  For
    X_phi=(a e^{-i phi}+a^dagger e^{i phi})/sqrt(2), each distribution has
    variance 1/2 and opposite means.  This is an ordinary-quantum
    counterexample to inferring beyond-quantum physics from equal reduced
    qubit states or equal cavity photon number.
    """

    amplitudes = np.asarray(alpha, dtype=complex)
    projected_amplitude = np.abs(
        np.real(amplitudes * np.exp(-1.0j * quadrature_phase))
    )
    return 2.0 * ndtr(2.0 * projected_amplitude) - 1.0


def candidate_fibre_variation(
    z: ArrayLike,
    measurement_scale: float,
    epsilon: float,
    probe_time: float,
) -> RealArray:
    """Maximal short-probe TV distance on a matched pure-state fibre slice.

    The optimization is over two relational holonomy memories with
    sin(Gamma_1)=+1 and sin(Gamma_2)=-1 while all declared standard
    environment/process variables and probe resources are fixed.  It is a
    property of VD-Hopf-1, not a claim about Nature.
    """

    z_array = np.asarray(z, dtype=float)
    maximal_mean_difference = (
        2.0
        * measurement_scale
        * abs(epsilon)
        * (1.0 - z_array**2)
        * probe_time
    )
    return equal_variance_gaussian_tv(
        maximal_mean_difference, variance=probe_time
    )


def candidate_global_fibre_variation(
    measurement_scale: float, epsilon: float, probe_time: float
) -> float:
    """Protocol-fixed V_star for VD-Hopf-1 in the Gaussian limit."""

    return float(
        candidate_fibre_variation(
            z=0.0,
            measurement_scale=measurement_scale,
            epsilon=epsilon,
            probe_time=probe_time,
        )
    )


def shots_per_orientation(
    epsilon: ArrayLike,
    measurement_scale: float,
    probe_time: float,
    solid_angle: float = np.pi,
    z0: float = 0.0,
    gamma_rel_offset: float = 0.0,
    sigma_threshold: float = 5.0,
    power: float = 0.90,
) -> RealArray:
    """Gaussian short-probe shot budget per loop orientation.

    For equal sample counts N in the two arms, the difference-of-means
    standard error is sqrt(2/(N T)).  The normal approximation requires
    E[Z] >= sigma_threshold + Phi^{-1}(power).
    """

    from scipy.stats import norm

    eps = np.asarray(epsilon, dtype=float)
    delta_rate = np.abs(
        target_rate_difference(
            solid_angle,
            measurement_scale=measurement_scale,
            epsilon=eps,
            z0=z0,
            gamma_rel_offset=gamma_rel_offset,
        )
    )
    required_z = sigma_threshold + norm.ppf(power)
    with np.errstate(divide="ignore"):
        return 2.0 * required_z**2 / (delta_rate**2 * probe_time)


def binary_gaussian_mutual_information(
    mean: ArrayLike, variance: float, quadrature_order: int = 96
) -> RealArray:
    """I(H;Y) in bits for H=+/- equiprobable and Y|H~N(+/-mean,var).

    Gauss-Hermite quadrature evaluates
    I = 1 - E_Z log2(1 + exp(-2 m (m + sqrt(v) Z) / v)).
    """

    means = np.atleast_1d(np.asarray(mean, dtype=float))
    nodes, weights = np.polynomial.hermite.hermgauss(quadrature_order)
    standard_normal = np.sqrt(2.0) * nodes
    normalized_weights = weights / np.sqrt(np.pi)
    result = []
    for m in means:
        log_likelihood_ratio = 2.0 * m * (
            m + np.sqrt(variance) * standard_normal
        ) / variance
        conditional_entropy = np.sum(
            normalized_weights
            * np.logaddexp(0.0, -log_likelihood_ratio)
            / np.log(2.0)
        )
        result.append(max(0.0, 1.0 - float(conditional_entropy)))
    output = np.asarray(result)
    return output if np.ndim(mean) else output[0]


@dataclass(frozen=True)
class TargetConfig:
    measurement_scale: float = 1.0
    epsilon: float = 0.02
    probe_time: float = 0.05
    time_step: float = 0.001
    z0: float = 0.0
    solid_angle: float = np.pi


def simulate_target_records(
    config: TargetConfig,
    paths_per_orientation: int,
    seed: int = 20260821,
    return_final_z: bool = False,
) -> tuple[RealArray, ...]:
    """Euler-Maruyama simulation of the candidate QND target.

    The projected state follows the standard pure efficient QND martingale at
    leading order.  The extra drift appears only in the phenomenological record
    law; this is not a complete quantum instrument.  The relational holonomy
    memory Gamma_rel is fixed during the short readout.
    """

    if config.probe_time <= 0.0 or config.time_step <= 0.0:
        raise ValueError("probe_time and time_step must be positive")
    steps = int(round(config.probe_time / config.time_step))
    if not np.isclose(steps * config.time_step, config.probe_time):
        raise ValueError("probe_time must be an integer multiple of time_step")

    rng = np.random.default_rng(seed)
    outputs: list[RealArray] = []
    final_z: list[RealArray] = []
    for orientation in (1.0, -1.0):
        omega = orientation * config.solid_angle
        gamma_rel = float(geometric_fibre_shift(omega))
        s = np.full(paths_per_orientation, np.arctanh(config.z0))
        y = np.zeros(paths_per_orientation)
        for _ in range(steps):
            z = np.tanh(s)
            fibre_term = deformation_term(z, gamma_rel)
            drift = epr_record_drift(
                z, gamma_rel, config.measurement_scale, config.epsilon
            )
            dy = (
                drift * config.time_step
                + np.sqrt(config.time_step) * rng.standard_normal(y.size)
            )
            # Remove the candidate-only record drift from the state update so
            # that dz=a(1-z^2)dW and E[z_t]=z_0 remain QND at this order.
            s += config.measurement_scale * (
                dy
                - config.measurement_scale
                * config.epsilon
                * fibre_term
                * config.time_step
            )
            y += dy
        outputs.append(y)
        final_z.append(np.tanh(s))
    if return_final_z:
        return outputs[0], outputs[1], final_z[0], final_z[1]
    return outputs[0], outputs[1]


def fit_frozen_epsilon(
    solid_angles: ArrayLike,
    observed_rate_differences: ArrayLike,
    measurement_scale: float,
    z0: float = 0.0,
    weights: ArrayLike | None = None,
) -> tuple[float, float]:
    """Weighted one-parameter fit used only on the calibration split.

    Returns (epsilon_hat, standard_error).  Held-out predictions are obtained
    by substituting epsilon_hat into target_rate_difference without refitting.
    """

    omega = np.asarray(solid_angles, dtype=float)
    y = np.asarray(observed_rate_differences, dtype=float)
    design = target_rate_difference(
        omega, measurement_scale=measurement_scale, epsilon=1.0, z0=z0
    )
    supplied_weights = weights is not None
    w = np.ones_like(y) if weights is None else np.asarray(weights, dtype=float)
    information = float(np.sum(w * design**2))
    if information <= 0.0:
        raise ValueError("calibration design has zero information")
    epsilon_hat = float(np.sum(w * design * y) / information)
    if supplied_weights:
        standard_error = float(1.0 / np.sqrt(information))
    else:
        residual = y - epsilon_hat * design
        dof = max(1, y.size - 1)
        residual_variance = float(np.sum(w * residual**2) / dof)
        standard_error = float(np.sqrt(residual_variance / information))
    return epsilon_hat, standard_error


def bernoulli_log_likelihood(logit: ArrayLike, outcomes: ArrayLike) -> float:
    """Stable Bernoulli log likelihood, useful for model checks."""

    eta = np.asarray(logit, dtype=float)
    y = np.asarray(outcomes, dtype=float)
    p = expit(eta)
    return float(
        np.sum(
            y * np.log(np.clip(p, 1e-300, 1.0))
            + (1.0 - y) * np.log(np.clip(1.0 - p, 1e-300, 1.0))
        )
    )
