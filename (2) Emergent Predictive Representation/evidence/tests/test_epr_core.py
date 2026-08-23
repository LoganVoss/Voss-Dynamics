import numpy as np

from evidence.epr_core import (
    PAULIS,
    TargetConfig,
    bloch_vector,
    born_probability,
    candidate_fibre_variation,
    candidate_global_fibre_variation,
    chsh_value,
    coherent_memory_homodyne_tv,
    compatible_peirce_degrees,
    deformation_term,
    density_from_spinor,
    geometric_fibre_shift,
    hopf_map,
    qnd_coordinate_inverse,
    qnd_coordinates,
    relational_density,
    shots_per_orientation,
    simulate_target_records,
    singlet_joint_probability,
    spinor,
    target_rate_difference,
)


def test_relational_map_is_a_density_matrix() -> None:
    for s in (-4.0, -0.7, 0.0, 1.2, 4.0):
        for theta in np.linspace(0.0, 2.0 * np.pi, 9):
            for radius in (0.0, 0.4, 1.0):
                rho = relational_density(s, theta, radius)
                assert np.isclose(np.trace(rho), 1.0)
                assert np.allclose(rho, rho.conj().T)
                assert np.linalg.eigvalsh(rho).min() >= -1e-12


def test_hopf_projection_forgets_global_phase_only() -> None:
    psi = spinor(0.73, 1.17, 0.0)
    rho = density_from_spinor(psi)
    for chi in np.linspace(-3.0 * np.pi, 3.0 * np.pi, 17):
        shifted = np.exp(1.0j * chi) * psi
        assert np.allclose(density_from_spinor(shifted), rho, atol=1e-13)
        assert np.allclose(hopf_map(shifted), hopf_map(psi), atol=1e-13)


def test_qnd_chart_is_orthogonal_and_invertible() -> None:
    values = [
        (-2.1, -1.7, -0.8),
        (0.0, 0.0, 0.0),
        (0.4, 2.7, 1.3),
    ]
    for s, theta, phi in values:
        x_phi, c_phi = qnd_coordinates(s, theta, phi)
        s_back, theta_back = qnd_coordinate_inverse(x_phi, c_phi, phi)
        assert np.allclose([s_back, theta_back], [s, theta])
        assert np.isclose(x_phi**2 + c_phi**2, s**2 + theta**2)


def test_scalar_selection_integer_solution_is_complex_degree() -> None:
    assert compatible_peirce_degrees() == [2]


def test_born_probability_matches_trace_rule() -> None:
    rng = np.random.default_rng(11)
    for _ in range(100):
        r = rng.normal(size=3)
        r /= max(1.0, np.linalg.norm(r))
        axis = rng.normal(size=3)
        axis /= np.linalg.norm(axis)
        rho = 0.5 * (
            np.eye(2, dtype=complex) + np.tensordot(r, PAULIS, axes=1)
        )
        effect = 0.5 * (
            np.eye(2, dtype=complex) + np.tensordot(axis, PAULIS, axes=1)
        )
        assert np.isclose(
            born_probability(r, axis), np.real(np.trace(rho @ effect))
        )


def test_singlet_is_no_signalling_and_saturates_tsirelson() -> None:
    rng = np.random.default_rng(17)
    for _ in range(20):
        axis_a = rng.normal(size=3)
        axis_b = rng.normal(size=3)
        for outcome_a in (-1, 1):
            marginal = sum(
                singlet_joint_probability(outcome_a, outcome_b, axis_a, axis_b)
                for outcome_b in (-1, 1)
            )
            assert np.isclose(marginal, 0.5)
    assert np.isclose(abs(chsh_value()), 2.0 * np.sqrt(2.0))


def test_deformation_obeys_target_symmetries() -> None:
    z = np.linspace(-1.0, 1.0, 21)
    gamma_rel = 0.73
    assert np.allclose(
        deformation_term(z, -gamma_rel), -deformation_term(z, gamma_rel)
    )
    assert np.allclose(
        deformation_term(-z, gamma_rel), deformation_term(z, gamma_rel)
    )
    assert np.isclose(deformation_term(1.0, gamma_rel), 0.0)
    assert np.isclose(deformation_term(-1.0, gamma_rel), 0.0)


def test_geometric_loop_pair_has_same_rho_but_opposite_fibre() -> None:
    psi = spinor(0.0, 0.0)
    rho = density_from_spinor(psi)
    omega = np.pi
    chi_plus = float(geometric_fibre_shift(omega))
    chi_minus = float(geometric_fibre_shift(-omega))
    assert np.allclose(density_from_spinor(np.exp(1j * chi_plus) * psi), rho)
    assert np.allclose(density_from_spinor(np.exp(1j * chi_minus) * psi), rho)
    assert np.isclose(chi_plus, -chi_minus)


def test_target_signal_and_shot_budget_scaling() -> None:
    delta = target_rate_difference(
        np.pi, measurement_scale=1.0, epsilon=0.02, z0=0.0
    )
    assert np.isclose(delta, -0.04)
    assert np.isclose(
        target_rate_difference(
            np.pi,
            measurement_scale=1.0,
            epsilon=0.02,
            z0=0.0,
            gamma_rel_offset=np.pi / 2.0,
        ),
        0.0,
        atol=1e-15,
    )
    n1 = shots_per_orientation(0.02, 1.0, 0.05)
    n2 = shots_per_orientation(0.01, 1.0, 0.05)
    assert np.isclose(n2 / n1, 4.0)


def test_fibre_variation_is_zero_iff_candidate_coupling_is_zero() -> None:
    assert candidate_global_fibre_variation(1.0, 0.0, 0.05) == 0.0
    positive = candidate_global_fibre_variation(1.0, 0.02, 0.05)
    assert positive > 0.0
    z = np.linspace(-1.0, 1.0, 101)
    profile = candidate_fibre_variation(z, 1.0, 0.02, 0.05)
    assert np.argmax(profile) == 50
    assert np.isclose(profile[0], 0.0)
    assert np.isclose(profile[-1], 0.0)


def test_standard_quantum_coherent_memory_can_mimic_history_signal() -> None:
    alpha = 0.3
    assert np.isclose(abs(alpha) ** 2, abs(-alpha) ** 2)
    assert coherent_memory_homodyne_tv(0.0) == 0.0
    assert coherent_memory_homodyne_tv(alpha) > 0.0
    assert np.isclose(
        coherent_memory_homodyne_tv(1j * alpha, quadrature_phase=0.0),
        0.0,
    )
    assert coherent_memory_homodyne_tv(
        1j * alpha, quadrature_phase=np.pi / 2.0
    ) > 0.0


def test_candidate_state_update_preserves_qnd_population_martingale() -> None:
    config = TargetConfig(epsilon=0.2, probe_time=0.05, time_step=0.001, z0=0.0)
    _, _, z_plus, z_minus = simulate_target_records(
        config,
        paths_per_orientation=20_000,
        seed=1234,
        return_final_z=True,
    )
    assert abs(float(np.mean(z_plus)) - config.z0) < 0.01
    assert abs(float(np.mean(z_minus)) - config.z0) < 0.01


def test_bloch_round_trip() -> None:
    rho = relational_density(0.44, -0.91, 0.72)
    r = bloch_vector(rho)
    assert np.linalg.norm(r) <= 1.0 + 1e-12
