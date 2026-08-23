#!/usr/bin/env python3
"""
Blind QND coordinate / observability discovery experiment.

The discovery routines are NOT given:
  chi = theta - tan(phi) atanh(z)
  q   = atanh(z)
  the optimal Hamiltonian axis

They are given only the supplied Ito coefficients and an objective:
  1) find a first integral of the stochastic flow,
  2) find a quotient coordinate whose increment is proportional to dY,
  3) maximize a control observability criterion and construct a decoder.

Requires: numpy, scipy, sympy
"""

import math
import numpy as np
import sympy as sp
from scipy.linalg import expm

# ---------------------------------------------------------------------
# Stage 1: blind invariant discovery from stochastic characteristics
# ---------------------------------------------------------------------

z = sp.symbols("z", real=True)
a, phi = sp.symbols("a phi", positive=True, finite=True)

bz = a * sp.cos(phi) * (1 - z**2)
btheta = a * sp.sin(phi)
muz = sp.Integer(0)
mutheta = a**2 * sp.sin(phi) * sp.cos(phi) * z

# Search for a stochastic first integral F by solving the characteristic
# condition b . grad F = 0.  Normalize dF/dtheta = 1.
#
# Then F = theta + g(z) and g'(z) is forced by the SDE coefficients.
gprime = sp.simplify(-btheta / bz)
g = sp.integrate(gprime, z)

# Verify BOTH Ito terms vanish:
# drift(F) = F_z mu_z + F_theta mu_theta + 1/2 F_zz bz^2
gpp = sp.diff(gprime, z)
noise_F = sp.simplify(btheta + gprime * bz)
drift_F = sp.simplify(mutheta + gprime * muz + sp.Rational(1, 2) * gpp * bz**2)

# ---------------------------------------------------------------------
# Stage 2: blind observable quotient discovery
# ---------------------------------------------------------------------

# Measurement:
# dY = kappa*z dt + dW,  kappa = a cos(phi)
kappa = a * sp.cos(phi)

# Seek q=f(z) and a constant lambda so dq=lambda*dY.
# Fix the irrelevant scale by requiring f'(0)=1.
lam = sp.simplify(bz.subs(z, 0))
fprime = sp.simplify(lam / bz)
f = sp.integrate(fprime, z)

fpp = sp.diff(fprime, z)
diff_q = sp.simplify(fprime * bz)
drift_q = sp.simplify(fprime * muz + sp.Rational(1, 2) * fpp * bz**2)

diff_residual = sp.simplify(diff_q - lam)
drift_residual = sp.simplify(drift_q - lam * kappa * z)

# ---------------------------------------------------------------------
# Stage 3: blind Hamiltonian-axis discovery
# ---------------------------------------------------------------------

def observability_matrix(n, omega=1.0):
    nx, ny, nz = n
    return np.array([
        [0.0, 0.0, 1.0],
        [-omega * ny, omega * nx, 0.0],
        [omega**2 * nz * nx,
         omega**2 * nz * ny,
         omega**2 * (nz**2 - 1.0)]
    ])

rng = np.random.default_rng(20260820)
N_AXES = 200_000

axes = rng.normal(size=(N_AXES, 3))
axes /= np.linalg.norm(axes, axis=1)[:, None]

# D-optimal local observability score: |det O|.
scores = np.empty(N_AXES)
for i, n in enumerate(axes):
    scores[i] = abs(np.linalg.det(observability_matrix(n)))

best_idx = int(np.argmax(scores))
best_axis = axes[best_idx]
best_score = float(scores[best_idx])
best_angle_deg = float(np.degrees(np.arccos(abs(best_axis[2]))))

# Post-discovery analytic optimum for comparison only.
analytic_optimum_score = 2.0 / (3.0 * math.sqrt(3.0))
analytic_optimum_abs_nz = 1.0 / math.sqrt(3.0)

# ---------------------------------------------------------------------
# Construct R(Phi) using integrated mean measurement histories.
#
# Hamiltonian-only weak-probe forward map:
#    rdot = A r
#    E[dY]/dt = kappa z(t)
#
# Set kappa=1 for reconstruction scale.  For sample times t_j:
#    E[Y(t_j)] = D_j r0
#
# R is the Moore-Penrose inverse of D.
# ---------------------------------------------------------------------

n = best_axis
nx, ny, nz = n
A = np.array([
    [0.0, -nz, ny],
    [nz, 0.0, -nx],
    [-ny, nx, 0.0]
])

ez = np.array([0.0, 0.0, 1.0])

def integral_expm(A, t):
    dim = A.shape[0]
    block = np.block([
        [A, np.eye(dim)],
        [np.zeros((dim, dim)), np.zeros((dim, dim))]
    ])
    E = expm(block * t)
    return E[:dim, dim:]

sample_times = np.array([0.43, 1.10, 2.20, 3.00])
D = np.vstack([ez @ integral_expm(A, t) for t in sample_times])

rank_D = int(np.linalg.matrix_rank(D, tol=1e-12))
singular_values_D = np.linalg.svd(D, compute_uv=False)
condition_D = float(singular_values_D[0] / singular_values_D[-1])

R = np.linalg.pinv(D)

# Withheld validation states, uniformly distributed in the Bloch ball.
N_TEST = 10_000
states = rng.normal(size=(N_TEST, 3))
states /= np.linalg.norm(states, axis=1)[:, None]
states *= rng.random(N_TEST)[:, None] ** (1.0 / 3.0)

mean_records = states @ D.T
reconstructed = mean_records @ R.T
errors = np.linalg.norm(reconstructed - states, axis=1)

# Controls known to be degenerate are evaluated only after the blind search.
rank_z_axis = np.linalg.matrix_rank(observability_matrix(np.array([0.0, 0.0, 1.0])))
rank_x_axis = np.linalg.matrix_rank(observability_matrix(np.array([1.0, 0.0, 0.0])))

print("STAGE 1: invariant discovery")
print("  discovered g'(z) =", sp.simplify(gprime))
print("  discovered integral g(z) =", g)
print("  noise residual =", noise_F)
print("  Ito drift residual =", drift_F)
print()

print("STAGE 2: observable quotient discovery")
print("  discovered f'(z) =", sp.simplify(fprime))
print("  discovered integral f(z) =", f)
print("  discovered lambda =", lam)
print("  diffusion residual dq-lambda*dY =", diff_residual)
print("  drift residual dq-lambda*dY =", drift_residual)
print()

print("STAGE 3: blind control discovery")
print("  sampled axes =", N_AXES)
print("  best axis =", best_axis)
print("  |n_z| =", abs(best_axis[2]))
print("  angle from z (deg) =", best_angle_deg)
print("  score |det O| =", best_score)
print("  post-hoc analytic |n_z| optimum =", analytic_optimum_abs_nz)
print("  post-hoc analytic score optimum =", analytic_optimum_score)
print("  rank(z-axis control) =", rank_z_axis)
print("  rank(x-axis control) =", rank_x_axis)
print()

print("DECODER / HELD-OUT VERIFICATION")
print("  sample times =", sample_times)
print("  forward map rank =", rank_D)
print("  forward singular values =", singular_values_D)
print("  forward condition number =", condition_D)
print("  withheld states =", N_TEST)
print("  mean reconstruction error =", float(np.mean(errors)))
print("  max reconstruction error =", float(np.max(errors)))
