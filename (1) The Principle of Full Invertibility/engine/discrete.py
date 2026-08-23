"""
Exact finite analogue of the (r, v) skeleton of Ψ.

The continuous mechanical map is

    r' = r + v Δt
    v' = γ v + (Δt/m) F(r')

Its finite modular analog on Z/MZ is

    r' = r + v          (mod M)
    v' = g v + F(r')    (mod M)

If gcd(g, M) = 1 this map is a bijection of the finite set, regardless
of F: every state has exactly one predecessor.
"""

from __future__ import annotations

import math

import numpy as np


def mod_inv(a: int, m: int) -> int:
    return pow(int(a), -1, int(m))


def force_1d(r: int, m: int) -> int:
    # arbitrary nonlinear force — invertibility must not depend on this
    return (3 * r * r + 7 * r + 11) % m


def forward(r: int, v: int, m: int, g: int) -> tuple[int, int]:
    rp = (r + v) % m
    vp = (g * v + force_1d(rp, m)) % m
    return rp, vp


def inverse(rp: int, vp: int, m: int, g: int) -> tuple[int, int]:
    ginv = mod_inv(g, m)
    v = (ginv * (vp - force_1d(rp, m))) % m
    r = (rp - v) % m
    return r, v


def enumerate_bijection(m: int = 251, g: int = 3) -> dict:
    if math.gcd(g, m) != 1:
        raise ValueError("g must be coprime to M or the analog is not invertible")
    n = m * m
    image = np.full(n, -1, dtype=np.int32)
    collisions = 0
    for r in range(m):
        for v in range(m):
            rp, vp = forward(r, v, m, g)
            key = rp * m + vp
            if image[key] != -1:
                collisions += 1
            image[key] = r * m + v
    missing = int(np.sum(image < 0))
    # verify inverse on every state
    bad_inv = 0
    for r in range(m):
        for v in range(m):
            rp, vp = forward(r, v, m, g)
            rr, vv = inverse(rp, vp, m, g)
            if (rr, vv) != (r, v):
                bad_inv += 1
    return {
        "M": m,
        "g": g,
        "state_space": n,
        "collisions": collisions,
        "missing_images": missing,
        "bad_inverses": bad_inv,
        "is_bijection": collisions == 0 and missing == 0 and bad_inv == 0,
        "entropy_bits": float(np.log2(n)),
    }


def jacobian_block_det_numeric(
    dim_v: int, gamma: float, dt: float, DF: np.ndarray
) -> dict:
    """
    Build the (r, v) Jacobian of the Euler map in any dimension d = dim_v

        [ I              , dt I                      ]
        [ (dt) DF        , γ I + (dt²) DF            ]

    and confirm det = γ^d independently of DF.
    """
    d = dim_v
    I = np.eye(d)
    A = np.block(
        [
            [I, dt * I],
            [dt * DF, gamma * I + (dt**2) * DF],
        ]
    )
    det = float(np.linalg.det(A))
    return {
        "det": det,
        "gamma_pow": float(gamma**d),
        "rel_err": abs(det - gamma**d) / max(abs(gamma**d), 1e-15),
    }
