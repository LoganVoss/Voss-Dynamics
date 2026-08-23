#!/usr/bin/env python3
"""Equal-ω phase memory, Lyapunov of Ψ, and γ=0 collision census."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from engine.discrete import enumerate_bijection
from engine.psi import FieldParams, euler_step, kuramoto_R, make_field, wrap


def _jsonable(x):
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, (np.floating, float)):
        v = float(x)
        return None if math.isnan(v) or math.isinf(v) else v
    if isinstance(x, (np.integer, int)):
        return int(x)
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    return x


def dphi(a, b):
    return float(abs(np.arctan2(np.sin(a - b), np.cos(a - b))))


def phase_memory(equal_omega: bool) -> dict:
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
    f = make_field(n=2, seed=2, params=prm)
    f.r = np.array([[-0.12, 0.0], [0.12, 0.0]])
    f.v[:] = 0
    f.phi = np.array([0.4, 2.9])
    if equal_omega:
        f.omega[:] = 0.05
    else:
        f.omega = np.array([0.03, 0.08])
    for _ in range(100):
        euler_step(f)
        f.r = np.array([[-0.12, 0.0], [0.12, 0.0]])
        f.v[:] = 0
    locked = dphi(f.phi[0], f.phi[1])
    f.r = np.array([[-15.0, 0.0], [15.0, 0.0]])
    f.v[:] = 0
    trace = []
    for t in range(300):
        euler_step(f)
        f.r = np.array([[-15.0, 0.0], [15.0, 0.0]])
        f.v[:] = 0
        if t in (0, 24, 49, 99, 199, 299):
            trace.append({"t": t + 1, "dphi": dphi(f.phi[0], f.phi[1]), "R": kuramoto_R(f)})
    return {
        "equal_omega": equal_omega,
        "omega": f.omega.tolist(),
        "dphi_locked": locked,
        "trace": trace,
        "dphi_final": trace[-1]["dphi"],
        "drift": trace[-1]["dphi"] - trace[0]["dphi"],
    }


def lyapunov() -> dict:
    prm = FieldParams(gamma=1.0, theta_drift=0.0, wall_bounce=0.0)
    a = make_field(n=8, seed=5, params=prm)
    b = a.copy()
    eps0 = 1e-8
    b.r = b.r + eps0
    logs = []
    d0 = float(np.linalg.norm(a.flatten() - b.flatten()))
    for t in range(1, 181):
        euler_step(a)
        euler_step(b)
        d = float(np.linalg.norm(a.flatten() - b.flatten()))
        if t in (20, 40, 80, 120, 160, 180):
            logs.append(
                {
                    "t": t,
                    "sep": d,
                    "lambda_proxy": math.log(d / d0) / t,
                }
            )
    return {"d0": d0, "trace": logs, "lambda_180": logs[-1]["lambda_proxy"]}


def gamma_zero_census() -> dict:
    # finite analog with g=0: every v collapses
    m = 64
    seen = set()
    for r in range(m):
        for v in range(m):
            rp = (r + v) % m
            vp = (0 * v + (3 * rp * rp + 7 * rp + 11)) % m
            seen.add((rp, vp))
    coprime = enumerate_bijection(m=127, g=5)
    return {
        "M": m,
        "g": 0,
        "state_space": m * m,
        "unique_images": len(seen),
        "lost_states": m * m - len(seen),
        "coprime_127_5_is_bijection": coprime["is_bijection"],
        "coprime_state_space": coprime["state_space"],
    }


def main():
    results = {
        "equal_omega": phase_memory(True),
        "unequal_omega": phase_memory(False),
        "lyapunov": lyapunov(),
        "gamma_zero_census": gamma_zero_census(),
    }
    out = ROOT / "results" / "omega.json"
    out.write_text(json.dumps(_jsonable(results), indent=2))
    print(json.dumps(_jsonable(results), indent=2))
    print("wrote", out)


if __name__ == "__main__":
    main()
