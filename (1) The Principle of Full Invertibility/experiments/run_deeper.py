#!/usr/bin/env python3
"""Second-wave tests. Follow the fractures the first battery opened."""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from engine.observe import observe, promote, reset_tokens
from engine.psi import (
    FieldParams,
    euler_inverse_step,
    euler_step,
    kuramoto_R,
    make_field,
    run,
    run_inverse,
    state_distance,
    wrap,
)


def _jsonable(x):
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, (np.floating, float)):
        v = float(x)
        if math.isnan(v) or math.isinf(v):
            return None
        return v
    if isinstance(x, (np.integer, int)):
        return int(x)
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if x is None:
        return None
    return x


def experiment_chunked_inverse() -> dict:
    """
    Compare one-shot inversion with independently verified checkpoint windows.

    Stored checkpoints do not improve a single final-state reconstruction.
    They allow each short interval to be inverted from its own exact endpoint,
    which separates local window conditioning from accumulated long-horizon
    error.
    """
    rows = []
    for steps, chunk, gamma in [
        (200, 20, 0.982),
        (200, 10, 0.982),
        (400, 20, 0.982),
        (250, 10, 1.0),
        (400, 20, 1.0),
        (640, 16, 1.0),
    ]:
        prm = FieldParams(gamma=gamma, theta_drift=0.0, wall_bounce=0.0)
        f = make_field(n=8, seed=21, params=prm)
        origin = f.copy()
        checkpoints = [origin.copy()]
        for _ in range(0, steps, chunk):
            run(f, chunk, kind="euler")
            checkpoints.append(f.copy())
        # one-shot from the end
        oneshot = checkpoints[-1].copy()
        run_inverse(oneshot, steps, kind="euler")
        one = state_distance(origin, oneshot)
        # Invert every checkpoint window independently. Do not pretend that
        # resetting to a stored endpoint is one final-state decode.
        window_err = []
        for k in range(len(checkpoints) - 1, 0, -1):
            reconstructed = checkpoints[k].copy()
            run_inverse(reconstructed, chunk, kind="euler")
            window_err.append(
                state_distance(checkpoints[k - 1], reconstructed)["l2"]
            )
        rows.append(
            {
                "steps": steps,
                "chunk": chunk,
                "gamma": gamma,
                "oneshot_l2": one["l2"],
                "window_max_l2": max(window_err) if window_err else None,
                "window_median_l2": float(np.median(window_err))
                if window_err
                else None,
                "window_success_fraction_1e4": float(
                    np.mean(np.asarray(window_err) < 1e-4)
                )
                if window_err
                else None,
                "oneshot_ok": one["l2"] < 1e-4,
            }
        )
    return {"rows": rows}


def experiment_true_single_state_decode() -> dict:
    """
    Decode using ONLY the final microstate — no checkpoints.
    Sweep T until float64 loses the past. This is the epistemic horizon
    of a finite-precision observer who nevertheless holds the full z.
    """
    rows = []
    for gamma in (1.0, 0.982):
        for steps in (20, 40, 60, 80, 100, 120, 160, 200):
            prm = FieldParams(gamma=gamma, theta_drift=0.0, wall_bounce=0.0)
            f = make_field(n=8, seed=21, params=prm)
            origin = f.copy()
            run(f, steps, kind="euler")
            run_inverse(f, steps, kind="euler")
            d = state_distance(origin, f)
            rows.append(
                {
                    "gamma": gamma,
                    "steps": steps,
                    "l2": d["l2"],
                    "max_dr": d["max_dr"],
                    "ok_1e6": d["l2"] < 1e-6,
                }
            )
    return {"rows": rows}


def experiment_walled_attractor() -> dict:
    """Original visualization used walls. Does meaning-compression appear then?"""
    rows = []
    for wall_bounce, gamma, label in [
        (0.0, 1.0, "free_reversible"),
        (0.0, 0.982, "free_damped"),
        (0.4, 0.982, "walled_damped"),
        (0.4, 1.0, "walled_reversible"),
    ]:
        prm = FieldParams(
            gamma=gamma, theta_drift=0.0, wall_bounce=wall_bounce, wall=3.6
        )
        Ls, vars_, Rs, locked = [], [], [], []
        for s in range(30):
            f = make_field(n=14, seed=70 + s, params=prm)
            start = observe(f)
            run(f, 100, kind="euler")
            end = observe(f)
            Ls.append(end["L_total"] - start["L_total"])
            vars_.append(end["residual_var"])
            Rs.append(kuramoto_R(f))
            # fraction of tokens inside the lock radius of a declared anchor
            d = np.linalg.norm(f.r[:, None, :] - f.p[:6][None, :, :], axis=-1)
            locked.append(float(np.mean(d.min(axis=1) < 1.15)))
        rows.append(
            {
                "label": label,
                "mean_dL": float(np.mean(Ls)),
                "frac_end_shorter": float(np.mean(np.array(Ls) < 0)),
                "mean_residual_var": float(np.mean(vars_)),
                "mean_R": float(np.mean(Rs)),
                "mean_lock_frac": float(np.mean(locked)),
            }
        )
    return {"rows": rows}


def experiment_frozen_entanglement() -> dict:
    """Lock, then freeze positions far apart so they cannot drift back."""
    prm = FieldParams(
        gamma=1.0,
        beta=0.4,
        dt=0.25,
        alpha=0.0,
        wave_amp=0.0,
        mutual_rep=0.0,
        mutual_att=0.0,
        theta_drift=0.0,
        wall_bounce=0.0,
    )
    f = make_field(n=2, seed=1, params=prm, spread=0.2)
    f.r = np.array([[-0.15, 0.0], [0.15, 0.0]])
    f.v[:] = 0
    f.phi = np.array([0.3, 2.7])
    for _ in range(80):
        euler_step(f)
        f.v[:] = 0
        f.r = np.array([[-0.15, 0.0], [0.15, 0.0]])
    dphi_lock = float(
        abs(np.arctan2(np.sin(f.phi[0] - f.phi[1]), np.cos(f.phi[0] - f.phi[1])))
    )
    R_lock = kuramoto_R(f)
    # teleport far and freeze positions; only phases may run
    f.r = np.array([[-12.0, 0.0], [12.0, 0.0]])
    f.v[:] = 0
    trace = []
    for t in range(200):
        euler_step(f)
        f.r = np.array([[-12.0, 0.0], [12.0, 0.0]])
        f.v[:] = 0
        if t in (0, 19, 49, 99, 199):
            trace.append(
                {
                    "t": t + 1,
                    "dphi": float(
                        abs(
                            np.arctan2(
                                np.sin(f.phi[0] - f.phi[1]),
                                np.cos(f.phi[0] - f.phi[1]),
                            )
                        )
                    ),
                    "R": kuramoto_R(f),
                }
            )
    pre1 = float(f.phi[1])
    f.phi[0] = 0.0
    instant = float(abs(np.sin(f.phi[1] - pre1)))
    euler_step(f)
    f.r = np.array([[-12.0, 0.0], [12.0, 0.0]])
    after = float(
        abs(np.arctan2(np.sin(f.phi[1] - pre1), np.cos(f.phi[1] - pre1)))
    )
    # control: same snap at small distance
    g = make_field(n=2, seed=1, params=prm, spread=0.2)
    g.r = np.array([[-0.15, 0.0], [0.15, 0.0]])
    g.v[:] = 0
    g.phi = np.array([0.3, 2.7])
    for _ in range(80):
        euler_step(g)
        g.v[:] = 0
        g.r = np.array([[-0.15, 0.0], [0.15, 0.0]])
    pre = float(g.phi[1])
    g.phi[0] = 0.0
    euler_step(g)
    near_after = float(abs(np.arctan2(np.sin(g.phi[1] - pre), np.cos(g.phi[1] - pre))))
    return {
        "dphi_after_lock": dphi_lock,
        "R_after_lock": R_lock,
        "separated_trace": trace,
        "dphi_drift_over_200": trace[-1]["dphi"] - trace[0]["dphi"],
        "instantaneous_other_phase_change": instant,
        "far_step_influence": after,
        "near_step_influence": near_after,
        "influence_ratio_near_over_far": near_after / max(after, 1e-15),
    }


def experiment_gamma_zero() -> dict:
    """Construct a translational collision at γ=0."""
    prm = FieldParams(gamma=0.0, theta_drift=0.0, wall_bounce=0.0)
    f = make_field(n=6, seed=3, params=prm)
    a = f.copy()
    b = f.copy()
    b.v += 0.4
    euler_step(a)
    euler_step(b)
    # after one step, velocity is determined only by F(r'), and r' = r + v dt
    # so different v still changes r' — γ=0 does not kill injectivity of the
    # full step unless dt=0. Kill dt-coupling: compare two states with same
    # r, different v, after imagining the map at dt=0 — skip.
    # Better: two states that reach the same r' with γ=0:
    # r' = r + v dt. Choose (r1,v1), (r2,v2) with r1+v1 dt = r2+v2 dt
    # and same phases. Then F same, v' = 0 + (dt/m)F same. Collision.
    g = make_field(n=6, seed=3, params=prm)
    h = g.copy()
    h.v = g.v + 0.7
    h.r = g.r - 0.7 * g.params.dt
    # now r+v dt equal
    r1 = g.r + g.v * g.params.dt
    r2 = h.r + h.v * h.params.dt
    euler_step(g)
    euler_step(h)
    return {
        "constructed_same_r_prime": float(np.max(np.abs(r1 - r2))),
        "final_r_l2": float(np.linalg.norm(g.r - h.r)),
        "final_v_l2": float(np.linalg.norm(g.v - h.v)),
        "final_phi_l2": float(np.linalg.norm(np.sin(g.phi - h.phi))),
        "collision": float(np.linalg.norm(g.flatten() - h.flatten())) < 1e-10,
        "note": "γ=0 sends every velocity to (Δt/m)F(r'). Distinct (r,v) with the same r+vΔt become the same future.",
    }


def experiment_recursion() -> dict:
    prm = FieldParams(gamma=0.982, theta_drift=0.0, wall_bounce=0.4, wall=3.6)
    f = make_field(n=16, seed=17, params=prm)
    layers = []
    for ell in range(25):
        run(f, 50, kind="euler")
        obs = observe(f)
        layers.append(
            {
                "layer": ell + 1,
                "phrase": obs["phrase"],
                "L_total": obs["L_total"],
                "L_data": obs["L_data"],
                "L_phrase": obs["L_phrase"],
                "n_anchors": int(len(f.Q)),
                "residual_var": obs["residual_var"],
                "R": kuramoto_R(f),
                "n_clusters": obs["n_clusters"],
            }
        )
        promote(f, obs["clusters"], charge_scale=1.12)
        reset_tokens(f, seed=1000 + ell, spread=3.0)
    phrases = [L["phrase"] for L in layers]
    return {
        "layers": layers,
        "unique_phrases": len(set(phrases)),
        "L_data_first": layers[0]["L_data"],
        "L_data_last": layers[-1]["L_data"],
        "residual_first": layers[0]["residual_var"],
        "residual_last": layers[-1]["residual_var"],
        "anchors_last": layers[-1]["n_anchors"],
        "self_similar": len(set(phrases)) < 10,
    }


def experiment_quantized_shadow() -> dict:
    """
    Discrete shadow of the continuous law: quantize full Z at several
    precisions. Count collisions among 400 ICs as a function of time.
    If the law is injective, collisions appear only from rounding.
    """
    prm = FieldParams(gamma=0.982, theta_drift=0.0, wall_bounce=0.0)
    ens = 200
    fields = [make_field(n=6, seed=5000 + i, params=prm) for i in range(ens)]
    rows = []
    for t in (0, 20, 60, 120):
        if t > 0:
            # advance from last
            prev = 0
            # we rebuild from start for cleanliness
            fields = [make_field(n=6, seed=5000 + i, params=prm) for i in range(ens)]
            for f in fields:
                run(f, t, kind="euler")
        flats = np.stack([f.flatten() for f in fields])
        rec = {"t": t, "ens": ens}
        for dec in (6, 4, 3, 2):
            u = len(np.unique(np.round(flats, dec), axis=0))
            rec[f"unique_{dec}dec"] = int(u)
            rec[f"collisions_{dec}dec"] = int(ens - u)
        # positions only
        pos = np.stack([f.r.ravel() for f in fields])
        rec["unique_pos_3dec"] = int(len(np.unique(np.round(pos, 3), axis=0)))
        rows.append(rec)
    return {"rows": rows}


def experiment_three_body_phase() -> dict:
    """Frustration: three tokens, equally spaced, can they all lock?"""
    prm = FieldParams(
        gamma=1.0,
        beta=0.35,
        alpha=0.0,
        wave_amp=0.0,
        mutual_rep=0.0,
        mutual_att=0.0,
        theta_drift=0.0,
        wall_bounce=0.0,
        dt=0.2,
    )
    f = make_field(n=3, seed=0, params=prm)
    f.r = np.array([[1.0, 0.0], [-0.5, 0.866], [-0.5, -0.866]])
    f.v[:] = 0
    f.phi = np.array([0.0, 2.0, 4.0])
    hist = []
    for t in range(200):
        euler_step(f)
        f.r = np.array([[1.0, 0.0], [-0.5, 0.866], [-0.5, -0.866]])
        f.v[:] = 0
        if t in (19, 49, 99, 199):
            hist.append(
                {
                    "t": t + 1,
                    "R": kuramoto_R(f),
                    "phi": wrap(f.phi).tolist(),
                }
            )
    return {
        "trace": hist,
        "final_R": hist[-1]["R"],
        "locks_to_global_sync": hist[-1]["R"] > 0.9,
    }


def experiment_beta_threshold() -> dict:
    rows = []
    for beta in np.linspace(0.05, 0.8, 16):
        prm = FieldParams(gamma=0.982, beta=float(beta), theta_drift=0.0, wall_bounce=0.0)
        f = make_field(n=10, seed=6, params=prm)
        origin = f.copy()
        run(f, 30, kind="euler")
        res = run_inverse(f, 30, kind="euler")
        d = state_distance(origin, f)
        rows.append(
            {
                "beta": float(beta),
                "l2": d["l2"],
                "phase_residual": max(res),
                "ok": d["l2"] < 1e-6,
            }
        )
    ok_betas = [r["beta"] for r in rows if r["ok"]]
    return {
        "rows": rows,
        "max_safe_beta": max(ok_betas) if ok_betas else None,
        "min_fail_beta": min((r["beta"] for r in rows if not r["ok"]), default=None),
    }


def experiment_phase_only_information() -> dict:
    """Holding r,v,θ fixed, test the phase map at the default β."""
    prm = FieldParams(beta=0.08, dt=0.35, theta_drift=0.0, wall_bounce=0.0)
    f = make_field(n=12, seed=9, params=prm)
    # freeze r,v — only step phases via the same update used in Ψ
    from engine.psi import _advance_phases

    origin_phi = f.phi.copy()
    r0, v0 = f.r.copy(), f.v.copy()
    for _ in range(40):
        _advance_phases(f)
        f.r, f.v = r0, v0
    # invert 40 phase-only steps
    from engine.psi import pairwise

    def invert_phase_once(field):
        _, dist = pairwise(field.r, field.params.eps)
        phi_old = wrap(field.phi - field.omega * field.params.dt)
        for _ in range(40):
            s = np.sin(phi_old[None, :] - phi_old[:, None])
            couple = np.sum(s / dist, axis=1)
            phi_old = wrap(
                field.phi
                - field.omega * field.params.dt
                - field.params.beta * field.params.dt * couple
            )
        field.phi = phi_old

    for _ in range(40):
        invert_phase_once(f)
        f.r, f.v = r0, v0
    err = float(
        np.max(np.abs(np.arctan2(np.sin(f.phi - origin_phi), np.cos(f.phi - origin_phi))))
    )
    return {"phase_only_40step_recovery_max_abs": err, "ok": err < 1e-8}


def main() -> dict:
    t0 = time.time()
    results = {}
    jobs = [
        ("chunked_inverse", experiment_chunked_inverse),
        ("single_state_decode", experiment_true_single_state_decode),
        ("walled_attractor", experiment_walled_attractor),
        ("frozen_entanglement", experiment_frozen_entanglement),
        ("gamma_zero", experiment_gamma_zero),
        ("recursion_25", experiment_recursion),
        ("quantized_shadow", experiment_quantized_shadow),
        ("three_body_phase", experiment_three_body_phase),
        ("beta_threshold", experiment_beta_threshold),
        ("phase_only_information", experiment_phase_only_information),
    ]
    for name, fn in jobs:
        print(f"[run] {name}", flush=True)
        t1 = time.time()
        try:
            results[name] = fn()
            results[name]["_seconds"] = time.time() - t1
            print(f"      ok  {results[name]['_seconds']:.2f}s", flush=True)
        except Exception as e:
            results[name] = {"error": repr(e)}
            print(f"      FAIL {e!r}", flush=True)
    results["_meta"] = {"elapsed_s": time.time() - t0}
    out = ROOT / "results" / "deeper.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(_jsonable(results), indent=2))
    print(f"wrote {out}", flush=True)
    return results


if __name__ == "__main__":
    main()
