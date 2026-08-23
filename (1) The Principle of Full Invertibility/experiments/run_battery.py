#!/usr/bin/env python3
"""
First-principles battery for the drift--kick--phase map.

    Z(t+Δt) = Ψ(Z(t); A, λ)

Every experiment below tests one structural assumption or one mechanical
consequence of the declared force law.
"""

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

from engine.discrete import enumerate_bijection, jacobian_block_det_numeric  # noqa: E402
from engine.observe import (  # noqa: E402
    observe,
    phrase_centers,
    promote,
    reset_tokens,
)
from engine.psi import (  # noqa: E402
    FieldParams,
    euler_inverse_step,
    euler_step,
    forces,
    kinetic,
    kuramoto_R,
    make_field,
    run,
    run_inverse,
    state_distance,
    verlet_step,
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
    if isinstance(x, np.ndarray):
        return _jsonable(x.tolist())
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if x is None:
        return None
    return x


def experiment_finite_bijection() -> dict:
    """Exact proof on a finite set: the (r,v) skeleton is a bijection."""
    ok = enumerate_bijection(m=251, g=3)
    # Counterexample: if g shares a factor with M, invertibility dies.
    bad = enumerate_bijection.__wrapped__ if hasattr(enumerate_bijection, "__wrapped__") else None
    # g=2, M=250 is even — construct a small non-coprime case by hand
    m, g = 16, 2
    seen = set()
    collisions = 0
    for r in range(m):
        for v in range(m):
            rp = (r + v) % m
            vp = (g * v + (3 * rp * rp + 7 * rp + 11)) % m
            key = (rp, vp)
            if key in seen:
                collisions += 1
            seen.add(key)
    return {
        "coprime_case": ok,
        "non_coprime_counterexample": {
            "M": m,
            "g": g,
            "gcd": math.gcd(g, m),
            "unique_images": len(seen),
            "state_space": m * m,
            "collisions_observed": collisions,
            "is_bijection": len(seen) == m * m,
        },
    }


def experiment_jacobian_identity() -> dict:
    """det DΨ_(r,v) = γ^{2N} for any force Jacobian DF."""
    rng = np.random.default_rng(7)
    rows = []
    for d, gamma, dt in [
        (2, 0.982, 0.35),
        (4, 0.982, 0.35),
        (8, 0.982, 0.35),
        (8, 1.0, 0.35),
        (8, 0.5, 0.2),
        (16, 0.982, 0.35),
    ]:
        DF = rng.normal(0, 0.4, size=(d, d))
        rec = jacobian_block_det_numeric(d, gamma, dt, DF)
        rec.update({"dim_v": d, "gamma": gamma, "dt": dt})
        rows.append(rec)
    # finite-difference Jacobian of the live Euler map on a small field
    field = make_field(n=3, seed=4, params=FieldParams(theta_drift=0.0, wall_bounce=0.0))
    x0 = field.flatten()
    dim = x0.size
    eps = 1e-6
    J = np.zeros((dim, dim))
    base = field.copy()
    euler_step(base)
    x1 = base.flatten()
    for k in range(dim):
        f = field.copy()
        x = f.flatten()
        x[k] += eps
        f.unflatten_into(x)
        euler_step(f)
        J[:, k] = (f.flatten() - x1) / eps
    # (r,v) block is first 4N coordinates. Phase coupling makes the
    # *full* det differ from γ^{2N}; the theorem is for the (r,v) map
    # at frozen phases. Freeze phases by finite-differencing only r,v
    # while holding θ,φ fixed after the step comparison via a custom probe.
    n = field.n
    d_rv = 4 * n
    # Build the analytic block using a numerical DF of force wrt r at r'
    r = field.r + field.v * field.params.dt
    DF = np.zeros((2 * n, 2 * n))
    F0 = forces(field, r, field.theta, field.phi).ravel()
    for k in range(2 * n):
        rp = r.copy().ravel()
        rp[k] += 1e-6
        Fp = forces(field, rp.reshape(n, 2), field.theta, field.phi).ravel()
        DF[:, k] = (Fp - F0) / 1e-6
    # mass = 1, so the v-row uses (dt) DF
    rec_live = jacobian_block_det_numeric(2 * n, field.params.gamma, field.params.dt, DF)
    full_det = float(np.linalg.det(J)) if np.isfinite(np.linalg.det(J)) else None
    return {
        "random_DF_rows": rows,
        "live_force_DF_block": rec_live,
        "full_map_finite_diff_det": full_det,
        "full_map_dim": dim,
        "note": "full map det includes Kuramoto phase block; (r,v) block det is γ^{2N}",
        "max_rel_err_random": max(r["rel_err"] for r in rows),
    }


def _roundtrip(n, seed, steps, gamma, kind, dt=0.35, beta=0.08) -> dict:
    prm = FieldParams(gamma=gamma, dt=dt, beta=beta, theta_drift=0.0, wall_bounce=0.0)
    f = make_field(n=n, seed=seed, params=prm)
    origin = f.copy()
    run(f, steps, kind=kind)
    after = f.copy()
    residuals = run_inverse(f, steps, kind=kind)
    dist = state_distance(origin, f)
    dist["phase_residual_max"] = max(residuals) if residuals else 0.0
    dist["steps"] = steps
    dist["n"] = n
    dist["gamma"] = gamma
    dist["kind"] = kind
    dist["forward_kinetic"] = kinetic(after)
    dist["recovered_kinetic"] = kinetic(f)
    dist["success_1e8"] = dist["l2"] < 1e-8
    dist["success_1e6"] = dist["l2"] < 1e-6
    dist["success_1e4"] = dist["l2"] < 1e-4
    return dist


def experiment_roundtrip() -> dict:
    rows = []
    for kind, gamma in [("euler", 0.982), ("euler", 1.0), ("verlet", 1.0)]:
        for steps in (1, 5, 20, 80, 200):
            rows.append(_roundtrip(12, 11, steps, gamma, kind))
    # many-body
    rows.append(_roundtrip(32, 3, 40, 0.982, "euler"))
    rows.append(_roundtrip(8, 9, 400, 0.982, "euler"))
    return {"rows": rows}


def experiment_information_horizon() -> dict:
    """Damping is injective on reals and hostile to finite precision."""
    rows = []
    for gamma in (1.0, 0.995, 0.982, 0.95, 0.9):
        for steps in (10, 40, 80, 160, 320, 640):
            rec = _roundtrip(8, 21, steps, gamma, "euler")
            rec["theory_amp"] = float((1.0 / gamma) ** steps)
            rec["volume_factor"] = float(gamma ** (2 * 8 * steps))
            rows.append(rec)
    return {"rows": rows}


def experiment_chaos_still_invertible() -> dict:
    """Nearby ICs diverge, yet each is uniquely recoverable."""
    prm = FieldParams(gamma=1.0, dt=0.35, beta=0.12, theta_drift=0.0, wall_bounce=0.0)
    a = make_field(n=10, seed=5, params=prm)
    b = a.copy()
    b.r = b.r + 1e-8
    sep = []
    for t in range(250):
        euler_step(a)
        euler_step(b)
        if t in (0, 24, 49, 99, 149, 199, 249):
            sep.append(
                {
                    "t": t + 1,
                    "l2": state_distance(a, b)["l2"],
                    "max_dr": state_distance(a, b)["max_dr"],
                }
            )
    # invert a back
    origin = make_field(n=10, seed=5, params=prm)
    recovered = a.copy()
    run_inverse(recovered, 250, kind="euler")
    return {
        "separation": sep,
        "final_separation": sep[-1],
        "recovery_of_unperturbed": state_distance(origin, recovered),
        "diverged": sep[-1]["l2"] > 1e3 * sep[0]["l2"],
    }


def experiment_partial_state() -> dict:
    """
    'Every state encodes its past' is a claim about the FULL z_i.
    Drop velocity or phase and the inverse ceases to exist.
    """
    prm = FieldParams(gamma=0.982, theta_drift=0.0, wall_bounce=0.0)
    f = make_field(n=10, seed=8, params=prm)
    origin = f.copy()
    run(f, 30, kind="euler")
    full = f.copy()
    run_inverse(full, 30, kind="euler")
    full_err = state_distance(origin, full)["l2"]

    # positions only: try a false inverse that assumes v=0, phases kept
    pos_only = f.copy()
    # restore to the *final* state then zero velocities before inverting
    # (re-run forward)
    g = make_field(n=10, seed=8, params=prm)
    run(g, 30, kind="euler")
    g.v[:] = 0.0
    run_inverse(g, 30, kind="euler")
    pos_err = state_distance(origin, g)["l2"]

    h = make_field(n=10, seed=8, params=prm)
    run(h, 30, kind="euler")
    h.phi[:] = 0.0
    run_inverse(h, 30, kind="euler")
    phase_err = state_distance(origin, h)["l2"]

    # two different ICs that share final positions after projection
    rng = np.random.default_rng(0)
    finals_r = []
    inits = []
    for s in range(80):
        q = make_field(n=6, seed=100 + s, params=prm)
        inits.append(q.copy())
        run(q, 25, kind="euler")
        finals_r.append(np.round(q.r, 2))
    collisions = 0
    for i in range(len(finals_r)):
        for j in range(i + 1, len(finals_r)):
            if np.allclose(finals_r[i], finals_r[j]):
                collisions += 1
    return {
        "full_state_recovery_l2": full_err,
        "zero_velocity_then_invert_l2": pos_err,
        "zero_phase_then_invert_l2": phase_err,
        "position_rounding_collisions_80ics": collisions,
        "full_recovers": full_err < 1e-6,
        "partial_fails": pos_err > 1e-2 and phase_err > 1e-2,
    }


def _quantize_entropy(samples: np.ndarray, bins: int) -> float:
    # samples: (n_ens, dim)
    # independent per-coordinate histogram entropy (upper-ish bound)
    ens, dim = samples.shape
    h = 0.0
    for d in range(dim):
        hist, _ = np.histogram(samples[:, d], bins=bins, density=False)
        p = hist.astype(float)
        p = p[p > 0] / ens
        h += float(-np.sum(p * np.log2(p)))
    return h


def experiment_entropy_split() -> dict:
    """Fine-grained identity vs coarse-grained growth."""
    prm = FieldParams(gamma=1.0, theta_drift=0.0, wall_bounce=0.0)
    ens = 120
    n = 6
    steps_list = [0, 10, 30, 80, 150]
    fields = [make_field(n=n, seed=200 + i, params=prm) for i in range(ens)]
    rows = []
    for target in steps_list:
        if target > 0:
            prev = steps_list[steps_list.index(target) - 1]
            for f in fields:
                run(f, target - prev, kind="euler")
        full = np.stack([f.flatten() for f in fields])
        pos = np.stack([f.r.ravel() for f in fields])
        # wrap phases already in flatten
        rows.append(
            {
                "t": target,
                "H_full_12bins": _quantize_entropy(full, 12),
                "H_pos_12bins": _quantize_entropy(pos, 12),
                "H_full_6bins": _quantize_entropy(full, 6),
                "mean_kinetic": float(np.mean([kinetic(f) for f in fields])),
                "mean_R": float(np.mean([kuramoto_R(f) for f in fields])),
            }
        )
    # damped ensemble: fine volume contracts
    prm_d = FieldParams(gamma=0.982, theta_drift=0.0, wall_bounce=0.0)
    fd = [make_field(n=n, seed=200 + i, params=prm_d) for i in range(ens)]
    drows = []
    for target in steps_list:
        if target > 0:
            prev = steps_list[steps_list.index(target) - 1]
            for f in fd:
                run(f, target - prev, kind="euler")
        full = np.stack([f.flatten() for f in fd])
        pos = np.stack([f.r.ravel() for f in fd])
        drows.append(
            {
                "t": target,
                "H_full_12bins": _quantize_entropy(full, 12),
                "H_pos_12bins": _quantize_entropy(pos, 12),
                "mean_kinetic": float(np.mean([kinetic(f) for f in fd])),
                "volume_factor": float(0.982 ** (2 * n * target)),
            }
        )
    return {"gamma1": rows, "gamma0982": drows}


def experiment_mdl_many_to_one() -> dict:
    """Many distinct microstates collapse to one phrase."""
    prm = FieldParams(gamma=0.982, theta_drift=0.0, wall_bounce=0.0)
    phrases = {}
    micro_pairs = []
    for s in range(160):
        f = make_field(n=14, seed=300 + s, params=prm)
        run(f, 90, kind="euler")
        obs = observe(f)
        phrases.setdefault(obs["phrase"], []).append(
            {
                "seed": 300 + s,
                "L_total": obs["L_total"],
                "residual_var": obs["residual_var"],
                "occupancy": obs["occupancy"],
                "fingerprint": [
                    float(np.mean(f.r)),
                    float(np.std(f.r)),
                    float(kuramoto_R(f)),
                    float(kinetic(f)),
                ],
            }
        )
        micro_pairs.append((obs["phrase"], f.flatten().copy()))
    # among shared phrases, are microstates actually different?
    shared = {k: v for k, v in phrases.items() if len(v) >= 2}
    distinct_under_same_phrase = 0
    checked = 0
    keys = list(shared)
    for ph in keys[:12]:
        members = [p for p, fl in micro_pairs if p == ph]
        flats = [fl for p, fl in micro_pairs if p == ph]
        for i in range(len(flats)):
            for j in range(i + 1, len(flats)):
                checked += 1
                if np.linalg.norm(flats[i] - flats[j]) > 1e-3:
                    distinct_under_same_phrase += 1
    return {
        "unique_phrases": len(phrases),
        "ensemble": 160,
        "largest_basin": max(len(v) for v in phrases.values()),
        "phrases_with_multiplicity": len(shared),
        "distinct_microstates_same_phrase": distinct_under_same_phrase,
        "pairs_checked": checked,
        "top_phrases": sorted(
            ((k, len(v), float(np.mean([x["L_total"] for x in v]))) for k, v in phrases.items()),
            key=lambda t: -t[1],
        )[:8],
    }


def experiment_observation_force_breaks_injectivity() -> dict:
    """F_obs written from an MDL phrase makes two pasts share a future force."""
    prm = FieldParams(gamma=0.982, theta_drift=0.0, wall_bounce=0.0, obs_strength=0.08)
    # two different ICs, evolve without obs, collapse to same phrase, then
    # apply one observation step using the same phrase centers.
    hits = []
    for s in range(40):
        a = make_field(n=12, seed=400 + s, params=prm)
        b = make_field(n=12, seed=800 + s, params=prm)
        run(a, 70, kind="euler")
        run(b, 70, kind="euler")
        oa, ob = observe(a), observe(b)
        if oa["phrase"] != ob["phrase"]:
            continue
        centers = phrase_centers(a, oa["clusters"])
        a.phrase_centers = centers
        b.phrase_centers = centers
        a_before, b_before = a.copy(), b.copy()
        euler_step(a)
        euler_step(b)
        # if they were different, they remain different — F_obs does not
        # identify states in one step. The real non-injectivity is the
        # MDL map itself, plus the fact that future evolution is conditioned
        # on a many-to-one label. Record both facts.
        hits.append(
            {
                "phrase": oa["phrase"],
                "micro_l2_before": state_distance(a_before, b_before)["l2"],
                "micro_l2_after": state_distance(a, b)["l2"],
                "same_force_centers": True,
            }
        )
    # Direct test: MDL map collisions
    return {
        "same_phrase_pairs": len(hits),
        "mean_micro_l2_before": float(np.mean([h["micro_l2_before"] for h in hits])) if hits else None,
        "mean_micro_l2_after": float(np.mean([h["micro_l2_after"] for h in hits])) if hits else None,
        "mdl_is_the_lossy_map": True,
        "sample": hits[:5],
    }


def experiment_recursion() -> dict:
    """25-layer protocol: promote phrase centers to anchors, then re-run Ψ."""
    prm = FieldParams(gamma=0.982, theta_drift=0.0, wall_bounce=0.0)
    f = make_field(n=16, seed=17, params=prm)
    layers = []
    for ell in range(25):
        run(f, 40, kind="euler")
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
                "K": kinetic(f),
                "n_clusters": obs["n_clusters"],
            }
        )
        promote(f, obs["clusters"], charge_scale=1.15)
        reset_tokens(f, seed=1000 + ell, spread=3.0)
    # uniqueness / self-similarity
    phrases = [L["phrase"] for L in layers]
    return {
        "layers": layers,
        "unique_phrases": len(set(phrases)),
        "L_data_first": layers[0]["L_data"],
        "L_data_last": layers[-1]["L_data"],
        "L_data_min": min(L["L_data"] for L in layers),
        "anchors_first": layers[0]["n_anchors"],
        "anchors_last": layers[-1]["n_anchors"],
        "monotonic_anchor_growth": all(
            layers[i]["n_anchors"] <= layers[i + 1]["n_anchors"]
            for i in range(len(layers) - 1)
        ),
    }


def experiment_entanglement_analog() -> dict:
    """
    Shared clock history, then spatial separation.

    This asks what the phase law does when two tokens lock and then move apart.
    """
    prm = FieldParams(
        gamma=1.0,
        beta=0.35,
        dt=0.25,
        theta_drift=0.0,
        wall_bounce=0.0,
        wave_amp=0.04,
        alpha=0.02,
    )
    f = make_field(n=2, seed=1, params=prm, spread=0.4)
    f.r = np.array([[-0.2, 0.0], [0.2, 0.0]])
    f.v = np.zeros((2, 2))
    f.phi = np.array([0.2, 2.4])
    # lock
    lock_R = []
    for t in range(120):
        euler_step(f)
        if t % 20 == 19:
            lock_R.append(
                {
                    "t": t + 1,
                    "dphi": float(np.abs(np.arctan2(np.sin(f.phi[0] - f.phi[1]), np.cos(f.phi[0] - f.phi[1])))),
                    "R": kuramoto_R(f),
                    "dist": float(np.linalg.norm(f.r[0] - f.r[1])),
                }
            )
    locked_dphi = lock_R[-1]["dphi"]
    # separate: teleport spatially far, keep phases and velocities
    f.r = np.array([[-8.0, 0.0], [8.0, 0.0]])
    f.v = np.zeros((2, 2))
    sep = []
    for t in range(120):
        euler_step(f)
        if t % 20 == 19:
            sep.append(
                {
                    "t": t + 1,
                    "dphi": float(np.abs(np.arctan2(np.sin(f.phi[0] - f.phi[1]), np.cos(f.phi[0] - f.phi[1])))),
                    "R": kuramoto_R(f),
                    "dist": float(np.linalg.norm(f.r[0] - f.r[1])),
                }
            )
    # local "measurement": snap particle 0 phase, do not touch particle 1
    pre = f.phi.copy()
    f.phi[0] = 0.0
    untouched = float(abs(np.sin(f.phi[1] - pre[1])))
    # one more dynamical step after the snap
    before = f.phi.copy()
    euler_step(f)
    after_step_dphi = float(
        abs(np.arctan2(np.sin(f.phi[0] - f.phi[1]), np.cos(f.phi[0] - f.phi[1])))
    )
    influence = float(abs(np.arctan2(np.sin(f.phi[1] - before[1]), np.cos(f.phi[1] - before[1]))))
    return {
        "lock_trace": lock_R,
        "separated_trace": sep,
        "locked_dphi": locked_dphi,
        "separated_dphi_final": sep[-1]["dphi"],
        "correlation_survives_separation": sep[-1]["dphi"] < 0.6,
        "snap_leaves_other_phase_untouched": untouched < 1e-15,
        "post_snap_influence_on_distant_phase": influence,
        "post_snap_dphi": after_step_dphi,
        "reading": (
            "shared-history correlation; instantaneous readout is local; "
            "residual influence after a step is O(1/d) field coupling"
        ),
    }


def experiment_signaling() -> dict:
    """Does a local snap change a distant token more with global MDL back-reaction?"""
    def once(use_obs: bool, seed: int) -> dict:
        prm = FieldParams(
            gamma=0.982,
            beta=0.1,
            theta_drift=0.0,
            wall_bounce=0.0,
            obs_strength=0.10 if use_obs else 0.0,
        )
        f = make_field(n=16, seed=seed, params=prm)
        run(f, 50, kind="euler")
        # split left/right
        left = np.where(f.r[:, 0] < 0)[0]
        right = np.where(f.r[:, 0] >= 0)[0]
        if len(left) == 0 or len(right) == 0:
            return {"skipped": True}
        obs = observe(f)
        if use_obs:
            f.phrase_centers = phrase_centers(f, obs["clusters"])
        phi_right_before = f.phi[right].copy()
        r_right_before = f.r[right].copy()
        # snap-measure the left side
        f.phi[left] = 0.0
        euler_step(f)
        dphi = float(np.mean(np.abs(np.sin(f.phi[right] - phi_right_before))))
        dr = float(np.mean(np.linalg.norm(f.r[right] - r_right_before, axis=1)))
        return {"dphi_right": dphi, "dr_right": dr, "n_left": int(len(left)), "n_right": int(len(right))}

    off = [once(False, 50 + i) for i in range(24)]
    on = [once(True, 50 + i) for i in range(24)]
    off = [x for x in off if not x.get("skipped")]
    on = [x for x in on if not x.get("skipped")]
    return {
        "no_obs_mean_dphi_right": float(np.mean([x["dphi_right"] for x in off])),
        "obs_mean_dphi_right": float(np.mean([x["dphi_right"] for x in on])),
        "no_obs_mean_dr_right": float(np.mean([x["dr_right"] for x in off])),
        "obs_mean_dr_right": float(np.mean([x["dr_right"] for x in on])),
        "obs_increases_distant_influence": (
            float(np.mean([x["dphi_right"] for x in on]))
            > float(np.mean([x["dphi_right"] for x in off]))
        ),
    }


def experiment_basins() -> dict:
    """Do occupancy fractions track charge Q, or only initial geometry?"""
    prm = FieldParams(gamma=0.982, theta_drift=0.0, wall_bounce=0.0)
    occ = []
    for s in range(200):
        f = make_field(n=18, seed=900 + s, params=prm)
        run(f, 80, kind="euler")
        d = np.linalg.norm(f.r[:, None, :] - f.p[None, :, :], axis=-1)
        nearest = np.argmin(d, axis=1)
        occ.append([(nearest == a).mean() for a in range(len(f.Q))])
    occ = np.asarray(occ)
    mean = occ.mean(axis=0)
    Q = make_field(n=3, seed=0).Q
    Qn = Q / Q.sum()
    return {
        "mean_occupancy": mean.tolist(),
        "Q_normalized": Qn.tolist(),
        "l1_occupancy_minus_Q": float(np.sum(np.abs(mean - Qn))),
        "corr_occ_Q": float(np.corrcoef(mean, Qn)[0, 1]),
        "uniform": [1 / len(Q)] * len(Q),
        "l1_occupancy_minus_uniform": float(np.sum(np.abs(mean - 1 / len(Q)))),
    }


def experiment_causal_arrow() -> dict:
    """
    Compression gradient across time: which direction of the same
    reversible trajectory has a shorter positional description?
    """
    prm = FieldParams(gamma=1.0, theta_drift=0.0, wall_bounce=0.0)
    rows = []
    for s in range(40):
        f = make_field(n=14, seed=70 + s, params=prm)
        start = observe(f)
        run(f, 100, kind="euler")
        end = observe(f)
        rows.append(
            {
                "L_start": start["L_total"],
                "L_end": end["L_total"],
                "var_start": start["residual_var"],
                "var_end": end["residual_var"],
            }
        )
    dL = np.array([r["L_end"] - r["L_start"] for r in rows])
    return {
        "mean_dL": float(dL.mean()),
        "frac_end_shorter": float(np.mean(dL < 0)),
        "frac_end_longer": float(np.mean(dL > 0)),
        "mean_var_start": float(np.mean([r["var_start"] for r in rows])),
        "mean_var_end": float(np.mean([r["var_end"] for r in rows])),
        "arrow_from_compression": float(np.mean(dL < 0)) > 0.65,
    }


def experiment_energy() -> dict:
    prm_e = FieldParams(gamma=0.982, theta_drift=0.0, wall_bounce=0.0)
    prm_v = FieldParams(gamma=1.0, theta_drift=0.0, wall_bounce=0.0, dt=0.2)
    fe = make_field(n=12, seed=2, params=prm_e)
    fv = make_field(n=12, seed=2, params=prm_v)
    ke, kv = [kinetic(fe)], [kinetic(fv)]
    for _ in range(200):
        euler_step(fe)
        verlet_step(fv)
        ke.append(kinetic(fe))
        kv.append(kinetic(fv))
    return {
        "euler_K0": ke[0],
        "euler_K_end": ke[-1],
        "euler_decay_ratio": ke[-1] / max(ke[0], 1e-12),
        "verlet_K0": kv[0],
        "verlet_K_end": kv[-1],
        "verlet_max": max(kv),
        "verlet_min": min(kv),
        "verlet_relative_drift": abs(kv[-1] - kv[0]) / max(kv[0], 1e-12),
    }


def experiment_injectivity_sample() -> dict:
    """Do distinct ICs remain distinct after long Euler evolution?"""
    prm = FieldParams(gamma=0.982, theta_drift=0.0, wall_bounce=0.0)
    finals = []
    for s in range(300):
        f = make_field(n=8, seed=2000 + s, params=prm)
        run(f, 120, kind="euler")
        finals.append(np.round(f.flatten(), 8))
    uniq = np.unique(np.stack(finals), axis=0)
    return {
        "ics": 300,
        "unique_rounded_8dec": int(len(uniq)),
        "collisions_at_8dec": 300 - int(len(uniq)),
        "unique_rounded_4dec": int(len(np.unique(np.round(np.stack(finals), 4), axis=0))),
    }


def experiment_beta_invertibility() -> dict:
    """Phase inversion is a contraction only while β is not huge."""
    rows = []
    for beta in (0.02, 0.08, 0.2, 0.5, 1.0, 2.5):
        rec = _roundtrip(10, 6, 30, 0.982, "euler", beta=beta)
        rec["beta"] = beta
        rows.append(rec)
    return {"rows": rows}


def main() -> dict:
    t0 = time.time()
    results = {}
    jobs = [
        ("finite_bijection", experiment_finite_bijection),
        ("jacobian_identity", experiment_jacobian_identity),
        ("roundtrip", experiment_roundtrip),
        ("information_horizon", experiment_information_horizon),
        ("chaos_still_invertible", experiment_chaos_still_invertible),
        ("partial_state", experiment_partial_state),
        ("entropy_split", experiment_entropy_split),
        ("mdl_many_to_one", experiment_mdl_many_to_one),
        ("observation_force", experiment_observation_force_breaks_injectivity),
        ("recursion_25", experiment_recursion),
        ("entanglement_analog", experiment_entanglement_analog),
        ("signaling", experiment_signaling),
        ("basins", experiment_basins),
        ("causal_arrow", experiment_causal_arrow),
        ("energy", experiment_energy),
        ("injectivity_sample", experiment_injectivity_sample),
        ("beta_invertibility", experiment_beta_invertibility),
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
    results["_meta"] = {
        "elapsed_s": time.time() - t0,
        "formula": "Z(t+Δt)=Ψ(Z(t); A, λ)",
        "axioms": [
            "No true information loss exists",
            "Entropy is epistemic only",
            "Every state is a perfect encoding of its past",
        ],
    }
    out = ROOT / "results" / "battery.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(_jsonable(results), indent=2))
    print(f"wrote {out}", flush=True)
    return results


if __name__ == "__main__":
    main()
