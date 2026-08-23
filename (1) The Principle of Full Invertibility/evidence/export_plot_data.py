#!/usr/bin/env python3
"""Export manuscript plot tables from the stored JSON ledgers."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "evidence"


def load(path: Path):
    return json.loads(path.read_text())


def write_csv(name: str, fieldnames: list[str], rows: list[dict]) -> None:
    path = EVIDENCE / name
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {path}")


def main() -> None:
    battery = load(ROOT / "results" / "battery.json")
    deeper = load(ROOT / "results" / "deeper.json")
    omega = load(ROOT / "results" / "omega.json")
    verification = load(EVIDENCE / "thesis_verification.json")

    corrected_clock = verification["decoupled_clock_memory_and_readout"]
    eq = {row["t"]: row for row in corrected_clock["equal_omega_trace"]}
    uq = {row["t"]: row for row in corrected_clock["unequal_omega_trace"]}
    clock_rows = [
        {
            "t": t,
            "equal_delta_rad": abs(eq[t]["delta_actual"]),
            "unequal_delta_rad": abs(uq[t]["delta_actual"]),
            "equal_R": abs(math.cos(eq[t]["delta_actual"] / 2.0)),
            "unequal_R": abs(math.cos(uq[t]["delta_actual"] / 2.0)),
        }
        for t in sorted(eq)
    ]
    write_csv(
        "clock_memory.csv",
        ["t", "equal_delta_rad", "unequal_delta_rad", "equal_R", "unequal_R"],
        clock_rows,
    )
    reset_rows = corrected_clock["controlled_reset_rows"]
    write_csv(
        "clock_readout.csv",
        [
            "distance",
            "softened_distance",
            "reset_vs_control_effect",
            "analytic_prediction",
            "absolute_residual",
        ],
        reset_rows,
    )

    horizon_rows = []
    for row in verification["replicated_information_horizon"]["rows"]:
        horizon_rows.append(
            {
                "gamma": row["gamma"],
                "steps": row["steps"],
                "certified_seed_count": row["certified_seed_count"],
                "median_log10_l2": row["certified_median_log10_l2"],
                "p10_log10_l2": row["certified_p10_log10_l2"],
                "p90_log10_l2": row["certified_p90_log10_l2"],
                "fraction_below_1e6": row["certified_fraction_below_1e6"],
                "fraction_below_1e4": row["certified_fraction_below_1e4"],
            }
        )
    write_csv(
        "information_horizon.csv",
        [
            "gamma",
            "steps",
            "certified_seed_count",
            "median_log10_l2",
            "p10_log10_l2",
            "p90_log10_l2",
            "fraction_below_1e6",
            "fraction_below_1e4",
        ],
        horizon_rows,
    )
    for gamma, suffix in ((1.0, "gamma1"), (0.982, "gamma0982")):
        selected = [row for row in horizon_rows if row["gamma"] == gamma]
        write_csv(
            f"information_horizon_{suffix}.csv",
            [
                "gamma",
                "steps",
                "certified_seed_count",
                "median_log10_l2",
                "p10_log10_l2",
                "p90_log10_l2",
                "fraction_below_1e6",
                "fraction_below_1e4",
            ],
            selected,
        )

    recursion_rows = [
        {
            "layer": row["layer"],
            "anchors": row["n_anchors"],
            "residual_variance": row["residual_var"],
            "L_data": row["L_data"],
            "L_phrase": row["L_phrase"],
            "L_total": row["L_total"],
            "R": row["R"],
        }
        for row in deeper["recursion_25"]["layers"]
    ]
    write_csv(
        "recursion.csv",
        [
            "layer",
            "anchors",
            "residual_variance",
            "L_data",
            "L_phrase",
            "L_total",
            "R",
        ],
        recursion_rows,
    )

    entropy_rows = []
    fixed_entropy = verification["fixed_partition_marginal_entropy"]["rows"]
    gamma_one = {row["t"]: row for row in fixed_entropy if row["gamma"] == 1.0}
    gamma_damped = {
        row["t"]: row for row in fixed_entropy if row["gamma"] == 0.982
    }
    for t in sorted(gamma_one):
        entropy_rows.append(
            {
                "t": t,
                "full_gamma1_bits": gamma_one[t]["full_marginal_sum_bits"],
                "position_gamma1_bits": gamma_one[t][
                    "position_marginal_sum_bits"
                ],
                "full_gamma0982_bits": gamma_damped[t]["full_marginal_sum_bits"],
                "position_gamma0982_bits": gamma_damped[t][
                    "position_marginal_sum_bits"
                ],
            }
        )
    write_csv(
        "quantized_entropy.csv",
        [
            "t",
            "full_gamma1_bits",
            "position_gamma1_bits",
            "full_gamma0982_bits",
            "position_gamma0982_bits",
        ],
        entropy_rows,
    )

    lyapunov_rows = [
        {
            "t": row["t"],
            "separation": row["sep"],
            "lambda_proxy": row["lambda_proxy"],
        }
        for row in omega["lyapunov"]["trace"]
    ]
    write_csv(
        "lyapunov.csv",
        ["t", "separation", "lambda_proxy"],
        lyapunov_rows,
    )

    beta_rows = [
        {
            "beta": row["beta"],
            "fixed_point_l2": row["l2"],
            "phase_residual": row["phase_residual"],
            "solver_success": row["ok"],
        }
        for row in deeper["beta_threshold"]["rows"]
    ]
    write_csv(
        "fixed_point_beta_sweep.csv",
        ["beta", "fixed_point_l2", "phase_residual", "solver_success"],
        beta_rows,
    )

    wall_rows = [
        {
            "condition": row["label"],
            "mean_delta_L_bits": row["mean_dL"],
            "fraction_shorter": row["frac_end_shorter"],
            "mean_residual_variance": row["mean_residual_var"],
            "mean_R": row["mean_R"],
            "mean_lock_fraction": row["mean_lock_frac"],
        }
        for row in deeper["walled_attractor"]["rows"]
    ]
    write_csv(
        "boundary_conditions.csv",
        [
            "condition",
            "mean_delta_L_bits",
            "fraction_shorter",
            "mean_residual_variance",
            "mean_R",
            "mean_lock_fraction",
        ],
        wall_rows,
    )

    determinant_rows = [
        {
            "n": row["n"],
            "seed": row["seed"],
            "det_full": row["det_full_finite_difference"],
            "det_factorized": row["factorized_prediction"],
            "relative_error": row["relative_error"],
            "phase_q_bound": row["phase_contraction_q_bound"],
        }
        for row in verification["determinant_factorization"]["rows"]
    ]
    write_csv(
        "determinant_factorization.csv",
        [
            "n",
            "seed",
            "det_full",
            "det_factorized",
            "relative_error",
            "phase_q_bound",
        ],
        determinant_rows,
    )


if __name__ == "__main__":
    main()
