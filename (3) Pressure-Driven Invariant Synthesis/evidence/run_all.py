#!/usr/bin/env python3
"""Generate every Paper III metric, audit artifact, and figure."""

from __future__ import annotations

from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pdis.canonical import (  # noqa: E402
    canonical_features,
    normalized_laplacian_eigenvalues,
    recurrence_adjacency,
)
from pdis.datasets import (  # noqa: E402
    Trajectory,
    kuramoto_domain,
    logistic_domain,
    lorenz_domain,
    noisy_transform,
    sealed_ucr_domain,
    similarity_transform,
    synthetic_domains,
)
from pdis.kuramoto_posthoc import run_posthoc_diagnostic  # noqa: E402
from pdis.programs import Program, build_grammar  # noqa: E402
from pdis.statistics import collision_indicators, stable_seed  # noqa: E402
from pdis.synthesis import (  # noqa: E402
    evaluate_domain,
    paired_domain_reduction,
    representation,
    synthesize,
)


SEED = 20260822
THRESHOLD = 0.22
MAXIMUM_ROUNDS = 3
BASE_NAMES = ("radial_lag1", "path_tortuosity")
OUTPUT = ROOT / "evidence" / "outputs"
FIGURES = ROOT / "thesis" / "figures"
RAW = ROOT / "evidence" / "data" / "raw"


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def feature_table(records: list[Trajectory], mode: str = "original") -> dict[str, dict[str, float]]:
    table: dict[str, dict[str, float]] = {}
    for index, item in enumerate(records):
        if mode == "original":
            values = item.values
        elif mode == "similarity":
            values = similarity_transform(item.values, SEED + index)
        elif mode == "noise":
            values = noisy_transform(item.values, SEED + 100_000 + index)
        else:
            raise ValueError(mode)
        table[item.identifier] = canonical_features(values)
        if (index + 1) % 40 == 0 or index + 1 == len(records):
            print(f"features[{mode}] {index + 1}/{len(records)}", flush=True)
    return table


def merge_tables(*tables: dict[str, dict[str, float]]) -> dict[str, dict[str, float]]:
    result: dict[str, dict[str, float]] = {}
    for table in tables:
        overlap = set(result).intersection(table)
        if overlap:
            raise ValueError(f"duplicate feature records: {sorted(overlap)[:3]}")
        result.update(table)
    return result


def selected_union(folds: dict[str, dict[str, object]]) -> list[str]:
    return sorted(
        {
            row["name"]
            for fold in folds.values()
            for row in fold["synthesis"]["selected"]
        }
    )


def admissible_programs(
    programs: list[Program],
    records: list[Trajectory],
    table: dict[str, dict[str, float]],
    similarity_table: dict[str, dict[str, float]],
    noise_table: dict[str, dict[str, float]],
) -> list[Program]:
    accepted: list[Program] = []
    for program in programs:
        original = np.asarray([program.evaluate(table[item.identifier]) for item in records])
        similar = np.asarray([program.evaluate(similarity_table[item.identifier]) for item in records])
        noisy = np.asarray([program.evaluate(noise_table[item.identifier]) for item in records])
        invariance = np.max(np.abs(original - similar) / (1.0 + np.abs(original)))
        noise = np.quantile(np.abs(original - noisy) / (1.0 + np.abs(original)), 0.95)
        if invariance <= 5.0e-6 and noise <= 0.20:
            accepted.append(program)
    return accepted


def random_program_null(
    programs: list[Program],
    audit_records: list[Trajectory],
    audit_table: dict[str, dict[str, float]],
    observed: float,
    *,
    repetitions: int = 999,
) -> dict[str, object]:
    rng = np.random.default_rng(stable_seed(audit_records[0].domain + "-random-program"))
    sampled = rng.choice(len(programs), size=repetitions, replace=True)
    reductions = []
    for index in sampled:
        result = paired_domain_reduction(
            audit_records,
            audit_table,
            BASE_NAMES,
            (programs[int(index)],),
            threshold=THRESHOLD,
        )
        reductions.append(float(result["absolute_reduction"]))
    null = np.asarray(reductions)
    return {
        "repetitions": repetitions,
        "admissible_program_count": len(programs),
        "p_value": (1.0 + float(np.sum(null >= observed))) / (repetitions + 1.0),
        "null_mean": float(np.mean(null)),
        "null_95": [float(np.quantile(null, 0.025)), float(np.quantile(null, 0.975))],
    }


def full_pipeline_label_null(
    sources: dict[str, list[Trajectory]],
    audit_records: list[Trajectory],
    source_table: dict[str, dict[str, float]],
    similarity_table: dict[str, dict[str, float]],
    noise_table: dict[str, dict[str, float]],
    audit_table: dict[str, dict[str, float]],
    grammar: list[Program],
    observed: float,
    *,
    repetitions: int = 99,
) -> dict[str, object]:
    rng = np.random.default_rng(stable_seed(audit_records[0].domain + "-pipeline-null"))
    reductions = []
    selected_counts = []
    for _ in range(repetitions):
        permuted_sources: dict[str, list[Trajectory]] = {}
        for name, records in sources.items():
            shuffled = rng.permutation([item.label for item in records])
            permuted_sources[name] = [replace(item, label=int(label)) for item, label in zip(records, shuffled)]
        result = synthesize(
            permuted_sources,
            source_table,
            similarity_table,
            noise_table,
            grammar,
            base_names=BASE_NAMES,
            threshold=THRESHOLD,
            maximum_rounds=MAXIMUM_ROUNDS,
        )
        selected_counts.append(len(result.selected))
        audit = paired_domain_reduction(
            audit_records,
            audit_table,
            BASE_NAMES,
            result.selected,
            threshold=THRESHOLD,
        )
        reductions.append(float(audit["absolute_reduction"]))
    null = np.asarray(reductions)
    return {
        "repetitions": repetitions,
        "maximum_rounds": MAXIMUM_ROUNDS,
        "p_value": (1.0 + float(np.sum(null >= observed))) / (repetitions + 1.0),
        "null_mean": float(np.mean(null)),
        "null_95": [float(np.quantile(null, 0.025)), float(np.quantile(null, 0.975))],
        "fraction_selecting_any_program": float(np.mean(np.asarray(selected_counts) > 0)),
    }


def grammar_ablations(grammar: list[Program]) -> dict[str, list[Program]]:
    def operands(program: Program) -> tuple[str, ...]:
        return (program.left,) if program.right is None else (program.left, program.right)

    return {
        "primitive_only": [program for program in grammar if program.complexity == 1],
        "recurrence_only": [
            program for program in grammar if all(name.startswith("recurrence_") for name in operands(program))
        ],
        "topology_only": [
            program for program in grammar if all(name.startswith(("h0_", "h1_")) for name in operands(program))
        ],
    }


def null_p_value(
    records: list[Trajectory],
    table: dict[str, dict[str, float]],
    base_names: tuple[str, ...],
    selected: tuple[Program, ...],
    observed: float,
    repetitions: int = 999,
) -> dict[str, object]:
    labels = np.asarray([item.label for item in records], dtype=int)
    base_rep = representation(records, table, base_names, ())
    final_rep = representation(records, table, base_names, selected)
    rng = np.random.default_rng(stable_seed(records[0].domain + "-null"))
    null = []
    for index in range(repetitions):
        permuted = rng.permutation(labels)
        seed = stable_seed(records[0].domain) + index + 10_000
        before = collision_indicators(base_rep, permuted, threshold=THRESHOLD, seed=seed)
        after = collision_indicators(final_rep, permuted, threshold=THRESHOLD, seed=seed)
        null.append(float(np.mean(before) - np.mean(after)))
    null_array = np.asarray(null)
    p_value = (1.0 + float(np.sum(null_array >= observed))) / (repetitions + 1.0)
    return {
        "repetitions": repetitions,
        "p_value": p_value,
        "null_mean": float(np.mean(null_array)),
        "null_95": [float(np.quantile(null_array, 0.025)), float(np.quantile(null_array, 0.975))],
    }


def theorem_checks(
    records: list[Trajectory],
    table: dict[str, dict[str, float]],
    similarity_table: dict[str, dict[str, float]],
    selected: tuple[Program, ...],
) -> dict[str, object]:
    adjacency_equal = []
    for index, item in enumerate(records[:24]):
        original, _, _ = recurrence_adjacency(item.values)
        transformed, _, _ = recurrence_adjacency(similarity_transform(item.values, SEED + index))
        adjacency_equal.append(bool(np.array_equal(original, transformed)))

    program_errors: dict[str, float] = {}
    for program in selected:
        errors = []
        for item in records:
            a = program.evaluate(table[item.identifier])
            b = program.evaluate(similarity_table[item.identifier])
            errors.append(abs(a - b) / (1.0 + abs(a)))
        program_errors[program.name] = float(max(errors))

    # Known graph checks: P4 has lambda_2=0.5; 2P2 is disconnected so lambda_2=0.
    p4 = np.zeros((4, 4))
    for left, right in ((0, 1), (1, 2), (2, 3)):
        p4[left, right] = p4[right, left] = 1.0
    two_p2 = np.zeros((4, 4))
    two_p2[0, 1] = two_p2[1, 0] = 1.0
    two_p2[2, 3] = two_p2[3, 2] = 1.0
    p4_eig = normalized_laplacian_eigenvalues(p4)
    two_p2_eig = normalized_laplacian_eigenvalues(two_p2)

    return {
        "similarity_adjacency_trials": len(adjacency_equal),
        "similarity_adjacency_exact": int(sum(adjacency_equal)),
        "selected_program_max_relative_error": program_errors,
        "lambda2_path4": float(p4_eig[1]),
        "lambda2_disconnected_2p2": float(two_p2_eig[1]),
        "lambda2_indexing_pass": bool(abs(p4_eig[1] - 0.5) < 1.0e-12 and abs(two_p2_eig[1]) < 1.0e-12),
    }


def noncompleteness_check(item: Trajectory) -> dict[str, object]:
    original, _, _ = recurrence_adjacency(item.values)
    reversed_adjacency, _, _ = recurrence_adjacency(item.values[::-1])
    # Reversal changes vertex order. Spectra remain equal although oriented histories differ.
    a = normalized_laplacian_eigenvalues(original)
    b = normalized_laplacian_eigenvalues(reversed_adjacency)
    return {
        "trajectory": item.identifier,
        "histories_equal": bool(np.array_equal(item.values, item.values[::-1])),
        "recurrence_spectrum_max_abs_difference": float(np.max(np.abs(a - b))),
        "interpretation": "time orientation is not identified by the recurrence spectrum",
    }


def plot_transfer(results: pd.DataFrame) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    order = list(results["domain"])
    x = np.arange(len(order))
    width = 0.36
    fig, ax = plt.subplots(figsize=(8.0, 3.8))
    ax.bar(x - width / 2, results["risk_before"], width, label="base representation", color="#69747D")
    ax.bar(x + width / 2, results["risk_after"], width, label="after frozen synthesis", color="#1D5F8A")
    ax.set_ylabel("cross-witness collision risk")
    ax.set_xticks(x, order, rotation=18, ha="right")
    ax.set_ylim(0, max(0.05, float(results[["risk_before", "risk_after"]].to_numpy().max()) * 1.18))
    ax.grid(axis="y", color="#D5DDE1", linewidth=0.7)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, ncol=2, loc="upper right")
    fig.tight_layout()
    fig.savefig(FIGURES / "frozen_transfer.pdf")
    fig.savefig(FIGURES / "frozen_transfer.png", dpi=220)
    plt.close(fig)


def plot_programs(result: dict[str, object]) -> None:
    selected = result["synthesis_all_synthetic"]["selected"]
    names = [row["name"] for row in selected]
    if not names:
        return
    real = result["sealed_real_domains"]
    domains = list(real)
    matrix = np.asarray([[real[domain]["program_auc"].get(name, 0.5) for name in names] for domain in domains])
    fig, ax = plt.subplots(figsize=(7.4, 2.8 + 0.35 * len(domains)))
    image = ax.imshow(matrix, vmin=0.5, vmax=1.0, cmap="Blues", aspect="auto")
    ax.set_xticks(np.arange(len(names)), [name.replace("normalized_difference", "ndiff") for name in names], rotation=20, ha="right")
    ax.set_yticks(np.arange(len(domains)), domains)
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            ax.text(column, row, f"{matrix[row, column]:.2f}", ha="center", va="center", color="#17232D")
    fig.colorbar(image, ax=ax, label="orientation-free AUC")
    fig.tight_layout()
    fig.savefig(FIGURES / "program_auc.pdf")
    fig.savefig(FIGURES / "program_auc.png", dpi=220)
    plt.close(fig)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    discovery = synthetic_domains(n_per_class=28)
    discovery_records = [item for records in discovery.values() for item in records]
    archives = {
        "ECGFiveDays": RAW / "ECGFiveDays.zip",
        "Earthquakes": RAW / "Earthquakes.zip",
    }
    for name, path in archives.items():
        if not path.exists():
            raise FileNotFoundError(f"missing sealed archive {name}: {path}")

    # Only source/discovery records are materialized before program freeze.
    print(f"discovery records: {len(discovery_records)}", flush=True)
    discovery_table = feature_table(discovery_records, "original")
    discovery_similarity = feature_table(discovery_records, "similarity")
    discovery_noise = feature_table(discovery_records, "noise")

    feature_names = sorted(next(iter(discovery_table.values())))
    grammar = build_grammar(feature_names, BASE_NAMES)
    grammar_manifest = {
        "status": "internal_freeze_not_externally_timestamped",
        "seed": SEED,
        "base_names": list(BASE_NAMES),
        "collision_threshold": THRESHOLD,
        "maximum_rounds": MAXIMUM_ROUNDS,
        "minimum_audit_reduction": 0.15,
        "minimum_audit_auc": 0.75,
        "minimum_selection_auc": 0.70,
        "nuisance_group": "translation x positive global scale x orthogonal channel action",
        "program_count": len(grammar),
        "programs": [program.to_dict() for program in grammar],
    }
    write_json(OUTPUT / "grammar_manifest.json", grammar_manifest)
    declared_split_manifest = {
        "discovery_synthetic": {
            name: {"count": len(records), "identifiers": [item.identifier for item in records]}
            for name, records in discovery.items()
        },
        "independent_synthetic_audit": {
            "per_class": 120,
            "seeds": {"logistic": SEED + 191, "lorenz": SEED + 192, "kuramoto": SEED + 193},
        },
        "sealed_real": {
            "archives": ["ECGFiveDays", "Earthquakes"],
            "split": "TEST only",
            "balanced_per_class": 32,
        },
        "archives": {name: {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)} for name, path in archives.items()},
        "independence_unit": "one independently seeded trajectory or one non-overlapping archive case",
    }
    write_json(OUTPUT / "declared_split_manifest.json", declared_split_manifest)

    folds: dict[str, dict[str, object]] = {}
    fold_objects: dict[str, object] = {}
    for holdout in discovery:
        sources = {name: records for name, records in discovery.items() if name != holdout}
        synthesis = synthesize(
            sources,
            discovery_table,
            discovery_similarity,
            discovery_noise,
            grammar,
            base_names=BASE_NAMES,
            threshold=THRESHOLD,
            maximum_rounds=MAXIMUM_ROUNDS,
        )
        fold_objects[holdout] = synthesis
        folds[holdout] = {
            "source_domains": sorted(sources),
            "holdout_domain": holdout,
            "synthesis": synthesis.to_dict(),
        }

    synthesis_all = synthesize(
        discovery,
        discovery_table,
        discovery_similarity,
        discovery_noise,
        grammar,
        base_names=BASE_NAMES,
        threshold=THRESHOLD,
        maximum_rounds=MAXIMUM_ROUNDS,
    )

    frozen_programs = {
        "status": "frozen_before_audit_or_real-domain_loading_within_this_run",
        "external_timestamp": None,
        "grammar_manifest_sha256": sha256(OUTPUT / "grammar_manifest.json"),
        "declared_split_manifest_sha256": sha256(OUTPUT / "declared_split_manifest.json"),
        "all_source_domains": synthesis_all.to_dict(),
        "outer_folds": {name: fold.to_dict() for name, fold in fold_objects.items()},
    }
    write_json(OUTPUT / "frozen_programs.json", frozen_programs)
    print("programs frozen; loading audit and real test domains", flush=True)

    independent_audit = {
        "logistic": [
            replace(item, identifier=f"audit-{item.identifier}")
            for item in logistic_domain(n_per_class=120, seed=SEED + 191)
        ],
        "lorenz": [
            replace(item, identifier=f"audit-{item.identifier}")
            for item in lorenz_domain(n_per_class=120, seed=SEED + 192)
        ],
        "kuramoto": [
            replace(item, identifier=f"audit-{item.identifier}")
            for item in kuramoto_domain(n_per_class=120, seed=SEED + 193)
        ],
    }
    sealed = {
        "ecgfivedays": sealed_ucr_domain("ECGFiveDays", archives["ECGFiveDays"], per_class=32, seed=SEED + 41),
        "earthquakes": sealed_ucr_domain("Earthquakes", archives["Earthquakes"], per_class=32, seed=SEED + 42),
    }
    audit_records = [item for records in independent_audit.values() for item in records]
    real_records = [item for records in sealed.values() for item in records]
    postfreeze_records = audit_records + real_records
    postfreeze_table = feature_table(postfreeze_records, "original")
    postfreeze_similarity = feature_table(postfreeze_records, "similarity")
    postfreeze_noise = feature_table(postfreeze_records, "noise")
    table = merge_tables(discovery_table, postfreeze_table)
    similarity_table = merge_tables(discovery_similarity, postfreeze_similarity)
    noise_table = merge_tables(discovery_noise, postfreeze_noise)
    all_records = discovery_records + postfreeze_records

    transfer_rows: list[dict[str, object]] = []
    ablation_grammars = grammar_ablations(grammar)
    for holdout, synthesis in fold_objects.items():
        heldout = independent_audit[holdout]
        reduction = paired_domain_reduction(
            heldout, table, synthesis.base_names, synthesis.selected, threshold=THRESHOLD
        )
        heldout_metrics = evaluate_domain(
            heldout, table, synthesis.base_names, synthesis.selected, threshold=THRESHOLD
        )
        fixed_null = null_p_value(
            heldout,
            table,
            synthesis.base_names,
            synthesis.selected,
            observed=float(reduction["absolute_reduction"]),
        )
        sources = {name: records for name, records in discovery.items() if name != holdout}
        admissible = admissible_programs(
            grammar,
            [item for records in sources.values() for item in records],
            discovery_table,
            discovery_similarity,
            discovery_noise,
        )
        random_null = random_program_null(
            admissible,
            heldout,
            table,
            observed=float(reduction["absolute_reduction"]),
        )
        pipeline_null = full_pipeline_label_null(
            sources,
            heldout,
            discovery_table,
            discovery_similarity,
            discovery_noise,
            table,
            grammar,
            observed=float(reduction["absolute_reduction"]),
        )
        ablations: dict[str, object] = {}
        for ablation_name, ablation_grammar in ablation_grammars.items():
            ablation = synthesize(
                sources,
                discovery_table,
                discovery_similarity,
                discovery_noise,
                ablation_grammar,
                base_names=BASE_NAMES,
                threshold=THRESHOLD,
                maximum_rounds=1,
            )
            ablations[ablation_name] = {
                "selected": [program.to_dict() for program in ablation.selected],
                "audit": paired_domain_reduction(
                    heldout,
                    table,
                    BASE_NAMES,
                    ablation.selected,
                    threshold=THRESHOLD,
                ),
            }
        folds[holdout]["independent_audit"] = {
            **reduction,
            "metrics": heldout_metrics,
            "fixed_program_label_permutation": fixed_null,
            "admissible_random_program": random_null,
            "full_pipeline_source_label_permutation": pipeline_null,
            "ablations": ablations,
        }
        transfer_rows.append({"domain": holdout, "track": "independent-synthetic-audit", **reduction})
        print(f"independent audit {holdout}: {reduction}", flush=True)

    # This interpretation audit occurs only after the source-selected program
    # has been frozen and evaluated.  It cannot affect selection or upgrade the
    # internally held-out result to prospective evidence.
    kuramoto_posthoc = run_posthoc_diagnostic(
        independent_audit["kuramoto"],
        table,
        fold_objects["kuramoto"].base_names,
        fold_objects["kuramoto"].selected,
        threshold=THRESHOLD,
        output_dir=OUTPUT,
        figure_dir=FIGURES,
    )

    sealed_results: dict[str, dict[str, object]] = {}
    for name, records in sealed.items():
        reduction = paired_domain_reduction(
            records, table, synthesis_all.base_names, synthesis_all.selected, threshold=THRESHOLD
        )
        metrics = evaluate_domain(
            records, table, synthesis_all.base_names, synthesis_all.selected, threshold=THRESHOLD
        )
        null = null_p_value(
            records,
            table,
            synthesis_all.base_names,
            synthesis_all.selected,
            observed=float(reduction["absolute_reduction"]),
        )
        sealed_results[name] = {**reduction, **metrics, "permutation_null": null}
        transfer_rows.append({"domain": name, "track": "sealed-real", **reduction})
        print(f"sealed {name}: {reduction}", flush=True)

    selected = synthesis_all.selected
    checks = theorem_checks(all_records, table, similarity_table, selected)
    checks["time_reversal_noncompleteness"] = noncompleteness_check(discovery["logistic"][0])

    noise_errors: dict[str, dict[str, float]] = {}
    for program in selected:
        values = []
        for item in all_records:
            a = program.evaluate(table[item.identifier])
            b = program.evaluate(noise_table[item.identifier])
            values.append(abs(a - b) / (1.0 + abs(a)))
        noise_errors[program.name] = {
            "median": float(np.median(values)),
            "p95": float(np.quantile(values, 0.95)),
            "maximum": float(np.max(values)),
        }

    transfer_frame = pd.DataFrame(transfer_rows)
    transfer_frame.to_csv(OUTPUT / "transfer_metrics.csv", index=False)
    result = {
        "seed": SEED,
        "collision_threshold": THRESHOLD,
        "claim_boundary": "internal frozen holdouts; no prospective registration or population injectivity claim",
        "synthetic_leave_one_domain_out": folds,
        "synthesis_all_synthetic": synthesis_all.to_dict(),
        "sealed_real_domains": sealed_results,
        "post_hoc_diagnostics": {"kuramoto": kuramoto_posthoc},
        "theorem_checks": checks,
        "noise_robustness": noise_errors,
        "legacy_seismic_status": "historical exploratory pilot only; excluded from confirmatory metrics",
    }
    write_json(OUTPUT / "metrics.json", result)

    evaluated_split_manifest = {
        "independent_synthetic_audit": {
            name: {"count": len(records), "identifiers": [item.identifier for item in records]}
            for name, records in independent_audit.items()
        },
        "sealed_real": {
            name: {"split": "TEST only", "count": len(records), "identifiers": [item.identifier for item in records]}
            for name, records in sealed.items()
        },
    }
    write_json(OUTPUT / "evaluated_split_manifest.json", evaluated_split_manifest)
    plot_transfer(transfer_frame)
    plot_programs(result)

    manifest_paths = sorted(
        [path for path in OUTPUT.glob("*") if path.is_file() and path.name != "manifest.sha256"]
        + sorted(FIGURES.glob("*.pdf"))
        + sorted(FIGURES.glob("*.png"))
    )
    manifest_lines = [f"{sha256(path)}  {path.relative_to(ROOT)}" for path in manifest_paths]
    (OUTPUT / "manifest.sha256").write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")
    print(f"wrote {OUTPUT / 'metrics.json'}", flush=True)


if __name__ == "__main__":
    main()
