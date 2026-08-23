# Paper III Development Audit

## Purpose

This document preserves the chronology of the Paper III evidence pipeline, including a failed null test and the corrective rerun. It is not an external registration record. The frozen-program artifact has `external_timestamp: null`; its exact status is “frozen before audit or real-domain loading within this run.”

The current claim boundary is therefore:

> Internal frozen holdouts; no prospective registration or population injectivity claim.

The legacy seismic pilot is excluded throughout.

## Stage 1: replacement of the recovered legacy engine

Forensic testing of Voss-Codex 0.15.0 identified non-reinserted features, unequal before/after thresholds, mixed unary and relational function types, incorrect normalized-Laplacian indexing, a one-perturbation stability gate, auxiliary-domain leakage, and noncausal seismic preprocessing. The exact findings are recorded in `evidence/LEGACY_ENGINE_AUDIT.md`.

Paper III was therefore implemented as a new, typed unary program system rather than as a numerical continuation of the legacy engine. The replacement introduced:

- similarity-normalized trajectory coordinates;
- a typed unary grammar;
- actual program-column reinsertion;
- one fixed collision threshold before and after insertion;
- source-only candidate construction;
- disjoint, independently seeded outer synthetic audit trajectories;
- TEST-only balanced subsets of two real UCR archives;
- fixed, source-independent coordinate squashing;
- invariance and perturbation checks;
- explicit split, grammar, program, output, and archive manifests.

## Stage 2: the first full-pipeline label-null failed

The first complete statistical audit was deliberately treated as a development run. Although its observed synthetic reductions were large, the full-pipeline source-label permutation exposed an invalid selection behavior: every one of 99 null repetitions selected at least one program in every leave-one-domain-out fold.

These values are superseded and are not the contents of the current `evidence/outputs/` directory. They are retained here because the failed diagnostic materially changed the algorithm.

| Outer synthetic domain | Observed reduction in failed run | Bootstrap 95% interval | Full-pipeline null selected any program | Full-pipeline p-value | Admissible random-program comparator p-value |
|---|---:|---:|---:|---:|---:|
| Kuramoto | 0.9250 | [0.8750, 0.9667] | 99/99 | 0.22 | 0.039 |
| Logistic | 0.9917 | [0.9750, 1.0000] | 99/99 | 0.38 | 0.124 |
| Lorenz | 0.4000 | [0.3167, 0.4917] | 99/99 | 1.00 | 0.183 |

The problem was not repaired by reinterpretation. Candidate ranking admitted programs that passed only similarity and noise checks; under permuted labels, the discoverer could still select a candidate on every run. The full-pipeline null therefore showed that the selection rule itself lacked a meaningful effect gate.

## Stage 3: pre-ranking candidate gates were added

The `synthesize()` selection rule was changed so that a candidate is rejected before ranking unless it satisfies all three source-inner criteria:

| Gate | Required value |
|---|---:|
| Minimum collision-risk reduction on every source inner-audit partition | 0.15 |
| Minimum orientation-free AUC on every source inner-audit partition | 0.75 |
| Minimum orientation-free AUC on every source selection partition | 0.70 |

These are source-domain selection gates. They do not inspect the outer synthetic audit trajectories or the real archive cases. The invariance tolerance (`5e-6`) and source noise tolerance (`0.20` at the 95th percentile) remain separate admissibility requirements.

The thresholds are encoded in both `src/pdis/synthesis.py` and `evidence/outputs/grammar_manifest.json`.

## Stage 4: independent final rerun after the fix

After the gates were added, the complete evidence pipeline was rerun from the declared seed `20260822`. Within the run, programs were selected and written to `frozen_programs.json` before the independent synthetic audit trajectories or real TEST cases were loaded and featurized.

The independent synthetic audit used 120 trajectories per class for each outer domain, with seeds declared before generation in `declared_split_manifest.json`. None of those audit identifiers appears in the discovery partitions.

### Final synthetic results

| Outer domain | Frozen program | Reduction | Bootstrap 95% interval | Program AUC | Fixed-program label p-value | Admissible random-program comparator p-value | Full-pipeline null selected any program | Full-pipeline p-value |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Kuramoto | `log_ratio(increment_cv,radial_permutation_entropy)` | 0.9167 | [0.8667, 0.9667] | 0.9976 | 0.001 | 0.034 | 0/99 | 0.01 |
| Logistic | `log_ratio(increment_cv,radial_permutation_entropy)` | 0.9917 | [0.9750, 1.0000] | 1.0000 | 0.001 | 0.119 | 0/99 | 0.01 |
| Lorenz | `log_ratio(covariance_entropy,h1_max_persistence)` | 0.4500 | [0.3583, 0.5333] | 0.8278 | 0.027 | 0.166 | 0/99 | 0.01 |

The repaired pipeline selected no program in any of the 297 source-label-permutation runs. This resolves the specific failure observed in Stage 2. The 99-repetition design limits the minimum attainable full-pipeline p-value to `0.01`.

The admissible random-program comparator is intentionally named narrowly. It samples one program that satisfies the source invariance and noise conditions and applies that program to the outer audit; it is not a rerun of the complete source-selection pipeline. The selected program exceeds this comparator at `p < 0.05` for Kuramoto, but not for Logistic or Lorenz.

### Final sealed-real results

The all-synthetic program was frozen before loading the real TEST cases. Each archive contributes a reproducibly sampled, balanced subset of 32 cases per class.

| Real TEST domain | Reduction | Bootstrap 95% interval | Program AUC | Label-permutation p-value |
|---|---:|---:|---:|---:|
| Earthquakes | 0.2500 | [0.1250, 0.4063] | 0.7314 | 0.41 |
| ECGFiveDays | 0.0000 | [0.0000, 0.0000] | 0.5020 | 1.00 |

These real-domain results do not confirm transfer. Earthquakes has a positive point estimate but does not beat its label-permutation null, and ECGFiveDays shows no reduction. They are reported as mixed or null transfer evidence, not as validation across real domains.

### Verification state

The completed final rerun produced the following SHA-256 values:

| Artifact | SHA-256 |
|---|---|
| `evidence/outputs/metrics.json` | `331516062988c4043c25509566ae9c051318a59b7fa2149d34a68d0780ed630c` |
| `evidence/outputs/frozen_programs.json` | `eaa7bb650838115831a8ef352720077243006929f0aa90e1cbd106ffcaa3795d` |
| `evidence/outputs/grammar_manifest.json` | `f569a8e557b6aae147c1323cdef4550635f394739428ca0b3173932ea516b86c` |
| `evidence/outputs/manifest.sha256` | `4808263fef263bd2c4ba55fc51e034c0059cf91a8eac97b4ab211d88a2fe42ee` |

Every entry in `manifest.sha256` verified after the final run, including the split manifest that records both raw-archive hashes, every JSON/CSV output, and all generated PDF/PNG figures. The core test suite passed 8/8 tests. The third-party archives themselves are acquired through the checksum-enforcing downloader and are not redistributed in the public repository.

Reproduction and verification commands from the Paper III directory are:

```sh
uv sync --frozen
uv run python -m evidence.run_all
shasum -a 256 -c evidence/outputs/manifest.sha256
uv run pytest
```

On systems with GNU Coreutils, `sha256sum -c evidence/outputs/manifest.sha256` is equivalent.

## Remaining claim limits

- The successful synthetic result is an internally frozen, seeded benchmark result, not an externally registered prospective experiment.
- The outer unit is one independently seeded trajectory; the collision statistic uses one disjoint cross-label pairing per domain.
- The full-pipeline null has 99 repetitions.
- Superiority over the admissible random-program comparator is established for only one of three synthetic outer domains.
- The source-domain noise gate is a selection-time condition. The all-record diagnostic for the all-synthetic program has a 95th-percentile relative noise error of approximately `0.265`, above the source selection threshold, and should not be described as universal robustness.
- Neither real archive supplies confirmatory transfer evidence.
- No empirical result establishes injectivity over a population or completeness of the representation.
- The legacy seismic pilot is historical exploratory work and is excluded from all confirmatory totals.
