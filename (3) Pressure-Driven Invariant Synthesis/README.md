# Pressure-Driven Invariant Synthesis

> **Claim boundary:** internal frozen holdouts; no prospective registration or population injectivity claim

This is Paper III of Voss Dynamics. Papers I and II study how distinctions are preserved, collapsed, and tested for predictive significance under fixed representations. Paper III makes the representation dynamic: observed collisions create pressure for a typed unary program, the program must pass source-domain invariance, noise, separation, and AUC gates, and an accepted program is appended to the representation before evaluation on untouched records.

The package contains conditional finite-sample theorems, an executable synthesis engine, a complete evidence runner, exact split and artifact manifests, and mixed/null real-domain transfer results. It does not claim recovery of a physical law, general injectivity, or prospectively registered validation.

## Package map

| Component | Location |
|---|---|
| Manuscript source | [`thesis/main.tex`](thesis/main.tex) |
| Bibliography | [`thesis/references.bib`](thesis/references.bib) |
| Canonical representation | [`src/pdis/canonical.py`](src/pdis/canonical.py) |
| Typed program grammar | [`src/pdis/programs.py`](src/pdis/programs.py) |
| Synthesis and reinsertion | [`src/pdis/synthesis.py`](src/pdis/synthesis.py) |
| Collision statistics | [`src/pdis/statistics.py`](src/pdis/statistics.py) |
| Complete evidence runner | [`evidence/run_all.py`](evidence/run_all.py) |
| Machine-readable final metrics | [`evidence/outputs/metrics.json`](evidence/outputs/metrics.json) |
| Post-hoc Kuramoto interpretation | [`evidence/outputs/kuramoto_posthoc.json`](evidence/outputs/kuramoto_posthoc.json) · [record table](evidence/outputs/kuramoto_posthoc_records.csv) · [channel-subset table](evidence/outputs/kuramoto_channel_subsets.csv) |
| Frozen programs and source audits | [`evidence/outputs/frozen_programs.json`](evidence/outputs/frozen_programs.json) |
| Declared and evaluated splits | [`evidence/outputs/declared_split_manifest.json`](evidence/outputs/declared_split_manifest.json), [`evidence/outputs/evaluated_split_manifest.json`](evidence/outputs/evaluated_split_manifest.json) |
| Development chronology, including the failed null | [`evidence/DEVELOPMENT_AUDIT.md`](evidence/DEVELOPMENT_AUDIT.md) |
| Recovered legacy-engine audit | [`evidence/LEGACY_ENGINE_AUDIT.md`](evidence/LEGACY_ENGINE_AUDIT.md) |
| Detailed reproduction protocol | [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) |

## What the implementation establishes

The mathematical layer proves finite-sample statements under explicit assumptions: append-only coordinates refine representation fibres; a successful separating grammar terminates on a finite sample; and max-metric collision witnesses cannot be introduced by appending a fixed coordinate while retaining the metric, squash, and threshold.

The computational layer checks the declared similarity action, correct normalized-Laplacian zero multiplicity, actual feature reinsertion, fixed-threshold nesting, deterministic generation, and TEST-only archive loading. These are implementation checks, not population-level identification results.

## Final evidence

The final seed is `20260822`, and the fixed collision threshold is `0.22`. Each outer synthetic audit uses 120 independently generated trajectories per class and one disjoint cross-label matching, giving 120 paired collision trials.

### Independent synthetic audits

| Outer domain | Frozen program | Absolute collision-risk reduction | Bootstrap 95% interval | Program AUC | Full-pipeline source-label null | Admissible random-program comparator p-value |
|---|---|---:|---:|---:|---:|---:|
| Kuramoto | `log_ratio(increment_cv,radial_permutation_entropy)` | 0.9167 | [0.8667, 0.9667] | 0.9976 | 0/99 selections; p = 0.01 | 0.034 |
| Logistic | `log_ratio(increment_cv,radial_permutation_entropy)` | 0.9917 | [0.9750, 1.0000] | 1.0000 | 0/99 selections; p = 0.01 | 0.119 |
| Lorenz | `log_ratio(covariance_entropy,h1_max_persistence)` | 0.4500 | [0.3583, 0.5333] | 0.8278 | 0/99 selections; p = 0.01 | 0.166 |

The repaired selection pipeline passes the full-pipeline label-null diagnostic. Superiority to the narrower admissible random-program comparator reaches `p < 0.05` for Kuramoto only.

### Post-hoc Kuramoto interpretation

Because Kuramoto is the only fold that beats the admissible random-program comparator, the exact audit was regenerated after the result was known while retaining the simulator’s latent phases. The observed sine-channel arrays match the released audit bit-for-bit. This is an explicitly **computed, post-hoc diagnostic**; it does not upgrade the internally held-out claim.

| Diagnostic | Result |
|---|---:|
| Frozen scalar label AUC | 0.997569 |
| Latent mean first- / second-harmonic order label AUC | 0.998889 / 0.998889 |
| Spearman: frozen scalar vs mean first- / second-harmonic order | 0.860826 / 0.906991 |
| Stratified record-bootstrap 95% intervals | [0.832641, 0.885370] / [0.884555, 0.926683] |
| Weak-class Spearman: scalar vs mean second-harmonic order | 0.944031 [0.908270, 0.963261] |
| Mean subset AUC, 1 / 2 / 3 / 4 / 5 channels | 0.522 / 0.950 / 0.992 / 0.998 / 0.998 |
| All-pairs post-repair collision risk | 0.09090 |

The key interpretation is relational: a single sine channel is null on average, while recomputing subset-specific canonicalization and the full scalar on two or more channels rapidly restores separation. For the noiseless sensor vector, the exact identity ‖sin θ‖² = N/2 · (1 − Re Z₂) connects radial trajectory structure to generalized phase order. The released array adds Gaussian measurement noise and is then centered and globally scaled, so the identity does not hold exactly for the computed feature. Within the weak-coupling class, the ratio correlates more strongly with mean |Z₂| than either component alone (0.944 versus 0.888 for increment CV and −0.864 for radial permutation entropy); paired record-bootstrap intervals for both improvements exclude zero. This is a post-hoc within-regime pattern, not confirmation. Across the coarse, widely separated labels, increment CV alone has slightly higher AUC than the ratio. The frozen scalar is therefore an **order-sensitive multichannel surrogate inside this generator**—not the Kuramoto order parameter itself, a conserved law, a critical-coupling estimator, a universal order coordinate, or a uniquely necessary separator.

### Frozen real TEST audits

| Domain | Absolute collision-risk reduction | Bootstrap 95% interval | Program AUC | Label-permutation p-value |
|---|---:|---:|---:|---:|
| Earthquakes | 0.2500 | [0.1250, 0.4063] | 0.7314 | 0.41 |
| ECGFiveDays | 0.0000 | [0.0000, 0.0000] | 0.5020 | 1.00 |

The real-domain results do not confirm transfer. Earthquakes has a positive point estimate that does not exceed its permutation null; ECGFiveDays shows no reduction.

## Freeze semantics

“Frozen” has one precise meaning in this package: within the final evidence run, the grammar, source splits, selected programs, and source diagnostics are serialized before independent synthetic audit trajectories and real TEST cases are loaded and featurized. `frozen_programs.json` has `external_timestamp: null`. This is an internal execution-order boundary, not an external registration or prospective commitment.

The [development audit](evidence/DEVELOPMENT_AUDIT.md) records the first failed full-pipeline label-null, the pre-ranking gates added in response, and the independent final rerun. The failed run is not represented as final evidence.

## Quick reproduction

Install [uv](https://docs.astral.sh/uv/), then run from this directory:

```sh
sh evidence/data/download_ucr.sh
uv sync --frozen
uv run pytest -q
uv run python evidence/run_all.py
shasum -a 256 -c evidence/outputs/manifest.sha256
```

Use `sha256sum -c evidence/outputs/manifest.sha256` instead of `shasum` on systems with GNU Coreutils. The evidence run regenerates JSON, CSV, and the six figure files covered by the manifest. See [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) for environment details, expected values, integrity hashes, and failure handling.

## Data and license provenance

The real-data inputs are the official UCR/UEA Time Series Classification archives `ECGFiveDays.zip` and `Earthquakes.zip`. Their official URLs, byte sizes, and SHA-256 hashes are recorded in [`evidence/data/README.md`](evidence/data/README.md). Only each archive's `*_TEST.ts` member is read. The public repository ships the checksum-enforcing downloader and exact provenance records, not the third-party archives themselves. The archives remain governed by the source repository and original contributors' terms.

Python dependencies remain governed by their respective upstream licenses. Under the repository-level [`LICENSE`](../LICENSE), original software is MIT-licensed and original manuscript text and documentation are CC BY 4.0; third-party materials remain excluded.

## Legacy exclusion

Voss-Codex 0.15.0 and the earlier seismic benchmark were forensic inputs to the redesign, not empirical inputs to the final results. The recovered engine did not reinsert synthesized features into evaluated vectors, mixed unary and relational APIs, misindexed algebraic connectivity, and used leaking seismic evaluation paths. The seismic pilot is classified as historical exploratory work and is excluded from all Paper III confirmatory totals. Exact findings and source identities are in [`evidence/LEGACY_ENGINE_AUDIT.md`](evidence/LEGACY_ENGINE_AUDIT.md).

## Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). Cite the preferred manuscript citation together with the exact repository revision used for computation.
