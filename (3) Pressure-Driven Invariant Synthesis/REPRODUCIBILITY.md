# Reproducibility

> **Claim boundary:** internal frozen holdouts; no prospective registration or population injectivity claim

This protocol reproduces the Paper III tests, seeded evidence tables, statistical diagnostics, manifests, and figures. It does not reproduce or incorporate the excluded legacy seismic pilot.

## Release environment

The final committed evidence was generated with:

| Component | Release environment |
|---|---|
| Operating system | macOS 26.5.2, Apple silicon |
| Python | 3.14.6 |
| uv | 0.12.3 |
| Tectonic, for the manuscript | 0.17.0 |
| Evidence seed | `20260822` |
| Collision threshold | `0.22` |

`pyproject.toml` declares Python 3.11 or newer. `uv.lock` is the authoritative numerical dependency lock and contains platform/Python markers. The exact release lock has SHA-256:

```text
ef100a23d0696711e8234cc2975f633ac5830606d67a2d3caf9322442150625e  uv.lock
```

For the closest replication, use Python 3.14 on a 64-bit platform and do not update the lock before running the audit.

## Inputs

The synthetic discovery and audit records are generated deterministically by `src/pdis/datasets.py`. The independent audit seeds are declared in `evidence/outputs/declared_split_manifest.json`.

The only real-data inputs are the following externally acquired archives; they are not redistributed in the public repository:

| File | SHA-256 |
|---|---|
| `evidence/data/raw/ECGFiveDays.zip` | `11457a2d590711598eac1a0ab5c58d43c5e3c5c2c86521809d69f7c0b6b3edd1` |
| `evidence/data/raw/Earthquakes.zip` | `927cfb732988055850a74efa169dfd633bc7f578095c35319bebb83453501bf1` |

If either archive is absent, retrieve and verify it with:

```sh
sh evidence/data/download_ucr.sh
```

The download helper uses the official UCR/UEA Time Series Classification URLs, verifies SHA-256 before installing a file, and refuses to overwrite a mismatched existing archive. See `evidence/data/README.md` for source links and provenance. The runner reads only `ECGFiveDays_TEST.ts` and `Earthquakes_TEST.ts`, then draws a deterministic balanced subset of 32 cases per class.

## One-command sequence

Run the following from the Paper III directory:

```sh
uv sync --frozen
uv run pytest -q
uv run python evidence/run_all.py
shasum -a 256 -c evidence/outputs/manifest.sha256
```

On GNU/Linux, use:

```sh
sha256sum -c evidence/outputs/manifest.sha256
```

Expected test result:

```text
..........                                                               [100%]
```

The evidence runner is intentionally a write-producing command. It regenerates:

- `evidence/outputs/*.json`;
- `evidence/outputs/*.csv`;
- `thesis/figures/frozen_transfer.{pdf,png}`;
- `thesis/figures/program_auc.{pdf,png}`;
- `thesis/figures/kuramoto_posthoc.{pdf,png}`;
- `evidence/outputs/manifest.sha256`.

Do not hand-edit a generated artifact to make a checksum pass. A missing or mismatched manifest entry means the evidence bundle must be regenerated or the release packaging corrected before citation.

## Optional manuscript build

The manuscript source is `thesis/main.tex`, with bibliography `thesis/references.bib`. Build it separately with:

```sh
tectonic -X compile thesis/main.tex --outdir thesis/build
```

PDF creation metadata can vary by backend or version. The evidence manifest covers the committed evidence figures, not a bitwise-reproducible manuscript PDF.

## Freeze and loading order

The final runner enforces this chronology:

1. Generate the three synthetic discovery domains, 28 trajectories per class.
2. Compute discovery, similarity-transform, and noise-transform feature tables.
3. Write the grammar and declared split manifests.
4. Select each leave-one-domain-out program and the all-synthetic program from source records only.
5. Serialize `evidence/outputs/frozen_programs.json`.
6. Generate the independent synthetic audit records with different seeds.
7. Load and featurize the real UCR TEST cases.
8. Evaluate fixed programs, nulls, ablations, theorem checks, and primary figures.
9. Only after the Kuramoto result is fixed, regenerate that exact audit with latent phases and write the explicitly post-hoc interpretation diagnostic.
10. Write the evaluated split manifest and release-integrity manifest.

The freeze is an in-process ordering guarantee. `frozen_programs.json` records `external_timestamp: null`; it is not evidence of an external registration or prospective commitment.

## Expected final results

### Independent synthetic audits

| Outer domain | Reduction | Bootstrap 95% interval | AUC | Fixed-program label p-value | Admissible random-program comparator p-value | Full-pipeline null |
|---|---:|---:|---:|---:|---:|---:|
| Kuramoto | 0.9167 | [0.8667, 0.9667] | 0.9976 | 0.001 | 0.034 | 0/99 selected; p = 0.01 |
| Logistic | 0.9917 | [0.9750, 1.0000] | 1.0000 | 0.001 | 0.119 | 0/99 selected; p = 0.01 |
| Lorenz | 0.4500 | [0.3583, 0.5333] | 0.8278 | 0.027 | 0.166 | 0/99 selected; p = 0.01 |

The first development version failed the full-pipeline source-label null by selecting a program in 99/99 permutations for every fold. Candidate effect/AUC gates were added before ranking, and the entire audit was rerun. The superseded run and repair are documented in `evidence/DEVELOPMENT_AUDIT.md`.

### Frozen real TEST audits

| Domain | Reduction | Bootstrap 95% interval | AUC | Label-permutation p-value |
|---|---:|---:|---:|---:|
| Earthquakes | 0.2500 | [0.1250, 0.4063] | 0.7314 | 0.41 |
| ECGFiveDays | 0.0000 | [0.0000, 0.0000] | 0.5020 | 1.00 |

These values are mixed/null evidence. They do not confirm transfer to real domains.

### Post-hoc Kuramoto interpretation

The runner then reproduces a diagnostic that was designed after the Kuramoto result was known. The observed audit arrays must regenerate with maximum absolute discrepancy `0.0`. Expected values are:

| Diagnostic | Expected value |
|---|---:|
| Frozen scalar label AUC | 0.997569 |
| Manual scalar vs frozen `Program.evaluate` maximum error | 0.0 |
| Latent mean first- / second-harmonic order label AUC | 0.998889 / 0.998889 |
| Spearman correlation with mean first- / second-harmonic order | 0.860826 / 0.906991 |
| Stratified record-bootstrap 95% intervals | [0.832641, 0.885370] / [0.884555, 0.926683] |
| Weak-class correlation with mean second-harmonic order | 0.944031 [0.908270, 0.963261] |
| Mean channel-subset AUC, 1 / 2 / 3 / 4 / 5 channels | 0.5216 / 0.9499 / 0.9915 / 0.9981 / 0.9976 |
| All-pairs post-repair collision risk | 0.090903 |

For the first-harmonic interval, seed `9922` independently resamples 120 records with replacement inside each label stratum for 10,000 repetitions. The corresponding second-harmonic pooled and weak-class intervals each reset an independent generator to seed `9923`; the weak interval resamples only the 120 weak-coupling records. A dependent-correlation bootstrap independently resets to seed `9924`, resamples the same weak records for all three correlations in each of 10,000 repetitions, and gives scalar-minus-increment and scalar-minus-oriented-entropy differences `0.056469 [0.028756, 0.094744]` and `0.079945 [0.042654, 0.127708]`. Pairing sensitivity uses seed `555` and 20,000 random perfect cross-label matchings; its reported range is descriptive, not a confidence interval. All seeds, methods, and unrounded values are serialized in `kuramoto_posthoc.json`.

Every channel-count mean averages all `C(5,m)` subsets. For each subset the runner starts from the corresponding raw observed columns, reruns centering and global RMS canonicalization, then recomputes both primitives and the frozen ratio. This output is `COMPUTED_POST_HOC_NO_CLAIM_UPGRADE`: an interpretation audit on the same synthetic records, not an additional held-out experiment.

### Deterministic theorem and implementation checks

The final metrics record:

- exact similarity adjacency in 24/24 trials;
- normalized-Laplacian lambda-2 equal to `0.5` for the four-vertex path;
- lambda-2 equal to `0.0` for two disconnected two-vertex paths;
- maximum selected-program similarity-transform relative error approximately `2.58e-14`;
- an explicit time-reversal noncompleteness witness.

## Release-integrity anchors

The final evidence artifacts before repository packaging have these SHA-256 anchors:

| Artifact | SHA-256 |
|---|---|
| `evidence/outputs/metrics.json` | `9408907d336c189be8e027b9e24d420aaac6876255a447fc5d1a2f0cfd7525b8` |
| `evidence/outputs/frozen_programs.json` | `eaa7bb650838115831a8ef352720077243006929f0aa90e1cbd106ffcaa3795d` |
| `evidence/outputs/grammar_manifest.json` | `f569a8e557b6aae147c1323cdef4550635f394739428ca0b3173932ea516b86c` |
| `evidence/outputs/manifest.sha256` | `2a13383bcfa765acf9386355db6ef63fa1487b1908e633f03b2522ee121b91e7` |

The manifest is an integrity record for the generated bundle, not an external timestamp. Repository metadata files are outside that generated-artifact manifest.

## Interpretation limits

- The statistical audits are internally frozen and seeded, not externally registered.
- The synthetic result does not establish population injectivity or recovery of physical laws.
- The admissible random-program comparator is not a complete rerun of source selection; only the source-label permutation is the full-pipeline null.
- Only Kuramoto exceeds the admissible random-program comparator at `p < 0.05`.
- The all-record 95th-percentile relative noise diagnostic is approximately `0.265`, so the source-domain noise gate is not evidence of universal noise robustness.
- Neither real archive provides confirmatory transfer evidence.
- The Kuramoto interpretation diagnostic was designed after the result was seen. Its permitted reading is an order-sensitive multichannel surrogate inside the named generator, with a descriptive second-harmonic association inside the weak class—not a new conserved quantity, universal order parameter, or claim that composition is required for the coarse target separation.

## License and legacy boundary

The UCR/UEA archives and Python dependencies retain their upstream licenses and terms. They are not relicensed by this package. Under the repository-level `LICENSE`, original software is MIT-licensed and original manuscript text and documentation are CC BY 4.0; third-party materials remain excluded.

The recovered Voss-Codex 0.15.0 engine and the historical seismic pilot are excluded from this reproduction. Their defects, unavailable raw inputs, and leaking evaluation paths are documented in `evidence/LEGACY_ENGINE_AUDIT.md`. No legacy collision count, Watchtower result, stability score, transfer proxy, or seismic detection rate enters `metrics.json`.
