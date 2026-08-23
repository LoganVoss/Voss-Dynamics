# Paper III claim ledger

This ledger separates mathematical results, deterministic implementation checks, internal held-out evidence, and unresolved claims. The authoritative numerical record is `evidence/outputs/metrics.json`; dataset membership and hashes are in `evidence/outputs/declared_split_manifest.json`, `evidence/outputs/evaluated_split_manifest.json`, and `evidence/outputs/manifest.sha256`.

## Status vocabulary

- **PROVED**: follows analytically under the stated assumptions; computation is only a check.
- **COMPUTED**: deterministic result of the archived implementation and seed `20260822`.
- **HELD-OUT**: measured on records unavailable to the corresponding synthesis step. These are internal frozen holdouts, not prospectively registered experiments.
- **OPEN**: not established by this paper and prohibited as an affirmative conclusion.

## Proved under stated assumptions

### P3-P1 — Exact observation refinement

**Claim.** For any unary observable `f`, appending it as `O_(k+1)(x) = (O_k(x), f(x))` refines the equivalence relation induced by `O_k`. The refinement is strict exactly when `f` is nonconstant on at least one `O_k` fibre.

**Assumptions.** Equality is exact and `f` is a function of one trajectory, not an unanchored pairwise score.

**Executable witness.** `evidence/tests/test_pdis_core.py::test_synthesized_program_reenters_unary_representation`.

### P3-P2 — Conditional finite-sample completion

**Claim.** If a finite sample has `n` elements and its initial observation induces `c_0` classes, and if the candidate class contains a separator for an unresolved nonsingleton fibre at every successful round, then repeated strict refinement reaches sample injectivity in at most `n - c_0` successful additions.

**Boundary.** This is conditional on separator existence and successful search. It says nothing about population injectivity, global observability, or finite-noise identity.

**Executable witness.** `evidence/tests/test_pdis_core.py::test_finite_sample_synthesis_splits_a_collision_class` checks one strict split, not the general proof.

### P3-P3 — Threshold-collision nesting

**Claim.** With the threshold fixed, old coordinates left unchanged, and representation distance equal to the coordinatewise maximum after the fixed monotone squash, appending a coordinate cannot create a new threshold collision.

**Boundary.** Re-estimating coordinate weights, rescaling old coordinates, changing the squash, or changing the threshold can break nesting.

**Executable witness.** `evidence/tests/test_pdis_core.py::test_threshold_collisions_are_nested_after_append`.

### P3-P4 — Quantile-recurrence similarity invariance

**Claim.** Under `x_t -> a Q x_t + b`, where `a > 0` and `Q` is orthogonal, all Euclidean interpoint distances and their fixed empirical quantile threshold scale by `a`; with deterministic tie handling, recurrence adjacency is unchanged.

**Boundary.** The result does not cover anisotropic scaling, arbitrary nonlinear observation maps, resampling, missingness, or general time warps. It also does not make the recurrence spectrum a complete descriptor.

**Executable witness.** `evidence/tests/test_pdis_core.py::test_quantile_recurrence_similarity_invariance`.

### P3-P5 — Exact binomial collision-risk bound

**Claim.** For `K` independent Bernoulli collision witnesses and zero observed collisions, the exact one-sided `1 - alpha` Clopper--Pearson upper bound is `1 - alpha^(1/K)`.

**Boundary.** The result does not license treating all `n(n-1)/2` dependent pairs as independent. In this implementation evaluation uses one disjoint cross-label matching; source selection alone uses the dependent all-pairs statistic descriptively.

**Executable witness.** `evidence/tests/test_pdis_core.py::test_zero_collision_bound_is_exact_rule_of_three_generalization` checks that 0/299 gives an upper 95% bound below 0.0101.

## Computed implementation facts

### P3-C1 — Frozen grammar and selected all-source program

The manifest contains 271 expressions, two base coordinates (`radial_lag1`, `path_tortuosity`), threshold `0.22`, at most three rounds, and nuisance group “translation x positive global scale x orthogonal channel action.” Synthesis on all three discovery domains selected one complexity-3 program:

`log_ratio(increment_cv, radial_permutation_entropy)`.

The selected expression is a nuisance-invariant label separator. It is not claimed to be a dynamical first integral or conservation law.

**Evidence.** `evidence/outputs/grammar_manifest.json`, `evidence/outputs/frozen_programs.json`.

### P3-C2 — Similarity and graph-indexing checks

- Recurrence adjacency agreed exactly in 24/24 tested similarity transforms.
- The all-source selected program had maximum relative similarity error `2.58e-14` on the recorded check.
- Normalized-Laplacian indexing returned `lambda_2(P_4) = 0.5` and `lambda_2(2P_2) = 0`, so algebraic connectivity uses the second-smallest eigenvalue including zero multiplicity.

**Evidence.** `metrics.json -> theorem_checks`.

### P3-C3 — Spectral non-completeness witness

For trajectory `logistic-0-000`, reversing temporal order produced a different oriented history while the recurrence normalized-Laplacian spectra agreed to maximum absolute difference `3.83e-15`. Thus the spectrum does not identify time orientation.

**Evidence.** `metrics.json -> theorem_checks.time_reversal_noncompleteness`.

### P3-C4 — Source-fit diagnostics are not transfer evidence

After fitting on all three discovery domains, collision counts on 28 disjoint cross-label pairs per domain changed as follows: Kuramoto `28 -> 2`, logistic `28 -> 0`, and Lorenz `24 -> 0`.

**Boundary.** These are source-fit diagnostics. They must not be described as held-out or cross-domain results.

**Evidence.** `metrics.json -> synthesis_all_synthetic.source_before/source_after`.

### P3-C5 — Recorded noise sensitivity

For the all-source selected expression under the archived noise perturbation, relative deviation had median `0.0340`, 95th percentile `0.2651`, and maximum `0.4148`.

**Boundary.** This is one prescribed synthetic perturbation distribution. It does not establish general robustness, and its recorded 95th percentile exceeds `0.20` on the broader post-freeze check.

**Evidence.** `metrics.json -> noise_robustness`.

## Internal held-out evidence

### P3-H1 — Leakage-controlled synthetic leave-one-domain-out transfer

Each fold synthesized on two discovery domains and was evaluated on 240 independently seeded audit trajectories from the omitted domain, yielding 120 disjoint cross-label witnesses:

| Held-out domain | Frozen source-selected program | Collisions before -> after | Absolute reduction | Paired bootstrap 95% interval | Final risk upper 95% | Orientation-free AUC |
|---|---|---:|---:|---:|---:|---:|
| Logistic | `log_ratio(increment_cv, radial_permutation_entropy)` | 119 -> 0 / 120 | 0.9917 | [0.9750, 1.0000] | 0.0247 | 1.0000 |
| Lorenz | `log_ratio(covariance_entropy, h1_max_persistence)` | 110 -> 56 / 120 | 0.4500 | [0.3583, 0.5333] | 0.5457 | 0.8278 |
| Kuramoto | `log_ratio(increment_cv, radial_permutation_entropy)` | 120 -> 10 / 120 | 0.9167 | [0.8667, 0.9667] | 0.1373 | 0.9976 |

**Permitted conclusion.** Under the frozen grammar, threshold, generators, labels, and nuisance group, source-selected programs reduced collision risk on independent trajectories from each of the three named held-out synthetic systems.

**Prohibited conclusion.** This does not establish transfer to arbitrary dynamics or a population of physical domains.

**Evidence.** `metrics.json -> synthetic_leave_one_domain_out`, `transfer_metrics.csv`, and the two split manifests.

### P3-H2 — Synthetic null tests

- Fixed-program held-out label permutations used 999 repetitions and gave `p = 0.001` for logistic, `p = 0.027` for Lorenz, and `p = 0.001` for Kuramoto.
- Full source-pipeline label permutations used 99 repetitions. No permuted fold selected any program; the reported `p = 0.01` is the minimum attainable corrected value and should be described as resolution-limited.
- Against one independently sampled program from the nuisance-admissible pool per repetition, the winner exceeded the random-program comparator for Kuramoto (`p = 0.034`) but not for logistic (`p = 0.119`) or Lorenz (`p = 0.166`). This is not an equal-compute rerun of source selection.

**Boundary.** The random-program results do not support a uniform claim that guided synthesis beats random search in every domain.

**Evidence.** `metrics.json -> synthetic_leave_one_domain_out.*.independent_audit`.

### P3-H3 — Ablations

Primitive-only, recurrence-only, and topology-only searches did not reproduce the full program's reduction uniformly. The strongest topology-only reductions were `0.2917` on Kuramoto and `0.3833` on Lorenz; recurrence-only achieved `0.1333` on logistic and zero on the other two held-out systems; primitive-only achieved `0.1333` on Lorenz and zero on logistic and Kuramoto.

**Permitted conclusion.** No restricted family was sufficient across all three synthetic holdouts under this protocol.

**Boundary.** The experiment does not prove that recurrence or persistent-homology primitives are individually necessary, nor does it compare against full external symbolic-regression or learned-representation baselines.

**Evidence.** `metrics.json -> synthetic_leave_one_domain_out.*.independent_audit.ablations`.

### P3-H4 — Frozen transfer to two UCR archive test sets

The all-source synthetic program and threshold were frozen before reading balanced subsets of official UCR TEST members, 32 cases per class:

| UCR domain | Collisions before -> after | Absolute reduction | Paired bootstrap 95% interval | Label-permutation p | Program AUC |
|---|---:|---:|---:|---:|---:|
| Earthquakes | 31 -> 23 / 32 | 0.2500 | [0.1250, 0.4063] | 0.41 | 0.7314 |
| ECGFiveDays | 32 -> 32 / 32 | 0.0000 | [0, 0] | 1.00 | 0.5020 |

**Permitted conclusion.** A frozen synthetic-selected expression produced a descriptive reduction on the sampled Earthquakes TEST cases and none on ECGFiveDays.

**Prohibited conclusion.** The real-data results do not establish statistically significant transfer: the Earthquakes fixed-program label-permutation test was nonsignificant, and ECGFiveDays was a clear null. UCR series may also share acquisition sources, so treating matched cases as independent population draws is an additional modeling assumption; the reported binomial bounds and pair bootstrap are not cluster-robust physical-domain inference.

**Evidence.** `metrics.json -> sealed_real_domains`, `transfer_metrics.csv`, declared archive SHA-256 values, and `test_real_loader_reads_only_declared_test_member`.

### P3-H5 — Reproducibility boundary

The declared archive hashes, declared/evaluated record identifiers, grammar, frozen programs, metrics, figures, and split manifests are checksummed in `evidence/outputs/manifest.sha256`. The third-party archives are acquired separately through the checksum-enforcing downloader rather than redistributed. The claim boundary recorded by the pipeline is: “internal frozen holdouts; no prospective registration or population injectivity claim.” The legacy Yellowstone seismic analysis is explicitly excluded from confirmatory metrics and retained only as a historical exploratory pilot.

## Open or explicitly unsupported

### P3-O1 — Population injectivity and global observability

Not established. Finite-sample collision reduction cannot prove state identity, global observability, or injectivity outside the sampled task-relevant distinction relation.

### P3-O2 — Universal or arbitrary-domain transfer

Not established. Three synthetic system families and two small archive evaluations support claims only about those named generators, records, grammar, nuisance group, threshold, and binary distinction tasks.

### P3-O3 — Real-domain efficacy

Not established. One real archive domain was null and the other did not reject its fixed-program label-permutation null. Independent event-, subject-, or machine-grouped external replication remains required.

### P3-O4 — Conserved-law discovery

Not established. The selected programs are nuisance-invariant discriminative coordinates; no recovered expression was shown constant along an orbit or matched to an independently known first integral.

### P3-O5 — Causal or mechanistic interpretation

Not established. Label separation, recurrence summaries, persistence features, or Koopman-style coordinates do not by themselves identify a physical mechanism or causal variable.

### P3-O6 — Complete descriptor or canonical form

Not established. Recurrence spectra and persistence summaries are stable/useful descriptors but are not complete invariants; the time-reversal witness demonstrates one concrete ambiguity.

### P3-O7 — Comprehensive baseline superiority

Not established. The archived evidence contains restricted-family ablations and random-program nulls, not equal-budget implementations of PySR/SRBench methods, SINDy, DMD/EDMD/HAVOK, catch22, ROCKET, TS2Vec, InceptionTime, or HIVE-COTE 2.0.

### P3-O8 — Prospective confirmation

Not established. A prospectively timestamped protocol, larger permutation budgets, cluster-aware real-domain sampling, and independent replication are future work.
