# Legacy Engine Forensic Audit

## Status and scope

This document records the forensic examination of the earlier Voss-Codex engine, version 0.15.0. It is a development-history artifact, not evidence produced by the current Paper III implementation.

The legacy code and its seismic pilot are excluded from the confirmatory metrics in `evidence/outputs/metrics.json`. No legacy performance number is used to support the empirical claims of Paper III.

## Recovered source identity

The working copy of the old package had been removed, but its complete source was recovered from the local Git repository `Voss Attractor System` at commit:

```text
4132c63910f19727738122faed6bf0c7a38d2592
Initial snapshot: Voss Attractor System (V16 / H29)
2026-08-05T15:59:17-07:00
```

The package root at that commit is `Voss-Codex/`. Relevant Git objects are:

| Recovered path | Git blob |
|---|---|
| `Voss-Codex/src/voss_codex/canonical.py` | `e01c956e3f54b46b308df8e9f5212135bb41c2cd` |
| `Voss-Codex/src/voss_codex/engine.py` | `780388573b303644687ac8544d90bbca09322864` |
| `Voss-Codex/src/voss_codex/invariants.py` | `d91cf9d8e998c6e4cde49ecf36822cc750b9b9cd` |
| `Voss-Codex/src/voss_codex/persistence.py` | `4e39103ea29e3a3bbe8a3ebd44ffc1818e48ab12` |
| `Voss-Codex/src/voss_codex/provenance.py` | `dc85c36696b5e1df2d739914367a77fd9a387e1f` |
| `Voss-Codex/src/voss_codex/scoring.py` | `328b9eda64821f29bc18ac7bf23c35babe1253f6` |
| `Voss-Codex/src/voss_codex/synthesis.py` | `c7dd30ae20628eb849d3ba465a066c93462cd510` |
| `Voss-Codex/src/voss_codex/transfer.py` | `79fd9699141c5160ee32e5cd000294285e053de6` |
| `Voss-Codex/src/voss_codex/validation.py` | `a2a044009a2c15ed750b5b0be73c68835d776a84` |
| `Voss-Codex/src/voss_codex/watchtower.py` | `4e412079bd5cb399f2f177f9c94d50132e5f7160` |
| `Voss-Codex/benchmarks/seismic_prospective.py` | `a46d00605b6ca994b0f17fb9a0619740220c40e1` |
| `Voss-Codex/benchmarks/fast_yellowstone_benchmark.py` | `6f8bf0963b48f384b6bedb3a3715dd423d1d3de8` |

A deterministic `git archive` of `Voss-Codex/` at that commit has SHA-256:

```text
f58a8d4867265d96326db12448e61c4291a5b4355e664e8f8718a2b9c0bd51cb
```

The forensic checks below were performed against a temporary extraction of those exact Git objects. The recovered repository was not edited.

## Confirmed defects

### 1. Synthesized features were not reinserted into the evaluated representation

`engine.py` constructs the trajectory vectors once from the base invariant functions. It later adds synthesized callables to `inv_funcs`, but the final collision computation reuses the original vectors. The `vectors` argument supplied to the legacy synthesis routine is itself unused.

Consequently, the feature namespace can report a synthesized observable that never contributes a column to the representation whose collisions are counted.

The apparent before/after improvement could also be generated solely by unequal thresholds. In an isolated reproduction with two one-dimensional vectors at distance `0.08`:

| Configuration | Initial collisions | Final collisions | Synthesized features |
|---|---:|---:|---:|
| Legacy defaults: initial `0.09`, final `0.065` | 1 | 0 | 0 |
| Equal thresholds: initial `0.09`, final `0.09` | 1 | 1 | 0 |

This is a threshold artifact, not synthesized resolution. A separate boundary test with `max_rounds=0` and an initial collision raised `UnboundLocalError` because `complexity_level` was never assigned.

### 2. Unary and relational functions shared an incompatible namespace

The base engine contract is unary: `f(trajectory) -> scalar`. Candidate scoring instead constructs a trajectory difference and evaluates relational functions of `trajectory_a - trajectory_b`. The legacy engine inserts those relational functions into the same registry as unary invariants. `watchtower.py` then evaluates the returned registry on a window-reference difference.

These are different mathematical types with different domains. Treating them as interchangeable makes the meaning of the returned representation depend on the caller.

### 3. Algebraic connectivity used the wrong eigenvalue

`canonical.py` removed zero eigenvalues and then selected `positive[1]`. That usually returns the third eigenvalue of a connected graph and can report a positive value for a disconnected graph. Direct checks gave:

| Graph | Legacy value | Correct normalized-Laplacian lambda-2 |
|---|---:|---:|
| Three-vertex path | 2.0 | 1.0 |
| Four-vertex path | 1.5 | 0.5 |
| Two disconnected two-vertex paths | 2.0 | 0.0 |

The correct implementation retains zero multiplicity, sorts the full spectrum, and reads `eig[1]`.

### 4. The advertised canonical field was not globally similarity invariant

For a deterministic `180 x 3` trajectory and transformed copy `y = 3.7 x Q + b`, 13 of 30 legacy canonical coordinates changed. Examples included:

| Coordinate | Original | Transformed |
|---|---:|---:|
| `geometry.mean_distance` | 2.89793 | 10.7223 |
| `persistence_b1_max_lifetime` | 1.30596 | 4.83206 |
| `temporal.autocorr_lag1` | -0.032658 | -0.387756 |

Only the recurrence-adjacency-derived subset built from distance quantiles had the stated similarity behavior. A claim about the entire field was therefore unsupported.

### 5. Some quantile-graph summaries were construction constants

Across 20 independent continuous random trajectories of shape `80 x 3`, the legacy fast construction produced exactly the same recurrence density (`0.147569444444`) and mean degree (`7.083333333333`) on every run. Those fields mainly encode the selected recurrence quantile rather than trajectory-specific structure.

### 6. The stability gate could not detect instability

The legacy budget fixed `MAX_PERTURBATIONS=1`. Stability was calculated from the variance across perturbations; the variance of one value is zero. Repeated tests gave a stability score of `1.0` for a norm, a first-coordinate projection, and an explosive exponential function alike.

### 7. The candidate search was narrower than the discovery language suggested

The level-one priority list was sorted and capped at five candidates, yielding only predeclared persistence summaries: `b0_count`, `b0_final`, `b0_max_lifetime`, `b0_mean_lifetime`, and `b1_count`. Reported persistence “discoveries” were therefore selections from a small fixed shortlist, not open-ended invariant synthesis.

### 8. Auxiliary-domain transfer was not held out

The synthesis budget allowed one scored auxiliary domain, and the implementation sliced the supplied list to its first member. Provenance could nevertheless list every supplied auxiliary. More importantly, the source and auxiliary domains were combined in the separation and transfer objective. Performance on those same auxiliaries was training performance, not held-out transfer.

The recorded provenance generally used wildcard feature fields and generic names such as `aux_0`, which is insufficient to reconstruct the exact audit boundary.

### 9. The persistence calculation and benchmark proxy were overstated

The recovered persistence module performs a graph component/cyclomatic sweep. It is not persistent first homology of a clique filtration. In addition, one benchmark harness generated `transfer_retention_proxy` from a fixed formula rather than measuring transfer on a held-out domain.

The recovered package contained no committed test directory.

## Legacy seismic pilot

The associated seismic work was also inspected. Its principal scripts were:

- `waveform_voss_benchmark.py`
- `fast_yellowstone_benchmark.py`
- `global_waveform_spotcheck.py`
- `analyze_voss_merge.py`

The original `merge.csv` and `merge.hdf5` inputs were no longer present, so the results cannot be rerun from the retained artifacts. The retained summaries also disagree across tuning stages: one run reports 10/11 pre-P detections, a later run reports 5/11, and the tuned fast run again reports 10/11. In the earlier 10/11 result, the robust baseline also detected 10/11.

Confirmed sources of look-ahead or evaluation leakage include:

- median/MAD normalization over an entire waveform before scoring its pre-arrival prefix;
- full-trace normalization in the fast benchmark;
- retrospective, hard-coded selection of known swarm targets;
- preference for manually reviewed traces and larger arrival indices;
- reuse of baseline/noise traces as negative controls;
- use of the same purported controls for weight learning, threshold calibration, and final control scoring;
- iterative post-hoc tuning on the same eleven targets;
- standardization of a complete chronological catalog before defining an early baseline.

Accordingly, the seismic analysis is classified as a **historical exploratory pilot only**. It is excluded from the Paper III confirmatory metrics and is not described as prospective evidence.

## Reuse boundary

The current Paper III implementation reuses only the conceptual development path: collision pressure, candidate construction, validation, and actual reinsertion. It does not reuse the legacy engine's collision counts, stability scores, transfer claims, Watchtower performance, persistence labels, or seismic performance numbers.

