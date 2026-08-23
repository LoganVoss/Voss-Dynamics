# Emergent Predictive Representation

This repository contains a conditional mathematical thesis, reproducible internal checks, and a draft bounded-model experiment for the Voss Dynamics program.

The central result is deliberately conditional:

- a normalized nontrivial stationary real pointer has a two-dimensional irreducible cyclic span;
- a binary Euclidean-Jordan face of Peirce degree `d` has an unnormalized QND-filter algebra of dimension `1 + d(d - 1)/2`;
- an explicit accessibility plus injective-surjective response axiom forces `d = 2`, the complex Hermitian family;
- a pure qubit admits the redundant Hopf lift `Q(ψ) = ψψ†`; the theorem yields the ray/projector, not a physical point of `S³`;
- the QND transverse invariant is encoded in `ρ`;
- VD-Hopf-1 separately postulates a physical closed-history holonomy memory and one orientation-odd record term.

No physical deviation from quantum mechanics is claimed. Synthetic data exercise only the internal identities and Gaussian calibration/holdout software. A complete instrument, mixed/composite theory, bounded process tensor, qutrit–cavity null, and path score remain open.

## Main deliverables

- `Emergent-Predictive-Representation.pdf` — finished thesis PDF
- `thesis/main.tex` — LaTeX source used to build the paper
- `thesis/references.bib` — primary-source bibliography
- `evidence/` — symbolic, numerical, statistical, and plotting code
- `experiment/PREREGISTRATION.md` — non-executable draft protocol and missing fields
- `evidence/outputs/frozen_prediction.json` — machine-readable candidate prediction
- `evidence/outputs/manifest.sha256` — local release-integrity hashes (not an external timestamp)

## Reproduce

```bash
uv sync --frozen --all-groups
uv run pytest
uv run python -m evidence.run_all
tectonic -X compile thesis/main.tex --outdir thesis/build
```

The evidence seed is `20260821`. Generated target records are marked `synthetic_only: true`.

## Evidence boundary

| Status | Meaning |
|---|---|
| Proved | Follows from stated mathematical assumptions |
| Computed | Reproduced by deterministic code or seeded simulation |
| Candidate | Explicit falsifiable law, not inferred from Nature |
| Open | Requires hardware, stronger axioms, or independent replication |

The proposed coupling parameter is `ε`. The null `ε = 0` is the **memoryless reduced-qubit null**, not all of standard quantum theory.
