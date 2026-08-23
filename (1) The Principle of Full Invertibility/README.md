# The Principle of Full Invertibility

Paper I of Voss Dynamics asks a deceptively simple question: when two possible states become indistinguishable, **which map merged the distinction?**

The paper separates five operations that are often conflated: dynamical evolution, observation, numerical representation, inverse selection, and boundary handling. The result is a provenance calculus for information loss. A transition can be locally nonsingular but globally many-to-one; an observation can merge states even when the dynamics are bijective; finite precision can hide a distinction without destroying it in the exact model; and a boundary clamp can erase it explicitly.

## Main results

- On the admissible drift–kick stratum,
  $$
  \det D\Psi=\gamma^{2N}\det D_\phi G_{\mathbf r^+},
  $$
  with the force Jacobian cancelling from the determinant.
- An explicit target has three distinct complete phase preimages, separating local nonsingularity from global invertibility.
- A finite modular realization is bijective exactly when its velocity multiplier is a unit modulo the state-space modulus.
- A computed census maps 160 complete states to 103 descriptions while all 68 checked shared-description pairs remain dynamically distinct.
- Transition folds, observation fibres, numerical access, and boundary clipping are assigned different mathematical types and different recovery obligations.

Read the finished [paper](The_Principle_of_Full_Invertibility.pdf), then use the manuscript’s Appendix C and the [evidence guide](evidence/README.md) as the claim-to-artifact ledger.

## Reproduce

~~~bash
uv sync --frozen
uv run python evidence/thesis_verification.py
uv run python experiments/run_battery.py
shasum -a 256 -c evidence/SHA256SUMS.txt
~~~

The committed numbers are internal computations of the stated models, not measurements of an external physical system. The manifest verifies the canonical release artifacts; regenerating an artifact can change its hash when environment metadata changes.

## Package map

- `The_Principle_of_Full_Invertibility.pdf` — canonical paper
- `The_Principle_of_Full_Invertibility.tex` and `sections/` — source
- `engine/` — exact-state dynamics, observation, and reconstruction primitives
- `experiments/` and `results/` — seeded computational studies
- `evidence/` — verification program, machine-readable outputs, plot tables, and checksums
- `uv.lock` — pinned Python dependency resolution

The paper establishes a method for locating and classifying a lost distinction. It does **not** claim a universal physical conservation law for information.
