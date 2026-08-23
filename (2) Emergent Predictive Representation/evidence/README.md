# Evidence package

This directory separates exact identities, internal numerical regressions, and a synthetic Gaussian protocol smoke test.

## Files

- `epr_core.py` — state maps, Hopf geometry, QND/Jordan calculations, quantum-limit checks, candidate law, fibre variation, power, and simulation
- `run_all.py` — deterministic end-to-end run and figure generation
- `tests/test_epr_core.py` — mathematical and numerical regression tests
- `outputs/metrics.json` — machine-readable results
- `outputs/frozen_prediction.json` — frozen synthetic validation specification
- `outputs/frozen_predictions.csv` — calibration and synthetic holdout cells
- `outputs/power_table.csv` — sample-size and information calculations
- `outputs/manifest.sha256` — local release-integrity hashes, not an external timestamped commitment

## Run

```bash
uv run pytest
uv run python -m evidence.run_all
```

## Interpretation

The scalar-selection result is conditional on Euclidean Jordan geometry and an accessibility plus injective-surjective map onto an unnormalized filter algebra. The code checks the dimension arithmetic, not the physical axiom or Jordan classification.

VD-Hopf-1 is a pure-state phenomenological record ansatz, not a complete instrument. For a fixed protocol and matched bounded standard-process state, its short-probe objective is

\[
V^\star_{{\rm p},\eta_0}
\simeq2\Phi(a|\varepsilon|\sqrt{T})-1.
\]

Positive variation in a simulator containing `ε ≠ 0` proves only that the ansatz is internally distinguishable. Zero from a finite search is not a no-go theorem. The current protocol file explicitly lists the theory, model, statistics, and commitment work still required before hardware.

The suite also encodes an ordinary-quantum countermodel: coherent cavity memories `|+α⟩` and `|-α⟩` have equal photon number and the same reduced qubit state but different homodyne laws. This is why a positive reduced-state witness is not automatically beyond-quantum evidence.
