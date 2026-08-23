# Evidence package — The Principle of Full Invertibility

This directory contains the executable tests, machine-readable ledgers, and
plot data used by the paper. Its purpose is to make every quantitative claim
traceable to code and a reproducible artifact.

## Reproduce

Run from the Paper I directory:

```bash
uv sync --frozen
uv run python experiments/run_battery.py
uv run python experiments/run_deeper.py
uv run python experiments/run_omega.py
uv run python evidence/thesis_verification.py
uv run python evidence/export_plot_data.py
tectonic -X compile The_Principle_of_Full_Invertibility.tex --outdir build
```

The final run used Python 3.9.6, NumPy 2.0.2, Tectonic 0.17.0, and arm64 macOS.

## Core evidence

1. **Full determinant factorization**

   \[
   \det D\Psi=\gamma^{2N}\det D_\phi G_{\mathbf r^+}.
   \]

   Fifteen live-map finite-difference probes agree with the factorized value;
   maximum relative discrepancy is approximately \(2.01\times10^{-8}\).

2. **Exact two-clock law**

   \[
   \delta^+=\delta+\Delta\omega\,\Delta t
   -\frac{2\beta\Delta t}{D}\sin\delta \pmod{2\pi}.
   \]

   Two hundred random probes have maximum angular residual
   \(1.55\times10^{-15}\).

3. **Global phase branch structure**

   At seed 22, step 9, one phase target has three refined roots. Combining each
   root with the algebraic translational inverse produces three complete
   pre-states. The alternative states lie \(7.386\) and \(6.738\) from the
   generating pre-state and return to the target within
   \(1.6\times10^{-13}\).

4. **Certified float64 round-trip recovery**

   Twelve seeds, two values of \(\gamma\), and six horizons produce 144
   endpoint inversions. Displayed quantiles include only trajectories with
   \(q<1\) at every step and wrapped phase residual below \(10^{-12}\).

5. **Description fibres**

   Across 160 complete states, the geometric phrase readout produces 103
   phrases. Every one of 68 checked shared-phrase pairs remains a pair of
   distinct complete states.

6. **Discrete entropy and fixed observation partitions**

   A nonuniform distribution preserves Shannon entropy under the modular
   permutation to \(3.6\times10^{-15}\) bits. The fixed-partition ensemble
   diagnostic reports sums of coordinate marginals separately from exact
   distinction cardinality.

7. **Relational phase transport**

   After coupling is set exactly to zero, equal-frequency phase difference is
   conserved and unequal-frequency drift follows the exact linear law. A
   reset/no-reset control matches the analytic \(1/\widetilde D\) response.

8. **Boundary provenance**

   Two in-domain one-token states at distance \(0.2\) become identical under
   the wall clamp, locating the distinction loss at the boundary operation.

## Plot tables

- `clock_memory.csv`
- `clock_readout.csv`
- `information_horizon.csv`
- `information_horizon_gamma1.csv`
- `information_horizon_gamma0982.csv`
- `quantized_entropy.csv`
- `determinant_factorization.csv`

Additional CSVs preserve exploratory diagnostics. All plot tables are generated
from `results/*.json` and `thesis_verification.json` by
`export_plot_data.py`.

## Public repository

The paper directory is the complete publication package: manuscript PDF and
source, engine, experiments, raw JSON ledgers, generated CSVs, checksums, and
reproduction instructions. No separate submission archive is required.
