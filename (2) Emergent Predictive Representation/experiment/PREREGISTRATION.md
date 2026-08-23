# VD-Hopf-1 draft protocol

Version: `VD-Hopf-1`  
Status: **not executable as a confirmatory preregistration**  
Implemented today: synthetic cell-level Gaussian smoke test  
Still required: mixed-state/instrument theory, bounded qutrit–cavity and process-tensor nulls, numerical thresholds, path score, and external commitment

## Scope

Question: do oppositely oriented closed control histories that return a qubit to the same **reduced** density matrix produce different later records after a declared bounded reset and memory model?

A positive result first establishes insufficiency of that reduced-state/process model. It cannot exclude standard quantum mechanics with a larger reference, cavity, controller, leakage manifold, or environment.

## Candidate history variable

For a periodic pure-state lift,

\[
\Gamma_{\rm rel}[C]
=
\oint_C \mathcal A\pmod{2\pi},
\qquad
\mathcal A=+i\langle\psi|d\psi\rangle.
\]

With the orientation convention used here,

\[
\Gamma_{\rm rel}[C]=-\Omega[C]/2.
\]

For a nonperiodic lift, the endpoint Pancharatnam term must be included. Dynamical phase and the physical reference path must be calculated explicitly.

**Holonomy-memory postulate:** some physical carrier retains `Γrel` through the declared delay/reset while every variable represented in the frozen bounded standard-quantum process state `η` is matched between arms within specified confidence margins. Unrepresented modes remain a limitation. Geometry does not prove this postulate.

The current ansatz is restricted to conditioned pure, efficient trajectories. A mixed-state/purification extension or a frozen near-purity/no-jump selection rule with quantified bias is required before hardware.

## Models

Memoryless reduced-qubit null:

\[
P(Y_{0:T}\mid C_+,\rho_0)
=
P(Y_{0:T}\mid C_-,\rho_0).
\]

VD-Hopf-1 phenomenological record law:

\[
dY_t
=
a\left[
z_t+\varepsilon(1-z_t^2)\sin\Gamma_{{\rm rel},t}
\right]dt+dW_t,
\]

\[
ds_t
=
a\,dY_t
-a^2\varepsilon(1-z_t^2)\sin\Gamma_{{\rm rel},t}\,dt,
\qquad
z_t=\tanh s_t.
\]

This keeps the declared projected population as the standard QND martingale,

\[
dz_t=a(1-z_t^2)dW_t,
\]

but it is not a derived mixed-state, completely positive, or composable instrument.

For mirror loops with initial offset `Γ0`,

\[
\Delta_\star(\Omega)
=
-2a\varepsilon\cos\Gamma_0\sin(\Omega/2)
\left[
\frac1T\int_0^T\mathbb E(1-z_t^2)\,dt
\right].
\]

The synthetic code fixes `Γ0 = 0` and uses the short-probe limit.

## Required explicit control construction

Before registration, supply:

- control Hamiltonians `H+(t)` and `H-(t)`;
- cyclic-state and realized-path calculations;
- total, dynamic, and geometric phase calculations;
- the physical reference path;
- reset and delay transformation laws for `Γrel` and every standard memory mode;
- independently compiled physical implementations of orientation reversal.

The histories must match, within numeric confidence margins:

- initial and final reduced qubit states;
- qutrit leakage populations;
- duration, pulse count, spectra, and integrated power;
- relaxation and dephasing exposure;
- residual dynamical phase;
- complex cavity amplitude and covariance, not photon number alone;
- measured controller/reference summaries represented in the bounded null;
- measured readout-chain summaries and calibration, each with a confidence margin.

Sham loops use matched resources and zero enclosed area.

## Experimental blocks

1. Prepare the frozen near-pure state and log the purity/selection criterion.
2. Randomly assign `C+`, `C-`, sham, or hardware-control arms in short interleaved drift blocks.
3. Execute the explicit loop.
4. On a non-target subset, run SPAM-robust tomography or GST and a one-sided equivalence test.
5. Characterize qutrit leakage and complex cavity residue.
6. Apply the fixed reset/ringdown procedure.
7. Acquire a short nominally `Z`-QND record with no Rabi drive.
8. Preserve raw IQ data and all timing, control, temperature, and calibration metadata.

Randomized benchmarking is auxiliary; it cannot establish final-state equality.

## Bounded null models

Freeze before target reveal:

1. memoryless reduced-qubit model, `ε = 0`;
2. qutrit plus measured cavity-mode model, including level-dependent IQ response and measurement-induced transitions;
3. completely positive process tensor with fixed:
   - time grid,
   - intervention basis,
   - system/leakage dimension,
   - memory length or MPO bond dimension,
   - regularization,
   - stationarity assumptions,
   - confidence construction,
   - out-of-distribution prediction rule;
4. pulse-spectrum, IQ imbalance, transfer-function, and signed-leakage models;
5. frozen multi-mode ringdown/coherent-spectator model;
6. bounded low-rank history baselines.

An unrestricted history model or process tensor can reproduce the sine law and is not a falsifiable null.

## Data partition and commitment

Before acquisition, assign records and solid-angle cells to:

- apparatus calibration;
- coefficient calibration;
- held-out validation;
- diagnostics excluded from the primary score.

Only coefficient calibration may estimate `ε`. Preprocessing, exclusions, covariance, model capacity, stopping, and missing-data rules freeze before validation reveal.

`evidence/outputs/frozen_prediction.json` is a synthetic software artifact. Its local SHA-256 manifest is a release-integrity snapshot, not an external timestamped commitment. A physical study needs an independent registry or cryptographic timestamp.

## Implemented synthetic fit

At `Γ0 = 0`, calibration uses

\[
\widehat\varepsilon
=
\arg\min_\varepsilon
\sum_{j\in{\rm cal}}
\frac{
[\widehat\Delta_j
+2a\varepsilon(1-z_0^2)\sin(\Omega_j/2)]^2
}{\sigma_j^2}.
\]

Held-out cells use

\[
\Delta_j^{\rm frozen}
=
-2a\widehat\varepsilon(1-z_0^2)
\sin(\Omega_j^{\rm hold}/2).
\]

The code implements a Gaussian summary log score and propagates the shared calibration variance. It does **not** implement path likelihood or process-tensor comparison.

## Protocol-level fibre variation

Fix protocol `p = (a,T,instrument)` and all bounded standard process variables `η0`:

\[
\mathcal F_{\rho,\eta_0}
=
\{(\rho,\Gamma_{\rm rel},\eta_0):
\Gamma_{\rm rel}\in S^1\}.
\]

\[
V_{\rm p}(\rho;\eta_0)
=
\sup_{\Lambda_1,\Lambda_2\in\mathcal F_{\rho,\eta_0}}
D_{\rm TV}\!\left[
P_{\rm p}(\cdot\mid\Lambda_1),
P_{\rm p}(\cdot\mid\Lambda_2)
\right].
\]

In the short Gaussian ansatz,

\[
V^\star_{{\rm p},\eta_0}
\simeq
2\Phi(a|\varepsilon|\sqrt T)-1.
\]

An experiment with two loops estimates a lower bound under state-equivalence assumptions; it cannot empirically optimize over all fibres. A finite numerical search returning zero is not a no-go theorem.

## Oracle power only

For independent white Gaussian records and known nuisance parameters,

\[
\operatorname{SE}(\widehat\Delta)=\sqrt{\frac{2}{NT}},
\]

\[
N
\ge
\frac{2(z_\alpha+z_{1-\beta})^2}
{4a^2\varepsilon^2(1-z_0^2)^2T
\cos^2\Gamma_0\sin^2(\Omega/2)}.
\]

This omits drift/autocorrelation, calibration uncertainty, nuisance fitting, multiplicity, exclusions, tomography, and process-model comparison. Final power must use blockwise simulation of the frozen qutrit–cavity and process models.

## Numeric fields that must be frozen before execution

- minimum effect of interest;
- total sample count and allocation;
- one- or two-sided `α`, multiplicity, and power;
- block length and effective-sample/autocorrelation rule;
- trace-distance and leakage equivalence margins;
- largest null record difference allowed by those margins;
- path-log-score model-comparison threshold;
- process-tensor confidence set and capacity;
- negative-control margins;
- exclusion, missing-data, and stopping rules;
- replication criterion.

## Bounded-model stage pass

After all numeric fields exist, every condition must hold:

- final reduced-state and leakage equivalence margins pass;
- sham/zero-area controls are null;
- sign reverses under loop reversal;
- equal-area distinct paths agree within tolerance;
- area, `z0`, `Γ0`, and delay scaling match the frozen law;
- ringdown, leakage, pulse, and IQ controls satisfy their margins;
- VD-Hopf-1 beats the frozen reduced-qubit, qutrit–cavity, and finite-memory models by the threshold;
- the bounded process-tensor prediction fails outside its confidence set.

This stage establishes unexplained history dependence relative to the listed models. Independent replication is a separate confirmation stage. It is not evidence against every possible quantum environment.

## Failure

Any of the following rejects VD-Hopf-1 at the tested scale:

- `ε` is consistent with zero at the achieved sensitivity;
- final reduced states or leakage differ enough to explain the record;
- a bounded ordinary quantum model predicts validation records;
- signal scaling fails;
- analysis changes after reveal;
- confirmation fails.

Failure leaves separate branches: conditional reconstruction, modified dynamics of a sufficient `ρ`, or searches in temporal/composite regimes. It is not positive evidence for those branches.
