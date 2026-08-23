# Contributing

Contributions that improve correctness, reproducibility, or claim discipline are welcome.

## High-value contributions

- independent reproduction on a clean environment;
- a minimal failing case for a theorem statement or implementation;
- leakage audits and genuinely independent split designs;
- prospectively frozen replication on a new domain;
- compute-matched baselines;
- corrections to bibliography, provenance, or licensing;
- PDF accessibility and layout fixes.

## Before opening a pull request

1. Identify the claim status affected: proved, computed, candidate/internally held out, or open.
2. Add or update a test when code behavior changes.
3. Do not update a checksum manifest until the numerical difference is explained.
4. Keep external data out of the repository unless redistribution is clearly permitted and the file is small.
5. State whether any target, audit, or holdout data influenced method selection.
6. Preserve superseded failures in the development audit when they materially shaped the method.

No pull request should convert a synthetic result into a physical claim, an internal freeze into preregistration, or a named-domain result into universal transfer.

