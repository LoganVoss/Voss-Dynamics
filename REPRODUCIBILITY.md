# Reproducibility

## Release boundary

This repository is a source-and-evidence release for three manuscripts. Paper I and Paper II report proofs and internal computations; Paper II additionally specifies a candidate protocol. Paper III reports proofs, deterministic checks, internally frozen synthetic audits, and two null/mixed real-data stress tests.

No artifact in this repository is an external registration timestamp. Git history and SHA-256 manifests establish release integrity after publication, not prospective commitment before project development.

## Requirements

- macOS or Linux
- Python 3.11 or later for Papers II and III
- [uv](https://docs.astral.sh/uv/) for locked environments
- [Tectonic](https://tectonic-typesetting.github.io/) for PDF builds
- shasum or sha256sum for integrity verification

The committed PDFs are the canonical release builds. TeX engine and plotting-backend metadata can change binary hashes even when page content is unchanged.

## Paper I

~~~bash
cd "(1) The Principle of Full Invertibility"
uv sync --frozen
uv run python evidence/thesis_verification.py
uv run python experiments/run_battery.py
tectonic The_Principle_of_Full_Invertibility.tex
~~~

Expected manuscript filename: The_Principle_of_Full_Invertibility.pdf.

Paper I’s numerical outputs are simulated model computations. They are not measurements of an external physical system.

## Paper II

~~~bash
cd "(2) Emergent Predictive Representation"
uv sync --frozen --all-groups
uv run pytest -q
uv run python -m evidence.run_all
cd thesis
tectonic main.tex
cp main.pdf ../Emergent-Predictive-Representation.pdf
cd ..
~~~

Expected manuscript filename: Emergent-Predictive-Representation.pdf.

The evidence suite should pass 13 tests. Synthetic target records exercise algebra, baselines, and the calibration–freeze–reveal software path; they do not validate VD-Hopf-1 in Nature.

## Paper III

~~~bash
cd "(3) Pressure-Driven Invariant Synthesis"
sh evidence/data/download_ucr.sh
uv sync --frozen
uv run pytest -q
uv run python evidence/run_all.py
shasum -a 256 -c evidence/outputs/manifest.sha256
cd thesis
tectonic main.tex
cp main.pdf ../Pressure-Driven-Invariant-Synthesis.pdf
cd ..
~~~

Expected manuscript filename: Pressure-Driven-Invariant-Synthesis.pdf.

The evidence suite should pass 10 tests. The full run is deterministic under the pinned environment and declared seed, but it is substantially slower than the unit suite because it computes recurrence and persistent-homology features across source, audit, and real records. After the frozen Kuramoto result is evaluated, the runner also regenerates the exact audit with latent phases retained and writes an explicitly post-hoc mechanistic diagnostic; that step does not participate in program selection.

### Paper III freeze semantics

Within the final runner:

1. source synthetic trajectories are generated and featurized;
2. grammar and declared split manifests are written;
3. programs are selected and written to frozen_programs.json;
4. only then are independent synthetic audit trajectories and UCR TEST cases loaded and featurized.

The external timestamp is null. The [development audit](./%283%29%20Pressure-Driven%20Invariant%20Synthesis/evidence/DEVELOPMENT_AUDIT.md) records an earlier full-pipeline null failure and the gate redesign that followed. Final audit seeds were changed after that redesign, but the project is not prospectively registered.

## Unified volume

The series source builds the requested two-page front matter. A pinned PDF assembler then appends the three canonical papers while preserving their link annotations and importing their outlines:

~~~bash
cd series
tectonic main.tex
uv run --with pypdf==6.10.0 python build_volume.py
~~~

Expected output: Voss-Dynamics-Information-Representation-and-Discovery.pdf.

## Integrity and provenance

- Paper-specific evidence manifests cover committed numerical outputs and figures.
- Paper III’s raw UCR archive hashes and selected TEST identifiers are machine-readable.
- Root release hashes are generated at publication time and identify the exact public bundle.
- Raw or third-party data remain governed by upstream terms.
- The historical seismic pilot is excluded from Paper III’s empirical result totals; see its forensic audit.

## Reporting a mismatch

Open an issue with:

- operating system and architecture;
- Python, uv, Tectonic, and dependency versions;
- exact command and full traceback/log;
- repository commit;
- the first mismatching checksum or metric.

Do not silently regenerate and commit changed evidence. A numerical difference should be explained before any manifest is updated.
