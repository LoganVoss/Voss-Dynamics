# Evidence Data

Paper III uses two public univariate archives from the official UCR/UEA Time Series Classification repository. The public Voss Dynamics repository does not redistribute the archives; the checksum-enforcing downloader acquires them from the official URLs. The evidence runner reads only each archive's `*_TEST.ts` member; it does not load the TRAIN split. From each TEST split it selects a deterministic balanced subset of at most 32 cases per class using the seeds declared in `evidence/outputs/declared_split_manifest.json`.

## Archive provenance

| Archive | Official source | Bytes | SHA-256 |
|---|---|---:|---|
| `ECGFiveDays.zip` | [UCR/UEA archive](https://www.timeseriesclassification.com/aeon-toolkit/ECGFiveDays.zip) | 1,523,014 | `11457a2d590711598eac1a0ab5c58d43c5e3c5c2c86521809d69f7c0b6b3edd1` |
| `Earthquakes.zip` | [UCR/UEA archive](https://www.timeseriesclassification.com/aeon-toolkit/Earthquakes.zip) | 781,477 | `927cfb732988055850a74efa169dfd633bc7f578095c35319bebb83453501bf1` |

The archive catalog and dataset documentation are maintained at [timeseriesclassification.com](https://www.timeseriesclassification.com/dataset.php).

## Download and verification

From the Paper III directory, run:

```sh
sh evidence/data/download_ucr.sh
```

The script downloads into `evidence/data/raw/`, verifies the expected SHA-256 before installing each file, skips an already-correct archive, and refuses to overwrite an existing archive whose checksum differs.

To verify the archives independently on macOS:

```sh
shasum -a 256 evidence/data/raw/ECGFiveDays.zip evidence/data/raw/Earthquakes.zip
```

On systems with GNU Coreutils:

```sh
sha256sum evidence/data/raw/ECGFiveDays.zip evidence/data/raw/Earthquakes.zip
```

The expected hashes are also recorded in `evidence/outputs/declared_split_manifest.json` and covered by `evidence/outputs/manifest.sha256`.

## Evaluation boundary

These archives are the only real-data inputs to the current Paper III evidence pipeline. They were loaded after the all-synthetic program was frozen within the run. This is an internal frozen-holdout design, not an externally registered experiment. The legacy seismic pilot and its unavailable raw waveform/catalog inputs are excluded.
