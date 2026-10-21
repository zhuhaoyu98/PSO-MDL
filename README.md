# PSO-MDL

[中文说明](README_zh.md)

## Overview

PSO-MDL selects molecular-descriptor subspaces by particle swarm optimization
and ranks candidates by their equal-weight mean rank across independent models.
The study concerns PBP2-related antibacterial screening. A rank is a relative
descriptor-space priority, not a prediction of antibacterial activity or binding.

## Requirements

MATLAB R2024a or later, Statistics and Machine Learning Toolbox, the default
JVM, and working PNG/PDF graphics export. Dragon is not required: the four
input groups below already contain the computed descriptors and identities.
Allow sufficient disk space for approximately 10.1 GB of inputs and the newly
generated results. Formal training is computationally intensive.

## Usage

Download **all** repository files, including all 100 descriptor parts. Extract
the project to any writable folder and set MATLAB's Current Folder to that
project root. No absolute-path edits, setup script, historical models, or
precomputed intermediate files are needed for a first run:

```matlab
step00_prepare_dataset
step01_preprocess_descriptors
step02_train_models
step03_rank_candidates
step04_generate_figures
```

| Step | Reads | Produces |
|---|---|---|
| 00 | The four input groups listed below | Identity audits, frozen background/holdout partition, descriptor missingness |
| 01 | Step 00 outputs and raw descriptors | Reference/background matrices and fitted preprocessing |
| 02 | Step 01 matrices and preprocessing metadata | 64 newly trained models, optimization histories and statistics |
| 03 | Steps 00–02 outputs and raw candidate descriptors | `candidate_ranking.tsv`, single-model ranks and stability statistics |
| 04 | This run's Step 02/03 outputs | Figure 3A–C and Figure 4A–E, each as PNG and PDF |

All generated files go to `results/research/`, relative to the project root.
Each Step reads existing upstream outputs; none silently loads historical
models or substitutes simulated inputs. Step 04 plots the current run.
For a separate repeat run, change only `cfg.Output.Directory` in
`functions/pso_mdl_config.m` to a new directory. Existing outputs are protected
by default. Other input paths and research parameters are configured in that
same file if needed.

The formal settings remain 16 references, 20,000 background molecules,
306,922 screening candidates, 2,781 retained descriptors, and 64 seeds × 3,000
PSO generations. P is the mean non-medoid proximity percentile (0–100), with
J = −RFP. Internally, PSO compares the unscaled key
K = −RF sum(p_t) and exports J = K/(nRef−1), preserving the historical
floating-point comparison path without changing the mathematical objective.

All 64 models and the complete candidate ranking were reproduced exactly
from scratch on NIRO using MATLAB R2024a, Intel MKL 2023.2 (AVX2), and 28
computational threads. The default configuration requests the same thread
limit and records the actual environment. Different CPU, MATLAB or BLAS
environments can run the workflow but are not guaranteed bitwise-identical
optimization trajectories. A BLAS-branch mismatch is reported as a warning.

## Data

Only the four original input groups are retained under `data/private/`:

```text
data/private/
  candidate_descriptor_parts/   # 100 numbered TXT parts + descriptor_parts.json
  candidate_manifest.tsv
  reference_descriptors.txt
  reference_manifest.tsv
```

Keep every numbered part and its JSON manifest unchanged. The code streams
the parts as one table; no merged TXT or external data directory is needed.
The folder name `private` is retained for compatibility; these exact input
paths are now allowlisted in `.gitignore`. Other data and generated outputs
remain ignored. No saved models, previous rankings, test inputs, or historical
paper-result bundle are included.

Licensed Dragon 7 calculated all 2D descriptors, with 3D disabled and default
save-all export; no descriptor exclusion/filtering was applied in Dragon.
Both reference and candidate tables contain 3,855 descriptors. Subsequent
MATLAB quality control removes 18 descriptors exceeding 5% candidate
missingness and 1,056 with near-zero/invalid scales, leaving 2,781.
See [data/README.md](data/README.md) for input and output details.

The input payload is approximately 10.1 GB; its files are below 100 MiB but
exceed the browser-upload limit. Splitting does not remove GitHub's push-size
or repository-size constraints. Use an appropriate Git-based data distribution
plan and confirm redistribution permission before publishing third-party
derived data. See [GitHub's large-file guidance](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github)
and [repository limits](https://docs.github.com/en/repositories/creating-and-managing-repositories/repository-limits).
The local release has not been uploaded automatically.

## Citation and License

Please cite Zhu Haoyu, **PSO-MDL, version 1.0.0**, using
[CITATION.cff](CITATION.cff). The [MIT License](LICENSE) covers the author's
original code, not third-party molecular data, MATLAB, Dragon, or other software.
