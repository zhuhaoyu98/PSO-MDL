# PSO-MDL

## Overview

PSO-MDL selects molecular-descriptor subspaces by particle swarm optimization
and ranks candidates by their equal-weight mean rank across independent models.
The study concerns PBP2-related antibacterial screening.

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
Each Step reads existing upstream outputs. Step 04 plots the current run.
For a separate repeat run, change only `cfg.Output.Directory` in
`functions/pso_mdl_config.m` to a new directory. Existing outputs are protected
by default. Other input paths and research parameters are configured in that
same file if needed.

The formal settings remain 16 references, 20,000 background molecules,
306,922 screening candidates, 2,781 retained descriptors, and 64 seeds × 3,000
PSO generations. P is the mean non-medoid proximity percentile (0–100), with
J = −RFP. Internally, PSO compares the unscaled key
K = −RF sum(p_t) and exports J = K/(nRef−1), without changing the
mathematical objective.

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

Licensed Dragon 7 calculated all 2D descriptors, with 3D disabled and default
save-all export; no descriptor exclusion/filtering was applied in Dragon.
Both reference and candidate tables contain 3,855 descriptors. Subsequent
MATLAB quality control removes 18 descriptors exceeding 5% candidate
missingness and 1,056 with near-zero/invalid scales, leaving 2,781.
See [data/README.md](data/README.md) for input and output details.

## Helper functions

The table lists every MATLAB helper file under `functions/`.

| File (relative to `functions/`) | Purpose |
|---|---|
| [pso_mdl_config.m](functions/pso_mdl_config.m) | Define input/output paths, preprocessing settings, PSO parameters and figure options. |
| [resolve_pso_mdl_config.m](functions/resolve_pso_mdl_config.m) | Resolve configuration and project-relative paths, initialize helper folders and check runtime dependencies. |
| [io/assert_output_available.m](functions/io/assert_output_available.m) | Prevent unintended overwriting of existing outputs unless explicitly enabled. |
| [io/audit_descriptor_missingness.m](functions/io/audit_descriptor_missingness.m) | Count descriptor-wise NaN and Inf values while checking streamed table alignment. |
| [io/compute_file_sha256.m](functions/io/compute_file_sha256.m) | Calculate a file's SHA-256 digest in memory-efficient blocks. |
| [io/count_text_data_rows.m](functions/io/count_text_data_rows.m) | Count nonempty data rows after the header of a text table. |
| [io/descriptor_table_stream.m](functions/io/descriptor_table_stream.m) | Stream a single TXT table or ordered descriptor parts, validating their structure. |
| [io/load_paper_results.m](functions/io/load_paper_results.m) | Validate and load a supplied saved-result bundle; this optional bundle is not included. |
| [io/match_input_files.m](functions/io/match_input_files.m) | Match configured file patterns, including empty optional input paths. |
| [io/pso_mdl_paths.m](functions/io/pso_mdl_paths.m) | Create run-specific results, cache, audit, figure and status directories. |
| [io/read_descriptor_columns.m](functions/io/read_descriptor_columns.m) | Read requested descriptor columns for all molecules, preserving column order. |
| [io/read_descriptor_header.m](functions/io/read_descriptor_header.m) | Validate table headers and return ordered identifier and descriptor names. |
| [io/read_descriptor_rows.m](functions/io/read_descriptor_rows.m) | Read requested molecule rows with all raw descriptors, preserving row order. |
| [io/resolve_input_file.m](functions/io/resolve_input_file.m) | Resolve a required or optional input file and reject ambiguous matches. |
| [io/write_json.m](functions/io/write_json.m) | Write structured metadata as indented UTF-8 JSON. |
| [preprocessing/deduplicate_candidate_manifest.m](functions/preprocessing/deduplicate_candidate_manifest.m) | Apply identity-based deduplication and reference-overlap rules, recording exclusion reasons. |
| [preprocessing/select_background_partition.m](functions/preprocessing/select_background_partition.m) | Select background molecules and separate screening/excluded records with reproducible membership. |
| [preprocessing/weighted_descriptor_indices.m](functions/preprocessing/weighted_descriptor_indices.m) | Identify descriptors on the name-based emphasis list for Step 01 scaling. |
| [scoring/build_fitness_context.m](functions/scoring/build_fitness_context.m) | Fix the reference medoid and cache per-descriptor squared distance components. |
| [scoring/column_percentile.m](functions/scoring/column_percentile.m) | Calculate a linearly interpolated percentile for each column of a numeric matrix. |
| [scoring/evaluate_pso_mdl_fitness.m](functions/scoring/evaluate_pso_mdl_fitness.m) | Evaluate J = -RFP and return objective components and the PSO comparison key. |
| [optimization/configure_numerical_execution.m](functions/optimization/configure_numerical_execution.m) | Set the computational thread limit and record MATLAB/BLAS environment details. |
| [optimization/run_sign_pso.m](functions/optimization/run_sign_pso.m) | Optimize a descriptor mask using seeded, sign-encoded PSO and record optimization histories. |
| [ranking/aggregate_candidate_ranks.m](functions/ranking/aggregate_candidate_ranks.m) | Combine equal-weight model ranks and summarize cross-model rank uncertainty. |
| [ranking/score_candidates_streaming.m](functions/ranking/score_candidates_streaming.m) | Stream descriptors, apply frozen preprocessing and compute per-model distances to the medoid. |
| [visualization/export_publication_figure.m](functions/visualization/export_publication_figure.m) | Export figures as resolution-controlled PNG and vector PDF files. |
| [visualization/plot_saved_paper_figures.m](functions/visualization/plot_saved_paper_figures.m) | Redraw Figures 3-4 from a supplied saved-result bundle without rerunning PSO. |

## Citation and License

Please cite Zhu Haoyu, **PSO-MDL, version 1.0.0**, using
[CITATION.cff](CITATION.cff). The [MIT License](LICENSE) covers the author's
original code, not third-party molecular data, MATLAB, Dragon, or other software.
