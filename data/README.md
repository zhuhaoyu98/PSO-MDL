# Input data and generated outputs

## Four input groups

Default input paths are relative to this project, not to a particular computer:

| Path under data/private/ | Contents |
|---|---|
| candidate_descriptor_parts/ | 100 ordered TXT parts and descriptor_parts.json; 328,841 candidate records × 3,855 descriptors |
| candidate_manifest.tsv | Candidate identities in exactly the descriptor-record order |
| reference_descriptors.txt | 16 reference records × the same 3,855 descriptors |
| reference_manifest.tsv | Reference identities in exactly the descriptor-record order |

Descriptor files are UTF-8, tab-delimited, with `No.`, `NAME`, then descriptor
columns. Candidate and reference descriptor names and order must match exactly.
Missing values are not zero. Identity manifests use `INCHIKEY`; candidate
`DragonName` must match `NAME`, and `IDNUMBER` is used for tie-breaking.

Keep all `candidates_part_001.txt` through `candidates_part_100.txt` and
`descriptor_parts.json`. Each part repeats the unchanged header; `No.` remains
global across parts. The manifest records ordering, row ranges, byte counts and
SHA-256 values. The reader checks inventory, headers, sizes and row continuity;
Step 00 also checks all part hashes. Do not modify or renumber individual parts.
Git attributes preserve the original bytes of all four input groups on checkout.

## Descriptor calculation and preprocessing

Licensed Dragon 7 calculated all 2D descriptors, with 3D disabled, using default
save-all export and no descriptor exclusion/filtering during calculation or
export. Both raw exports contain 3,855 descriptors. Subsequent **MATLAB** quality
control removes 18 descriptors exceeding 5% candidate missingness and 1,056
with near-zero/invalid scales, retaining 2,781.

Step 00 deduplicates identities, excludes reference overlaps and fixes the
background/holdout split. Background matching uses MW, MLOGP, TPSA(Tot), nHDon,
nHAcc and RBN; the study selects 20,000 background molecules and retains
306,922 screening candidates after 1,919 other exclusions. Descriptor
availability is audited on the complete candidate table; imputation and scales
are fitted only on the frozen background. No per-candidate dimensions are dropped.

Scales use `mad(X,0,1)`, the mean absolute deviation about the mean.
Name-matched scale multipliers are 0.2; nonfinite scales or scales below 1e-4
are excluded. Coordinates are `X./scale`, without centering. Nonfinite retained
values are replaced with background medians. The medoid is fixed in the full
retained reference space. P is the mean non-medoid background-relative proximity
percentile (0–100); J = −RFP. All original research parameters are unchanged.

## Run in order

Run the five root Step files in MATLAB. They create `results/research/`:

| Step | Reads | Creates under results/research/ |
|---|---|---|
| 00 | Four input groups above | cache/dataset_metadata.mat, dataset_partition.mat, screening_holdout.tsv, descriptor_missingness_dataset.tsv; audit tables |
| 01 | Step 00 outputs and raw descriptors | cache/data01_model.mat, data02_background.mat, preprocess_metadata.mat |
| 02 | Step 01 matrices and metadata | models.mat, optimization histories and summaries; 64 new runs × 3,000 generations |
| 03 | Steps 00–02 outputs and raw candidate descriptors | candidate_ranking.tsv/.mat and rank-stability statistics |
| 04 | Current Step 02/03 outputs | figures/FIG3A–C and FIG4A–E, as PNG and PDF |

Intermediate files are generated locally, not required as initial inputs.
To repeat the workflow without overwriting an existing run, select a new
`cfg.Output.Directory` in `functions/pso_mdl_config.m`. A relocated project
uses its own local inputs and outputs; no original workstation path is required.
