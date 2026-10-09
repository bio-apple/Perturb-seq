# CLI ↔ YAML mapping

Precedence for `python -m perturbseq run`: **CLI overrides YAML**. If `--config` is omitted and `configs/default.yaml` exists, that file is loaded automatically (`repro.default_config_path`).

Sources of truth: `configs/default.yaml`, `PipelineConfig` in `src/perturbseq/pipeline.py`, `_cli_overrides` in `src/perturbseq/cli.py`.

## `run` — PipelineConfig keys

| YAML key (`configs/default.yaml`) | CLI flag | Type / notes |
| --- | --- | --- |
| `input_dir` | `--input-dir` | Path; **required** (YAML and/or CLI) |
| `output_dir` | `--output-dir` | Path; **required** (YAML and/or CLI) |
| `sample_id` | `--sample-id` | string; default `sample1` |
| `control` | `--control` | string; default `NT` |
| `replicate_col` | `--replicate-col` | string or `null`; obs column for biological replicates (≥2 → inferential PyDESeq2) |
| `de_prefer_pseudobulk` | `--de-prefer-pseudobulk` | bool; default `true` — prefer sample_id/pseudo-rep pseudobulk even without bio reps |
| `n_pseudo_replicates` | `--pseudo-replicates` | int or `null`; N≥2 technical pseudo-replicates for single-sample exploratory PyDESeq2 |
| `singlet_only` | `--keep-multiplets` | YAML `true`/`false`. CLI flag **sets `singlet_only: false`** when present (no “force singlet” flag) |
| `n_mads` | `--n-mads` | float; MAD QC threshold |
| `perturbation_aware_qc` | `--perturbation-aware-qc` | bool / `store_true`; fit MAD on NT/control only, apply to all (preferred for Perturb-seq) |
| `n_perms` | `--n-perms` | int; E-test permutations |
| `min_cells_per_pert` | `--min-cells-per-pert` | int; min cells per perturbation for distance / DE grouping |
| `etest_power_min_cells` | `--etest-power-min-cells` | int; default 50; below → `low_power` / suppress `significant_adj_reported` |
| `secondary_distance_metrics` | `--secondary-distance-metrics` | YAML list or CLI comma-separated (e.g. `mmd,wasserstein`); graceful skip if deps missing |
| `n_bootstrap` | `--n-bootstrap` | int; default 100; cell bootstrap for E-distance CI (`0` = skip) |
| `skip_mixscape` | `--skip-mixscape` | bool / `store_true` |
| `force_mixscape` | `--force-mixscape` | bool / `store_true`; run Mixscape even if target count &gt; `mixscape_max_targets` |
| `mixscape_max_targets` | `--mixscape-max-targets` | int; auto-skip Mixscape above this unless forced / subset |
| `mixscape_mode` | `--mixscape-mode` | `auto` \| `skip` \| `force` \| `subset` |
| `mixscape_targets` | `--mixscape-targets` | YAML list or CLI comma-separated genes; implies subset |
| `mixscape_top_n` | `--mixscape-top-n` | int; top-N by pre-Mixscape E-distance; implies subset |
| `de_top_n` | `--de-top-n` | int |
| `report_plot_top_n` | `--report-plot-top-n` | int; default 15; cap guide-consistency / target-validation PNGs |
| `n_jobs` | `--n-jobs` | int; `1` = sequential, `-1` = all CPUs (DE / E-test outer loops) |
| `skip_cell_annotation` | `--skip-cell-annotation` | bool / `store_true` |
| `skip_distance` | `--skip-distance` | bool / `store_true` |
| `skip_de` | `--skip-de` | bool / `store_true` |
| `perturbation_type` | `--perturbation-type` | string; Mixscape label e.g. `KO` / `KD` |
| `perturbation_space` | `--perturbation-space` | `pca_silhouette` (default) \| `kmeans` \| `lr_classifier` — perturbation embedding + clusters |
| `control_patterns` | `--control-patterns` | YAML list of regexes; CLI = **comma-separated** regexes |
| `random_state` | `--random-state` | int; PCA / neighbors / UMAP seed |
| `min_cells` | *(YAML only)* | int; gene filter min cells (QC) |
| `n_top_genes` | *(YAML only)* | int; HVG count |
| `n_pcs` | *(YAML only)* | int |
| `leiden_resolution` | *(YAML only)* | float |
| `guide_merge` | `--guide-merge` | `none` \| `equal` \| `umi` \| `confidence` \| `umi_confidence` |
| `on_target_lfc_cutoff` | `--on-target-lfc-cutoff` | float; Guide QC on-target \|log2FC\| cutoff (default `0.25`) |
| `on_target_min_fail_guides` | `--on-target-min-fail-guides` | int; ≥N failing adequate guides → `potential_low_efficiency` (default `2`) |
| `guide_reassign` | `--guide-reassign` | `off` \| `compare` \| `apply_max` \| `apply_gmm` |
| `guide_reassign_min_umi` | `--guide-reassign-min-umi` | float; min CRISPR UMI for max/GMM |
| `de_covariates` | `--de-covariates` | `true` / `false` / YAML list / CLI comma-separated obs cols |

Unknown YAML keys are stored on `PipelineConfig.extra` (not validated as pipeline fields).

## `run` — CLI-only (not PipelineConfig / not in default.yaml)

| CLI flag | Effect |
| --- | --- |
| `--config` | Path to YAML; default search `configs/default.yaml` |
| `--resume` | Skip stages that already succeeded with unchanged inputs/params |
| `--dry-run` | Validate DRAGEN inputs; print planned stages + resolved params; **no analysis** |

## `config` — generate a commented template

Writes a fully commented YAML covering `configs/default.yaml` values plus all `PipelineConfig` fields:

```bash
python -m perturbseq config --generate
python -m perturbseq config generate -o my_pipeline.yaml
# equivalent: perturbseq config --generate -o configs/pipeline.template.yaml
```

## Other subcommands (no YAML merge)

These do **not** use `configs/default.yaml` / `PipelineConfig`.

### `write-demo`

| CLI | Default |
| --- | --- |
| `--output-dir` | `data/demo` |

### `annotate`

| CLI | Notes |
| --- | --- |
| `--h5ad` | required |
| `--output-dir` | required |
| `--inplace` | overwrite input h5ad |

### `report`

| CLI | Notes |
| --- | --- |
| `--output-dir` | required; existing results dir |
| `--control` | optional override (else from `report.json`) |
| `--min-cells` | optional; report helper, **not** the same as YAML `min_cells` QC |
| `--html-detail-limit` | optional cap on detailed HTML sections |
| `--h5ad` | optional; expression for guide-consistency / target-validation PNGs (auto-discovers `*.tertiary.h5ad`) |
| `--report-plot-top-n` | optional; max genes for those PNGs (YAML `report_plot_top_n`, default 15) |

### `guide-qc`

| CLI | Default |
| --- | --- |
| `--h5ad` | required |
| `--output-dir` | required |
| `--control` | `NT` |
| `--min-cells` | `10` |
| `--min-median-umi` | `5.0` |
| `--min-detection-rate` | `0.5` |

## Example: same config two ways

YAML (`configs/default.yaml` excerpt):

```yaml
skip_mixscape: true
n_jobs: 4
control_patterns:
  - "^nt$"
  - "^negctrl"
```

Equivalent CLI overrides on top of defaults:

```bash
python -m perturbseq run \
  --input-dir data/demo \
  --output-dir results/demo \
  --skip-mixscape \
  --n-jobs 4 \
  --control-patterns '^nt$,^negctrl'
```

Dry-run (validate + plan only):

```bash
python -m perturbseq run \
  --input-dir data/demo \
  --output-dir results/demo \
  --skip-mixscape \
  --dry-run
```
