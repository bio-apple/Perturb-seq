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
| `replicate_col` | `--replicate-col` | string or `null`; obs column for biological replicates |
| `singlet_only` | `--keep-multiplets` | YAML `true`/`false`. CLI flag **sets `singlet_only: false`** when present (no “force singlet” flag) |
| `n_mads` | `--n-mads` | float; MAD QC threshold |
| `n_perms` | `--n-perms` | int; E-test permutations |
| `skip_mixscape` | `--skip-mixscape` | bool / `store_true` |
| `force_mixscape` | `--force-mixscape` | bool / `store_true`; run Mixscape even if target count &gt; `mixscape_max_targets` |
| `mixscape_max_targets` | `--mixscape-max-targets` | int; skip Mixscape above this unless forced |
| `de_top_n` | `--de-top-n` | int |
| `n_jobs` | `--n-jobs` | int; `1` = sequential, `-1` = all CPUs (DE / E-test outer loops) |
| `skip_cell_annotation` | `--skip-cell-annotation` | bool / `store_true` |
| `skip_distance` | `--skip-distance` | bool / `store_true` |
| `skip_de` | `--skip-de` | bool / `store_true` |
| `perturbation_type` | `--perturbation-type` | string; Mixscape label e.g. `KO` / `KD` |
| `control_patterns` | `--control-patterns` | YAML list of regexes; CLI = **comma-separated** regexes |
| `random_state` | `--random-state` | int; PCA / neighbors / UMAP seed |
| `min_cells` | *(YAML only)* | int; gene filter min cells (QC) |
| `n_top_genes` | *(YAML only)* | int; HVG count |
| `n_pcs` | *(YAML only)* | int |
| `leiden_resolution` | *(YAML only)* | float |
| `min_cells_per_pert` | *(YAML only)* | int; min cells per perturbation for distance / DE grouping |

Unknown YAML keys are stored on `PipelineConfig.extra` (not validated as pipeline fields).

## `run` — CLI-only (not PipelineConfig / not in default.yaml)

| CLI flag | Effect |
| --- | --- |
| `--config` | Path to YAML; default search `configs/default.yaml` |
| `--resume` | Skip stages that already succeeded with unchanged inputs/params |

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
