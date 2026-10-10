# Perturb-seq tertiary analysis

[English](README.md) | [简体中文](README.zh-CN.md)

Tertiary analysis on **DRAGEN scRNA CRISPR / Perturb-seq secondary outputs**, using **scanpy + pertpy** (default Python stack). The [Seurat Mixscape vignette](https://satijalab.org/seurat/articles/mixscape_vignette) is conceptual parity / an optional R path — see [MIXSCAPE_SEURAT_MAP.md](docs/MIXSCAPE_SEURAT_MAP.md). Aligned with [sc-best-practices · perturbation modeling](https://www.sc-best-practices.org/conditions/perturbation-modeling/), [pertpy](https://pertpy.scverse.org/en/latest/index.html) ([Heumos et al. 2026](https://www.nature.com/articles/s41592-025-02909-7)), and the [Illumina Connected Multiomics Perturb-seq walkthrough](https://help.multiomics.illumina.com/icm/analyses/walkthroughs/perturb-seq).

**Docs:** [tutorial notebook](docs/tutorial.ipynb) · [tutorial.md](docs/tutorial.md) · [CLI ↔ YAML](docs/CLI_YAML.md) · [Statistical caveats](docs/STATISTICAL_CAVEATS.md) · [Step dependencies](docs/ANALYSIS_DEPENDENCIES.md) · [sc-best-practices map + module checklist](docs/SC_BEST_PRACTICES_MAP.md#tertiary-module-coverage-checklist) · [Mixscape ↔ Seurat vignette](docs/MIXSCAPE_SEURAT_MAP.md)

> Chinese readers: [README.zh-CN.md](README.zh-CN.md) covers the same install / run / outputs / caveats in 中文 (full detail). This English file is the maintainable primary README for packaging (`pyproject.toml` → `readme`).

## Inputs (DRAGEN secondary)

Per sample (e.g. `sample1`):

| File | Role |
| --- | --- |
| `{sample}.scRNA.filtered.matrix.mtx.gz` | Filtered cell × feature counts (genes + CRISPR) |
| `{sample}.scRNA.filtered.barcodes.tsv.gz` | Cell barcodes |
| `{sample}.scRNA.filtered.features.tsv.gz` | Gene/gRNA IDs, names, feature types |
| `{sample}.scRNA.feature_barcode_reference.csv` | gRNA reference |
| `{sample}.scRNA.positive_cell_guide_assignments.csv` | DRAGEN GMM guide calling |

Field definitions: [DRAGEN v4.5 CRISPR mode](https://help.dragen.illumina.com/dragen-v4.5/product-guides/dragen-v4.5/dragen-single-cell-pipeline/dragen-scrna-illumina#crispr-mode). Place files in a directory such as `data/raw/`.

## Pipeline (look at the data first)

```
DRAGEN MEX + guide assignment
  → split Gene Expression / CRISPR Direct Capture
  → annotate guide_id, gene_target, NT vs perturbed
  → QC (n_counts / n_genes / %MT, MAD; prefer --perturbation-aware-qc)  [filtering]
  → default: keep num_features == 1 singlets               [filtering]
  → normalize / log1p (keep layers['counts']) + HVG/PCA/Leiden  [descriptive]
  → UMAP                                                   [visualization only]
  → cell annotation: cell cycle + state scores             [descriptive; skipped by default for sample_type=cell_line]
  → Mixscape: KO/NP; optional KO+NT subset                 [filtering/classification]
  → E-distance / E-test (pre/post Mixscape)                [inferential]
  → DE: prefer pseudobulk + PyDESeq2 (bio reps → inferential; else exploratory); Wilcoxon only if <2 pseudobulk units
  → perturbation space (`--perturbation-space`: pca_silhouette | kmeans | lr_classifier)
  → Guide QC + report plots (guide-consistency / target-validation PNGs)
```

Controls default to patterns like `NegCtrl*`, `NT`/`NTC`, `non-targeting` → `NT`. Override with `--control-patterns`.

## Install

**Recommended** (includes `pertpy[de]`: Mixscape / E-distance / E-test / PyDESeq2):

```bash
pip install -e ".[de,dev]" -c constraints.txt
# or
conda env create -f environment.yml
conda activate perturbseq-tertiary
```

**Core only** (no pertpy): `pip install -e ".[dev]" -c constraints.txt` — then always pass `--skip-mixscape --skip-distance`.

Entry points: `python -m perturbseq …` or `perturbseq` / `perturbseq-tertiary`.

### Tested versions / Compatibility

This pipeline was end-to-end tested with **pertpy==1.4.0** and **scanpy==1.12.4** (conda env `perturbseq-tertiary`, Python 3.12.15; versions from `importlib.metadata`). Supported ranges below match `pyproject.toml` / `constraints.txt` / `environment.yml` — prefer those files over a drifted local env when debugging installs.

| Package | Tested (E2E) | Supported range |
| --- | --- | --- |
| Python | 3.12.15 | `>=3.10` |
| scanpy | 1.12.4 | `>=1.10,<2` |
| pertpy | 1.4.0 | `>=1.3,<2` (`pertpy` / `de` / `all` extras) |
| anndata | 0.13.4 | `>=0.10,<0.14` |

### Without pertpy

| Capability | `.[de]` | Core only |
| --- | --- | --- |
| Input validation, guide annotate, QC/MAD, singlet | ✅ | ✅ |
| Normalize / HVG / PCA / Leiden / UMAP, cell annotation | ✅ | ✅ |
| Guide QC | ✅ | ✅ |
| Exploratory DE (Wilcoxon fallback / no-bio-rep pseudobulk) | ✅ (PyDESeq2 needs pertpy) | ✅ Wilcoxon only |
| Mixscape / E-distance / E-test / PyDESeq2 / pert. clusters | ✅ | ❌ → explicit `skipped` + `reason` in `report.json` |

## Run

Interactive walkthrough (synthetic demo, offline): **[docs/tutorial.ipynb](docs/tutorial.ipynb)**.

```bash
jupyter notebook docs/tutorial.ipynb
# or: jupyter lab docs/tutorial.ipynb
```

CLI smoke test:

```bash
python -m perturbseq write-demo --output-dir data/demo
python -m perturbseq run --input-dir data/demo --output-dir results/demo --sample-id sample1
```

YAML + resume (flag map: [docs/CLI_YAML.md](docs/CLI_YAML.md)):

```bash
python -m perturbseq run \
  --config configs/default.yaml \
  --input-dir data/demo \
  --output-dir results/demo \
  --sample-id sample1 \
  --resume
```

Real data: `--input-dir data/raw --output-dir results/sample1 --sample-id sample1 --control NT`.

Useful flags (full map: [CLI_YAML.md](docs/CLI_YAML.md)): `--config`, `--resume`, `--dry-run`, `--random-state`, `--n-jobs`, `--keep-multiplets`, `--perturbation-aware-qc`, `--replicate-col`, `--de-prefer-pseudobulk`, `--pseudo-replicates`, `--n-bootstrap` (E-distance CI; `0` = skip), `--etest-power-min-cells`, `--sample-type` (`cell_line` default → skip full cell annotation; `primary`/`mixed` → run; override with `--skip-cell-annotation` / `--run-cell-annotation`), `--skip-mixscape` / `--skip-distance` / `--skip-de`, `--mixscape-mode` / `--mixscape-targets` / `--mixscape-top-n`, `--perturbation-type`, `--perturbation-space`, `--guide-merge`, `--on-target-lfc-cutoff` / `--on-target-min-fail-guides`, `--guide-reassign`, `--de-covariates`, `--report-plot-top-n`. KD/CRISPRi caveats: [STATISTICAL_CAVEATS](docs/STATISTICAL_CAVEATS.md).

Standalone: `annotate`, `guide-qc`, `report` (rebuild HTML/JSON + optional report PNGs without re-running analysis).

## Outputs

```
results/sample1/
  sample1.tertiary.h5ad
  report.json / report.html          # verdict + per-perturbation cards; evidence_level on DE
  run_manifest.json / run.log
  stages/<stage>/status.json         # inputs/params_hash for --resume
  figures/                           # qc_violin/scatter/cell_counts, guide_*, umap*, composition_*, volcano_*, de_heatmap, …
  tables/
    composition_audit.csv
    composition_by_perturbation.csv  # cluster × gene_target counts/fractions
    guide_assignment_summary.csv
    qc_cell_counts.csv
    qc_filter_by_guide.csv
    edistance.csv, etest.csv, distances.csv, optional distance_mmd.csv, de_*.csv
    guide_qc.csv, gene_guide_consistency.csv, qc_warnings.csv
    gene_guide_weighted.csv          # when guide_merge ≠ none
    perturbation_clusters.csv        # (+ optional perturbation_space_embeddings.csv)
    pre_mixscape/ / post_mixscape/   # if Mixscape ran
```

Stages: `1_input_validation` → `2_preprocessing_qc` → `3a_perturbation_modeling` → `3b_statistical_inference` → `4_robustness` → `5_report`. Module checklist + book-chapter map: [SC_BEST_PRACTICES_MAP.md](docs/SC_BEST_PRACTICES_MAP.md#tertiary-module-coverage-checklist). Matrix layers / pre–post Mixscape: [ANALYSIS_DEPENDENCIES.md](docs/ANALYSIS_DEPENDENCIES.md).

Key `obs`: `guide_id`, `gene_target`, `perturbation`, `num_features`, `guide_umi`, `mixscape_class*`. `layers['counts']` = raw UMI; `X` = log-norm. Guide QC layers + on-target / `potential_low_efficiency` / optional `guide_merge` / `guide_reassign`: see [STATISTICAL_CAVEATS](docs/STATISTICAL_CAVEATS.md) (not re-expanded here).

## 常见坑 / Common pitfalls

Full discussion + citations: **[docs/STATISTICAL_CAVEATS.md](docs/STATISTICAL_CAVEATS.md)** · also in [tutorial](docs/tutorial.md#常见坑--common-pitfalls).

| Pitfall | Short why |
| --- | --- |
| Multi-guide inconsistency | Same gene, discordant guide effects — do not pool into one “gene KO” without review (`gene_guide_consistency.csv`) |
| Low on-target / multi-guide fail | ≥2 adequate guides fail expected target log2FC (`potential_low_efficiency`) — do **not** conclude “no phenotype” |
| Mixscape labels many NP | Transcriptomic resemblance to NT ≠ proven failed edit |
| No FDR claims without replicates | Cells ≠ biological replicates; no-bio-rep DE (Wilcoxon or exploratory PyDESeq2) is exploratory only (Squair 2021) |
| UMAP/Leiden as “effect” | Global structure is not perturbation evidence |
| Filtering changes hits | QC / singlet / Mixscape KO drop cells — check `composition_audit.csv`; prefer `--perturbation-aware-qc`; see `qc_filter_by_guide.csv` |
| Inherited DRAGEN assignment | This pipeline does not re-call guides by default |

## Method choices (why these)

| Step | Choice | Not default |
| --- | --- | --- |
| Guide assignment | DRAGEN GMM | Re-calling (e.g. crispat) |
| Non-perturbed cells | Mixscape (Papalexi 2021; stronger for KO than KD/CRISPRi); [Seurat vignette map](docs/MIXSCAPE_SEURAT_MAP.md) | Treat all targeting cells as KO |
| Effect size | E-distance / E-test (Peidli 2024) | UMAP separation alone |
| DE | Pseudobulk + PyDESeq2 first (bio reps → inferential; else exploratory); Wilcoxon only as fallback | Silent cell-level Wilcoxon as “the” DE / population FDR |
| Unseen perturbation prediction | Not a default step | scGen / foundation models as primary |

## Project layout

```
src/perturbseq/   # io, guides, qc, preprocessing, perturbation, statistics, …
configs/default.yaml
docs/tutorial.ipynb   # interactive demo (offline)
docs/tutorial.md      # CLI twin + pitfalls
docs/STATISTICAL_CAVEATS.md
tests/
```

Lint (core API):

```bash
pip install -e ".[dev]" -c constraints.txt
ruff check src tests
python -m mypy src/perturbseq/qc.py src/perturbseq/perturbation.py src/perturbseq/statistics.py src/perturbseq/_deps.py
```
