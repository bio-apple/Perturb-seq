# Tutorial: demo data → HTML report

End-to-end walkthrough for new users. Structure follows [sc-best-practices · perturbation modeling](https://www.sc-best-practices.org/conditions/perturbation-modeling/): motivation → steps → inspection → caveats.

Related: [CLI ↔ YAML](CLI_YAML.md) · [Statistical caveats](STATISTICAL_CAVEATS.md) · [Step dependencies](ANALYSIS_DEPENDENCIES.md)

## Motivation

You have DRAGEN CRISPR-mode secondary outputs (MEX + guide assignments) and want a tertiary pass: QC, Mixscape, E-distance, DE, and a readable report. This tutorial uses **synthetic** DRAGEN-like files so you can verify the install without real sequencing data.

Demo data is small (~120 cells, few guides). Treat results as a **pipeline smoke test**, not biology.

## Prerequisites

```bash
conda env create -f environment.yml
conda activate perturbseq-tertiary
# from repo root, with package importable (editable install or PYTHONPATH=src)
```

## Step 1 — Write demo inputs

```bash
python -m perturbseq write-demo --output-dir data/demo
```

**Why:** Creates a minimal DRAGEN-shaped bundle (`sample1.scRNA.*`) with NegCtrl + targeting guides, singletons, and a few multiplets / empty assignments. Same file names the real pipeline expects.

**Inspect:**

```bash
ls data/demo/
# sample1.scRNA.filtered.matrix.mtx.gz
# sample1.scRNA.filtered.barcodes.tsv.gz
# sample1.scRNA.filtered.features.tsv.gz
# sample1.scRNA.feature_barcode_reference.csv
# sample1.scRNA.positive_cell_guide_assignments.csv
```

## Step 2 — Run the tertiary pipeline

```bash
python -m perturbseq run \
  --config configs/default.yaml \
  --input-dir data/demo \
  --output-dir results/demo \
  --sample-id sample1
```

**Parameter rationale (demo defaults):**

| Choice | Default | Why here |
| --- | --- | --- |
| `control: NT` | YAML / `--control` | Demo NegCtrl guides map to `NT` via `control_patterns` |
| `singlet_only: true` | (omit `--keep-multiplets`) | Cleaner assignment for a first run |
| `n_mads: 5` | YAML | Adaptive QC; demo is tiny so few cells drop |
| `n_jobs: 1` | YAML | Reproducible sequential DE / E-test |
| `random_state: 0` | YAML | Stable PCA / neighbors / UMAP |
| Mixscape on | (omit `--skip-mixscape`) | Demo has few targets (&lt; `mixscape_max_targets`) so Mixscape runs |
| `n_perms: 200` | YAML | E-test permutations; raise for publication runs |

Faster smoke test (skip Mixscape + distance):

```bash
python -m perturbseq run \
  --config configs/default.yaml \
  --input-dir data/demo \
  --output-dir results/demo_fast \
  --sample-id sample1 \
  --skip-mixscape --skip-distance
```

Resume after interrupt (same inputs/params → skip finished stages):

```bash
python -m perturbseq run \
  --config configs/default.yaml \
  --input-dir data/demo \
  --output-dir results/demo \
  --sample-id sample1 \
  --resume
```

## Step 3 — Inspect intermediate outputs

```
results/demo/
  sample1.tertiary.h5ad
  report.html / report.json
  run_manifest.json / run.log
  stages/<stage>/status.json   # + optional adata checkpoints
  figures/
  tables/
```

**Suggested order:**

1. **`run_manifest.json` / `run.log`** — Did every stage succeed? Check versions and checksums.
2. **`tables/composition_audit.csv`** — How many cells remain after QC / singlet / Mixscape filters.
3. **`tables/edistance.csv`** (and `pre_mixscape/` / `post_mixscape/` if Mixscape ran) — Effect sizes vs NT.
4. **`tables/etest.csv`** — Permutation p-values / padj for distances.
5. **`tables/de_*.csv`** — Without replicates, expect exploratory Wilcoxon (see caveats).
6. **`tables/guide_qc.csv`**, **`gene_guide_consistency.csv`**, **`qc_warnings.csv`** — Guide-level QC.
7. **`figures/`** — QC violin/scatter, UMAP (viz only), volcano, E-distance plots.
8. **`sample1.tertiary.h5ad`** — `obs`: `guide_id`, `gene_target`, `perturbation`, `mixscape_class*`; `layers['counts']` raw; `X` log-norm.

Quick Python peek:

```python
import scanpy as sc
adata = sc.read_h5ad("results/demo/sample1.tertiary.h5ad")
print(adata)
print(adata.obs[["guide_id", "gene_target", "perturbation"]].head())
```

## Step 4 — Open the report

```bash
open results/demo/report.html   # macOS; or open in any browser
```

Or rebuild HTML/JSON from an existing results directory without re-running analysis:

```bash
python -m perturbseq report --output-dir results/demo
```

`report.html` is organized around questions per perturbation (cell counts, effect size, consistency, top genes, warnings). Prefer CSV / h5ad / `report.json` for downstream scripting.

## Step 5 — Optional follow-ons

```bash
# Cell-cycle / state annotation only
python -m perturbseq annotate \
  --h5ad results/demo/sample1.tertiary.h5ad \
  --output-dir results/demo \
  --inplace

# Guide QC only
python -m perturbseq guide-qc \
  --h5ad results/demo/sample1.tertiary.h5ad \
  --output-dir results/demo \
  --control NT
```

## Caveats (demo and beyond)

- **Cells ≠ biological replicates.** Demo has no `replicate_col` → Wilcoxon DE is exploratory.
- **UMAP / Leiden** describe structure; they are not perturbation evidence.
- **Mixscape NP** ≠ proven biological null; see [STATISTICAL_CAVEATS.md](STATISTICAL_CAVEATS.md).
- Demo NegCtrl / targeting modules are synthetic — do not over-interpret gene names.

When ready for real data, point `--input-dir` at a DRAGEN sample folder and keep the same inspection order. Map all flags in [CLI_YAML.md](CLI_YAML.md).
