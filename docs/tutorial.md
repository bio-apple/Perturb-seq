# Tutorial: demo data → HTML report

End-to-end walkthrough for new users. Structure follows [sc-best-practices · perturbation modeling](https://www.sc-best-practices.org/conditions/perturbation-modeling/): motivation → steps → inspection → caveats.

**Interactive notebook (recommended):** [tutorial.ipynb](tutorial.ipynb) — synthetic demo runs offline; optional Papalexi 2021 cell needs network + pertpy.

```bash
# from repo root, after install
jupyter notebook docs/tutorial.ipynb
# or: jupyter lab docs/tutorial.ipynb
```

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

**Why:** Creates a minimal DRAGEN-shaped bundle (`sample1.scRNA.*`) with NegCtrl + targeting guides, singletons, and a few multiplets / empty assignments. Same file names the real pipeline expects. The demo also plants **IFNGR2** with two guides that have **conflicting target-gene effects** so guide QC can flag multi-guide inconsistency (see Common pitfalls below).

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
| `de_prefer_pseudobulk: true` | YAML | Prefer PyDESeq2 when ≥2 `sample_id` or `--pseudo-replicates`; single-sample demo → Wilcoxon fallback |
| `n_jobs: 1` | YAML | Reproducible sequential DE / E-test |
| `random_state: 0` | YAML | Stable PCA / neighbors / UMAP |
| Mixscape on | (omit `--skip-mixscape`) | Demo has few targets (&lt; `mixscape_max_targets`) so Mixscape runs |
| `n_perms: 200` | YAML | E-test permutations; raise for publication runs |
| `n_bootstrap: 100` | YAML / `--n-bootstrap` | E-distance percentile CI (`0` = skip) |
| `etest_power_min_cells: 50` | YAML / `--etest-power-min-cells` | Below → `low_power`; Peidli ~50–100 safer, ~200 more stable |
| `secondary_distance_metrics` | YAML / CLI | Default `mmd,wasserstein`; skip gracefully if JAX missing |
| `perturbation_space` | `pca_silhouette` | Descriptive clusters; also `kmeans` / `lr_classifier` |
| `report_plot_top_n: 15` | YAML | Cap guide-consistency / target-validation PNGs |

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

1. **`run_manifest.json` / `run.log`** — Did every stage succeed? Check versions and checksums (`params_hash` for `--resume`).
2. **`tables/composition_audit.csv`** — How many cells remain after QC / singlet / Mixscape filters.
3. **`tables/edistance.csv`** / **`distances.csv`** (and `pre_mixscape/` / `post_mixscape/` if Mixscape ran) — Effect sizes vs NT (+ bootstrap CI columns when `n_bootstrap` &gt; 0).
4. **`tables/etest.csv`** — Permutation p-values / padj; check `low_power` and `significant_adj_reported`.
5. **`tables/de_*.csv`** — Check `evidence_level` / method: single-sample demo → exploratory Wilcoxon; ≥2 samples or `--pseudo-replicates` → exploratory PyDESeq2; true bio reps → inferential (see caveats).
6. **`tables/guide_qc.csv`**, **`gene_guide_consistency.csv`**, **`qc_warnings.csv`** — Guide-level QC. Optional **`gene_guide_weighted.csv`** when `guide_merge` ≠ `none`.
7. **`figures/`** — QC, UMAP (viz only), volcano, E-distance, plus `guide_consistency_*` / `target_validation_*` report PNGs.
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

`report.html` opens with an **Analysis verdict** card (one sentence from guide QC + E-test + DE — also in `report.json` as `verdict` / `summary_sentence`). Below that, each perturbation answers cell counts, effect size, consistency, top genes, warnings. Cards surface guide `interpretation` / `qc_warnings` (and optional weighted summary when `guide_merge` ≠ `none`). Prefer CSV / h5ad / `report.json` for downstream scripting. Skipped optional steps appear as objects like `"mixscape": {"skipped": true, "reason": "user_skip", ...}` — not missing keys or bare `null`.

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

## Optional — public dataset (Papalexi 2021)

Needs **network** (first download) and `pertpy`. Not DRAGEN MEX — useful for comparing with Mixscape literature. Full interactive cell: [tutorial.ipynb](tutorial.ipynb).

```python
# optional
import pertpy as pt
mdata = pt.data.papalexi_2021()
print(mdata)
```

## 常见坑 / Common pitfalls

Full text + literature: **[STATISTICAL_CAVEATS.md](STATISTICAL_CAVEATS.md)**.

| Pitfall | Why | What to do |
| --- | --- | --- |
| **Multi-guide inconsistency** | Same gene, different guides → discordant directions (assignment, UMI, incomplete edit, off-target, cytotoxicity) | Inspect `gene_guide_consistency.csv` / `qc_warnings.csv`; do **not** pool inconsistent guides into one “gene KO” claim |
| **Low on-target / multi-guide fail** | ≥2 adequate guides fail target log2FC in expected direction (`potential_low_efficiency`) | Flag low-efficiency guides; do **not** conclude “no phenotype” for the gene |
| **Mixscape labels many NP** | Targeting barcode but transcriptome still resembles NT; pool-dependent; **not** proof the guide failed | Read composition audit; NP ≠ proven biological null |
| **No FDR claims without replicates** | Cells ≠ biological replicates; no-bio-rep DE (exploratory PyDESeq2 or Wilcoxon fallback) is not population FDR | Hypothesis ranking only; ≥2 true replicates → `evidence_level=inferential` |
| **UMAP / Leiden as effect** | Global structure ≠ perturbation evidence | Prefer E-distance / E-test (+ replicate-aware DE) |
| **Perturbation-space clusters as mechanism** | `pca_silhouette` / `kmeans` / `lr_classifier` are descriptive embeddings | Compare methods via `--perturbation-space`; do not equate co-cluster with shared pathway |
| **Filtering changes composition** | MAD QC, singlet, Mixscape KO+NT drop cells | Compare `tables/composition_audit.csv`; prefer `--perturbation-aware-qc` so strong phenotypes do not set MAD cutoffs; check `qc_filter_by_guide.csv` / `cytotoxicity_qc_depletion` |
| **Inherited guide assignment** | Tertiary pipeline does not re-call guides by default | Assignment quality bounds everything downstream |

### Hands-on: multi-guide inconsistency (demo)

Synthetic `write-demo` plants **IFNGR2** with two guides: `IFNGR2_1` knocks down the target gene; `IFNGR2_2` has a conflicting (upregulated) target effect. After Step 2:

```bash
python - <<'PY'
import pandas as pd
c = pd.read_csv("results/demo/tables/gene_guide_consistency.csv")
print(c.loc[c["gene_target"] == "IFNGR2", ["gene_target", "directions", "guides_consistent"]])
w = pd.read_csv("results/demo/tables/qc_warnings.csv")
print(w.loc[w["warning"] == "inconsistent_guides"])
PY
```

Expect `guides_consistent == False` and an `inconsistent_guides` warning. The HTML verdict / IFNGR2 card should surface the same flag. Do **not** pool those guides into one KO claim.

Demo NegCtrl / targeting modules are synthetic — do not over-interpret gene names.

When ready for real data, point `--input-dir` at a DRAGEN sample folder and keep the same inspection order. Map all flags in [CLI_YAML.md](CLI_YAML.md).
