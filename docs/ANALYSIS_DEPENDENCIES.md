# Analysis step dependencies

This document separates **descriptive** structure discovery from steps that **filter cells** or change **inference scope**. Global HVG / PCA / neighbors / Leiden / UMAP are **not** perturbation-effect evidence by default.

Module map: `preprocessing.py` (normalize/PCA/UMAP), `perturbation.py` (Mixscape/E-distance), `statistics.py` (DE), `composition.py` (filter composition audit + cluster×perturbation), orchestrated by `pipeline.py`. See also `docs/structure.md`. Pipeline step ↔ sc-best-practices chapter sections + biologist module checklist (DRAGEN / ICM / sc-best-practices): [SC_BEST_PRACTICES_MAP.md](SC_BEST_PRACTICES_MAP.md#tertiary-module-coverage-checklist).

## Matrix layers

| Slot | Contents | Typical use |
| --- | --- | --- |
| `layers['counts']` | Raw UMI counts (post QC/singlet gene filter) | Pseudobulk DE (PyDESeq2); never discarded |
| `adata.X` | `normalize_total` + `log1p` | Wilcoxon exploratory DE; HVG/PCA input |
| `layers['X_pert']` | Mixscape perturbation signature (if run) | Mixscape classification; optional post-Mixscape PCA |
| `obsm['X_pca']` | Low-dimensional embedding | E-distance / clustering input; **not** a count matrix |
| `obsm['X_umap']` | 2D layout for plots | **Visualization only** |

## Dependency table

| Stage | Kind | Input matrix / labels | Output | Changes inference scope? | Caveats |
| --- | --- | --- | --- | --- | --- |
| Load + guide annotate | descriptive | DRAGEN MEX + assignments | `guide_id`, `gene_target`, `perturbation` | No | Composition already reflects DRAGEN calling |
| Guide reassignment (optional) | QC / optional filtering | CRISPR count matrix + DRAGEN calls | comparison CSV; optional override | **Yes** if `apply_*` | Default `off`; CatchR / Cell Ranger FB not installed |
| QC (MAD + min cells/gene; optional perturbation-aware) | filtering | raw `X` (+ NT mask when aware) | filtered AnnData + QC log + `qc_filter_by_guide.csv` | **Yes** — drops cells/genes | Global MAD can drop strong phenotypes; prefer `perturbation_aware_qc` |
| Singlet filter | filtering | `num_features` | singlet cells | **Yes** — drops 0/≥2 guide cells | Can deplete rare guides or doublets carrying real biology |
| Normalize / log1p | transform | `layers['counts']` → `X` | log-norm expression | No (representation) | Downstream on `X` is not count-scale |
| HVG / PCA / neighbors | descriptive | log-norm `X` (HVG) | `X_pca`, graph | No for perturbation claims | Global structure; confounded by cell cycle, ambient, guide load |
| UMAP | descriptive (viz) | neighbor graph | `X_umap` | No | **Do not** treat separation as KO evidence |
| Leiden | descriptive | neighbor graph | `obs['leiden']` | No | Cluster labels ≠ perturbation classes |
| Cell annotation | descriptive | log-norm scores / markers | phase, cell_state, … | No (unless you later filter on them) | Gated by `sample_type` (default `cell_line` → skip); cell-line “states”, not tissue taxonomy |
| Mixscape | filtering + classification | signature on expression | `mixscape_class*`, `X_pert`; optional `uns['mixscape_lda']` | **Yes** — NP vs KO changes who counts as perturbed | NP ≠ proven null; depends on NT pool size; Seurat vignette map: [MIXSCAPE_SEURAT_MAP.md](MIXSCAPE_SEURAT_MAP.md) |
| Post-Mixscape KO+NT subset | filtering | `mixscape_class_global` | analysis object (default inferential subset) | **Yes** | Before/after in `mixscape_ko_filter` + `tables/mixscape_ko_filter_composition.csv` |
| E-distance / E-test | inferential (effect size) | `X_pca` from `X_pert` when Mixscape ran | `edistance.csv`, `etest.csv`, `distances.csv`, optional `distance_mmd.csv` | Uses KO+NT set | Embedding choice changes ranks; bootstrap CI `edistance_ci_low`/`edistance_ci_high` (`n_bootstrap`, default 100); `low_power` when `n_cells < etest_power_min_cells`; secondary metrics may skip (`missing_jax` / `failed`) |
| DE (Wilcoxon fallback) | exploratory | log-norm `X` | `de_*.csv` | Uses `de_group` | Only when &lt;2 pseudobulk units; cell-level p-values inflate significance; **covariates not modeled**; `evidence_level=exploratory` |
| DE (PyDESeq2, bio reps) | inferential | `layers['counts']` pseudobulk | `de_*.csv` | `de_group` + `replicate_col` + optional covariates | `replicate_col` ≥2 levels; `method=pydeseq2_pseudobulk`; `evidence_level=inferential` |
| DE (PyDESeq2, no bio reps) | exploratory | `layers['counts']` pseudobulk | `de_*.csv` | `de_group` + `sample_id` or technical pseudo-reps | `method=pydeseq2_pseudobulk_no_bio_reps`; not a substitute for bio reps (Squair 2021) |
| Perturbation clustering | descriptive | method-dependent (`perturbation_space`: mean `X_pca` / KMeans / LR coeffs) | `perturbation_clusters.csv`, optional `perturbation_space_embeddings.csv` | Uses current cell set; method in provenance | Pathway-like grouping, not proof of mechanism |
| Guide QC | descriptive / QC flags | log-norm `X` + obs labels (+ optional `qc_filter_by_guide`) | `guide_qc.csv`, `gene_guide_consistency.csv`, `qc_warnings.csv` | No cell drop by itself | On-target proxy / `potential_low_efficiency` ≠ “no phenotype”; optional `guide_merge` is review-only |
| Report plots | visualization | ranked perturbations + expression h5ad | `figures/guide_consistency_*.png`, `figures/target_validation_*.png` | No | Capped by `report_plot_top_n` (default 15); rebuild via `perturbseq report` |

\*Only true multi-rep `replicate_col` yields `evidence_level=inferential`. No-bio-rep pseudobulk and Wilcoxon stay `exploratory`. Machine-readable `evidence_level` is on each DE CSV row, `report.json` → `de` / `experimental_design`, and per-perturbation summaries.

## Pre- vs post-Mixscape paths

When Mixscape **runs successfully**:

| Path | Cells | Embedding for E-distance | Tables |
| --- | --- | --- | --- |
| `pre_mixscape` | All cells after QC/singlet (includes NP) | `X_pca` from log-norm HVG | `tables/pre_mixscape/` |
| `post_mixscape` | `mixscape_class_global` ∈ {control, KO} | `X_pca` recomputed from `layers['X_pert']` | `tables/post_mixscape/` and primary `tables/edistance.csv` |

DE contrasts use Mixscape labels (`mixscape_class`) when Mixscape succeeds — **downstream DE depends on that classification** (KO class vs control; NP cells are not in KO contrasts). The KO+NT filtered object is the default E-distance / clustering subset; composition before vs after NP removal is recorded under `mixscape.ko_filter`.

**Perturbation embedding:** `layers['X_pert']` is the pipeline’s built-in Mixscape signature (see `report.perturbation_embedding`). Optional LDA (`Mixscape.lda` → `uns['mixscape_lda']`, `report.mixscape.lda`) matches the Seurat vignette’s MixscapeLDA viz path — not the primary E-distance embedding. Seurat step map: [MIXSCAPE_SEURAT_MAP.md](MIXSCAPE_SEURAT_MAP.md). SCEPTRE / MIMOSCA / PerturbNet are external advanced options — see [STATISTICAL_CAVEATS.md](STATISTICAL_CAVEATS.md).

When Mixscape is **skipped** or fails, `report.json` records an explicit status object (never a missing key or bare `null`):

```json
"mixscape": {
  "skipped": true,
  "reason": "user_skip|too_many_targets|subset_no_targets|missing_pertpy|failed",
  "detail": "...",
  "n_targets": 785,
  "mixscape_max_targets": 40,
  "mixscape_mode": "auto",
  "estimate": {
    "note": "Heuristic cost ∝ O(n_cells × n_targets); not a wall-clock benchmark. ...",
    "n_cells": 50000,
    "n_targets": 785,
    "approx_work_units": 39250000,
    "approx_relative_to_threshold": 19.62,
    "approx_memory_hint_gb": 15.7
  }
}
"edistance": {"skipped": true, "reason": "user_skip|missing_pertpy|failed", "detail": "..."}
"de": {"skipped": false, "reason": null, "n_contrasts": 3, "scope": "gene_target (...)", ...}
```

Success shape: `{"skipped": false, "reason": null, ...rich fields...}`. Subset success may include `"subset": true` and `"selected_targets": [...]`. E-distance/DE then use `gene_target` on the post-QC object when Mixscape was skipped; no KO/NP filter is applied.

## Composition sensitivity

Each filter (QC, singlet, Mixscape KO-only, optional state/guide filters) can change counts of cells per `gene_target` / guide / state. The pipeline writes `tables/composition_audit.csv` with snapshots at each stage. If a strong E-distance or DE hit disappears after a filter, treat the conclusion as **preprocessing-sensitive**.

## Report language

- UMAP / Leiden: visualization or descriptive structure only.
- Perturbation claims: prefer E-distance / E-test and (when replicates exist) pseudobulk DE, with matrix provenance from `report.json` → `matrix_provenance` (always includes `pca_source`, `n_hvg`, `n_pcs`). Prefer `significant_adj_reported` over raw `significant_adj` when `low_power` is set.
