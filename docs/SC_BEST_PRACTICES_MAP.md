# Pipeline ↔ sc-best-practices mapping

Methodological basis for each tertiary-analysis stage relative to the [sc-best-practices · perturbation modeling](https://www.sc-best-practices.org/conditions/perturbation-modeling/) chapter (Heumos / Lotfollahi / Ji; pertpy-based). Stage names match `pipeline.py` / `stages/*/`: `1_input_validation` → `5_report`.

Chapter URL base: `https://www.sc-best-practices.org/conditions/perturbation-modeling/`  
Anchors below are the live HTML `id`s on that page (some H2s use `conditions-perturbation-modeling-key-takeaway-*` rather than a slugified title).

For matrix layers, pre/post-Mixscape paths, and which steps change inference scope, see [ANALYSIS_DEPENDENCIES.md](ANALYSIS_DEPENDENCIES.md). Module layout: [structure.md](structure.md).

## Mapping table

| Pipeline step | Purpose | sc-best-practices section | Notes |
| --- | --- | --- | --- |
| **1_input_validation** — DRAGEN MEX + barcode/feature ingest | Load filtered count matrix and CRISPR features | *pipeline-specific / not in chapter* | Chapter starts from a prepared MuData / AnnData artifact, not Illumina DRAGEN secondary outputs. |
| **1 / 2** — Guide annotation from DRAGEN assignments | Attach `guide_id`, `gene_target`, `perturbation`, NT vs targeting | [Assigning guides](https://www.sc-best-practices.org/conditions/perturbation-modeling/#assigning-guides) | **Default** keeps DRAGEN GMM calls (`positive_cell_guide_assignments.csv`). Chapter assigns the max-count guide with a count threshold via `pertpy.pp.GuideAssignment`. Same goal; different caller by default. |
| **2** — Optional `guide_reassign` (compare / apply max-UMI or GMM) | Tertiary check or override of DRAGEN calls | [Assigning guides](https://www.sc-best-practices.org/conditions/perturbation-modeling/#assigning-guides) (related) | Off by default. Closest to the chapter’s max-guide logic when `apply_max` is used; still pipeline-specific plumbing. |
| **2_preprocessing_qc** — QC (n_counts / n_genes / %MT, MAD; optional perturbation-aware) + singlet (`num_features == 1`) | Drop low-quality / multi-guide cells | *pipeline-specific / not in chapter* | Prefer `--perturbation-aware-qc` (MAD on NT only). Per-guide QC loss → `qc_filter_by_guide.csv` / cytotoxicity flags. See `composition_audit.csv`. |
| **2** — `normalize_total` + `log1p`, keep `layers['counts']`; HVG / PCA / neighbors | Expression representation + embedding for downstream tools | Prep shown inside [Identifying perturbed cells](https://www.sc-best-practices.org/conditions/perturbation-modeling/#conditions-perturbation-modeling-key-takeaway-2) and [Comparing perturbations at scale](https://www.sc-best-practices.org/conditions/perturbation-modeling/#conditions-perturbation-modeling-key-takeaway-3) | Same scanpy pattern the chapter uses before Mixscape and before E-distance. Not itself evidence of a KO effect. |
| **2** — UMAP / Leiden | Structure visualization / descriptive clusters | *not claimed as perturbation evidence* (chapter focus is Mixscape / E-distance / DE / clustering of perturbations) | Aligns with book takeaways: do not treat global UMAP separation as KO proof. See [STATISTICAL_CAVEATS.md](STATISTICAL_CAVEATS.md). |
| **2** — Cell cycle / state annotation | Descriptive labels on cells | *pipeline-specific / not in chapter* | Optional; not used as default perturbation evidence. |
| **3a_perturbation_modeling** — Mixscape (`X_pert`, KO/NP) | Identify cells where the perturbation “worked” | [Identifying perturbed cells](https://www.sc-best-practices.org/conditions/perturbation-modeling/#conditions-perturbation-modeling-key-takeaway-2) | Direct match (`pertpy.tl.Mixscape`). NP ≠ proven failed edit. Optional KO+NT subset mirrors chapter’s focus on perturbed vs control. Seurat vignette parity: [MIXSCAPE_SEURAT_MAP.md](MIXSCAPE_SEURAT_MAP.md). |
| **3a** — E-distance / E-test (+ CI / power) | Quantify and test transcriptome shift vs control | [Effect size](https://www.sc-best-practices.org/conditions/perturbation-modeling/#effect-size) (under Comparing perturbations at scale) | Direct match (`pertpy.tl.Distance` / `DistanceTest` on `X_pca`). Pre/post-Mixscape embeddings; `n_bootstrap` → `edistance_ci_*`; `etest_power_min_cells` → `low_power` / `significant_adj_reported`. |
| **3a** — Perturbation-space clustering (`--perturbation-space`) | Group perturbations with similar effects | [Perturbations with similar effects](https://www.sc-best-practices.org/conditions/perturbation-modeling/#perturbations-with-similar-effects) | Same goal. Methods: `pca_silhouette` (default; mean PCA + Leiden), `kmeans`, `lr_classifier`. Descriptive only — see [STATISTICAL_CAVEATS.md](STATISTICAL_CAVEATS.md). |
| **3a** — Guide QC / multi-guide consistency / on-target proxy | Per-guide assignment, target effect, cytotoxicity; gene-level agreement | *pipeline-specific / not in chapter* | `guide_qc.csv` / `gene_guide_consistency.csv` / `potential_low_efficiency`; HTML cards + optional report PNGs. |
| **3b_statistical_inference** — Pseudobulk + PyDESeq2 (preferred default) | Genes affected by each perturbation | [Differential expression per perturbation](https://www.sc-best-practices.org/conditions/perturbation-modeling/#differential-expression-per-perturbation) | Direct match (`PseudobulkSpace` + `PyDESeq2`). Bio reps → `inferential`; else `sample_id` / `--pseudo-replicates` → exploratory `pydeseq2_pseudobulk_no_bio_reps`. |
| **3b** — Wilcoxon on log-norm cells (&lt;2 pseudobulk units) | Exploratory DE only | *partial — chapter recommends pseudobulk DE* | Loud fallback when single-sample and no pseudo-reps; not a substitute for replicate-aware design (Squair 2021). |
| **4_robustness** — Sensitivity / effect consistency / confidence flags | Flag preprocessing- or guide-sensitive claims | *pipeline-specific / not in chapter* | Operational hardening around the same scientific claims the chapter motivates. |
| **5_report** — H5AD / CSV / JSON / HTML / provenance / report plots | Deliverables and machine-readable summaries | *pipeline-specific / not in chapter* | Includes `guide_consistency_*` / `target_validation_*` PNGs (`report_plot_top_n`); not part of the book’s notebook narrative. |
| *(not a pipeline stage)* — Unseen-perturbation prediction | Predict responses never measured | [Predicting unseen perturbations](https://www.sc-best-practices.org/conditions/perturbation-modeling/#predicting-unseen-perturbations) | **Out of scope.** Documented as external / non-default (scGen, foundation models); chapter cautions that simple baselines often win. |

## Chapter sections (reference)

| Section title (as published) | Anchor |
| --- | --- |
| Motivation | [`#conditions-perturbation-modeling-key-takeaway-1`](https://www.sc-best-practices.org/conditions/perturbation-modeling/#conditions-perturbation-modeling-key-takeaway-1) |
| Dataset | [`#dataset`](https://www.sc-best-practices.org/conditions/perturbation-modeling/#dataset) |
| Assigning guides | [`#assigning-guides`](https://www.sc-best-practices.org/conditions/perturbation-modeling/#assigning-guides) |
| Identifying perturbed cells | [`#conditions-perturbation-modeling-key-takeaway-2`](https://www.sc-best-practices.org/conditions/perturbation-modeling/#conditions-perturbation-modeling-key-takeaway-2) |
| Differential expression per perturbation | [`#differential-expression-per-perturbation`](https://www.sc-best-practices.org/conditions/perturbation-modeling/#differential-expression-per-perturbation) |
| Comparing perturbations at scale | [`#conditions-perturbation-modeling-key-takeaway-3`](https://www.sc-best-practices.org/conditions/perturbation-modeling/#conditions-perturbation-modeling-key-takeaway-3) |
| Effect size | [`#effect-size`](https://www.sc-best-practices.org/conditions/perturbation-modeling/#effect-size) |
| Perturbations with similar effects | [`#perturbations-with-similar-effects`](https://www.sc-best-practices.org/conditions/perturbation-modeling/#perturbations-with-similar-effects) |
| Predicting unseen perturbations | [`#predicting-unseen-perturbations`](https://www.sc-best-practices.org/conditions/perturbation-modeling/#predicting-unseen-perturbations) |

## Honest boundary

This pipeline is **aligned with** the chapter’s path (assign guides → identify perturbed cells → DE on pseudobulks → E-distance / cluster similar perturbations), not a line-by-line reproduction. DRAGEN ingest, default GMM assignment, MAD/singlet QC, guide-level QC, robustness flags, and the HTML/JSON report are **pipeline-specific**. Unseen-perturbation prediction is **intentionally omitted**.

## Tertiary module coverage checklist

Single status table for the biologist-facing modules (DRAGEN inputs → scanpy/pertpy tertiary; ICM walkthrough-style deliverables). Status: **done** / **partial** / **gap**. Seurat is **not** the default runtime — see [MIXSCAPE_SEURAT_MAP.md](MIXSCAPE_SEURAT_MAP.md).

External refs:

- [DRAGEN v4.5 CRISPR mode](https://help.dragen.illumina.com/dragen-v4.5/product-guides/dragen-v4.5/dragen-single-cell-pipeline/dragen-scrna-illumina#crispr-mode) — secondary inputs (`*.filtered.*.gz`, `feature_barcode_reference.csv`, `positive_cell_guide_assignments.csv`)
- [sc-best-practices · perturbation modeling](https://www.sc-best-practices.org/conditions/perturbation-modeling/)
- [Illumina Connected Multiomics Perturb-seq walkthrough](https://help.multiomics.illumina.com/icm/analyses/walkthroughs/perturb-seq)

| Module | Recommended outputs | Priority | Status | Pipeline artifacts |
| --- | --- | --- | --- | --- |
| 细胞 QC | QC violin, scatter, before/after cell counts | 必做 | **done** | `figures/qc_violin.png`, `qc_scatter_*.png`, `qc_histograms.png`, `qc_cell_counts.png` + `tables/qc_cell_counts.csv` |
| gRNA QC | guide cell-count barplot, guides-per-cell dist, assignment summary | 必做 | **done** | `figures/guide_cell_counts.png`, `guide_composition.png` (bar labels), `tables/guide_assignment_summary.csv` (+ `guide_qc.csv`) |
| UMAP / clustering | cluster UMAP, marker table, cell-state annotation | 必做 | **done** | `figures/umap.png`, `umap_cell_annotation.png`, `tables/cluster_markers.csv`, `cell_annotations.csv` |
| Perturbation 映射 | UMAP by guide, target gene, NTC | 必做 | **done** | `figures/umap_perturbation.png` (`guide_id` / `gene_target` / `is_ntc`) |
| 细胞组成变化 | per-perturbation proportions/counts in clusters | 推荐 | **done** | `tables/composition_by_perturbation.csv`, `figures/composition_by_perturbation.png` (+ stage `composition_audit.csv`) |
| 差异表达 | DEG vs NTC per perturbation, volcano, heatmap | when suitable control | **done** | `tables/de_*.csv`, `figures/volcano_*.png`, `figures/de_heatmap.png` (skipped if no contrasts) |
| 通路分析 | GO/GSEA, enrichment plots | when reliable DEG/ranking | **gap** | Not run by default; DE gene lists in HTML/CSV only. Optional external enrichment on `de_*.csv` |

HTML embedding: all `figures/*.png` are linked from `report.html`; per-perturbation cards prefer `umap_perturbation.png`, volcano, and shared `de_heatmap.png` when present.
