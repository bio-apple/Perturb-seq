# Statistical caveats

Limitations that apply when interpreting this pipeline’s outputs. Short “常见坑 / Common pitfalls” tables: [tutorial.md](tutorial.md#常见坑--common-pitfalls) · [tutorial.ipynb](tutorial.ipynb). Step-level matrix provenance: [ANALYSIS_DEPENDENCIES.md](ANALYSIS_DEPENDENCIES.md). Method choices: [README](../README.md) / [README.zh-CN](../README.zh-CN.md).

## Cells are not biological replicates

Treating single cells as independent samples for DE / testing **inflates** significance and confuses technical with biological variation (Squair et al. 2021). Default DE path prefers **pseudobulk + PyDESeq2**:

| Design | Method label | `evidence_level` |
| --- | --- | --- |
| `replicate_col` with ≥2 true bio-rep levels | `pydeseq2_pseudobulk` | `inferential` |
| No bio reps, but ≥2 `sample_id` levels (or `--pseudo-replicates N`) | `pydeseq2_pseudobulk_no_bio_reps` | `exploratory` |
| Single sample and no pseudo-replicates | `wilcoxon_cell_level_exploratory` (loud fallback) | `exploratory` |

Pseudobulk without true biological replicates (including technical pseudo-replicates) is **still not** a substitute for replicate-aware inference — shared dispersion can improve ranking, but population FDR claims require real bio reps.

- Squair et al., *Nature Communications* (2021): confronting false discoveries in single-cell differential expression — [doi:10.1038/s41467-021-25960-2](https://doi.org/10.1038/s41467-021-25960-2)

## Mixscape: NP ≠ proven biological null

Mixscape classifies cells that received a targeting guide but still resemble controls as **NP** (non-perturbed), vs **KO** (perturbed). That is a **transcriptomic** classification relative to the NT pool:

- NP does **not** prove the guide failed, the edit was absent, or the gene is non-essential in another assay.
- KO filtering changes who enters E-distance / DE (see composition audit).
- Large libraries: default **auto** skip when `#targets > mixscape_max_targets` unless `--force-mixscape` / `mixscape_mode=force` — genome-scale Mixscape is expensive and pool-dependent. Recorded as `report.json` → `mixscape.reason = "too_many_targets"` with `n_targets`, `mixscape_max_targets`, and a rough `estimate` (`approx_work_units` ∝ O(n_cells × n_targets), not a benchmark).
- Partial runs: `mixscape_mode=subset` with `--mixscape-targets GENE1,GENE2` and/or `--mixscape-top-n` (top E-distance genes after the cheap pre-Mixscape ranking). Results cover **only** control + selected targets — do not generalize Mixscape labels to the rest of the library. Per-batch / replicate splitting still uses `replicate_col` → Mixscape `split_by`.

**Canonical Seurat vignette** (CalcPerturbSig / RunMixscape / MixscapeLDA ↔ our steps): [MIXSCAPE_SEURAT_MAP.md](MIXSCAPE_SEURAT_MAP.md) → [satijalab.org Mixscape vignette](https://satijalab.org/seurat/articles/mixscape_vignette).

### CRISPRi/a / knockdown (`perturbation_type=KD`)

Mixscape was developed and validated primarily in a **CRISPR KO** setting (binary-ish transcriptomic shift vs NT). For **knockdown / CRISPRi/a** (`--perturbation-type KD` or similar), assumptions are **weaker**: partial repression produces graded effects, so NP vs perturbed labels are less decisive. Prefer E-distance, target-gene knockdown, and multi-guide consistency; treat Mixscape classes as exploratory filters. The HTML report surfaces this under **Statistical caveats** when `perturbation_type` is KD-like.

- Papalexi et al., *Nature Genetics* (2021): Mixscape / ECCITE-seq — [doi:10.1038/s41588-021-00778-2](https://doi.org/10.1038/s41588-021-00778-2)
- Heumos et al., *Nature Methods* (2026): pertpy framework — [doi:10.1038/s41592-025-02909-7](https://doi.org/10.1038/s41592-025-02909-7)

## E-distance / E-test

E-distance quantifies multivariate shift (typically in PCA space) vs control; E-test uses permutations (`n_perms`). Embedding choice (pre- vs post-Mixscape PCA, `X_pert`) changes ranks — check `report.json` → `matrix_provenance` (`pca_source`, `n_hvg`, `n_pcs`). Distance ≠ mechanism.

**Sample size / power.** Peidli et al. show E-test power is unstable at small per-group *n*; prefer roughly **≥50–100 cells per perturbation** for reporting significance, with **~200** more stable. This pipeline keeps computing distances for groups above `min_cells_per_pert`, but marks E-test rows with `low_power: true` when `n_cells < etest_power_min_cells` (default 50). Use `significant_adj_reported` (not raw `significant_adj`) for claims — low-power hits are never silently reported as significant.

**Bootstrap CI.** With `n_bootstrap` (default 100; `0` = skip), `edistance.csv` / `etest.csv` gain percentile `edistance_ci_low` / `edistance_ci_high` (cell resampling). HTML effect cards prefer this CI when present.

**Multi-metric sensitivity.** Alongside `edistance.csv`, the pipeline may write `distance_mmd.csv` / combined `distances.csv` (and attempt Wasserstein if OTT-JAX is installed). Skipped secondary metrics are recorded with `skipped` + `reason` (e.g. `missing_jax`).

- Peidli et al., *Nature Methods* (2024): scPerturb / E-distance — [doi:10.1038/s41592-024-02296-5](https://doi.org/10.1038/s41592-024-02296-5)
- pertpy Distance / DistanceTest: Heumos et al. 2026 (above)

## DE without true biological replicates stays exploratory

Default (`de_prefer_pseudobulk: true`): aggregate by `sample_id` + perturbation and run exploratory PyDESeq2 when ≥2 samples exist, or when the user sets `n_pseudo_replicates` / `--pseudo-replicates N` (technical splits — documented risk). Only if fewer than two pseudobulk units can be formed does the pipeline fall back to cell-level Wilcoxon, with a loud `evidence_level=exploratory` warning. Use all no-bio-rep DE for **ranking / hypothesis generation**, not population FDR. True multi-rep `replicate_col` → `evidence_level=inferential` (Squair 2021).

## UMAP / Leiden are not perturbation evidence

Global HVG → PCA → neighbors → Leiden / UMAP describe **dataset structure** (batch, cell cycle, ambient RNA, guide load, …). Separation on UMAP or Leiden membership is **not** evidence of a KO effect. Prefer E-distance / E-test and (when available) replicate-aware DE.

Aligned with [sc-best-practices · perturbation modeling](https://www.sc-best-practices.org/conditions/perturbation-modeling/) and [ANALYSIS_DEPENDENCIES.md](ANALYSIS_DEPENDENCIES.md).

## Perturbation-space method choice is descriptive

`--perturbation-space` / YAML `perturbation_space` selects how perturbations are embedded and clustered (`pca_silhouette` default | `kmeans` | `lr_classifier`). All paths write descriptive `perturbation_clusters.csv` (+ optional `perturbation_space_embeddings.csv`); they are **pathway-like groupings**, not proof of shared mechanism. Method + backend are recorded in `report.json` → `perturbation_clusters` / `matrix_provenance.perturbation_clusters`.

| Method | Embedding | Typical use |
| --- | --- | --- |
| `pca_silhouette` | Mean `X_pca` per target + Leiden (pertpy PseudobulkSpace when available) | Global effect comparison (default) |
| `kmeans` | Mean `X_pca` → sklearn KMeans (silhouette-selected *k*) | Response-pattern clusters (pertpy KMeansSpace-like intent; local) |
| `lr_classifier` | Logistic-regression coefficients (pertpy `LRClassifierSpace` when available, else local) | Marker-like / interpretable similarity space |

Do not treat cluster co-membership alone as genetic interaction evidence; use E-distance / DE / multi-guide consistency for claims.

## Filtering changes composition

QC (MAD), singlet filter (`num_features == 1`), and Mixscape KO+NT subset all **drop cells** and can deplete rare guides or alter state mix. Compare stages in `tables/composition_audit.csv`. Preprocessing-sensitive hits should be treated cautiously.

**Perturbation-aware QC (preferred for Perturb-seq):** global MAD on `n_counts` / `n_genes` / `%MT` can systematically remove the strongest phenotypes (essential-gene KO, apoptosis, arrest) because those cells define the outlier tails — a false-negative trap. Prefer `--perturbation-aware-qc` / `perturbation_aware_qc: true`: fit MAD thresholds on **NT/control cells only**, then apply to all cells. Per-guide / gene_target QC removal fractions are written to `tables/qc_filter_by_guide.csv`; guides with ≥90% cells removed get a `cytotoxicity_qc_depletion` warning in `guide_qc.csv` / `qc_warnings.csv` (not a silent discard).

## Guide inconsistency vs biology

`gene_guide_consistency.csv` / guide QC flag same-gene guides with discordant effect directions. Do **not** pool inconsistent guides into one “gene KO” conclusion without inspecting assignment quality, UMI, cytotoxicity proxies, and possible off-target / incomplete edits. Inconsistency can be technical **or** biological; the tables separate “weak assignment” from “assigned but no target effect” via `interpretation` / warnings — still requires human review.

**On-target efficiency proxy.** Per-guide `on_target_score` / `on_target_pass` check that target-gene log2FC moves in the expected direction for `perturbation_type` (KO/KD/CRISPRi → downregulation; CRISPRa → upregulation). If ≥ `on_target_min_fail_guides` (default 2) adequate guides fail this proxy, gene-level `potential_low_efficiency` / `interpretation=potential_low_efficiency_guides` marks **potential low-efficiency guides** — do **not** conclude “no phenotype” for the gene.

Optional `guide_merge` (`none` default | `equal` | `umi` | `confidence` | `umi_confidence`) adds a **descriptive** weighted summary (weighted mean target log2FC + weight-voted direction) for multi-guide genes. It does **not** drop per-guide rows, clear `guides_consistent=False`, or replace gene-level KO calls. Treat weighted columns as a review aid, not a silent merge of discordant guides.

## Guide assignment is inherited from DRAGEN

This tertiary pipeline does **not** re-call guides by default. Assignment quality (GMM on CRISPR capture) bounds everything downstream.

Optional `--guide-reassign compare|apply_max|apply_gmm` compares DRAGEN calls to a simple max-UMI baseline and a lightweight per-guide 2-component GMM on CRISPR counts (sklearn). Default remains **off** — DRAGEN secondary GMM is preferred. External alternatives (not installed here): Cell Ranger Feature Barcode, CatchR, pertpy `GuideAssignment`.

- [DRAGEN scRNA CRISPR mode](https://help.dragen.illumina.com/dragen-v4.5/product-guides/dragen-v4.5/dragen-single-cell-pipeline/dragen-scrna-illumina#crispr-mode)
- [Illumina Connected Multiomics · Perturb-seq walkthrough](https://help.multiomics.illumina.com/icm/analyses/walkthroughs/perturb-seq)

## DE covariates and Wilcoxon limits

When `de_covariates: true` (default), available obs columns among `phase`, `pct_counts_mt`, and `log_n_counts` (from `n_counts`) are included in the **PyDESeq2** design formula. **Wilcoxon** fallback cannot regress covariates via `scanpy.tl.rank_genes_groups`; requested covariates are recorded in `report.json` / table notes as ignored — still exploratory only.

## Advanced perturbation methods (not bundled)

| Method | Role | In this repo |
| --- | --- | --- |
| Mixscape `layers['X_pert']` | Perturbation signature; KO/NP; post-Mixscape E-distance PCA | **Implemented** (pertpy); Seurat vignette map: [MIXSCAPE_SEURAT_MAP.md](MIXSCAPE_SEURAT_MAP.md) |
| SCEPTRE | Causal / resampling inference for CRISPR screens | Documented only — R/Bioconductor ecosystem |
| MIMOSCA | Covariate-aware linear perturbation model | Documented only — optional external |
| PerturbNet / deep autoencoders | Learned low-dim phenotype embeddings | Documented only — do **not** pull heavy DL stacks into the default path |

Prefer Mixscape `X_pert` + E-distance here; treat SCEPTRE/MIMOSCA/PerturbNet as advanced follow-ups when you need causal identification or custom embeddings.

## Practical checklist

| Claim you want | Prefer | Avoid treating as proof |
| --- | --- | --- |
| Guide worked / cell perturbed | Mixscape KO + target knockdown / multi-guide consistency | UMAP blob alone |
| Perturbation shifts transcriptome | E-distance + E-test (tagged embedding) | Leiden cluster label |
| Gene-level DE between conditions | Pseudobulk + true bio replicates | Cell-level Wilcoxon / no-bio-rep pseudobulk as population FDR |
| Gene is validated KO | Orthogonal assay / consistent guides | Single NP-rich guide merged with KO cells |
