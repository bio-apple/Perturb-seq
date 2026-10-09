# Statistical caveats

Limitations that apply when interpreting this pipeline’s outputs. Step-level matrix provenance: [ANALYSIS_DEPENDENCIES.md](ANALYSIS_DEPENDENCIES.md). Method choices: README «方法选型».

## Cells are not biological replicates

Treating single cells as independent samples for DE / testing **inflates** significance and confuses technical with biological variation. Prefer **pseudobulk across true replicates** when `replicate_col` has ≥2 levels (PyDESeq2 path). Without replicates, cell-level Wilcoxon results are **exploratory only** — not between-condition inferential claims.

- Squair et al., *Nature Communications* (2021): confronting false discoveries in single-cell differential expression — [doi:10.1038/s41467-021-25960-2](https://doi.org/10.1038/s41467-021-25960-2)

## Mixscape: NP ≠ proven biological null

Mixscape classifies cells that received a targeting guide but still resemble controls as **NP** (non-perturbed), vs **KO** (perturbed). That is a **transcriptomic** classification relative to the NT pool:

- NP does **not** prove the guide failed, the edit was absent, or the gene is non-essential in another assay.
- KO filtering changes who enters E-distance / DE (see composition audit).
- Large libraries: default skip when `#targets > mixscape_max_targets` unless `--force-mixscape` — genome-scale Mixscape is expensive and pool-dependent.

- Papalexi et al., *Nature Genetics* (2021): Mixscape / ECCITE-seq — [doi:10.1038/s41588-021-00778-2](https://doi.org/10.1038/s41588-021-00778-2)
- Heumos et al., *Nature Methods* (2026): pertpy framework — [doi:10.1038/s41592-025-02909-7](https://doi.org/10.1038/s41592-025-02909-7)

## E-distance / E-test

E-distance quantifies multivariate shift (typically in PCA space) vs control; E-test uses permutations (`n_perms`). Embedding choice (pre- vs post-Mixscape PCA, `X_pert`) changes ranks — check `report.json` → `matrix_provenance` / pca_source. Distance ≠ mechanism.

- Peidli et al., *Nature Methods* (2024): scPerturb / E-distance — [doi:10.1038/s41592-024-02296-5](https://doi.org/10.1038/s41592-024-02296-5)
- pertpy Distance / DistanceTest: Heumos et al. 2026 (above)

## Wilcoxon without replicates is exploratory

If there is no usable biological replicate column, the pipeline may still emit Wilcoxon DE tables. Use them for **ranking / hypothesis generation**, not for FDR-style claims about population effects. With replicates → pseudobulk + PyDESeq2 (Squair 2021).

## UMAP / Leiden are not perturbation evidence

Global HVG → PCA → neighbors → Leiden / UMAP describe **dataset structure** (batch, cell cycle, ambient RNA, guide load, …). Separation on UMAP or Leiden membership is **not** evidence of a KO effect. Prefer E-distance / E-test and (when available) replicate-aware DE.

Aligned with [sc-best-practices · perturbation modeling](https://www.sc-best-practices.org/conditions/perturbation-modeling/) and [ANALYSIS_DEPENDENCIES.md](ANALYSIS_DEPENDENCIES.md).

## Filtering changes composition

QC (MAD), singlet filter (`num_features == 1`), and Mixscape KO+NT subset all **drop cells** and can deplete rare guides or alter state mix. Compare stages in `tables/composition_audit.csv`. Preprocessing-sensitive hits should be treated cautiously.

## Guide inconsistency vs biology

`gene_guide_consistency.csv` / guide QC flag same-gene guides with discordant effect directions. Do **not** pool inconsistent guides into one “gene KO” conclusion without inspecting assignment quality, UMI, cytotoxicity proxies, and possible off-target / incomplete edits. Inconsistency can be technical **or** biological; the tables separate “weak assignment” from “assigned but no target effect” via `interpretation` / warnings — still requires human review.

## Guide assignment is inherited from DRAGEN

This tertiary pipeline does **not** re-call guides by default. Assignment quality (GMM on CRISPR capture) bounds everything downstream.

- [DRAGEN scRNA CRISPR mode](https://help.dragen.illumina.com/dragen-v4.5/product-guides/dragen-v4.5/dragen-single-cell-pipeline/dragen-scrna-illumina#crispr-mode)
- [Illumina Connected Multiomics · Perturb-seq walkthrough](https://help.multiomics.illumina.com/icm/analyses/walkthroughs/perturb-seq)

## Practical checklist

| Claim you want | Prefer | Avoid treating as proof |
| --- | --- | --- |
| Guide worked / cell perturbed | Mixscape KO + target knockdown / multi-guide consistency | UMAP blob alone |
| Perturbation shifts transcriptome | E-distance + E-test (tagged embedding) | Leiden cluster label |
| Gene-level DE between conditions | Pseudobulk + replicates | Cell-level Wilcoxon p-values alone |
| Gene is validated KO | Orthogonal assay / consistent guides | Single NP-rich guide merged with KO cells |
