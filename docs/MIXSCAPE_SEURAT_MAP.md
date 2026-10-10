# Pipeline ↔ Seurat Mixscape vignette

Canonical external reference for Mixscape concepts used by this repo:

**[Seurat Mixscape vignette](https://satijalab.org/seurat/articles/mixscape_vignette)** (Papalexi / Satija lab; ECCITE-seq THP-1 demo)

Paper: Papalexi et al., *Nature Genetics* (2021) — [doi:10.1038/s41588-021-00778-2](https://doi.org/10.1038/s41588-021-00778-2).  
Python implementation here: **pertpy** `pt.tl.Mixscape` (Heumos et al. 2026), wired in `perturbation.run_mixscape` / stage `3a_perturbation_modeling`.

Caveats (NP ≠ proven null, KO vs KD, auto-skip): [STATISTICAL_CAVEATS.md](STATISTICAL_CAVEATS.md). Matrix layers / pre–post paths: [ANALYSIS_DEPENDENCIES.md](ANALYSIS_DEPENDENCIES.md).

## Step mapping (Seurat → this pipeline)

| Seurat vignette step | Seurat API (vignette) | This pipeline / pertpy | Status |
| --- | --- | --- | --- |
| Local perturbation signature (PRTB) | `CalcPerturbSig` → assay `"PRTB"` | `Mixscape.perturbation_signature` → `layers['X_pert']` | **Implemented** (`n_neighbors` default 20; optional `split_by` = `replicate_col`) |
| Classify KO vs NP (vs NT) | `RunMixscape` → `mixscape_class*` | `Mixscape.mixscape` → `obs['mixscape_class']`, `mixscape_class_global`, `mixscape_class_p_ko` (name follows `perturbation_type`) | **Implemented**; default inferential subset keeps control + KO |
| Inspect perturbation-score densities | `PlotPerturbScore` | pertpy `Mixscape.plot_perturbscore` (needs a target gene) | **Gap** — not auto-exported; scores live under `adata.uns['mixscape']` when Mixscape ran; call pertpy interactively |
| Posterior P(KO) violin | `VlnPlot(..., "mixscape_class_p_ko")` | `obs['mixscape_class_p_ko']` (+ pertpy `plot_violin`) | **Partial** — column written; no default figure |
| Mixscape DE heatmap | `MixscapeHeatmap` | pertpy `Mixscape.plot_heatmap` | **Gap** — use pipeline DE tables / HTML instead |
| LDA embedding of KO+NT | `MixscapeLDA` → LDA UMAP | `Mixscape.lda` → `uns['mixscape_lda']`; `figures/mixscape/lda_umap.png` when LDA succeeds | **Implemented** (best-effort; `report.json` → `mixscape.lda`) |
| % KO / NP / NT by guide | stacked bar of classes | composition audit + Mixscape UMAP by `mixscape_class_global` | **Partial** — counts in `report.mixscape_global` / `composition_audit.csv`; no per-sgRNA Seurat-style facet |

## Stack role: Seurat vs scanpy + pertpy

| Role | Tooling |
| --- | --- |
| **Default tertiary runtime** | Python **scanpy + pertpy** (this repo’s CLI / stages) |
| **Conceptual / vignette parity** | [Seurat Mixscape vignette](https://satijalab.org/seurat/articles/mixscape_vignette) — step map above; implemented via pertpy, not `Seurat::` R calls |
| **Optional R path** | Not shipped. Export `*.h5ad` / DE CSVs and run Seurat Mixscape interactively if you need vignette-exact plots (`PlotPerturbScore`, `MixscapeHeatmap`, etc.) |

Required biologist modules (QC, gRNA QC, UMAP, perturbation map, composition, DE) are covered in **Python**; Seurat is documentation + optional external R, not a second pipeline. Module checklist: [SC_BEST_PRACTICES_MAP.md § Tertiary module coverage](SC_BEST_PRACTICES_MAP.md#tertiary-module-coverage-checklist).

## Honest gaps (intentionally not a Seurat port)

- No Seurat object / R glue in-tree; names follow **pertpy** (`X_pert`, snake_case `obs` columns).
- Vignette viz helpers that need a **chosen target gene** (`PlotPerturbScore`, heatmap, P(KO) violin for one class) stay interactive via pertpy — the tertiary CLI does not invent per-gene figure APIs.
- Primary effect-size path after Mixscape is **E-distance on PCA of `X_pert`**, not LDA UMAP (LDA is visualization / embedding, as in the vignette).
- Genome-scale libraries may **auto-skip** Mixscape (`mixscape_max_targets`) — see caveats; the Seurat vignette uses a small ECCITE panel.

## Quick interactive checks (after a successful Mixscape run)

```python
import pertpy as pt
# adata = ... tertiary h5ad or in-memory AnnData with mixscape_* columns
ms = pt.tl.Mixscape()
# Seurat PlotPerturbScore analogue — pick a gene present in the library:
# ms.plot_perturbscore(adata, pert_key="gene_target", target_gene="IFNGR2")
# Seurat MixscapeLDA viz (also written by the pipeline when lda succeeds):
# ms.plot_lda(adata, control="NT", return_fig=True)
print(adata.obs["mixscape_class_global"].value_counts())
print("lda" in (adata.uns.get("_perturbseq_mixscape") or {}), "mixscape_lda" in adata.uns)
```

Public ECCITE object for literature comparison: `pertpy` `papalexi_2021()` — see [tutorial.md](tutorial.md#optional--public-dataset-papalexi-2021).
