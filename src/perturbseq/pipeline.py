from __future__ import annotations

import importlib.metadata
import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

import pandas as pd
from anndata import AnnData

from perturbseq.analysis import cluster_perturbations, run_deseq2_or_wilcoxon, run_edistance, run_mixscape
from perturbseq.annotate import DEFAULT_CONTROL_PATTERNS, annotate_guides, filter_singlets
from perturbseq.io import load_dragen_sample, write_h5ad
from perturbseq.plots import plot_edistance, plot_guide_composition, plot_qc, plot_umap, plot_volcano
from perturbseq.preprocess import preprocess_rna
from perturbseq.qc import add_qc_metrics, filter_cells


@dataclass
class PipelineConfig:
    input_dir: Path
    output_dir: Path
    sample_id: str = "sample1"
    singlet_only: bool = True
    n_mads: float = 5.0
    min_cells: int = 3
    n_top_genes: int = 2000
    n_pcs: int = 30
    leiden_resolution: float = 0.5
    control: str = "NT"
    replicate_col: str | None = None
    skip_mixscape: bool = False
    skip_distance: bool = False
    skip_de: bool = False
    n_perms: int = 200
    min_cells_per_pert: int = 10
    mixscape_max_targets: int = 40
    force_mixscape: bool = False
    de_top_n: int = 10
    perturbation_type: str = "KO"
    random_state: int = 0
    control_patterns: tuple[str, ...] = DEFAULT_CONTROL_PATTERNS
    extra: dict = field(default_factory=dict)


def _ko_label(gene: str, perturbation_type: str) -> str:
    return f"{gene} {perturbation_type}"


def _mixscape_subset(adata: AnnData, control: str) -> AnnData:
    if "mixscape_class_global" not in adata.obs:
        return adata
    keep = adata.obs["mixscape_class_global"].astype(str).isin([control, "KO"])
    return adata[keep].copy()


def run_pipeline(config: PipelineConfig) -> dict:
    config.output_dir.mkdir(parents=True, exist_ok=True)
    figures = config.output_dir / "figures"
    tables = config.output_dir / "tables"
    figures.mkdir(exist_ok=True)
    tables.mkdir(exist_ok=True)
    versions = {}
    for pkg in ("perturbseq", "scanpy", "anndata", "pandas", "numpy", "pertpy"):
        try:
            versions[pkg] = importlib.metadata.version(pkg)
        except importlib.metadata.PackageNotFoundError:
            versions[pkg] = None
    report: dict = {"sample_id": config.sample_id, "versions": versions, "steps": []}

    print(f"Loading {config.sample_id} from {config.input_dir}", flush=True)
    rna, _crispr, assignments, feature_ref = load_dragen_sample(config.input_dir, config.sample_id)
    print(f"Loaded {rna.n_obs} cells × {rna.n_vars} genes", flush=True)
    rna = annotate_guides(rna, assignments, feature_ref, config.control_patterns)
    rna = add_qc_metrics(rna)
    plot_qc(rna, figures)
    plot_guide_composition(rna, figures)
    report["n_cells_loaded"] = int(rna.n_obs)
    report["n_genes_loaded"] = int(rna.n_vars)
    report["guide_counts"] = rna.obs["num_features"].value_counts().sort_index().to_dict()
    report["gene_target_counts"] = rna.obs["gene_target"].value_counts().to_dict()

    print("Running QC filters", flush=True)
    rna, qc_log = filter_cells(rna, n_mads=config.n_mads, min_cells=config.min_cells)
    report["qc"] = qc_log
    n_before_singlet = int(rna.n_obs)
    rna = filter_singlets(rna, singlet_only=config.singlet_only)
    report["n_after_singlet"] = int(rna.n_obs)
    report["n_removed_nonsinglet"] = n_before_singlet - int(rna.n_obs)
    if rna.n_obs < 20:
        raise ValueError(f"Only {rna.n_obs} cells remain after QC/singlet filters")

    print(f"Preprocessing {rna.n_obs} singlets", flush=True)
    rna = preprocess_rna(
        rna,
        n_top_genes=config.n_top_genes,
        n_pcs=config.n_pcs,
        leiden_resolution=config.leiden_resolution,
        random_state=config.random_state,
    )
    umap_color = ["leiden", "perturbation"]
    if rna.obs["gene_target"].nunique() <= 40:
        umap_color.insert(1, "gene_target")
    plot_umap(rna, figures, color=umap_color)

    n_targets = int(
        rna.obs.loc[~rna.obs["gene_target"].isin([config.control, "unassigned"]), "gene_target"].nunique()
    )
    report["n_gene_targets"] = n_targets
    mixscape_ok = False
    skip_mixscape = config.skip_mixscape
    if not skip_mixscape and n_targets > config.mixscape_max_targets and not config.force_mixscape:
        skip_mixscape = True
        report["mixscape_skipped"] = (
            f"{n_targets} gene targets > mixscape_max_targets={config.mixscape_max_targets}; "
            "use --force-mixscape to override. Genome-scale screens are quantified with E-distance."
        )
    if not skip_mixscape:
        try:
            run_mixscape(
                rna,
                control=config.control,
                split_by=config.replicate_col,
                perturbation_type=config.perturbation_type,
            )
            mixscape_ok = True
            report["mixscape_global"] = rna.obs["mixscape_class_global"].value_counts().to_dict()
            plot_umap(rna, figures / "mixscape", color=["mixscape_class_global", "perturbation"])
        except Exception as exc:  # noqa: BLE001
            report["mixscape_error"] = str(exc)
            report["steps"].append("mixscape_failed")

    analysis_obj = _mixscape_subset(rna, config.control) if mixscape_ok else rna
    if "X_pert" in analysis_obj.layers:
        analysis_pca = analysis_obj.copy()
        analysis_pca.X = analysis_pca.layers["X_pert"]
        import scanpy as sc

        sc.pp.pca(analysis_pca, n_comps=min(config.n_pcs, analysis_pca.n_obs - 1, analysis_pca.n_vars - 1))
        analysis_obj.obsm["X_pca"] = analysis_pca.obsm["X_pca"]

    edistances = pd.DataFrame()
    if not config.skip_distance:
        try:
            print("Computing E-distance", flush=True)
            edistances, etest = run_edistance(
                analysis_obj,
                groupby="gene_target",
                contrast=config.control,
                n_perms=config.n_perms,
                min_cells=config.min_cells_per_pert,
                random_state=config.random_state,
            )
            edistances.to_csv(tables / "edistance.csv")
            plot_edistance(edistances, figures)
            if etest is not None:
                etest.to_csv(tables / "etest.csv")
            report["edistance_top"] = edistances.head(10).to_dict()
        except Exception as exc:  # noqa: BLE001
            report["distance_error"] = str(exc)

    de_tables: list[pd.DataFrame] = []
    if not config.skip_de:
        counts = rna.copy()
        if mixscape_ok:
            counts.obs["de_group"] = counts.obs["mixscape_class"].astype(str)
            groups = [
                g
                for g in counts.obs["de_group"].unique()
                if str(g).endswith(f" {config.perturbation_type}")
            ]
            reference = config.control
        else:
            counts.obs["de_group"] = counts.obs["gene_target"].astype(str)
            groups = [g for g in counts.obs["de_group"].unique() if g not in {config.control, "unassigned"}]
            reference = config.control
        if not edistances.empty:
            ranked = [g for g in edistances.index if g in set(map(str, groups))]
            groups = ranked[: config.de_top_n] if ranked else list(groups)[: config.de_top_n]
        else:
            groups = (
                counts.obs["de_group"].astype(str).value_counts().reindex(groups).sort_values(ascending=False).index.tolist()
            )[: config.de_top_n]
        report["de_groups"] = [str(g) for g in groups]
        for group in groups:
            n_group = int((counts.obs["de_group"] == group).sum())
            if n_group < config.min_cells_per_pert:
                continue
            try:
                table = run_deseq2_or_wilcoxon(
                    counts,
                    group=group,
                    reference=reference,
                    replicate_col=config.replicate_col,
                    groupby="de_group",
                )
            except Exception as exc:  # noqa: BLE001
                report.setdefault("de_errors", {})[str(group)] = str(exc)
                continue
            if table.empty:
                continue
            safe = str(group).replace(" ", "_")
            table.to_csv(tables / f"de_{safe}.csv", index=False)
            plot_volcano(table, figures, safe)
            de_tables.append(table.head(50))
        if de_tables:
            pd.concat(de_tables, ignore_index=True).to_csv(tables / "de_top50_concat.csv", index=False)
        report["n_de_contrasts"] = len(de_tables)
        if config.replicate_col is None:
            report["de_note"] = (
                "Single-sample run: DE is cell-level Wilcoxon and is exploratory. "
                "Do not treat cell-level p-values as biological-replicate inference."
            )

    try:
        pert_clusters = cluster_perturbations(analysis_obj, groupby="gene_target")
        pert_clusters.to_csv(tables / "perturbation_clusters.csv")
    except Exception as exc:  # noqa: BLE001
        report["cluster_error"] = str(exc)

    write_h5ad(rna, config.output_dir / f"{config.sample_id}.tertiary.h5ad")
    report["output_h5ad"] = str(config.output_dir / f"{config.sample_id}.tertiary.h5ad")
    report["config"] = {k: (str(v) if isinstance(v, Path) else v) for k, v in asdict(config).items() if k != "extra"}
    (config.output_dir / "report.json").write_text(json.dumps(report, indent=2, default=str))
    return report
