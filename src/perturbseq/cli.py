from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from anndata import read_h5ad

from perturbseq.annotation_policy import SAMPLE_TYPES
from perturbseq.cell_annotation import annotate_cells, annotation_summary, write_annotation_tables
from perturbseq.demo import write_demo_dragen
from perturbseq.guide_qc import GUIDE_MERGE_MODES, run_guide_qc, write_guide_qc_tables
from perturbseq.guide_reassignment import GUIDE_REASSIGN_MODES
from perturbseq.io import write_h5ad
from perturbseq.pipeline import (
    PipelineConfig,
    dry_run_plan,
    pipeline_config_from_mapping,
    run_pipeline,
)
from perturbseq.plots import plot_cell_annotation, plot_guide_qc
from perturbseq.report import write_analysis_report
from perturbseq.repro import (
    deep_merge,
    default_config_path,
    load_yaml_config,
    write_config_template,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Tertiary Perturb-seq analysis from DRAGEN CRISPR-mode outputs (pertpy-first)."
    )
    sub = parser.add_subparsers(dest="command", required=True)

    demo = sub.add_parser("write-demo", help="Write a tiny DRAGEN-like input bundle")
    demo.add_argument("--output-dir", type=Path, default=Path("data/demo"))

    run = sub.add_parser("run", help="Run the tertiary pipeline")
    run.add_argument(
        "--config",
        type=Path,
        default=None,
        help=f"YAML config (default search: {default_config_path()})",
    )
    run.add_argument("--input-dir", type=Path, default=None)
    run.add_argument("--output-dir", type=Path, default=None)
    run.add_argument("--sample-id", default=None)
    run.add_argument("--control", default=None)
    run.add_argument("--replicate-col", default=None)
    run.add_argument("--keep-multiplets", action="store_true")
    run.add_argument("--n-mads", type=float, default=None)
    run.add_argument(
        "--perturbation-aware-qc",
        action="store_true",
        help=(
            "Fit MAD QC thresholds on NT/control cells only, then apply to all cells. "
            "Preferred for Perturb-seq so strong phenotypes are not used to set cutoffs that drop them."
        ),
    )
    run.add_argument("--n-perms", type=int, default=None)
    run.add_argument(
        "--min-cells-per-pert",
        type=int,
        default=None,
        help="Min cells per perturbation for distance / DE grouping",
    )
    run.add_argument(
        "--etest-power-min-cells",
        type=int,
        default=None,
        help="Below this n_cells, mark E-test as low_power (default 50; Peidli ~200 more stable)",
    )
    run.add_argument(
        "--secondary-distance-metrics",
        default=None,
        help="Comma-separated secondary metrics after E-distance (e.g. mmd,wasserstein)",
    )
    run.add_argument(
        "--n-bootstrap",
        type=int,
        default=None,
        help="Cell bootstrap replicates for E-distance CI (default 100; 0 = skip)",
    )
    run.add_argument("--skip-mixscape", action="store_true")
    run.add_argument("--force-mixscape", action="store_true")
    run.add_argument("--mixscape-max-targets", type=int, default=None)
    run.add_argument(
        "--mixscape-mode",
        default=None,
        choices=["auto", "skip", "force", "subset"],
        help="auto (default: skip if n_targets>max), skip, force (all targets), or subset",
    )
    run.add_argument(
        "--mixscape-targets",
        default=None,
        help="Comma-separated gene targets for subset Mixscape (implies mixscape_mode=subset)",
    )
    run.add_argument(
        "--mixscape-top-n",
        type=int,
        default=None,
        help="Run Mixscape on top-N E-distance genes (+ control); implies subset mode",
    )
    run.add_argument("--de-top-n", type=int, default=None)
    run.add_argument(
        "--report-plot-top-n",
        type=int,
        default=None,
        help="Max genes for report guide-consistency / target-validation PNGs (default 15)",
    )
    run.add_argument(
        "--n-jobs",
        type=int,
        default=None,
        help="Parallel workers for per-group DE / E-test (default 1; -1 = all CPUs)",
    )
    run.add_argument(
        "--sample-type",
        default=None,
        choices=list(SAMPLE_TYPES),
        help=(
            "Sample biology for annotation policy: cell_line (default; skip full annotation), "
            "primary/mixed (run), unknown (skip with note). "
            "Overridden by --skip-cell-annotation / --run-cell-annotation. "
            "Clustering/UMAP always run."
        ),
    )
    ann_grp = run.add_mutually_exclusive_group()
    ann_grp.add_argument(
        "--skip-cell-annotation",
        action="store_true",
        help="Force skip cell-cycle / state annotation (overrides --sample-type)",
    )
    ann_grp.add_argument(
        "--run-cell-annotation",
        action="store_true",
        help="Force run cell-cycle / state annotation (overrides --sample-type)",
    )
    run.add_argument("--skip-distance", action="store_true")
    run.add_argument("--skip-de", action="store_true")
    run.add_argument("--perturbation-type", default=None, help="Mixscape label, e.g. KO / KD / perturbation")
    run.add_argument(
        "--perturbation-space",
        default=None,
        choices=["pca_silhouette", "kmeans", "lr_classifier"],
        help=(
            "Perturbation-space construction: pca_silhouette (default; mean PCA + Leiden), "
            "kmeans (mean-PCA KMeans), lr_classifier (LR-coefficient embedding)"
        ),
    )
    run.add_argument(
        "--control-patterns",
        default=None,
        help="Comma-separated regexes matched against guide/target names",
    )
    run.add_argument(
        "--resume",
        action="store_true",
        help="Skip stages that already succeeded with unchanged inputs/params",
    )
    run.add_argument("--random-state", type=int, default=None, help="Random seed for PCA/neighbors/UMAP/etc.")
    run.add_argument(
        "--guide-merge",
        default=None,
        choices=list(GUIDE_MERGE_MODES),
        help="Multi-guide weighted summary: none (default, flag only) | equal | umi | confidence | umi_confidence",
    )
    run.add_argument(
        "--on-target-lfc-cutoff",
        type=float,
        default=None,
        help="Guide QC on-target |log2FC| cutoff vs control (default 0.25; KO/KD/i ↓, CRISPRa ↑)",
    )
    run.add_argument(
        "--on-target-min-fail-guides",
        type=int,
        default=None,
        help="≥N adequate guides failing on-target → potential_low_efficiency (default 2)",
    )
    run.add_argument(
        "--guide-reassign",
        default=None,
        choices=list(GUIDE_REASSIGN_MODES),
        help=(
            "Optional tertiary guide re-call vs DRAGEN: off (default) | compare | apply_max | apply_gmm. "
            "compare writes agreement QC; apply_* overrides guide_id (keeps *_dragen backups)."
        ),
    )
    run.add_argument(
        "--guide-reassign-min-umi",
        type=float,
        default=None,
        help="Min CRISPR UMI for max/GMM reassignment (default 1.0)",
    )
    run.add_argument(
        "--de-covariates",
        default=None,
        help=(
            "DE covariates: 'true' (default candidates: phase,pct_counts_mt,log_n_counts), "
            "'false' to disable, or comma-separated obs columns for PyDESeq2 design"
        ),
    )
    run.add_argument(
        "--de-prefer-pseudobulk",
        default=None,
        help=(
            "Prefer sample_id(+perturbation) pseudobulk + PyDESeq2 even without bio reps "
            "(default true). 'false' forces cell-level Wilcoxon fallback."
        ),
    )
    run.add_argument(
        "--pseudo-replicates",
        type=int,
        default=None,
        dest="n_pseudo_replicates",
        help=(
            "Split cells into N≥2 technical pseudo-replicates for exploratory single-sample "
            "PyDESeq2 (not biological replicates; Squair et al. 2021)"
        ),
    )
    run.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate inputs and print planned stages + resolved params (no analysis)",
    )

    config_cmd = sub.add_parser("config", help="Config helpers (generate a commented YAML template)")
    config_cmd.add_argument(
        "action",
        nargs="?",
        choices=["generate"],
        default=None,
        help="Optional action; 'generate' is the same as --generate",
    )
    config_cmd.add_argument(
        "--generate",
        action="store_true",
        help="Write a fully commented YAML template (default.yaml + PipelineConfig fields)",
    )
    config_cmd.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("configs/pipeline.template.yaml"),
        help="Output path for --generate (default: configs/pipeline.template.yaml)",
    )

    annotate = sub.add_parser("annotate", help="Add cell annotation to an existing tertiary h5ad")
    annotate.add_argument("--h5ad", type=Path, required=True)
    annotate.add_argument("--output-dir", type=Path, required=True)
    annotate.add_argument("--inplace", action="store_true", help="Overwrite the input h5ad")

    report_cmd = sub.add_parser(
        "report",
        help="Rebuild HTML/JSON analysis report from an existing results directory (no full re-run)",
    )
    report_cmd.add_argument("--output-dir", type=Path, required=True)
    report_cmd.add_argument("--control", default=None, help="Override control label (default: from report.json)")
    report_cmd.add_argument("--min-cells", type=int, default=None)
    report_cmd.add_argument(
        "--html-detail-limit",
        type=int,
        default=None,
        help="Optional cap on detailed HTML sections (DE targets always included)",
    )
    report_cmd.add_argument(
        "--h5ad",
        type=Path,
        default=None,
        help="Optional expression h5ad for guide-consistency / target-validation plots "
        "(auto-discovers *.tertiary.h5ad in output-dir when omitted)",
    )
    report_cmd.add_argument(
        "--report-plot-top-n",
        type=int,
        default=None,
        help="Max genes for guide-consistency / target-validation PNGs (default: from config / 15)",
    )

    guide_qc = sub.add_parser("guide-qc", help="Guide-level QC / efficacy on an existing tertiary h5ad")
    guide_qc.add_argument("--h5ad", type=Path, required=True)
    guide_qc.add_argument("--output-dir", type=Path, required=True)
    guide_qc.add_argument("--control", default="NT")
    guide_qc.add_argument("--min-cells", type=int, default=10)
    guide_qc.add_argument("--min-median-umi", type=float, default=5.0)
    guide_qc.add_argument("--min-detection-rate", type=float, default=0.5)
    guide_qc.add_argument(
        "--guide-merge",
        default="none",
        choices=list(GUIDE_MERGE_MODES),
        help="Weighted gene-level guide summary (default none = flag inconsistency only)",
    )
    guide_qc.add_argument(
        "--perturbation-type",
        default="KO",
        help="Perturbation modality for on-target direction (KO/KD/CRISPRi → ↓, CRISPRa → ↑)",
    )
    guide_qc.add_argument(
        "--on-target-lfc-cutoff",
        type=float,
        default=0.25,
        help="|log2FC| vs control required for on_target_pass (default 0.25)",
    )
    guide_qc.add_argument(
        "--on-target-min-fail-guides",
        type=int,
        default=2,
        help="≥N adequate guides failing on-target → gene potential_low_efficiency (default 2)",
    )
    return parser


def _cli_overrides(args: argparse.Namespace) -> dict:
    """Map CLI args to PipelineConfig fields; omit unset optional values."""
    data: dict = {}
    if args.input_dir is not None:
        data["input_dir"] = args.input_dir
    if args.output_dir is not None:
        data["output_dir"] = args.output_dir
    if args.sample_id is not None:
        data["sample_id"] = args.sample_id
    if args.control is not None:
        data["control"] = args.control
    if args.replicate_col is not None:
        data["replicate_col"] = args.replicate_col
    if args.n_mads is not None:
        data["n_mads"] = args.n_mads
    if getattr(args, "perturbation_aware_qc", False):
        data["perturbation_aware_qc"] = True
    if args.n_perms is not None:
        data["n_perms"] = args.n_perms
    if getattr(args, "min_cells_per_pert", None) is not None:
        data["min_cells_per_pert"] = args.min_cells_per_pert
    if getattr(args, "etest_power_min_cells", None) is not None:
        data["etest_power_min_cells"] = args.etest_power_min_cells
    if getattr(args, "secondary_distance_metrics", None) is not None:
        data["secondary_distance_metrics"] = tuple(
            p.strip() for p in args.secondary_distance_metrics.split(",") if p.strip()
        )
    if getattr(args, "n_bootstrap", None) is not None:
        data["n_bootstrap"] = args.n_bootstrap
    if args.mixscape_max_targets is not None:
        data["mixscape_max_targets"] = args.mixscape_max_targets
    if getattr(args, "mixscape_mode", None) is not None:
        data["mixscape_mode"] = args.mixscape_mode
    if getattr(args, "mixscape_targets", None) is not None:
        data["mixscape_targets"] = tuple(
            t.strip() for t in args.mixscape_targets.split(",") if t.strip()
        )
    if getattr(args, "mixscape_top_n", None) is not None:
        data["mixscape_top_n"] = args.mixscape_top_n
    if args.de_top_n is not None:
        data["de_top_n"] = args.de_top_n
    if getattr(args, "report_plot_top_n", None) is not None:
        data["report_plot_top_n"] = args.report_plot_top_n
    if getattr(args, "n_jobs", None) is not None:
        data["n_jobs"] = args.n_jobs
    if args.perturbation_type is not None:
        data["perturbation_type"] = args.perturbation_type
    if getattr(args, "perturbation_space", None) is not None:
        data["perturbation_space"] = args.perturbation_space
    if args.random_state is not None:
        data["random_state"] = args.random_state
    if getattr(args, "guide_merge", None) is not None:
        data["guide_merge"] = args.guide_merge
    if getattr(args, "on_target_lfc_cutoff", None) is not None:
        data["on_target_lfc_cutoff"] = args.on_target_lfc_cutoff
    if getattr(args, "on_target_min_fail_guides", None) is not None:
        data["on_target_min_fail_guides"] = args.on_target_min_fail_guides
    if getattr(args, "guide_reassign", None) is not None:
        data["guide_reassign"] = args.guide_reassign
    if getattr(args, "guide_reassign_min_umi", None) is not None:
        data["guide_reassign_min_umi"] = args.guide_reassign_min_umi
    if getattr(args, "de_covariates", None) is not None:
        raw = str(args.de_covariates).strip()
        low = raw.lower()
        if low in {"true", "yes", "1"}:
            data["de_covariates"] = True
        elif low in {"false", "no", "0", "none"}:
            data["de_covariates"] = False
        else:
            data["de_covariates"] = tuple(p.strip() for p in raw.split(",") if p.strip())
    if getattr(args, "de_prefer_pseudobulk", None) is not None:
        raw = str(args.de_prefer_pseudobulk).strip().lower()
        data["de_prefer_pseudobulk"] = raw in {"true", "yes", "1"}
    if getattr(args, "n_pseudo_replicates", None) is not None:
        data["n_pseudo_replicates"] = args.n_pseudo_replicates
    if args.control_patterns is not None:
        data["control_patterns"] = tuple(
            p.strip() for p in args.control_patterns.split(",") if p.strip()
        )
    # Booleans: only override when flag present
    if args.keep_multiplets:
        data["singlet_only"] = False
    if args.skip_mixscape:
        data["skip_mixscape"] = True
    if args.force_mixscape:
        data["force_mixscape"] = True
    if getattr(args, "sample_type", None) is not None:
        data["sample_type"] = args.sample_type
    if args.skip_cell_annotation:
        data["skip_cell_annotation"] = True
    elif getattr(args, "run_cell_annotation", False):
        data["skip_cell_annotation"] = False
    if args.skip_distance:
        data["skip_distance"] = True
    if args.skip_de:
        data["skip_de"] = True
    return data


def build_run_config(args: argparse.Namespace) -> tuple[PipelineConfig, Path | None]:
    config_path: Path | None = args.config
    base: dict = {}
    if config_path is None and default_config_path().is_file():
        config_path = default_config_path()
    if config_path is not None:
        base = load_yaml_config(config_path)
    merged = deep_merge(base, _cli_overrides(args))
    return pipeline_config_from_mapping(merged), config_path


def _print_dry_run(config: PipelineConfig, config_path: Path | None) -> int:
    plan = dry_run_plan(config)
    print("Dry-run OK — inputs validated (no analysis run).", flush=True)
    if config_path is not None:
        print(f"Config: {config_path}", flush=True)
    print(f"Sample: {plan['sample_id']}", flush=True)
    print(f"Input:  {plan['input_dir']}", flush=True)
    print(f"Output: {plan['output_dir']}", flush=True)
    print("\nInput files:", flush=True)
    for key, path in plan["input_files"].items():
        print(f"  {key}: {path}", flush=True)
    print("\nPlanned stages:", flush=True)
    for entry in plan["stages"]:
        notes = f"  ({', '.join(entry['notes'])})" if entry["notes"] else ""
        print(f"  {entry['stage']:28s} {entry['action']}{notes}", flush=True)
    print("\nResolved parameters:", flush=True)
    for key, value in sorted(plan["params"].items()):
        print(f"  {key}: {value}", flush=True)
    return 0


def _run_annotate(h5ad: Path, output_dir: Path, inplace: bool) -> int:
    output_dir.mkdir(parents=True, exist_ok=True)
    figures = output_dir / "figures"
    tables = output_dir / "tables"
    figures.mkdir(exist_ok=True)
    tables.mkdir(exist_ok=True)
    print(f"Loading {h5ad}", flush=True)
    adata = read_h5ad(h5ad)
    print(f"Annotating {adata.n_obs} cells", flush=True)
    adata, markers = annotate_cells(adata)
    write_annotation_tables(adata, markers, tables)
    plot_cell_annotation(adata, figures)
    summary = annotation_summary(adata)
    out_h5ad = h5ad if inplace else output_dir / h5ad.name.replace(".h5ad", ".annotated.h5ad")
    write_h5ad(adata, out_h5ad)
    (output_dir / "cell_annotation_report.json").write_text(json.dumps(summary, indent=2, default=str))
    print(f"Wrote annotations to {tables}", flush=True)
    print(f"Wrote {out_h5ad}", flush=True)
    print(f"phase: {summary.get('phase')}", flush=True)
    print(f"cell_state: {summary.get('cell_state')}", flush=True)
    return 0


def _run_guide_qc(
    h5ad: Path,
    output_dir: Path,
    control: str,
    min_cells: int,
    min_median_umi: float,
    min_detection_rate: float,
    guide_merge: str = "none",
    perturbation_type: str = "KO",
    on_target_lfc_cutoff: float = 0.25,
    on_target_min_fail_guides: int = 2,
) -> int:
    output_dir.mkdir(parents=True, exist_ok=True)
    figures = output_dir / "figures"
    tables = output_dir / "tables"
    figures.mkdir(exist_ok=True)
    tables.mkdir(exist_ok=True)
    print(f"Loading {h5ad}", flush=True)
    adata = read_h5ad(h5ad)
    print(
        f"Guide QC on {adata.n_obs} cells "
        f"(guide_merge={guide_merge}, perturbation_type={perturbation_type})",
        flush=True,
    )
    guide_df, consistency_df, warnings_df, summary = run_guide_qc(
        adata,
        control=control,
        min_cells=min_cells,
        min_median_umi=min_median_umi,
        min_detection_rate=min_detection_rate,
        guide_merge=guide_merge,
        perturbation_type=perturbation_type,
        on_target_lfc_cutoff=on_target_lfc_cutoff,
        on_target_min_fail_guides=on_target_min_fail_guides,
    )
    write_guide_qc_tables(guide_df, consistency_df, warnings_df, tables)
    plot_guide_qc(guide_df, consistency_df, figures)
    (output_dir / "guide_qc_report.json").write_text(json.dumps(summary, indent=2, default=str))
    main_report = output_dir / "report.json"
    if main_report.exists():
        try:
            payload = json.loads(main_report.read_text())
        except json.JSONDecodeError:
            payload = {}
        payload["guide_qc"] = summary
        steps = payload.setdefault("steps", [])
        if "guide_qc" not in steps:
            steps.append("guide_qc")
        main_report.write_text(json.dumps(payload, indent=2, default=str))
    print(f"Wrote {tables / 'guide_qc.csv'}", flush=True)
    print(f"Wrote {tables / 'gene_guide_consistency.csv'}", flush=True)
    print(f"Wrote {tables / 'qc_warnings.csv'}", flush=True)
    if (tables / "gene_guide_weighted.csv").exists():
        print(f"Wrote {tables / 'gene_guide_weighted.csv'}", flush=True)
    print(
        f"warnings={summary.get('n_warnings')} "
        f"inconsistent_genes={summary.get('n_genes_inconsistent_guides')}",
        flush=True,
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "write-demo":
        path = write_demo_dragen(args.output_dir)
        print(f"Wrote demo DRAGEN files to {path}")
        return 0
    if args.command == "config":
        if args.generate or args.action == "generate":
            path = write_config_template(args.output)
            print(f"Wrote commented YAML template to {path}")
            return 0
        parser.error("config: pass --generate (or 'generate') to write a template")
    if args.command == "annotate":
        return _run_annotate(args.h5ad, args.output_dir, args.inplace)
    if args.command == "guide-qc":
        return _run_guide_qc(
            args.h5ad,
            args.output_dir,
            args.control,
            args.min_cells,
            args.min_median_umi,
            args.min_detection_rate,
            args.guide_merge,
            args.perturbation_type,
            args.on_target_lfc_cutoff,
            args.on_target_min_fail_guides,
        )
    if args.command == "report":
        report = write_analysis_report(
            args.output_dir,
            control=args.control,
            min_cells=args.min_cells,
            html_detail_limit=args.html_detail_limit,
            h5ad=args.h5ad,
            report_plot_top_n=args.report_plot_top_n,
        )
        n = report.get("perturbation_summary", {}).get("n_perturbations", 0)
        print(f"Wrote {args.output_dir / 'report.html'} ({n} perturbations)")
        print(f"Updated {args.output_dir / 'report.json'} → perturbations")
        plots = report.get("report_plots") or {}
        n_fig = len(plots.get("generated") or [])
        n_skip = len(plots.get("skipped_plots") or [])
        if n_fig or n_skip:
            print(f"Report plots: {n_fig} generated, {n_skip} skipped", flush=True)
        return 0
    try:
        config, config_path = build_run_config(args)
        if args.dry_run:
            return _print_dry_run(config, config_path)
        report = run_pipeline(config, resume=args.resume, config_path=config_path)
    except (FileNotFoundError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    print(f"Finished. Report: {config.output_dir / 'report.json'}")
    print(f"HTML report: {config.output_dir / 'report.html'}")
    print(f"Manifest: {config.output_dir / 'run_manifest.json'}")
    print(f"Run log: {config.output_dir / 'run.log'}")
    print(f"Cells after singlet filter: {report.get('n_after_singlet')}")
    cell_ann = report.get("cell_annotation")
    if isinstance(cell_ann, dict):
        if cell_ann.get("skipped"):
            print(f"Cell annotation: skipped ({cell_ann.get('reason')})")
        elif cell_ann.get("cell_state") is not None:
            print(f"Cell annotation: {cell_ann.get('cell_state')}")
    mix = report.get("mixscape")
    if isinstance(mix, dict):
        if mix.get("skipped"):
            print(f"Mixscape: skipped ({mix.get('reason')})")
        elif "mixscape_global" in report:
            print(f"Mixscape: {report['mixscape_global']}")
        elif mix.get("global_counts"):
            print(f"Mixscape: {mix['global_counts']}")
    elif "mixscape_global" in report:
        print(f"Mixscape: {report['mixscape_global']}")
    edist = report.get("edistance")
    if isinstance(edist, dict) and edist.get("skipped"):
        print(f"E-distance: skipped ({edist.get('reason')})")
    if "guide_qc" in report:
        gq = report["guide_qc"]
        print(
            f"Guide QC: warnings={gq.get('n_warnings')} "
            f"consistent={gq.get('n_genes_consistent_guides')} "
            f"inconsistent={gq.get('n_genes_inconsistent_guides')}"
        )
    if "perturbation_summary" in report:
        print(f"Perturbations summarized: {report['perturbation_summary'].get('n_perturbations')}")
    de = report.get("de") if isinstance(report.get("de"), dict) else {}
    de_note = de.get("detail") or report.get("de_note")
    if de_note:
        print(de_note)
    elif de.get("skipped"):
        print(f"DE: skipped ({de.get('reason')})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
