from __future__ import annotations

import argparse
import json
from pathlib import Path

from anndata import read_h5ad

from perturbseq.cell_annotation import annotate_cells, annotation_summary, write_annotation_tables
from perturbseq.demo import write_demo_dragen
from perturbseq.guide_qc import run_guide_qc, write_guide_qc_tables
from perturbseq.io import write_h5ad
from perturbseq.pipeline import PipelineConfig, pipeline_config_from_mapping, run_pipeline
from perturbseq.plots import plot_cell_annotation, plot_guide_qc
from perturbseq.report import write_analysis_report
from perturbseq.repro import deep_merge, default_config_path, load_yaml_config


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
    run.add_argument("--n-perms", type=int, default=None)
    run.add_argument("--skip-mixscape", action="store_true")
    run.add_argument("--force-mixscape", action="store_true")
    run.add_argument("--mixscape-max-targets", type=int, default=None)
    run.add_argument("--de-top-n", type=int, default=None)
    run.add_argument("--skip-cell-annotation", action="store_true")
    run.add_argument("--skip-distance", action="store_true")
    run.add_argument("--skip-de", action="store_true")
    run.add_argument("--perturbation-type", default=None, help="Mixscape label, e.g. KO / KD / perturbation")
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

    guide_qc = sub.add_parser("guide-qc", help="Guide-level QC / efficacy on an existing tertiary h5ad")
    guide_qc.add_argument("--h5ad", type=Path, required=True)
    guide_qc.add_argument("--output-dir", type=Path, required=True)
    guide_qc.add_argument("--control", default="NT")
    guide_qc.add_argument("--min-cells", type=int, default=10)
    guide_qc.add_argument("--min-median-umi", type=float, default=5.0)
    guide_qc.add_argument("--min-detection-rate", type=float, default=0.5)
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
    if args.n_perms is not None:
        data["n_perms"] = args.n_perms
    if args.mixscape_max_targets is not None:
        data["mixscape_max_targets"] = args.mixscape_max_targets
    if args.de_top_n is not None:
        data["de_top_n"] = args.de_top_n
    if args.perturbation_type is not None:
        data["perturbation_type"] = args.perturbation_type
    if args.random_state is not None:
        data["random_state"] = args.random_state
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
    if args.skip_cell_annotation:
        data["skip_cell_annotation"] = True
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
) -> int:
    output_dir.mkdir(parents=True, exist_ok=True)
    figures = output_dir / "figures"
    tables = output_dir / "tables"
    figures.mkdir(exist_ok=True)
    tables.mkdir(exist_ok=True)
    print(f"Loading {h5ad}", flush=True)
    adata = read_h5ad(h5ad)
    print(f"Guide QC on {adata.n_obs} cells", flush=True)
    guide_df, consistency_df, warnings_df, summary = run_guide_qc(
        adata,
        control=control,
        min_cells=min_cells,
        min_median_umi=min_median_umi,
        min_detection_rate=min_detection_rate,
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
        )
    if args.command == "report":
        report = write_analysis_report(
            args.output_dir,
            control=args.control,
            min_cells=args.min_cells,
            html_detail_limit=args.html_detail_limit,
        )
        n = report.get("perturbation_summary", {}).get("n_perturbations", 0)
        print(f"Wrote {args.output_dir / 'report.html'} ({n} perturbations)")
        print(f"Updated {args.output_dir / 'report.json'} → perturbations")
        return 0
    config, config_path = build_run_config(args)
    report = run_pipeline(config, resume=args.resume, config_path=config_path)
    print(f"Finished. Report: {config.output_dir / 'report.json'}")
    print(f"HTML report: {config.output_dir / 'report.html'}")
    print(f"Manifest: {config.output_dir / 'run_manifest.json'}")
    print(f"Run log: {config.output_dir / 'run.log'}")
    print(f"Cells after singlet filter: {report.get('n_after_singlet')}")
    if "cell_annotation" in report:
        print(f"Cell annotation: {report['cell_annotation'].get('cell_state')}")
    if "mixscape_global" in report:
        print(f"Mixscape: {report['mixscape_global']}")
    if "guide_qc" in report:
        gq = report["guide_qc"]
        print(
            f"Guide QC: warnings={gq.get('n_warnings')} "
            f"consistent={gq.get('n_genes_consistent_guides')} "
            f"inconsistent={gq.get('n_genes_inconsistent_guides')}"
        )
    if "perturbation_summary" in report:
        print(f"Perturbations summarized: {report['perturbation_summary'].get('n_perturbations')}")
    if "de_note" in report:
        print(report["de_note"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
