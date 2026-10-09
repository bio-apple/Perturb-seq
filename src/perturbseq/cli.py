from __future__ import annotations

import argparse
from pathlib import Path

from perturbseq.annotate import DEFAULT_CONTROL_PATTERNS
from perturbseq.demo import write_demo_dragen
from perturbseq.pipeline import PipelineConfig, run_pipeline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Tertiary Perturb-seq analysis from DRAGEN CRISPR-mode outputs (pertpy-first)."
    )
    sub = parser.add_subparsers(dest="command", required=True)

    demo = sub.add_parser("write-demo", help="Write a tiny DRAGEN-like input bundle")
    demo.add_argument("--output-dir", type=Path, default=Path("data/demo"))

    run = sub.add_parser("run", help="Run the tertiary pipeline")
    run.add_argument("--input-dir", type=Path, required=True)
    run.add_argument("--output-dir", type=Path, required=True)
    run.add_argument("--sample-id", default="sample1")
    run.add_argument("--control", default="NT")
    run.add_argument("--replicate-col", default=None)
    run.add_argument("--keep-multiplets", action="store_true")
    run.add_argument("--n-mads", type=float, default=5.0)
    run.add_argument("--n-perms", type=int, default=200)
    run.add_argument("--skip-mixscape", action="store_true")
    run.add_argument("--force-mixscape", action="store_true")
    run.add_argument("--mixscape-max-targets", type=int, default=40)
    run.add_argument("--de-top-n", type=int, default=10)
    run.add_argument("--skip-distance", action="store_true")
    run.add_argument("--skip-de", action="store_true")
    run.add_argument("--perturbation-type", default="KO", help="Mixscape label, e.g. KO / KD / perturbation")
    run.add_argument(
        "--control-patterns",
        default=",".join(DEFAULT_CONTROL_PATTERNS),
        help="Comma-separated regexes matched against guide/target names",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "write-demo":
        path = write_demo_dragen(args.output_dir)
        print(f"Wrote demo DRAGEN files to {path}")
        return 0
    config = PipelineConfig(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        sample_id=args.sample_id,
        singlet_only=not args.keep_multiplets,
        n_mads=args.n_mads,
        control=args.control,
        replicate_col=args.replicate_col,
        skip_mixscape=args.skip_mixscape,
        force_mixscape=args.force_mixscape,
        mixscape_max_targets=args.mixscape_max_targets,
        de_top_n=args.de_top_n,
        skip_distance=args.skip_distance,
        skip_de=args.skip_de,
        n_perms=args.n_perms,
        perturbation_type=args.perturbation_type,
        control_patterns=tuple(p.strip() for p in args.control_patterns.split(",") if p.strip()),
    )
    report = run_pipeline(config)
    print(f"Finished. Report: {config.output_dir / 'report.json'}")
    print(f"Cells after singlet filter: {report.get('n_after_singlet')}")
    if "mixscape_global" in report:
        print(f"Mixscape: {report['mixscape_global']}")
    if "de_note" in report:
        print(report["de_note"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
