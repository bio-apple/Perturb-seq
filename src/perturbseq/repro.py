"""Reproducibility helpers: checksums, config, stage status, resume, manifest."""

from __future__ import annotations

import hashlib
import json
import logging
import os
import platform
import sys
from collections.abc import Iterable
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml

CHUNK_SIZE = 1 << 20  # 1 MiB

# Align with pipeline.PIPELINE_STAGES (diagram stages 1–5; 3a∥3b).
STAGE_ORDER = (
    "1_input_validation",
    "2_preprocessing_qc",
    "3a_perturbation_modeling",
    "3b_statistical_inference",
    "4_robustness",
    "5_report",
)

MANIFEST_NAME = "run_manifest.json"
LOG_NAME = "run.log"
STATUS_NAME = "status.json"
CHECKPOINT_NAME = "adata.h5ad"


def sha256_file(path: Path | str, chunk_size: int = CHUNK_SIZE) -> str:
    """Stream file content and return hex sha256 digest."""
    path = Path(path)
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def sha256_paths(paths: Iterable[Path | str]) -> dict[str, str]:
    """Checksum existing paths; missing paths map to empty string."""
    out: dict[str, str] = {}
    for path in paths:
        p = Path(path)
        key = str(p)
        out[key] = sha256_file(p) if p.is_file() else ""
    return out


def collect_versions(packages: Iterable[str] | None = None) -> dict[str, str | None]:
    import importlib.metadata

    packages = tuple(
        packages
        or (
            "perturbseq",
            "scanpy",
            "anndata",
            "pandas",
            "numpy",
            "scipy",
            "pertpy",
            "pyyaml",
        )
    )
    versions: dict[str, str | None] = {}
    for pkg in packages:
        try:
            versions[pkg] = importlib.metadata.version(pkg)
        except importlib.metadata.PackageNotFoundError:
            versions[pkg] = None
    return versions


def collect_environment() -> dict[str, Any]:
    return {
        "python": sys.version,
        "executable": sys.executable,
        "platform": platform.platform(),
        "system": platform.system(),
        "machine": platform.machine(),
        "hostname": platform.node(),
        "cwd": os.getcwd(),
        "utc_now": datetime.now(timezone.utc).isoformat(),
    }


def _jsonable(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, tuple):
        return list(value)
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_jsonable(v) for v in value]
    return value


def config_params_dict(raw: dict[str, Any]) -> dict[str, Any]:
    """Normalize a config mapping for manifest / status JSON."""
    return _jsonable(raw)


def load_yaml_config(path: Path | str) -> dict[str, Any]:
    path = Path(path)
    with path.open() as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Config root must be a mapping: {path}")
    return data


def deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    merged = dict(base)
    for key, value in override.items():
        if value is None:
            continue
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def default_config_path() -> Path:
    # src/perturbseq/repro.py → repo root / configs/default.yaml
    return Path(__file__).resolve().parents[2] / "configs" / "default.yaml"


# Comments for PipelineConfig fields. Keep in sync with pipeline.PipelineConfig
# and configs/default.yaml.
_CONFIG_FIELD_COMMENTS: dict[str, str] = {
    "input_dir": "Required path to DRAGEN sample directory (or pass --input-dir)",
    "output_dir": "Required results directory (or pass --output-dir)",
    "sample_id": "Sample prefix used in DRAGEN filenames",
    "singlet_only": "Keep only cells with num_features == 1; CLI --keep-multiplets sets false",
    "n_mads": "MAD multiplier for cell QC outlier filter",
    "min_cells": "Drop genes expressed in fewer than this many cells",
    "n_top_genes": "Highly variable gene count for PCA/neighbors",
    "n_pcs": "Number of PCA components",
    "leiden_resolution": "Leiden clustering resolution (descriptive)",
    "control": "Control / NT label in gene_target / perturbation",
    "replicate_col": "Optional obs column for biological replicates (null = none)",
    "skip_mixscape": "Skip Mixscape classification",
    "skip_cell_annotation": "Skip cell-cycle / cell-state annotation",
    "skip_distance": "Skip E-distance / E-test",
    "skip_de": "Skip differential expression contrasts",
    "n_perms": "E-test permutations",
    "min_cells_per_pert": "Min cells per perturbation for distance / DE grouping",
    "etest_power_min_cells": "Below this n_cells, mark E-test as low_power (default 50)",
    "secondary_distance_metrics": "Secondary metrics after E-distance (e.g. mmd, wasserstein)",
    "mixscape_max_targets": "Skip Mixscape when unique targets exceed this (unless force)",
    "force_mixscape": "Run Mixscape even above mixscape_max_targets",
    "mixscape_mode": "auto | skip | force | subset (CLI --skip/--force still work)",
    "mixscape_targets": "Target list when mixscape_mode=subset",
    "mixscape_top_n": "Optional top-N targets for Mixscape subset mode",
    "de_top_n": "Top DE genes retained per contrast in summaries",
    "n_jobs": "Parallel workers for per-group DE / E-test (1 = sequential; -1 = all CPUs)",
    "perturbation_type": "Mixscape label suffix, e.g. KO / KD / perturbation",
    "random_state": "Seed for PCA / neighbors / UMAP / etc.",
    "control_patterns": "Regexes matched against guide/target names to label NT/control",
    "guide_merge": "Multi-guide summary: none | equal | umi | confidence | umi_confidence",
}


def _yaml_block(key: str, value: Any) -> str:
    if isinstance(value, tuple):
        value = list(value)
    return yaml.dump({key: value}, default_flow_style=False, sort_keys=False).rstrip()


def generate_config_template() -> str:
    """Fully commented YAML template from default.yaml + PipelineConfig fields."""
    from dataclasses import MISSING, fields

    from perturbseq.pipeline import PipelineConfig

    loaded: dict[str, Any] = {}
    defaults_path = default_config_path()
    if defaults_path.is_file():
        loaded = load_yaml_config(defaults_path)

    lines = [
        "# Fully commented template for the tertiary Perturb-seq pipeline.",
        "# Generate: python -m perturbseq config --generate",
        "#           python -m perturbseq config generate -o my.yaml",
        "# Precedence: CLI flags override YAML. See docs/CLI_YAML.md.",
        "#",
        "# Required: set input_dir and output_dir here and/or via CLI.",
        "",
    ]

    field_order = [f.name for f in fields(PipelineConfig) if f.name != "extra"]
    # Prefer paths first for readability.
    ordered = [k for k in ("input_dir", "output_dir") if k in field_order]
    ordered += [k for k in field_order if k not in ordered]

    for key in ordered:
        lines.append(f"# {_CONFIG_FIELD_COMMENTS.get(key, key)}")
        if key in loaded:
            lines.append(_yaml_block(key, loaded[key]))
        elif key == "input_dir":
            lines.append("# input_dir: data/demo")
        elif key == "output_dir":
            lines.append("# output_dir: results/demo")
        else:
            field_obj = next(f for f in fields(PipelineConfig) if f.name == key)
            if field_obj.default is not MISSING:
                val = field_obj.default
            elif field_obj.default_factory is not MISSING:  # type: ignore[comparison-overlap]
                val = field_obj.default_factory()
            else:
                val = None
            lines.append(_yaml_block(key, val))
        lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def write_config_template(path: Path | str) -> Path:
    """Write generate_config_template() to ``path``; return the path."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(generate_config_template(), encoding="utf-8")
    return path


class RunTracker:
    """Persist per-stage status, checkpoints, run log, and machine-readable manifest."""

    def __init__(self, output_dir: Path | str, *, resume: bool = False, seed: int = 0) -> None:
        self.output_dir = Path(output_dir)
        self.resume = resume
        self.seed = seed
        self.stages_dir = self.output_dir / "stages"
        self.manifest_path = self.output_dir / MANIFEST_NAME
        self.log_path = self.output_dir / LOG_NAME
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.stages_dir.mkdir(parents=True, exist_ok=True)
        self._setup_logging()
        self.manifest: dict[str, Any] = self._load_or_init_manifest()

    def _setup_logging(self) -> None:
        self.logger = logging.getLogger(f"perturbseq.run.{id(self)}")
        self.logger.setLevel(logging.INFO)
        self.logger.handlers.clear()
        self.logger.propagate = False
        formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
        file_handler = logging.FileHandler(self.log_path, mode="a", encoding="utf-8")
        file_handler.setFormatter(formatter)
        stream_handler = logging.StreamHandler(sys.stdout)
        stream_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)
        self.logger.addHandler(stream_handler)

    def _load_or_init_manifest(self) -> dict[str, Any]:
        if self.resume and self.manifest_path.is_file():
            return json.loads(self.manifest_path.read_text())
        return {
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "resume": self.resume,
            "random_seed": self.seed,
            "environment": collect_environment(),
            "versions": collect_versions(),
            "inputs": {},
            "params": {},
            "stages": {},
            "status": "running",
        }

    def set_run_metadata(
        self,
        *,
        inputs: dict[str, str] | None = None,
        params: dict[str, Any] | None = None,
        config_path: str | None = None,
    ) -> None:
        self.manifest["updated_utc"] = datetime.now(timezone.utc).isoformat()
        self.manifest["resume"] = self.resume
        self.manifest["random_seed"] = self.seed
        self.manifest["environment"] = collect_environment()
        self.manifest["versions"] = collect_versions()
        if inputs is not None:
            self.manifest["inputs"] = inputs
        if params is not None:
            self.manifest["params"] = params
        if config_path is not None:
            self.manifest["config_path"] = config_path
            cfg = Path(config_path)
            if cfg.is_file():
                self.manifest["config_sha256"] = sha256_file(cfg)
        self.write_manifest()

    def stage_dir(self, stage: str) -> Path:
        path = self.stages_dir / stage
        path.mkdir(parents=True, exist_ok=True)
        return path

    def status_path(self, stage: str) -> Path:
        return self.stage_dir(stage) / STATUS_NAME

    def checkpoint_path(self, stage: str) -> Path:
        return self.stage_dir(stage) / CHECKPOINT_NAME

    def load_status(self, stage: str) -> dict[str, Any] | None:
        path = self.status_path(stage)
        if not path.is_file():
            return None
        return json.loads(path.read_text())

    def write_status(self, stage: str, payload: dict[str, Any]) -> None:
        path = self.status_path(stage)
        path.write_text(json.dumps(_jsonable(payload), indent=2, default=str))
        self.manifest["stages"][stage] = {
            "status": payload.get("status"),
            "status_path": str(path),
            "started_utc": payload.get("started_utc"),
            "finished_utc": payload.get("finished_utc"),
            "error": payload.get("error"),
        }
        self.write_manifest()

    def write_manifest(self) -> None:
        self.manifest["updated_utc"] = datetime.now(timezone.utc).isoformat()
        self.manifest_path.write_text(json.dumps(_jsonable(self.manifest), indent=2, default=str))

    def invalidate_from(self, stage: str) -> None:
        """Remove status/checkpoints for stage and all downstream stages."""
        if stage not in STAGE_ORDER:
            raise ValueError(f"Unknown stage: {stage}")
        start = STAGE_ORDER.index(stage)
        for name in STAGE_ORDER[start:]:
            sdir = self.stages_dir / name
            for artifact in (STATUS_NAME, CHECKPOINT_NAME):
                path = sdir / artifact
                if path.exists():
                    path.unlink()
            self.manifest.get("stages", {}).pop(name, None)
            self.logger.info("Invalidated stage %s", name)
        self.write_manifest()

    def should_skip(
        self,
        stage: str,
        *,
        params: dict[str, Any],
        input_checksums: dict[str, str],
    ) -> bool:
        if not self.resume:
            return False
        status = self.load_status(stage)
        if status is None:
            return False
        if status.get("status") not in {"success", "skipped"}:
            return False
        prev_params = status.get("params") or {}
        prev_inputs = status.get("input_checksums") or {}
        if prev_params != params or prev_inputs != input_checksums:
            self.logger.info(
                "Stage %s inputs/params changed; invalidating from here",
                stage,
            )
            self.invalidate_from(stage)
            return False
        outputs = status.get("outputs") or {}
        for out_path in outputs.values():
            if out_path and not Path(out_path).exists():
                self.logger.info("Stage %s output missing (%s); will rerun", stage, out_path)
                self.invalidate_from(stage)
                return False
        self.logger.info("Resuming: skip stage %s (unchanged)", stage)
        return True

    def last_checkpoint_before(self, stage: str) -> Path | None:
        """Return newest existing checkpoint among stages strictly before `stage`."""
        if stage not in STAGE_ORDER:
            return None
        idx = STAGE_ORDER.index(stage)
        for name in reversed(STAGE_ORDER[:idx]):
            path = self.checkpoint_path(name)
            if path.is_file():
                return path
        return None

    def begin_stage(
        self,
        stage: str,
        *,
        inputs: dict[str, str],
        params: dict[str, Any],
    ) -> dict[str, Any]:
        payload = {
            "stage": stage,
            "status": "running",
            "started_utc": datetime.now(timezone.utc).isoformat(),
            "finished_utc": None,
            "inputs": list(inputs.keys()),
            "input_checksums": inputs,
            "outputs": {},
            "output_checksums": {},
            "params": params,
            "error": None,
        }
        self.write_status(stage, payload)
        self.logger.info("START stage=%s", stage)
        return payload

    def finish_stage(
        self,
        stage: str,
        payload: dict[str, Any],
        *,
        status: str,
        outputs: dict[str, str] | None = None,
        error: str | None = None,
        extras: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        outputs = outputs or {}
        payload["status"] = status
        payload["finished_utc"] = datetime.now(timezone.utc).isoformat()
        payload["outputs"] = outputs
        payload["output_checksums"] = {
            key: (sha256_file(path) if path and Path(path).is_file() else "")
            for key, path in outputs.items()
        }
        payload["error"] = error
        if extras:
            payload["extras"] = extras
        self.write_status(stage, payload)
        self.logger.info("END stage=%s status=%s", stage, status)
        return payload

    def finalize(self, *, status: str = "success", report_path: str | None = None) -> None:
        self.manifest["status"] = status
        if report_path:
            self.manifest["report_path"] = report_path
        self.write_manifest()
        self.logger.info("Run finished status=%s", status)
