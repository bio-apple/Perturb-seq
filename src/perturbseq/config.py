"""Load the run configuration and stop when required customer fields are missing."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import yaml

CRISPR_ALIASES = {
    "KO": "CRISPR_KO",
    "CRISPR_KO": "CRISPR_KO",
    "i": "CRISPRi",
    "CRISPRi": "CRISPRi",
    "a": "CRISPRa",
    "CRISPRa": "CRISPRa",
    "other": "other",
}
ORGANISMS = ("human", "mouse")
_ORGANISM_ERROR = "experiment.organism 只能填 human 或 mouse。"
DIRECTIONS = ("down", "up", "bidirectional", "unknown")
POLICIES = ("drop",)
MITO_PREFIX_BY_ORGANISM = {"human": "MT-", "mouse": "mt-"}
ZERO_GUIDE_LABEL = "zero_guide"
# Default GEX marker panels for detection heatmap / dotplot when plotting.markers is unset.
DEFAULT_MARKERS_BY_ORGANISM = {
    "human": (
        "ACTB",
        "GAPDH",
        "B2M",
        "MALAT1",
        "NEAT1",
        "UBC",
        "PPIA",
        "EEF1A1",
    ),
    "mouse": (
        "Actb",
        "Gapdh",
        "B2m",
        "Malat1",
        "Neat1",
        "Ubc",
        "Ppia",
        "Eef1a1",
    ),
}
CONTROL_SOURCE_ZERO = "zero_guide_cells"
CONTROL_SOURCE_DESIGNED = "designed_nt"
# Shared width for preprocess PCA (n_pcs) and Mixscape CalcPerturbSig (n_dims).
DEFAULT_EMBEDDING_NDIMS = 30
DEFAULT_SIGNATURE_NDIMS = DEFAULT_EMBEDDING_NDIMS
TEMPORARY_CONTROL_WARNING = (
    "当前对照是临时 zero-guide cells，这不是设计好的非靶向对照。"
    "zero-guide 细胞可能在 UMI、基因数、线粒体比例、环境 RNA 和细胞状态上与靶向细胞不同。"
    "后续 log2FC、Mixscape、Mixscale、E-distance 都是相对这组临时对照。"
    "UMI、基因数和线粒体比例的对比只作描述，分布差异不会中止分析。"
)
TEMPORARY_CONTROL_NOTE = TEMPORARY_CONTROL_WARNING
_OTHER_ORGANELLE_SYMBOLS = (
    ("pt-", "叶绿体"),
    ("PT-", "叶绿体"),
    ("Cp-", "叶绿体"),
    ("CP-", "叶绿体"),
)
_OTHER_ORGANELLE_IDS = (("ATCG", "叶绿体"),)


class ConfigError(RuntimeError):
    """The run is not allowed to start."""


# UMI absolute cutoffs: both set on the CLI/config, or both omitted (filter skipped).
_QC_CUTOFFS = ("min_umi", "max_umi")
# Gene bounds and mito fraction are software defaults (optional CLI overrides).
_GENE_CUTOFFS = ("min_genes", "max_genes")


def software_config_path() -> Path:
    """Default analysis-parameter file. It does not name a sample or a directory."""
    candidates = [
        Path(__file__).resolve().parents[2] / "config" / "config.yaml",
        Path.cwd() / "config" / "config.yaml",
    ]
    for path in candidates:
        if path.is_file():
            return path
    raise ConfigError("找不到分析软件默认参数文件。请用 --config 指定。")


def read_config(path: str | Path) -> dict[str, Any]:
    """Load a YAML mapping without checking customer fields."""
    config_path = Path(path)
    with config_path.open() as handle:
        loaded = yaml.safe_load(handle)
    if not isinstance(loaded, dict):
        raise ConfigError("配置文件必须是一个 YAML 映射。")
    loaded["_config_path"] = str(config_path.resolve())
    return loaded


def apply_customer_arguments(config: dict[str, Any], customer: dict[str, Any]) -> dict[str, Any]:
    """Replace project fields with the values the customer passed on the command line.

    Paths, sample name, organism, guide separator, and the designed NT name in
    the YAML file are not used.
    """
    config["experiment"] = {
        "sample_id": customer.get("sample_id"),
        "organism": customer.get("organism"),
        "crispr_mode": customer.get("crispr_mode"),
        "crispr_mode_description": customer.get("crispr_mode_description"),
        "expected_target_direction": customer.get("expected_target_direction"),
    }
    inputs = config.setdefault("input", {})
    if not isinstance(inputs, dict):
        inputs = {}
        config["input"] = inputs
    inputs["input_dir"] = customer.get("input_dir")
    inputs["guide_separator"] = customer.get("guide_separator")
    for key in ("matrix", "barcodes", "features", "feature_reference", "assignments"):
        inputs.pop(key, None)

    config["output"] = {"output_dir": customer.get("output_dir")}
    run = config.get("run")
    if isinstance(run, dict):
        run.pop("output_dir", None)

    labels = config.setdefault("labels", {})
    if not isinstance(labels, dict):
        labels = {}
        config["labels"] = labels
    labels["nt_control"] = _optional_text(customer.get("nt_control"))

    qc = config.setdefault("qc", {})
    if not isinstance(qc, dict):
        qc = {}
        config["qc"] = qc
    for key in (
        *_QC_CUTOFFS,
        "thresholds_confirmed",
        "mito_gene_prefix",
        "mito_prefix",
        "max_ribo_frac",
    ):
        qc.pop(key, None)
    provided = {key: customer.get(key) for key in _QC_CUTOFFS}
    if any(value is not None for value in provided.values()):
        qc.update(provided)
    for key in _GENE_CUTOFFS:
        if customer.get(key) is not None:
            qc[key] = customer[key]
    if customer.get("max_mito_frac") is not None:
        qc["max_mito_frac"] = customer["max_mito_frac"]

    cell_type_key = _blank_to_none(customer.get("cell_type_key"))
    batch_key = _blank_to_none(customer.get("batch_key"))
    preprocess = config.setdefault("preprocess", {})
    if isinstance(preprocess, dict):
        preprocess["cell_type_key"] = cell_type_key
        preprocess["batch_key"] = batch_key
    perturb = config.setdefault("perturb", {})
    if isinstance(perturb, dict):
        perturb["cell_type_key"] = cell_type_key
    return config


# Leaf values, including lists, are not walked. ``plotting.palette`` is a color map:
# any name under it is read by lookup, so those names are not unknown keys.
_ANY_KEYS = object()
_ALLOWED_KEYS: dict[str, Any] = {
    "_config_path": None,
    "document_version": None,
    "experiment": {
        "sample_id": None,
        "organism": None,
        "crispr_mode": None,
        "crispr_mode_description": None,
        "expected_target_direction": None,
        "resolved_expected_direction": None,
        "biological_sample_count": None,
        "resolved_organism": None,
        "resolved_mito_prefix": None,
        "species_check": None,
        "organism_evidence": None,
        "other_organelle_prefixes": None,
    },
    "input": {
        "barcode_policy": None,
        "input_dir": None,
        "guide_separator": None,
        "matrix": None,
        "barcodes": None,
        "features": None,
        "feature_reference": None,
        "assignments": None,
        "cell_annotations": None,
    },
    "output": {"output_dir": None},
    "labels": {
        "feature_type_crispr": None,
        "feature_type_gex": None,
        "multi_guide_policy": None,
        "targets_per_cell": None,
        # Rejected in validate_config with a rename message.
        "guides_per_cell": None,
        "nt_control": None,
        "flag_guides_above": None,
    },
    # Filled by assign via resolve_control_spec; single source of truth for contrasts.
    "control": {
        "control_source": None,
        "control_label": None,
        "temporary": None,
        "control_phrase": None,
        "note": None,
        "fallback_from": None,
    },
    "qc": {
        "min_cells_per_target": None,
        "min_cells_per_guide": None,
        "min_guide_umi": None,
        "min_umi": None,
        "max_umi": None,
        "min_genes": None,
        "max_genes": None,
        "max_mito_frac": None,
        # Rejected in validate_config (unsupported key).
        "max_ribo_frac": None,
        "mito_gene_prefix": None,
        "mito_prefix": None,
    },
    "preprocess": {
        "normalize_target_sum": None,
        "n_hvg": None,
        "n_pcs": None,
        "n_neighbors": None,
        "leiden_resolution": None,
        "batch_key": None,
        "cell_type_key": None,
    },
    "mixscape": {
        "enabled": None,
        "min_de_genes": None,
        "logfc_threshold": None,
        "pval_cutoff": None,
        "max_de_genes": None,
        "n_dims": None,
        "n_control_neighbors": None,
        "iter_num": None,
        "cell_type_key": None,
    },
    "mixscale": {
        "min_de_genes": None,
        "logfc_threshold": None,
        "pval_cutoff": None,
        "max_de_genes": None,
    },
    "statistics": {
        "mode": None,
        "aggregation": None,
        "pseudocount_cpm": None,
        "bootstrap_n": None,
        "bootstrap_ci": None,
        "min_nt_cells": None,
        "fdr_threshold": None,
        "allow_random_split_replicates": None,
    },
    "perturb": {
        "mixscape_min_de_genes": None,
        "mixscape_logfc": None,
        "mixscape_pval": None,
        "cell_type_key": None,
    },
    "target_expression": {
        "min_nt_detection_fraction": None,
        "min_nt_mean_cpm": None,
        "direction_log2fc_tolerance": None,
    },
    "singler": {"enabled": None},
    "plotting": {
        "formats": None,
        "mixscape_heatmap_max_genes": None,
        "mixscape_heatmap_max_cells": None,
        "dpi": None,
        "width": None,
        "height": None,
        "theme": None,
        "max_legend_categories": None,
        "scatter_max_points": None,
        "volcano_lfc": None,
        "heatmap_max_genes": None,
        "save_plot_data": None,
        "hist_bins": None,
        "markers": None,
        "palette": _ANY_KEYS,
    },
    "run": {"random_seed": None, "output_dir": None},
    "runtime": {"seed": None},
}


def _reject_unknown_keys(config: dict[str, Any]) -> None:
    """Stop when a mapping contains a key the analysis does not read."""
    unknown = _unknown_key_paths(config, _ALLOWED_KEYS, "")
    if unknown:
        raise ConfigError("未识别的配置键 " + "、".join(unknown) + "。")


def _unknown_key_paths(node: Any, allowed: Any, prefix: str) -> list[str]:
    if not isinstance(node, dict):
        return []
    if allowed is _ANY_KEYS:
        found: list[str] = []
        for key, value in node.items():
            path = f"{prefix}.{key}" if prefix else str(key)
            if isinstance(value, dict):
                found.extend(_unknown_key_paths(value, {}, path))
        return found
    if not isinstance(allowed, dict):
        return [f"{prefix}.{key}" if prefix else str(key) for key in node]
    found = []
    for key, value in node.items():
        path = f"{prefix}.{key}" if prefix else str(key)
        if key not in allowed:
            found.append(path)
            continue
        found.extend(_unknown_key_paths(value, allowed[key], path))
    return found


def validate_config(config: dict[str, Any]) -> None:
    _reject_unknown_keys(config)
    experiment = _section(config, "experiment")
    sample_id = _required_text(experiment, "sample_id", "experiment.sample_id")
    if any(separator in sample_id for separator in ("/", "\\", "\n")):
        raise ConfigError("experiment.sample_id 只能是一个样本名。")

    _organism(experiment)
    mode = _crispr_mode(experiment)
    description = experiment.get("crispr_mode_description")
    direction = experiment.get("expected_target_direction")
    experiment["resolved_expected_direction"] = _resolve_direction(mode, direction, description)
    experiment["biological_sample_count"] = 1

    labels = _section(config, "labels")
    _nt_control(labels)
    policy = labels.get("multi_guide_policy", "drop")
    if policy == "keep_first":
        raise ConfigError("labels.multi_guide_policy=keep_first 已不再支持，必须改为 drop。")
    if policy not in POLICIES:
        raise ConfigError("labels.multi_guide_policy 只能是 drop。")
    if labels.get("flag_guides_above") is not None:
        raise ConfigError(
            "labels.flag_guides_above 已停用，不能作为剔除阈值。请使用 labels.targets_per_cell。"
        )
    if "guides_per_cell" in labels:
        raise ConfigError(
            "labels.guides_per_cell 已改为 labels.targets_per_cell："
            "按每个细胞的唯一靶基因数判定，不再按 guide 条数。"
        )
    targets_per_cell = labels.get("targets_per_cell", 1)
    if (
        not isinstance(targets_per_cell, int)
        or isinstance(targets_per_cell, bool)
        or targets_per_cell < 1
    ):
        raise ConfigError(
            "labels.targets_per_cell 必须是正整数。默认 1，表示每个细胞恰好一个靶基因"
            "（同一靶基因的多条 guide 仍算单靶）。"
        )
    labels["targets_per_cell"] = targets_per_cell

    qc = _section(config, "qc")
    if "max_ribo_frac" in qc:
        raise ConfigError(
            "不支持 qc.max_ribo_frac：请删除该键（及命令行 --max-ribo-frac）。QC 指标为 nCount / nFeature / percent.mt。"
        )
    for key in ("min_cells_per_target", "min_cells_per_guide"):
        value = qc.get(key)
        if not isinstance(value, int) or isinstance(value, bool) or value < 1:
            raise ConfigError(f"qc.{key} 必须是正整数。")
    _require_confirmed_qc(qc)
    mito_prefix(config)

    inputs = _section(config, "input")
    separator = inputs.get("guide_separator")
    if not isinstance(separator, str) or separator == "":
        raise ConfigError("多 guide 分隔符未配置。input.guide_separator 必须是已核实的非空字符串。")
    barcode_policy = inputs.get("barcode_policy", "strict")
    if barcode_policy != "strict":
        raise ConfigError("当前只实现 barcode_policy=strict，不能静默改用交集。")
    _normalize_runtime(config)
    _section(config, "preprocess")
    _section(config, "perturb")
    _section(config, "mixscape")
    _section(config, "mixscale")
    _validate_embedding_ndims(config)
    _validate_mixscape_fallback(config)
    statistics = config.setdefault("statistics", {"mode": "descriptive"})
    _validate_statistics(statistics)
    _section(config, "plotting")
    seed = config["runtime"].get("seed", 0)
    if not isinstance(seed, int) or isinstance(seed, bool):
        raise ConfigError("随机种子必须是整数。")


def resolve_input_paths(config: dict[str, Any]) -> dict[str, Path]:
    """Return the five DRAGEN files for the single configured sample."""
    sample_id = str(config["experiment"]["sample_id"])
    input_cfg = config["input"]
    input_dir = Path(input_cfg.get("input_dir", "."))
    explicit = {
        "matrix": input_cfg.get("matrix"),
        "barcodes": input_cfg.get("barcodes"),
        "features": input_cfg.get("features"),
        "feature_reference": input_cfg.get("feature_reference"),
        "assignments": input_cfg.get("assignments"),
    }
    if any(value is not None for value in explicit.values()) and not all(
        value is not None for value in explicit.values()
    ):
        raise ConfigError("input 中的矩阵、barcode、feature、参考表和分配文件必须同时填写。")

    if explicit["matrix"] is None:
        paths = {
            "matrix": input_dir / f"{sample_id}.scRNA.filtered.matrix.mtx.gz",
            "barcodes": input_dir / f"{sample_id}.scRNA.filtered.barcodes.tsv.gz",
            "features": input_dir / f"{sample_id}.scRNA.filtered.features.tsv.gz",
            "feature_reference": input_dir / f"{sample_id}.scRNA.feature_barcode_reference.csv",
            "assignments": input_dir / f"{sample_id}.scRNA.positive_cell_guide_assignments.csv",
        }
    else:
        paths = {key: Path(str(value)) for key, value in explicit.items()}

    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        raise ConfigError("缺少输入文件:\n" + "\n".join(missing))

    assignment_files = sorted(input_dir.glob("*.positive_cell_guide_assignments.csv"))
    extras = [
        path.name for path in assignment_files if path.resolve() != paths["assignments"].resolve()
    ]
    if extras:
        raise ConfigError("检测到多样本输入，流程终止。目录中还有: " + ", ".join(extras))
    return paths


def nt_control_name(config: dict[str, Any]) -> str:
    """Configured non-targeting label. Empty means this run has no designed NT."""
    labels = config.get("labels")
    if not isinstance(labels, dict):
        return ""
    return _optional_text(labels.get("nt_control"))


def _designed_control_spec(designed: str) -> dict[str, Any]:
    return {
        "control_source": CONTROL_SOURCE_DESIGNED,
        "control_label": designed,
        "temporary": False,
        "control_phrase": designed,
        "note": "",
        "fallback_from": "",
    }


def _temporary_control_spec(*, fallback_from: str = "") -> dict[str, Any]:
    note = TEMPORARY_CONTROL_NOTE
    if fallback_from:
        note = (
            f"{TEMPORARY_CONTROL_NOTE} "
            f"已配置 labels.nt_control={fallback_from!r}，但数据中无可用设计 NT 细胞，"
            "已回退到临时 zero_guide。"
        )
    return {
        "control_source": CONTROL_SOURCE_ZERO,
        "control_label": ZERO_GUIDE_LABEL,
        "temporary": True,
        "control_phrase": "zero-guide cells",
        "note": note,
        "fallback_from": fallback_from,
    }


def bind_control_spec(config: dict[str, Any], spec: dict[str, Any]) -> dict[str, Any]:
    """Persist the resolved contrast so every step reads the same control_spec."""
    block = {
        "control_source": str(spec["control_source"]),
        "control_label": str(spec["control_label"]),
        "temporary": bool(spec["temporary"]),
        "control_phrase": str(spec.get("control_phrase", spec["control_label"])),
        "note": str(spec.get("note", "")),
        "fallback_from": str(spec.get("fallback_from", "")),
    }
    config["control"] = block
    return dict(block)


def control_spec(config: dict[str, Any]) -> dict[str, Any]:
    """Contrast for this run: designed NT when available, else temporary zero_guide.

    After ``assign``, ``config['control']`` is the single source of truth (including
    data-level fallback when ``labels.nt_control`` is set but no NT cells exist).
    Before assign, this is provisional from ``labels.nt_control`` only.
    """
    block = config.get("control")
    if isinstance(block, dict) and block.get("control_label") not in (None, ""):
        return {
            "control_source": str(block.get("control_source", CONTROL_SOURCE_ZERO)),
            "control_label": str(block["control_label"]),
            "temporary": bool(block.get("temporary", False)),
            "control_phrase": str(
                block.get("control_phrase", block["control_label"])
            ),
            "note": str(block.get("note", "")),
            "fallback_from": str(block.get("fallback_from", "")),
        }
    designed = nt_control_name(config)
    if designed:
        return _designed_control_spec(designed)
    return _temporary_control_spec()


def resolve_control_spec(config: dict[str, Any], obs) -> dict[str, Any]:
    """Finalize the contrast from config + assignment ledger (call in assign).

    Designed NT wins only when usable NT cells exist. Otherwise zero-guide cells
    become an explicit temporary control.
    """
    from perturbseq.assign import control_cell_mask

    designed = nt_control_name(config)
    if designed:
        designed_spec = _designed_control_spec(designed)
        mask = np.asarray(control_cell_mask(obs, designed_spec), dtype=bool)
        if "qc_pass" in getattr(obs, "columns", []):
            qc = obs["qc_pass"]
            if getattr(qc, "dtype", None) is not None and str(qc.dtype) == "bool":
                qc_ok = qc.fillna(False).to_numpy(dtype=bool)
            else:
                qc_ok = qc.astype(str).str.lower().isin(["true", "1"]).to_numpy()
            mask = mask & np.asarray(qc_ok, dtype=bool)
        if int(mask.sum()) > 0:
            return bind_control_spec(config, designed_spec)
        return bind_control_spec(config, _temporary_control_spec(fallback_from=designed))
    return bind_control_spec(config, _temporary_control_spec())


def restore_control_spec(config: dict[str, Any], adata) -> dict[str, Any] | None:
    """Reload a previously resolved contrast from ``adata.uns`` (resume path)."""
    uns = getattr(adata, "uns", None)
    if not isinstance(uns, dict):
        return None
    label = uns.get("control_label")
    source = uns.get("control_source")
    if label in (None, "") or source in (None, ""):
        return None
    temporary = uns.get("control_temporary")
    if temporary is None:
        temporary = str(source) == CONTROL_SOURCE_ZERO
    return bind_control_spec(
        config,
        {
            "control_source": str(source),
            "control_label": str(label),
            "temporary": bool(temporary),
            "control_phrase": str(uns.get("control_phrase", label)),
            "note": str(uns.get("control_note", "")),
            "fallback_from": str(uns.get("control_fallback_from", "")),
        },
    )


def verify_species_consistency(
    config: dict[str, Any],
    gene_ids,
    gene_symbols,
) -> dict[str, Any]:
    """Check gene ids and symbols against the organism the user already set.

    The organism is not inferred. Human must not meet ENSMUSG or mt-. Mouse
    must not meet ENSG or MT-. Chloroplast prefixes are recorded and are not
    mitochondrial.
    """
    experiment = config["experiment"]
    organism = str(experiment.get("organism", ""))
    if organism not in MITO_PREFIX_BY_ORGANISM:
        raise ConfigError(_ORGANISM_ERROR)
    prefix = mito_prefix(config)
    ids = [str(gene_id) for gene_id in gene_ids]
    symbols = [str(symbol) for symbol in gene_symbols]
    human_ids, human_id_examples = _matching(ids, "ENSG")
    mouse_ids, mouse_id_examples = _matching(ids, "ENSMUSG")
    human_symbols, human_symbol_examples = _matching(symbols, "MT-")
    mouse_symbols, mouse_symbol_examples = _matching(symbols, "mt-")
    evidence = [
        f"输入物种 {organism}",
        f"线粒体前缀 {prefix}",
        _marker_text("基因号 ENSG", human_ids, human_id_examples),
        _marker_text("基因号 ENSMUSG", mouse_ids, mouse_id_examples),
        _marker_text("基因符号 MT-", human_symbols, human_symbol_examples),
        _marker_text("基因符号 mt-", mouse_symbols, mouse_symbol_examples),
    ]
    problems: list[str] = []
    if human_ids and mouse_ids:
        problems.append(
            "基因号互相矛盾："
            f"{_marker_text('ENSG', human_ids, human_id_examples)}；"
            f"{_marker_text('ENSMUSG', mouse_ids, mouse_id_examples)}"
        )
    elif organism == "human" and mouse_ids:
        problems.append(
            f"{_marker_text('基因号 ENSMUSG', mouse_ids, mouse_id_examples)}，与输入物种 human 不一致"
        )
    elif organism == "mouse" and human_ids:
        problems.append(
            f"{_marker_text('基因号 ENSG', human_ids, human_id_examples)}，与输入物种 mouse 不一致"
        )
    if human_symbols and mouse_symbols:
        problems.append(
            "线粒体符号互相矛盾："
            f"{_marker_text('MT-', human_symbols, human_symbol_examples)}；"
            f"{_marker_text('mt-', mouse_symbols, mouse_symbol_examples)}"
        )
    elif organism == "human" and mouse_symbols:
        problems.append(
            f"{_marker_text('基因符号 mt-', mouse_symbols, mouse_symbol_examples)}，"
            "与输入物种 human 不一致"
        )
    elif organism == "mouse" and human_symbols:
        problems.append(
            f"{_marker_text('基因符号 MT-', human_symbols, human_symbol_examples)}，"
            "与输入物种 mouse 不一致"
        )
    experiment["resolved_organism"] = organism
    experiment["resolved_mito_prefix"] = prefix
    experiment["organism_evidence"] = evidence + problems
    experiment["other_organelle_prefixes"] = _other_organelles(ids, symbols)
    experiment["species_check"] = "failed" if problems else "passed"
    if problems:
        raise ConfigError("物种一致性验证失败：" + "；".join(problems) + "。")
    return config


def default_markers(config: dict[str, Any]) -> list[str]:
    """Marker gene symbols for GEX detection heatmap / dotplot.

    Prefer ``plotting.markers`` when set; otherwise organism-specific housekeeping
    defaults from ``DEFAULT_MARKERS_BY_ORGANISM``.
    """
    configured = config.get("plotting", {}).get("markers")
    if configured is not None:
        if not isinstance(configured, list | tuple) or not configured:
            raise ConfigError("plotting.markers 必须是非空基因符号列表。")
        markers = [str(gene).strip() for gene in configured if str(gene).strip()]
        if not markers:
            raise ConfigError("plotting.markers 必须是非空基因符号列表。")
        return markers
    organism = str(config["experiment"].get("organism", ""))
    defaults = DEFAULT_MARKERS_BY_ORGANISM.get(organism)
    if defaults is None:
        raise ConfigError(_ORGANISM_ERROR)
    return list(defaults)


def mito_prefix(config: dict[str, Any]) -> str:
    """Mitochondrial symbol prefix. Human is MT-, mouse is mt-."""
    organism = str(config["experiment"].get("organism", ""))
    expected = MITO_PREFIX_BY_ORGANISM.get(organism)
    if expected is None:
        raise ConfigError(_ORGANISM_ERROR)
    configured = _optional_text(
        config["qc"].get("mito_gene_prefix", config["qc"].get("mito_prefix"))
    )
    if configured and configured != expected:
        raise ConfigError(f"experiment.organism={organism} 时线粒体基因前缀必须是 {expected}。")
    config["qc"]["mito_gene_prefix"] = expected
    return expected


def expected_expression(config: dict[str, Any]) -> tuple[str, str]:
    """Return the direction code and the sentence used in figure captions."""
    mode = config["experiment"]["crispr_mode"]
    if mode == "CRISPRi":
        return "down", "i 预期靶基因表达下降。"
    if mode == "CRISPRa":
        return "up", "a 预期靶基因表达升高。"
    if mode == "CRISPR_KO":
        return "not_required", "KO 不以 RNA 表达下降作为唯一成功标准。"
    direction = str(config["experiment"]["expected_target_direction"])
    description = str(config["experiment"]["crispr_mode_description"]).strip()
    sentence = {
        "down": "客户预期靶基因表达下降。",
        "up": "客户预期靶基因表达升高。",
        "bidirectional": "客户预期靶基因表达双向变化。",
        "unknown": "客户未指定靶基因表达方向。",
    }[direction]
    return direction, f"{description} {sentence}"


def output_dir(config: dict[str, Any]) -> Path:
    if "output" in config and isinstance(config["output"], dict):
        return Path(config["output"]["output_dir"])
    return Path(config["run"]["output_dir"])


def _organism(experiment: dict[str, Any]) -> str:
    raw = experiment.get("organism")
    if not isinstance(raw, str) or not raw.strip():
        raise ConfigError(_ORGANISM_ERROR)
    normalized = raw.strip().lower()
    if normalized not in ORGANISMS:
        raise ConfigError(_ORGANISM_ERROR)
    experiment["organism"] = normalized
    return normalized


def _matching(values: list[str], prefix: str) -> tuple[int, list[str]]:
    examples: list[str] = []
    count = 0
    for value in values:
        if not value.startswith(prefix):
            continue
        count += 1
        if value not in examples and len(examples) < 3:
            examples.append(value)
    return count, examples


def _marker_text(label: str, count: int, examples: list[str]) -> str:
    if not examples:
        return f"{label} {count} 个"
    return f"{label} {count} 个，例如 {'、'.join(examples)}"


def _other_organelles(gene_ids: list[str], gene_symbols: list[str]) -> list[str]:
    found: list[str] = []
    for prefix, name in _OTHER_ORGANELLE_SYMBOLS:
        if any(symbol.startswith(prefix) for symbol in gene_symbols):
            label = f"{name}基因符号 {prefix}"
            if label not in found:
                found.append(label)
    for prefix, name in _OTHER_ORGANELLE_IDS:
        if any(gene_id.startswith(prefix) for gene_id in gene_ids):
            label = f"{name}基因号 {prefix}"
            if label not in found:
                found.append(label)
    return found


def _crispr_mode(experiment: dict[str, Any]) -> str:
    raw = experiment.get("crispr_mode")
    if raw is None or (isinstance(raw, str) and not raw.strip()):
        mode = "CRISPR_KO"
    elif isinstance(raw, str) and raw.strip() in CRISPR_ALIASES:
        mode = CRISPR_ALIASES[raw.strip()]
    else:
        raise ConfigError("experiment.crispr_mode 必须是 KO、i、a 或 other。")
    experiment["crispr_mode"] = mode
    return mode


def _nt_control(labels: dict[str, Any]) -> str:
    if "nt_control" not in labels or labels.get("nt_control") is None:
        labels["nt_control"] = ""
        return ""
    value = labels.get("nt_control")
    if not isinstance(value, str):
        raise ConfigError("labels.nt_control 必须是字符串或留空。")
    labels["nt_control"] = value.strip()
    return labels["nt_control"]


def _resolve_direction(mode: str, direction: Any, description: Any) -> str:
    if mode == "CRISPRi":
        if direction not in (None, "down"):
            raise ConfigError(
                "experiment.crispr_mode=i 时 expected_target_direction 只能是 null 或 down。"
            )
        return "down"
    if mode == "CRISPRa":
        if direction not in (None, "up"):
            raise ConfigError(
                "experiment.crispr_mode=a 时 expected_target_direction 只能是 null 或 up。"
            )
        return "up"
    if mode == "CRISPR_KO":
        if direction not in (None, "not_required"):
            raise ConfigError(
                "experiment.crispr_mode=KO 时 expected_target_direction 只能是 null 或 not_required。"
            )
        return "not_required"
    if not isinstance(description, str) or not description.strip():
        raise ConfigError("experiment.crispr_mode 为 other 时必须填写 crispr_mode_description。")
    if direction not in DIRECTIONS:
        raise ConfigError(
            "experiment.crispr_mode 为 other 时 expected_target_direction 必须是 "
            "down、up、bidirectional 或 unknown。"
        )
    return str(direction)


def _require_confirmed_qc(qc: dict[str, Any]) -> None:
    _require_gene_cutoffs(qc)
    _require_mito_cutoff(qc)
    present = [key for key in _QC_CUTOFFS if qc.get(key) is not None]
    if not present:
        return
    if len(present) != len(_QC_CUTOFFS):
        raise ConfigError(
            "命令行上的 --min-umi、--max-umi 必须同时填写，"
            "或者全部留空（留空则不做 UMI 绝对阈值过滤）。"
            "线粒体比例上限由 qc.max_mito_frac 单独设置（默认 0.2）。"
            "表达基因数上下限由 qc.min_genes / qc.max_genes 单独设置。"
        )
    for key in _QC_CUTOFFS:
        value = qc.get(key)
        if isinstance(value, bool) or not isinstance(value, int | float):
            raise ConfigError(f"qc.{key} 必须是数值。")
    if float(qc["min_umi"]) >= float(qc["max_umi"]):
        raise ConfigError("qc.min_umi 必须小于 qc.max_umi。")


def _require_gene_cutoffs(qc: dict[str, Any]) -> None:
    """Keep cells with min_genes < n_genes < max_genes (drop ≤ min or ≥ max)."""
    for key in _GENE_CUTOFFS:
        value = qc.get(key)
        if value is None or isinstance(value, bool) or not isinstance(value, int | float):
            raise ConfigError(f"qc.{key} 必须是数值。默认 min_genes=200、max_genes=6000。")
        if float(value) < 0:
            raise ConfigError(f"qc.{key} 不能为负。")
        qc[key] = float(value)
    if float(qc["min_genes"]) >= float(qc["max_genes"]):
        raise ConfigError("qc.min_genes 必须小于 qc.max_genes。")


def _require_mito_cutoff(qc: dict[str, Any]) -> None:
    """Drop cells with mito_frac ≥ max_mito_frac (default 0.2 = 20%)."""
    value = qc.get("max_mito_frac")
    if value is None or isinstance(value, bool) or not isinstance(value, int | float):
        raise ConfigError(
            "qc.max_mito_frac 必须是数值。默认 0.2（排除线粒体 reads ≥ 20% 的细胞）。"
        )
    if not 0 <= float(value) <= 1:
        raise ConfigError("qc.max_mito_frac 必须在 0 到 1 之间。")
    qc["max_mito_frac"] = float(value)


def _blank_to_none(value: Any) -> str | None:
    text = _optional_text(value)
    return text or None


def _validate_embedding_ndims(config: dict[str, Any]) -> None:
    """Keep preprocess.n_pcs and mixscape.n_dims aligned (default 30)."""
    preprocess = config["preprocess"]
    mixscape = config["mixscape"]
    n_pcs = preprocess.get("n_pcs", DEFAULT_EMBEDDING_NDIMS)
    n_dims = mixscape.get("n_dims", DEFAULT_EMBEDDING_NDIMS)
    for name, value in (("preprocess.n_pcs", n_pcs), ("mixscape.n_dims", n_dims)):
        if not isinstance(value, int) or isinstance(value, bool) or value < 2:
            raise ConfigError(f"{name} 必须是 ≥2 的整数。")
    if n_pcs != n_dims:
        raise ConfigError(
            f"preprocess.n_pcs={n_pcs} 与 mixscape.n_dims={n_dims} 不一致；"
            f"二者须相同（默认 {DEFAULT_EMBEDDING_NDIMS}）。"
        )
    preprocess["n_pcs"] = n_pcs
    mixscape["n_dims"] = n_dims


def _validate_mixscape_fallback(config: dict[str, Any]) -> None:
    """Keep mixscape.logfc_threshold and perturb.mixscape_logfc from drifting apart.

    The primary key is ``mixscape.logfc_threshold`` (Seurat default 0.25).
    ``perturb.mixscape_logfc`` is only a fallback when the primary key is absent.
    """
    mix = config.get("mixscape")
    perturb = config.get("perturb")
    if not isinstance(mix, dict) or not isinstance(perturb, dict):
        return
    primary = mix.get("logfc_threshold")
    fallback = perturb.get("mixscape_logfc")
    if primary is None or fallback is None:
        return
    if isinstance(primary, bool) or isinstance(fallback, bool):
        return
    if not isinstance(primary, int | float) or not isinstance(fallback, int | float):
        return
    if float(primary) != float(fallback):
        raise ConfigError(
            "mixscape.logfc_threshold 与 perturb.mixscape_logfc 不一致："
            f"{primary} vs {fallback}。"
            "请只改 mixscape.logfc_threshold，并把回退键改成相同数值"
            "（Mixscape 默认 0.25，勿与 mixscale.logfc_threshold=0.2 混用）。"
        )


def _validate_statistics(statistics: dict[str, Any]) -> None:
    if not isinstance(statistics, dict):
        raise ConfigError("statistics 必须是映射。")
    mode = statistics.get("mode", "descriptive")
    if mode == "guide_exploratory":
        raise ConfigError("unsupported_inference_backend")
    if mode != "descriptive":
        raise ConfigError("statistics.mode 当前只定义了 descriptive。")
    if statistics.get("allow_random_split_replicates") is True:
        raise ConfigError("禁止用随机拆分细胞制造生物学重复。")
    alpha = statistics.get("pseudocount_cpm", 1.0)
    if isinstance(alpha, bool) or not isinstance(alpha, int | float) or float(alpha) <= 0:
        raise ConfigError("statistics.pseudocount_cpm 必须是正数。")
    n_boot = statistics.get("bootstrap_n", 200)
    if isinstance(n_boot, bool) or not isinstance(n_boot, int) or n_boot < 1:
        raise ConfigError("statistics.bootstrap_n 必须是正整数。")
    statistics["bootstrap_n"] = n_boot
    level = statistics.get("bootstrap_ci", 0.95)
    if isinstance(level, bool) or not isinstance(level, int | float) or not (0 < float(level) < 1):
        raise ConfigError("statistics.bootstrap_ci 必须是 (0, 1) 内的数。")
    statistics["bootstrap_ci"] = float(level)


def _normalize_runtime(config: dict[str, Any]) -> None:
    run = config.get("run") if isinstance(config.get("run"), dict) else {}
    output = config.get("output") if isinstance(config.get("output"), dict) else {}
    runtime = config.get("runtime") if isinstance(config.get("runtime"), dict) else {}
    output_path = output.get("output_dir", run.get("output_dir"))
    if not isinstance(output_path, str) or not output_path.strip():
        raise ConfigError("必须填写 output.output_dir 或 run.output_dir。")
    seed = runtime.get("seed", run.get("random_seed", 0))
    config["output"] = {**output, "output_dir": output_path}
    config["runtime"] = {**runtime, "seed": seed}


def _section(config: dict[str, Any], name: str) -> dict[str, Any]:
    section = config.get(name)
    if not isinstance(section, dict):
        raise ConfigError(f"配置缺少 {name} 段。")
    return section


def _required_text(section: dict[str, Any], key: str, label: str) -> str:
    value = section.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ConfigError(f"{label} 为必填。")
    return value.strip()


def _optional_text(value: Any) -> str:
    if not isinstance(value, str):
        return ""
    return value.strip()
