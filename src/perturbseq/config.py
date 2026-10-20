"""Load the run configuration and stop when required customer fields are missing."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

CRISPR_MODES = ("CRISPR_KO", "CRISPRi", "CRISPRa", "other")
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
NO_NT_CONTROL = "no_nt_control"
ZERO_GUIDE_LABEL = "zero_guide"
CONTROL_SOURCE_ZERO = "zero_guide_cells"
CONTROL_SOURCE_DESIGNED = "designed_nt"
TEMPORARY_CONTROL_NOTE = (
    "临时对照，不是设计好的 Non-Targeting。"
    "对照细胞是没有检测到 guide RNA 的 zero-guide cells。"
)
_OTHER_ORGANELLE_SYMBOLS = (
    ("pt-", "叶绿体"),
    ("PT-", "叶绿体"),
    ("Cp-", "叶绿体"),
    ("CP-", "叶绿体"),
)
_OTHER_ORGANELLE_IDS = (("ATCG", "叶绿体"),)


class ConfigError(RuntimeError):
    """The run is not allowed to start."""


_QC_CUTOFFS = ("min_umi", "max_umi", "min_genes", "max_mito_frac")


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


def load_config(path: str | Path) -> dict[str, Any]:
    loaded = read_config(path)
    validate_config(loaded)
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
    for key in (*_QC_CUTOFFS, "thresholds_confirmed", "mito_gene_prefix", "mito_prefix"):
        qc.pop(key, None)
    provided = {key: customer.get(key) for key in _QC_CUTOFFS}
    if any(value is not None for value in provided.values()):
        qc.update(provided)

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


def validate_config(config: dict[str, Any]) -> None:
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
            "labels.flag_guides_above 已停用，不能作为剔除阈值。请使用 labels.guides_per_cell。"
        )
    guides_per_cell = labels.get("guides_per_cell", 1)
    if (
        not isinstance(guides_per_cell, int)
        or isinstance(guides_per_cell, bool)
        or guides_per_cell < 1
    ):
        raise ConfigError(
            "labels.guides_per_cell 必须是正整数。默认 1，表示每个细胞恰好一条 guide。"
        )
    labels["guides_per_cell"] = guides_per_cell

    qc = _section(config, "qc")
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
    statistics = config.setdefault("statistics", {"mode": "descriptive"})
    _validate_statistics(statistics)
    _section(config, "augur")
    _section(config, "plotting")
    seed = config["runtime"].get("seed", 0)
    if not isinstance(seed, int) or isinstance(seed, bool):
        raise ConfigError("随机种子必须是整数。")

    augur = config["augur"]
    classifier = augur.get("classifier", "random_forest")
    if classifier not in {"random_forest", "logistic_regression"}:
        raise ConfigError("augur.classifier 只能是 random_forest 或 logistic_regression。")
    folds = augur.get("n_folds", 5)
    if not isinstance(folds, int) or isinstance(folds, bool) or folds < 2:
        raise ConfigError("augur.n_folds 必须是不小于 2 的整数。")


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


def control_spec(config: dict[str, Any]) -> dict[str, Any]:
    """Resolve the contrast. A designed NT is never replaced by zero-guide cells.

    Without ``labels.nt_control``, cells that have no detected guide RNA are a
    temporary control. The label is ``zero_guide``, not a Non-Targeting guide name.
    """
    designed = nt_control_name(config)
    if designed:
        return {
            "control_source": CONTROL_SOURCE_DESIGNED,
            "control_label": designed,
            "temporary": False,
            "control_phrase": designed,
            "note": "",
        }
    return {
        "control_source": CONTROL_SOURCE_ZERO,
        "control_label": ZERO_GUIDE_LABEL,
        "temporary": True,
        "control_phrase": "zero-guide cells",
        "note": TEMPORARY_CONTROL_NOTE,
    }


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
        return "down", "CRISPRi 预期靶基因表达下降。"
    if mode == "CRISPRa":
        return "up", "CRISPRa 预期靶基因表达升高。"
    if mode == "CRISPR_KO":
        return "not_required", "CRISPR_KO 不以 RNA 表达下降作为唯一成功标准。"
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
        raise ConfigError(
            "experiment.crispr_mode 必须是 KO / CRISPR_KO、i / CRISPRi、a / CRISPRa 或 other。"
        )
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
            raise ConfigError("CRISPRi 的 expected_target_direction 只能是 null 或 down。")
        return "down"
    if mode == "CRISPRa":
        if direction not in (None, "up"):
            raise ConfigError("CRISPRa 的 expected_target_direction 只能是 null 或 up。")
        return "up"
    if mode == "CRISPR_KO":
        if direction not in (None, "not_required"):
            raise ConfigError(
                "CRISPR_KO 的 expected_target_direction 只能是 null 或 not_required。"
            )
        return "not_required"
    if not isinstance(description, str) or not description.strip():
        raise ConfigError("crispr_mode 为 other 时必须填写 crispr_mode_description。")
    if direction not in DIRECTIONS:
        raise ConfigError(
            "crispr_mode 为 other 时 expected_target_direction 必须是 "
            "down、up、bidirectional 或 unknown。"
        )
    return str(direction)


def _require_confirmed_qc(qc: dict[str, Any]) -> None:
    n_mad = qc.get("n_mad", 3)
    if isinstance(n_mad, bool) or not isinstance(n_mad, int | float) or float(n_mad) <= 0:
        raise ConfigError("qc.n_mad 必须是正数。")
    qc["n_mad"] = float(n_mad)
    present = [key for key in _QC_CUTOFFS if qc.get(key) is not None]
    if not present:
        return
    if len(present) != len(_QC_CUTOFFS):
        raise ConfigError(
            "命令行上的 --min-umi、--max-umi、--min-genes、--max-mito-frac 必须同时填写，"
            "或者全部留空，留空时用 median ± qc.n_mad。"
        )
    for key in _QC_CUTOFFS:
        value = qc.get(key)
        if isinstance(value, bool) or not isinstance(value, int | float):
            raise ConfigError(f"qc.{key} 必须是数值。")
    if float(qc["min_umi"]) >= float(qc["max_umi"]):
        raise ConfigError("qc.min_umi 必须小于 qc.max_umi。")
    mito = float(qc["max_mito_frac"])
    if not 0 <= mito <= 1:
        raise ConfigError("qc.max_mito_frac 必须在 0 到 1 之间。")


def _blank_to_none(value: Any) -> str | None:
    text = _optional_text(value)
    return text or None


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


