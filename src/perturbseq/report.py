"""Report writing: provenance helpers + per-perturbation analysis HTML/JSON."""

from __future__ import annotations

import html
import json
import math
import re
from dataclasses import asdict
from pathlib import Path
from typing import Any

import pandas as pd
from anndata import AnnData

from perturbseq.io import write_h5ad

CONTROL_DEFAULT = "NT"
SKIP_TARGETS = frozenset({"unassigned", ""})
NA = "N/A"

# Machine-safe analysis status blocks in report.json (always present for optional stages).
ANALYSIS_STATUS_KEYS = ("mixscape", "edistance", "de", "cell_annotation", "perturbation_clusters")

# Core checklist mirrored from docs/STATISTICAL_CAVEATS.md (kept in-package for report HTML).
CORE_STATISTICAL_CAVEATS: tuple[str, ...] = (
    "Cells are not biological replicates — cell-level tests inflate significance; prefer pseudobulk + PyDESeq2 (even without bio reps stays exploratory).",
    "Mixscape NP ≠ proven biological null — NP is transcriptomic resemblance to controls, not proof the guide failed.",
    "E-distance / E-test quantify multivariate shift, not mechanism; embedding choice changes ranks.",
    "Wilcoxon / no-bio-rep pseudobulk DE is exploratory only — not for FDR-style population claims (Squair 2021).",
    "UMAP / Leiden describe dataset structure, not perturbation-effect evidence.",
    "Filtering (QC / singlet / Mixscape KO) changes composition — check composition_audit.csv; prefer perturbation-aware QC.",
    "Guide inconsistency: do not pool discordant same-gene guides into one KO conclusion without review.",
    "Guide assignment is inherited from DRAGEN — this pipeline does not re-call guides by default.",
)


def skipped_status(reason: str, detail: str | None = None, **extra: Any) -> dict[str, Any]:
    """Explicit skip/degrade payload for report.json (never bare null / missing key)."""
    payload: dict[str, Any] = {"skipped": True, "reason": reason}
    if detail is not None:
        payload["detail"] = detail
    if extra:
        payload.update(extra)
    return payload


def completed_status(**extra: Any) -> dict[str, Any]:
    """Success payload: skipped=false + optional rich fields."""
    return {"skipped": False, "reason": None, **extra}


def status_note(block: Any) -> str | None:
    """Human-readable one-liner from an analysis status object."""
    if not isinstance(block, dict):
        return None
    detail = block.get("detail")
    reason = block.get("reason")
    if block.get("skipped"):
        if detail and reason:
            note = f"{reason}: {detail}"
        else:
            note = str(detail or reason or "skipped")
        est = block.get("estimate")
        if isinstance(est, dict) and est.get("note"):
            note = f"{note} | estimate: {est['note']}"
        return note
    return str(detail) if detail else None


def mixscape_kd_caveat(perturbation_type: str | None) -> str | None:
    """Caveat when Mixscape labels are used for knockdown / CRISPRi/a (weaker than KO)."""
    if perturbation_type is None:
        return None
    token = str(perturbation_type).strip().upper()
    if token in {"KO", "KNOCKOUT"}:
        return None
    if token in {"KD", "CRISPRI", "CRISPRA"} or "KD" in token or "CRISPRI" in token or "CRISPRA" in token:
        return (
            "Mixscape assumptions are weaker for knockdown / CRISPRi/a than for CRISPR KO: "
            "partial repression yields graded (not binary) transcriptomic shifts, so NP vs "
            f"perturbed ({perturbation_type}) labels are less decisive. Prefer E-distance / "
            "target-gene knockdown and multi-guide consistency; see docs/STATISTICAL_CAVEATS.md."
        )
    return None


# ---------------------------------------------------------------------------
# Provenance / pipeline helpers (shared with resume / finalize paths)
# ---------------------------------------------------------------------------


def collect_package_versions(packages: tuple[str, ...] | None = None) -> dict[str, str | None]:
    import importlib.metadata

    packages = packages or ("perturbseq", "scanpy", "anndata", "pandas", "numpy", "pertpy")
    versions: dict[str, str | None] = {}
    for pkg in packages:
        try:
            versions[pkg] = importlib.metadata.version(pkg)
        except importlib.metadata.PackageNotFoundError:
            versions[pkg] = None
    return versions


def init_report(sample_id: str, config: Any | None = None) -> dict:
    report: dict = {
        "sample_id": sample_id,
        "versions": collect_package_versions(),
        "steps": [],
        "stages": [],
    }
    if config is not None:
        report["config"] = {
            k: (str(v) if isinstance(v, Path) else v) for k, v in asdict(config).items() if k != "extra"
        }
    return report


def write_report_json(report: dict, output_dir: Path, filename: str = "report.json") -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / filename
    path.write_text(json.dumps(report, indent=2, default=str))
    return path


def write_provenance_manifest(
    output_dir: Path,
    report: dict,
    *,
    inputs: dict[str, str] | None = None,
) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "provenance_manifest.json"
    payload = {
        "sample_id": report.get("sample_id"),
        "versions": report.get("versions"),
        "config": report.get("config"),
        "stages": report.get("stages"),
        "steps": report.get("steps"),
        "confidence_flags": report.get("confidence_flags"),
        "inputs": inputs or {},
        "outputs": {
            "h5ad": report.get("output_h5ad"),
            "report_json": str(output_dir / "report.json"),
            "report_html": str(output_dir / "report.html"),
            "summary_table_csv": str(output_dir / "summary_table.csv"),
            "summary_table_md": str(output_dir / "summary_table.md"),
        },
    }
    path.write_text(json.dumps(payload, indent=2, default=str))
    return path


# ---------------------------------------------------------------------------
# Table loading & per-perturbation summaries
# ---------------------------------------------------------------------------


def _read_csv(path: Path, **kwargs) -> pd.DataFrame | None:
    if not path.exists() or path.stat().st_size == 0:
        return None
    try:
        return pd.read_csv(path, **kwargs)
    except pd.errors.EmptyDataError:
        return None


def _safe_float(value: Any) -> float | None:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(out) or math.isinf(out):
        return None
    return out


def _safe_int(value: Any) -> int | None:
    f = _safe_float(value)
    return int(f) if f is not None else None


def _na(value: Any = None) -> Any:
    return NA if value is None else value


def _perturbation_from_de_file(path: Path, df: pd.DataFrame) -> str | None:
    if "contrast" in df.columns and len(df):
        contrast = str(df["contrast"].iloc[0])
        m = re.match(r"(.+)_vs_", contrast)
        if m:
            return m.group(1)
    stem = path.stem.removeprefix("de_")
    return stem if stem else None


def _index_lookup(df: pd.DataFrame | None, key: str) -> pd.Series | None:
    if df is None or df.empty:
        return None
    if key in df.index:
        return df.loc[key]
    return None


def _row_lookup(df: pd.DataFrame | None, col: str, key: str) -> pd.Series | None:
    if df is None or df.empty or col not in df.columns:
        return None
    hit = df[df[col].astype(str) == key]
    if hit.empty:
        return None
    return hit.iloc[0]


def load_result_tables(output_dir: Path) -> dict[str, Any]:
    """Load optional result tables from an existing pipeline output directory."""
    tables = output_dir / "tables"
    report_path = output_dir / "report.json"
    report = json.loads(report_path.read_text()) if report_path.exists() else {}

    edistance = _read_csv(tables / "edistance.csv")
    if edistance is not None and "gene_target" not in edistance.columns and edistance.shape[1] >= 1:
        edistance = edistance.rename(columns={edistance.columns[0]: "gene_target"})

    etest = _read_csv(tables / "etest.csv", index_col=0)
    distances = _read_csv(tables / "distances.csv")
    if distances is not None and "gene_target" not in distances.columns and distances.shape[1] >= 1:
        distances = distances.rename(columns={distances.columns[0]: "gene_target"})
    guide_qc = _read_csv(tables / "guide_qc.csv")
    guide_consistency = _read_csv(tables / "gene_guide_consistency.csv")
    qc_warnings = _read_csv(tables / "qc_warnings.csv")
    guide_weighted = _read_csv(tables / "gene_guide_weighted.csv")
    pert_clusters = _read_csv(tables / "perturbation_clusters.csv")
    cell_ann = _read_csv(tables / "cell_annotations.csv")
    composition = _read_csv(tables / "composition_by_perturbation.csv")
    if composition is None:
        composition = _read_csv(tables / "perturbation_composition.csv")
    audit = _read_csv(tables / "composition_audit.csv")
    if audit is None:
        audit = _read_csv(tables / "analysis_audit.csv")

    de_by_pert: dict[str, pd.DataFrame] = {}
    if tables.exists():
        for path in sorted(tables.glob("de_*.csv")):
            if path.name == "de_top50_concat.csv":
                continue
            df = pd.read_csv(path)
            if df.empty:
                continue
            pert = _perturbation_from_de_file(path, df)
            if pert:
                de_by_pert[pert] = df

    return {
        "report": report,
        "edistance": edistance,
        "etest": etest,
        "distances": distances,
        "guide_qc": guide_qc,
        "guide_consistency": guide_consistency,
        "qc_warnings": qc_warnings,
        "guide_weighted": guide_weighted,
        "pert_clusters": pert_clusters,
        "cell_annotations": cell_ann,
        "composition": composition,
        "audit": audit,
        "de_by_pert": de_by_pert,
        "figures_dir": output_dir / "figures",
    }


def list_perturbations(
    tables: dict[str, Any],
    control: str = CONTROL_DEFAULT,
    min_cells: int = 10,
) -> list[str]:
    """Union of perturbation IDs from analysis tables (prefer post-QC evidence)."""
    names: set[str] = set()
    ed = tables.get("edistance")
    if ed is not None and "gene_target" in ed.columns:
        names.update(ed["gene_target"].astype(str).tolist())
    et = tables.get("etest")
    if et is not None:
        names.update(map(str, et.index.tolist()))
    for pert in tables.get("de_by_pert") or {}:
        names.add(pert)
    gc = tables.get("guide_consistency")
    if gc is not None and "gene_target" in gc.columns:
        names.update(gc["gene_target"].astype(str).tolist())
    gq = tables.get("guide_qc")
    if gq is not None and "gene_target" in gq.columns:
        if "is_control" in gq.columns:
            mask = ~gq["is_control"].astype(bool)
            names.update(gq.loc[mask, "gene_target"].astype(str).tolist())
        else:
            names.update(gq["gene_target"].astype(str).tolist())

    if not names:
        report = tables.get("report") or {}
        for key, n in (report.get("gene_target_counts") or {}).items():
            if key not in SKIP_TARGETS and key != control and int(n) >= min_cells:
                names.add(str(key))

    names.discard(control)
    names -= SKIP_TARGETS
    return sorted(names)


def _qc_status(n_cells: int | None, min_cells: int, guide_status: str | None) -> str:
    if guide_status:
        return str(guide_status)
    if n_cells is None:
        return NA
    if n_cells < min_cells:
        return f"low_n (<{min_cells} cells)"
    return "ok"


def _top_de_genes(de: pd.DataFrame | None, n: int = 10) -> list[dict[str, Any]]:
    if de is None or de.empty:
        return []
    df = de.copy()
    gene_col = next((c for c in ("names", "gene", "gene_id", "symbol") if c in df.columns), None)
    if gene_col is None:
        return []
    lfc_col = next(
        (c for c in ("logfoldchanges", "log2FoldChange", "lfc", "logFC") if c in df.columns),
        None,
    )
    padj_col = next((c for c in ("pvals_adj", "padj", "p_adj", "FDR") if c in df.columns), None)
    pval_col = next((c for c in ("pvals", "pvalue", "p_value") if c in df.columns), None)
    if padj_col:
        ascending = [True, False] if lfc_col else [True]
        cols = [padj_col, lfc_col] if lfc_col else [padj_col]
        df = df.sort_values(cols, ascending=ascending)
    elif lfc_col:
        df = df.reindex(df[lfc_col].abs().sort_values(ascending=False).index)
    out: list[dict[str, Any]] = []
    for _, row in df.head(n).iterrows():
        entry: dict[str, Any] = {"gene": str(row[gene_col])}
        if lfc_col:
            entry["logfoldchange"] = _safe_float(row[lfc_col])
        if padj_col:
            entry["pvalue_adj"] = _safe_float(row[padj_col])
        if pval_col:
            entry["pvalue"] = _safe_float(row[pval_col])
        if "lfcSE" in df.columns and lfc_col and pd.notna(row.get("lfcSE")):
            lfc = _safe_float(row[lfc_col])
            se = _safe_float(row["lfcSE"])
            if lfc is not None and se is not None:
                entry["ci_low"] = lfc - 1.96 * se
                entry["ci_high"] = lfc + 1.96 * se
        out.append(entry)
    return out


def _composition_for_pert(
    cell_ann: pd.DataFrame | None,
    composition: pd.DataFrame | None,
    pert: str,
    control: str,
    cache: dict[str, dict[str, Any]] | None = None,
) -> dict[str, Any]:
    if composition is not None and not composition.empty:
        row = _row_lookup(composition, "gene_target", pert) or _row_lookup(composition, "perturbation", pert)
        if row is not None:
            return {
                "source": "composition_table",
                "values": {k: (_safe_float(v) if _safe_float(v) is not None else v) for k, v in row.items()},
            }

    if cache is not None:
        if pert in cache:
            return cache[pert]
        return {"source": NA, "note": "No cell_annotations.csv or composition table"}

    if cell_ann is None or cell_ann.empty or "gene_target" not in cell_ann.columns:
        return {"source": NA, "note": "No cell_annotations.csv or composition table"}

    sub = cell_ann[cell_ann["gene_target"].astype(str) == pert]
    if sub.empty:
        return {"source": "cell_annotations", "n_cells": 0, "note": "No annotated cells for this perturbation"}

    result: dict[str, Any] = {"source": "cell_annotations", "n_cells": int(len(sub))}
    for col in ("cell_state", "phase", "cluster_annotation", "leiden"):
        if col in sub.columns:
            result[f"{col}_counts"] = sub[col].astype(str).value_counts().head(8).to_dict()
    ctrl = cell_ann[cell_ann["gene_target"].astype(str) == control]
    if not ctrl.empty and "cell_state" in cell_ann.columns:
        pert_frac = sub["cell_state"].astype(str).value_counts(normalize=True)
        ctrl_frac = ctrl["cell_state"].astype(str).value_counts(normalize=True)
        delta = (pert_frac - ctrl_frac).dropna().sort_values(key=abs, ascending=False).head(5)
        result["cell_state_delta_vs_control"] = {k: round(float(v), 4) for k, v in delta.items()}
    return result


def summarize_perturbation(
    pert: str,
    tables: dict[str, Any],
    *,
    control: str = CONTROL_DEFAULT,
    min_cells: int = 10,
    top_genes: int = 10,
) -> dict[str, Any]:
    """Build the five-question summary for one perturbation."""
    report = tables.get("report") or {}
    config = report.get("config") or {}
    replicate_col = config.get("replicate_col")
    gene_counts = report.get("gene_target_counts") or {}

    ed = tables.get("edistance")
    ed_row = _row_lookup(ed, "gene_target", pert) if ed is not None else None
    n_cells = _safe_int(ed_row["n_cells"]) if ed_row is not None and "n_cells" in ed_row.index else None
    if n_cells is None and pert in gene_counts:
        n_cells = _safe_int(gene_counts[pert])

    gc_row = _row_lookup(tables.get("guide_consistency"), "gene_target", pert)
    n_guides = _safe_int(gc_row["n_guides"]) if gc_row is not None and "n_guides" in gc_row.index else None
    guide_status = None
    if gc_row is not None and "potential_low_efficiency" in gc_row.index:
        ple = gc_row["potential_low_efficiency"]
        if pd.notna(ple) and bool(ple):
            guide_status = "potential_low_efficiency_guides"
    if guide_status is None and gc_row is not None and "guides_consistent" in gc_row.index:
        cons_val = gc_row["guides_consistent"]
        if pd.isna(cons_val):
            guide_status = "consistency_unknown"
        elif bool(cons_val):
            guide_status = "guides_consistent"
        else:
            guide_status = "guides_inconsistent"

    gq = tables.get("guide_qc")
    if n_guides is None and gq is not None and "gene_target" in gq.columns:
        n_from_gq = int((gq["gene_target"].astype(str) == pert).sum())
        if n_from_gq > 0:
            n_guides = n_from_gq

    if replicate_col is None:
        n_replicates_display: Any = 1
    else:
        n_replicates_display = NA
        if gc_row is not None:
            for col in ("n_replicates", "n_samples", "n_reps"):
                if col in gc_row.index:
                    n_replicates_display = _na(_safe_int(gc_row[col]))
                    break

    qc_status = _qc_status(n_cells, min_cells, guide_status)

    et_row = _index_lookup(tables.get("etest"), pert)
    edistance = _safe_float(ed_row["edistance"]) if ed_row is not None and "edistance" in ed_row.index else None
    low_power = None
    if et_row is not None and "low_power" in et_row.index:
        low_power = bool(et_row["low_power"])
    elif ed_row is not None and "n_cells" in ed_row.index:
        power_min = int(
            (report.get("config") or {}).get("etest_power_min_cells")
            or (report.get("edistance") or {}).get("etest_power_min_cells")
            or 50
        )
        n_for_power = _safe_int(ed_row["n_cells"])
        if n_for_power is not None:
            low_power = n_for_power < power_min
    effect: dict[str, Any] = {
        "metric": "edistance" if edistance is not None else NA,
        "effect_size": _na(edistance),
        "confidence_interval": NA,
        "pvalue": NA,
        "pvalue_adj": NA,
        "significant_adj": NA,
        "significant_adj_reported": NA,
        "low_power": _na(low_power) if low_power is not None else NA,
        "note": None,
    }
    # Prefer E-distance bootstrap CI when present (matches metric=edistance).
    ci_src = None
    if ed_row is not None and "edistance_ci_low" in ed_row.index:
        ci_src = ed_row
    elif et_row is not None and "edistance_ci_low" in et_row.index:
        ci_src = et_row
    if ci_src is not None and "edistance_ci_high" in ci_src.index:
        ci_lo = _safe_float(ci_src.get("edistance_ci_low"))
        ci_hi = _safe_float(ci_src.get("edistance_ci_high"))
        if ci_lo is not None and ci_hi is not None:
            ci_level = _safe_float(ci_src.get("ci_level")) if "ci_level" in ci_src.index else 0.95
            n_boot = _safe_int(ci_src.get("n_bootstrap")) if "n_bootstrap" in ci_src.index else None
            effect["confidence_interval"] = {
                "low": ci_lo,
                "high": ci_hi,
                "level": ci_level if ci_level is not None else 0.95,
                "source": "edistance_bootstrap_cells",
                "n_bootstrap": _na(n_boot),
            }
    dist_row = _row_lookup(tables.get("distances"), "gene_target", pert)
    secondary_metrics: dict[str, Any] = {}
    if dist_row is not None:
        for col in dist_row.index:
            if col in {
                "gene_target",
                "edistance",
                "n_cells",
                "embedding",
                "pca_source",
                "metric",
                "index",
                "edistance_ci_low",
                "edistance_ci_high",
                "n_bootstrap",
                "ci_level",
            }:
                continue
            val = _safe_float(dist_row[col])
            if val is not None:
                secondary_metrics[str(col)] = val
    effect["secondary_metrics"] = secondary_metrics if secondary_metrics else NA
    if et_row is not None:
        effect["pvalue"] = _na(_safe_float(et_row.get("pvalue")))
        effect["pvalue_adj"] = _na(_safe_float(et_row.get("pvalue_adj")))
        if "significant_adj" in et_row.index:
            effect["significant_adj"] = bool(et_row["significant_adj"])
        elif isinstance(effect["pvalue_adj"], float):
            effect["significant_adj"] = effect["pvalue_adj"] < 0.05
        if "significant_adj_reported" in et_row.index:
            effect["significant_adj_reported"] = bool(et_row["significant_adj_reported"])
        elif isinstance(effect["significant_adj"], bool):
            effect["significant_adj_reported"] = bool(effect["significant_adj"]) and not bool(
                low_power
            )
        if "low_power" in et_row.index:
            effect["low_power"] = bool(et_row["low_power"])
        if edistance is None and "distance" in et_row.index:
            effect["effect_size"] = _na(_safe_float(et_row["distance"]))
            effect["metric"] = "edistance"
        if low_power:
            note = "low_power: E-test significance not reported (n_cells below etest_power_min_cells)"
            effect["note"] = f"{effect['note']}; {note}" if effect.get("note") else note
    elif edistance is not None:
        effect["note"] = "E-distance available; E-test not run or not tested for this target"

    de = (tables.get("de_by_pert") or {}).get(pert)
    de_padj_min = None
    de_method = None
    de_evidence_level = None
    if de is not None and not de.empty:
        de_method = str(de["method"].iloc[0]) if "method" in de.columns else NA
        if "evidence_level" in de.columns:
            de_evidence_level = str(de["evidence_level"].iloc[0])
        elif de_method in {
            "wilcoxon_cell_level_exploratory",
            "pydeseq2_pseudobulk_no_bio_reps",
        }:
            de_evidence_level = "exploratory"
        elif de_method == "pydeseq2_pseudobulk":
            de_evidence_level = "inferential"
        top1 = _top_de_genes(de, n=1)
        # Only fall back to DESeq2 gene-level CI when E-distance bootstrap CI is absent.
        if effect["confidence_interval"] == NA and top1 and "ci_low" in top1[0]:
            effect["confidence_interval"] = {
                "low": top1[0]["ci_low"],
                "high": top1[0]["ci_high"],
                "level": 0.95,
                "source": "DESeq2 lfcSE (top DE gene)",
            }
        padj_col = next((c for c in ("pvals_adj", "padj") if c in de.columns), None)
        if padj_col:
            de_padj_min = _safe_float(de[padj_col].min())
        if effect["confidence_interval"] == NA and (
            replicate_col is None
            or de_method
            in {"wilcoxon_cell_level_exploratory", "pydeseq2_pseudobulk_no_bio_reps"}
        ):
            note = "No biological-replicate CI (exploratory DE path)"
            effect["note"] = f"{effect['note']}; {note}" if effect.get("note") else note
    if de_evidence_level is None:
        de_block_el = (report.get("de") if isinstance(report.get("de"), dict) else {}) or {}
        if de_block_el.get("evidence_level"):
            de_evidence_level = str(de_block_el["evidence_level"])
        else:
            design = report.get("experimental_design") or {}
            if design.get("evidence_level"):
                de_evidence_level = str(design["evidence_level"])

    top_genes_list = _top_de_genes(de, n=top_genes)
    composition = _composition_for_pert(
        tables.get("cell_annotations"),
        tables.get("composition"),
        pert,
        control,
        cache=tables.get("composition_cache"),
    )

    consistency: dict[str, Any] = {
        "across_guides": NA,
        "across_samples": NA,
        "across_cell_states": composition if composition.get("source") != NA else NA,
    }
    if gc_row is not None:
        across_guides: dict[str, Any] = {
            "n_guides": n_guides,
            "guides_consistent": None if pd.isna(gc_row.get("guides_consistent")) else bool(gc_row["guides_consistent"]),
            "majority_direction": gc_row.get("majority_direction", NA),
            "median_target_log2fc": _safe_float(gc_row.get("median_target_log2fc")),
            "directions": gc_row.get("directions", NA),
        }
        if "guide_merge" in gc_row.index and pd.notna(gc_row.get("guide_merge")):
            across_guides["guide_merge"] = str(gc_row["guide_merge"])
        if "n_low_efficiency_guides" in gc_row.index:
            across_guides["n_low_efficiency_guides"] = _safe_int(gc_row.get("n_low_efficiency_guides"))
        if "n_on_target_pass" in gc_row.index:
            across_guides["n_on_target_pass"] = _safe_int(gc_row.get("n_on_target_pass"))
        if "potential_low_efficiency" in gc_row.index and pd.notna(gc_row.get("potential_low_efficiency")):
            across_guides["potential_low_efficiency"] = bool(gc_row["potential_low_efficiency"])
        if "interpretation" in gc_row.index and pd.notna(gc_row.get("interpretation")):
            across_guides["interpretation"] = str(gc_row["interpretation"])
        if "low_efficiency_guides" in gc_row.index and pd.notna(gc_row.get("low_efficiency_guides")):
            across_guides["low_efficiency_guides"] = str(gc_row["low_efficiency_guides"])
        consistency["across_guides"] = across_guides
    else:
        consistency["across_guides"] = "N/A — gene_guide_consistency.csv not present"

    if replicate_col is None:
        consistency["across_samples"] = "N/A — single-sample run (no biological replicates)"
    else:
        consistency["across_samples"] = {
            "replicate_col": replicate_col,
            "n_replicates": n_replicates_display,
            "note": "Cross-replicate consistency metric not computed; see DESeq2 if available",
        }

    # Per-guide interpretations + structured QC warnings for HTML cards / report.json
    guide_interpretations: list[dict[str, Any]] = []
    if gq is not None and not gq.empty and "gene_target" in gq.columns:
        gq_pert = gq[gq["gene_target"].astype(str) == pert]
        for _, grow in gq_pert.iterrows():
            if grow.get("is_control"):
                continue
            entry: dict[str, Any] = {
                "guide_id": str(grow.get("guide_id", "")),
                "interpretation": str(grow.get("interpretation", NA)),
                "assignment_confidence": str(grow.get("assignment_confidence", NA)),
                "target_effect_direction": str(grow.get("target_effect_direction", NA)),
                "n_cells": _safe_int(grow.get("n_cells")),
            }
            if "target_expr_log2fc_vs_control" in grow.index:
                entry["target_log2fc"] = _safe_float(grow.get("target_expr_log2fc_vs_control"))
            if "on_target_score" in grow.index:
                entry["on_target_score"] = _safe_float(grow.get("on_target_score"))
            if "on_target_pass" in grow.index and pd.notna(grow.get("on_target_pass")):
                entry["on_target_pass"] = bool(grow["on_target_pass"])
            if "low_efficiency_guide" in grow.index and pd.notna(grow.get("low_efficiency_guide")):
                entry["low_efficiency_guide"] = bool(grow["low_efficiency_guide"])
            guide_interpretations.append(entry)

    interp_by_guide = {str(i["guide_id"]): i["interpretation"] for i in guide_interpretations}
    gene_level_interp = None
    if gc_row is not None and "interpretation" in gc_row.index and pd.notna(gc_row.get("interpretation")):
        gene_level_interp = str(gc_row["interpretation"])
    if gene_level_interp is None and guide_interpretations:
        # Prefer low-efficiency / knockdown / inconsistency labels when present.
        for pref in (
            "potential_low_efficiency_guide",
            "target_knockdown_detected",
            "target_upregulated",
            "inconsistent_guides",
        ):
            hit = next((i for i in guide_interpretations if i["interpretation"] == pref), None)
            if hit:
                gene_level_interp = hit["interpretation"]
                break
        if gene_level_interp is None:
            gene_level_interp = guide_interpretations[0]["interpretation"]

    qc_warning_rows: list[dict[str, Any]] = []
    qw = tables.get("qc_warnings")
    if qw is not None and not qw.empty and "gene_target" in qw.columns:
        pert_warns = qw[qw["gene_target"].astype(str) == pert]
        for _, wrow in pert_warns.iterrows():
            entity = str(wrow.get("entity", ""))
            level = str(wrow.get("level", ""))
            linked = interp_by_guide.get(entity) if level == "guide" else gene_level_interp
            qc_warning_rows.append(
                {
                    "level": level,
                    "entity": entity,
                    "warning": str(wrow.get("warning", "qc")),
                    "detail": str(wrow.get("detail", "")),
                    "interpretation": linked or NA,
                }
            )

    weighted_summary: Any = NA
    w_src = tables.get("guide_weighted")
    w_row = _row_lookup(w_src, "gene_target", pert) if w_src is not None else None
    if w_row is None and gc_row is not None and "weighted_target_log2fc" in gc_row.index:
        w_row = gc_row
    if w_row is not None and "weighted_target_log2fc" in w_row.index and pd.notna(w_row.get("weighted_target_log2fc")):
        weighted_summary = {
            "guide_merge": str(w_row.get("guide_merge", gc_row.get("guide_merge") if gc_row is not None else NA)),
            "weighted_target_log2fc": _safe_float(w_row.get("weighted_target_log2fc")),
            "weighted_direction": str(w_row.get("weighted_direction", NA)),
            "weighted_direction_vote": str(w_row.get("weighted_direction_vote", NA)),
            "weight_sum": _safe_float(w_row.get("weight_sum")),
            "n_guides_weighted": _safe_int(w_row.get("n_guides_weighted")),
            "guide_weights": w_row.get("guide_weights", NA),
            "note": "Weighted summary is descriptive; inconsistency flags are retained",
        }

    limitations: list[str] = []
    warnings: list[str] = []
    mix_note = status_note(report.get("mixscape")) or (
        str(report["mixscape_skipped"]) if report.get("mixscape_skipped") else None
    )
    if mix_note:
        limitations.append(mix_note)
    de_block = report.get("de") if isinstance(report.get("de"), dict) else {}
    de_note = de_block.get("detail") or de_block.get("note") or report.get("de_note")
    if de_note:
        limitations.append(str(de_note))
    if replicate_col is None:
        limitations.append("Single-sample: no replicate-level confidence intervals from DESeq2")
    if de is None:
        limitations.append("No DE table for this perturbation (only top-N by E-distance are tested by default)")
    elif de_method == "wilcoxon_cell_level_exploratory":
        warnings.append("DE p-values are cell-level Wilcoxon (exploratory); not biological-replicate inference")
    elif de_method == "pydeseq2_pseudobulk_no_bio_reps":
        warnings.append(
            "DE uses exploratory PyDESeq2 without true biological replicates "
            "(sample_id or technical pseudo-replicates); not population inference"
        )
    if n_cells is not None and n_cells < min_cells:
        warnings.append(f"Low cell count ({n_cells} < {min_cells})")
    if edistance is None and de is None:
        warnings.append("No E-distance or DE evidence — human review recommended before claiming an effect")
    if et_row is not None and effect.get("significant_adj") is False and edistance is not None and edistance > 0:
        warnings.append("E-distance ranked but not significant after multiple-testing adjustment")
    if low_power:
        warnings.append(
            "low_power: cell count below etest_power_min_cells — do not treat E-test significance as conclusive"
        )
    if (
        et_row is not None
        and effect.get("significant_adj") is True
        and effect.get("significant_adj_reported") is False
    ):
        warnings.append("Raw significant_adj=True but significant_adj_reported=False due to low_power")
    if guide_status == "potential_low_efficiency_guides":
        warnings.append(
            "Potential low-efficiency guides (on-target proxy failed for ≥2 adequate guides) — "
            "do not conclude no phenotype; review gene_guide_consistency.csv"
        )
    if guide_status == "guides_inconsistent":
        warnings.append("Guide effects inconsistent across guides targeting this gene — review gene_guide_consistency.csv")
    for wrow in qc_warning_rows[:8]:
        warnings.append(f"{wrow['warning']}: {wrow['detail']}")
    limitations.append("Pathway enrichment not run (gene lists from DE tables only)")

    guide_qc_block: dict[str, Any] = {
        "status": qc_status,
        "warnings": qc_warning_rows[:8] if qc_warning_rows else [],
        "interpretations": guide_interpretations,
        "weighted_summary": weighted_summary,
    }

    figures = _resolve_pert_figures(tables.get("figures_dir"), pert)

    cluster = None
    pc = tables.get("pert_clusters")
    if pc is not None:
        crow = _row_lookup(pc, "gene_target", pert)
        if crow is not None and "pert_cluster" in crow.index:
            cluster = crow["pert_cluster"]

    return {
        "perturbation": pert,
        "counts": {
            "n_cells": _na(n_cells),
            "n_guides": _na(n_guides),
            "n_replicates": n_replicates_display,
            "qc_status": qc_status,
        },
        "guide_qc": guide_qc_block,
        "effect_vs_control": {
            **effect,
            "control": control,
            "de_min_pvalue_adj": _na(de_padj_min),
            "de_method": _na(de_method),
            "evidence_level": _na(de_evidence_level),
        },
        "consistency": consistency,
        "top_affected": {
            "genes": top_genes_list if top_genes_list else NA,
            "pathways": "N/A — pathway enrichment not run (optional; DE gene list above)",
            "evidence": de_method or ("edistance" if edistance is not None else NA),
            "evidence_level": _na(de_evidence_level),
            "perturbation_cluster": _na(cluster),
        },
        "limitations": limitations,
        "warnings": warnings,
        "figures": figures,
    }


def _resolve_pert_figures(figures_dir: Path | None, pert: str) -> dict[str, str | None]:
    """Relative paths under output_dir for card thumbnails (None if missing)."""
    out: dict[str, str | None] = {
        "volcano": None,
        "umap": None,
        "edistance": None,
        "guide_consistency": None,
        "target_validation": None,
        "effect_summary": None,
    }
    if figures_dir is None or not Path(figures_dir).exists():
        return out
    figures_dir = Path(figures_dir)
    safe = re.sub(r"[^\w.\-]+", "_", str(pert).strip()) or pert
    for candidate in (
        figures_dir / f"volcano_{pert}.png",
        figures_dir / f"volcano_{safe}.png",
        figures_dir / f"volcano_{pert.replace(' ', '_')}.png",
    ):
        if candidate.exists():
            out["volcano"] = f"figures/{candidate.name}"
            break
    for name in ("umap.png", "umap_cell_annotation.png"):
        if (figures_dir / name).exists():
            out["umap"] = f"figures/{name}"
            break
    if (figures_dir / "edistance.png").exists():
        out["edistance"] = "figures/edistance.png"
    if (figures_dir / "perturbation_effect_summary.png").exists():
        out["effect_summary"] = "figures/perturbation_effect_summary.png"
    for candidate in (
        figures_dir / f"guide_consistency_{pert}.png",
        figures_dir / f"guide_consistency_{safe}.png",
    ):
        if candidate.exists():
            out["guide_consistency"] = f"figures/{candidate.name}"
            break
    for candidate in (
        figures_dir / f"target_validation_{pert}.png",
        figures_dir / f"target_validation_{safe}.png",
    ):
        if candidate.exists():
            out["target_validation"] = f"figures/{candidate.name}"
            break
    return out


def build_analysis_verdict(
    summaries: list[dict[str, Any]],
    report: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Auto-generate a skim-friendly analysis verdict from guide QC + E-test + DE.

    Returns structured counts plus ``summary_sentence`` for report.json / HTML.
    """
    report = report or {}
    n = len(summaries)
    n_sig = 0
    n_low_power = 0
    n_inconsistent = 0
    for s in summaries:
        eff = s.get("effect_vs_control") or {}
        if eff.get("low_power") is True:
            n_low_power += 1
        if eff.get("significant_adj_reported") is True:
            n_sig += 1
        guides = (s.get("consistency") or {}).get("across_guides")
        if isinstance(guides, dict) and guides.get("guides_consistent") is False:
            n_inconsistent += 1
        elif (s.get("counts") or {}).get("qc_status") == "guides_inconsistent":
            n_inconsistent += 1

    caveats: list[str] = []
    config = report.get("config") or {}
    has_exploratory_de = any(
        (s.get("effect_vs_control") or {}).get("de_method")
        in {"wilcoxon_cell_level_exploratory", "pydeseq2_pseudobulk_no_bio_reps"}
        for s in summaries
    )
    de_block = report.get("de") if isinstance(report.get("de"), dict) else {}
    de_text = str(de_block.get("detail") or report.get("de_note") or "")
    if (
        config.get("replicate_col") is None
        or has_exploratory_de
        or "exploratory" in de_text.lower()
        or "Wilcoxon" in de_text
    ):
        caveats.append(
            "DE without true biological replicates is exploratory only "
            "(cells ≠ bio reps; Squair 2021)"
        )
    if n_low_power:
        caveats.append(
            f"{n_low_power} perturbation(s) marked low_power (n_cells below etest_power_min_cells)"
        )
    if n_inconsistent:
        caveats.append(
            f"{n_inconsistent} gene(s) with inconsistent multi-guide effects — do not pool without review"
        )

    sentence = (
        f"Among {n} perturbation(s), {n_sig} have E-test significance reported after power "
        f"filtering (non-low-power); {n_inconsistent} gene(s) show inconsistent multi-guide effects."
    )
    if caveats:
        # Avoid repeating the inconsistency clause already in the lead sentence.
        lead_caveats = [c for c in caveats if "inconsistent multi-guide" not in c]
        if lead_caveats:
            sentence += " Key caveats: " + "; ".join(lead_caveats[:3]) + "."

    return {
        "n_perturbations": n,
        "n_significant_reported": n_sig,
        "n_low_power": n_low_power,
        "n_guides_inconsistent": n_inconsistent,
        "caveats": caveats,
        "summary_sentence": sentence,
    }


def build_summary_table(perturbations: list[dict[str, Any]]) -> pd.DataFrame:
    """Publication-ready one-row-per-perturbation summary."""
    rows: list[dict[str, Any]] = []
    for s in perturbations:
        eff = s["effect_vs_control"]
        genes = s["top_affected"]["genes"]
        top = genes[0] if isinstance(genes, list) and genes else {}
        guides = s["consistency"].get("across_guides")
        guides_consistent: Any = NA
        if isinstance(guides, dict):
            gc = guides.get("guides_consistent")
            guides_consistent = NA if gc is None else gc
        warns = s.get("warnings") or []
        rows.append(
            {
                "perturbation": s["perturbation"],
                "n_cells": s["counts"]["n_cells"],
                "n_guides": s["counts"]["n_guides"],
                "n_replicates": s["counts"]["n_replicates"],
                "qc_status": s["counts"]["qc_status"],
                "edistance": eff.get("effect_size", NA),
                "edistance_ci_low": (
                    (eff.get("confidence_interval") or {}).get("low", NA)
                    if isinstance(eff.get("confidence_interval"), dict)
                    and (eff.get("confidence_interval") or {}).get("source") == "edistance_bootstrap_cells"
                    else NA
                ),
                "edistance_ci_high": (
                    (eff.get("confidence_interval") or {}).get("high", NA)
                    if isinstance(eff.get("confidence_interval"), dict)
                    and (eff.get("confidence_interval") or {}).get("source") == "edistance_bootstrap_cells"
                    else NA
                ),
                "pvalue_adj": eff.get("pvalue_adj", NA),
                "significant_adj": eff.get("significant_adj", NA),
                "significant_adj_reported": eff.get("significant_adj_reported", NA),
                "low_power": eff.get("low_power", NA),
                "de_method": eff.get("de_method", NA),
                "evidence_level": eff.get("evidence_level", NA),
                "de_min_pvalue_adj": eff.get("de_min_pvalue_adj", NA),
                "top_gene": top.get("gene", NA) if top else NA,
                "top_gene_lfc": top.get("logfoldchange", NA) if top else NA,
                "guides_consistent": guides_consistent,
                "n_warnings": len(warns),
                "warnings": "; ".join(warns) if warns else NA,
            }
        )
    return pd.DataFrame(rows)


def _dataframe_to_markdown(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No perturbations summarized._\n"
    cols = [str(c) for c in df.columns]
    header = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join("---" for _ in cols) + " |"
    lines = [header, sep]
    for _, row in df.iterrows():
        cells = []
        for c in cols:
            val = row[c]
            if val is None or (isinstance(val, float) and math.isnan(val)):
                text = NA
            else:
                text = str(val).replace("|", "\\|").replace("\n", " ")
            cells.append(text)
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines) + "\n"


def write_summary_tables(
    output_dir: Path | str,
    perturbations: list[dict[str, Any]],
) -> tuple[Path, Path]:
    """Write summary_table.csv and summary_table.md next to report.html."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    df = build_summary_table(perturbations)
    csv_path = output_dir / "summary_table.csv"
    md_path = output_dir / "summary_table.md"
    df.to_csv(csv_path, index=False)
    md_path.write_text(
        "# Perturbation summary table\n\n"
        "Publication-ready overview generated with the HTML report.\n\n"
        + _dataframe_to_markdown(df)
    )
    return csv_path, md_path


def _precompute_state_composition(cell_ann: pd.DataFrame | None, control: str) -> dict[str, dict[str, Any]]:
    """One-pass cell_state / phase counts per gene_target (avoids O(perts × cells) rescans)."""
    if cell_ann is None or cell_ann.empty or "gene_target" not in cell_ann.columns:
        return {}
    out: dict[str, dict[str, Any]] = {}
    ctrl_frac = None
    if "cell_state" in cell_ann.columns:
        ctrl = cell_ann[cell_ann["gene_target"].astype(str) == control]
        if not ctrl.empty:
            ctrl_frac = ctrl["cell_state"].astype(str).value_counts(normalize=True)
    for gene, sub in cell_ann.groupby(cell_ann["gene_target"].astype(str), sort=False):
        entry: dict[str, Any] = {"source": "cell_annotations", "n_cells": int(len(sub))}
        for col in ("cell_state", "phase", "cluster_annotation"):
            if col in sub.columns:
                entry[f"{col}_counts"] = sub[col].astype(str).value_counts().head(8).to_dict()
        if ctrl_frac is not None and "cell_state" in sub.columns:
            pert_frac = sub["cell_state"].astype(str).value_counts(normalize=True)
            delta = (pert_frac - ctrl_frac).dropna().sort_values(key=abs, ascending=False).head(5)
            entry["cell_state_delta_vs_control"] = {k: round(float(v), 4) for k, v in delta.items()}
        out[str(gene)] = entry
    return out


def build_perturbation_summaries(
    tables: dict[str, Any],
    *,
    control: str = CONTROL_DEFAULT,
    min_cells: int = 10,
    top_genes: int = 10,
) -> list[dict[str, Any]]:
    # Cache composition lookups so summarize_perturbation does not rescan cell_annotations.
    if "composition_cache" not in tables:
        tables = dict(tables)
        tables["composition_cache"] = _precompute_state_composition(tables.get("cell_annotations"), control)
        # Drop heavy frame after caching (optional composition table still used).
        tables["cell_annotations"] = None
    perts = list_perturbations(tables, control=control, min_cells=min_cells)
    summaries = [
        summarize_perturbation(p, tables, control=control, min_cells=min_cells, top_genes=top_genes) for p in perts
    ]

    def sort_key(s: dict[str, Any]) -> tuple:
        eff = s["effect_vs_control"]["effect_size"]
        eff_f = eff if isinstance(eff, (int, float)) else -1.0
        n = s["counts"]["n_cells"]
        n_f = n if isinstance(n, int) else -1
        return (-eff_f, -n_f, s["perturbation"])

    summaries.sort(key=sort_key)
    return summaries


def _fmt(value: Any) -> str:
    if value is None or value == NA:
        return NA
    if isinstance(value, float):
        if abs(value) >= 0.01 or value == 0:
            return f"{value:.4g}"
        return f"{value:.2e}"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, (list, dict)):
        return html.escape(json.dumps(value, indent=2, default=str))
    return html.escape(str(value))


def _fig_thumb(rel: str | None, label: str, *, max_width: int = 220) -> str:
    if not rel:
        return ""
    safe = html.escape(rel)
    return (
        f'<a class="fig-thumb" href="{safe}" title="{html.escape(label)}">'
        f'<img src="{safe}" alt="{html.escape(label)}" style="max-width:{max_width}px"/>'
        f"<span>{html.escape(label)}</span></a>"
    )


def _caveats_checklist_html(
    report: dict[str, Any] | None = None,
    *,
    caveats: tuple[str, ...] | list[str] | None = None,
) -> str:
    """Collapsible checklist from STATISTICAL_CAVEATS.md cores + run-specific notes."""
    items: list[str] = list(caveats if caveats is not None else CORE_STATISTICAL_CAVEATS)
    if report:
        ptype = (report.get("config") or {}).get("perturbation_type")
        kd_note = mixscape_kd_caveat(ptype if isinstance(ptype, str) else None)
        if kd_note and kd_note not in items:
            items.insert(0, kd_note)
        for c in report.get("statistical_caveats") or []:
            text = str(c)
            if text and text not in items:
                items.append(text)
    lis = "".join(
        f'<li><label><input type="checkbox"/> {html.escape(text)}</label></li>' for text in items
    )
    return f"""
<details class="caveats" id="statistical-caveats" open>
  <summary>Statistical caveats checklist (read before interpreting results)</summary>
  <ul class="checklist">{lis}</ul>
  <p class="meta">Full text: <code>docs/STATISTICAL_CAVEATS.md</code></p>
</details>
"""


def _key_viz_section(report: dict[str, Any]) -> str:
    """Top-of-page biologist figures (PNG links, not base64)."""
    plots = report.get("report_plots") if isinstance(report.get("report_plots"), dict) else {}
    effect = plots.get("effect_summary")
    thumbs = []
    if effect:
        thumbs.append(_fig_thumb(str(effect), "Perturbation effect summary", max_width=520))
    gc = plots.get("guide_consistency") or {}
    tv = plots.get("target_validation") or {}
    link_bits = []
    if isinstance(gc, dict) and gc:
        link_bits.append(
            "<li>Guide consistency: "
            + ", ".join(
                f'<a href="{html.escape(v)}">{html.escape(k)}</a>' for k, v in sorted(gc.items())
            )
            + "</li>"
        )
    if isinstance(tv, dict) and tv:
        link_bits.append(
            "<li>Target validation: "
            + ", ".join(
                f'<a href="{html.escape(v)}">{html.escape(k)}</a>' for k, v in sorted(tv.items())
            )
            + "</li>"
        )
    skips = plots.get("skipped_plots") or []
    skip_html = ""
    if skips:
        skip_html = (
            "<ul class=\"meta\">"
            + "".join(
                f"<li>Skipped <code>{html.escape(str(s.get('plot')))}</code>: "
                f"{html.escape(str(s.get('detail') or s.get('reason') or 'skipped'))}</li>"
                for s in skips[:12]
            )
            + "</ul>"
        )
    if not thumbs and not link_bits and not skip_html:
        return ""
    body = ""
    if thumbs:
        body += f'<div class="fig-row key-viz">{"".join(thumbs)}</div>'
    if link_bits:
        body += "<ul>" + "".join(link_bits) + "</ul>"
    body += skip_html
    return f"""
<section class="key-viz" id="key-visualizations">
  <h2>Key visualizations</h2>
  <p class="meta">Effect summary ranks multivariate shift; guide consistency tests multi-guide agreement;
  target validation checks KO/KD of the perturbed gene itself.</p>
  {body}
</section>
"""


def render_html_report(
    report: dict[str, Any],
    perturbations: list[dict[str, Any]],
    *,
    sample_id: str,
    figure_links: list[str] | None = None,
    caveats: tuple[str, ...] | list[str] | None = None,
) -> str:
    """HTML with relative paths into figures/ and per-perturbation Q&A sections."""
    figure_links = figure_links or []
    n = len(perturbations)
    rows = []
    for s in perturbations:
        eff = s["effect_vs_control"]
        lp = eff.get("low_power")
        lp_badge = ' <span class="badge">low_power</span>' if lp is True else ""
        rows.append(
            "<tr>"
            f"<td><a href='#pert-{html.escape(s['perturbation'])}'>{html.escape(s['perturbation'])}</a>"
            f"{lp_badge}</td>"
            f"<td>{_fmt(s['counts']['n_cells'])}</td>"
            f"<td>{_fmt(eff['effect_size'])}</td>"
            f"<td>{_fmt(eff['pvalue_adj'])}</td>"
            f"<td>{_fmt(eff.get('significant_adj_reported', eff.get('significant_adj')))}</td>"
            f"<td>{_fmt(s['counts']['qc_status'])}</td>"
            "</tr>"
        )

    sections = []
    for s in perturbations:
        pid = html.escape(s["perturbation"])
        genes = s["top_affected"]["genes"]
        if genes == NA:
            gene_html = f"<p>{NA}</p>"
        else:
            items = [
                f"<li><code>{html.escape(str(g.get('gene')))}</code> "
                f"LFC={_fmt(g.get('logfoldchange'))} padj={_fmt(g.get('pvalue_adj'))}</li>"
                for g in genes[:10]
            ]
            gene_html = "<ol>" + "".join(items) + "</ol>"
        figs = s.get("figures") or {}
        thumbs = "".join(
            [
                _fig_thumb(figs.get("volcano"), f"Volcano {s['perturbation']}", max_width=280),
                _fig_thumb(figs.get("guide_consistency"), "Guide consistency", max_width=220),
                _fig_thumb(figs.get("target_validation"), "Target validation", max_width=220),
                _fig_thumb(figs.get("umap"), "UMAP (viz only)", max_width=180),
                _fig_thumb(figs.get("edistance"), "E-distance bar", max_width=180),
            ]
        )
        fig_block = f'<div class="fig-row">{thumbs}</div>' if thumbs else f"<p>{NA}</p>"
        warn = s.get("warnings") or []
        lim = s.get("limitations") or []
        cons = s["consistency"]
        gqc = s.get("guide_qc") or {}
        interp = gqc.get("interpretations") or []
        if interp:
            interp_items = []
            for i in interp[:8]:
                extras = (
                    f"conf={html.escape(str(i.get('assignment_confidence')))}, "
                    f"dir={html.escape(str(i.get('target_effect_direction')))}, "
                    f"n={_fmt(i.get('n_cells'))}"
                )
                if i.get("on_target_score") is not None:
                    extras += f", on_target_score={_fmt(i.get('on_target_score'))}"
                if "on_target_pass" in i:
                    extras += f", on_target_pass={_fmt(i.get('on_target_pass'))}"
                if i.get("low_efficiency_guide"):
                    extras += ", <strong>low_efficiency_guide</strong>"
                interp_items.append(
                    f"<li><code>{html.escape(str(i.get('guide_id')))}</code>: "
                    f"<em>{html.escape(str(i.get('interpretation')))}</em> ({extras})</li>"
                )
            interp_html = "<ul>" + "".join(interp_items) + "</ul>"
        else:
            interp_html = f"<p>{NA}</p>"
        across = cons.get("across_guides") if isinstance(cons.get("across_guides"), dict) else {}
        if across.get("potential_low_efficiency"):
            low_eff_note = (
                f'<p class="warn"><strong>Potential low-efficiency guides</strong> '
                f"(n_low_efficiency={_fmt(across.get('n_low_efficiency_guides'))}/"
                f"{_fmt(across.get('n_guides'))}; "
                f"pass={_fmt(across.get('n_on_target_pass'))}; "
                f"interp={_fmt(across.get('interpretation'))}) — "
                f"do <em>not</em> conclude no phenotype for this gene.</p>"
            )
        else:
            low_eff_note = ""
        qc_warns = gqc.get("warnings") or []
        if qc_warns:
            warn_items = []
            for w in qc_warns[:8]:
                interp = w.get("interpretation")
                interp_s = (
                    f" — <em>{html.escape(str(interp))}</em>"
                    if interp and interp != NA
                    else ""
                )
                warn_items.append(
                    f'<li class="warn"><code>{html.escape(str(w.get("warning")))}</code>'
                    f" [{html.escape(str(w.get('level')))}/{html.escape(str(w.get('entity')))}]: "
                    f"{html.escape(str(w.get('detail')))}{interp_s}</li>"
                )
            qc_warn_html = "<ul>" + "".join(warn_items) + "</ul>"
        else:
            qc_warn_html = f"<p>{NA}</p>"
        wsum = gqc.get("weighted_summary")
        if isinstance(wsum, dict):
            weighted_html = (
                "<ul>"
                f"<li>Mode: {_fmt(wsum.get('guide_merge'))}</li>"
                f"<li>Weighted target log2FC: {_fmt(wsum.get('weighted_target_log2fc'))}</li>"
                f"<li>Weighted direction: {_fmt(wsum.get('weighted_direction'))} "
                f"(vote: {_fmt(wsum.get('weighted_direction_vote'))})</li>"
                f"<li>Weights: {_fmt(wsum.get('guide_weights'))}</li>"
                f"<li>Note: {_fmt(wsum.get('note'))}</li>"
                "</ul>"
            )
        else:
            weighted_html = f"<p>{NA} — set <code>guide_merge</code> to equal|umi|confidence|umi_confidence</p>"
        sections.append(
            f"""
<section class="pert" id="pert-{pid}">
  <h2>{pid}</h2>
  <h3>1. Counts &amp; QC</h3>
  <ul>
    <li>Cell count: {_fmt(s['counts']['n_cells'])}</li>
    <li>Guide count: {_fmt(s['counts']['n_guides'])}</li>
    <li>Replicate count: {_fmt(s['counts']['n_replicates'])}</li>
    <li>QC status: {_fmt(s['counts']['qc_status'])}</li>
  </ul>
  <h4>Guide interpretations</h4>
  {low_eff_note}
  {interp_html}
  <h4>QC warnings</h4>
  {qc_warn_html}
  <h4>Weighted guide summary</h4>
  {weighted_html}
  <h3>2. Effect vs {html.escape(str(s['effect_vs_control'].get('control', 'NT')))}</h3>
  <ul>
    <li>Metric: {_fmt(s['effect_vs_control']['metric'])}</li>
    <li>Effect size: {_fmt(s['effect_vs_control']['effect_size'])}</li>
    <li>Confidence interval: {_fmt(s['effect_vs_control']['confidence_interval'])}</li>
    <li>Adjusted p-value (E-test): {_fmt(s['effect_vs_control']['pvalue_adj'])}</li>
    <li>Significant (adj, raw): {_fmt(s['effect_vs_control']['significant_adj'])}</li>
    <li>Significant (adj, reported): {_fmt(s['effect_vs_control'].get('significant_adj_reported'))}</li>
    <li>low_power: {_fmt(s['effect_vs_control'].get('low_power'))}</li>
    <li>Secondary metrics: {_fmt(s['effect_vs_control'].get('secondary_metrics'))}</li>
    <li>DE method: {_fmt(s['effect_vs_control']['de_method'])}</li>
    <li>DE evidence_level: {_fmt(s['effect_vs_control'].get('evidence_level'))}</li>
    <li>Note: {_fmt(s['effect_vs_control'].get('note'))}</li>
  </ul>
  <h3>3. Consistency</h3>
  <ul>
    <li>Across guides: {_fmt(cons['across_guides'])}</li>
    <li>Across samples: {_fmt(cons['across_samples'])}</li>
    <li>Across cell types/states: {_fmt(cons['across_cell_states'])}</li>
  </ul>
  <h3>4. Top affected genes / pathways</h3>
  {gene_html}
  <p>Pathways: {_fmt(s['top_affected']['pathways'])}</p>
  <p>Evidence: {_fmt(s['top_affected']['evidence'])}; cluster: {_fmt(s['top_affected']['perturbation_cluster'])}</p>
  <h3>5. Limitations &amp; warnings</h3>
  <ul>{''.join(f'<li class="warn">{html.escape(w)}</li>' for w in warn) or f'<li>{NA}</li>'}</ul>
  <ul>{''.join(f'<li>{html.escape(x)}</li>' for x in lim)}</ul>
  <h3>Figures</h3>
  {fig_block}
</section>
"""
        )

    fig_list = "".join(f'<li><a href="{html.escape(f)}">{html.escape(f)}</a></li>' for f in figure_links)
    verdict = report.get("verdict") if isinstance(report.get("verdict"), dict) else {}
    summary_sentence = report.get("summary_sentence") or verdict.get("summary_sentence")
    if summary_sentence:
        verdict_html = f"""
<section class="verdict" id="analysis-verdict">
  <h2>Analysis verdict</h2>
  <p>{html.escape(str(summary_sentence))}</p>
</section>
"""
    else:
        verdict_html = ""
    global_notes = []
    for key in ANALYSIS_STATUS_KEYS:
        block = report.get(key)
        note = status_note(block)
        if note:
            label = html.escape(key)
            skipped = isinstance(block, dict) and block.get("skipped")
            cls = ' class="warn"' if skipped else ""
            global_notes.append(f"<li{cls}><code>{label}</code>: {html.escape(note)}</li>")
        elif isinstance(block, dict) and block.get("skipped") is False:
            global_notes.append(f"<li><code>{html.escape(key)}</code>: ran</li>")
        # Surface structured Mixscape skip estimate when present
        if key == "mixscape" and isinstance(block, dict):
            est = block.get("estimate")
            if isinstance(est, dict):
                bits = []
                if est.get("approx_work_units") is not None:
                    bits.append(f"approx_work_units={est['approx_work_units']}")
                if est.get("approx_relative_to_threshold") is not None:
                    bits.append(f"relative_to_threshold={est['approx_relative_to_threshold']}")
                if est.get("approx_memory_hint_gb") is not None:
                    bits.append(f"approx_memory_hint_gb={est['approx_memory_hint_gb']}")
                if bits:
                    global_notes.append(
                        f'<li class="warn"><code>mixscape.estimate</code>: {html.escape("; ".join(bits))}</li>'
                    )
            if block.get("subset"):
                targets = block.get("selected_targets") or []
                global_notes.append(
                    '<li class="warn"><code>mixscape.subset</code>: results cover only '
                    f"{html.escape(str(len(targets)))} selected targets + control — not the full library.</li>"
                )
    # Legacy flat notes (older report.json without structured status blocks)
    legacy_map = {
        "mixscape_skipped": "mixscape",
        "mixscape_error": "mixscape",
        "distance_error": "edistance",
        "de_note": "de",
    }
    for key, structured in legacy_map.items():
        if report.get(key) and not isinstance(report.get(structured), dict):
            global_notes.append(f"<li>{html.escape(str(report[key]))}</li>")

    caveats_html = _caveats_checklist_html(report, caveats=caveats)
    key_viz_html = _key_viz_section(report)

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<title>Perturb-seq report — {html.escape(sample_id)}</title>
<style>
  body {{ font-family: Georgia, 'Times New Roman', serif; margin: 2rem auto; max-width: 960px;
         line-height: 1.45; color: #1a1a1a; background: #faf9f7; padding: 0 1rem; }}
  h1,h2,h3,h4 {{ font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 0.9rem; margin: 1rem 0; }}
  th, td {{ border: 1px solid #ccc; padding: 0.35rem 0.5rem; text-align: left; }}
  th {{ background: #eee; }}
  .pert {{ border-top: 2px solid #333; margin-top: 2rem; padding-top: 0.5rem; }}
  .warn {{ color: #8a4b08; }}
  .badge {{ background: #f0e6d8; color: #8a4b08; font-size: 0.75rem; padding: 0.1rem 0.35rem;
            margin-left: 0.35rem; border: 1px solid #d9c7a8; }}
  code {{ background: #eee; padding: 0 0.2rem; }}
  .meta {{ color: #555; font-size: 0.95rem; }}
  img {{ border: 1px solid #ddd; }}
  .downloads a {{ margin-right: 1rem; }}
  .caveats {{ background: #f3efe6; border: 1px solid #d4cbb8; padding: 0.6rem 0.9rem; margin: 1rem 0; }}
  .caveats summary {{ cursor: pointer; font-weight: 600; font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; }}
  .verdict {{ background: #e8f0e9; border: 1px solid #b7c9b9; padding: 0.75rem 1rem; margin: 1rem 0; }}
  .verdict h2 {{ margin: 0 0 0.35rem; font-size: 1.1rem; }}
  .verdict p {{ margin: 0; }}
  .key-viz {{ background: #eef2f6; border: 1px solid #c5d0db; padding: 0.75rem 1rem; margin: 1rem 0; }}
  .key-viz h2 {{ margin: 0 0 0.35rem; font-size: 1.1rem; }}
  .checklist {{ list-style: none; padding-left: 0; }}
  .checklist li {{ margin: 0.35rem 0; }}
  .fig-row {{ display: flex; flex-wrap: wrap; gap: 0.75rem; align-items: flex-start; }}
  .fig-thumb {{ display: inline-flex; flex-direction: column; align-items: center; text-decoration: none;
                color: #333; font-size: 0.8rem; max-width: 280px; }}
  .fig-thumb img {{ display: block; }}
  .key-viz .fig-thumb {{ max-width: 560px; }}
</style>
</head>
<body>
<h1>Perturb-seq analysis report</h1>
<p class="meta">Sample: <strong>{html.escape(sample_id)}</strong> · {n} perturbations summarized</p>
<p class="meta">Machine-readable: <code>report.json</code> → <code>perturbations</code> / <code>verdict</code>. CSV / <code>.h5ad</code> remain downstream interfaces.</p>
{verdict_html}
{key_viz_html}
{caveats_html}
<p class="downloads">Publication-ready summary:
  <a href="summary_table.csv" download>summary_table.csv</a>
  <a href="summary_table.md" download>summary_table.md</a>
</p>
<h2>Run notes</h2>
<ul>{''.join(global_notes) or f'<li>{NA}</li>'}</ul>
<h2>Figures</h2>
<ul>{fig_list or f'<li>{NA}</li>'}</ul>
<h2>Perturbation index</h2>
<table>
<thead><tr><th>Perturbation</th><th>n_cells</th><th>E-distance</th><th>padj</th><th>sig_adj_reported</th><th>QC</th></tr></thead>
<tbody>
{''.join(rows)}
</tbody>
</table>
{''.join(sections)}
</body>
</html>
"""


def write_analysis_report(
    output_dir: Path | str,
    *,
    control: str | None = None,
    min_cells: int | None = None,
    top_genes: int = 10,
    html_detail_limit: int | None = None,
    adata: AnnData | None = None,
    h5ad: Path | str | None = None,
    report_plot_top_n: int | None = None,
) -> dict[str, Any]:
    """Enrich report.json with perturbations and write report.html."""
    output_dir = Path(output_dir)
    tables = load_result_tables(output_dir)
    report = dict(tables.get("report") or {})
    control = control or (report.get("config") or {}).get("control") or CONTROL_DEFAULT
    min_cells = min_cells or int((report.get("config") or {}).get("min_cells_per_pert") or 10)
    cfg = report.get("config") or {}
    top_n_plots = report_plot_top_n
    if top_n_plots is None:
        top_n_plots = int(cfg.get("report_plot_top_n") or 15)

    # Biologist figures before summaries so card thumbnails resolve.
    from perturbseq.report_plots import generate_report_plots

    plot_status = generate_report_plots(
        output_dir,
        tables=tables,
        adata=adata,
        h5ad=h5ad,
        control=control,
        top_n=top_n_plots,
        min_guide_cells=min_cells,
    )
    report["report_plots"] = plot_status

    summaries = build_perturbation_summaries(
        tables, control=control, min_cells=min_cells, top_genes=top_genes
    )
    report["perturbations"] = {s["perturbation"]: s for s in summaries}
    report["perturbation_summary"] = {
        "n_perturbations": len(summaries),
        "control": control,
        "n_with_edistance": sum(
            1 for s in summaries if isinstance(s["effect_vs_control"]["effect_size"], (int, float))
        ),
        "n_with_de": sum(1 for s in summaries if s["top_affected"]["genes"] != NA),
        "n_warnings": sum(len(s.get("warnings") or []) for s in summaries),
    }
    verdict = build_analysis_verdict(summaries, report)
    report["verdict"] = verdict
    report["summary_sentence"] = verdict["summary_sentence"]

    figures_dir = output_dir / "figures"
    figure_links: list[str] = []
    if figures_dir.exists():
        figure_links = sorted(f"figures/{p.name}" for p in figures_dir.glob("*.png"))

    html_perts = summaries if html_detail_limit is None else summaries[:html_detail_limit]
    if html_detail_limit is not None:
        de_names = set(tables.get("de_by_pert") or {})
        have = {s["perturbation"] for s in html_perts}
        for s in summaries:
            if s["perturbation"] in de_names and s["perturbation"] not in have:
                html_perts.append(s)
                have.add(s["perturbation"])

    sample_id = str(report.get("sample_id") or output_dir.name)
    html_text = render_html_report(report, html_perts, sample_id=sample_id, figure_links=figure_links)

    write_summary_tables(output_dir, summaries)
    write_report_json(report, output_dir)
    (output_dir / "report.html").write_text(html_text)
    write_provenance_manifest(output_dir, report)
    steps = list(report.get("steps") or [])
    if "analysis_report" not in steps:
        steps.append("analysis_report")
        report["steps"] = steps
        write_report_json(report, output_dir)
    return report


def write_html_summary(
    report: dict,
    output_dir: Path,
    filename: str = "report.html",
    *,
    adata: AnnData | None = None,
    h5ad: Path | str | None = None,
) -> Path:
    """Write HTML; prefer rich per-perturbation report when tables exist."""
    output_dir = Path(output_dir)
    tables_dir = output_dir / "tables"
    if tables_dir.exists() and any(tables_dir.glob("*.csv")):
        # Persist current report first so the builder can merge.
        write_report_json(report, output_dir)
        write_analysis_report(output_dir, adata=adata, h5ad=h5ad)
        path = output_dir / "report.html"
        if filename != "report.html":
            target = output_dir / filename
            target.write_text(path.read_text())
            return target
        return path
    # Minimal fallback (no analysis tables yet)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / filename
    sample = html.escape(str(report.get("sample_id", "")))
    steps = report.get("steps") or []
    stages = report.get("stages") or {}
    stage_items = stages.items() if isinstance(stages, dict) else [(s, "") for s in stages]
    rows = []
    for key in (
        "n_cells_loaded",
        "n_after_singlet",
        "n_gene_targets",
        "n_de_contrasts",
        "de_scope",
        "composition_audit",
    ):
        if key in report:
            rows.append(f"<tr><td>{html.escape(key)}</td><td>{html.escape(str(report[key]))}</td></tr>")
    for key in ANALYSIS_STATUS_KEYS:
        block = report.get(key)
        if isinstance(block, dict) and "skipped" in block:
            reason = block.get("reason")
            label = "skipped" if block.get("skipped") else "ran"
            text = f"{label}" + (f" ({reason})" if reason else "")
            if block.get("detail"):
                text += f" — {block['detail']}"
            rows.append(f"<tr><td>{html.escape(key)}</td><td>{html.escape(text)}</td></tr>")
        elif key == "cell_annotation" and key in report:
            rows.append(f"<tr><td>{html.escape(key)}</td><td>{html.escape(str(report[key]))}</td></tr>")
    prov = report.get("matrix_provenance") or {}
    if prov:
        rows.append(
            "<tr><td>matrix_provenance</td><td><code>"
            + html.escape(json.dumps(prov, default=str)[:800])
            + "</code></td></tr>"
        )
    notes = report.get("composition_notes") or []
    if notes:
        note_html = "; ".join(html.escape(str(n)) for n in notes[:8])
        rows.append(f"<tr><td>composition_notes</td><td>{note_html}</td></tr>")
    path.write_text(
        f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><title>Perturb-seq report — {sample}</title></head>
<body>
<h1>Perturb-seq tertiary report</h1>
<p>Sample: <strong>{sample}</strong></p>
<table><tbody>{''.join(rows)}</tbody></table>
<h2>Stages</h2>
<ul>{''.join(f'<li>{html.escape(str(k))}: {html.escape(str(v))}</li>' for k, v in stage_items)}</ul>
<h2>Steps</h2>
<ul>{''.join(f'<li>{html.escape(str(s))}</li>' for s in steps)}</ul>
<p>Machine-readable JSON: <code>report.json</code></p>
</body></html>
"""
    )
    return path


def finalize_outputs(
    adata: AnnData,
    report: dict,
    output_dir: Path,
    sample_id: str,
    *,
    inputs: dict[str, str] | None = None,
) -> dict:
    """Write H5AD + JSON/HTML reports + provenance manifest."""
    h5ad_path = output_dir / f"{sample_id}.tertiary.h5ad"
    write_h5ad(adata, h5ad_path)
    report["output_h5ad"] = str(h5ad_path)
    write_report_json(report, output_dir)
    # Pass in-memory adata so report plots skip a second disk read.
    write_html_summary(report, output_dir, adata=adata, h5ad=h5ad_path)
    # Reload enriched report (perturbations section) written by write_analysis_report.
    enriched_path = output_dir / "report.json"
    if enriched_path.exists():
        try:
            report = json.loads(enriched_path.read_text())
        except json.JSONDecodeError:
            pass
    write_provenance_manifest(output_dir, report, inputs=inputs)
    steps = list(report.get("steps") or [])
    if "report" not in steps:
        steps.append("report")
    report["steps"] = steps
    write_report_json(report, output_dir)
    return report


__all__ = [
    "ANALYSIS_STATUS_KEYS",
    "CORE_STATISTICAL_CAVEATS",
    "build_analysis_verdict",
    "build_perturbation_summaries",
    "build_summary_table",
    "collect_package_versions",
    "completed_status",
    "finalize_outputs",
    "init_report",
    "list_perturbations",
    "load_result_tables",
    "mixscape_kd_caveat",
    "render_html_report",
    "skipped_status",
    "status_note",
    "summarize_perturbation",
    "write_analysis_report",
    "write_html_summary",
    "write_provenance_manifest",
    "write_report_json",
    "write_summary_tables",
]
