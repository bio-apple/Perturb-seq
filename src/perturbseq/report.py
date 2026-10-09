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
    guide_qc = _read_csv(tables / "guide_qc.csv")
    guide_consistency = _read_csv(tables / "gene_guide_consistency.csv")
    qc_warnings = _read_csv(tables / "qc_warnings.csv")
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
        "guide_qc": guide_qc,
        "guide_consistency": guide_consistency,
        "qc_warnings": qc_warnings,
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
    if gc_row is not None and "guides_consistent" in gc_row.index:
        cons_val = gc_row["guides_consistent"]
        if pd.isna(cons_val):
            guide_status = "consistency_unknown"
        elif bool(cons_val):
            guide_status = "guides_consistent"
        else:
            guide_status = "guides_inconsistent"

    gq = tables.get("guide_qc")
    if n_guides is None and gq is not None and "gene_target" in gq.columns:
        n_guides = int((gq["gene_target"].astype(str) == pert).sum())

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
    effect: dict[str, Any] = {
        "metric": "edistance" if edistance is not None else NA,
        "effect_size": _na(edistance),
        "confidence_interval": NA,
        "pvalue": NA,
        "pvalue_adj": NA,
        "significant_adj": NA,
        "note": None,
    }
    if et_row is not None:
        effect["pvalue"] = _na(_safe_float(et_row.get("pvalue")))
        effect["pvalue_adj"] = _na(_safe_float(et_row.get("pvalue_adj")))
        if "significant_adj" in et_row.index:
            effect["significant_adj"] = bool(et_row["significant_adj"])
        elif isinstance(effect["pvalue_adj"], float):
            effect["significant_adj"] = effect["pvalue_adj"] < 0.05
        if edistance is None and "distance" in et_row.index:
            effect["effect_size"] = _na(_safe_float(et_row["distance"]))
            effect["metric"] = "edistance"
    elif edistance is not None:
        effect["note"] = "E-distance available; E-test not run or not tested for this target"

    de = (tables.get("de_by_pert") or {}).get(pert)
    de_padj_min = None
    de_method = None
    if de is not None and not de.empty:
        de_method = str(de["method"].iloc[0]) if "method" in de.columns else NA
        top1 = _top_de_genes(de, n=1)
        if top1 and "ci_low" in top1[0]:
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
            replicate_col is None or de_method == "wilcoxon_cell_level_exploratory"
        ):
            note = "No replicate CI (single-sample / exploratory Wilcoxon)"
            effect["note"] = f"{effect['note']}; {note}" if effect.get("note") else note

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
        consistency["across_guides"] = {
            "n_guides": n_guides,
            "guides_consistent": None if pd.isna(gc_row.get("guides_consistent")) else bool(gc_row["guides_consistent"]),
            "majority_direction": gc_row.get("majority_direction", NA),
            "median_target_log2fc": _safe_float(gc_row.get("median_target_log2fc")),
            "directions": gc_row.get("directions", NA),
        }
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

    limitations: list[str] = []
    warnings: list[str] = []
    if report.get("mixscape_skipped"):
        limitations.append(str(report["mixscape_skipped"]))
    if report.get("de_note"):
        limitations.append(str(report["de_note"]))
    if replicate_col is None:
        limitations.append("Single-sample: no replicate-level confidence intervals from DESeq2")
    if de is None:
        limitations.append("No DE table for this perturbation (only top-N by E-distance are tested by default)")
    elif de_method == "wilcoxon_cell_level_exploratory":
        warnings.append("DE p-values are cell-level Wilcoxon (exploratory); not biological-replicate inference")
    if n_cells is not None and n_cells < min_cells:
        warnings.append(f"Low cell count ({n_cells} < {min_cells})")
    if edistance is None and de is None:
        warnings.append("No E-distance or DE evidence — human review recommended before claiming an effect")
    if et_row is not None and effect.get("significant_adj") is False and edistance is not None and edistance > 0:
        warnings.append("E-distance ranked but not significant after multiple-testing adjustment")
    if guide_status == "guides_inconsistent":
        warnings.append("Guide effects inconsistent across guides targeting this gene — review gene_guide_consistency.csv")
    qw = tables.get("qc_warnings")
    if qw is not None and not qw.empty and "gene_target" in qw.columns:
        pert_warns = qw[qw["gene_target"].astype(str) == pert]
        for _, wrow in pert_warns.head(5).iterrows():
            warnings.append(f"{wrow.get('warning', 'qc')}: {wrow.get('detail', wrow.get('entity', ''))}")
    limitations.append("Pathway enrichment not run (gene lists from DE tables only)")

    figures = tables.get("figures_dir")
    volcano_rel = None
    if figures is not None:
        for candidate in (figures / f"volcano_{pert}.png", figures / f"volcano_{pert.replace(' ', '_')}.png"):
            if candidate.exists():
                volcano_rel = f"figures/{candidate.name}"
                break

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
        "effect_vs_control": {
            **effect,
            "control": control,
            "de_min_pvalue_adj": _na(de_padj_min),
            "de_method": _na(de_method),
        },
        "consistency": consistency,
        "top_affected": {
            "genes": top_genes_list if top_genes_list else NA,
            "pathways": "N/A — pathway enrichment not run (optional; DE gene list above)",
            "evidence": de_method or ("edistance" if edistance is not None else NA),
            "perturbation_cluster": _na(cluster),
        },
        "limitations": limitations,
        "warnings": warnings,
        "figures": {"volcano": volcano_rel},
    }


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


def render_html_report(
    report: dict[str, Any],
    perturbations: list[dict[str, Any]],
    *,
    sample_id: str,
    figure_links: list[str] | None = None,
) -> str:
    """HTML with relative paths into figures/ and per-perturbation Q&A sections."""
    figure_links = figure_links or []
    n = len(perturbations)
    rows = []
    for s in perturbations:
        eff = s["effect_vs_control"]
        rows.append(
            "<tr>"
            f"<td><a href='#pert-{html.escape(s['perturbation'])}'>{html.escape(s['perturbation'])}</a></td>"
            f"<td>{_fmt(s['counts']['n_cells'])}</td>"
            f"<td>{_fmt(eff['effect_size'])}</td>"
            f"<td>{_fmt(eff['pvalue_adj'])}</td>"
            f"<td>{_fmt(eff['significant_adj'])}</td>"
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
        volcano = s.get("figures", {}).get("volcano")
        vol_html = (
            f'<p><img src="{html.escape(volcano)}" alt="volcano {pid}" style="max-width:420px"/></p>'
            if volcano
            else f"<p>{NA}</p>"
        )
        warn = s.get("warnings") or []
        lim = s.get("limitations") or []
        cons = s["consistency"]
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
  <h3>2. Effect vs {html.escape(str(s['effect_vs_control'].get('control', 'NT')))}</h3>
  <ul>
    <li>Metric: {_fmt(s['effect_vs_control']['metric'])}</li>
    <li>Effect size: {_fmt(s['effect_vs_control']['effect_size'])}</li>
    <li>Confidence interval: {_fmt(s['effect_vs_control']['confidence_interval'])}</li>
    <li>Adjusted p-value (E-test): {_fmt(s['effect_vs_control']['pvalue_adj'])}</li>
    <li>Significant (adj): {_fmt(s['effect_vs_control']['significant_adj'])}</li>
    <li>DE method: {_fmt(s['effect_vs_control']['de_method'])}</li>
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
  <h3>Volcano</h3>
  {vol_html}
</section>
"""
        )

    fig_list = "".join(f'<li><a href="{html.escape(f)}">{html.escape(f)}</a></li>' for f in figure_links)
    global_notes = []
    for key in ("mixscape_skipped", "de_note", "distance_error", "mixscape_error"):
        if report.get(key):
            global_notes.append(f"<li>{html.escape(str(report[key]))}</li>")

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<title>Perturb-seq report — {html.escape(sample_id)}</title>
<style>
  body {{ font-family: Georgia, 'Times New Roman', serif; margin: 2rem auto; max-width: 960px;
         line-height: 1.45; color: #1a1a1a; background: #faf9f7; padding: 0 1rem; }}
  h1,h2,h3 {{ font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 0.9rem; margin: 1rem 0; }}
  th, td {{ border: 1px solid #ccc; padding: 0.35rem 0.5rem; text-align: left; }}
  th {{ background: #eee; }}
  .pert {{ border-top: 2px solid #333; margin-top: 2rem; padding-top: 0.5rem; }}
  .warn {{ color: #8a4b08; }}
  code {{ background: #eee; padding: 0 0.2rem; }}
  .meta {{ color: #555; font-size: 0.95rem; }}
  img {{ border: 1px solid #ddd; }}
</style>
</head>
<body>
<h1>Perturb-seq analysis report</h1>
<p class="meta">Sample: <strong>{html.escape(sample_id)}</strong> · {n} perturbations summarized</p>
<p class="meta">Machine-readable: <code>report.json</code> → <code>perturbations</code>. CSV / <code>.h5ad</code> remain downstream interfaces.</p>
<h2>Run notes</h2>
<ul>{''.join(global_notes) or f'<li>{NA}</li>'}</ul>
<h2>Figures</h2>
<ul>{fig_list or f'<li>{NA}</li>'}</ul>
<h2>Perturbation index</h2>
<table>
<thead><tr><th>Perturbation</th><th>n_cells</th><th>E-distance</th><th>padj</th><th>sig_adj</th><th>QC</th></tr></thead>
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
) -> dict[str, Any]:
    """Enrich report.json with perturbations and write report.html."""
    output_dir = Path(output_dir)
    tables = load_result_tables(output_dir)
    report = dict(tables.get("report") or {})
    control = control or (report.get("config") or {}).get("control") or CONTROL_DEFAULT
    min_cells = min_cells or int((report.get("config") or {}).get("min_cells_per_pert") or 10)

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

    write_report_json(report, output_dir)
    (output_dir / "report.html").write_text(html_text)
    write_provenance_manifest(output_dir, report)
    steps = list(report.get("steps") or [])
    if "analysis_report" not in steps:
        steps.append("analysis_report")
        report["steps"] = steps
        write_report_json(report, output_dir)
    return report


def write_html_summary(report: dict, output_dir: Path, filename: str = "report.html") -> Path:
    """Write HTML; prefer rich per-perturbation report when tables exist."""
    output_dir = Path(output_dir)
    tables_dir = output_dir / "tables"
    if tables_dir.exists() and any(tables_dir.glob("*.csv")):
        # Persist current report first so the builder can merge.
        write_report_json(report, output_dir)
        write_analysis_report(output_dir)
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
        "mixscape_status",
        "de_scope",
        "composition_audit",
    ):
        if key in report:
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
    write_html_summary(report, output_dir)
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
    "build_perturbation_summaries",
    "collect_package_versions",
    "finalize_outputs",
    "init_report",
    "list_perturbations",
    "load_result_tables",
    "render_html_report",
    "summarize_perturbation",
    "write_analysis_report",
    "write_html_summary",
    "write_provenance_manifest",
    "write_report_json",
]
