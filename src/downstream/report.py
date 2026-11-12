"""
The contact sheet: one self-contained HTML page holding every figure.

Eyeballing a run means looking at forty PNGs. Opening them one by one in a file
browser loses the grouping and the context; this page keeps them in analysis
order, with the QC section first (thresholds and field flags come before any
figure that depends on them) and the run's headline numbers at the top.

Images are referenced by relative path rather than embedded, so the page stays
small and reflects re-rendered figures on reload: it lives alongside the
``figures/`` tree it points at.
"""

from __future__ import annotations

import html
import logging
from datetime import datetime
from pathlib import Path

import pandas as pd

import provenance

logger = logging.getLogger(__name__)

# Section order and the human-readable blurb shown under each heading.
SECTIONS = [
    ("qc", "Quality control",
     "Segmentation yield, staining consistency and the evidence behind every "
     "positivity threshold. Read this before trusting anything below."),
    ("signaling", "Signaling rings (QC)",
     "Where the nuclear:cytoplasmic ratio was measured for signaling markers "
     "(β-catenin, SMAD2/3, ...): nucleus cores (orange) and the cytoplasmic ring "
     "inside cell bodies (sky blue). Rings should hug the nuclei."),
    ("distributions", "Distributions",
     "Per-marker intensity by condition. Violins show cells; the black points are "
     "field medians. Signaling markers are shown relative to their reference wells."),
    ("composition", "Composition",
     "Positive fractions with confidence intervals across replicates (see the note "
     "at the top for what the replicate is), the graded neg/low/mid/high breakdown, "
     "and marker co-expression."),
    ("joint", "Joint marker views",
     "Two-marker densities and per-cell scatters with the positivity gates "
     "(dashed) and tertile edges (dotted) drawn on. With three or more markers "
     "the scatters are colored by each remaining marker."),
    ("plate", "Plate layout & uniformity",
     "Well statistics in plate geometry (edge effects, pipetting slips) and the "
     "field-to-field spread within each well."),
    ("trends", "Trends & factor effects",
     "Readouts across the configured time/dose axis, and the marginal effect of "
     "each experimental factor."),
    ("effects", "Effects, dose trends and power",
     "From every per-field number (positive fraction, median). Trends: replicate "
     "values against a numeric condition column with the fitted slope, its bootstrap "
     "95 % CI and a permutation p. Differences from the reference condition ± 95 % "
     "CI. Power: the smallest difference two conditions could be shown to differ by, "
     "against the number of independent experiments."),
]

_CSS = """
:root { color-scheme: light dark; }
body { font-family: -apple-system, "Segoe UI", Roboto, sans-serif; margin: 0;
       background: #fbfbfc; color: #1a1a1a; }
header { background: #0072B2; color: white; padding: 1.4rem 2rem; }
header h1 { margin: 0 0 .3rem; font-size: 1.35rem; }
header p { margin: 0; opacity: .9; font-size: .85rem; }
header p.version { margin-top: .3rem; opacity: .7; font-size: .75rem; }
header a { color: white; }
main { padding: 1.2rem 2rem 3rem; max-width: 1500px; margin: 0 auto; }
section { margin-bottom: 2.4rem; }
h2 { font-size: 1.05rem; border-bottom: 2px solid #0072B2; padding-bottom: .3rem;
     margin-bottom: .3rem; }
.blurb { font-size: .82rem; color: #555; margin: .3rem 0 1rem; }
.grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(340px, 1fr));
        gap: 1rem; }
figure { margin: 0; background: white; border: 1px solid #e2e2e6; border-radius: 6px;
         padding: .6rem; }
figure img { width: 100%; height: auto; display: block; cursor: zoom-in; }
figcaption { font-size: .72rem; color: #666; margin-top: .4rem; word-break: break-all; }
table { border-collapse: collapse; font-size: .78rem; background: white; }
th, td { border: 1px solid #e2e2e6; padding: .3rem .55rem; text-align: left; }
th { background: #f0f2f5; }
.stats { display: flex; gap: 1.6rem; flex-wrap: wrap; margin: 1rem 0 1.6rem; }
.stat { background: white; border: 1px solid #e2e2e6; border-radius: 6px;
        padding: .7rem 1.1rem; min-width: 110px; }
.stat .value { font-size: 1.35rem; font-weight: 600; color: #0072B2; }
.stat .label { font-size: .72rem; color: #666; text-transform: uppercase;
               letter-spacing: .04em; }
.warn { background: #fff4e5; border-left: 4px solid #E69F00; padding: .7rem 1rem;
        font-size: .82rem; border-radius: 4px; margin-bottom: 1rem; }
a { color: #0072B2; }
@media (prefers-color-scheme: dark) {
  body { background: #16181c; color: #e8e8ea; }
  figure, table, .stat { background: #1e2126; border-color: #2c3038; }
  th { background: #24272e; }
  .blurb, figcaption, .stat .label { color: #a0a4ab; }
  .warn { background: #2e2517; }
}
"""


def _figure_cards(directory: Path, root: Path) -> str:
    if not directory.exists():
        return ""
    cards = []
    for png in sorted(directory.glob("*.png")):
        rel = png.relative_to(root).as_posix()
        cards.append(
            f'<figure><a href="{html.escape(rel)}" target="_blank">'
            f'<img src="{html.escape(rel)}" loading="lazy" alt="{html.escape(png.stem)}"></a>'
            f'<figcaption>{html.escape(png.name)}</figcaption></figure>'
        )
    return f'<div class="grid">{"".join(cards)}</div>' if cards else ""


def _stat_card(value, label) -> str:
    return f'<div class="stat"><div class="value">{value}</div><div class="label">{label}</div></div>'


def _table_html(df: pd.DataFrame, max_rows: int = 40) -> str:
    if df is None or df.empty:
        return "<p><em>none</em></p>"
    shown = df.head(max_rows)
    head = "".join(f"<th>{html.escape(str(c))}</th>" for c in shown.columns)
    body = "".join(
        "<tr>" + "".join(f"<td>{html.escape(_fmt(v))}</td>" for v in row) + "</tr>"
        for row in shown.itertuples(index=False)
    )
    more = (f"<p><em>… {len(df) - max_rows} more row(s) in the CSV</em></p>"
            if len(df) > max_rows else "")
    return f"<table><tr>{head}</tr>{body}</table>{more}"


def _fmt(value) -> str:
    if isinstance(value, float):
        return "-" if pd.isna(value) else (f"{value:.4g}")
    return str(value)


def build_report(ctx) -> Path:
    """Render ``report.html`` at the analysis root."""
    root = ctx.out
    parts: list[str] = []

    # Headline numbers.
    stats_html = ""
    try:
        cells = ctx.data.cells
        stats_html = "".join([
            _stat_card(f"{len(cells):,}", "cells"),
            _stat_card(f"{cells['field_id'].nunique():,}", "fields"),
            _stat_card(f"{cells['well'].nunique():,}" if "well" in cells else "-", "wells"),
            _stat_card(f"{len(ctx.markers)}", "markers"),
            _stat_card(f"{cells['dataset'].nunique()}", "datasets"),
        ])
    except Exception:
        logger.debug("report: cell tables unavailable for the summary bar", exc_info=True)

    # Anything the run refused to gate is stated up front, not buried.
    warnings_html = ""
    unit = "field"
    try:
        import json
        with open(root / "tables" / "replicate_unit.json") as f:
            rec = json.load(f)
        unit = rec.get("unit", unit)
        warnings_html += f'<div class="warn">{html.escape(rec["note"])}</div>'
    except Exception:
        pass
    try:
        sat = pd.read_csv(root / "tables" / "qc_saturated_nuclei.csv")
        if len(sat):
            warnings_html += (f'<div class="warn">{sat["cell_uid"].nunique():,} nuclei sat at '
                              'the detector ceiling and were left out (bright debris; '
                              '<code>tables/qc_saturated_nuclei.csv</code>).</div>')
    except Exception:
        pass
    try:
        thr = ctx.thresholds
        graded = thr[thr["verdict"] != "ok"]
        if len(graded):
            names = ", ".join(sorted(set(graded["marker"].astype(str))))
            warnings_html += (
                f'<div class="warn"><strong>{len(graded)} marker/scope combination(s) '
                f'have no positivity gate</strong> ({html.escape(names)}): no positive '
                f'fraction is reported for them. See the threshold evidence panel and '
                f'<code>tables/thresholds.csv</code>.</div>'
            )
    except Exception:
        thr = None

    # Known sections in analysis order, then any group a custom step invented.
    known = {name for name, _, _ in SECTIONS}
    figures_root = root / "figures"
    extra = sorted(d.name for d in figures_root.iterdir()
                   if d.is_dir() and d.name not in known) if figures_root.exists() else []
    sections = list(SECTIONS) + [(name, name.replace("_", " ").title(),
                                  "Produced by a custom analysis step.") for name in extra]

    for folder, title, blurb in sections:
        cards = _figure_cards(figures_root / folder, root)
        if not cards:
            continue
        parts.append(f'<section id="{folder}"><h2>{html.escape(title)}</h2>'
                     f'<p class="blurb">{html.escape(blurb)}</p>{cards}</section>')

    # Key tables inline.
    tables_html = ""
    for name, caption in (("thresholds.csv", "Thresholds and the evidence for each"),
                          ("control_assignment.csv",
                           "Which wells supplied each marker's negative control"),
                          ("conditions.csv", "Condition summary (statistics across replicates)"),
                          ("comparisons_positive_fraction.csv",
                           "Each condition against its marker's control: positive fraction")):
        path = root / "tables" / name
        if path.exists():
            tables_html += (f"<h3>{html.escape(caption)} "
                            f'<a href="tables/{name}">({name})</a></h3>'
                            + _table_html(pd.read_csv(path)))
    if tables_html:
        parts.append(f'<section id="tables"><h2>Tables</h2>{tables_html}</section>')

    nav = " · ".join(f'<a href="#{f}">{html.escape(t)}</a>'
                     for f, t, _ in SECTIONS if (root / "figures" / f).exists())
    unit_text = {"field": "fields", "well": "wells", "experiment": "experiments"}.get(unit, unit)

    page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Downstream analysis: {html.escape(root.name)}</title>
<style>{_CSS}</style></head>
<body>
<header>
  <h1>Downstream analysis: {html.escape(root.parent.name)}</h1>
  <p>Generated {datetime.now():%Y-%m-%d %H:%M} · replicates: {html.escape(unit_text)} ·
     {nav}</p>
  <p class="version">{html.escape(provenance.summary())}</p>
</header>
<main>
<div class="stats">{stats_html}</div>
{warnings_html}
{"".join(parts)}
</main></body></html>"""

    out = root / "report.html"
    out.write_text(page, encoding="utf-8")
    logger.info("report → %s", out)
    return out


__all__ = ["build_report"]
