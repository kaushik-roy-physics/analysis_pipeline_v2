"""
Example: an experiment-specific analysis that plugs into the generic pipeline.

Copy this next to your ``config.yaml`` (e.g. as ``downstream_analysis/custom.py``)
and point the config at it::

    downstream:
      extra_modules: [downstream_analysis.custom]
      steps: [assemble, qc, thresholds, aggregate, gating, distributions,
              composition, fate_calls, report]          # your step, in order

Nothing in ``src/downstream`` needs to change. Your step gets the same context,
the same checkpointing (``--force-step fate_calls``), and the same
"one failure doesn't kill the run" handling as the built-ins.

The step below is a real, if simple, example: a three-way fate call built from
the band assignments the generic ``gating`` step already computed. Swap the
marker names and the rules for your experiment's biology: the plumbing stays.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from downstream import register
from downstream.plots import save_figure, shorten
from downstream import style

logger = logging.getLogger(__name__)

# --- your experiment's biology, in one place -------------------------------
PROGENITOR_MARKER = "ISL1"
DIFFERENTIATED_MARKER = "NKX2.5"

FATES = {
    "Progenitor": "high progenitor marker, low differentiation marker",
    "Transitioning": "both present",
    "Differentiated": "high differentiation marker, progenitor marker low/absent",
    "Other": "neither",
}
FATE_COLORS = {
    "Progenitor": style.OKABE_ITO["green"],
    "Transitioning": style.OKABE_ITO["yellow"],
    "Differentiated": style.OKABE_ITO["vermillion"],
    "Other": style.GREY,
}


@register("fate_calls", after=["gating"], config_keys=["gating", "conditions"],
          description="Example: per-cell fate classification from the graded bands")
def step_fate_calls(ctx) -> None:
    ctx.reset_figures("fate")      # each step owns one figure group
    banded = ctx.banded
    markers = set(banded["marker"].unique())
    missing = {PROGENITOR_MARKER, DIFFERENTIATED_MARKER} - markers
    if missing:
        logger.warning("fate_calls: %s absent from this experiment; skipping.", missing)
        return

    # Bands are per (cell × marker); pivot to one row per cell.
    wide = banded.pivot_table(index="cell_uid", columns="marker", values="band",
                             observed=True, aggfunc="first")
    wide = wide[[PROGENITOR_MARKER, DIFFERENTIATED_MARKER]].dropna(how="any")
    if wide.empty:
        logger.warning("fate_calls: no cells measured for both markers; skipping.")
        return

    prog = wide[PROGENITOR_MARKER].astype(str)
    diff = wide[DIFFERENTIATED_MARKER].astype(str)

    fate = np.select(
        [
            (prog == "high") & (diff == "neg"),
            (prog != "neg") & (diff != "neg"),
            (prog == "neg") & (diff.isin(["mid", "high"])),
        ],
        ["Progenitor", "Transitioning", "Differentiated"],
        default="Other",
    )

    cells = pd.DataFrame({"cell_uid": wide.index, "fate": fate})
    cells["field_id"] = ctx.data.cell_column("field_id", cells["cell_uid"].to_numpy())
    for key in ctx.condition_keys:
        cells[key] = ctx.data.cell_column(key, cells["cell_uid"].to_numpy())
    from downstream.aggregate import condition_label
    cells["condition"] = condition_label(cells, ctx.condition_keys)

    # Fractions per field first, then averaged: fields are the replicates.
    per_field = (cells.groupby(["condition", "field_id", "fate"], observed=True)
                 .size().rename("n").reset_index())
    totals = per_field.groupby(["condition", "field_id"], observed=True)["n"].transform("sum")
    per_field["fraction"] = per_field["n"] / totals

    table = per_field.groupby(["condition", "fate"], observed=True).agg(
        n_cells=("n", "sum"),
        n_fields=("field_id", "nunique"),
        fraction=("fraction", "mean"),
        fraction_sd=("fraction", "std"),
    ).reset_index()
    ctx.write_table(table, "fate_proportions.csv")

    _plot(ctx, table)


def _plot(ctx, table: pd.DataFrame) -> None:
    import matplotlib.pyplot as plt
    from downstream.plots import ordered_conditions

    order = ordered_conditions(table["condition"].unique(), ctx.cfg, ctx.control_conditions)
    pivot = (table.pivot_table(index="condition", columns="fate", values="fraction",
                               observed=True, aggfunc="first")
             .reindex(order).fillna(0.0))

    fig, ax = plt.subplots(figsize=(max(6.0, 0.6 * len(pivot) + 2.4), 4.4))
    bottom = np.zeros(len(pivot))
    x = np.arange(len(pivot))
    for fate in FATES:
        if fate not in pivot.columns:
            continue
        vals = 100 * pivot[fate].to_numpy(float)
        ax.bar(x, vals, bottom=bottom, color=FATE_COLORS[fate], label=fate,
               edgecolor="white", linewidth=0.4)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels([shorten(c) for c in pivot.index], rotation=40, ha="right")
    ax.set_ylabel("Cells (%)")
    ax.set_ylim(0, 100)
    ax.set_title(f"Fate composition ({PROGENITOR_MARKER} × {DIFFERENTIATED_MARKER})")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
    save_figure(fig, ctx.figure("fate", "fate_composition.png"))


__all__ = ["step_fate_calls"]
