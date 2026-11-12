"""
Composition figures: how much of a population is where.

* :func:`positive_fraction_bars`: % positive per condition, error bars from the
  **fields** (never the cells), with the field values themselves overplotted.
* :func:`band_composition_bars`: the graded version: neg / low / mid / high
  stacked per condition, which shows a shift *within* the positive population
  that a binary fraction cannot.
* :func:`coexpression_grid`: the 4×4 band matrix per condition, i.e. what
  fraction of cells co-express A and B and in which part of each distribution.
* :func:`quadrant_bars`: the same data collapsed to the familiar
  double-positive / single / double-negative quadrants.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .. import style
from ..gating import BAND_ORDER, HIGH, LOW, MID, NEG
from .base import (plural, annotate_n, condition_palette, grid_shape, no_data,
                   ordered_conditions, save_figure, shorten, shorten_labels)

logger = logging.getLogger(__name__)

# Ordered bands share a sequential ramp so "more positive" reads as "darker":
# in color and in grayscale. Sampled from the sequential map rather than picked
# by hand, because a hand-picked qualitative triple is not monotone in luminance.
_RAMP = style.band_ramp(3)
BAND_COLORS = {
    NEG: style.NEUTRAL_FILL,
    LOW: _RAMP[0],
    MID: _RAMP[1],
    HIGH: _RAMP[2],
}


def positive_fraction_bars(conditions: pd.DataFrame, fields: pd.DataFrame, marker: str,
                           cfg: dict, out_path, controls=(), facet: str = ""):
    """% positive per condition: field-level CI, individual fields overplotted."""
    sub = conditions[conditions["marker"] == marker]
    if sub.empty or not np.isfinite(sub["positive_fraction"]).any():
        return None
    order = ordered_conditions(sub["condition"].unique(), cfg, controls)
    sub = sub.set_index("condition").reindex(order).dropna(subset=["positive_fraction"])
    if sub.empty:
        return None
    palette = condition_palette(sub.index, controls)

    values = 100 * sub["positive_fraction"].to_numpy(float)
    lo = 100 * sub["positive_fraction_ci_lo"].to_numpy(float)
    hi = 100 * sub["positive_fraction_ci_hi"].to_numpy(float)
    err = np.vstack([np.clip(values - lo, 0, None), np.clip(hi - values, 0, None)])
    err = np.nan_to_num(err, nan=0.0)

    fig, ax = plt.subplots(figsize=(max(6.0, 0.6 * len(sub) + 2.0), 4.4))
    x = np.arange(len(sub))
    ax.bar(x, values, color=[palette[c] for c in sub.index], alpha=0.85,
           yerr=err, error_kw={"lw": 1.0, "ecolor": style.ANNOTATION})

    fsub = fields[fields["marker"] == marker]
    rng = np.random.default_rng(0)
    for i, cond in enumerate(sub.index):
        vals = 100 * fsub.loc[fsub["condition"] == cond, "positive_fraction"].to_numpy(float)
        vals = vals[np.isfinite(vals)]
        if vals.size:
            ax.scatter(np.full(vals.size, i) + rng.uniform(-0.16, 0.16, vals.size), vals,
                       s=9, color="black", alpha=0.75, zorder=3, linewidths=0)

    ax.set_xticks(x)
    ax.set_xticklabels(shorten_labels(sub.index), rotation=40, ha="right")
    ax.set_ylabel(f"{marker}-positive cells (%)")
    ax.set_ylim(0, min(100, max(105, float(np.nanmax(values)) * 1.25)))
    unit = str(sub["replicate_unit"].iloc[0]) if "replicate_unit" in sub and len(sub) else "field"
    ncol = "n_replicates" if "n_replicates" in sub else "n_fields"
    ax.set_title(f"{marker}-positive fraction (error bars: 95% CI over {plural(unit)})"
                 + (f"\n{facet}" if facet else ""), pad=16)
    annotate_n(ax, sub[ncol].fillna(0).astype(int).to_list(), offset=0, unit=unit)
    return save_figure(fig, out_path)


def band_composition_bars(bands: pd.DataFrame, marker: str, cfg: dict, out_path,
                          controls=()):
    """Stacked neg/low/mid/high composition per condition."""
    sub = bands[bands["marker"] == marker]
    if sub.empty:
        return None
    order = ordered_conditions(sub["condition"].unique(), cfg, controls)
    pivot = (sub.pivot_table(index="condition", columns="band", values="fraction",
                             observed=True, aggfunc="first")
             .reindex(order).fillna(0.0))
    present = [b for b in BAND_ORDER if b in pivot.columns]
    if not present:
        return None

    fig, ax = plt.subplots(figsize=(max(6.0, 0.6 * len(pivot) + 2.4), 4.4))
    bottom = np.zeros(len(pivot))
    x = np.arange(len(pivot))
    for band in present:
        vals = 100 * pivot[band].to_numpy(float)
        ax.bar(x, vals, bottom=bottom, color=BAND_COLORS.get(band, style.GREY),
               label=band, edgecolor="white", linewidth=0.4)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(shorten_labels(pivot.index), rotation=40, ha="right")
    ax.set_ylabel("Cells (%)")
    ax.set_ylim(0, 100)
    ax.set_title(f"{marker}: graded composition (tertiles within the positive population)")
    ax.legend(title="band", bbox_to_anchor=(1.02, 1), loc="upper left")
    return save_figure(fig, out_path)


def coexpression_grid(coex: pd.DataFrame, marker_a: str, marker_b: str,
                      cfg: dict, out_path, controls=()):
    """One band × band heatmap per condition (the graded co-expression readout)."""
    sub = coex[(coex["marker_a"] == marker_a) & (coex["marker_b"] == marker_b)]
    if sub.empty:
        return None
    order = ordered_conditions(sub["condition"].unique(), cfg, controls)
    bands_a = [b for b in BAND_ORDER if b in set(sub["band_a"].astype(str))]
    bands_b = [b for b in BAND_ORDER if b in set(sub["band_b"].astype(str))]

    nrows, ncols = grid_shape(len(order), max_cols=int((cfg.get("visualization") or {})
                                                      .get("panel_cols", 4)))
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.7 * ncols, 2.55 * nrows),
                             squeeze=False)
    vmax = float(np.nanpercentile(100 * sub["fraction"].to_numpy(float), 99)) or 1.0
    # Panel titles are the only label these subplots carry, so they must stay
    # distinct from one another; shorten_labels widens until they are.
    titles = dict(zip(order, shorten_labels(order, 24)))

    for ax, cond in zip(axes.ravel(), order):
        cell = sub[sub["condition"] == cond]
        if cell.empty:
            no_data(ax, titles[cond])
            continue
        mat = (cell.pivot_table(index="band_b", columns="band_a", values="fraction",
                                observed=True, aggfunc="first")
               .reindex(index=bands_b, columns=bands_a).fillna(0.0).to_numpy() * 100)
        im = ax.imshow(mat, cmap=style.CMAP_SEQUENTIAL, vmin=0, vmax=vmax, origin="lower")
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                ax.text(j, i, f"{mat[i, j]:.0f}", ha="center", va="center", fontsize=6.5,
                        color=style.text_on_cmap(mat[i, j], 0, vmax, style.CMAP_SEQUENTIAL))
        ax.set_xticks(range(len(bands_a)), bands_a, fontsize=6.5, rotation=40, ha="right")
        ax.set_yticks(range(len(bands_b)), bands_b, fontsize=6.5)
        # Two lines: long condition names otherwise run into the neighboring panel.
        ax.set_title(f"{titles[cond]}\nn={int(cell['n_cells'].sum()):,}", fontsize=7.5)

    for ax in axes.ravel()[len(order):]:
        ax.axis("off")
    for row in axes:
        row[0].set_ylabel(marker_b, fontsize=8.5)
    # Label the bottom-most *used* panel of each column, not just the last row:
    # a partly-filled final row would otherwise leave columns unlabeled.
    for col in range(ncols):
        for row in range(nrows - 1, -1, -1):
            if row * ncols + col < len(order):
                axes[row][col].set_xlabel(marker_a, fontsize=8.5)
                break
    fig.colorbar(im, ax=axes, shrink=0.6, label="% of cells")
    fig.suptitle(f"{marker_a} × {marker_b}: graded co-expression by condition", fontsize=11)
    return save_figure(fig, out_path)


def quadrant_bars(coex: pd.DataFrame, marker_a: str, marker_b: str, cfg: dict,
                  out_path, controls=()):
    """Double-positive / single-positive / double-negative stacked per condition."""
    sub = coex[(coex["marker_a"] == marker_a) & (coex["marker_b"] == marker_b)].copy()
    if sub.empty:
        return None
    a_pos = sub["band_a"].astype(str) != NEG
    b_pos = sub["band_b"].astype(str) != NEG
    labels = np.select(
        [a_pos & b_pos, a_pos & ~b_pos, ~a_pos & b_pos],
        [f"{marker_a}+ {marker_b}+", f"{marker_a}+ {marker_b}−", f"{marker_a}− {marker_b}+"],
        default=f"{marker_a}− {marker_b}−",
    )
    sub["quadrant"] = labels
    agg = sub.groupby(["condition", "quadrant"], observed=True)["fraction"].sum().reset_index()
    order = ordered_conditions(agg["condition"].unique(), cfg, controls)
    pivot = agg.pivot_table(index="condition", columns="quadrant", values="fraction",
                            observed=True, aggfunc="sum").reindex(order).fillna(0.0)

    quad_order = [f"{marker_a}+ {marker_b}+", f"{marker_a}+ {marker_b}−",
                  f"{marker_a}− {marker_b}+", f"{marker_a}− {marker_b}−"]
    colors = [style.OKABE_ITO["purple"], style.OKABE_ITO["blue"],
              style.OKABE_ITO["orange"], style.NEUTRAL_FILL]

    fig, ax = plt.subplots(figsize=(max(6.0, 0.6 * len(pivot) + 2.6), 4.4))
    bottom = np.zeros(len(pivot))
    x = np.arange(len(pivot))
    for quad, color in zip(quad_order, colors):
        if quad not in pivot.columns:
            continue
        vals = 100 * pivot[quad].to_numpy(float)
        ax.bar(x, vals, bottom=bottom, color=color, label=quad,
               edgecolor="white", linewidth=0.4)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(shorten_labels(pivot.index), rotation=40, ha="right")
    ax.set_ylabel("Cells (%)")
    ax.set_ylim(0, 100)
    ax.set_title(f"{marker_a} × {marker_b}: quadrant composition")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
    return save_figure(fig, out_path)


__all__ = ["positive_fraction_bars", "band_composition_bars", "coexpression_grid",
           "quadrant_bars", "BAND_COLORS"]
