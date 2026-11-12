"""
Two-marker joint views.

:func:`pair_density` draws one small-multiple panel per condition: a 2-D
histogram of marker A against marker B with the positivity gates and the tertile
band edges drawn on top, so the co-expression matrix computed in
:mod:`downstream.gating` can be read straight off the density.

Density rather than a raw scatter is the default because a well holds thousands
of cells and an alpha-blended scatter saturates long before it shows structure.
:func:`pair_scatter_colored` draws the same grid as a cell-level scatter; when
the panel has more than two markers it colors every cell by a third marker's
intensity, one figure per third marker, so a three-way relationship can be read
from a single figure. :func:`pair_scatter` is kept for the small-n case and uses
a shared axis range so panels remain comparable.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LogNorm

from .. import style
from .base import (grid_shape, no_data, ordered_conditions, save_figure, shorten,
                   shorten_labels)

logger = logging.getLogger(__name__)


def _limits(values: np.ndarray, cap_pct: float) -> tuple[float, float]:
    vals = values[np.isfinite(values)]
    if vals.size == 0:
        return (0.0, 1.0)
    hi = float(np.nanpercentile(vals, cap_pct))
    lo = float(max(0.0, np.nanmin(vals)))
    return (lo, hi if hi > lo else lo + 1.0)


def _gate_lines(ax, thr, edges, axis: str, over_density: bool = False):
    """Draw the positivity gate (dashed) and tertile edges (dotted) on one axis.

    Over a density map the line lands on colors ranging from near-white to
    near-black, so a flat color is guaranteed to disappear somewhere in the
    panel. *over_density* adds a contrasting outline that keeps the gate (the
    single most important annotation on the figure) visible over all of them.
    """
    draw = ax.axvline if axis == "x" else ax.axhline
    effects = {"path_effects": style.halo(2.6)} if over_density else {}
    if thr is not None and np.isfinite(thr):
        draw(thr, color=style.THRESHOLD, lw=1.4, ls="--", zorder=4, **effects)
    for edge in (edges or []):
        if edge is not None and np.isfinite(edge):
            draw(edge, color=style.ANNOTATION, lw=1.0, ls=":", zorder=4,
                 **({"path_effects": style.halo(2.2)} if over_density else {}))


def pair_density(wide: pd.DataFrame, marker_a: str, marker_b: str,
                 conditions: pd.Series, gates: dict, cfg: dict, out_path):
    """2-D histogram of A vs B per condition, with gates and band edges overlaid.

    ``gates`` carries ``{'thr_a', 'thr_b', 'edges_a', 'edges_b'}``.
    """
    vis = cfg.get("visualization") or {}
    cap_pct = float(vis.get("cap_percentile", 99.5))
    bins = int(vis.get("density_bins", 60))

    frame = wide.copy()
    frame["condition"] = conditions.reindex(frame.index).to_numpy()
    frame = frame.dropna(subset=[marker_a, marker_b, "condition"])
    if frame.empty:
        return None

    order = ordered_conditions(frame["condition"].unique(), cfg)
    xlim = _limits(frame[marker_a].to_numpy(float), cap_pct)
    ylim = _limits(frame[marker_b].to_numpy(float), cap_pct)

    nrows, ncols = grid_shape(len(order), max_cols=int(vis.get("panel_cols", 4)))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.0 * ncols, 2.85 * nrows),
                             squeeze=False, sharex=True, sharey=True)

    xbins = np.linspace(*xlim, bins)
    ybins = np.linspace(*ylim, bins)
    titles = dict(zip(order, shorten_labels(order, 24)))
    # One shared color scale across panels: per-panel autoscaling would make a
    # sparse well look as dense as a crowded one.
    counts_max = 1.0
    for cond in order:
        sub = frame[frame["condition"] == cond]
        pair = np.isfinite(sub[marker_a].to_numpy(float)) & np.isfinite(sub[marker_b].to_numpy(float))
        if pair.sum() < 10:
            continue
        hist, _, _ = np.histogram2d(sub[marker_a].to_numpy(float)[pair],
                                    sub[marker_b].to_numpy(float)[pair],
                                    bins=[xbins, ybins])
        counts_max = max(counts_max, float(hist.max()))
    norm = LogNorm(vmin=1, vmax=max(counts_max, 2.0))

    mesh = None
    for ax, cond in zip(axes.ravel(), order):
        sub = frame[frame["condition"] == cond]
        x = sub[marker_a].to_numpy(float)
        y = sub[marker_b].to_numpy(float)
        keep = np.isfinite(x) & np.isfinite(y)
        if keep.sum() < 10:
            no_data(ax, titles[cond])
            continue
        *_, mesh = ax.hist2d(x[keep], y[keep], bins=[xbins, ybins],
                             cmap=style.CMAP_DENSITY, norm=norm, rasterized=True)
        _gate_lines(ax, gates.get("thr_a"), gates.get("edges_a"), "x", over_density=True)
        _gate_lines(ax, gates.get("thr_b"), gates.get("edges_b"), "y", over_density=True)
        ax.set_title(f"{titles[cond]}\nn={int(keep.sum()):,}", fontsize=7.5)

    for ax in axes.ravel()[len(order):]:
        ax.axis("off")
    for row in axes:
        row[0].set_ylabel(marker_b)
    # Bottom-most used panel per column (the final row may be partly empty).
    for col in range(ncols):
        for row in range(nrows - 1, -1, -1):
            if row * ncols + col < len(order):
                axes[row][col].set_xlabel(marker_a)
                break
    fig.suptitle(f"{marker_a} vs {marker_b}: per-cell density "
                 f"(dashed = positivity gate, dotted = tertile edges)", fontsize=11)
    fig.tight_layout()
    if mesh is not None:
        fig.colorbar(mesh, ax=axes, shrink=0.55, label="cells per bin")
    return save_figure(fig, out_path)


def pair_scatter_colored(wide: pd.DataFrame, marker_a: str, marker_b: str,
                         conditions: pd.Series, gates: dict, cfg: dict, out_path,
                         color_by: str | None = None, thr_c: float | None = None):
    """Per-condition scatter of A vs B, optionally colored by a third marker C.

    Same panel grid, axis ranges and gate overlays as :func:`pair_density`, so the
    two figures can be laid side by side. With *color_by* each cell takes its C
    intensity on the sequential map, on ONE color scale shared by every panel
    (1st percentile → ``cap_percentile``); the colorbar carries C's positivity
    gate as a dashed line. Cells are drawn in ascending C so the C-high cells are
    not buried under the C-low bulk. Only cells measured for all plotted markers
    are shown. Without *color_by* (a two-marker panel) the points are one color.
    """
    vis = cfg.get("visualization") or {}
    cap_pct = float(vis.get("cap_percentile", 99.5))
    size = float(vis.get("joint_scatter_size", 2))
    alpha = float(vis.get("joint_scatter_alpha", 0.8))

    markers = [marker_a, marker_b] + ([color_by] if color_by else [])
    frame = wide.copy()
    frame["condition"] = conditions.reindex(frame.index).to_numpy()
    frame = frame.dropna(subset=markers + ["condition"])
    frame = frame[np.isfinite(frame[markers].to_numpy(float)).all(axis=1)]
    if len(frame) < 10:
        return None

    order = ordered_conditions(frame["condition"].unique(), cfg)
    xlim = _limits(frame[marker_a].to_numpy(float), cap_pct)
    ylim = _limits(frame[marker_b].to_numpy(float), cap_pct)
    norm = None
    if color_by:
        cvals = frame[color_by].to_numpy(float)
        lo = float(np.nanpercentile(cvals, 1))
        hi = float(np.nanpercentile(cvals, cap_pct))
        norm = plt.Normalize(vmin=lo, vmax=hi if hi > lo else lo + 1.0)

    nrows, ncols = grid_shape(len(order), max_cols=int(vis.get("panel_cols", 4)))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.0 * ncols, 2.85 * nrows),
                             squeeze=False, sharex=True, sharey=True)
    titles = dict(zip(order, shorten_labels(order, 24)))

    points = None
    for ax, cond in zip(axes.ravel(), order):
        sub = frame[frame["condition"] == cond]
        if len(sub) < 10:
            no_data(ax, titles[cond])
            continue
        if color_by:
            sub = sub.sort_values(color_by, kind="stable")
            points = ax.scatter(sub[marker_a], sub[marker_b], c=sub[color_by],
                                cmap=style.CMAP_SEQUENTIAL, norm=norm, s=size,
                                alpha=alpha, linewidths=0, rasterized=True)
        else:
            ax.scatter(sub[marker_a], sub[marker_b], color=style.color(0), s=size,
                       alpha=alpha, linewidths=0, rasterized=True)
        _gate_lines(ax, gates.get("thr_a"), gates.get("edges_a"), "x", over_density=True)
        _gate_lines(ax, gates.get("thr_b"), gates.get("edges_b"), "y", over_density=True)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_title(f"{titles[cond]}\nn={len(sub):,}", fontsize=7.5)

    for ax in axes.ravel()[len(order):]:
        ax.axis("off")
    for row in axes:
        row[0].set_ylabel(marker_b)
    for col in range(ncols):
        for row in range(nrows - 1, -1, -1):
            if row * ncols + col < len(order):
                axes[row][col].set_xlabel(marker_a)
                break
    head = f"{marker_a} vs {marker_b}"
    if color_by:
        head += f", colored by {color_by} (high-{color_by} cells drawn on top)"
    fig.suptitle(f"{head}\ndashed = positivity gate, dotted = tertile edges",
                 fontsize=11)
    fig.tight_layout()
    if points is not None:
        extend = ("both" if (cvals < norm.vmin).any() and (cvals > norm.vmax).any()
                  else "max" if (cvals > norm.vmax).any()
                  else "min" if (cvals < norm.vmin).any() else "neither")
        # An opaque mappable, so the bar (and its extend triangles) show the true
        # colors rather than the points' alpha.
        opaque = ScalarMappable(norm=norm, cmap=style.CMAP_SEQUENTIAL)
        cbar = fig.colorbar(opaque, ax=axes, shrink=0.55, extend=extend,
                            label=f"{color_by} intensity")
        if thr_c is not None and np.isfinite(thr_c) and norm.vmin <= thr_c <= norm.vmax:
            cbar.ax.axhline(thr_c, color=style.THRESHOLD, lw=1.6, ls="--",
                            path_effects=style.halo(2.6))
    return save_figure(fig, out_path)


def pair_scatter(wide: pd.DataFrame, marker_a: str, marker_b: str,
                 conditions: pd.Series, gates: dict, cfg: dict, out_path):
    """Single-panel scatter colored by condition: for few conditions / few cells."""
    vis = cfg.get("visualization") or {}
    cap_pct = float(vis.get("cap_percentile", 99.5))
    size = float(vis.get("scatter_size", 4))
    alpha = float(vis.get("scatter_alpha", 0.35))

    frame = wide.copy()
    frame["condition"] = conditions.reindex(frame.index).to_numpy()
    frame = frame.dropna(subset=[marker_a, marker_b, "condition"])
    if frame.empty:
        return None
    order = ordered_conditions(frame["condition"].unique(), cfg)

    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    legend_labels = shorten_labels(order, 24)
    for i, cond in enumerate(order):
        sub = frame[frame["condition"] == cond]
        ax.scatter(sub[marker_a], sub[marker_b], s=size, alpha=alpha,
                   color=style.color(i), marker=style.MARKERS[i % len(style.MARKERS)],
                   linewidths=0, label=legend_labels[i], rasterized=True)
    _gate_lines(ax, gates.get("thr_a"), gates.get("edges_a"), "x")
    _gate_lines(ax, gates.get("thr_b"), gates.get("edges_b"), "y")
    ax.set_xlim(*_limits(frame[marker_a].to_numpy(float), cap_pct))
    ax.set_ylim(*_limits(frame[marker_b].to_numpy(float), cap_pct))
    ax.set_xlabel(marker_a)
    ax.set_ylabel(marker_b)
    ax.set_title(f"{marker_a} vs {marker_b}")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", markerscale=2.5)
    return save_figure(fig, out_path)


def marker_correlation(wide: pd.DataFrame, cfg: dict, out_path, title: str = ""):
    """Spearman correlation matrix of the markers in one panel."""
    if wide.shape[1] < 2:
        return None
    corr = wide.corr(method="spearman")
    values = corr.to_numpy(dtype=float).copy()

    # The diagonal is 1.0 by construction. Left in, it owns both the color scale
    # and the eye; the off-diagonal, the only part carrying information, is
    # squeezed into a sliver of the ramp. Mask it out and scale to what is left.
    off = values[~np.eye(len(corr), dtype=bool)]
    off = off[np.isfinite(off)]
    # Symmetric about zero: for a correlation, "no relationship" must stay the
    # neutral midpoint of a diverging map, so only the SPAN is rescaled.
    limit = float(np.nanmax(np.abs(off))) if off.size else 1.0
    limit = max(limit, 0.1)
    shown = values.copy()
    np.fill_diagonal(shown, np.nan)

    cmap = plt.get_cmap(style.CMAP_DIVERGING).copy()
    cmap.set_bad(style.MASKED)

    fig, ax = plt.subplots(figsize=(1.1 * len(corr) + 2.4, 1.0 * len(corr) + 2.0))
    im = ax.imshow(shown, cmap=cmap, vmin=-limit, vmax=limit)
    ax.set_xticks(range(len(corr)))
    ax.set_xticklabels(corr.columns, rotation=40, ha="right")
    ax.set_yticks(range(len(corr)))
    ax.set_yticklabels(corr.index)
    for i in range(len(corr)):
        for j in range(len(corr)):
            if i == j:
                continue
            val = values[i, j]
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7.5,
                    color=style.text_on_cmap(val, -limit, limit, style.CMAP_DIVERGING))
    fig.colorbar(im, ax=ax, shrink=0.8, label="Spearman ρ")
    ax.set_title(title or "Marker correlation (per cell)\ndiagonal masked; "
                          f"scale ±{limit:.2f}", fontsize=10)
    return save_figure(fig, out_path)


__all__ = ["pair_density", "pair_scatter_colored", "pair_scatter", "marker_correlation"]
