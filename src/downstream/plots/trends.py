"""
Trends across an ordered axis, and factor effects in combinatorial designs.

:func:`timecourse` plots any readout against any numeric plate-map column
(``day``, ``timepoint``, ``dose_uM``, ``reagent_hours``), with one line per
series column. This is what makes "different treatments on different days, with
several fixed timepoints" a config setting rather than new code::

    trends:
      x: timepoint
      series: treatment
      facet: day

:func:`factor_effects_plot` is for factorial screens where conditions are encoded
as independent boolean columns (``bmp``, ``noggin``, ``fgf8``, …): it shows each
factor's average marginal effect with its confidence interval, so the drivers
are visible without reading 47 individual bars.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .. import style
from .base import grid_shape, no_data, save_figure, shorten, shorten_labels

logger = logging.getLogger(__name__)


def timecourse(fields: pd.DataFrame, marker: str, value_col: str, x_col: str,
               series_col: str | None, cfg: dict, out_path, facet_col: str | None = None,
               percent: bool = False):
    """Field-level *value_col* against *x_col*, one line per series, optional facets."""
    sub = fields[fields["marker"] == marker].copy()
    if sub.empty or x_col not in sub.columns:
        return None
    sub[x_col] = pd.to_numeric(sub[x_col], errors="coerce")
    sub = sub.dropna(subset=[x_col, value_col])
    if sub.empty or sub[x_col].nunique() < 2:
        return None

    facets = ([None] if not facet_col or facet_col not in sub.columns
              else list(dict.fromkeys(sub[facet_col].astype(str))))
    nrows, ncols = grid_shape(len(facets), max_cols=3)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.6 * ncols, 3.6 * nrows),
                             squeeze=False, sharey=True)
    scale = 100.0 if percent else 1.0

    # A numeric x-axis is what makes an ordered categorical readable as a trend,
    # but 0/1/2/3 means nothing to a reader. `trends.x_tick_labels` puts the
    # category names back on the axis without giving up the ordering.
    tick_labels = (cfg.get("trends") or {}).get("x_tick_labels") or {}

    for ax, facet in zip(axes.ravel(), facets):
        part = sub if facet is None else sub[sub[facet_col].astype(str) == facet]
        series = ([None] if not series_col or series_col not in part.columns
                  else list(dict.fromkeys(part[series_col].astype(str))))
        names = shorten_labels([s for s in series if s is not None], 20)
        legend = dict(zip([s for s in series if s is not None], names))
        for i, name in enumerate(series):
            block = part if name is None else part[part[series_col].astype(str) == name]
            grouped = block.groupby(x_col, observed=True)[value_col]
            mean = scale * grouped.mean()
            sem = scale * grouped.std(ddof=1) / np.sqrt(grouped.count().clip(lower=1))
            sty = style.series_style(i)
            ax.errorbar(mean.index.to_numpy(dtype=float), mean.to_numpy(),
                        yerr=sem.to_numpy(), capsize=2.5, markersize=4.5,
                        label=legend.get(name, marker), **sty)
        if tick_labels:
            ticks = sorted(float(v) for v in part[x_col].dropna().unique())
            ax.set_xticks(ticks)
            ax.set_xticklabels(
                [str(tick_labels.get(int(t) if float(t).is_integer() else t,
                                     tick_labels.get(t, t))) for t in ticks],
                rotation=30, ha="right")
        ax.set_xlabel(x_col)
        if facet is not None:
            ax.set_title(f"{facet_col} = {facet}", fontsize=9)
    for ax in axes.ravel()[len(facets):]:
        ax.axis("off")
    for row in axes:
        row[0].set_ylabel(f"{marker} {value_col}" + (" (%)" if percent else ""))
    axes.ravel()[0].legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7.5)
    fig.suptitle(f"{marker}: {value_col} across {x_col} "
                 f"(mean ± SEM over fields)", fontsize=11)
    fig.tight_layout()
    return save_figure(fig, out_path)


def factor_effects_plot(effects: pd.DataFrame, marker: str, value_col: str,
                        cfg: dict, out_path):
    """Marginal effect of each factor on one readout, sorted by magnitude."""
    sub = effects[(effects["marker"] == marker) & (effects["value"] == value_col)].copy()
    sub = sub.dropna(subset=["effect"])
    if sub.empty:
        return None
    sub = sub.sort_values("effect")

    fig, ax = plt.subplots(figsize=(6.6, max(2.6, 0.42 * len(sub) + 1.6)))
    y = np.arange(len(sub))
    significant = sub["q_value"].to_numpy(dtype=float) < 0.05
    colors = [style.OKABE_ITO["blue"] if s else style.GREY for s in significant]
    ax.barh(y, sub["effect"].to_numpy(float), color=colors, alpha=0.9)
    ax.axvline(0, color="black", lw=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{f} ({lo}→{hi})" for f, lo, hi in
                        zip(sub["factor"], sub["level_low"], sub["level_high"])], fontsize=8)
    ax.set_xlabel(f"Δ {marker} {value_col} (factor present − absent)")
    ax.set_title(f"{marker}: marginal factor effects\n"
                 "colored = BH q < 0.05, field-level test", fontsize=10)
    for i, (d, q) in enumerate(zip(sub["cohens_d"], sub["q_value"])):
        if np.isfinite(d):
            ax.text(sub["effect"].iloc[i], i, f"  d={d:.2f}", va="center", fontsize=6.5,
                    color=style.ANNOTATION)
    return save_figure(fig, out_path)


def interaction_heatmap(fields: pd.DataFrame, marker: str, value_col: str,
                        factor_a: str, factor_b: str, cfg: dict, out_path,
                        percent: bool = False):
    """Mean readout for every combination of two factors: the interaction view."""
    sub = fields[fields["marker"] == marker]
    if sub.empty or factor_a not in sub.columns or factor_b not in sub.columns:
        return None
    scale = 100.0 if percent else 1.0
    pivot = sub.pivot_table(index=factor_b, columns=factor_a, values=value_col,
                            observed=True, aggfunc="mean")
    counts = sub.pivot_table(index=factor_b, columns=factor_a, values="field_id",
                             observed=True, aggfunc="nunique")
    if pivot.empty:
        return None
    # pivot_table can hand back object dtype holding pd.NA; force a float matrix
    # so the NaN-aware maths below stays well defined.
    pivot = pivot.apply(pd.to_numeric, errors="coerce") * scale
    values = pivot.to_numpy(dtype=float)
    if not np.isfinite(values).any():
        return None

    vmin = float(np.nanmin(values))
    vmax = float(np.nanmax(values))
    fig, ax = plt.subplots(figsize=(1.5 * len(pivot.columns) + 2.6,
                                    1.2 * len(pivot.index) + 2.2))
    im = ax.imshow(values, cmap=style.CMAP_SEQUENTIAL, vmin=vmin, vmax=vmax)
    ax.set_xticks(range(len(pivot.columns)), [str(c) for c in pivot.columns])
    ax.set_yticks(range(len(pivot.index)), [str(i) for i in pivot.index])
    ax.set_xlabel(factor_a)
    ax.set_ylabel(factor_b)
    count_values = counts.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            val = values[i, j]
            if np.isfinite(val):
                n = count_values[i, j]
                label = f"{val:.1f}" + (f"\nn={int(n)}" if np.isfinite(n) else "")
                ax.text(j, i, label, ha="center", va="center", fontsize=7,
                        color=style.text_on_cmap(val, vmin, vmax, style.CMAP_SEQUENTIAL))
    fig.colorbar(im, ax=ax, shrink=0.8,
                 label=f"{marker} {value_col}" + (" (%)" if percent else ""))
    ax.set_title(f"{marker}: {factor_a} × {factor_b}", fontsize=10)
    return save_figure(fig, out_path)


__all__ = ["timecourse", "factor_effects_plot", "interaction_heatmap"]
