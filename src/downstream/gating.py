"""
Graded gating: negative / low / mid / high bands, and marker co-expression.

A binary positive/negative call throws away most of what an intensity
distribution says. This module splits the **positive** population of each marker
into tertiles, giving four ordered bands per marker::

    neg   below the positivity threshold
    low   bottom third of the positive population
    mid   middle third
    high  top third

Band edges are the 33rd/67th percentiles of the positive cells **pooled over the
test wells of that marker's scope** (that marker's *own* controls excluded, so the
cut is not dragged by a population that is negative by construction). Pooling
makes the edges a single absolute reference shared by every condition, which is
what makes bands comparable across treatments: per-condition tertiles would
compare each condition only against itself. ``split: median`` gives a 50/50
low/high variant.

Both qualifiers matter. The scope is the marker's own (:mod:`downstream.scope`),
so a plate-scoped marker gets one set of edges across every panel it appears in
rather than two incomparable sets. And the exclusion is per marker
(:mod:`downstream.controls`): a pluripotent well that is SOX2's *positive* control
must stay in SOX2's positive pool: dropping it would compute SOX2's tertiles
after discarding the brightest cells in the experiment and silently bias every
band downward.

Crossing two markers' bands yields the 4×4 co-expression matrix: for every
condition, the fraction of cells in each (marker A band × marker B band) cell.
That is the "what fraction co-express, and in which part of the distribution"
readout, and it degrades gracefully to the familiar 2×2 quadrant view by
collapsing the three positive bands.
"""

from __future__ import annotations

import itertools
import logging

import numpy as np
import pandas as pd

from . import scope as scope_mod
from .aggregate import condition_columns, condition_label
from .controls import Controls, control_mask
from .data import Data
from .thresholds import OK, SCOPE_COLUMNS

logger = logging.getLogger(__name__)

NEG = "neg"
LOW = "low"
MID = "mid"
HIGH = "high"
BAND_ORDER = [NEG, LOW, MID, HIGH]
BAND_ORDER_MEDIAN = [NEG, LOW, HIGH]


def positive_edges(gated: pd.DataFrame, controls, cfg: dict) -> pd.DataFrame:
    """Tertile (or median) edges of each (scope × marker) positive population.

    *controls* is a :class:`downstream.controls.Controls`; a plain boolean array is
    still accepted and treated as one marker-agnostic mask.
    """
    gate_cfg = cfg.get("gating") or {}
    split = str(gate_cfg.get("split", "tertile"))
    min_positive = int(gate_cfg.get("min_positive_cells", 50))

    available = [c for c in SCOPE_COLUMNS if c in gated.columns]
    markers = list(pd.unique(gated["marker"]))
    union_cols = scope_mod.union_scope_columns(markers, cfg, available)

    rows: list[dict] = []
    for marker in markers:
        marker = str(marker)
        scope_cols = scope_mod.scope_columns_for(marker, cfg, available)
        block = gated[gated["marker"] == marker].copy()
        block["_is_control"] = _mask_for(controls, marker)[block["cell_uid"].to_numpy()]

        groups = (block.groupby(scope_cols, observed=True) if scope_cols
                  else [((), block)])
        for keys, grp in groups:
            keys = keys if isinstance(keys, tuple) else (keys,)
            scope = dict(zip(scope_cols, keys))
            rec = {
                **scope_mod.fill_scope_row(scope, union_cols),
                "marker": marker,
                "scope_key": ("|".join(str(scope[c]) for c in scope_cols)
                              if scope_cols else ""),
            }
            # Test wells only: this marker's own negative control contributes no
            # meaningful positives. A well that is some *other* marker's control
            # stays in: it may be this marker's brightest population.
            test = grp[~grp["_is_control"]]
            pos = test.loc[test["positive"].to_numpy(dtype=float) > 0, "value"]
            pos = pos.to_numpy(dtype=float)
            pos = pos[np.isfinite(pos)]
            rec["n_positive"] = int(pos.size)
            rec["split"] = split
            if pos.size < min_positive:
                rec.update(edge_lo=np.nan, edge_hi=np.nan, median=np.nan,
                           note=f"only {pos.size} positive cells (<{min_positive})")
            elif split == "median":
                med = float(np.median(pos))
                rec.update(edge_lo=med, edge_hi=med, median=med, note="median split")
            else:
                rec.update(edge_lo=float(np.quantile(pos, 1 / 3)),
                           edge_hi=float(np.quantile(pos, 2 / 3)),
                           median=float(np.median(pos)), note="tertile split")
            rows.append(rec)
    out = pd.DataFrame(rows)
    logger.info("gating: edges for %d (scope × marker) combination(s), split=%s",
                len(out), split)
    return out


def _mask_for(controls, marker: str) -> np.ndarray:
    """Per-marker control mask, tolerating a plain boolean array."""
    if isinstance(controls, Controls):
        return controls.mask_for(marker)
    return np.asarray(controls, dtype=bool)


def assign_bands(gated: pd.DataFrame, edges: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Add an ordered ``band`` column to the gated long table (vectorized)."""
    gate_cfg = cfg.get("gating") or {}
    split = str(gate_cfg.get("split", "tertile"))
    available = [c for c in SCOPE_COLUMNS if c in gated.columns]

    # Same (marker, scope_key) join as apply_thresholds, so band edges follow the
    # marker's declared scope rather than a fixed column tuple.
    gated_keys = scope_mod.scope_keys(gated, cfg, available)
    if "scope_key" in edges.columns:
        edge_keys = edges["scope_key"].fillna("").astype(str)
    else:
        edge_keys = scope_mod.scope_keys(
            edges, cfg, [c for c in SCOPE_COLUMNS if c in edges.columns])

    index = pd.MultiIndex.from_arrays([edges["marker"].astype(str), edge_keys])
    keys = pd.MultiIndex.from_arrays([gated["marker"].astype(str), gated_keys])
    lo = pd.Series(edges["edge_lo"].to_numpy(), index=index).reindex(keys).to_numpy(float)
    hi = pd.Series(edges["edge_hi"].to_numpy(), index=index).reindex(keys).to_numpy(float)

    values = gated["value"].to_numpy(dtype=float)
    positive = gated["positive"].to_numpy(dtype=float)
    called = np.isfinite(positive)          # NaN = marker has no gate → no band
    is_pos = called & (positive > 0)

    band = np.full(len(gated), None, dtype=object)
    band[called & ~is_pos] = NEG
    if split == "median":
        band[is_pos & (values < lo)] = LOW
        band[is_pos & (values >= lo)] = HIGH
    else:
        band[is_pos & (values < lo)] = LOW
        band[is_pos & (values >= lo) & (values < hi)] = MID
        band[is_pos & (values >= hi)] = HIGH
    # Positive cells in a scope with no usable edges stay uncategorized (NaN).
    order = BAND_ORDER_MEDIAN if split == "median" else BAND_ORDER
    out = gated.copy()
    out["band"] = pd.Categorical(band, categories=order, ordered=True)
    return out


def marker_pairs(data: Data, banded: pd.DataFrame, cfg: dict) -> list[tuple[str, str]]:
    """Marker pairs to cross: explicit from config, or every within-panel pair."""
    gate_cfg = cfg.get("gating") or {}
    spec = gate_cfg.get("pairs", "auto")
    if isinstance(spec, list) and spec:
        return [(str(a), str(b)) for a, b in spec]

    pairs: set[tuple[str, str]] = set()
    for _, grp in banded.groupby("stain_set", observed=True):
        markers = sorted(m for m in grp["marker"].unique())
        pairs.update(itertools.combinations(markers, 2))
    return sorted(pairs)


def coexpression(data: Data, banded: pd.DataFrame, cfg: dict,
                 pairs: list[tuple[str, str]] | None = None) -> pd.DataFrame:
    """Per-condition 4×4 band co-occurrence for each marker pair.

    Returns a long table: one row per (condition × pair × band_a × band_b) with
    the cell count, the fraction of that condition's paired cells, and the
    field-level mean fraction ± SD (the honest error bar).
    """
    cond_cfg = cfg.get("conditions") or {}
    keys = condition_columns(data.cells, cond_cfg)
    pairs = pairs if pairs is not None else marker_pairs(data, banded, cfg)
    if not pairs:
        return pd.DataFrame()

    # Each cell belongs to exactly one plate × panel, so pairing is scope-safe by
    # construction; the scope travels along so figures can be split by panel.
    scope_cols = [c for c in SCOPE_COLUMNS if c in data.cells.columns]
    uid = banded["cell_uid"].to_numpy()
    meta = pd.DataFrame({
        "cell_uid": uid,
        "field_id": data.cell_column("field_id", uid),
    })
    for col in scope_cols:
        meta[col] = data.cell_column(col, uid)
    for k in keys:
        meta[k] = data.cell_column(k, uid)
    meta["condition"] = condition_label(meta, keys)
    meta = meta.drop_duplicates("cell_uid").set_index("cell_uid")

    rows: list[pd.DataFrame] = []
    for marker_a, marker_b in pairs:
        sub = banded[banded["marker"].isin([marker_a, marker_b])]
        wide = sub.pivot_table(index="cell_uid", columns="marker", values="band",
                               observed=True, aggfunc="first")
        if marker_a not in wide.columns or marker_b not in wide.columns:
            continue
        wide = wide[[marker_a, marker_b]].dropna(how="any")
        if wide.empty:
            continue
        joined = wide.join(meta, how="inner")

        counts = (joined.groupby(scope_cols + ["condition", "field_id", marker_a, marker_b],
                                 observed=True).size().rename("n").reset_index())
        totals = counts.groupby(scope_cols + ["condition", "field_id"],
                                observed=True)["n"].transform("sum")
        counts["field_fraction"] = counts["n"] / totals

        agg = counts.groupby(scope_cols + ["condition", marker_a, marker_b],
                             observed=True).agg(
            n_cells=("n", "sum"),
            n_fields=("field_id", "nunique"),
            field_fraction=("field_fraction", "mean"),
            field_fraction_sd=("field_fraction", "std"),
        ).reset_index()
        agg = agg.rename(columns={marker_a: "band_a", marker_b: "band_b"})
        cond_totals = agg.groupby(scope_cols + ["condition"],
                                  observed=True)["n_cells"].transform("sum")
        agg["fraction"] = agg["n_cells"] / cond_totals
        agg.insert(0, "marker_a", marker_a)
        agg.insert(1, "marker_b", marker_b)
        rows.append(agg)

    out = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    logger.info("gating: co-expression for %d pair(s) → %d rows", len(pairs), len(out))
    return out


def band_composition(data: Data, banded: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Per (condition × marker) fraction of cells in each band, with field error bars."""
    cond_cfg = cfg.get("conditions") or {}
    keys = condition_columns(data.cells, cond_cfg)
    uid = banded["cell_uid"].to_numpy()

    scope_cols = [c for c in SCOPE_COLUMNS if c in banded.columns]
    df = banded[["marker", "band"] + scope_cols].copy()
    df["field_id"] = data.cell_column("field_id", uid)
    for k in keys:
        df[k] = data.cell_column(k, uid)
    df["condition"] = condition_label(df, keys)

    # Scope is part of the key: a marker measured in two panels keeps two
    # separate compositions rather than being averaged into a meaningless one.
    counts = (df.groupby(scope_cols + ["condition", "field_id", "marker", "band"],
                         observed=True).size().rename("n").reset_index())
    totals = counts.groupby(scope_cols + ["condition", "field_id", "marker"],
                            observed=True)["n"].transform("sum")
    counts["field_fraction"] = counts["n"] / totals

    out = counts.groupby(scope_cols + ["condition", "marker", "band"], observed=True).agg(
        n_cells=("n", "sum"),
        n_fields=("field_id", "nunique"),
        fraction=("field_fraction", "mean"),
        fraction_sd=("field_fraction", "std"),
    ).reset_index()
    return out


def gate_summary(gated: pd.DataFrame, thresholds: pd.DataFrame) -> pd.DataFrame:
    """Which markers actually carry a usable gate: a quick trust check."""
    ok = thresholds[thresholds["verdict"] == OK]
    return (ok.groupby("marker", observed=True)
            .agg(n_scopes=("threshold", "size"),
                 threshold_min=("threshold", "min"),
                 threshold_max=("threshold", "max"))
            .reset_index())


__all__ = ["positive_edges", "assign_bands", "coexpression", "band_composition",
           "marker_pairs", "gate_summary", "BAND_ORDER", "BAND_ORDER_MEDIAN",
           "NEG", "LOW", "MID", "HIGH", "control_mask"]
