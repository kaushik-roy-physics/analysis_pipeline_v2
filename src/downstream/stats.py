"""
Statistics on field-level values.

Everything here takes an array of **per-field** observations, never per-cell
values: see :mod:`downstream.aggregate` for why. All routines are vectorized
over the bootstrap dimension (no Python resampling loop).
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def bootstrap_ci(values, n_boot: int = 2000, alpha: float = 0.05,
                 seed: int = 0, statistic=np.mean) -> tuple[float, float]:
    """Percentile bootstrap CI of *statistic* over field-level *values*.

    Resampling is done as one ``(n_boot, n)`` index matrix, so the whole bootstrap
    is a single vectorized reduction.
    """
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    n = values.size
    if n == 0:
        return (np.nan, np.nan)
    if n == 1:
        return (float(values[0]), float(values[0]))
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n, size=(n_boot, n))
    stats = statistic(values[idx], axis=1)
    lo, hi = np.percentile(stats, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return (float(lo), float(hi))


def cohens_d(a, b) -> float:
    """Standardized difference between two field-level samples."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    if a.size < 2 or b.size < 2:
        return np.nan
    pooled = np.sqrt(((a.size - 1) * a.var(ddof=1) + (b.size - 1) * b.var(ddof=1))
                     / (a.size + b.size - 2))
    return float((a.mean() - b.mean()) / pooled) if pooled > 0 else np.nan


def mann_whitney(a, b) -> float:
    """Two-sided Mann-Whitney U p-value; NaN when either sample is too small."""
    from scipy.stats import mannwhitneyu
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    if a.size < 3 or b.size < 3:
        return np.nan
    try:
        return float(mannwhitneyu(a, b, alternative="two-sided").pvalue)
    except ValueError:
        return np.nan


def benjamini_hochberg(pvalues) -> np.ndarray:
    """BH-adjusted q-values, NaNs preserved in place."""
    p = np.asarray(pvalues, dtype=float)
    out = np.full(p.shape, np.nan)
    finite = np.isfinite(p)
    if not finite.any():
        return out
    vals = p[finite]
    order = np.argsort(vals)
    ranked = vals[order]
    n = ranked.size
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    adjusted = np.empty_like(q)
    adjusted[order] = np.clip(q, 0, 1)
    out[finite] = adjusted
    return out


def smallest_p(n_a: int, n_b: int, exact: bool, n_perm: int) -> float:
    """The smallest p a permutation test of these sizes can return (a property of
    the design, not of the data): 5 vs 5 cannot go below 2/252 ≈ 0.008."""
    from math import comb
    if n_a < 1 or n_b < 1:
        return np.nan
    if not exact:
        return 1.0 / (n_perm + 1)
    total = comb(n_a + n_b, n_a)
    return (2.0 if n_a == n_b else 1.0) / total


def _replicates(frame, value_col, unit, cfg, need=2):
    """Replicate values at *unit*, stepping down (experiment → well → field) until
    every group has ≥ *need*. Returns ``(values, unit_used)``."""
    from .aggregate import UNIT_LADDER, replicate_values
    groups = [c for c in ("_group",) if c in frame.columns]
    start = UNIT_LADDER.index(unit) if unit in UNIT_LADDER else len(UNIT_LADDER) - 1
    for level in UNIT_LADDER[start:]:
        if level == "experiment" and not (cfg.get("stats") or {}).get("experiment_column"):
            continue
        try:
            vals = replicate_values(frame, groups, level, cfg, value_col)
        except Exception:                                   # noqa: BLE001
            continue
        if len(vals) and vals.groupby("_group").size().min() >= need \
                and vals["_group"].nunique() == frame["_group"].nunique():
            return vals, level
    return replicate_values(frame, groups, "field", cfg, value_col), "field"


def _two_sample(frame_a, frame_b, value_col, unit, cfg, n_perm, rng) -> dict:
    """Effect of a vs b at the replicate unit (Welch CI, Hedges' g, permutation p)."""
    from .effects import effect
    both = pd.concat([frame_a.assign(_group="a"), frame_b.assign(_group="b")],
                     ignore_index=True)
    vals, used = _replicates(both, value_col, unit, cfg)
    a = vals.loc[vals["_group"] == "a", value_col].to_numpy(float)
    b = vals.loc[vals["_group"] == "b", value_col].to_numpy(float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    res = effect(a, b, n_perm, rng)
    return {"unit": used, "pseudo_replicated": used != unit, "n": int(a.size),
            "n_reference": int(b.size), "difference": res["diff"], "ci_lo": res["ci_lo"],
            "ci_hi": res["ci_hi"], "hedges_g": res["hedges_g"], "cohens_d": cohens_d(a, b),
            "p_value": res["p_perm"], "exact": bool(res["exact"]),
            "smallest_possible_p": smallest_p(a.size, b.size, bool(res["exact"]), n_perm)}


def compare_to_reference(fields: pd.DataFrame, value_col: str, reference: str,
                         group_col: str = "condition", extra_keys=("marker",),
                         seed: int = 0, unit: str = "field", cfg: dict | None = None,
                         n_perm: int = 5000) -> pd.DataFrame:
    """Effect size + permutation test of every condition against a reference.

    Replicates are at *unit* (field, well or experiment; see
    :func:`downstream.aggregate.choose_unit`); a side with fewer than two
    replicates at that unit steps down a level and the row says so
    (``pseudo_replicated``). p-values are BH-adjusted within each marker.
    """
    cfg = cfg or {}
    rng = np.random.default_rng(seed)
    keys = [k for k in extra_keys if k in fields.columns]
    rows: list[dict] = []
    for key_vals, grp in fields.groupby(keys, observed=True) if keys else [((), fields)]:
        key_vals = key_vals if isinstance(key_vals, tuple) else (key_vals,)
        ref = grp[grp[group_col].astype(str) == str(reference)]
        for name, sub in grp.groupby(group_col, observed=True):
            if str(name) == str(reference) or ref.empty:
                continue
            vals = sub[value_col].to_numpy(float)
            refv = ref[value_col].to_numpy(float)
            rows.append({
                **dict(zip(keys, key_vals)),
                group_col: name, "value": value_col,
                "n_fields": int(np.isfinite(vals).sum()),
                "n_fields_reference": int(np.isfinite(refv).sum()),
                "mean": float(np.nanmean(vals)) if np.isfinite(vals).any() else np.nan,
                "mean_reference": float(np.nanmean(refv)) if np.isfinite(refv).any() else np.nan,
                **_two_sample(sub, ref, value_col, unit, cfg, n_perm, rng),
            })
    out = pd.DataFrame(rows)
    if len(out):
        out["q_value"] = np.nan
        for _, idx in out.groupby(keys, observed=True).groups.items() if keys else [(None, out.index)]:
            out.loc[idx, "q_value"] = benjamini_hochberg(out.loc[idx, "p_value"])
    return out


def factor_effects(fields: pd.DataFrame, value_col: str, factors: list[str],
                   extra_keys=("marker",), seed: int = 0, unit: str = "field",
                   cfg: dict | None = None, n_perm: int = 5000) -> pd.DataFrame:
    """Marginal effect of each 2-level factor on a readout, across replicates at *unit*.

    For a combinatorial screen (``bmp``, ``noggin``, ``fgf8``, … as separate
    plate-map columns) this is the natural first look: the average difference a
    factor makes across every background it appears in.
    """
    cfg = cfg or {}
    rng = np.random.default_rng(seed)
    keys = [k for k in extra_keys if k in fields.columns]
    rows: list[dict] = []
    groups = fields.groupby(keys, observed=True) if keys else [((), fields)]
    for key_vals, grp in groups:
        key_vals = key_vals if isinstance(key_vals, tuple) else (key_vals,)
        for factor in factors:
            if factor not in grp.columns:
                continue
            levels = pd.unique(grp[factor].dropna())
            if len(levels) != 2:
                continue
            lo, hi = sorted(levels, key=_level_key)
            a, b = grp[grp[factor] == hi], grp[grp[factor] == lo]
            res = _two_sample(a, b, value_col, unit, cfg, n_perm, rng)
            rows.append({
                **dict(zip(keys, key_vals)),
                "factor": factor, "level_high": str(hi), "level_low": str(lo),
                "value": value_col,
                "n_fields_high": int(np.isfinite(a[value_col].to_numpy(float)).sum()),
                "n_fields_low": int(np.isfinite(b[value_col].to_numpy(float)).sum()),
                "mean_high": float(np.nanmean(a[value_col])) if len(a) else np.nan,
                "mean_low": float(np.nanmean(b[value_col])) if len(b) else np.nan,
                "effect": res.pop("difference"), **res,
            })
    out = pd.DataFrame(rows)
    if len(out):
        # Multiplicity is corrected within each (scope × marker) family (the set
        # of factors tested on one readout), not globally across unrelated ones.
        if keys:
            out["q_value"] = np.nan
            for _, idx in out.groupby(keys, observed=True).groups.items():
                out.loc[idx, "q_value"] = benjamini_hochberg(out.loc[idx, "p_value"])
        else:
            out["q_value"] = benjamini_hochberg(out["p_value"])
    return out


def _level_key(v):
    """Numbers sort as numbers (0 < 5 < 100), text alphabetically."""
    try:
        return (0, float(v), "")
    except (TypeError, ValueError):
        return (1, 0.0, str(v))


def detect_factors(cells: pd.DataFrame, candidates: list[str] | None = None,
                   max_levels: int = 2) -> list[str]:
    """Plate-map columns that look like experimental factors (exactly 2 levels)."""
    from .data import ID_COLUMNS
    cols = candidates or [c for c in cells.columns if c not in ID_COLUMNS]
    out = []
    for col in cols:
        if col not in cells.columns:
            continue
        levels = pd.unique(cells[col].dropna())
        if 2 <= len(levels) <= max_levels:
            out.append(col)
    return out


__all__ = ["bootstrap_ci", "cohens_d", "mann_whitney", "benjamini_hochberg",
           "compare_to_reference", "factor_effects", "detect_factors", "smallest_p"]
