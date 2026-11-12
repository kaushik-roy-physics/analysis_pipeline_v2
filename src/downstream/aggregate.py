"""
Rollups: cell → field → condition.

**Cells are never the replicate.** A well holds thousands of cells, so treating
cells as replicates makes every comparison "significant" no matter how small the
effect. Each field is first reduced to one value per marker; which level is then
the replicate (``n``) is chosen by :func:`choose_unit`:

* ``experiment`` when ``stats.experiment_column`` names independent
  differentiations and every condition has ≥ ``stats.min_experiments_for_unit``
  (2) of them: fields are averaged per well, wells per experiment;
* ``well`` when every condition has ≥ ``stats.min_wells_for_well_unit`` (3)
  wells: fields of one well share a seeding and a staining reaction, so they are
  averaged first;
* ``field`` otherwise, and every table and the report say what that means (the
  fields of one well are partly dependent; with one well per condition, well and
  condition are confounded).

``stats.replicate`` (auto | field | well | experiment) overrides the choice.
Cell-level distributions are still plotted (they are the biology), but they
never set an error bar.

Three tables come out of here and everything downstream reads them:

``fields.csv``      one row per (field × marker): n cells, median/mean intensity,
                    positive fraction, plus the field's well and conditions.
``conditions.csv``  one row per (condition × marker): mean over fields ± SD, SEM
                    and a bootstrap CI, the field count, and the pooled cell count.
``wells.csv``       the same at well level: the unit you pipette, and the level
                    at which plate-position artifacts show up.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from . import controls, scope
from .data import Data
from .stats import bootstrap_ci
from .effects import unit_values

logger = logging.getLogger(__name__)

# Field-level identity carried into every rollup.
FIELD_KEYS = ["field_id", "dataset", "plate_id", "seg_mode", "stain_set", "well", "row", "col"]


def condition_columns(cells: pd.DataFrame, cond_cfg: dict) -> list[str]:
    """The plate-map columns that define a condition (config-declared)."""
    keys = cond_cfg.get("condition_key") or ["treatment"]
    keys = [keys] if isinstance(keys, str) else list(keys)
    present = [k for k in keys if k in cells.columns]
    missing = [k for k in keys if k not in cells.columns]
    if missing:
        logger.warning("conditions.condition_key names absent column(s) %s: ignoring. "
                       "Available: %s", missing, sorted(set(cells.columns)))
    if not present:
        raise ValueError(
            f"None of conditions.condition_key={keys} exist in the manifest. "
            f"Set it to one or more of: {sorted(set(cells.columns))}"
        )
    return present


UNIT_LADDER = ["experiment", "well", "field"]


def add_well_key(frame: pd.DataFrame) -> pd.DataFrame:
    """A ``well_key`` naming one physical well (dataset / plate / well)."""
    if "well_key" in frame.columns:
        return frame
    cols = [c for c in ("dataset", "plate_id", "well") if c in frame.columns]
    out = frame.copy()
    out["well_key"] = (out[cols].astype(str).agg("|".join, axis=1) if cols
                       else out["field_id"].astype(str))
    return out


def experiment_column(frame: pd.DataFrame, cfg: dict) -> str | None:
    """``stats.experiment_column`` when set and present in *frame*, else None."""
    col = (cfg.get("stats") or {}).get("experiment_column")
    return str(col) if col and str(col) in frame.columns else None


def choose_unit(fields: pd.DataFrame, cfg: dict) -> tuple[str, str]:
    """``(unit, note)``: the replicate for every interval and test (see module doc)."""
    st = cfg.get("stats") or {}
    mode = str(st.get("replicate", "auto")).lower()
    fields = add_well_key(fields)
    per = fields.groupby("condition", observed=True)["well_key"].nunique()
    fewest = int(per.min()) if len(per) else 0
    exp_col = experiment_column(fields, cfg)
    fewest_e = 0
    if exp_col:
        per_e = fields.groupby("condition", observed=True)[exp_col].nunique()
        fewest_e = int(per_e.min()) if len(per_e) else 0
    min_wells = int(st.get("min_wells_for_well_unit", 3))
    min_exp = int(st.get("min_experiments_for_unit", 2))
    if mode in ("field", "well", "experiment"):
        unit = mode
        if unit == "experiment" and not exp_col:
            raise ValueError("stats.replicate: experiment needs stats.experiment_column naming "
                             "a plate-map column present in the manifest")
    elif exp_col and fewest_e >= min_exp:
        unit = "experiment"
    else:
        unit = "well" if fewest >= min_wells else "field"
    if unit == "experiment":
        return unit, (f"Replicate = experiment ('{exp_col}': every condition has ≥ {fewest_e} "
                      "independent differentiations); fields are averaged per well and wells "
                      "per experiment, so n = experiments.")
    if unit == "well":
        return unit, (f"Replicate = well (every condition has ≥ {fewest} wells); fields are "
                      "averaged within each well first, so n = wells.")
    if fewest <= 1:
        return unit, ("Replicate = field. At least one condition sits in a single well, so "
                      "well-to-well variation cannot be separated from the condition: a "
                      "difference shows these fields differ, not that it would replicate in "
                      "another well. Add wells (≥ 3 per condition) to make the well the "
                      "replicate.")
    return unit, (f"Replicate = field (some condition has only {fewest} wells, fewer than "
                  f"{min_wells}). Fields of one well share a seeding and a staining reaction, "
                  "so the intervals and p-values are optimistic.")


def replicate_values(frame: pd.DataFrame, by: list[str], unit: str, cfg: dict,
                     value: str) -> pd.DataFrame:
    """One row per replicate at *unit*: fields as they are, wells = mean of their
    fields, experiments = mean of their well means. ``by`` + unit_id + *value*."""
    frame = add_well_key(frame)
    if unit == "field":
        f = frame[np.isfinite(frame[value].to_numpy(np.float64))]
        out = f[by + [value]].copy()
        out["unit_id"] = f["field_id"].astype(str).to_numpy()
        out["n_fields"] = 1
        return out
    out = unit_values(frame, by, unit, experiment_column(frame, cfg), value=value)
    return out.rename(columns={"n_fields": "n_fields"})


def condition_label(df: pd.DataFrame, keys: list[str]) -> pd.Series:
    """Compose the condition label from one or several columns."""
    if len(keys) == 1:
        return df[keys[0]].astype(str)
    return df[keys].astype(str).agg(" | ".join, axis=1)


def field_table(data: Data, gated: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Per (field × marker) summary: the replicate-level table."""
    cond_cfg = cfg.get("conditions") or {}
    keys = condition_columns(data.cells, cond_cfg)

    long = gated.copy()
    long["field_id"] = data.cell_column("field_id", long["cell_uid"].to_numpy())

    grouped = long.groupby(["field_id", "marker"], observed=True)
    out = grouped.agg(
        n_cells=("value", "size"),
        median=("value", "median"),
        mean=("value", "mean"),
        p25=("value", lambda s: np.nanpercentile(s, 25) if s.notna().any() else np.nan),
        p75=("value", lambda s: np.nanpercentile(s, 75) if s.notna().any() else np.nan),
        positive_fraction=("positive", "mean"),
        # min_count=1 so "no gate for this marker" stays NaN instead of summing to 0.
        n_positive=("positive", lambda s: s.sum(min_count=1)),
    ).reset_index()

    # Field metadata: one row per field, taken from the first cell of that field.
    # The control-declaring columns travel along so figures can gray out controls
    # without re-reading the cell table.
    control_col = cond_cfg.get("control_column", "class")
    candidates = [control_col, controls.CONTROL_FOR, controls.POSITIVE_FOR]
    extra = [c for c in dict.fromkeys(candidates)
             if c in data.cells.columns and c not in keys]
    meta_cols = [c for c in FIELD_KEYS if c in data.cells.columns and c != "field_id"]
    meta = (data.cells.groupby("field_id", observed=True)[meta_cols + keys + extra]
            .first().reset_index())
    meta["condition"] = condition_label(meta, keys)

    out = out.merge(meta, on="field_id", how="left")

    # Each row's gate scope, as the single key every per-marker rollup and test
    # groups by. Markers may sit at different scope levels, so one shared tuple of
    # scope columns cannot express this: see :mod:`downstream.scope`.
    available = [c for c in scope.UNION_ORDER if c in out.columns]
    out["scope_key"] = scope.scope_keys(out, cfg, available)

    logger.info("fields: %d field × marker rows over %d field(s)",
                len(out), out["field_id"].nunique())
    return out


def _component_columns(fields: pd.DataFrame, cfg: dict) -> list[str]:
    """The plate-map columns the condition label was composed from, if present."""
    try:
        return condition_columns(fields, cfg.get("conditions") or {})
    except ValueError:          # a reloaded table may carry only `condition`
        return []


def well_table(fields: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Per (well × marker) summary: fields aggregated within each well."""
    keys = [c for c in ("dataset", "plate_id", "stain_set", "seg_mode", "scope_key",
                        "well", "row", "col") if c in fields.columns]
    keys += [c for c in _component_columns(fields, cfg) if c in fields.columns]
    keys += [c for c in ("condition", "marker") if c in fields.columns]
    grouped = fields.groupby(keys, observed=True)
    out = grouped.agg(
        n_fields=("field_id", "nunique"),
        n_cells=("n_cells", "sum"),
        median=("median", "mean"),
        median_sd=("median", "std"),
        positive_fraction=("positive_fraction", "mean"),
        positive_fraction_sd=("positive_fraction", "std"),
    ).reset_index()
    # Field-to-field heterogeneity within the well (item 6).
    out["median_cv"] = out["median_sd"] / out["median"].replace(0, np.nan)
    logger.info("wells: %d well × marker rows", len(out))
    return out


def condition_table(fields: pd.DataFrame, cfg: dict, unit: str = "field") -> pd.DataFrame:
    """Per (condition × marker) summary; error bars and CIs across replicates at *unit*."""
    stats_cfg = cfg.get("stats") or {}
    n_boot = int(stats_cfg.get("n_bootstrap", 2000))
    seed = int(stats_cfg.get("random_state", 0))
    min_fields = int(stats_cfg.get("min_fields", 2))

    # Group by each marker's OWN gate scope, not by the union of everyone's scope
    # columns. A marker pooled across panels (`markers: {ISL1: {scope: plate}}`)
    # has one threshold and must get one row per condition, aggregating every
    # field the gate was estimated over; grouping on the `stain_set` its rows
    # still carry would split that condition back into one row per panel, each
    # with a third of the fields behind its CI. The pooled columns are marked `*`,
    # the same convention thresholds.csv uses.
    fields = scope.mark_pooled(fields, cfg)

    # The condition's component columns ride along as their own columns. They are
    # functionally determined by `condition` (it is their join), so adding them to
    # the group keys cannot change the grouping, but it makes the table pivotable
    # by factor and lets figures facet on one of them.
    keys = [c for c in ("dataset", "plate_id", "stain_set", "seg_mode", "scope_key")
            if c in fields.columns]
    keys += [c for c in _component_columns(fields, cfg) if c in fields.columns]
    keys += [c for c in ("condition", "marker") if c in fields.columns]
    rows: list[dict] = []
    for key_vals, grp in fields.groupby(keys, observed=True):
        key_vals = key_vals if isinstance(key_vals, tuple) else (key_vals,)
        rec = dict(zip(keys, key_vals))
        med = replicate_values(grp, [], unit, cfg, "median")["median"].to_numpy(float)
        frac = replicate_values(grp, [], unit, cfg, "positive_fraction")[
            "positive_fraction"].to_numpy(float)
        n_rep = int(max(np.isfinite(med).sum(), np.isfinite(frac).sum()))
        rec.update(
            replicate_unit=unit, n_replicates=n_rep,
            n_fields=int(len(grp)),
            n_wells=int(grp["well"].nunique()) if "well" in grp else np.nan,
            n_cells=int(grp["n_cells"].sum()),
            median=float(np.nanmean(med)) if np.isfinite(med).any() else np.nan,
            median_sd=float(np.nanstd(med, ddof=1)) if np.isfinite(med).sum() > 1 else np.nan,
            positive_fraction=float(np.nanmean(frac)) if np.isfinite(frac).any() else np.nan,
            positive_fraction_sd=(float(np.nanstd(frac, ddof=1))
                                  if np.isfinite(frac).sum() > 1 else np.nan),
        )
        rec["median_sem"] = (rec["median_sd"] / np.sqrt(max(n_rep, 1))
                             if np.isfinite(rec["median_sd"]) else np.nan)
        rec["positive_fraction_sem"] = (rec["positive_fraction_sd"] / np.sqrt(max(n_rep, 1))
                                        if np.isfinite(rec["positive_fraction_sd"]) else np.nan)
        if n_rep >= min_fields:
            lo, hi = bootstrap_ci(frac, n_boot=n_boot, seed=seed)
            rec["positive_fraction_ci_lo"], rec["positive_fraction_ci_hi"] = lo, hi
            lo, hi = bootstrap_ci(med, n_boot=n_boot, seed=seed)
            rec["median_ci_lo"], rec["median_ci_hi"] = lo, hi
        else:
            rec.update(positive_fraction_ci_lo=np.nan, positive_fraction_ci_hi=np.nan,
                       median_ci_lo=np.nan, median_ci_hi=np.nan)
        rows.append(rec)

    out = pd.DataFrame(rows)
    logger.info("conditions: %d condition × marker rows (median %d fields each)",
                len(out), int(out["n_fields"].median()) if len(out) else 0)
    return out


__all__ = ["field_table", "well_table", "condition_table", "condition_columns",
           "condition_label", "FIELD_KEYS", "choose_unit", "replicate_values", "add_well_key",
           "experiment_column", "UNIT_LADDER"]
