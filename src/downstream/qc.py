"""
Quality control: run this before believing any figure.

Two questions, both answered at the level of the imaged field:

1. **Did segmentation and staining behave?** Cell count, nuclear size and DAPI
   level per field. A field with a tenth of its neighbors' nuclei is a focus or
   segmentation failure; one with double is usually a debris-laden or
   over-confluent field.
2. **Is the well uniform?** Field-to-field spread of every readout within a well
   (item 6). A high coefficient of variation means the well average is hiding
   real heterogeneity, and any single field would have told a different story.

3. **Did each marker channel behave, field by field?**
   :func:`marker_field_qc`. The two checks above look only at DAPI, count and
   nuclear size: they say nothing at all about a marker channel, so an artifact
   confined to one antibody (non-specific patches, a laser flare, a bubble)
   passes them untouched. This one compares each field's marker median against
   **the other fields of its own well**, which is the only comparison where
   disagreement is necessarily technical: one well is one staining reaction on
   one population, so its fields should agree, whereas across conditions they
   are *supposed* to differ. Where the extractor's channel-health scalars are
   available they ride along, adding a clipping and a flooding flag.

Outlier fields are **flagged, never silently dropped**: the flags are written to
``qc_fields.csv`` / ``qc_marker_fields.csv`` and shown in the QC figure, and
dropping them is an explicit opt-in (``qc.drop_outlier_fields: true``, which acts
on the field-level flags only; a marker-level problem is not a reason to discard
that field's other, unaffected markers). This matters because thresholds are
estimated from cells pooled over a whole plate: one blown-out field can shift a
cut that then propagates into every well, so QC deliberately runs *before*
thresholding in the step order.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from .aggregate import condition_columns, condition_label
from .data import DAPI_NORM, RAW, Data

logger = logging.getLogger(__name__)

# Columns identifying the pool a field is compared against: its own well, in its
# own panel, on its own plate.
_WELL_KEYS = ("dataset", "plate_id", "stain_set", "well")


def mad_z(values: np.ndarray) -> np.ndarray:
    """Robust z-score: (x − median) / (1.4826 · MAD). Zero MAD → all zeros."""
    values = np.asarray(values, dtype=float)
    median = np.nanmedian(values)
    mad = np.nanmedian(np.abs(values - median))
    scale = 1.4826 * mad
    if not np.isfinite(scale) or scale == 0:
        return np.zeros_like(values)
    return (values - median) / scale


def field_qc(data: Data, cfg: dict) -> pd.DataFrame:
    """Per-field segmentation/staining metrics with robust outlier flags."""
    qc_cfg = cfg.get("qc") or {}
    z_max = float(qc_cfg.get("outlier_z", 3.5))
    min_cells = int(qc_cfg.get("min_cells_per_field", 20))

    cond_cfg = cfg.get("conditions") or {}
    keys = condition_columns(data.cells, cond_cfg)
    group_cols = [c for c in ("dataset", "seg_mode", "plate_id", "stain_set", "well", "row", "col")
                  if c in data.cells.columns]

    grouped = data.cells.groupby("field_id", observed=True)
    qc = grouped.agg(
        n_cells=("cell_uid", "size"),
        dapi_median=("dapi", "median"),
        area_median=("area", "median"),
        area_cv=("area", lambda s: float(np.nanstd(s) / np.nanmean(s))
                 if np.nanmean(s) else np.nan),
    ).reset_index()

    meta = grouped[group_cols + keys].first().reset_index()
    qc = qc.merge(meta, on="field_id", how="left")
    qc["condition"] = condition_label(qc, keys)

    # Robust z within each condition: comparing a field to its own biological
    # peers, not to the whole plate (a genuinely dense condition is not an outlier).
    # Peers share the dataset and segmentation mode too: a 2-D MIP and a 3-D
    # z-stack of the same field differ several-fold in count and DAPI level, and
    # pooling them would inflate the spread and hide real outliers.
    peers = [c for c in ("dataset", "seg_mode") if c in qc.columns] + ["condition"]
    grp = qc[peers].astype(str).agg("|".join, axis=1)
    for col in ("n_cells", "dapi_median", "area_median"):
        qc[f"{col}_z"] = (qc.groupby(grp)[col]
                          .transform(lambda s: pd.Series(mad_z(s.to_numpy()), index=s.index)))

    qc["flag_low_count"] = qc["n_cells"] < min_cells
    qc["flag_outlier"] = (
        (qc["n_cells_z"].abs() > z_max)
        | (qc["dapi_median_z"].abs() > z_max)
        | (qc["area_median_z"].abs() > z_max)
        | qc["flag_low_count"]
    )

    n_flagged = int(qc["flag_outlier"].sum())
    logger.info("qc: %d field(s), %d flagged (%.1f%%)", len(qc), n_flagged,
                100.0 * n_flagged / max(len(qc), 1))
    if n_flagged:
        worst = qc.loc[qc["flag_outlier"], ["field_id", "well", "condition", "n_cells"]].head(10)
        logger.info("qc: first flagged fields:\n%s", worst.to_string(index=False))
    return qc


def marker_field_qc(data: Data, cfg: dict,
                    channel_qc: pd.DataFrame | None = None) -> pd.DataFrame:
    """Per (field × marker) health, judged against the other fields of the well.

    Three independent flags, none of which changes a single analyzed number:

    ``flag_marker_outlier``
        The field's median for this marker is a robust-z outlier among the other
        fields of the same well × panel **and** deviates from the well's median
        by more than ``qc.marker_min_rel_dev``. Catches localized artifacts
        (patches of non-specific antibody, debris, a flare) that leave DAPI and
        the nucleus count perfectly normal and so pass :func:`field_qc` untouched.

        Both conditions are required because a robust z alone is unusable at the
        3-6 fields per well that most plates carry: the MAD of a handful of
        tightly agreeing numbers collapses toward zero, and then a field sitting
        3 % off its neighbors scores a z of 8. The relative floor is what keeps
        the flag pointed at artifacts rather than at ordinary agreement.
    ``flag_saturated``
        More than ``qc.max_saturated_fraction`` of that channel's pixels sat at
        the detector ceiling.
    ``flag_low_signal``
        Nuclei under ``qc.min_signal_to_background`` × their own background, i.e.
        indistinguishable from what they sit on.

    A note on what this can and cannot see: the robust z needs peers, so a well
    imaged at 2 fields cannot produce an outlier flag (the MAD of two numbers
    gives every point a z of 0) and one imaged at 3 will only catch the gross
    cases. That is a property of the design, not a failure of the check: with
    few fields per well the saturation and signal-to-background flags, which
    judge each field on its own, are what is doing the work.
    """
    qc_cfg = cfg.get("qc") or {}
    z_max = float(qc_cfg.get("marker_outlier_z", qc_cfg.get("outlier_z", 3.5)))
    min_rel = float(qc_cfg.get("marker_min_rel_dev", 0.25))
    max_sat = float(qc_cfg.get("max_saturated_fraction", 0.001))
    min_sb = float(qc_cfg.get("min_signal_to_background", 1.05))

    keys = [c for c in _WELL_KEYS if c in data.cells.columns]
    long = data.attach(keys + ["field_id"])

    # Each marker is judged in the space it is actually analyzed in, so the flag
    # and the figure it would explain are talking about the same numbers.
    spaces = {str(m): data.space_for(str(m), cfg) for m in data.markers}
    space_per_row = long["marker"].astype(str).map(spaces).to_numpy()
    values = np.where(space_per_row == DAPI_NORM,
                      long[DAPI_NORM].to_numpy(dtype=float),
                      long[RAW].to_numpy(dtype=float))
    long = long.assign(_value=values)

    per = (long.groupby(keys + ["field_id", "marker"], observed=True)
           .agg(n_cells=("cell_uid", "size"), median=("_value", "median")).reset_index())
    per["space"] = per["marker"].astype(str).map(spaces)

    per["median_z"] = (per.groupby(keys + ["marker"], observed=True)["median"]
                       .transform(lambda s: pd.Series(mad_z(s.to_numpy()), index=s.index)))
    per["n_fields_in_well"] = (per.groupby(keys + ["marker"], observed=True)["field_id"]
                               .transform("nunique"))
    well_median = (per.groupby(keys + ["marker"], observed=True)["median"]
                   .transform("median"))
    per["rel_dev"] = ((per["median"] - well_median).abs()
                      / well_median.abs().replace(0, np.nan))
    per["flag_marker_outlier"] = (per["median_z"].abs() > z_max) & (per["rel_dev"] > min_rel)

    if channel_qc is not None and len(channel_qc):
        cols = [c for c in ("field_id", "marker", "saturated_fraction",
                            "nuclei_saturated_fraction", "background",
                            "signal_to_background") if c in channel_qc.columns]
        merged = channel_qc[cols].copy()
        merged["marker"] = merged["marker"].astype(str)
        per["marker"] = per["marker"].astype(str)
        per = per.merge(merged, on=["field_id", "marker"], how="left")
        per["flag_saturated"] = per.get("saturated_fraction", np.nan) > max_sat
        per["flag_low_signal"] = per.get("signal_to_background", np.nan) < min_sb
    else:
        per["flag_saturated"] = False
        per["flag_low_signal"] = False

    per["flag_any"] = (per["flag_marker_outlier"] | per["flag_saturated"]
                       | per["flag_low_signal"])

    n_flag = int(per["flag_any"].sum())
    logger.info("qc: %d (field × marker) row(s), %d flagged (%.1f%%)",
                len(per), n_flag, 100.0 * n_flag / max(len(per), 1))
    if n_flag:
        show = [c for c in ("field_id", "well", "marker", "median", "median_z",
                            "saturated_fraction", "signal_to_background") if c in per.columns]
        worst = (per[per["flag_any"]].sort_values("median_z", key=abs, ascending=False)
                 .head(10)[show])
        logger.warning("qc: marker-channel problems, worst first:\n%s",
                       worst.to_string(index=False))
        for marker, grp in per[per["flag_any"]].groupby("marker", observed=True):
            logger.warning("qc: %s; %d suspect field(s); check them before trusting "
                           "its gate, especially any in a control well.",
                           marker, len(grp))
    return per


def well_uniformity(fields: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Field-to-field spread within each well: the heterogeneity readout (item 6)."""
    keys = [c for c in ("dataset", "plate_id", "stain_set", "well", "condition", "marker")
            if c in fields.columns]
    out = fields.groupby(keys, observed=True).agg(
        n_fields=("field_id", "nunique"),
        n_cells=("n_cells", "sum"),
        median_mean=("median", "mean"),
        median_sd=("median", "std"),
        median_min=("median", "min"),
        median_max=("median", "max"),
        positive_fraction_mean=("positive_fraction", "mean"),
        positive_fraction_sd=("positive_fraction", "std"),
    ).reset_index()
    out["median_cv"] = out["median_sd"] / out["median_mean"].replace(0, np.nan)
    out["positive_fraction_range"] = (
        fields.groupby(keys, observed=True)["positive_fraction"]
        .agg(lambda s: float(np.nanmax(s) - np.nanmin(s)) if s.notna().any() else np.nan)
        .to_numpy()
    )
    high = out["median_cv"] > float((cfg.get("qc") or {}).get("cv_warn", 0.5))
    if high.any():
        logger.warning("qc: %d well × marker combination(s) with field CV > %.2f; "
                       "well means are hiding real heterogeneity.",
                       int(high.sum()), float((cfg.get("qc") or {}).get("cv_warn", 0.5)))
    return out


def outlier_field_ids(qc: pd.DataFrame, cfg: dict) -> set:
    """Fields left out of the analysis: flagged ones, when ``qc.drop_outlier_fields``."""
    if not (cfg.get("qc") or {}).get("drop_outlier_fields", False) or "flag_outlier" not in qc:
        return set()
    return set(qc.loc[qc["flag_outlier"].astype(bool), "field_id"])


def apply_field_filter(data: Data, qc: pd.DataFrame, cfg: dict) -> Data:
    """Optionally drop flagged fields from the working tables (off by default)."""
    if not (cfg.get("qc") or {}).get("drop_outlier_fields", False):
        return data
    bad = set(qc.loc[qc["flag_outlier"], "field_id"])
    if not bad:
        return data
    keep_cells = ~data.cells["field_id"].isin(bad)
    cells = data.cells[keep_cells].reset_index(drop=True)
    kept_uids = set(cells["cell_uid"])
    intens = data.intensities[data.intensities["cell_uid"].isin(kept_uids)].reset_index(drop=True)
    # Re-index so cell_uid stays positional (the invariant Data.attach relies on).
    remap = pd.Series(np.arange(len(cells), dtype=np.int64), index=cells["cell_uid"].to_numpy())
    cells["cell_uid"] = remap.to_numpy()
    intens["cell_uid"] = remap.reindex(intens["cell_uid"].to_numpy()).to_numpy()
    logger.warning("qc: dropped %d flagged field(s) → %d cells remain", len(bad), len(cells))
    return Data(cells=cells, intensities=intens)


def saturated_nuclei(data: Data, channel_qc: pd.DataFrame | None, cfg: dict) -> pd.DataFrame:
    """Nuclei whose median intensity sits at the detector ceiling, in any channel.

    A nucleus whose *median* pixel is clipped is not a cell reading high: it is
    bright debris or a flare, and it stretches scales and gates. A nucleus is
    flagged when its raw median is ≥ ``qc.saturation_fraction`` (default 1.0) ×
    the ceiling of that channel in that image: ``qc.saturation_ceiling`` when set
    (e.g. 4095 for a 12-bit camera writing 16-bit files), else the extraction's
    per-channel ``saturation_value`` (the dtype maximum). The nuclear channel
    (the ``dapi`` column) counts too. Float images have no ceiling and are never
    flagged.

    Returns one row per flagged nucleus × channel: ``cell_uid``, ``field_id``,
    ``channel_name``, ``value``, ``ceiling``. Empty when nothing is clipped.
    """
    cols = ["cell_uid", "field_id", "channel_name", "value", "ceiling"]
    qc_cfg = cfg.get("qc") or {}
    frac = float(qc_cfg.get("saturation_fraction", 1.0))
    override = qc_cfg.get("saturation_ceiling")
    cells = data.cells
    field = cells["field_id"].to_numpy()
    markers = {str(m) for m in data.markers}
    have_qc = channel_qc is not None and len(channel_qc) and \
        {"field_id", "marker", "saturation_value"} <= set(channel_qc.columns)
    if override is None and not have_qc:
        return pd.DataFrame(columns=cols)
    if have_qc:
        cq = channel_qc[["field_id", "marker", "saturation_value"]].copy()
        cq["marker"] = cq["marker"].astype(str)
        ceil_of = cq.drop_duplicates(["field_id", "marker"]).set_index(["field_id", "marker"])[
            "saturation_value"]
        # The nuclear channel is the channel-QC row that is not a quantified marker.
        nuc = cq[~cq["marker"].isin(markers)].groupby("field_id")["saturation_value"].min()
    frames = []
    inten = data.intensities
    if RAW in inten:
        uid = inten["cell_uid"].to_numpy()
        names = inten["marker"].astype(str).to_numpy()
        if override is not None:
            ceil = np.full(len(inten), float(override))
        else:
            idx = pd.MultiIndex.from_arrays([field[uid], names])
            ceil = ceil_of.reindex(idx).to_numpy(np.float64)
        value = inten[RAW].to_numpy(np.float64)
        with np.errstate(invalid="ignore"):
            hit = np.isfinite(ceil) & (value >= frac * ceil)
        frames.append(pd.DataFrame({"cell_uid": uid[hit], "field_id": field[uid][hit],
                                    "channel_name": names[hit], "value": value[hit],
                                    "ceiling": ceil[hit]}))
    if "dapi" in cells:
        value = cells["dapi"].to_numpy(np.float64)
        ceil = (np.full(len(cells), float(override)) if override is not None
                else pd.Series(field).map(nuc).to_numpy(np.float64))
        with np.errstate(invalid="ignore"):
            hit = np.isfinite(ceil) & (value >= frac * ceil)
        frames.append(pd.DataFrame({"cell_uid": cells["cell_uid"].to_numpy()[hit],
                                    "field_id": field[hit], "channel_name": "nuclear",
                                    "value": value[hit], "ceiling": ceil[hit]}))
    out = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=cols)
    return out[cols]


__all__ = ["field_qc", "marker_field_qc", "well_uniformity", "apply_field_filter",
           "outlier_field_ids", "saturated_nuclei", "mad_z"]
