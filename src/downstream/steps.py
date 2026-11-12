"""
The built-in analysis steps.

Each is a small orchestration function: pull what it needs off the context, call
into the computational modules, write tables and figures. The science lives in
:mod:`downstream.thresholds`, :mod:`downstream.gating`, :mod:`downstream.aggregate`
and :mod:`downstream.stats`; the drawing lives in :mod:`downstream.plots`.

Order (dependencies are declared, not assumed)::

    signaling → assemble → qc → thresholds → aggregate → gating
                                    ↘ distributions, composition, joint, plate, trends,
                                      effects → report

``signaling`` (per image) does nothing unless a stain set holds a signaling
marker (β-catenin, SMAD2/3, ...; :mod:`downstream.signaling`). Optional steps
(``downstream.skip_steps``) are listed in :data:`OPTIONAL_STEPS`.

``qc`` deliberately precedes ``thresholds``: a cut is estimated from cells pooled
over the whole plate, so a blown-out field must be visible (and optionally
excluded) *before* it can move a threshold that then propagates everywhere.
"""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd

from . import aggregate as agg_mod
from . import captions as cap_mod
from . import effects as ef_mod
from . import signaling as sig_mod
from . import assemble as assemble_mod
from . import controls as controls_mod
from . import gating as gating_mod
from . import plots
from . import qc as qc_mod
from . import stats as stats_mod
from . import thresholds as thr_mod
from .aggregate import condition_label
from .assemble import assemble as assemble_cells
from .data import ID_COLUMNS
from .registry import register

logger = logging.getLogger(__name__)


def _noted(func):
    """Draw the normalization caption under every figure the step saves."""
    import functools

    @functools.wraps(func)
    def wrapper(ctx):
        try:
            note = _caption(ctx)
        except Exception:                                   # noqa: BLE001
            logger.debug("caption unavailable", exc_info=True)
            note = ""
        with plots.figure_note(note):
            return func(ctx)
    return wrapper


# ---------------------------------------------------------------------------
# 1. Assembly
# ---------------------------------------------------------------------------

def _n_jobs(cfg: dict) -> int:
    """Worker count: ``signaling.n_jobs``, else the LSF slot count, else up to 8 cores."""
    want = (cfg.get("signaling") or {}).get("n_jobs")
    if want:
        return max(1, int(want))
    lsf = os.environ.get("LSB_DJOB_NUMPROC")
    if lsf and lsf.isdigit():
        return max(1, int(lsf))
    return max(1, min(8, os.cpu_count() or 1))


def _pixel_size(result_dir: Path, cfg: dict) -> float:
    """µm per pixel: ``signaling.pixel_size_um``, else the image metadata."""
    override = (cfg.get("signaling") or {}).get("pixel_size_um")
    if override:
        return float(override)
    meta = Path(result_dir) / "image_metadata.json"
    if meta.exists():
        px = json.loads(meta.read_text()).get("pixel_size_x_um")
        if px and float(px) > 0:
            return float(px)
    raise ValueError(f"{Path(result_dir).name}: no pixel size in the image metadata; set "
                     "downstream.signaling.pixel_size_um (µm per pixel).")


SIG_FIG = "signaling"


def _signaling_plan(ctx):
    """``(jobs, matches)``: per image, the signaling channels and their readouts."""
    aliases = {str(k).strip().lower(): str(v)
               for k, v in (ctx.cfg.get("marker_aliases") or {}).items()}
    exclude = {str(m).strip().lower() for m in (ctx.cfg.get("exclude_markers") or [])} | {"dapi"}
    jobs, names = [], []
    for ds, man in assemble_mod.iter_manifests(ctx.datasets, ctx.cfg):
        dapi = int((ds.cfg.get("channels") or {}).get("dapi_index", 0))
        measured = {int(c) for c in (ds.cfg.get("channels") or {}).get("marker_channels", [])}
        for _, row in man.iterrows():
            ch2marker = assemble_mod._marker_map(row, aliases)
            channels = {}
            for ch, marker in ch2marker.items():
                if ch == dapi or marker.strip().lower() in exclude:
                    continue
                if measured and ch not in measured:
                    continue
                names.append(marker)
                match = sig_mod.classify(marker, ctx.cfg)
                if match is not None:
                    channels[int(ch)] = match.readout
            if channels:
                panel = " | ".join(f"{c}:{m}" for c, m in sorted(ch2marker.items()))
                jobs.append((ds.name, ds.directory / str(row["result_dir"]), channels, panel))
    matches = sig_mod.log_resolution(sorted(set(names)), ctx.cfg)
    return jobs, matches


@register("signaling",
          config_keys=["signaling", "marker_aliases", "exclude_markers", "manifest_basename",
                       "_selection"], always=True,
          description="Signaling markers: nuclear core, cytoplasmic ring and background per "
                      "nucleus (only when a stain set holds one)")
def step_signaling(ctx) -> bool:
    from .signaling import measure as sig_measure
    jobs, matches = _signaling_plan(ctx)
    if not matches:
        logger.info("signaling: no signaling marker in any stain set; nothing to measure.")
        return False
    jobs = [j for j in jobs if (j[1] / "nuclear_masks.tif").exists()]
    panels = {path: panel for _, path, _, panel in jobs}
    jobs = [j[:3] for j in jobs]
    if not jobs:
        logger.warning("signaling: no segmented image holds a signaling marker yet.")
        return False
    params = sig_mod.settings(ctx.cfg)
    force = ctx.force or getattr(ctx, "step_forced", False)
    todo, pixel = [], {}
    for name, path, channels in jobs:
        pixel[path] = _pixel_size(path, ctx.cfg)
        digest = sig_measure.settings_hash(params, channels)
        if force or not sig_measure.is_current(path, digest):
            todo.append((name, path, channels))
    n_jobs = 1 if len(todo) < 8 else min(_n_jobs(ctx.cfg), len(todo))
    logger.info("signaling: %d image(s) with signaling channels, %d to (re)measure on %d "
                "worker(s)", len(jobs), len(todo), n_jobs)
    if n_jobs > 1:
        from joblib import Parallel, delayed
        results = Parallel(n_jobs=n_jobs, backend="loky")(
            delayed(sig_measure.worker)(str(path), ch, params, pixel[path], force)
            for _, path, ch in todo)
    else:
        results = [sig_measure.worker(str(path), ch, params, pixel[path], force)
                   for _, path, ch in todo]
    failures = [(name, path, err) for (name, path, _), (_, err) in zip(todo, results) if err]
    for name, path, err in failures:
        logger.error("[%s] signaling failed for %s: %s", name, path.name, err)
    if failures and len(failures) == len(todo):
        raise RuntimeError("signaling: every image failed; see the errors above.")

    bg = sig_measure.background_table((name, path, panels[path]) for name, path, _ in jobs)
    ctx.write_table(bg, "signaling_background.csv")
    for r in bg[bg["outlier"]].itertuples(index=False):
        logger.warning("signaling: %s channel %d background %.4g is far from the other images "
                       "of %s with the same stain set; check that image.", r.result_dir, r.channel,
                       r.background, r.dataset)
    ctx.reset_figures(SIG_FIG)
    entries = [(path / sig_measure.QC_NAME, plots.short_name(path.name), "")
               for _, path, _ in sorted(jobs, key=lambda j: str(j[1]))
               if (path / sig_measure.QC_NAME).exists()]
    per_page = 12
    pages = [entries[i:i + per_page] for i in range(0, len(entries), per_page)]
    for k, page in enumerate(pages, 1):
        plots.overlay_gallery(page, ctx.figure(SIG_FIG, f"ring_qc_{k:02d}.png"), ncols=3,
                              thumb_px=900,
                              title=f"Signaling rings: page {k}/{len(pages)} (orange = nucleus "
                                    f"core, sky blue = cytoplasmic ring inside cell bodies)")
    return bool(todo)


def _invalidate(ctx) -> None:
    for key in ("data", "all_data", "gated", "controls", "thresholds", "saturated_mask"):
        ctx._cache.pop(key, None)


@register("assemble", after=["signaling"],
          config_keys=["marker_aliases", "exclude_markers", "manifest_basename", "_selection",
                       "_manifests", "signaling"],
          description="Join manifests + per-image joblibs into cells/intensities tables")
def step_assemble(ctx) -> None:
    assemble_cells(ctx.datasets, ctx.cfg, ctx.out)
    _invalidate(ctx)


# ---------------------------------------------------------------------------
# 2. Quality control
# ---------------------------------------------------------------------------

@register("qc", after=["assemble"],
          config_keys=["qc", "conditions", "markers", "default_space", "signaling",
                       "nuclear_stain"],
          description="Per-field segmentation/staining QC, per-marker channel QC, outlier "
                      "flags, saturated nuclei")
def step_qc(ctx) -> None:
    _invalidate(ctx)
    full = ctx.all_data
    qc = qc_mod.field_qc(full, ctx.cfg)
    ctx.write_table(qc, "qc_fields.csv")
    ctx.set("qc", qc)
    plots.count_qc(qc, ctx.cfg, ctx.figure("qc", "field_counts.png"))

    # Per (field × marker) health. Separate from the table above because it is a
    # different question at a different grain: `qc_fields` asks whether the field
    # was imaged and segmented properly, this asks whether each antibody channel
    # in it can be believed. Runs here, before `thresholds`, so a suspect control
    # field is reported before it can move a gate.
    channel_qc = None
    channel_qc_path = ctx.out / assemble_mod.CHANNEL_QC_NAME
    if channel_qc_path.exists():
        channel_qc = pd.read_csv(channel_qc_path)
    marker_qc = qc_mod.marker_field_qc(full, ctx.cfg, channel_qc)
    ctx.write_table(marker_qc, "qc_marker_fields.csv")
    ctx.set("marker_qc", marker_qc)

    # A histone nuclear label (H2B-GFP) is not a DNA stain: dividing by it adds noise.
    names = _nuclear_names(ctx)
    kind = cap_mod.nuclear_stain_kind(names, ctx.cfg)
    risky = cap_mod.histone_denominator_markers(
        {str(m): ctx.space_for(str(m)) for m in full.markers}, kind)
    if risky:
        logger.warning("qc: %s are divided by the nuclear channel (%s), a histone label, not a "
                       "DNA stain: its level varies from cell to cell, and if it is a reporter "
                       "this divides by the signal. Prefer space: raw (or set nuclear_stain: dna "
                       "if the name misleads).", ", ".join(risky), "/".join(names))

    # Nuclei clipped at the detector ceiling: listed here, left out of everything
    # after (Context.data). Outlier fields too, when qc.drop_outlier_fields.
    ctx.saturated_mask(full)
    hits = ctx._cache.get("saturated_hits")
    if hits is not None and len(hits):
        ctx.write_table(hits, "qc_saturated_nuclei.csv")
    elif ctx.has_table("qc_saturated_nuclei.csv"):
        ctx.table("qc_saturated_nuclei.csv").unlink()
    ctx._cache.pop("data", None)


def _nuclear_names(ctx) -> list[str]:
    """Names of the nuclear (segmentation) channel, from the extraction channel QC."""
    path = ctx.out / assemble_mod.CHANNEL_QC_NAME
    if not path.exists():
        return []
    qc = pd.read_csv(path)
    names = set()
    for ds in ctx.datasets:
        dapi = int((ds.cfg.get("channels") or {}).get("dapi_index", 0))
        sel = (qc["dataset"].astype(str) == ds.name) & (qc["channel"] == dapi)
        names |= set(qc.loc[sel, "marker"].astype(str))
    return sorted(names)


def _caption(ctx) -> str:
    """What each marker's value is normalized to, for every marker figure."""
    from .signaling import reference as sig_ref
    markers = [str(m) for m in ctx.markers]
    spaces = {m: ctx.space_for(m) for m in markers}
    ref = ctx.read_table(sig_ref.TABLE) if ctx.has_table(sig_ref.TABLE) else None
    signaling = {m: lab for m, (lab, _) in sig_ref.unit_labels(ref).items()}
    names = _nuclear_names(ctx)
    return cap_mod.normalization_caption(
        spaces, {}, signaling, nuclear_names=names,
        kind=cap_mod.nuclear_stain_kind(names, ctx.cfg),
        note=(ctx.cfg.get("visualization") or {}).get("caption"), scaled=False)


# ---------------------------------------------------------------------------
# 3. Thresholds
# ---------------------------------------------------------------------------

@register("thresholds", after=["assemble", "qc"],
          config_keys=["thresholds", "conditions", "markers", "default_space", "signaling",
                       "qc"],
          description="Positivity thresholds per scope × marker, with evidence")
def step_thresholds(ctx) -> None:
    controls = ctx.controls
    table = thr_mod.compute_thresholds(ctx.data, ctx.cfg, controls)
    ctx.write_table(table, "thresholds.csv")
    ctx.set("thresholds", table)

    # Which wells fed each marker's gate. A per-marker or pooled control is only
    # trustworthy if it is auditable, so this table is always written.
    assignment = controls.summary(ctx.data.markers)
    ctx.write_table(assignment, "control_assignment.csv")
    if controls.is_marker_specific():
        for _, row in assignment.iterrows():
            logger.info("control %-10s ← %s (%d cells, %s wells)", row["marker"],
                        row["control_source"], row["n_control_cells"],
                        row.get("n_control_wells", "?"))

    _signaling_reference(ctx, table)
    plots.threshold_qc(ctx.gated, table, ctx.data.cells, controls, ctx.cfg,
                       ctx.figure("qc", "threshold_evidence.png"))

    graded = table[table["verdict"] != thr_mod.OK]
    if len(graded):
        logger.warning("%d (scope × marker) combination(s) have NO gate and are "
                       "reported as graded-only: %s", len(graded),
                       ", ".join(sorted(set(graded["marker"].astype(str)))))


def _signaling_reference(ctx, table: pd.DataFrame) -> None:
    """Reference level of each signaling marker per scope (signaling_reference.csv):
    every later table and figure reads signaling markers relative to it."""
    from .signaling import reference as sig_ref
    matches = sig_mod.classify_all([str(m) for m in ctx.data.markers], ctx.cfg)
    ctx._cache.pop("gated", None)
    if not matches:
        if ctx.has_table(sig_ref.TABLE):
            ctx.table(sig_ref.TABLE).unlink()
        return
    long = thr_mod.apply_thresholds(ctx.data, table, ctx.cfg)
    long["field_id"] = ctx.data.cell_column("field_id", long["cell_uid"].to_numpy())
    ref = sig_ref.reference_table(long, ctx.data.cells, ctx.cfg, matches)
    if len(ref):
        ctx.write_table(ref, sig_ref.TABLE)
    elif ctx.has_table(sig_ref.TABLE):
        ctx.table(sig_ref.TABLE).unlink()
    ctx._cache.pop("gated", None)


# ---------------------------------------------------------------------------
# 4. Rollups
# ---------------------------------------------------------------------------

@register("aggregate", after=["thresholds"], config_keys=["conditions", "stats", "qc"],
          description="Cell → field → well → condition rollups; replicate = field, well or "
                      "experiment")
def step_aggregate(ctx) -> None:
    fields = agg_mod.field_table(ctx.data, ctx.gated, ctx.cfg)
    exp_col = (ctx.cfg.get("stats") or {}).get("experiment_column")
    fields = agg_mod.add_well_key(_with_plate_columns(ctx, fields, [exp_col]))
    ctx.write_table(fields, "fields.csv")
    ctx.set("fields", fields)

    unit, note = agg_mod.choose_unit(fields, ctx.cfg)
    logger.info("aggregate: %s", note)
    with open(ctx.table("replicate_unit.json"), "w") as f:
        json.dump({"unit": unit, "note": note}, f, indent=2)
    ctx.set("unit", unit)

    wells = agg_mod.well_table(fields, ctx.cfg)
    ctx.write_table(wells, "wells.csv")
    ctx.set("wells", wells)

    conditions = agg_mod.condition_table(fields, ctx.cfg, unit)
    ctx.write_table(conditions, "conditions.csv")
    ctx.set("conditions", conditions)

    uniformity = qc_mod.well_uniformity(fields, ctx.cfg)
    ctx.write_table(uniformity, "well_uniformity.csv")

    _write_comparisons(ctx, fields, unit)


def _replicate_unit(ctx) -> str:
    try:
        with open(ctx.table("replicate_unit.json")) as f:
            return json.load(f)["unit"]
    except (OSError, ValueError, KeyError):
        return "field"


def _write_comparisons(ctx, fields: pd.DataFrame, unit: str = "field") -> None:
    """Compare every condition against each marker's **own** negative control.

    One global reference is wrong as soon as controls are marker-specific: the
    pluripotent well is the right baseline for BRA and the *positive* control for
    SOX2, so comparing SOX2 against it reports the effect backwards. An explicit
    ``stats.reference_condition`` still overrides everything.

    Grouping is by ``scope_key`` rather than by the union of scope columns, so a
    plate-scoped marker is compared against its reference pooled the same way its
    gate was: with the union, its two panels would be compared separately and a
    reference present in only one of them would drop the other entirely.
    """
    explicit = (ctx.cfg.get("stats") or {}).get("reference_condition")
    key_cols = tuple(c for c in ("scope_key", "marker") if c in fields.columns)

    per_value: dict[str, list[pd.DataFrame]] = {"positive_fraction": [], "median": []}
    chosen: dict[str, str] = {}
    for marker in ctx.markers:
        sub = fields[fields["marker"] == marker]
        if sub.empty:
            continue
        reference = explicit
        if reference is None:
            candidates = ctx.control_conditions_for(marker) & set(
                sub["condition"].astype(str).unique())
            reference = sorted(candidates)[0] if candidates else None
        if reference is None:
            logger.info("aggregate: %s has no control condition; no comparisons.", marker)
            continue
        chosen[marker] = reference
        for value_col in per_value:
            comp = stats_mod.compare_to_reference(
                sub, value_col, reference, extra_keys=key_cols, unit=unit, cfg=ctx.cfg,
                n_perm=int((ctx.cfg.get("stats") or {}).get("n_permutations", 5000)),
                seed=int((ctx.cfg.get("stats") or {}).get("random_state", 0)))
            if len(comp):
                comp.insert(0, "reference_condition", reference)
                per_value[value_col].append(comp)

    for value_col, blocks in per_value.items():
        if blocks:
            ctx.write_table(pd.concat(blocks, ignore_index=True),
                            f"comparisons_{value_col}.csv")
    if chosen:
        distinct = sorted(set(chosen.values()))
        if len(distinct) == 1:
            logger.info("aggregate: comparisons computed against reference %r", distinct[0])
        else:
            logger.info("aggregate: per-marker references; %s",
                        ", ".join(f"{m}→{r}" for m, r in sorted(chosen.items())))


# ---------------------------------------------------------------------------
# 5. Graded gating + co-expression
# ---------------------------------------------------------------------------

@register("gating", after=["thresholds"], config_keys=["gating", "conditions"],
          description="Tertile bands within each positive population + co-expression matrices")
def step_gating(ctx) -> None:
    edges = gating_mod.positive_edges(ctx.gated, ctx.controls, ctx.cfg)
    ctx.write_table(edges, "band_edges.csv")
    ctx.set("edges", edges)

    banded = gating_mod.assign_bands(ctx.gated, edges, ctx.cfg)
    ctx.set("banded", banded)

    composition = gating_mod.band_composition(ctx.data, banded, ctx.cfg)
    ctx.write_table(composition, "band_composition.csv")

    pairs = gating_mod.marker_pairs(ctx.data, banded, ctx.cfg)
    coex = gating_mod.coexpression(ctx.data, banded, ctx.cfg, pairs)
    if len(coex):
        ctx.write_table(coex, "coexpression.csv")


# ---------------------------------------------------------------------------
# 6. Figures
# ---------------------------------------------------------------------------

@register("distributions", after=["thresholds", "aggregate"],
          config_keys=["visualization", "conditions"],
          description="Per-marker violin (with field medians) and ECDF figures")
@_noted
def step_distributions(ctx) -> None:
    ctx.reset_figures("distributions")
    gated, fields = ctx.gated, ctx.fields
    for marker, scope, suffix in ctx.marker_scopes():
        thr = ctx.threshold_for(marker, scope)
        space = ctx.space_for(marker)
        sub_gated = ctx.scoped(gated, scope)
        sub_fields = ctx.scoped(fields, scope)
        # This marker's own controls: a pluripotent well is grayed out as the
        # control on the BRA panel and drawn normally where it is SOX2's positive.
        marker_controls = ctx.control_conditions_for(marker)
        # One pass when no facet is configured, one figure per facet value when
        # there is. The gate is a property of the marker, not of the facet, so the
        # same threshold line is drawn on every panel.
        for facet, (part_gated, part_fields), relabel in ctx.facets(sub_gated, sub_fields):
            tag = f"{_safe(marker)}{_safe(suffix)}{_safe(ctx.facet_suffix(facet))}"
            title = ctx.facet_title(facet)
            controls = {relabel.get(c, c) for c in marker_controls}
            plots.violin_by_condition(part_gated, part_fields, marker, space, thr, ctx.cfg,
                                      ctx.figure("distributions", f"violin_{tag}.png"),
                                      controls, facet=title)
            sub = part_gated[part_gated["marker"] == marker].copy()
            sub["is_control"] = sub["condition"].astype(str).isin(controls)
            plots.ecdf_by_condition(sub, marker, space, thr, ctx.cfg,
                                    ctx.figure("distributions", f"ecdf_{tag}.png"),
                                    facet=title)


@register("composition", after=["gating", "aggregate"],
          config_keys=["visualization", "conditions", "gating"],
          description="Positive fractions, graded band composition, co-expression grids")
@_noted
def step_composition(ctx) -> None:
    ctx.reset_figures("composition")
    conditions, fields = ctx.conditions, ctx.fields

    for marker, scope, suffix in ctx.marker_scopes():
        marker_controls = ctx.control_conditions_for(marker)
        for facet, (part_cond, part_fields), relabel in ctx.facets(
                ctx.scoped(conditions, scope), ctx.scoped(fields, scope)):
            tag = f"{_safe(marker)}{_safe(suffix)}{_safe(ctx.facet_suffix(facet))}"
            plots.positive_fraction_bars(
                part_cond, part_fields, marker, ctx.cfg,
                ctx.figure("composition", f"positive_fraction_{tag}.png"),
                {relabel.get(c, c) for c in marker_controls},
                facet=ctx.facet_title(facet))

    if ctx.has_table("band_composition.csv"):
        bands = ctx.read_table("band_composition.csv")
        for marker, scope, suffix in ctx.marker_scopes():
            plots.band_composition_bars(
                ctx.scoped(bands, scope), marker, ctx.cfg,
                ctx.figure("composition", f"bands_{_safe(marker)}{_safe(suffix)}.png"),
                ctx.control_conditions_for(marker))

    if ctx.has_table("coexpression.csv"):
        # A pair spans two markers whose controls may differ, so the marker-agnostic
        # set is the honest one to shade with here.
        controls = ctx.control_conditions
        coex = ctx.read_table("coexpression.csv")
        scope_cols = [c for c in ctx.scope_columns if c in coex.columns]
        group_cols = scope_cols + ["marker_a", "marker_b"]
        multi = coex[scope_cols].drop_duplicates().shape[0] > 1 if scope_cols else False
        for keys, part in coex.groupby(group_cols, observed=True):
            keys = keys if isinstance(keys, tuple) else (keys,)
            info = dict(zip(group_cols, keys))
            marker_a, marker_b = info["marker_a"], info["marker_b"]
            suffix = ("_" + "-".join(str(info[c]) for c in scope_cols)) if multi else ""
            tag = f"{_safe(marker_a)}_vs_{_safe(marker_b)}{_safe(suffix)}"
            plots.coexpression_grid(part, marker_a, marker_b, ctx.cfg,
                                    ctx.figure("composition", f"coexpression_{tag}.png"),
                                    controls)
            plots.quadrant_bars(part, marker_a, marker_b, ctx.cfg,
                                ctx.figure("composition", f"quadrants_{tag}.png"), controls)


@register("joint", after=["gating"], config_keys=["visualization", "gating", "conditions"],
          description="Two-marker density and scatter panels with gates (scatter colored by a third marker), and the correlation matrix")
@_noted
def step_joint(ctx) -> None:
    ctx.reset_figures("joint")
    gated = ctx.gated
    edges = ctx.edges if ctx.has_table("band_edges.csv") else pd.DataFrame()
    scope_cols = ctx.scope_columns
    scopes = (gated[scope_cols].drop_duplicates().to_dict("records") if scope_cols
              else [{}])
    multi = len(scopes) > 1

    # One panel's markers are only comparable within that panel, so densities and
    # correlations are computed per plate × panel.
    for scope in scopes:
        sub = ctx.scoped(gated, scope)
        if sub.empty:
            continue
        suffix = ("_" + "-".join(str(v) for v in scope.values())) if multi else ""
        wide = sub.pivot_table(index="cell_uid", columns="marker", values="value",
                               observed=True, aggfunc="first")
        conditions = sub.drop_duplicates("cell_uid").set_index("cell_uid")["condition"]
        scoped_edges = ctx.scoped(edges, scope) if len(edges) else edges

        for marker_a, marker_b in gating_mod.marker_pairs(ctx.data, sub, ctx.cfg):
            if marker_a not in wide.columns or marker_b not in wide.columns:
                continue
            gates = {
                "thr_a": ctx.threshold_for(marker_a, scope),
                "thr_b": ctx.threshold_for(marker_b, scope),
                "edges_a": _edges_for(scoped_edges, marker_a),
                "edges_b": _edges_for(scoped_edges, marker_b),
            }
            tag = f"{_safe(marker_a)}_vs_{_safe(marker_b)}{_safe(suffix)}"
            plots.pair_density(wide, marker_a, marker_b, conditions, gates, ctx.cfg,
                               ctx.figure("joint", f"density_{tag}.png"))

            # Scatter companions: one per remaining marker of the panel, each
            # coloring the cells by that marker; a two-marker panel gets one
            # plain scatter.
            others = [m for m in wide.columns if m not in (marker_a, marker_b)]
            if not others:
                plots.pair_scatter_colored(wide, marker_a, marker_b, conditions, gates,
                                           ctx.cfg, ctx.figure("joint", f"scatter_{tag}.png"))
            for marker_c in others:
                tag_c = (f"{_safe(marker_a)}_vs_{_safe(marker_b)}_by_{_safe(marker_c)}"
                         f"{_safe(suffix)}")
                plots.pair_scatter_colored(
                    wide, marker_a, marker_b, conditions, gates, ctx.cfg,
                    ctx.figure("joint", f"scatter_{tag_c}.png"),
                    color_by=marker_c, thr_c=ctx.threshold_for(marker_c, scope))

        if wide.shape[1] >= 2:
            plots.marker_correlation(
                wide.dropna(how="any"), ctx.cfg,
                ctx.figure("joint", f"marker_correlation{_safe(suffix)}.png"),
                title=f"Marker correlation (per cell){suffix and ': ' + suffix.lstrip('_')}")


@register("plate", after=["aggregate"], config_keys=["visualization", "conditions"],
          description="Plate-position heatmaps and per-well field uniformity")
@_noted
def step_plate(ctx) -> None:
    ctx.reset_figures("plate")
    wells, fields = ctx.wells, ctx.fields
    # Sequentially-numbered wells carry no plate coordinate, so the heatmaps have
    # no geometry to draw. Say so once: a silently missing figure looks like a bug.
    if "row" not in wells.columns or "col" not in wells.columns:
        logger.info("plate: no row/col in the manifest (wells are numbered by "
                    "acquisition order): layout heatmaps skipped, per-well "
                    "uniformity still drawn.")
    for marker, scope, suffix in ctx.marker_scopes():
        tag = f"{_safe(marker)}{_safe(suffix)}"
        controls = ctx.control_conditions_for(marker)
        sub_wells = ctx.scoped(wells, scope)
        sub_fields = ctx.scoped(fields, scope)
        plots.plate_heatmap(sub_wells, "positive_fraction", marker, ctx.cfg,
                            ctx.figure("plate", f"plate_positive_{tag}.png"),
                            title=f"{marker}-positive fraction by plate position", percent=True)
        plots.plate_heatmap(sub_wells, "median", marker, ctx.cfg,
                            ctx.figure("plate", f"plate_median_{tag}.png"),
                            title=f"{marker} median intensity by plate position")
        plots.field_uniformity(sub_fields, marker, "positive_fraction", ctx.cfg,
                               ctx.figure("plate", f"uniformity_{tag}.png"), controls)


@register("trends", after=["aggregate"],
          config_keys=["trends", "conditions", "visualization", "stats"],
          description="Time/dose courses and factorial marginal effects")
@_noted
def step_trends(ctx) -> None:
    ctx.reset_figures("trends")
    trends_cfg = ctx.cfg.get("trends") or {}
    x_col = trends_cfg.get("x")
    series_col = trends_cfg.get("series")
    facet_col = trends_cfg.get("facet")

    # `field_table` keeps only the condition keys, but `trends.x` is explicitly
    # meant to be ANY numeric plate-map column: dose, timepoint, a coordinate.
    # Without pulling those back from the cell table the timecourse finds no such
    # column and draws nothing at all, silently.
    fields = _with_plate_columns(ctx, ctx.fields, [x_col, series_col, facet_col])

    if x_col and x_col not in fields.columns:
        logger.warning("trends: x=%r is not a plate-map column (have: %s); "
                       "no time/dose courses drawn.", x_col,
                       ", ".join(sorted(set(ctx.data.cells.columns) - set(ID_COLUMNS))))
    elif x_col:
        for marker, scope, suffix in ctx.marker_scopes():
            tag = f"{_safe(marker)}{_safe(suffix)}"
            sub = ctx.scoped(fields, scope)
            plots.timecourse(sub, marker, "positive_fraction", x_col, series_col,
                             ctx.cfg, ctx.figure("trends", f"timecourse_positive_{tag}.png"),
                             facet_col, percent=True)
            plots.timecourse(sub, marker, "median", x_col, series_col, ctx.cfg,
                             ctx.figure("trends", f"timecourse_median_{tag}.png"),
                             facet_col)
    else:
        logger.info("trends: no 'trends.x' configured; skipping time/dose courses.")

    factors = trends_cfg.get("factors")
    if factors is None:
        factors = stats_mod.detect_factors(ctx.data.cells)
        # Never treat a control declaration as an experimental factor: neither the
        # class column nor the reserved per-marker control keys.
        control_col = (ctx.cfg.get("conditions") or {}).get("control_column", "class")
        reserved = {control_col, controls_mod.CONTROL_FOR, controls_mod.POSITIVE_FOR}
        factors = [f for f in factors if f not in reserved]
    if factors:
        logger.info("trends: factorial analysis over %s", factors)
        merged = _with_plate_columns(ctx, fields, factors)

        # Scoped by the marker's own gate scope: the same wells imaged as MIPs and
        # as z-stacks are two views of one experiment, not twice the replicates.
        # Pooling them would inflate n and average two differently-gated
        # measurements. ``scope_key`` carries that per marker.
        keys = tuple(c for c in ("scope_key", "marker") if c in merged.columns)
        for value_col in ("positive_fraction", "median"):
            effects = stats_mod.factor_effects(
                merged, value_col, factors, extra_keys=keys, unit=_replicate_unit(ctx),
                cfg=ctx.cfg,
                n_perm=int((ctx.cfg.get("stats") or {}).get("n_permutations", 5000)))
            if not len(effects):
                continue
            ctx.write_table(effects, f"factor_effects_{value_col}.csv")
            for marker, scope, suffix in ctx.marker_scopes():
                plots.factor_effects_plot(
                    ctx.scoped(effects, scope), marker, value_col, ctx.cfg,
                    ctx.figure("trends",
                               f"factors_{value_col}_{_safe(marker)}{_safe(suffix)}.png"))
        if len(factors) >= 2:
            a, b = factors[0], factors[1]
            for marker, scope, suffix in ctx.marker_scopes():
                plots.interaction_heatmap(
                    ctx.scoped(merged, scope), marker, "positive_fraction", a, b, ctx.cfg,
                    ctx.figure("trends", f"interaction_{_safe(a)}_{_safe(b)}_"
                                         f"{_safe(marker)}{_safe(suffix)}.png"),
                    percent=True)


# ---------------------------------------------------------------------------
# 6b. Effects, trends with uncertainty, variance components, power
# ---------------------------------------------------------------------------

EFFECTS_FIG = "effects"
EFFECT_TABLES = ("effects_variance.csv", "effects_power.csv", "effects_trends.csv",
                 "effects_vs_reference.csv")


def _remove_tables(ctx, *names) -> None:
    for name in names:
        if ctx.has_table(name):
            ctx.table(name).unlink()


def _numeric(frame: pd.DataFrame, col: str):
    vals = pd.to_numeric(frame[col], errors="coerce")
    present = frame[col].notna()
    return vals if present.any() and vals[present].notna().all() else None


def _key_frame(frame: pd.DataFrame, cols: list[str]) -> list[pd.Series]:
    return [frame[c].astype(object).where(frame[c].notna(), "").astype(str).rename(c)
            for c in cols]


def _units_at(part, by, unit, cfg, need: int = 2):
    """Replicate values at *unit*, stepping down until every condition has ≥ *need*."""
    start = agg_mod.UNIT_LADDER.index(unit)
    for level in agg_mod.UNIT_LADDER[start:]:
        if level == "experiment" and not agg_mod.experiment_column(part, cfg):
            continue
        vals = agg_mod.replicate_values(part, by, level, cfg, "value")
        if len(vals) and vals.groupby("condition").size().min() >= need:
            return vals, level
    return agg_mod.replicate_values(part, by, "field", cfg, "value"), "field"


def _field_scalars(ctx) -> pd.DataFrame:
    """Every per-field number, long: positive fraction and median per marker × scope."""
    fields = ctx.fields
    exp_col = (ctx.cfg.get("stats") or {}).get("experiment_column")
    keys = ctx.condition_keys
    fields = _with_plate_columns(ctx, fields, keys + [exp_col])
    fields = agg_mod.add_well_key(fields)
    id_cols = [c for c in dict.fromkeys(["field_id", "well_key", "well", "condition", "marker",
                                         "scope_key"] + keys + ([exp_col] if exp_col else []))
               if c in fields.columns]
    long = fields.melt(id_vars=id_cols, value_vars=["positive_fraction", "median"],
                       var_name="metric", value_name="value")
    long = long[np.isfinite(long["value"].to_numpy(np.float64))]
    multi = (long.groupby("marker", observed=True)["scope_key"].nunique() > 1
             if "scope_key" in long else {})
    subject = long["marker"].astype(str)
    if "scope_key" in long and len(long):
        subject = np.where(long["marker"].map(multi).fillna(False),
                           subject + " [" + long["scope_key"].astype(str) + "]", subject)
    long = long.assign(family="marker", subject=subject, readout="field")
    long["metric_id"] = long["subject"].astype(str) + "|" + long["metric"].astype(str)
    return long


@register("effects", after=["aggregate"],
          config_keys=["effects", "stats", "conditions"],
          description="Per-field numbers → variance components (field / well / experiment), "
                      "effects vs a reference with CI, dose/time trends with CI, power")
def step_effects(ctx) -> None:
    ctx.reset_figures(EFFECTS_FIG)
    prm = ef_mod.settings(ctx.cfg)
    if not bool(prm["enabled"]):
        _remove_tables(ctx, *EFFECT_TABLES)
        return
    sc = _field_scalars(ctx)
    if sc.empty:
        _remove_tables(ctx, *EFFECT_TABLES)
        return
    unit = _replicate_unit(ctx)
    exp_col = agg_mod.experiment_column(sc, ctx.cfg)
    cond_keys = ctx.condition_keys
    rng = np.random.default_rng(int(prm["random_state"]))
    ident = ["metric_id", "family", "subject", "readout", "metric"]

    var_rows, power_rows = [], []
    targets = prm["power"].get("targets") or []
    for key, part in sc.groupby(ident, observed=True, dropna=False):
        vc = ef_mod.variance_components(part, ["condition"], exp_col)
        per_well = part.groupby("well_key").size()
        wells = part.groupby(["condition"] + ([exp_col] if exp_col else []))["well_key"].nunique()
        row = dict(zip(ident, key), **vc, wells_per_experiment=float(wells.median()),
                   fields_per_well=float(per_well.median()))
        var_rows.append(row)
        ref_mean = float(part["value"].mean())
        for pr in ef_mod.power_table(vc, prm, row["wells_per_experiment"],
                                     row["fields_per_well"], ref_mean):
            out = dict(zip(ident, key), **pr)
            for t in ef_mod.match_targets(key[0], key[4], targets):
                delta = float(t["delta"]) if t.get("delta") is not None else \
                    float(t.get("delta_rel", np.nan)) * abs(ref_mean)
                var = ef_mod.group_mean_var(pr["sd_experiment_used"], vc["sd_well"],
                                            vc["sd_field"], pr["experiments"],
                                            pr["wells_per_experiment"], pr["fields_per_well"])
                out.update(target_delta=delta, power_at_target=ef_mod.power_for(
                    delta, var, pr["experiments"], float(prm["power"]["alpha"])))
            power_rows.append(out)
    ctx.write_table(pd.DataFrame(var_rows), "effects_variance.csv")
    ctx.write_table(pd.DataFrame(power_rows), "effects_power.csv")

    if prm["trend_columns"] == "auto":
        trend_cols = [c for c in cond_keys if c in sc and _numeric(sc, c) is not None
                      and sc[c].nunique() >= int(prm["min_levels"])]
    else:
        trend_cols = [c for c in (prm["trend_columns"] or []) if c in sc]
    trend_rows = []
    for col in trend_cols:
        others = [k for k in cond_keys if k != col and k in sc]
        for key, part in sc.groupby(_key_frame(sc, ident + others), observed=True):
            if part[col].nunique() < int(prm["min_levels"]):
                continue
            vals, used = _units_at(part, ["condition", col], unit, ctx.cfg)
            x = pd.to_numeric(vals[col], errors="coerce").to_numpy(float)
            scale = ef_mod.choose_scale(np.unique(x[np.isfinite(x)]), prm["scale"])
            res = ef_mod.trend(x, vals["value"].to_numpy(float), scale, int(prm["n_boot"]),
                               int(prm["n_perm"]), rng,
                               exponential=col in (prm["exponential_columns"] or []))
            group = " | ".join(f"{k}={v}" for k, v in zip(others, key[len(ident):])
                               if v != "") or "all"
            trend_rows.append(dict(zip(ident, key[:len(ident)]), column=col, group=group,
                                   unit=used, pseudo_replicated=used != unit, **res))
    trends = pd.DataFrame(trend_rows)
    if len(trends):
        ctx.write_table(trends, "effects_trends.csv")
    else:
        _remove_tables(ctx, "effects_trends.csv")

    eff_rows = []
    for key, part in sc.groupby(ident, observed=True, dropna=False):
        groups = [("all", part)]
        ref_of = {}
        if prm["reference"] == "auto" and trend_cols:
            col = trend_cols[0]
            others = [k for k in cond_keys if k != col and k in part]
            groups = list(part.groupby(_key_frame(part, others))) if others else [("all", part)]
            for gkey, g in groups:
                xs = pd.to_numeric(g[col], errors="coerce")
                ref_of[str(gkey)] = str(g.loc[xs.idxmin(), "condition"]) if xs.notna().any() \
                    else None
        elif prm["reference"] not in (None, "auto"):
            ref_of = {"all": str(prm["reference"])}
        for gkey, g in groups:
            ref = ref_of.get(str(gkey))
            if ref is None or ref not in set(g["condition"].astype(str)):
                continue
            for cond in sorted(set(g["condition"].astype(str)) - {ref},
                               key=plots.base.natural_key):
                pair = g[g["condition"].astype(str).isin([cond, ref])]
                vals, used = _units_at(pair, ["condition"], unit, ctx.cfg)
                a = vals.loc[vals["condition"].astype(str) == cond, "value"].to_numpy(float)
                r = vals.loc[vals["condition"].astype(str) == ref, "value"].to_numpy(float)
                eff_rows.append(dict(zip(ident, key), condition=cond, reference=ref, unit=used,
                                     pseudo_replicated=used != unit,
                                     **ef_mod.effect(a, r, int(prm["n_perm"]), rng)))
    effects = pd.DataFrame(eff_rows)
    if len(effects):
        ctx.write_table(effects, "effects_vs_reference.csv")
    else:
        _remove_tables(ctx, "effects_vs_reference.csv")

    values = set(prm.get("figure_values") or ["positive_fraction", "median"])
    pw_all = pd.DataFrame(power_rows)
    for subject, part in sc[sc["metric"].isin(values)].groupby("subject"):
        tag = _safe(subject)
        if len(trends):
            t = trends[trends["subject"] == subject]
            if len(t):
                plots.trend_panels(part, t, unit, exp_col,
                                   ctx.figure(EFFECTS_FIG, f"trend_{tag}.png"),
                                   title=f"{subject}: trend with "
                                         f"{', '.join(sorted(set(t['column'])))}")
        if len(effects):
            e = effects[effects["subject"] == subject]
            if len(e):
                plots.effect_forest(e, ctx.figure(EFFECTS_FIG, f"vs_reference_{tag}.png"),
                                    title=f"{subject}: difference from the reference condition")
        pw = pw_all[pw_all["subject"] == subject] if len(pw_all) else pw_all
        if len(pw):
            plots.power_curves(pw, ctx.figure(EFFECTS_FIG, f"power_{tag}.png"),
                               title=f"{subject}: minimum detectable difference "
                                     f"(α = {prm['power']['alpha']}, power = "
                                     f"{prm['power']['power']})")


# ---------------------------------------------------------------------------
# 7. Report
# ---------------------------------------------------------------------------

@register("report", after=["distributions", "composition", "joint", "plate", "trends",
                          "effects", "gating"],
          config_keys=["report", "skip_steps"],
          description="Single scrollable HTML contact sheet of every figure and table")
def step_report(ctx) -> None:
    from .report import build_report
    build_report(ctx)


# ---------------------------------------------------------------------------

def _with_plate_columns(ctx, fields: pd.DataFrame, columns) -> pd.DataFrame:
    """Attach plate-map *columns* to the field table, matched on ``field_id``.

    The field table is built from the condition keys only, so every other
    plate-map column (precisely what ``trends.x`` and the factorial
    analysis consume) has to be pulled back off the cell table. A field is
    wholly within one well, so ``first()`` is the whole value, not a summary.
    """
    out = fields
    for col in columns:
        if not col or col in out.columns or col not in ctx.data.cells.columns:
            continue
        if out is fields:
            out = fields.copy()
        mapping = ctx.data.cells.groupby("field_id", observed=True)[col].first()
        out[col] = out["field_id"].map(mapping)
    return out


def _safe(name: str) -> str:
    """Filename-safe marker name ('NKX2.5' → 'NKX2-5')."""
    return "".join(c if c.isalnum() or c in "-_" else "-" for c in str(name))


def _edges_for(edges: pd.DataFrame, marker: str) -> list[float]:
    if edges is None or not len(edges):
        return []
    sub = edges[edges["marker"] == marker]
    if sub.empty:
        return []
    return [float(sub["edge_lo"].median()), float(sub["edge_hi"].median())]


# ---------------------------------------------------------------------------
# Optional steps: what `downstream.skip_steps` may leave out
# ---------------------------------------------------------------------------
# The core (signaling → assemble → qc → thresholds → aggregate → report) always
# runs: every readout depends on it. These add readouts and figures on top. A
# skipped step's earlier outputs move to <analysis>/_skipped/<step>/ so the
# report never shows a stale figure.
OPTIONAL_STEPS = {
    # step: (figure folders, table name patterns)
    "gating": ([], ["band_edges.csv", "band_composition.csv", "coexpression.csv"]),
    "distributions": (["distributions"], []),
    "composition": (["composition"], []),
    "joint": (["joint"], []),
    "plate": (["plate"], []),
    "trends": (["trends"], ["factor_effects_*.csv"]),
    "effects": ([EFFECTS_FIG], ["effects_*.csv"]),
}
SKIP_DEPENDANTS: dict = {}


def resolve_skips(skip) -> set[str]:
    """Validate ``downstream.skip_steps``."""
    skip = {str(s).strip() for s in (skip or []) if str(s).strip()}
    bad = sorted(skip - set(OPTIONAL_STEPS))
    if bad:
        raise ValueError(f"downstream.skip_steps: {bad} cannot be skipped. Optional steps: "
                         f"{sorted(OPTIONAL_STEPS)} (the rest is the core of every readout).")
    return skip


def set_aside(out_dir: Path, step: str) -> None:
    """Move a skipped step's previous figures and tables to ``_skipped/<step>/``."""
    import shutil
    from .registry import checkpoint_path
    figs, tables = OPTIONAL_STEPS[step]
    dest = Path(out_dir) / "_skipped" / step
    moved = []
    for folder in figs:
        src = Path(out_dir) / "figures" / folder
        if src.exists():
            target = dest / "figures" / folder
            if target.exists():
                shutil.rmtree(target)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(src), str(target))
            moved.append(f"figures/{folder}")
    for pattern in tables:
        for src in sorted((Path(out_dir) / "tables").glob(pattern)):
            (dest / "tables").mkdir(parents=True, exist_ok=True)
            shutil.move(str(src), str(dest / "tables" / src.name))
            moved.append(f"tables/{src.name}")
    ck = checkpoint_path(out_dir, step)
    if ck.exists():
        ck.unlink()
    if moved:
        logger.info("[skip] %-14s skipped by downstream.skip_steps; its earlier outputs moved "
                    "to %s (%d files/folders)", step, dest, len(moved))


__all__ = ["OPTIONAL_STEPS", "resolve_skips", "set_aside", "step_signaling", "step_effects",
           "step_assemble", "step_qc", "step_thresholds", "step_aggregate",
           "step_gating", "step_distributions", "step_composition", "step_joint",
           "step_plate", "step_trends", "step_report"]
