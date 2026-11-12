"""
A synthetic confocal experiment with known answers, written as the image
pipeline would have left it (masks, channel TIFFs, metadata, per-nucleus joblib),
so ``check``, ``manifest`` and ``downstream`` run end to end without Cellpose.

Design: doses 0 / 1 / 10 (``dose_ng_ml``), 3 wells each, 2 fields per well,
predefined points (Seq order). Channels:

0 DAPI    nuclear stain, flat
1 SOX2    standard marker: the fraction of positive nuclei rises with dose
          (0 % → 45 % → 85 %); dose-0 wells are the negative control
2 b-CAT   signaling marker (nuclear readout): nuclear level 1× / 2× / 4× the
          dose-0 level, cytoplasm constant, on a background of 100
3 ECAD    present in the images, left out with ``exclude_markers``

One SOX2 nucleus in the first field sits at the 16-bit ceiling (65535).
"""

from __future__ import annotations

import json
from pathlib import Path

import joblib
import numpy as np
import tifffile
import yaml

PIXEL_UM = 0.65
DOSES = (0, 1, 10)
POS_FRAC = {0: 0.0, 1: 0.45, 10: 0.85}
BCAT = {0: 1.0, 1: 2.0, 10: 4.0}


def _field(rng, dose: int, saturate: bool = False, shape=(160, 160), r=5, step=16):
    ys, xs = np.meshgrid(np.arange(step, shape[0] - step + 1, step),
                         np.arange(step, shape[1] - step + 1, step), indexing="ij")
    centers = np.column_stack([ys.ravel(), xs.ravel()]) + rng.integers(-2, 3, (ys.size, 2))
    yy, xx = np.indices(shape)
    masks = np.zeros(shape, np.uint32)
    cyto = np.zeros(shape, bool)
    for i, (cy, cx) in enumerate(centers, 1):
        d2 = (yy - cy) ** 2 + (xx - cx) ** 2
        masks[d2 <= r * r] = i
        cyto |= (d2 > (r + 1) ** 2) & (d2 <= (r + 5) ** 2)
    labels = np.arange(1, len(centers) + 1)
    n = len(labels)
    pos = rng.random(n) < POS_FRAC[dose]
    sox2_level = np.where(pos, 3000.0, 400.0) * rng.normal(1, 0.05, n)
    if saturate:
        sox2_level[0] = 65535
    dapi = np.where(masks > 0, 2000, 100).astype(np.uint16)
    sox2 = np.full(shape, 100.0)
    bcat = np.full(shape, 100.0)
    ecad = np.full(shape, 150.0)
    bcat[cyto] = 600.0
    for i, lab in enumerate(labels):
        m = masks == lab
        sox2[m] = sox2_level[i]
        bcat[m] = 100.0 + 500.0 * BCAT[dose] * rng.normal(1, 0.05)
        ecad[m] = 800.0
    chans = [dapi, sox2, bcat, ecad]
    chans = [np.clip(c, 0, 65535).astype(np.uint16) for c in chans]
    cents = centers.astype(np.float32)
    return masks, chans, labels, cents


def _write_result(res: Path, masks, chans, labels, cents) -> None:
    from scipy.ndimage import median
    res.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(res / "nuclear_masks.tif", masks)
    for i, c in enumerate(chans):
        tifffile.imwrite(res / f"ch{i}.tif", c)
    (res / "image_metadata.json").write_text(json.dumps({
        "file_format": ".tif", "height_px": masks.shape[0], "width_px": masks.shape[1],
        "n_channels": 4, "n_z": 1, "is_zstack": False,
        "pixel_size_x_um": PIXEL_UM, "pixel_size_y_um": PIXEL_UM,
        "channel_names": ["DAPI", "SOX2", "b-CAT", "ECAD"]}))
    (res / "segmentation_metadata.json").write_text(json.dumps({"mode": "2d"}))
    raw = {ch: median(chans[ch].astype(np.float32), masks, labels).astype(np.float32)
           for ch in (1, 2, 3)}
    dapi = median(chans[0].astype(np.float32), masks, labels).astype(np.float32)
    joblib.dump({
        "n_cells": len(labels), "labels": labels.astype(np.int64), "mode": "2d",
        "is_zstack": False, "dapi_index": 0, "marker_channels": [1, 2, 3],
        "centroid_axes": "yx", "raw_intensities": raw,
        "dapi_normalized_intensities": {ch: v / dapi for ch, v in raw.items()},
        "dapi_per_cell": dapi, "mask_sizes": np.full(len(labels), 80.0),
        "centroids": cents,
        "global_stats": {ch: {"saturation_value": 65535.0} for ch in (0, 1, 2, 3)},
        "channel_names": {0: "DAPI", 1: "SOX2", 2: "b-CAT", 3: "ECAD"},
    }, res / "per_nuclei_intensities.joblib")
    for name in (".imgproc_output", ".loaded", ".segmented", ".extracted"):
        (res / name).touch()


def make_experiment(root: Path, wells_per_dose: int = 3, fields_per_well: int = 2,
                    downstream: dict | None = None, seed: int = 0) -> Path:
    """Write the experiment under *root*; returns the config path."""
    root = Path(root)
    data = root / "zmax"
    data.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    seq, well, wells = 0, 0, {}
    for dose in DOSES:
        for _ in range(wells_per_dose):
            well += 1
            wells[well] = {"dose_ng_ml": dose, **({"class": "control"} if dose == 0 else {})}
            for _ in range(fields_per_well):
                stem = f"Point{seq:04d}_Seq{seq:04d}"
                (data / f"{stem}.tif").touch()
                _write_result(data / stem, *_field(rng, dose, saturate=(seq == 0)))
                seq += 1
    pm = {"plate": "96well", "plate_id": "synthetic", "mode": "predefined",
          "points_per_well": fields_per_well,
          "stains": {"default": {0: "DAPI", 1: "SOX2", 2: "b-CAT", 3: "ECAD"}},
          "wells": wells}
    (data / "plate_map.yaml").write_text(yaml.safe_dump(pm, sort_keys=False))
    cfg = {
        "defaults": {"input": {"recursive": True, "file_types": [".tif"]},
                     "channels": {"dapi_index": 0, "marker_channels": [1, 2, 3]},
                     "segmentation": {"mode": "2d"}},
        "datasets": [{"name": "zmax", "directory": "./zmax"}],
        "downstream": {
            "results_dir": "auto",
            "exclude_markers": ["DAPI", "ECAD"],
            "conditions": {"condition_key": ["dose_ng_ml"], "control_column": "class",
                           "control_values": ["control"]},
            "default_space": "raw",
            "trends": {"x": "dose_ng_ml", "series": None},
            "signaling": {"reference_values": ["control"], "markers": {"b-CAT": "nc"},
                          "n_jobs": 1},
            "stats": {"n_permutations": 999},
            "effects": {"n_boot": 300, "n_perm": 300},
            "thresholds": {"min_control_cells": 50},
        },
    }
    if downstream:
        from pipeline_config import deep_merge
        cfg["downstream"] = deep_merge(cfg["downstream"], downstream)
    path = root / "config.yaml"
    path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    return path


__all__ = ["make_experiment", "DOSES", "POS_FRAC", "BCAT", "PIXEL_UM"]
