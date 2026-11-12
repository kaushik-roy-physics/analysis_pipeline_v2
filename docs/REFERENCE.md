# Confocal Image-Processing & Analysis Pipeline

A modular, format-agnostic pipeline for **multichannel confocal microscopy**. It
discovers image files under a directory, converts each raw acquisition into
analysis-ready arrays, segments nuclei with **Cellpose-SAM**, extracts
per-nucleus intensities, and then turns those measurements plus your plate map
into biology: positivity thresholds, graded co-expression, QC, publication
figures, and a single HTML contact sheet.

The pipeline is designed for shared lab use: point it at a folder, run it, and
every image gets a self-contained results directory; point the downstream stage
at the same config and every experiment gets a self-contained analysis directory.

> **This is the command-line reference manual.** Most people should start with
> the GUI described in the main [README](../README.md). New to plate maps and
> configs? Read [`TUTORIAL.md`](TUTORIAL.md) first.

---

## What it does

**Per image**: three independently runnable steps:

| Step        | Command           | Produces (per image)                                                                                |
| ----------- | ----------------- | --------------------------------------------------------------------------------------------------- |
| **load**    | `--steps load`    | `image_metadata.json`, `multichannel.tif`, `ch{n}.tif` (+ `*_MIP.tif` for z-stacks)                 |
| **segment** | `--steps segment` | `nuclear_masks.tif`, `segmentation_overlay.png`, `size_filter_qc.png`, `segmentation_metadata.json` |
| **extract** | `--steps extract` | `per_nuclei_intensities.joblib`                                                                     |

Each step writes a hidden checkpoint (`.loaded`, `.segmented`, `.extracted`) so
re-running the pipeline only redoes what is missing.

**Per experiment**: the `downstream` subcommand, eleven more checkpointed steps
from per-image joblibs to finished figures. See
[Downstream analysis](#downstream-analysis).

### Key features

- **Any file type**: `.nd2`, `.tif`/`.tiff`, `.ome.tif`, `.oir`, `.czi`,
  `.lif`, and anything else `bioio` can read. Backends are auto-selected per
  file; axis order is normalized automatically.
- **MIP vs z-stack auto-detection**: the loader inspects each image and writes
  the appropriate outputs. Z-stacks additionally get per-channel maximum
  intensity projections.
- **2-D or 3-D segmentation**: 2-D runs on the DAPI MIP (works for any input);
  3-D runs on the DAPI z-stack with anisotropy taken from the metadata.
- **Optional cytoplasm channel**: DAPI + cytoplasm can be passed together to
  Cellpose-SAM for joint segmentation.
- **Adaptive nucleus filtering**: removes debris/oversplit masks using a
  robust size (and optional DAPI) distribution estimator, with a QC figure.
- **Per-nucleus quantification**: raw and DAPI-normalized intensities plus
  nucleus **centroids** (2-D or 3-D) for locating saturated/outlier signal.
- **Measurement-health guardrails**: every channel is checked for clipping and
  for nuclei indistinguishable from their own background; every marker channel is
  compared field-to-field *within its own well*; and every gate is tested for
  how much of it rests on a single control field. All three are detectors:
  they warn and write columns, they never change a number. See
  [Guardrails](#guardrails-catching-a-broken-measurement).
- **Three acquisition schemes**: JOBS 'Generate Points' (`Well07`), 'Predefined
  Points' (`Seq0042`) and Spinning-Disk confocal (`WellC10`, alphanumeric), all
  selected per plate and all resolving to the same canonical well numbering.
- **Generated plate maps**: `init-plate-map` scans a folder and writes a
  commented starter `plate_map.yaml` pre-filled with the mode, wells and channels
  it actually found.
- **Per-marker gate scope and per-marker controls**: pool a shared antibody
  across stain panels, and declare that a well is the negative control for some
  markers and the positive control for others.

---

## Installation

```bash
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt
```

**GPU (recommended for segmentation):** install a CUDA build of PyTorch that
matches your driver *before* the rest, e.g.

```bash
pip install torch --index-url https://download.pytorch.org/whl/cu126
pip install -r requirements.txt
```

**Extra file formats:** the universal `bioio` backend needs a plugin per
format: e.g. `pip install bioio-czi` for `.czi`, `pip install bioio-bioformats`
(plus a Java runtime) for `.oir`. `.nd2` and TIFF variants work out of the box.

---

## Configuration: one file for the whole experiment

All behavior is driven by [`config.yaml`](../templates/config.yaml). `defaults:` holds what
is shared; `datasets:` lists the acquisition folders, and each entry is
**deep-merged over the defaults**, so a folder segmented in 2-D and one
segmented in 3-D, or one folder per antibody panel, live in the same file:

```yaml
defaults:
  input:
    exclude_dirs: [analysis, jpegs, pngs]   # folders that never hold input images
  channels:
    dapi_index: 0              # DAPI channel index (0-based)
    marker_channels: [1, 2, 3] # channels to quantify
  segmentation:
    mode: 2d                   # 2d | 3d | null  (null → 2-D on DAPI MIP)
    min_mask_size: 15
    adaptive_size_filter: true
  intensity:
    normalize_by_dapi: true

datasets:
  - {name: zmax_stain1, directory: ./zmax_stain1}
  - {name: zmax_stain2, directory: ./zmax_stain2}
  - {name: zstacks,     directory: ./zstacks, segmentation: {mode: 3d}}
```

A dataset states only what makes it different. Its `name` becomes the `dataset`
column in the manifest, which is what lets downstream pool or facet across 2-D
and 3-D or across panels. `--datasets <name> …` runs a subset.

> Configs with no `datasets:` list still work: the legacy top-level
> `input.directory` is treated as a single dataset, so existing configs and
> command lines are unchanged.

### One copy of the code, one config per experiment

**Every relative path in a config resolves against the directory holding that
config: never against the directory you happen to run from.** That single rule
lets the code live once, somewhere central, while each experiment keeps its own
config next to its own data:

```
KR/
├── project_structure/         # the code: one copy, nothing is written here
├── exp_033/
│   ├── config.yaml            # datasets: ./zmax_stain1, ./zstacks_stain1, …
│   ├── zmax_stain1/           # images, plate_map.yaml, per-image results
│   ├── zstacks_stain1/
│   ├── analysis/              # results_dir: ./analysis → lands here
│   └── pipeline.log           # the run log, beside the run it describes
└── exp_034/ …
```

```bash
PIPE=/path/to/KR/project_structure
EXP=/path/to/KR/exp_033

python "$PIPE/main.py" check      --config "$EXP/config.yaml"
python "$PIPE/main.py" process    --config "$EXP/config.yaml"
python "$PIPE/main.py" downstream --config "$EXP/config.yaml"
```

Starting the next experiment is `cp exp_033/config.yaml exp_034/` and editing
the dataset names. Keep the `directory:` values **relative** and the folder stays
portable: move or rename it and nothing breaks.

Two things that follow from the anchoring rule:

- Give sibling panel folders of the same physical plate an explicit, identical
  `plate_id:` in their plate maps. It defaults to the folder name, so
  `zmax_stain1` and `zmax_stain2` would otherwise count as different plates and
  cross-panel pooling (`scope: plate`) would silently do nothing. `check`
  reports which markers are poolable.
- `downstream --out /local/scratch/exp_033` redirects only the analysis tree, if
  writing figures back to a network share is slow.

---

## Usage

```bash
# Everything: load → segment → extract, for every dataset
python main.py

# Only convert/load images (e.g. to inspect channels first)
python main.py --steps load

# One dataset, segmentation + extraction (assumes load already ran)
python main.py process --datasets zmax_stain1 --steps segment extract

# Override the segmentation mode for the selected datasets
python main.py --seg-mode 3d

# Recompute a single step / everything
python main.py --force-step segment
python main.py --force

# Only some images: by Seq number (as in the filename) and/or by well.
# The same flags work for `check` and `downstream`; `input.select` in the
# config sets a default subset.
python main.py process --seq 0-9 15
python main.py process --wells 1-3          # or C1-C4 on the spinning disk
python main.py downstream --seq 0-9 15

# Author, validate, annotate, analyze
python main.py init-plate-map        # write starter plate maps, then edit them
python main.py check
python main.py manifest
python main.py downstream
```

Useful flags: `--config <path>`, `--log-level DEBUG`. A full run log is also
written to `pipeline.log` **in the config's own directory**, so each experiment
folder accumulates its own history no matter where the command was typed; if
that directory is not writable the run continues, logging to the working
directory instead. Failures on one image are logged with a full traceback and do
not stop the remaining images.

---

## Outputs

Each image `path/to/img.nd2` produces a sibling directory `path/to/img/`:

```
img/
├── image_metadata.json          # pixel sizes, channels, MIP-vs-stack, axes, ...
├── multichannel.tif             # CYX (MIP) or CZYX (z-stack)
├── ch0.tif, ch1.tif, ...        # per-channel images
├── ch0_MIP.tif, ...             # per-channel MIPs (z-stacks only)
├── multichannel_MIP.tif         # (z-stacks only)
├── nuclear_masks.tif            # label image (2-D, or z-stack for 3-D mode)
├── segmentation_overlay.png     # nuclei outlined on DAPI (+ _montage for 3-D)
├── size_filter_qc.png           # adaptive filter diagnostic
├── segmentation_metadata.json   # parameters + nucleus counts before/after
└── per_nuclei_intensities.joblib
```

### Reading the intensity results

```python
import joblib
r = joblib.load("path/to/img/per_nuclei_intensities.joblib")

r["n_cells"]                       # number of nuclei
r["labels"]                        # mask label id per nucleus
r["centroids"]                     # (n, 2) [y,x] or (n, 3) [z,y,x]
r["centroid_axes"]                 # "yx" or "zyx"
r["mask_sizes"]                    # pixels (2-D) / voxels (3-D) per nucleus
r["dapi_per_cell"]                 # per-nucleus DAPI level
r["raw_intensities"][2]            # per-nucleus median of channel 2
r["dapi_normalized_intensities"][2]  # channel 2 / DAPI per nucleus
r["channel_names"]                 # {index: name}
r["global_stats"][2]               # mean/median/percentiles for channel 2
r["global_stats"][2]["saturated_fraction"]      # pixels at the detector ceiling
r["global_stats"][2]["signal_to_background"]    # nuclear median ÷ image background
```

Arrays are aligned row-for-row by `labels`, so you can build a per-nucleus table
(intensities + centroids + sizes) directly for your downstream analysis.

---

## Organizing images by experimental condition

The image-processing steps above are deliberately **condition-agnostic**: they
never need to know which well or treatment an image came from. A separate,
fully decoupled layer (`main.py manifest`) maps each image to its **plate well**
and the **experimental condition** that well represents, and writes a single
`manifest.parquet` / `manifest.csv` you join to the per-image results downstream.

### 1. Author a `plate_map.yaml` per acquisition folder

**Start by generating one**, rather than writing it from scratch:

```bash
python main.py init-plate-map                    # one per dataset folder
python main.py init-plate-map --data-dir ./zmax_stain1 --plate-id exp033
python main.py init-plate-map --mode spinning_disk --force
```

This scans the folder and writes a heavily commented template pre-filled with
what it found: the acquisition mode detected from the filenames, the wells (or
Seq indices) actually present, the images-per-well count, and the channel names
from `image_metadata.json` when `process --steps load` has already run. Lines you
must edit are marked `<-- EDIT`. An existing plate map is never overwritten
without `--force`: it is hand-authored knowledge.

Then edit it. The schema in full:

```yaml
plate: 18well            # 8well | 18well | 96well  (or geometry: {rows, cols})
plate_id: exp033         # optional; defaults to the folder name (see the note below)
mode: predefined         # predefined | generate | spinning_disk

# Predefined/Seq only: map the 0-based global point index to a well:
points_per_well: 5       # constant N, wells filled in order from well 1
# seq_ranges:            #   …or explicit ranges for variable counts (inclusive)
#   1: [0, 4]
#   3: [10, 14]

stains:                  # channel index → marker. 'default' unless a well sets 'stain'.
  default: {0: DAPI, 1: GM130, 2: pS6, 3: LC3}
  setB:    {0: DAPI, 1: Tom20, 2: pS6, 3: LC3}

wells:                   # well → condition; range keys like "2-4" are expanded
  1:   {treatment: Control, dose_uM: 0, day: 1, class: control}
  2-4: {treatment: DrugX,   dose_uM: 10, day: 1, stain: setB}
```

Wells use the **NIKON Elements** row-major numbering you already see in JOBS:
8-well `01-04 / 05-08`, 18-well `01-06 / 07-12 / 13-18`, 96-well `01-12 / …`.

* **`generate` mode**: `Well07_*.nd2` carries the well number directly. Using
  only some wells needs nothing special: you just have no files for the rest.

* **`predefined` mode**: `*_Seq0042.nd2` carries the global point index, which
  is *acquisition order*. The filename says nothing about **where on the plate**
  the point came from, so wells here are numbered **sequentially**: well 1 is the
  first well you imaged, well 2 the second. Using 13 wells of an 18-well plate,
  in any order, is still wells 1-13: you never have to look up physical well
  numbers. Say how the points divide up with one of:
  
  * bare `points_per_well: N`; equal counts per well, in acquisition order. The
    usual case.
  * `points_per_well: {1: 5, 2: 8}` or a parallel list; variable counts.
  * `seq_ranges:`: explicit `well: [start, end]`; overrides everything.
  
  Because there is no plate position, **`row`/`col` are not written to the
  manifest** and the plate-layout heatmaps skip themselves (the step logs why).
  If your predefined keys really *are* plate wells and you want the layout and
  edge-effect figures, opt in with `well_numbering: plate` and key `wells:` (and
  `imaged_wells:`/`seq_ranges:`) by real well numbers.

* **`spinning_disk` mode**: the Spinning-Disk confocal's 'GeneratePoints' writes
  an **alphanumeric** well label, `WellC10_*.nd2` → row C, column 10. Key the
  `wells:` block by label instead of by number, and same-row column ranges work
  the same way:
  
  ```yaml
  mode: spinning_disk
  plate: 96well
  wells:
    C10:   {treatment: Control, class: control}
    C1-C4: {treatment: DrugX, dose_uM: 10}
  ```
  
  `WellE06` and `WellE6` are the same well. Labels resolve to canonical
  row-major well numbers (`C10` on a 96-well plate → 34), so everything
  downstream sees the same numeric `well` regardless of mode. Cross-row ranges
  (`B10-C2`) are rejected rather than guessed at.

* **Per-well stains**: wells default to the `default` stain set; tag the
  exceptions with `stain: <set>`. Each manifest row then carries its own
  `ch{n}_marker` columns, so a `groupby` never mixes markers.

> **Filenames carry no marker meaning.** The only things read out of a filename
> are the well number (`Well07`), the point index (`Seq0042`) or the well label
> (`WellC10`) that the acquisition software writes. Which markers an image
> carries is whatever `plate_map.yaml` says it is: never guessed from the name.

> **`plate_id` matters when panels live in sibling folders.** It defaults to the
> folder name, so `zmax_stain1/` and `zmax_stain2/` get *different* plate ids and
> are treated as different batches. If they are two antibody panels of the same
> physical plate, give both maps the same explicit `plate_id:`: that is what
> makes cross-panel gate pooling possible (see
> [Gate scope](#gate-scope-pooling-a-marker-across-panels)).

#### Reserved well keys

Four keys in a `wells:` entry are structural rather than experimental
conditions. Everything else you write becomes a free-form column you can group,
facet and plot by.

| key            | meaning                                                                                                                    |
| -------------- | -------------------------------------------------------------------------------------------------------------------------- |
| `stain`        | selects a named entry from `stains:`                                                                                       |
| `class`        | marks the well as the negative control (values configured in `config.yaml`)                                                |
| `positive_for` | markers this well is the **positive** control for: it stays the negative control for every other marker                   |
| `control_for`  | markers this well **is** the negative control for (`all` accepted); an explicit name beats a blanket declaration elsewhere |

The last two exist for the case where one condition is a negative control for
*some* of the panel and a positive control for the rest: see
[Per-marker negative controls](#per-marker-negative-controls).

### 1b. Several stain sets in one folder

Normally each stain set gets its own folder, and each folder is a `datasets:`
entry in `config.yaml`: that is the simplest arrangement and needs nothing
extra. If instead you acquired several stain sets into a **single** folder, the
per-well `stain:` key can't express it (the same well appears once per stain),
so declare the Seq blocks:

```yaml
stains:
  stain1: {0: DAPI, 1: ISL1, 2: NKX2.5}
  stain2: {0: DAPI, 1: HCN4, 2: MEF2C}

stain_sweeps:                       # non-overlapping, inclusive, 0-based
  stain1: {seq_range: [0, 107],   points_per_well: 6}
  stain2: {seq_range: [108, 215], points_per_well: 6}
```

Each block is that stain's own sweep of the plate, so **its Seq indices are
local to the block**: with the above, `Seq0108` is stain2's *first* point and
maps to well 1, not to well 19. Each block takes the same mapping options as the
plate as a whole (`points_per_well`, `imaged_wells`, `seq_ranges`), and
overlapping blocks are rejected outright.

Stain resolution order is `stain_sweeps` → the well's `stain:` key → `default`.

### 1c. Check before you run

```bash
python main.py check
```

Reads the config, the plate maps and the file listing (no images, no GPU, no
writes) and reports per dataset: how many images map to a well, how many wells
have conditions, which stain set each image gets with its channel→marker map,
and any Seq blocks. It fails loudly on images that map to no well, wells with
images but no condition entry, and stain sets with no channel map; it warns on
wells declared but never imaged. Run it first: these are the mistakes that
otherwise surface after an hour of segmentation.

It also answers the two questions that decide whether your thresholds will mean
anything:

* **does every marker have a negative control?** It lists the control wells, any
  `control_for` / `positive_for` redirection, and errors when a marker has been
  demoted on every candidate control well and has nowhere left to be gated from.
* **is any marker poolable across panels?** When a marker appears in two panels
  of one plate on the *same* channel, it names it as a candidate for
  `scope: plate`. On a different channel it says so too: there, separate gates
  are almost certainly correct.

### 2. Build the manifest

```bash
python main.py manifest                                   # auto-discover all plate maps
python main.py manifest --plate-map data/exp1/plate_map.yaml   # a single plate
python main.py manifest --data-dir /path/to/data --out run42_manifest
```

This writes `manifest.parquet` + `manifest.csv` at each dataset root, one row per
image: `dataset, plate_id, image_path, result_dir, well, seq_index,
<your condition columns…>, stain_set, ch0_marker …`, plus `row, col` when the
well number is a real plate coordinate (see `well_numbering` above).

### 3. Join it to the per-nucleus results downstream

```python
import joblib, pandas as pd
from pathlib import Path

root = Path("data")
m = pd.read_parquet(root / "manifest.parquet")

frames = []
for _, r in m.iterrows():
    res = joblib.load(root / r.result_dir / "per_nuclei_intensities.joblib")
    df = pd.DataFrame({"dapi": res["dapi_per_cell"], "ch2": res["raw_intensities"][2]})
    df = df.assign(well=r.well, treatment=r.treatment, ch2_marker=r.ch2_marker)
    frames.append(df)

cells = pd.concat(frames, ignore_index=True)
cells.groupby(["treatment", "ch2_marker"])["ch2"].median()   # per-condition summary
```

`manifest` and `process` never depend on each other: you can build the manifest
before, after, or without ever running segmentation.

---

## Downstream analysis

`python main.py downstream` turns the per-image joblibs + manifests into tables
and figures. It is generic: it knows nothing about your markers or treatments,
reading all of that from the plate map and config.

```bash
python main.py downstream                       # every step
python main.py downstream --steps assemble qc   # just these
python main.py downstream --force-step thresholds
python main.py downstream --list-steps
```

### Three ideas the design rests on

**1. A tidy cell table is the interface.** `assemble` writes one pair of tables
for the whole experiment:

```
analysis/cells.parquet        one row per nucleus: ids, geometry, DAPI,
                              + every condition column from the plate map
analysis/intensities.parquet  one row per (nucleus × marker): raw, dapi_norm
```

Everything after that is a `groupby`. Any plate-map column is therefore a
grouping, faceting or x-axis variable for free, which is why time courses, dose
responses and factorial screens need configuration rather than code. `assemble`
is the only IO-heavy step; the rest run in seconds off parquet.

**2. Cells are never the replicate.** A well holds thousands of cells, so treating
cells as replicates makes everything "significant". Each field is reduced to one
value per marker, and the replicate is chosen by `stats.replicate: auto`:
**experiment** when `stats.experiment_column` names independent differentiations
and every condition has ≥ `min_experiments_for_unit` (2) of them; else **well**
when every condition has ≥ `min_wells_for_well_unit` (3) wells (fields of one well
share a seeding and a staining reaction, so they are averaged per well first);
else **field**, with a caveat in the report. `tables/replicate_unit.json` and the
top of the report say which. Comparisons are permutation tests on the replicate
means (exact when the design is small; `smallest_possible_p` is listed), with a
Welch CI and Hedges' g; a side with fewer than two replicates at the chosen unit
steps down a level and the row says `pseudo_replicated`. Cell-level distributions
are still plotted (they are the biology), but they never set an error bar.

**3. Steps are plugins.** Analyses register themselves, so your experiment adds
its own without editing the package (see [Extending](#extending-the-pipeline)).

### The steps

| Step            | Produces                                                                                                                       |
| --------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| `signaling`     | only for signaling markers (β-catenin, SMAD2/3, …): per image, nucleus core, cytoplasmic ring and background per nucleus (`signaling_nuclei.parquet` next to each image), `signaling_background.csv`, ring QC pictures in `figures/signaling/` |
| `assemble`      | `cells.parquet`, `intensities.parquet` (+ `nuc_bgsub`, `nc_ratio` for signaling markers), `assemble_summary.csv`, `field_channel_qc.csv`; each image labeled with the mode its masks were really made in |
| `qc`            | `qc_fields.csv`, `qc_marker_fields.csv`, `qc_saturated_nuclei.csv` (nuclei at the detector ceiling, left out of everything after) + field count/DAPI QC figure; flags outlier fields and suspect marker channels; warns when a marker is divided by a histone nuclear label |
| `thresholds`    | `thresholds.csv` (incl. gate-stability columns), `control_assignment.csv`, `signaling_reference.csv` + the evidence panel behind every gate; signaling markers get no gate unless `signaling.gate: true` |
| `aggregate`     | `fields.csv`, `wells.csv`, `conditions.csv` (mean ± CI across replicates, `n_replicates`, `replicate_unit`), `replicate_unit.json`, `well_uniformity.csv`, `comparisons_*.csv` (permutation p, `smallest_possible_p`, CI, Hedges' g, q) |
| `gating`        | `band_edges.csv`, `band_composition.csv`, `coexpression.csv`                                                                   |
| `distributions` | violins (cells + field medians + gate) and ECDFs per marker                                                                    |
| `composition`   | positive-fraction bars, graded band bars, co-expression grids, quadrants                                                       |
| `joint`         | two-marker density and scatter panels with gates drawn on (scatters colored by a third marker, see below), correlation matrix |
| `plate`         | plate-geometry heatmaps, per-well field-uniformity strips                                                                      |
| `trends`        | time/dose courses, factorial marginal effects (at the replicate unit), interaction heatmaps                                    |
| `effects`       | from every per-field number: `effects_variance.csv` (field / well / experiment SDs), `effects_vs_reference.csv`, `effects_trends.csv` (slope per decade with bootstrap CI and permutation p), `effects_power.csv` (minimum detectable difference vs number of experiments), `figures/effects/` |
| `report`        | `report.html`: one scrollable contact sheet of everything                                                                     |

`gating`, `distributions`, `composition`, `joint`, `plate`, `trends` and
`effects` are optional: list them in `downstream.skip_steps` (GUI step 3,
*Analyses to run*) to leave them out; their earlier outputs move to
`<analysis>/_skipped/<step>/`. The rest is the core. Every marker figure carries a
line saying what its values are normalized to.

New options (2026-11): `results_dir: auto` (→ `analysis_2d` / `analysis_3d` /
`analysis_2d_3d`, after the segmentation results used), `seg_mode` (auto | 2d |
3d; also `--seg-mode`), `exclude_markers` (by name; DAPI always out),
`skip_steps`, `nuclear_stain`, `qc.drop_saturated_nuclei` /
`saturation_ceiling`, `stats.replicate` / `min_wells_for_well_unit` /
`experiment_column` / `n_permutations`, the `effects` block and the `signaling`
block. Every option is explained in `templates/config.yaml`; the statistics in
[STATISTICAL_METHODS.md](STATISTICAL_METHODS.md).

Each writes `.checkpoints/<step>.json` holding a hash of the config subtrees it
depends on, so a step re-runs only when something it actually cares about
changed. Change `thresholds.fpr` and the threshold-dependent steps redo
themselves while `assemble` does not.

### Joint marker views

For every marker pair A × B in a panel, `joint` writes into `figures/joint/`:

- `density_A_vs_B.png`: a 2-D histogram per condition (log color scale shared
  across panels).
- Scatter plots on the same per-condition grid, axis ranges and gate overlays
  (dashed = positivity gate, dotted = tertile edges):
  - **Three or more markers:** one figure per remaining marker C,
    `scatter_A_vs_B_by_C.png`. Each cell is colored by its C intensity (viridis),
    on one color scale shared by all conditions (1st percentile to
    `cap_percentile`; the colorbar's triangles mark values beyond that range).
    C's positivity gate is drawn as a dashed line on the colorbar. Cells are
    drawn in ascending C, so C-high cells sit on top of the C-low bulk. Only cells
    measured for A, B and C are shown. With three markers this gives three
    scatters, and with four markers twelve (six pairs × two colorings).
  - **Two markers:** one plain single-color scatter, `scatter_A_vs_B.png`.
- `marker_correlation.png`: the Spearman matrix.

Point size and opacity come from `visualization.joint_scatter_size` (default 2)
and `visualization.joint_scatter_alpha` (default 0.8).

### Positivity thresholds

One threshold per **scope × marker**, applied to every well in that scope. Never
per well, never per treatment.

The default method is the negative-control percentile: the threshold is the
`1 − fpr` quantile of that marker's own control values, i.e. an empirical
false-positive-rate cut (`fpr: 0.01` → the control's P99). This is Eq. (2) of
`Cell_fate_classification.pdf` with `fpr` in place of the paper's P95.

When a marker has no usable control, a 2-component Gaussian mixture is fitted to
the pooled distribution of every well in that scope, but only accepted if the
data are **demonstrably bimodal**: BIC gain over a 1-component fit, standardized
separation between the modes, and a minimum weight for the smaller component.
If any check fails, **no threshold is emitted**: the marker is marked
`graded_only`, its distributions are still plotted, and no positive fraction is
reported for it anywhere (NaN, never a silent 0 %). The reason is recorded in
`thresholds.csv` and stated at the top of `report.html`.

#### Gate scope: pooling a marker across panels

The **scope** is the pool a marker's threshold and band edges are estimated over.
It is declared per marker, and defaults to `panel`:

| level               | grouping columns                 | pools                   |
| ------------------- | -------------------------------- | ----------------------- |
| `panel` *(default)* | `dataset × plate_id × stain_set` | nothing                 |
| `plate`             | `plate_id × seg_mode`            | panels **and** datasets |
| `experiment`        | `seg_mode`                       | plates as well          |

The `panel` default encodes two claims: the **plate is the batch** (staining and
exposure differ between plates, so a cut is never shared across `plate_id`), and
a marker in **two panels gets two thresholds** because a different secondary,
filter or exposure makes the raw scales incomparable.

That second claim is a default, not a law. When the *same* antibody at the *same*
dilution appears in several panels across the wells of one plate (say
`[NKX2.5, ISL1, HAND1]` in some wells and `[MEF2C, ISL1, HAND2]` in others):
splitting ISL1 into two gates makes its positive fractions incomparable between
wells of one plate for no reason. Declare the wider scope:

```yaml
downstream:
  markers:
    ISL1: {scope: plate}          # one ISL1 gate, pooled over both panels
    TNNT2: {space: raw}           # unrelated: which intensity space it is gated in
```

ISL1 then gets a single gate estimated from every ISL1 cell on the plate, one set
of tertile edges, and one figure per marker instead of two suffixed ones.
`seg_mode` stays in the key at every level, so the same plate imaged as MIPs and
as z-stacks never shares a gate: those are two measurements of the same cells,
not two panels.

Two things to know:

* **Sibling panel folders become separate datasets.** `plate` scope drops
  `dataset` so it can pool them, but only if their plate maps declare the *same*
  explicit `plate_id`. Since `plate_id` defaults to the folder name, forgetting
  this makes pooling silently do nothing.
* **Nothing is pooled without you saying so.** A marker that *looks* poolable
  (two panels, one plate, same channel index) is reported by `check` and by the
  `thresholds` step as a candidate, with the marker name and the config line to
  add. Whether the antibody and dilution actually matched is a fact about the
  bench that no file records, so the pipeline names the opportunity and leaves
  the claim to you. On a *different* channel it says that too, and there the
  separate gates are almost certainly right.

Markers at different scope levels coexist in one run. `thresholds.csv` carries
the union of the grouping columns and writes `*` in the ones a marker was pooled
over, alongside `scope_level`, `scope_key` and `n_panels_pooled`.

#### Per-marker negative controls

The blanket `conditions.control_column` / `control_values` declaration marks a
well as the negative control **for every marker measured in it**. That breaks as
soon as one condition is a negative control for part of the panel and a positive
control for the rest: the normal situation in a pluripotency experiment, and
exactly Table S1 of `Cell_fate_classification.pdf`:

| marker | mTeSR (pluripotent)  | BMP ≥10 ng/ml        |
| ------ | -------------------- | -------------------- |
| SOX2   | **positive** control | **negative** control |
| BRA    | negative control     | positive control     |
| CDX2   | negative control     | positive control     |

Gate SOX2 off the pluripotent well and its threshold lands at the P99 of a
saturated population: nothing is called positive, and the failure is silent. The
same mistake hits the tertile bands, which exclude control wells from the
positive pool: SOX2's edges would be computed after discarding the brightest
cells in the experiment.

You fix it by writing down only the **exception**. Two lines in the plate map:

```yaml
wells:
  1-6: {treatment: mTeSR, class: pluripotent, positive_for: [SOX2]}
  7-9: {treatment: BMP50, control_for: [SOX2]}
```

`positive_for` demotes wells 1-6 for SOX2 alone: they stay the negative control
for every other marker. `control_for` names where SOX2's negative control
actually lives. Everything else on the plate keeps inheriting the pluripotent
control, and no marker name is hardcoded anywhere in the pipeline.

Resolution runs most-specific-first:

1. **plate map `control_for`** naming the marker explicitly;
2. **`conditions.marker_controls`** in `config.yaml`: the same override from the
   config side, for reanalysis without editing a plate map:
   
   ```yaml
   conditions:
     marker_controls:
       SOX2: {control_values: [bmp_high]}
   ```
   
   Skipped when the plate map already spoke, with a warning naming both;
3. **the blanket `class` declaration, minus any `positive_for`.**

A marker left with no control falls through to the GMM fallback as before:
reached honestly rather than by gating off the wrong population. Every marker's
resolved control is written to `control_assignment.csv` (cells, wells, conditions
and which rule fired) and to a `control_source` column in `thresholds.csv`.

This also changes the statistics: **each marker is compared against its own
negative control** in `comparisons_*.csv`, since comparing SOX2 to its positive
control would report the effect backwards. An explicit
`stats.reference_condition` still overrides everything.

### Guardrails: catching a broken measurement

A threshold is a quantile of a control population, so it inherits whatever is in
that population. The failure that follows is the quiet kind: contaminate part of
a **control** well (a patch of non-specific antibody, a bubble, one over-bright
z-slice in a projection) and the cut climbs, the marker reports *fewer*
positives, and it reads as biologically negative rather than obviously broken.
Nothing about the per-nucleus numbers looks wrong afterwards.

Three checks exist to make that visible, at the three places it can be seen.
**All of them are detectors**: they write columns and log warnings, and none of
them changes a single analyzed number. That is deliberate: a guardrail that can
quietly alter a result becomes something you have to reason about in every
future experiment.

**1. Channel health, at extraction.** Each channel gets three scalars in
`global_stats` (and in `analysis/field_channel_qc.csv`):

| field                  | meaning                                                              |
| ---------------------- | -------------------------------------------------------------------- |
| `saturated_fraction`   | pixels at the detector ceiling: information no analysis can recover |
| `background`           | median over the whole channel, i.e. the pedestal the nuclei sit on   |
| `signal_to_background` | nuclear median ÷ background: **the one number to look at**          |

A healthy channel has `signal_to_background` well above 1. At 1 the nuclei are
indistinguishable from what they sit on, which is the signature of a flooded
field, and specifically of a MIP containing a bright z-slice, where a maximum
projection *floors* the image, overwriting every pixel dimmer than the flare
rather than adding to it. Set `intensity.saturation_value` for a camera that
digitizes fewer bits than its container (a 12-bit sensor clips at 4095 and would
otherwise never register as saturated).

**2. Per-marker field QC, before thresholding.** `qc_fields.csv` looks only at
nucleus count, DAPI and nuclear size, so an artifact confined to one antibody
passes it untouched. `qc_marker_fields.csv` compares each field's marker median
against **the other fields of its own well**: the only comparison where
disagreement must be technical, since one well is one staining reaction while
different conditions are *supposed* to differ. A field is flagged when it fails
both a robust-z test and a relative-deviation floor (`qc.marker_min_rel_dev`);
requiring both matters, because at the 3-6 fields per well most plates carry, the
MAD of a few tightly agreeing numbers collapses toward zero and a field sitting
3 % off its neighbors would otherwise score a z of 8.

**3. Gate stability, at thresholding.** Every gate is recomputed with each control
field held out in turn, and the spread lands in `thresholds.csv` as
`gate_loo_min/max`, `gate_loo_spread_pct`, `gate_loo_worst_field` and
`n_control_fields`. Fields of one control well are replicates of one population,
so a healthy gate barely moves; one that jumps when a single field leaves is
describing that field, not the control. Crossing `thresholds.loo_warn_pct`
(default 25 %) logs a warning naming the field, and a gate resting on a single
control field always warns.

None of this drops anything. `qc.drop_outlier_fields` still acts on the
field-level flags only: one bad marker channel is not a reason to discard that
field's other, unaffected markers, and excluding a field remains an explicit
decision you make and record.

### Graded gating and co-expression

A binary call discards most of what a distribution says, so each marker's
**positive** population is split into tertiles, giving four ordered bands:
`neg / low / mid / high`. Edges are the 33rd/67th percentiles of the positive
cells pooled over the **test** wells of that plate × panel, so they are one
absolute reference shared by every condition (per-condition tertiles would
compare each condition only with itself). `gating.split: median` gives the 50/50
variant.

Crossing two markers' bands gives the 4×4 co-expression matrix per condition:
what fraction co-express, and in which part of each distribution, which
collapses to the familiar quadrant view in `quadrants_*.png`.

### Configuring it for your experiment

```yaml
downstream:
  conditions:
    condition_key: [treatment]      # or [treatment, dose_uM] for a composite
    facet_key: null                 # e.g. [treatment]: one figure per treatment
    control_column: class           # plate-map column declaring controls
    control_values: [pluripotent, control]
    marker_controls: {}             # exceptions only, e.g. {SOX2: {control_values: [bmp_high]}}

  default_space: dapi_norm          # nuclear markers
  markers:
    TNNT2: {space: raw}             # membrane/cytoplasmic markers stay raw
    ISL1:  {scope: plate}           # one gate pooled across the panels it appears in

  thresholds: {scope: panel}        # default gate scope for undeclared markers

  gating: {split: tertile, pairs: auto}

  trends:
    x: timepoint                    # any numeric plate-map column
    series: treatment
    facet: day
    factors: null                   # null → auto-detect two-level columns
```

#### Crossed designs: one figure per level instead of one crowded axis

When conditions are a product of two factors (three treatments × four progenitor
origins), putting all twelve on one axis buries the comparison you care about.
`facet_key` holds one factor fixed per figure and leaves the rest on the x-axis:

```yaml
conditions:
  condition_key: [treatment, ps_origin]   # a condition is still the full combination
  facet_key: [treatment]                  # …but each treatment gets its own figure
```

You get `violin_ISL1__Activin.png`, `violin_ISL1__Activin_Noggin.png`, … each with
`ps_origin` on the x-axis and the treatment named in the title. `ecdf_*` and
`positive_fraction_*` split the same way.

This is presentational only, and deliberately so: the facet column must be part
of `condition_key`, so every threshold, bootstrap CI and q-value is still computed
on the full condition. Faceting changes which panel a number is drawn in, never
how it was calculated. A `facet_key` naming a column outside `condition_key`, or
using up all of them, is ignored with a warning rather than silently splitting a
condition across panels.

The component columns also travel into `conditions.csv` and `wells.csv` as their
own columns, so the tables pivot by factor without re-parsing the joined label.

### Extending the pipeline

Put your experiment's own analysis next to the config and register it:

```python
# downstream_analysis/custom.py
from downstream import register

@register("fate_calls", after=["gating"], config_keys=["gating"])
def step_fate_calls(ctx):
    banded = ctx.banded                      # bands the generic step computed
    ...
    ctx.write_table(table, "fate_proportions.csv")
    ctx.reset_figures("fate")                # each step owns one figure group
```

```yaml
downstream:
  extra_modules: [downstream_analysis.custom]
```

Your step is then a first-class citizen: dependency-ordered, checkpointed,
selectable with `--steps` / `--force-step`, and its figures appear in
`report.html`. A complete worked example is in
[`examples/custom_analysis.py`](../examples/custom_analysis.py).

### Figures

All color comes from `src/downstream/style.py`: the **Okabe-Ito**
color-vision-deficiency-safe palette for categories, perceptually-uniform
viridis/cividis/magma for continuous data, and line style + marker varied
alongside color so figures survive grayscale printing. `tests/test_downstream.py`
fails if any hex literal appears elsewhere in the package, so the palette cannot
drift.

### Tests

```bash
python tests/test_downstream.py
python tests/test_experiment.py        # or: python -m pytest tests/
```

Synthetic fixtures only: no network, no image data, about a second. They cover
the threshold methods (including that a unimodal marker is *refused* a gate), the
no-call-propagates-as-NaN rule, band assignment, the bootstrap, config merging,
the palette rule, and, in `test_experiment.py`, Seq→well mapping, stain
resolution, stain sweeps being local to their block, image discovery pruning
(result dirs, exports, dot-dirs, sub-plates), and the guarantee that marker names
appearing in a filename never influence anything.

The guardrails are tested against a contaminated fixture in both directions: a
poisoned control field must raise the gate *and* be named by the leave-one-out
spread, a clean control must stay under the warning threshold, and the per-marker
field flag must fire on exactly the contaminated (field × marker) and on nothing
else: including the case that motivated the relative-deviation floor, where a
3 % deviation among tightly agreeing fields scores a robust z of 8.

---

## Project structure

```
.
├── main.py                       # CLI: process | init-plate-map | manifest | check | downstream
├── templates/config.yaml         # defaults + datasets + downstream settings (commented)
├── docs/                         # REFERENCE.md (this file), TUTORIAL.md
├── requirements.txt
├── examples/custom_analysis.py   # worked example of an experiment-specific step
├── tests/                        # fast synthetic tests (no data needed)
└── src/
    ├── pipeline_config.py        # defaults + datasets → N resolved dataset configs
    ├── image_processing/         # generic, condition-agnostic core
    │   ├── io.py                 # versatile multi-format reader → canonical CZYX
    │   ├── directory_walker.py   # recursive, output-aware image discovery
    │   ├── image_loader.py       # MIP/stack detection, channel TIFFs, metadata
    │   ├── segmentation.py       # Cellpose-SAM 2-D/3-D + adaptive filter + QC
    │   ├── intensity_extractor.py# per-nucleus intensities + centroids + channel health
    │   └── utils.py              # shared metadata / channel-path helpers
    ├── experiment/               # well/condition annotation layer (decoupled)
    │   ├── plate.py              # plate geometry + vectorized well↔(row,col)/label
    │   ├── parser.py             # filename → well number / seq index / well label
    │   ├── plate_map.py          # load + validate plate_map.yaml
    │   ├── template.py           # scan a folder → commented starter plate_map.yaml
    │   └── manifest.py           # walk data → join conditions → manifest.parquet/.csv
    └── downstream/               # generic analysis layer
        ├── data.py               # the cells/intensities star schema + fast joins
        ├── assemble.py           # manifests + joblibs → the cell tables
        ├── qc.py                 # field QC, per-marker channel QC, outlier flags, uniformity
        ├── scope.py              # per-marker gate scope (panel/plate/experiment)
        ├── controls.py           # per-marker negative-control resolution
        ├── thresholds.py         # control-percentile / validated-GMM gating + gate stability
        ├── gating.py             # tertile bands + co-expression matrices
        ├── aggregate.py          # cell → field → well → condition rollups
        ├── stats.py              # bootstrap CIs, effect sizes, factor effects
        ├── registry.py           # @register step plugins + checkpoint hashing
        ├── context.py            # lazily-loaded analysis products per step
        ├── runner.py             # ordering, checkpoints, keep-going execution
        ├── steps.py              # the built-in steps
        ├── report.py             # the HTML contact sheet
        ├── style.py              # THE palette (Okabe-Ito) + publication rcParams
        └── plots/                # figure builders (pure functions)
```

### Extending

- **New file format**: add a `_read_<fmt>` backend in `io.py` returning
  `(data, axes, meta)` and register its extension; everything downstream is
  format-agnostic.
- **New per-image step**: add a class with a `process(path, force=False)`
  method and a checkpoint file, then register it in `STEP_ORDER` /
  `build_steps()` in `main.py`.
- **New downstream analysis**: `@register(...)` in your own module; see
  [Extending the pipeline](#extending-the-pipeline) above.

The core helpers (`io.read_image`, vectorized label statistics, adaptive
filtering, overlays, the cell-table schema) are written to be reusable across any
image-analysis task.
