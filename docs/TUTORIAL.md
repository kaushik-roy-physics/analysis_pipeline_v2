> This tutorial explains the ideas behind plate maps and configs using the command
> line. The GUI (see the main [README](../README.md)) writes the same files for you,
> so everything here still applies.

# Tutorial: running the confocal pipeline on your own experiment

For someone who has never used this pipeline. It assumes you know your own
experiment, which wells got what, which antibody sat in which channel, and
nothing else. Work through parts 1-4 once; after that, parts 5-9 are a recipe
book you come back to when a new experiment looks different from the last one.

The README is the reference manual: every option, precisely defined. This is the
practical companion: what to do, in what order, and why the pipeline asks what
it asks.

---

## 1. The one idea to understand first

The pipeline is three layers that know nothing about each other:

```
   images  ──►  per-image results        (load → segment → extract)
                 "how bright is each nucleus"     knows nothing about your experiment

   plate_map.yaml ──►  manifest          (manifest)
                 "which well is which condition"  never touches an image

   both together  ──►  analysis          (downstream)
                 "what does it mean"              pure tables and figures
```

Two consequences that save you a lot of confusion:

**Filenames carry no meaning.** The only thing ever read out of a filename is the
well number or point index the microscope wrote. In the dataset this tutorial
draws on, files are named `..._561_nkx25_...` while that channel actually holds
HAND2: the acquisition template was stale. It made no difference to anything,
because what a channel contains is whatever `plate_map.yaml` says it contains. If
you ever catch yourself renaming files so the pipeline will understand them,
stop: edit the plate map instead.

**You can build the manifest before you have segmented anything.** `manifest` and
`process` never depend on each other. Author your plate map on day one, run
`check`, and fix your thinking about the experiment before spending an hour of
GPU time.

---

## 2. The whole workflow, six commands

The code lives in one place; each experiment keeps its own config next to its own
data. Every relative path in a config resolves against **the directory holding
that config**, never the directory you ran from, so this works from anywhere:

```bash
PIPE=/path/to/KR/project_structure          # the code, one copy
EXP=/path/to/KR/20260801_my_experiment/plate1   # your data + config

python "$PIPE/main.py" init-plate-map --config "$EXP/config.yaml"   # 1. scaffold
#   ... edit plate_map.yaml ...                                      # 2. think
python "$PIPE/main.py" check      --config "$EXP/config.yaml"       # 3. validate
python "$PIPE/main.py" process    --config "$EXP/config.yaml"       # 4. GPU work
python "$PIPE/main.py" manifest   --config "$EXP/config.yaml"       # 5. annotate
python "$PIPE/main.py" downstream --config "$EXP/config.yaml"       # 6. analyze
```

Steps 1-3 take five minutes and no GPU. **Do not skip step 3.** `check` reads the
config, the plate maps and the file listing, and tells you whether every image
maps to a well, whether every well has a condition, which markers each image gets
and whether every marker has a negative control to be gated against. These are
the mistakes that otherwise surface after an hour of segmentation.

Everything is checkpointed. Re-running only redoes what is missing; `--force-step
segment` or `downstream --force-step thresholds` redoes one thing on purpose.

---

## 3. Starting a new experiment

```bash
cp /path/to/previous_experiment/config.yaml "$EXP/config.yaml"
```

Copy from a *similar* experiment, then read every line you inherited. A config
carries opinions about the last experiment: its treatments, its reference
condition, its marker list. Inheriting those silently is the single most common
way to get a confident wrong answer out of this pipeline.

Then point it at your folders:

```yaml
datasets:
  - name: plate1_zmax
    directory: ./zmax
```

`name` becomes the `dataset` column, which is what lets downstream pool or facet
across folders. Keep `directory` **relative** and the whole experiment folder
stays portable.

Now scaffold the plate map:

```bash
python "$PIPE/main.py" init-plate-map --config "$EXP/config.yaml"
```

This scans the folder and writes a heavily commented `plate_map.yaml` pre-filled
with what it actually found: the acquisition mode detected from the filenames,
the Seq or well range present, the images-per-well count, and (if you have
already run `process --steps load`) the real channel names. Lines you must edit
are marked `<-- EDIT`. It never overwrites an existing plate map without
`--force`, because that file is hand-authored knowledge.

---

## 4. Anatomy of a plate map

Here is a complete, real one: three antibody panels across 15 wells of an
18-well plate, five fields each. Read it top to bottom; every block is explained
below.

```yaml
plate: 18well                       # 8well | 18well | 96well, or geometry: {rows, cols}
plate_id: 20260801_LM_PScoords      # the BATCH identity; see §7
mode: predefined                    # how filenames encode position

points_per_well: 5                  # 5 fields per well, wells filled in order

stains:                             # channel index → marker, per panel
  stain1: {0: DAPI, 1: ISL1, 2: HAND1, 3: NKX2-5}
  stain2: {0: DAPI, 1: ISL1, 2: HAND2, 3: MEF2C}
  stain3: {0: DAPI, 1: ISL1, 2: FOXF1, 3: TBX5}

wells:
  1:  {stain: stain1, treatment: LM, ps_origin: DEPLM, ps_axis: 0}
  2:  {stain: stain1, treatment: LM, ps_origin: APSLM, ps_axis: 1}
  3:  {stain: stain1, treatment: LM, ps_origin: MPSLM, ps_axis: 2}
  4:  {stain: stain1, treatment: LM, ps_origin: PPSLM, ps_axis: 3}
  5:  {stain: stain1, treatment: NegCtrl, ps_origin: pluripotent, class: control}
  6-10:  # …the same five conditions stained with stain2
  11-15: # …and with stain3
```

### Well keys: four are reserved, the rest are yours

Everything you invent under a well (`treatment`, `dose_uM`, `timepoint`,
`ps_origin`) becomes a column you can group, facet, plot and test by. Invent
freely; just spell things the same way everywhere. Four keys are structural
instead:

| key | meaning |
| --- | --- |
| `stain` | picks a named entry from `stains:` |
| `class` | marks this well as the negative control |
| `positive_for` | markers this well is the **positive** control for |
| `control_for` | markers this well **is** the negative control for |

`class` is the most important annotation in the file: positivity thresholds are
a percentile of the control population, so getting it wrong moves every number
downstream. The last two exist for one specific situation, covered in §7.

### Range keys

`2-4` means wells 2, 3 and 4. Use it freely; it keeps a 96-well map readable.

---

## 5. Which acquisition mode am I in?

Look at one filename and match it:

| filename looks like | mode | what the pipeline reads |
| --- | --- | --- |
| `Well07_....nd2` | `generate` | the well **number**, directly |
| `..._Seq0042.nd2` | `predefined` | a global 0-based **point index** |
| `WellC10_....nd2` | `spinning_disk` | an alphanumeric **well label** (row C, col 10) |

**`generate`** is the easy one: the filename says the well, and using only some
wells needs nothing special: you simply have no files for the rest.

**`predefined`** is the one that confuses everybody once. `Seq0042` is
*acquisition order*, not plate position. The filename says nothing about where on
the plate the point came from, so wells here are numbered **sequentially**: well
1 is the first well you imaged, well 2 the second. If you imaged 13 wells of an
18-well plate in a zig-zag, they are still wells 1-13 in this file, and you never
have to look up physical positions. Say how points divide up in one of three
ways:

```yaml
points_per_well: 5                    # equal counts, wells 1..K in order; the usual case
points_per_well: {1: 5, 2: 8, 3: 3}   # variable counts
seq_ranges: {1: [0, 4], 3: [10, 14]}  # explicit, inclusive, 0-based; overrides everything
```

Because there is no plate position, `row`/`col` are not written and the
plate-layout heatmaps skip themselves (the step logs why). If your keys really
*are* physical wells and you want those figures, opt in with `well_numbering:
plate` and key everything by real well numbers.

**`spinning_disk`** writes alphanumeric labels. Key `wells:` by label; same-row
ranges work (`C1-C4`), cross-row ranges are rejected rather than guessed at.
`WellE06` and `WellE6` are the same well.

---

## 6. Recipe book: five shapes an experiment can take

### A. One panel, one folder

The base case. Nothing special needed.

```yaml
# config.yaml
datasets:
  - {name: plate1, directory: ./zmax}

# plate_map.yaml
stains:
  default: {0: DAPI, 1: SOX2, 2: BRA, 3: CDX2}
wells:
  1-3: {treatment: mTeSR, class: control}
  4-6: {treatment: BMP50}
```

### B. Several panels, one folder per panel

Each panel is its own folder and its own dataset. **Give the plate maps the same
explicit `plate_id`**: it defaults to the folder name, so without this the two
panels count as different batches and cross-panel pooling silently does nothing.

```yaml
datasets:
  - {name: zmax_stain1, directory: ./zmax_stain1}
  - {name: zmax_stain2, directory: ./zmax_stain2}
```

```yaml
# in BOTH zmax_stain1/plate_map.yaml and zmax_stain2/plate_map.yaml
plate_id: exp033_plate1
```

### C. Several panels in ONE folder, in different wells

You stained 15 wells: five with panel 1, five with panel 2, five with panel 3.
Each well was stained **once**. Use a per-well `stain:` key: this is the case in
§4 above, and the plain `points_per_well: 5` mapping handles it because the wells
are simply visited in order.

### D. Several panels in ONE folder, over the *same* wells

You swept the whole plate with panel 1, then swept it again with panel 2. Now the
same well appears twice, so a per-well `stain:` key cannot express it. Declare
the Seq blocks instead:

```yaml
stain_sweeps:                          # non-overlapping, inclusive, 0-based
  stain1: {seq_range: [0, 107],   points_per_well: 6}
  stain2: {seq_range: [108, 215], points_per_well: 6}
```

Each block is that panel's **own sweep** of the plate, so its Seq indices are
local to the block: `Seq0108` is stain2's *first* point and maps to well 1, not
to well 19.

> **C or D?** Ask one question: *did the same physical well get imaged under two
> panels?* If no (C), the panels are different wells: use `stain:`. If yes (D),
> use `stain_sweeps:`. Getting this backwards is not cosmetic. Using sweeps for a
> case-C plate collapses three physically distinct wells into one "well 1"
> carrying fields from three different staining reactions, which corrupts the
> field→well rollup, the per-well uniformity QC, and the well count behind every
> error bar.

### E. The same plate as MIPs and as z-stacks

Two datasets, same `plate_id`, different segmentation mode:

```yaml
datasets:
  - {name: zmax,    directory: ./zmax}
  - {name: zstacks, directory: ./zstacks, segmentation: {mode: 3d, anisotropy: null}}
```

`seg_mode` stays in every gate key at every scope level, so the 2-D and 3-D
measurements of the same cells never share a threshold: they are two
measurements of one experiment, not twice the replicates.

---

## 7. Conditions, controls and the things that go wrong

### What is a condition?

```yaml
downstream:
  conditions:
    condition_key: [treatment]        # or [treatment, dose_uM] for a composite
    order: [Control, Low, High]       # explicit; the default is ALPHABETICAL
    control_column: class
    control_values: [pluripotent, control, negative_control]
```

Set `order:` whenever your conditions have a natural sequence. Alphabetical order
will otherwise scramble it: in the worked experiment, the anterior→posterior
streak axis DEP→APS→MPS→PPS sorts to APS, DEP, MPS, PPS, which silently reverses
the main independent variable in every bar chart.

Leave a factor out of `condition_key` if it is constant. With one treatment
across all arms, adding `treatment` puts the same word on every condition and
buys nothing but longer labels.

### Crossed designs: `facet_key`

Three treatments × four origins on one axis buries the comparison you care
about. Hold one factor fixed per figure:

```yaml
conditions:
  condition_key: [treatment, ps_origin]   # a condition is still the full combination
  facet_key: [treatment]                  # …but each treatment gets its own figure
```

This is **presentational only**, by design: the facet column must be part of
`condition_key`, so every threshold, bootstrap CI and q-value is still computed
on the full condition. Faceting changes which panel a number is drawn in, never
how it was calculated.

### Per-marker controls: the one that bites everybody

The blanket `class: control` declaration marks a well as the negative control for
**every marker measured in it**. That breaks the moment one condition is a
negative control for part of the panel and a positive control for the rest: the
normal situation in a pluripotency experiment:

| marker | mTeSR (pluripotent) | BMP ≥10 ng/ml |
| --- | --- | --- |
| SOX2 | **positive** control | **negative** control |
| BRA | negative control | positive control |

Gate SOX2 off the pluripotent well and its threshold lands at the P99 of a
saturated population: nothing is ever called positive, and the failure is silent.
Write down only the exception, in the plate map where it travels with the data:

```yaml
wells:
  1-6: {treatment: mTeSR, class: pluripotent, positive_for: [SOX2]}
  7-9: {treatment: BMP50, control_for: [SOX2]}
```

`positive_for` demotes wells 1-6 for SOX2 alone: they stay the negative control
for everything else. `control_for` says where SOX2's negative control actually
lives. `check` will tell you if a marker has been left with nowhere to be gated
from.

### One antibody in several panels

If the *same* antibody at the *same* dilution appears in several panels of one
plate, two separate gates make its positive fractions incomparable between wells
of one plate for no reason. Pool them:

```yaml
downstream:
  markers:
    ISL1: {scope: plate}       # one gate over every panel on this plate_id
```

`check` and the `thresholds` step name every marker that *could* be pooled (same
plate, several panels, same channel index) and leave the decision to you,
because whether the antibody and dilution actually matched is a fact about the
bench that no file records. On a *different* channel, separate gates are almost
certainly right. Verify it worked: `thresholds.csv` should show `scope_level=plate`
with `n_panels_pooled` > 1. If it says 1, your `plate_id`s differ.

---

## 8. Getting the analysis you actually want

### Trends need a *numeric* x

Any numeric plate-map column can be the x-axis of a trend figure. If your
independent variable is a named series rather than a number, add a numeric
companion column: this is the single most useful trick in the file:

```yaml
# plate_map.yaml: ps_origin is the name, ps_axis is the same thing as a number
wells:
  1: {ps_origin: DEPLM, ps_axis: 0}
  2: {ps_origin: APSLM, ps_axis: 1}
  3: {ps_origin: MPSLM, ps_axis: 2}
  4: {ps_origin: PPSLM, ps_axis: 3}
  5: {ps_origin: pluripotent}          # no coordinate → drops out of the trend
```

```yaml
# config.yaml
downstream:
  trends:
    x: ps_axis
    series: stain_set        # one line per panel
    x_tick_labels: {0: DEP-LM, 1: APS-LM, 2: MPS-LM, 3: PPS-LM}
    factors: []              # see below
```

Leaving the control's `ps_axis` blank is deliberate: it has no coordinate on that
axis, so it drops out of the trend figures instead of anchoring them at a
meaningless x.

Choosing `series:` well is worth a moment. With one shared antibody across three
panels, `series: stain_set` gives most markers a single line but gives the shared
one **three**: three independent staining reactions under one pooled gate, i.e.
a free technical-replication readout you did not have to plan.

### Factors

`factors: null` auto-detects any plate-map column with exactly two levels and
runs a factorial analysis on it. That is what you want for a crossed design
(`noggin: true/false`, `chir: true/false`). It is *not* what you want when your
only two-level column is something like `treatment: LM | NegCtrl`, where the
"marginal effect" is the tautology that differentiated cells differ from stem
cells. Set `factors: []` to turn it off.

### Per-marker settings

```yaml
downstream:
  default_space: dapi_norm     # nuclear markers
  markers:
    TNNT2: {space: raw}        # membrane/cytoplasmic; DAPI normalization is meaningless
    ISL1:  {scope: plate}      # pooled across panels
    NKX2-5: {fpr: 0.001}       # a stricter gate for one antibody
```

On `fpr`: it is the false-positive rate against the negative control, so `0.01`
puts the cut at the control's P99. A **smaller** fpr means a higher quantile, a
higher cut and fewer positives. Reach for a per-marker value when one antibody's
control is a much tighter spike than the rest of the panel's, so the shared P99
sits barely above background and calls nearly everything positive.

---

## 9. Reading the output

```
analysis/
├── cells.parquet          one row per nucleus + every plate-map column
├── intensities.parquet    one row per (nucleus × marker)
├── field_channel_qc.csv   measurement health per (field × channel)
├── report.html            ← start here
├── tables/
│   ├── qc_fields.csv          did each field image and segment properly
│   ├── qc_marker_fields.csv   did each marker channel behave, field by field
│   ├── thresholds.csv         every gate + the evidence + its stability
│   ├── control_assignment.csv which wells fed which marker's gate
│   ├── fields.csv / wells.csv / conditions.csv    the rollups
│   └── comparisons_*.csv      stats vs the reference condition
└── figures/
```

Which file answers which question:

| question | look at |
| --- | --- |
| "Is any of this trustworthy?" | `qc_fields.csv`, `qc_marker_fields.csv` |
| "Where did this threshold come from?" | `thresholds.csv` + the threshold evidence figure |
| "Which wells were the control?" | `control_assignment.csv` |
| "What is the number, with an error bar?" | `conditions.csv` |
| "Is the effect significant?" | `comparisons_*.csv` |
| "Show me everything" | `report.html` |

**n is never cells.** A well holds thousands of cells, so treating cells as
replicates makes everything "significant" no matter how small the effect. Each
field is reduced to one value per marker; then, with 3 or more wells per
condition, fields are averaged per well and **n = wells**; with fewer, n = fields
and the report says what that costs; with `stats.experiment_column` naming
independent differentiations, n = differentiations. Cell-level distributions are
still plotted (they are the biology), but they never set an error bar. When you
report an n, it is `n_replicates` in `conditions.csv`, and the unit is in
`replicate_unit.json`. Why, in plain words: `docs/STATISTICAL_METHODS.md`.

---

## 10. When a guardrail fires

Three checks watch the measurement itself. They only ever warn: they never
change a number, so it is on you to act.

**`... % of pixels at the detector ceiling`**: that channel is clipped and has
lost information no analysis can recover. Look at the image.

**`nuclei are only 1.00× their own background`**: the channel is flooded: the
nuclei are indistinguishable from what they sit on. In a maximum projection this
usually means one over-bright z-slice, which *floors* the projection, overwriting
every pixel dimmer than the flare. Background subtraction cannot recover it,
because the signal was overwritten rather than buried.

**`gate moves by 88 % when a single control field is dropped`**: the most
important one. Your threshold is describing one contaminated control field rather
than the control population. The warning names the field. Because a contaminated
control raises the cut, the marker will read as *negative* rather than obviously
broken: this warning may be the only sign anything is wrong.

What to do, in order:

1. **Look at the field.** `qc_marker_fields.csv` gives you `field_id`; the
   manifest maps it to a file. Open the image before deciding anything.
2. **If it is real, exclude the field, not the cells.** Artifacts are spatially
   clustered, and the field is already the replication unit, so dropping one
   costs a single replicate and removes the artifact wholesale. Never trim
   "extreme outlier cells": in a marker channel the extreme upper tail is where
   your real high expressers live, and trimming biases hardest in exactly the
   conditions with the strongest response.
3. **Record what you did.** Move the offending files into a subfolder, add that
   folder to `input.exclude_dirs`, and leave a note beside them saying what you
   measured and why. Nothing is deleted, the exclusion is one line, and the next
   person can undo it.
4. **Re-run and check the warning is gone.**

---

## 11. Troubleshooting

| symptom | cause | fix |
| --- | --- | --- |
| `check`: "N image(s) map to NO well" | Seq range wider than your `wells:` block, or wrong `points_per_well` | fix the mapping, or move genuinely unwanted images out of the folder |
| `check`: "well(s) X have images but no entry" | you imaged more wells than you described | add them to `wells:` |
| A marker reports 0 % everywhere | its gate is poisoned: check `gate_loo_spread_pct` and the control fields | exclude the bad control field |
| A marker reports ~100 % everywhere | control is a tight spike, P99 sits at background | per-marker `fpr: 0.001`, or read the graded bands instead |
| `graded_only` in `thresholds.csv` | no usable control and the data are not demonstrably bimodal | give it a control; the refusal is correct |
| `n_panels_pooled: 1` after setting `scope: plate` | the plate maps have different `plate_id`s | make them identical |
| Plate-layout figures missing | `predefined` mode has no plate positions | `well_numbering: plate`, if your keys really are wells |
| Conditions in the wrong order | the default is alphabetical | set `conditions.order:` |
| Trend figure empty | `trends.x` is not a numeric plate-map column | add a numeric companion column |

---

## 12. Exercises

Best done on a dataset that has already been analyzed, so you can check yourself.

1. **Read a plate map cold.** Take an experiment you did not set up. Without
   opening a config, write down: how many wells, how many fields each, which
   markers in which channels, which well is the control. Then run `check` and see
   if you were right.
2. **Break it on purpose.** Change `points_per_well` from 5 to 6 and run `check`.
   Read the error. Change it back.
3. **Find the control.** For one marker, open `control_assignment.csv` and
   `thresholds.csv` and write one sentence: "this gate is the Pxx of N cells from
   wells A, B, C." If you cannot, you do not yet know what the positive fractions
   mean.
4. **Move a gate.** Re-run `downstream --force-step thresholds` with `fpr: 0.001`
   for one marker and compare positive fractions. Decide which is more honest for
   that antibody, and be able to defend it.
5. **Add a numeric axis.** Take an experiment whose conditions have a natural
   order but no numeric column. Add one, set `trends.x`, re-run, and look at what
   the trend figure shows that the bar charts did not.
