# Confocal image analysis: with a point-and-click set-up

This pipeline takes a folder of confocal images (`.nd2`, `.tif`, `.czi`, …),
finds every nucleus, measures each marker in each nucleus, and turns the
numbers into thresholds, positive fractions, co-expression, QC and figures, all
collected in one `report.html`.

You describe your experiment once in a **GUI that runs on your laptop**
(Windows, macOS or Linux). The GPU-heavy segmentation then runs on the **HPC**,
and the analysis runs either on the HPC or on your laptop. You don't need to
write any code or YAML.

```
 your laptop (GUI)                      the HPC
 ─────────────────                      ───────
 1 pick the experiment folder
 2 fill in the plate map     ──┐
 3 choose the analysis         ├──►  files saved next to your data  ──►  GPU job: segmentation
 4 check everything            │     (plate_map.yaml, config.yaml,        CPU job: analysis + figures
 5 HPC settings              ──┘      hpc_job.env)                            │
                                                                            ▼
                                                  analysis_2d/report.html in your experiment folder
```

---

## Contents

1. [Quick start](#1-quick-start)
2. [Installing and starting the GUI](#2-installing-and-starting-the-gui)
3. [The five steps of the GUI](#3-the-five-steps-of-the-gui)
4. [Running on the HPC](#4-running-on-the-hpc)
5. [Running on only some images](#5-running-on-only-some-images)
6. [Where everything lives](#6-where-everything-lives)
7. [Reading the results](#7-reading-the-results)
8. [Troubleshooting](#8-troubleshooting)
9. [Command line (for people who want it)](#9-command-line-for-people-who-want-it)
10. [For maintainers](#10-for-maintainers)

---

## 1. Quick start

**One-time:** install Python on your laptop ([§2](#2-installing-and-starting-the-gui)).
The HPC environment builds itself on the first job ([§4.1](#41-the-python-environment-on-the-hpc)).

**Every experiment:**

1. Put the images on the lab share, one folder per acquisition, e.g.
   `KR/exp_033/zmax_stain1/`.
2. Double-click **`launch_gui.bat`** (Windows) or **`launch_gui.command`** (Mac).
   The GUI opens in your browser.
3. Work through the five steps in the left sidebar. Press **💾 Save** on each.
4. Press **🔎 Run check** in step 4 and fix anything it reports.
5. Log in to the HPC and paste the commands shown in step 5. The segmentation
   job starts, and the analysis job starts by itself when segmentation is done.
6. Open `exp_033/analysis_2d/report.html` (`analysis_3d` for 3-D results).

---

## 2. Installing and starting the GUI

### You need

- **Python 3.10 or newer.** Get it from <https://www.python.org/downloads/>.
  On Windows, tick **"Add python.exe to PATH"** in the installer.
- The lab share mounted on your laptop (the drive or folder where you see `KR/`).

### Starting it

| Computer    | Do this                                                                                                                                                                                       |
| ----------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Windows** | Double-click `launch_gui.bat` in this folder.                                                                                                                                                 |
| **macOS**   | Double-click `launch_gui.command`. If macOS refuses ("unidentified developer"), right-click → **Open**, or open Terminal and type `bash ` then drag the file into the window and press Enter. |
| **Linux**   | `bash /path/to/analysis_pipeline_v2/launch_gui.sh`                                                                                                                                            |

The **first start takes a few minutes**: it builds a private Python environment
in your home folder (`~/.confocal_pipeline_gui_venv`), so nothing else on your
computer is touched. After that it starts in seconds. It updates itself when the
pipeline's requirements change.

The GUI opens as a page in your web browser (`http://localhost:8501`). It runs
only on your computer and nobody else can reach it. **To stop it, close the
black terminal window.**

---

## 3. The five steps of the GUI

Every step saves an ordinary text file next to your data. You can close the GUI
at any time. When you open the experiment again, everything is reloaded: the GUI
reopens the last experiment by itself and every page starts from the files on
disk.

**Saving keeps your file.** Saving over an existing `config.yaml` or
`plate_map.yaml` asks first (the old version is kept as `.bak`) and warns louder
if the file was edited outside the GUI since it was opened. The file is saved
*into* the one on disk: an unchanged setting keeps its line exactly, a changed one
is replaced on its own line, and every comment typed by hand stays. The page lists
which settings a save will change.

### Step 1 · Experiment

Choose the **experiment folder**: the folder that holds either the images
themselves or one sub-folder per acquisition. The GUI lists every folder that
contains images, with how many images it holds and which acquisition mode the
filenames suggest.

- **use**: untick folders you don't want analyzed.
- **name**: how the folder is labeled in tables and figures.
- **segmentation**: `2d` segments the DAPI maximum projection (works for any
  image). `3d` segments the full z-stack (z-stack folders only, needs more GPU
  time).

### Step 2 · Plate maps: what is in each well

This is the only information the pipeline can't work out for itself. There is
one plate map per acquisition folder.

1. **Plate and acquisition.** Choose the plate type (8, 18 or 96 wells). The
   *acquisition mode* is detected from the filenames:
   
   | Filenames look like | Mode              | Wells are…                                                                |
   | ------------------- | ----------------- | ------------------------------------------------------------------------- |
   | `…_Seq0042.nd2`     | Predefined points | numbered in the **order you imaged them**: well 1 = the first well imaged |
   | `Well07_….nd2`      | Generate points   | the well number in the filename                                           |
   | `WellC10_….nd2`     | Spinning disk     | the well label in the filename (C10)                                      |
   
   **Plate ID:** if several folders are different stain panels of the **same
   physical plate**, give them the **same** plate ID.

2. **Predefined points only: which images belong to which well?** The wells
   table starts empty: add one row per well **in the order you imaged them** and
   type its number of images in the *images* column, or use *Add N wells with K
   images each*. *Same number of images in every well* and *explicit Seq ranges*
   are the alternatives. Nothing is pre-labeled as a control.

3. **Stain sets.** One column per stain set, one row per channel (channel 0 is
   the first channel saved). Type the **antibody/marker name** (SOX2), not the
   fluorophore (FITC). Add a second set only if some wells carry a different
   panel.

4. **Wells and conditions.** One row per well. Add a column for everything that
   differs between wells (`treatment`, `dose_uM`, `day`, …). Numbers stay
   numbers, so a dose or time column can later be a plot axis.
   
   - **class:** mark your **negative-control** wells as `control`. Positivity
     thresholds are computed from these wells, so this is the most important
     column.
   - **stain:** only for wells that use a stain set other than `default`.
   - **positive_for / control_for:** for the case where a control well is
     *positive* for one marker. Example: a pluripotent well is the negative
     control for BRA but SOX2-high. Write `SOX2` under *positive_for* for that
     well, and `SOX2` under *control_for* for the well that is SOX2-negative.
   - **✏️ Fill many wells at once:** type wells (`1-6, 9` or `C1-C4`), pick a
     column and a value, then press **Fill**.

5. **The plate picture** colors the wells by any column, so mistakes are easy
   to spot.

6. **Check and save.** The GUI loads your plate map with the pipeline's own code
   and maps every image to a well. Red messages must be fixed before running.
   Yellow messages are worth reading. **💾 Save plate map** writes
   `plate_map.yaml` into that folder and keeps the previous version as
   `plate_map.yaml.bak`.

### Step 3 · Analysis settings

Sensible defaults are already filled in. The **Essentials** tab covers what
most experiments need:

- **Channels:** which channel is DAPI and which channels to measure.
- **Which images:** leave empty for all, or test on a few. See [§5](#5-running-on-only-some-images).
- **Conditions:**
  - *Condition = combination of*: which plate-map columns define a condition,
    e.g. `treatment`, or `treatment` + `dose_uM`.
  - *One figure per value of*: splits a crowded crossed design into one figure
    per treatment.
  - *Control wells are marked in column … with value …*: normally `class` = `control`.
- **Markers** (one row per marker):
  - *analyze*: untick to leave a marker out of everything. It works by name, so
    E-CAD in some wells and FOXA2 in others, both on channel 1, can be kept or
    dropped separately (`downstream.exclude_markers`).
  - *readout*: *standard* (intensity and % positive) or, for signaling markers
    (β-catenin, SMAD2/3, SMAD4, … recognized by name), *signaling: nuclear* or
    *signaling: N:C*: read by where the protein is, relative to reference wells,
    never gated.
  - *intensity*: *÷ nuclear stain* for nuclear markers with a DNA stain; *raw* for
    membrane or cytoplasmic markers, and for every marker when the nuclear channel
    is not a DNA stain (H2B: the GUI warns).
  - *one gate across*: one threshold for a marker stained with the same antibody
    and dilution in several panels of one plate.
- **Trends:** choose a numeric column (dose, time, day) for dose-response and
  time-course plots.
- **What counts as a replicate:** leave on *auto* (wells when every condition has
  3 or more wells, else fields with a caveat in the report; independent
  differentiations when a column names them). See [§7](#7-reading-the-results).
- **Analyses to run:** everything by default. The core always runs (QC,
  thresholds, field / well / condition summaries and comparisons, the report);
  untick what the experiment does not need: *Distributions*, *Positive fractions
  and co-expression*, *Graded bands*, *Two-marker views*, *Plate layout*, *Dose /
  time courses*, *Effects, dose trends and power* (`downstream.skip_steps`). A
  skipped analysis's earlier results move to `_skipped/` in the analysis folder.
- **Segmentation results used, and the analysis folder:** what is on disk for each
  folder (*265 of 265 images in 2D*), with a warning when it is not the mode the
  folder is set to. *Analyze*: whatever exists, 2-D only, or 3-D only. The folder
  is named after what is used: `analysis_2d`, `analysis_3d` or `analysis_2d_3d`.

The other tabs hold every remaining option. Hover the ⓘ next to any of them for
an explanation. **💾 Save config.yaml** writes the file into the experiment
folder, with every option's explanation kept as a comment.

### Step 4 · Check & run

- **🔎 Run check** reads the config, the plate maps and the filenames. It opens
  no images and writes nothing. It reports images that map to no well, wells
  without a condition, missing stain sets, and markers without a negative
  control. **Run it before every HPC submission.**
- **Run downstream analysis** runs the analysis on your laptop once segmentation
  results exist. It takes minutes, with no GPU needed.
- *Segmentation on this computer* is only for quick tests. It needs Cellpose
  installed and ideally a GPU.

### Step 5 · HPC job

This step writes `hpc_job.env` (the job settings) into the experiment folder and
shows the commands to type on the HPC: `cd` to the experiment and one `bsub` line
per segmentation job, each with a note saying what it does, plus a table of
re-run commands ([§4.3](#43-watching-and-re-running)). Fill in, once:

- **Share on this computer / same share on the HPC.** For example, `Z:\` on
  your laptop is `/data/Warmflash` on the HPC. The GUI uses this pair to
  translate your folders into HPC paths, and remembers it.
- **Python environment on the HPC:** the folder from [§4.1](#41-the-python-environment-on-the-hpc).

Then choose which steps, folders and images the jobs should run.

---

## 4. Running on the HPC

There are three job files in `hpc/`. **They are generic: you never edit them.**

| Job file                 | Runs                                                   | Node |
| ------------------------ | ------------------------------------------------------ | ---- |
| `run_pipeline_2dseg.lsf` | load → segment (2-D, on the DAPI projection) → extract | GPU  |
| `run_pipeline_3dseg.lsf` | load → segment (3-D, on the z-stack) → extract         | GPU  |
| `run_downstream.lsf`     | check → manifest → downstream analysis + report        | CPU  |

### How a job finds your experiment

You submit **from your experiment folder**. LSF remembers that folder, and the
job reads `hpc_job.env` from it. No path is ever typed into a job file, so a
job can't run on the wrong data because of a stale path. The settings stay with
the experiment they describe: re-running an experiment months later needs
nothing but `cd` and `bsub`.

```
exp_033/
├── config.yaml        ← what to analyze          (GUI step 3)
├── hpc_job.env        ← how to run it on the HPC  (GUI step 5)
├── zmax_stain1/
│   ├── plate_map.yaml ← what is in each well     (GUI step 2)
│   └── *.nd2
└── analysis_2d/       ← results appear here (analysis_3d/ for 3-D results)
```

### 4.1 The Python environment on the HPC

Nothing to do by hand: every job checks its Python environment (**Python
environment on the HPC** in GUI step 5) before it starts, which takes a few
seconds, and installs only what is missing or wrong. The first job builds the
environment (several minutes); later jobs start at once.

The check (`hpc/check_env.py`) reads every package version against
`requirements.txt` and, for PyTorch, which CUDA the wheel was built for and
which GPU architectures it contains. The gpu-b6000 nodes are Blackwell GPUs
(sm_120) with CUDA 12.8: only the **cu128** build fits. An older build (cu126)
imports fine but has no sm_120 kernels; a newer one (the plain PyPI wheel is now
CUDA 13.0) needs a newer GPU driver. When PyTorch is wrong it is uninstalled and
force-reinstalled from the cu128 index; other packages are installed from
`requirements.txt` without touching PyTorch. Each segmentation job then runs a
real computation on the GPU and stops with a clear message if that fails.

If an environment is broken in a way the check does not repair, tick *Rebuild
the environment* in GUI step 5 (`FRESH_VENV="true"`): the next job deletes it and
builds it from scratch. Untick it afterwards. To build or check the
environment ahead of time on the login node:

```bash
bash /data/Warmflash/KR/analysis_pipeline_v2/hpc/setup_env.sh /scratch/$USER/confocal_v2_venv
```

> The paths in this README assume the lab share appears as `/data/Warmflash` on
> the HPC and that this pipeline lives in `KR/analysis_pipeline_v2`. If yours
> differ, set them in GUI step 5.

### 4.2 Every experiment

GUI step 5 prints these commands with your paths filled in:

```bash
cd '/data/Warmflash/KR/exp_033'                                            # the experiment folder
bsub < /data/Warmflash/KR/analysis_pipeline_v2/hpc/run_pipeline_2dseg.lsf    # segmentation on a GPU
```

That is all. When segmentation has **finished successfully**, the job submits the
analysis job (`run_downstream.lsf`, CPU) by itself. With a 2-D and a 3-D folder
there is one `bsub` line per job; whichever finishes last starts the analysis, and
a job that fails starts nothing (resubmit it and it will). `THEN_ANALYSIS="false"`
in `hpc_job.env` (a tick box in step 5) turns this off.

### 4.3 Watching and re-running

- `bjobs` lists your jobs; `bpeek <jobid>` shows a running job's output.
- Logs are written into the experiment folder: `lsf_seg2d.<jobid>.log`,
  `lsf_downstream.<jobid>.log`, and the pipeline's own `pipeline.log`.
- **Resubmitting is safe.** Every image and every analysis step is checkpointed,
  so a resubmitted job skips what is finished. A job that hit its time limit
  just continues where it stopped.
- **Re-running a part.** `hpc/run.sh`, typed from the experiment folder, submits
  the same job files with a one-off setting; `hpc_job.env` is not changed
  (`$PIPE` = `/data/Warmflash/KR/analysis_pipeline_v2`):

  | When                                               | Type (after `cd` to the experiment) |
  | -------------------------------------------------- | ----------------------------------- |
  | changed an analysis setting or a plate map         | `bsub < $PIPE/hpc/run_downstream.lsf` |
  | redo one analysis step anyway (and all after it)   | `bash $PIPE/hpc/run.sh redo thresholds` |
  | redo the whole analysis                            | `bash $PIPE/hpc/run.sh redo-all`    |
  | changed segmentation settings, or masks look wrong | `bash $PIPE/hpc/run.sh reseg`       |
  | try everything on a few images first               | `SEQ="0-4" bash $PIPE/hpc/run.sh`   |
  | only the 3-D segmentation                          | `SEG=3d bash $PIPE/hpc/run.sh seg`  |
  | see what a job would do, no GPU, nothing submitted | `bash $PIPE/hpc/run.sh test`        |
  | where are my jobs                                  | `bash $PIPE/hpc/run.sh status`      |

  `bash $PIPE/hpc/run.sh help` lists everything. A changed setting re-runs what
  depends on it by itself; new masks make that image's measurements stale, so
  they are redone too.
- **More time or memory:** add flags to `bsub` instead of editing the job file,
  e.g. `bsub -W 24:00 -M 96000 < …/run_pipeline_3dseg.lsf`.
- **Queue for the downstream job:** it asks for no queue, so the cluster's
  default CPU queue is used. Add `-q <queue>` if needed.

### 4.4 Testing settings without spending GPU time

```bash
cd '/data/Warmflash/KR/exp_033'
bash /data/Warmflash/KR/analysis_pipeline_v2/hpc/run.sh test
```

This prints the exact commands every job would run and stops. It works on the
login node (the job scripts are bash-4.2 safe); there, a torch import error in
the environment check is expected, since the login node is older than the
compute nodes.

### 4.5 Editing `hpc_job.env` by hand

It is a plain text file with one `NAME="value"` per line and an explanation
above each. A blank template is in `templates/hpc_job.env`. Copy it into the
experiment folder and edit it. Files edited on Windows are fine: the jobs
strip Windows line endings.

---

## 5. Running on only some images

By default every image is used. To test settings on a few fields first, or to
re-do part of a plate, select images by:

- **Seq number**, as written in the filename (`Seq0042` → `42`): `0-9, 15, 20-24`
- **well**: `1-3` (numbered wells) or `C1-C4, D5` (spinning disk). In
  predefined-points mode, wells are looked up through the plate map.

Giving both means an image must match both.

| Where                                                                   | Applies to                                   |
| ----------------------------------------------------------------------- | -------------------------------------------- |
| GUI step 3 → *Which images* (saved as `input.select` in `config.yaml`)  | everything: segmentation, check and analysis |
| GUI step 5 → *Which images* (saved as `SEQ` / `WELLS` in `hpc_job.env`) | the HPC jobs only; overrides the config      |
| Command line `--seq` / `--wells`                                        | that one command                             |

The analysis uses the same subset, so a test on five fields gives an analysis of
exactly those five fields. Clear the selection to go back to all images. The
analysis then rebuilds its tables automatically.

---

## 6. Where everything lives

```
analysis_pipeline_v2/              ← the code: one copy for the whole lab
├── launch_gui.bat / .command / .sh   start the GUI (Windows / macOS / Linux)
├── gui/                           the GUI (app.py) and its file writers (core.py)
├── hpc/                           generic job files, run.sh (re-runs), check_env.py, setup_env.sh
├── templates/                     config.yaml and hpc_job.env templates (fully commented)
├── main.py                        the command-line entry point
├── src/                           the pipeline itself
├── CHANGELOG.md · CONTRIBUTING.md  what changed in each version; how to change the code
├── docs/REFERENCE.md              full reference manual (every option, every output)
├── docs/TUTORIAL.md               the ideas behind plate maps and configs, with recipes
├── docs/STATISTICAL_METHODS.md    every statistical method, in plain words first
├── examples/                      adding your own analysis step
└── tests/                         fast tests, no image data needed
```

**Nothing is ever written into the code folder.** Everything belonging to an
experiment (settings, logs, results) is written into the experiment folder.
All paths inside `config.yaml` are relative to the experiment folder, so the
same files work from your laptop and from the HPC, and an experiment folder can
be moved or renamed without breaking anything.

Each image gets a results folder next to it (`img.nd2` → `img/`) with the
nuclear masks, a segmentation overlay to check by eye, and the per-nucleus
measurements. That folder holds **one** segmentation, 2-D or 3-D
(`segmentation_metadata.json` says which). To have both, keep MIPs and z-stacks
in two acquisition folders, one set to 2d and one to 3d. Results are always
labeled by the mode they were really made in.

The analysis folder is `analysis_2d/`, `analysis_3d/` or `analysis_2d_3d/`
(`downstream.results_dir: auto`), named after the segmentation results used.

---

## 7. Reading the results

The header of **`report.html`** and `pipeline_version.json` name the pipeline version
and commit that made the results.

Open **`report.html`** in the analysis folder. It is one scrollable page with
every figure, and the warnings come first, including what the replicate was. The
tables behind it are in `tables/` (`conditions.csv`, `wells.csv`,
`thresholds.csv`, `comparisons_*.csv`, `effects_*.csv`, …), and the per-nucleus
data are in `cells.parquet` / `intensities.parquet`.

Things to know when reading it:

- **n is never cells.** Each field is reduced to one value per marker. With 3 or
  more wells per condition the **well** is the replicate (fields averaged per
  well); with fewer, the field is, and the report says what that means (with one
  well per condition, well and condition are confounded). Name independent
  differentiations in `stats.experiment_column` to make them n.
- **Tests are permutation tests** on replicate means, exact when the design is
  small, and every row lists the smallest p the design can give (3 wells vs 3
  wells cannot go below 0.1). Plain-language explanations of every method are in
  [`docs/STATISTICAL_METHODS.md`](docs/STATISTICAL_METHODS.md).
- **Signaling markers** (β-catenin, SMADs, …) are shown relative to their
  reference wells (1 = the reference) and get no % positive.
- **Effects** (`figures/effects/`): each condition vs the reference with a CI,
  dose / time slopes with a CI, how variance splits between fields, wells and
  differentiations, and how many experiments a given effect would need.
- **Saturated nuclei** (median at the detector ceiling) are left out and listed
  in `tables/qc_saturated_nuclei.csv`.
- Every marker figure states its normalization underneath (raw, ÷ nuclear stain,
  or a signaling unit); a histone nuclear label (H2B) used as a denominator is
  warned about.
- A marker without a usable negative control and without two clear populations
  gets **no threshold** (*graded only*) rather than a made-up one. The reason is
  stated at the top of the report.

Every table and figure is described in [`docs/REFERENCE.md`](docs/REFERENCE.md).

---

## 8. Troubleshooting

| Problem                                         | What to do                                                                                                                                                                      |
| ----------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| *"Python 3.10 or newer is needed"*              | Install Python from python.org (Windows: tick "Add to PATH"), then start the launcher again.                                                                                    |
| The browser page does not open                  | Open <http://localhost:8501> yourself. If a second GUI is already running, it may be on 8502.                                                                                   |
| The GUI finds no images                         | Check that you picked the experiment folder, not a parent of it. Other image types can be added under step 3 → *Images & segmentation* → `file_types`.                          |
| *"N images map to no well"*                     | The Seq → well mapping does not cover every Seq number. Check *images per well* in step 2; the *Image → well table* shows what went where.                                      |
| *"wells … have images but no row"*              | Add those wells to the wells table (or press *Re-list wells from the images*).                                                                                                  |
| *"No well is marked as a negative control"*     | Set `class` = `control` on your control wells.                                                                                                                                  |
| Job log says *"no hpc_job.env in …"*            | Submit from the experiment folder (`cd` there first), and save step 5 in the GUI.                                                                                               |
| Job log says *"NEEDS FIXING: torch …"* | Nothing to do: the job reinstalls PyTorch from the cu128 index itself ([§4.1](#41-the-python-environment-on-the-hpc)). |
| Job log says *"still not right after reinstalling"* or *"the GPU test failed"* | Tick *Rebuild the environment* in GUI step 5, save, resubmit ([§4.1](#41-the-python-environment-on-the-hpc)). |
| *"No GPU visible to PyTorch"*                   | The segmentation job must go to a GPU queue. The job file asks for `gpu-b6000`; override with `bsub -q <gpu queue>` if needed.                                                  |
| A threshold looks wrong                         | Look at `thresholds.csv` and the *gate stability* warnings in the report. Usually a control well is mislabeled, or a control is positive for that marker (use *positive_for*). |

Still stuck? `pipeline.log` in the experiment folder has the full history,
including error tracebacks. Send it along when asking for help.

---

## 9. Command line (for people who want it)

**Running on your own computer.** Two scripts do what the HPC jobs do, without
LSF. Both take the experiment folder (the folder holding `config.yaml`) and
use the `cellpose_sam` conda environment unless you set `CONDA_ENV`:

```bash
# segmentation, 2-D (GPU if PyTorch sees one); test on a few images first
SEQ="0-4" DATASETS="zmax" bash run_local_2dseg.sh /path/to/exp_040
DATASETS="zmax" bash run_local_2dseg.sh /path/to/exp_040

# downstream analysis (CPU): manifest + every analysis step
DATASETS="zmax" bash run_local_downstream.sh /path/to/exp_040
```

Each script lists its options (`SEQ`, `WELLS`, `DATASETS`, `FORCE`, ...) in its
header. `DATASETS` takes the dataset names from `config.yaml`; leave out a
dataset that has not been segmented. GUI step 4 offers the same choice.

The GUI only writes files. The pipeline itself is the same command-line tool as
before, so everything can still be scripted:

```bash
PIPE=/data/Warmflash/KR/analysis_pipeline_v2
EXP=/data/Warmflash/KR/exp_033

python $PIPE/main.py check      --config $EXP/config.yaml
python $PIPE/main.py process    --config $EXP/config.yaml --seq 0-9     # load + segment + extract
python $PIPE/main.py manifest   --config $EXP/config.yaml
python $PIPE/main.py downstream --config $EXP/config.yaml
python $PIPE/main.py init-plate-map --config $EXP/config.yaml           # text-file plate-map starter
```

All options, outputs and the design behind them are in
[`docs/REFERENCE.md`](docs/REFERENCE.md). [`docs/TUTORIAL.md`](docs/TUTORIAL.md)
walks through plate maps and five common experiment layouts.

---

## 10. For maintainers

**Version control.** The code is a git repository; `main` is what the lab runs.
Changes are made on a branch in a separate clone, tested, written into
[`CHANGELOG.md`](CHANGELOG.md), merged, and released with a version number in
`VERSION` and a `vX.Y.Z` tag. [`CONTRIBUTING.md`](CONTRIBUTING.md) has the
workflow, the release rules, and how to rerun an old result with the code that
made it. `python main.py --version` prints the version and commit in use; every
results folder records it in `pipeline_version.json` and in the report header.

**Tests** (synthetic data only, a few seconds each):

```bash
python tests/test_experiment.py      # plate maps, Seq → well, stains
python tests/test_downstream.py      # thresholds, gating, statistics, palette rule
python tests/test_selection.py       # --seq / --wells / input.select
python tests/test_gui_core.py        # everything the GUI writes reads back exactly
python tests/test_pipeline_options.py  # end to end on a synthetic plate with known answers:
                                       # replicate unit, permutation p floor, signaling N:C,
                                       # saturated nuclei, excluded markers, skipped steps,
                                       # effects, GUI writers, HPC one-off settings (~1 min)
```

`tests/synthetic_fields.py` writes a 3-dose × 3-well plate as the image pipeline
would have left it (masks, channel TIFFs, joblibs), so `check`, `manifest` and
`downstream` run without Cellpose.

**How the GUI stays in step with the pipeline.** The GUI does not re-implement
any rule:

- folders are scanned with `experiment.template`;
- plate maps are validated by loading them with the pipeline's own `PlateMap`;
- the final word is always `main.py check`;
- `config.yaml` is written by filling `templates/config.yaml`, so a new option
  added to that template appears in the GUI's *Advanced* tabs with its comment
  as help text, without any GUI change.

**Changing cluster defaults** (GPU queue, memory, wall time): edit the `#BSUB`
lines at the top of the job files in `hpc/`. Modules and the CUDA build are set
in `templates/hpc_job.env` (`CUDA_MODULE`, `TORCH_INDEX`); the environment check
and repair are `hpc/check_env.py` and `ensure_env` in `hpc/common.sh`.

**Shared with `../micropattern_pipeline`** (2026-11): `src/image_processing/`,
`src/experiment/`, `src/selection.py`, `src/pipeline_config.py`,
`src/seg_status.py`, the signaling measurement (`src/downstream/signaling/`),
the effects module (colony → field), `gui/core.py`, the HPC job files and
`hpc/run.sh`. A fix in one belongs in the other.

**Relation to the previous version** (`../analysis_pipeline`):

- The analysis code is the same.
- Added: image subsets (`--seq`/`--wells`/`input.select`), the GUI, generic
  HPC jobs, and launchers.
- Existing `config.yaml` and `plate_map.yaml` files work unchanged and can be
  opened in the GUI.
