# Changelog

All notable changes to analysis_pipeline_v2. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the version
numbers follow [semantic versioning](https://semver.org) (see CONTRIBUTING.md).
Each entry says when a change alters numbers in existing results.

## [Unreleased]

## [2.0.0] - 2026-11-12

First version under version control. Version 2 of the confocal field pipeline,
including the 2026-11 changes listed in `KR/TODO_item8_analysis_pipeline_v2_proposals.md`.

### Pipeline
- Per-image steps `load`, `segment` (Cellpose-SAM, 2-D or 3-D), `extract`, each
  checkpointed; new masks make the extracted intensities stale.
- Downstream steps `assemble`, `qc`, `thresholds`, `gating`, `aggregate`,
  `distributions`, `composition`, `joint`, `plate`, `trends`, `effects`,
  `signaling`, `report`, with config-aware checkpoints.
- Positivity gates from each marker's own negative control (P99), or a checked
  two-component mixture; graded bands and co-expression.

### Added 2026-11
- Replicate unit chosen from the design (`stats.replicate`: field, well or
  experiment) and stated in every table and the report.
- Exact or Monte-Carlo permutation tests with their smallest possible p, Welch
  CI and Hedges' g; factor effects use the same tests.
- `effects` step: variance components, effects, dose trends, power.
- Signaling readouts (nuclear, N:C) relative to the reference wells; never gated
  unless asked.
- Saturated-nucleus filter, histone guard, normalization captions on figures.
- `exclude_markers` by name, `downstream.skip_steps`, `results_dir: auto`
  (`analysis_2d` / `analysis_3d`), comment-keeping GUI saves, `hpc/run.sh` with
  automatic downstream submission and `MP_<NAME>` one-off overrides.
- Plain-language `docs/STATISTICAL_METHODS.md`.

### Changed
- p-values now come from permutation tests on replicate means; with few
  replicates they are larger than before (and say their floor). Gates, positive
  fractions and medians are unchanged (checked on the 20261027 plate).
- Doses sort as numbers on every axis.
- Field QC compares a field only with fields of the same dataset, segmentation
  mode and condition.

### Version tracking
- `VERSION` file; `python main.py --version`; every run logs the version and git
  commit; downstream writes `pipeline_version.json` and shows the version in the
  report header.

### Fixed
- Report navigation links were printed as literal HTML.
- `check` assumed different default control classes than the analysis.
- Signaling background in fields with no cell-free pixels took its noise from
  the whole image (same fix as micropattern_pipeline 1.0.0).
