# Statistical methods

What each analysis step computes, why, what it assumes, and how to read the
numbers. The [README](../README.md) covers running the pipeline;
[REFERENCE.md](REFERENCE.md) lists every setting and file. The first section
explains every idea in plain words; the numbered sections give the exact recipe.

---

## The statistics in plain words

**Replicate (what counts as n).** A replicate is one independent repeat of the
comparison. A field with 800 cells is still one field; its cells grew, were
stained and were imaged together, so they are not 800 independent tries. Fields
of one well share a seeding and a staining reaction too. The pipeline uses cells
only to measure each field, and uses fields, wells or independent
differentiations as n (§1). Counting cells as n is called *pseudoreplication*: it
makes tiny, meaningless differences look "significant".

**Confidence interval (CI).** The error bar on a condition. A 95 % CI is a range
built so that, if the experiment were repeated many times, 95 % of such ranges
would contain the true value. A wide bar means the replicates disagree, or there
are few of them.

**p-value.** How surprising the data would be if the conditions truly did not
differ. p = 0.01 means a difference this large would turn up in about 1 of 100
experiments where nothing is going on. It is not the probability that the result
is wrong, and it says nothing about how large the effect is (that is the
difference and its CI).

**Permutation test (how every p-value here is made).** Write each replicate's
condition on a card, shuffle the cards, deal them back out and recompute the
difference between the conditions. Repeat thousands of times. If the condition
does not matter, the real labeling is just one more shuffle, so the real
difference should look ordinary among the shuffled ones. The p-value is the
fraction of shuffles giving a difference at least as large as the real one.
Nothing is assumed about bell curves or equal spreads.

**The smallest possible p.** Few replicates allow only a few shuffles. Three
wells against three can be shuffled 20 ways, so p can never go below 0.1 however
large the effect; five fields against five cannot go below 0.008. Every
comparison lists this floor (`smallest_possible_p`). More replicates, not more
cells, lower it.

**q-value (false discovery rate).** When many conditions or markers are tested,
some p-values come out small by luck. A q-value is a p-value adjusted for how
many tests were run (Benjamini-Hochberg): calling everything with q < 0.05
"different" means that on average at most 5 % of those calls are false.

**Bootstrap (how an error bar is made from a few replicates).** Draw the
replicates at random *with replacement* (some appear twice, some not at all),
recompute the mean, and repeat 2,000 times. The middle 95 % of those means is the
CI. For a dose trend the replicates are redrawn within each dose and the slope is
refitted each time.

**Positivity gate.** The intensity above which a cell is called positive: the
level that only 1 % of that marker's negative-control cells exceed. Without a
control, the pipeline looks for two clearly separate populations (bright and dim)
with a two-component mixture model; if it does not find them, the marker gets no
gate rather than a guess.

**Robust z and MAD (QC flags).** To find an odd field, each field is compared with
the median of its group in units of the typical spread (the median absolute
deviation, MAD). Medians are used because one bad field cannot drag them.

**Variance components.** Fields differ from each other, wells differ a bit more,
separate differentiations differ more still. The pipeline estimates how much of
the spread comes from each level, which says where extra effort pays: if most of
it is between differentiations, imaging more fields will not help.

**Power.** The chance that an experiment of a given size detects an effect of a
given size. From the measured variances the pipeline works out the smallest
difference 2, 3, 4 or 6 independent experiments could reliably detect.

---

## 0. From pixels to a value per cell

**Segmentation.** Cellpose-SAM, in 2-D on the DAPI maximum projection or in 3-D on
the z-stack; masks are size-filtered on a log scale to drop debris and merges.

**Per-cell intensity.** The median of each marker's pixels inside the nuclear
mask: robust to hot pixels and to bleed at the nuclear boundary. `dapi_norm` is
that median divided by the cell's own nuclear-stain median, which removes
nucleus-to-nucleus differences in DNA content and focus; it is only meaningful
when the nuclear channel is a DNA stain (`nuclear_stain` warns for a histone
label such as H2B, where dividing adds noise).

**Saturated cells.** A cell whose median pixel is at the detector ceiling in any
channel is clipped debris, not a bright cell. It is left out of every gate,
summary and test (`qc.drop_saturated_nuclei`) and listed in
`qc_saturated_nuclei.csv`.

## 1. The replicate

```
condition  ─►  (differentiation)  ─►  well  ─►  field  ─►  cell
```

Each field is first reduced to one number per marker: its positive fraction and
its median intensity (`fields.csv`). Then `stats.replicate: auto` chooses:

| Design | Replicate | What the report says |
| --- | --- | --- |
| `stats.experiment_column` names differentiations, ≥ 2 per condition | **experiment**: fields averaged per well, wells per experiment | conclusions generalize to new differentiations |
| every condition has ≥ 3 wells | **well**: fields averaged per well | conclusions generalize to new wells |
| some condition has 1 or 2 wells | **field** | the fields of one well are partly dependent, so intervals and p are optimistic; with one well per condition, well and condition are confounded |

The choice is written to `replicate_unit.json`, stated at the top of the report,
and carried in every table (`replicate_unit`, `n_replicates`).

## 2. Quality control

**Field QC.** Per field: cell count, median DAPI and median nuclear size, each
as a robust z within its condition; a field is flagged when |z| > `outlier_z`
(3.5), or when it has fewer than `min_cells_per_field` cells. Flags are reported;
`drop_outlier_fields: true` leaves flagged fields out.

**Marker QC.** Each field's marker median is compared with the other fields of
its own well (one staining reaction, so they should agree); a field is flagged
when it is beyond `marker_outlier_z` *and* at least `marker_min_rel_dev` (25 %)
from the well median. The two conditions together keep the flag on artifacts:
with 3 to 6 fields per well the MAD of a few agreeing numbers collapses, and a
field 3 % off would otherwise score a large z. Channel health from extraction
(saturated fraction, signal-to-background) rides along. Nothing here drops data.

**Uniformity.** `well_uniformity.csv`: field-to-field spread within each well;
`cv_warn` flags wells whose fields disagree.

## 3. Positivity gates

One gate per marker × scope (plate × panel by default), never per well. The gate
is the `1 − fpr` quantile (P99) of the marker's **own** negative-control cells
(`control_for` / `positive_for` make a well a control for some markers only).
Without a usable control (fewer than `min_control_cells`), a two-component
Gaussian mixture is fitted on log intensities and accepted only if it beats one
component by `min_bic_gain` in BIC, the two modes are `min_separation` standard
deviations apart, and the smaller holds at least `min_weight` of the cells.
Otherwise the marker is `graded_only`: no positive fraction anywhere.

**Gate stability.** Each gate is recomputed with each control field left out;
a spread above `loo_warn_pct` names the field that moves it (a contaminated
control field).

## 4. Graded bands and co-expression

Within each marker's positive population the cells are split into tertiles of
intensity (neg / low / mid / high), with edges per marker × scope.
`band_composition.csv` gives each condition's distribution over the bands and
`coexpression.csv` the 4 × 4 band co-occurrence for every marker pair of a panel,
as a fraction of cells and as the mean ± SD over fields.

## 5. Condition summaries

`conditions.csv`, per condition × marker: the mean of the replicate values
(§1), SD, SEM, and a 95 % percentile bootstrap CI (2,000 resamples of the
replicates; no CI with fewer than `min_fields` replicates).

## 6. Comparisons against each marker's control

`comparisons_positive_fraction.csv`, `comparisons_median.csv`. Each condition
against the control of *that marker* (or `stats.reference_condition`):

* replicate values at the chosen unit (§1); a side with fewer than two replicates
  steps down a level (experiment → well → field) and the row says
  `pseudo_replicated`;
* difference of means with a Welch 95 % CI, Hedges' g and Cohen's d;
* a two-sided permutation p on the difference of means: every relabeling is
  enumerated when there are at most `n_permutations` of them (`exact`), else
  `n_permutations` random ones; `smallest_possible_p` is listed;
* Benjamini-Hochberg q within each marker.

**Factor effects** (`trends` step, `factor_effects_*.csv`): for each two-level
plate-map column (± a ligand, ± an inhibitor), the mean difference it makes
across every background it appears in, tested the same way.

## 7. Effects, trends, variance components and power

The `effects` step takes every per-field number (each marker's positive fraction
and median) and answers four questions the same way for each.

**How much varies at each level?** Nested method-of-moments variance components
pooled over conditions: fields within a well (`sd_field`), wells beyond their
fields (`sd_well`), differentiations beyond their wells (`sd_experiment`). A
level that is not replicated (one well per condition, one differentiation) is not
estimable and reported as NaN, never as zero.

**How large is the effect?** Each condition against a reference (the lowest level
of the dose / time column within its group, or `effects.reference`): difference
of replicate means, Welch 95 % CI, Hedges' g, permutation p.

**Does it follow dose or time?** For numeric condition columns with at least
`min_levels` levels: the least-squares slope of replicate values against x, or
against log10 x when the positive levels span ten-fold or more (a zero level is
then left out of the fit and the table says so); a 95 % CI from a bootstrap of
replicates within each level; p from permuting the level labels; Spearman ρ.
Columns in `exponential_columns` also get a decay time constant τ with a CI.

**How many experiments are needed?** With E experiments, W wells per experiment
and F fields per well, the variance of a condition mean is
σ_E²/E + σ_W²/(E·W) + σ_F²/(E·W·F); the minimum detectable difference is
(t₁₋α/₂ + t_power)·√(2·var) with df = 2(E − 1). When σ_E cannot be estimated it is
taken as 0 and the row says the MDD is a lower bound;
`power.assumed_experiment_sd` adds an explicit assumption.

## 8. Signaling readouts

Signaling markers (β-catenin, SMAD2/3, SMAD4, SMAD1, … recognized by name) are
read by where the protein is, not by how much of it there is:

* **nuclear**: the median of a nucleus core (the mask eroded by a pixel) minus the
  image background;
* **N:C**: that, divided by the background-subtracted median of a cytoplasmic
  ring grown 1 to 4 µm outside the nucleus (split between neighbors, restricted
  to cell-body pixels that clear background + `body_k`·MAD).

The background is the median of cell-free pixels in a band just beyond the cells,
else of all cell-free pixels, else the camera offset (a field packed with cells
has no free pixels; its noise is then taken from the dimmest pixels). Values are
divided by the median of the reference wells of the same plate × panel
(`signaling.reference_values`, by default the control wells), so 1 = the
reference; with no reference they stay in their own units. Signaling markers get
no positivity gate unless `signaling.gate: true`. The ring QC pictures
(`figures/signaling/`) show where the cytoplasm was measured.

## 9. Assumptions and limitations

* **Exchangeability, not normality.** The permutation tests assume only that
  replicates are exchangeable under the null, which random assignment of
  conditions to wells provides.
* **Power comes from replicates.** Thousands of cells make each field precise;
  they cannot replace wells and differentiations.
* **One gate per scope** assumes the cells under it were stained and imaged alike.
* **Projections mix planes.** In a maximum projection, signal above and below a
  nucleus is assigned to it; this lowers N:C. Segment in 3-D when absolute N:C
  matters.
* **One differentiation is one replicate of the batch.** Between-experiment
  variance is then not estimable and power estimates are lower bounds.

## References

- Good P. *Permutation, Parametric, and Bootstrap Tests of Hypotheses*, 3rd ed.
  Springer, 2005.
- Efron B, Tibshirani RJ. *An Introduction to the Bootstrap.* Chapman & Hall, 1993.
- Benjamini Y, Hochberg Y. Controlling the false discovery rate. *J R Stat Soc B*
  1995;57:289-300.
- Lazic SE. The problem of pseudoreplication in neuroscientific studies.
  *BMC Neurosci* 2010;11:5.
- Pachitariu M, Rariden M, Stringer C. Cellpose-SAM: superhuman generalization
  for cellular segmentation. *bioRxiv* 2025.
