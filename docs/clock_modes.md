# Clock rate: which mode answers which question

The molecular clock is the one parameter where the simulated and the real
analysis can be made to disagree *by configuration rather than by biology*, so
it needs deciding deliberately rather than per run. This document is the
decision record.

The short version: **run both modes, on different arms, and report the three
headline numbers rather than two** -- plus, on simulated data, the
decomposition in 1a and 1b, which says how much of the headline is evolution
and how much is the importation regime.

---

## 1. The quantities

| symbol | what it is | how it is obtained | comparable to real data? |
|---|---|---|---|
| **mu_truth** | the simulation's actual root-to-tip rate | regression of Hamming distance from the reference on sampling date, straight off the painted FASTA. No tree, no inference. | **yes** |
| **mu_sim** | TreeTime's estimate from the simulated sequences | `augur refine` with no `--clock-rate`, on ncov's own `tree_raw.nwk` | yes |
| **mu_real** | TreeTime's estimate from real sequences | the same unconstrained refine on a real-data build | — |

### 1a. Decomposing `mu_truth`

`mu_truth` is a *pooled* slope, and pooling hides a split that matters. The
painter accumulates mutations along time since each chain's own importation,
not along calendar date, and each chain starts at its seed's divergence rather
than at zero. So:

| symbol | what it is | how obtained | comparable to real data? |
|---|---|---|---|
| **mu_lineage** | the painter's per-lineage rate | within-chain slope, free intercept per chain | no -- real data has no chain labels |
| **mu_between_chains** | the importation regime's contribution | chain means vs chain mean dates, tip-weighted | no |

and the decomposition is exact:

    mu_truth  =  w * mu_lineage  +  (1 - w) * mu_between_chains

with `w` the share of date variance lying within chains. On the Virginia Delta
run `w` is small -- most date variance is *between* chains, because chains are
young relative to the wave -- so the pooled slope is dominated by the
importation term and sits far below the per-lineage rate. On synthetic data
with a known 8.0e-4 clock, `mu_lineage` recovers 8.03e-4 while the pooled fit
returns 2.99e-4, a 2.7x attenuation, and the identity reconstructs the pooled
number exactly.

**Use pooled `mu_truth` for the `mu_sim` comparison** (both sides are rooted
on Wuhan and inherit the same mixture) **and `mu_lineage` when the question is
what the painter generated.** They are not interchangeable and neither is
wrong.

### 1b. What the seeds themselves contribute

The seed schedule is built from real collection dates, and real Delta genomes
grew more diverged through 2021, so seed divergence is *correlated with import
date*. That correlation enters the pooled slope directly, weighted by how many
tips each chain contributes -- which means the founder lottery reshapes it.

| symbol | what it is |
|---|---|
| **mu_seed_trend** | seed divergence vs import date, every importation, unweighted. The real-world trend in what arrived. |
| **mu_seed_trend_surviving** | the same, restricted to chains that took off. |
| **mu_seed_trend_weighted** | the same, weighted by tips per chain -- the version that actually enters `mu_truth`. |

The gap between the first and the last is the founder lottery's effect. It is
the main reason `mu_truth` need not equal what real surveillance of the same
wave would measure, even though the seeds are real sequences: real surveillance
applies the same "only lineages that grew get sequenced" filter, but draws its
own lottery. Two consequences worth stating plainly:

- **Non-overlapping imports.** The seeds are sequenced *local* genomes, so
  descendants of introductions carrying some local evolution, not the
  introductions themselves. Their divergence-vs-date trend can differ from the
  true introductions'.
- **Chains that never take off.** A seed whose chain fizzles contributes no
  tips and so drops out of the weighted trend entirely.

Neither affects `mu_sim` vs `mu_truth`, which are computed from the same tips
and the same chains. Both affect `mu_truth` vs `mu_real_observed`.

---

A further number would be worth recording but is **not yet implemented**, and is
**not** comparable to real data:

| **mu_placed** | substitutions the painter actually placed, summed along each transmission path | painter accounting -- *not built* | **no** |

`mu_placed` would count substitution *events*. `mu_truth` counts observable
*differences*, so it undercounts back-mutations and repeat hits at the same
site. That undercounting is the point: a real root-to-tip plot undercounts in
exactly the same way, which is what makes `mu_truth` comparable and
`mu_placed` not. Their ratio would be a saturation/homoplasy diagnostic -- if
it drifted far from 1, the simulated sequences are saturating.

Implementing it needs the painter to account for the mutations it places per
record, which it does not currently do, so row 5 of the table below is not
runnable yet.

### Why root-to-tip and not something cleverer

It is what TempEst and TreeTime regress, and what the ~8e-4 subs/site/year
figure in the SARS-CoV-2 literature is derived from. Using anything else makes
the comparison to published rates an apples-to-oranges exercise.

Hamming-distance-to-reference is a tree-free form of root-to-tip distance. It
is valid here because the painter begins from `genetic_painter.reference_fasta`
and `refine.root` is that same sequence, so distance-from-reference *is*
distance-from-root. If the root ever stops being the ancestral sequence, this
equivalence breaks and the measure needs revisiting.

Units: slope in substitutions per genome per day, divided by the number of
**comparable sites** and multiplied by 365, to land on
substitutions/site/year.

Comparable sites, not genome length. `data/reference/reference.fasta` is now
byte-identical to ncov's copy, which masks the first 100 and last 50 bases to
`N` -- standard practice for unreliable genome termini. `hamming_to_reference`
skips any position that is not an unambiguous base on both sides, so
divergence is measured over 29,750 sites rather than 29,903. Dividing by the
raw length would deflate the rate by 0.51%. The result records both numbers.

---

## 2. The two modes

### `operational` -- the de facto prior (ncov's 0.0008 / 0.0004)

Treats the pipeline as an opaque measurement instrument and asks what the data
look like *after* standard surveillance distortion.

- Both real and simulated data get the **identical** prior, which is the only
  way the resulting trees share a time axis.
- This is what `cfg/nextstrain/ncov/base.yaml` does today, by not overriding
  ncov's defaults.
- Keeps `clock_filter_iqd: 4`.

Use it for: anything comparing trees, clades, traits or Auspice output between
real and simulated. **It cannot support a claim about the simulation's own
evolutionary rate** -- the rate was assumed, not measured.

### `generative` -- unconstrained, rate learned from the data

Asks whether the simulation's emergent evolution resembles reality.

- No `--clock-rate`; TreeTime fits the rate by root-to-tip regression.
- Wants **`clock_filter_iqd: 0`**. The filter measures deviation from the
  fitted or prior regression line, so with a mismatched prior it prunes valid
  branches -- and prunes different amounts from different arms, which
  manufactures a difference between sampling strategies that is really an
  artefact of filtering.

Use it for: validating the generator, and for reporting a rate.

### Why this is not a switch inside the ncov config

`ncov/workflow/snakemake_rules/main_workflow.smk` reads
`config["refine"]["clock_rate"]` with `[]` and passes `--clock-rate`
unconditionally. Omitting the key leaves ncov's default in force; emptying it
breaks augur. So `generative` is a **separate PhyloGAS step** over ncov's own
intermediates, not a different ncov config:

```
results/<build>/tree_raw.nwk              <- refine's input tree
results/<build>/filtered.fasta            <- refine's input alignment
results/<build>/metadata_adjusted.tsv.xz  <- refine's input metadata
```

An in-ncov approximation exists -- widening `clock_std_dev` (0.0004 -> ~0.01)
flattens the prior so the data dominate -- but it is an approximation and
should be checked against the separate step before being relied on.

---

## 3. Which mode, which arm, which comparison

The arm matters as much as the mode. Root-to-tip regression depends on which
tips are present, so a rate fitted on a *sampled* arm is partly a property of
the sampler.

| # | Question | Mode | Arm | Comparison |
|---|---|---|---|---|
| 1 | Can TreeTime recover a clock we know the answer to? | generative | `all_infections` | `mu_sim` vs `mu_truth` |
| 2 | Does our simulation evolve like SARS-CoV-2 does? | generative | `all_infections` | `mu_truth` vs `mu_real` |
| 3 | How much does each sampling strategy bias the rate estimate? | generative | every `strategy` arm | `mu_sim(arm)` vs `mu_truth` |
| 4 | Would our simulated outbreak look right on an Auspice dashboard? | operational | `strategy` arms | trees, clades, traits vs a real build |
| 5 | Are the simulated sequences saturating? | n/a | `all_infections` | `mu_placed` vs `mu_truth` — **not implemented** |
| 6 | What rate did the painter actually generate? | generative | `all_infections` | `mu_lineage` vs the literature ~8e-4 |
| 7 | How much of the wave's apparent rate is importation? | generative | `all_infections` | `w`, with `mu_lineage` vs `mu_between_chains` |
| 8 | Would real surveillance of this wave measure the same? | generative | `all_infections` | `mu_seed_trend` vs `mu_seed_trend_weighted` |

Rows 1 and 2 are the decomposition that matters. Comparing `mu_sim` to
`mu_real` directly -- the obvious two-way comparison -- conflates them: a
discrepancy could be the simulation's physics *or* TreeTime's inference error,
and you cannot tell which. `mu_truth` separates them, and is available only
because this is a digital twin.

Rows 6 to 8 come out of `truth` mode in one pass, so they cost nothing extra
once the chain join is available. Like `mu_truth` they are properties of the
simulation rather than of an arm -- the FASTA they read is the full painted
set, so `--arm` only labels the row.

**Row 3 is a benchmark result in its own right**, and belongs with the other
sampling metrics. It is a question an empirical study cannot ask, because real
surveillance has exactly one realised sampling scheme.

### The trap

Running the generative mode on a `strategy` arm and reporting the result as
"the simulation's clock rate" is wrong -- that number is the rate *as seen
through that sampling strategy*. On the `all_infections` arm the rate is a
property of the simulation. Same command, different claim.

### Why `--arm` does not currently change `mu_truth`

In `truth` mode `--arm` is only a label on the output row. `_benchmark_clock`
defaults its FASTA to `genetic_painter.output_prefix`.sequences.fasta -- the
*full painted set* -- so `mu_truth` is simulation-wide whatever arm it is filed
under, and the Snakefile's per-arm loop writes the same value once per arm.

That is the right number, but it is right by default rather than by design, and
the trap above makes it look as though the arm were doing the work. Two
consequences worth knowing before anyone "fixes" it:

- Pointing that default at an arm's subset would silently turn every `mu_truth`
  row into the sampler-dependent quantity this section warns against.
- `mu_sim` (`inferred`) *does* depend on the arm, because it re-runs refine on
  that arm's build directory. So one row of the CSV varies by arm and the other
  does not, which is correct but reads as an inconsistency.

Row 1 of the table wants `all_infections` for `mu_sim`, not for `mu_truth`. If
`nextstrain.builds` has no `all_infections` entry -- it is commented out by
default -- then `mu_sim` is only ever measured on sampled arms, and row 1's
comparison cannot be made. `mu_truth` is unaffected.

---

## 4. Applying this to real surveillance

The operational comparison is only valid if both sides get the identical
distortion, which is why `base.yaml` holds deviations only and both arms share
it. ### How little real metadata is actually needed

Worth being precise about, because real line lists are disclosive and the
answer is "much less than you would think".

| comparison | real metadata required |
|---|---|
| `mu_real_observed` vs `mu_truth` (tree-free, rows 1-2) | **`strain` and `date`. Nothing else.** |
| `mu_real` via TreeTime | the same two, plus an alignment |
| operational tree / trait comparison (row 4) | whatever geography you want traits on |

Root-to-tip regression is divergence against date, so geography, demographics
and county never enter it. Two columns, at day resolution -- which is what
GenBank and Cov-Spectrum publish openly. **The generative comparison, the one
that validates the simulation, needs no sensitive data at all.**

Better still, run it with the *same estimator on both sides*:

```bash
phylogas benchmark clock --mode truth --kind simulated --arm all_infections
phylogas benchmark clock --mode truth --kind real --arm va_2021_observed \
    --fasta <aligned_real.fasta> --metadata <strain_date.tsv>
```

That yields `mu_truth` and `mu_real_observed` from identical arithmetic, with
no TreeTime on either side, so the difference is attributable to the data
rather than split between the data and two separate inference runs. It is the
most directly comparable pair available.

The one requirement is that real sequences be **aligned to the same
reference** -- Hamming distance needs equal lengths. Use ncov's
`results/<build>/filtered.fasta` from a real-data build, or `augur align`.
Unaligned records are skipped and counted rather than silently compared.

### Where sensitive metadata does become necessary

Only the operational comparison (row 4), and only to the resolution you choose
to run traits at. Three things make that easier than it first appears:

1. **`county` is not required on a real arm.** The county requirement comes
   from `phylogas benchmark mugration`, which reads the traits JSON as
   `models['county']` -- and that benchmark needs ABM ground truth, so it
   cannot run on real data anyway. A real arm can legitimately run traits at
   `division`, or at a coarser aggregate.
2. **Coarsen rather than withhold.** RUCC category or health district is
   already in the pipeline via `population.rucc_file`, and aggregating county
   up to one of those keeps the comparison meaningful while reducing
   disclosure risk in small counties.
3. **The same context inputs.** `refine.root: Wuhan/Hu-1/2019` must be
   satisfiable on both sides, and `data/phylogas/_shared/` already carries
   that record, so both arms root on the identical sequence.

### Metric scoping on real data

Everything with `needs_truth = True` in BeyondBaseline's `eval_names.py`
requires the ABM's hidden graph and cannot run on real data. What can: the
clock estimates, the Auspice outputs, and -- given a real line list to compare
the sequenced subset against -- `kl_targets` and `equity_<stratifier>`. Those
last two are arguably the most directly useful thing this framework offers a
health department, and they are also the ones that need the sensitive file.

`mu_real` should come from a period and geography matched to the simulated
wave, through the same unconstrained step. Otherwise differences in sampling
window, which dominate root-to-tip fits over short spans, get attributed to
biology.

---

## 5. Running it

```bash
# all three quantities for one arm
phylogas benchmark clock --config config.yaml --arm all_infections --mode all

# just the tree-free truth rate: no ncov run needed
phylogas benchmark clock --config config.yaml --mode truth

# the sampling-bias question, across arms
phylogas benchmark clock --config config.yaml --arm surs --mode inferred
```

`truth` mode also emits the decomposition of 1a and 1b whenever it can join
chain membership. It needs `component_id` keyed on `alias_pid`, which lives in
TwinSampler's all-events table -- the painted metadata has no chain column,
because the painter writes `real_strain` only for index cases. The path is
taken from the configured ascertainment output, or named explicitly:

```bash
phylogas benchmark clock --config config.yaml --mode truth \
    --components results/<project>/02_simulated_linelists/linelist_allevents.csv.xz
```

Without it the pooled `mu_truth` is still produced and `mu_lineage` reports
why it is absent; `--no-decompose` suppresses the extra rows entirely. The
join holds one entry per infection, so budget roughly a gigabyte of memory for
a 5M-infection run, or cap the sequence pass with `--max-records`.
`--takeoff-min` sets how many tips a chain needs before it counts as having
taken off, for the surviving/weighted seed rows.

Config:

```yaml
nextstrain:
  clock:
    mode: both          # operational | generative | both
    # max_records: 200000   # sampled across the file, not truncated
```

Results accumulate in `{benchmark.outdir}/clock_estimates.csv`, one row per
(arm, quantity), with the slope, the per-site-per-year rate, R^2, the tip
count and the date span. The decomposition rows add the chain count, the
largest and median chain size, and `within_date_variance_share` -- the `w` of
1a, which is what tells you whether the pooled number is dominated by
evolution or by importation. R^2 and the span are reported because a rate fitted
over a narrow window is unreliable regardless of mode, and a short sampling
window is the usual reason a clock estimate misbehaves.

For that to mean anything, the span has to describe the data the fit rests
on. `mu_truth` therefore bounds the window at `genetic_painter.start_date`
and drops the ids in `cfg/nextstrain/reference_id.txt`. Without both, a
single record -- the reference, carrying Wuhan-Hu-1's real 2019-12-26
collection date -- reported a 537-day span for a 70-day study, making the one
field you would check to judge the fit describe a period the fit barely
covers. `--date-min`/`--date-max` override the window and `--keep-context`
the exclusion; neither default applies to `--kind real`, where clipping to a
simulation's start date would be meaningless.

Low R^2 over a short window is expected, not a fault. Over 70 days at ~8e-4
the clock contributes only ~5 substitutions of divergence, while tips sampled
the same day differ by a comparable amount because divergence tracks chain
length rather than calendar time. On synthetic data with a known 8.0e-4 rate,
a 70-day window recovers 7.87e-4 at R^2 = 0.16, and a 400-day window recovers
7.94e-4 at R^2 = 0.86. The slope stays well determined because n is large;
it is the variance explained that falls.

---

## 6. What is assumed, and what would invalidate it

- **The painter has no clock parameter.** `rate_limited` is per replication
  cycle (`mutation_rate_per_cycle = 3.40e-6`, Poisson draws per transmission
  with burst sizes 10-1000). The per-site-per-year rate is emergent from the
  painter's parameters *and* EpiHiper's transmission timing -- the generation
  interval does half the work and lives outside PhyloGAS. So `mu_truth` must be
  measured; it cannot be looked up. A future model with an explicit time-based
  rate would change this.
- **Distance-from-reference does NOT equal distance-from-root, and that is
  all right.** This entry used to claim the equivalence held "while the
  painter starts from the configured reference". It does not start there.
  `genetic_painter.reference_fasta` is used only to emit the reference as the
  first output record; the config says so at that key -- "the ancestral
  genomes come from seed_fasta, not from here". Every chain begins at a real
  seed genome already carrying its own divergence from Wuhan: measured at mean
  41.6, sd 4.7, range 14-57 substitutions across the 3,322 seeds of the
  Virginia Delta run.

  What rescues the comparison is not that premise but a different one.
  `refine.root` is `Wuhan/Hu-1/2019`, so TreeTime's root-to-tip distances
  decompose exactly as the Hamming measure does: a shared Wuhan-to-Delta stem
  (constant, no effect on slope), plus that chain's seed offset, plus
  accumulation since its import. Both sides inherit the same seed offsets and
  the same import mixture, so `mu_sim` vs `mu_truth` remains apples to apples.
  The estimator is sound; the stated reason was wrong.

  What the old premise did get right is that `mu_truth` cannot be read as the
  *simulation's* rate. It is the wave's rate under this importation regime.
  `mu_lineage` is the painter's own rate -- see section 2a.
- **`mu_real` needs a matched period.** Unmatched windows make the comparison
  meaningless.
- **`max_records` samples, it does not truncate.** It used to `break` after
  the first N usable records. The painted FASTA is written in tick order, so
  that returned the earliest infections and a compressed window: the same data
  gave 8.5443e-04 over a 69-day span capped at 20,000, and 2.7506e-04 over the
  full 299 days uncapped -- a 3x discrepancy that was entirely an artefact of
  the cap. It is now a seeded reservoir sample, so the retained records span
  the whole window. Any number produced with a cap before 2026-10-05 should be
  recomputed.

  Capping is also no longer needed for speed. `hamming_to_reference` summed a
  generator over ~30k characters per record; vectorised, a 5.35M-record set
  costs about 3 minutes instead of 160.
- **Tip dates are onset dates, not collection dates, and not detection
  dates.** Three dates exist, all from `base_date + (tick - start_tick)`. The
  painter filters `exit_state` on `genetic_painter.painted_prefix` ("E2"), so
  the painted metadata is dated at *exposure*. TwinSampler filters on
  `prefix_override` (`A`, `P`, `I`, `dM`, `hM`) and dates on that event, so the
  line list -- and therefore everything ncov sees, via the sample CSV and
  `metadata_adjusted.tsv.xz` -- is dated at entry into a clinical state. A
  detection or report date does not exist anywhere in the chain: ascertainment
  draws a Bernoulli to decide inclusion and never shifts a date, and ncov's
  `adjust_regional_meta.py` does not touch dates either.

  A real GISAID `date` is specimen collection, which lags onset. Nothing here
  models that lag, so simulated tips sit earlier in the natural history than
  real ones. For the clock this is close to a constant offset -- it moves the
  intercept, not the slope -- but anything reading absolute dates (epi curves,
  time to detection, Rt) inherits a systematic shift.

- **`mu_truth` and `mu_sim` are dated differently, and it does not matter.**
  `mu_truth` reads the painted metadata (exposure dates); the tree is built
  from the line list (onset dates). Measured over 1,138,109 line-list rows the
  gap is 5.01 days, SD 0.477, with 89.2% at exactly 5 days. Treating the
  difference as errors-in-variables, the slope attenuates by 0.057% on a 69-day
  window and 0.004% on a 252-day one -- against a 6.8% gap between the measured
  8.5443e-04 and the 0.0008 prior, four orders of magnitude too small to
  affect the comparison. Recorded because the asymmetry is real and someone
  will notice it; it needs no correction.

  The *shape* is worth knowing separately: 3 days occurs 35,960 times against
  4 days' 4,498, non-monotonic around the mode. That is distinct routes through
  the state machine, not a dwell-time distribution being sampled -- a
  structural fact about `disease.json` recoverable from the line list without
  parsing the model.

- **The latent period is near-deterministic, which makes dating easier here
  than on real data.** Delta's incubation period is roughly 4 days median with
  an SD of 2-3 days, lognormal; here the exposure-to-onset SD is 0.48 days. So
  TreeTime is being asked to date tips whose true dates carry almost none of
  the noise real tips carry. This does not threaten `mu_sim` vs `mu_truth`,
  which is internal, but it bounds what that recovery validates: it shows
  TreeTime can recover a known clock on an easy dating problem, not that it
  will on a real one.
