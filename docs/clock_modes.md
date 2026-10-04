# Clock rate: which mode answers which question

The molecular clock is the one parameter where the simulated and the real
analysis can be made to disagree *by configuration rather than by biology*, so
it needs deciding deliberately rather than per run. This document is the
decision record.

The short version: **run both modes, on different arms, and report three
numbers rather than two.**

---

## 1. The three quantities

| symbol | what it is | how it is obtained | comparable to real data? |
|---|---|---|---|
| **mu_truth** | the simulation's actual root-to-tip rate | regression of Hamming distance from the reference on sampling date, straight off the painted FASTA. No tree, no inference. | **yes** |
| **mu_sim** | TreeTime's estimate from the simulated sequences | `augur refine` with no `--clock-rate`, on ncov's own `tree_raw.nwk` | yes |
| **mu_real** | TreeTime's estimate from real sequences | the same unconstrained refine on a real-data build | — |

A fourth number would be worth recording but is **not yet implemented**, and is
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

Rows 1 and 2 are the decomposition that matters. Comparing `mu_sim` to
`mu_real` directly -- the obvious two-way comparison -- conflates them: a
discrepancy could be the simulation's physics *or* TreeTime's inference error,
and you cannot tell which. `mu_truth` separates them, and is available only
because this is a digital twin.

**Row 3 is a benchmark result in its own right**, and belongs with the other
sampling metrics. It is a question an empirical study cannot ask, because real
surveillance has exactly one realised sampling scheme.

### The trap

Running the generative mode on a `strategy` arm and reporting the result as
"the simulation's clock rate" is wrong -- that number is the rate *as seen
through that sampling strategy*. On the `all_infections` arm the rate is a
property of the simulation. Same command, different claim.

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

Config:

```yaml
nextstrain:
  clock:
    mode: both          # operational | generative | both
```

Results accumulate in `{benchmark.outdir}/clock_estimates.csv`, one row per
(arm, quantity), with the slope, the per-site-per-year rate, R^2, the tip
count and the date span. R^2 and the span are reported because a rate fitted
over a narrow window is unreliable regardless of mode, and a short sampling
window is the usual reason a clock estimate misbehaves.

---

## 6. What is assumed, and what would invalidate it

- **The painter has no clock parameter.** `rate_limited` is per replication
  cycle (`mutation_rate_per_cycle = 3.40e-6`, Poisson draws per transmission
  with burst sizes 10-1000). The per-site-per-year rate is emergent from the
  painter's parameters *and* EpiHiper's transmission timing -- the generation
  interval does half the work and lives outside PhyloGAS. So `mu_truth` must be
  measured; it cannot be looked up. A future model with an explicit time-based
  rate would change this.
- **Distance-from-reference equals distance-from-root.** True while the painter
  starts from the configured reference and ncov roots on the same sequence.
- **`mu_real` needs a matched period.** Unmatched windows make the comparison
  meaningless.
