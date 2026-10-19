# What the ABM contributes to the twin's molecular clock

`clock_modes.md` defines the clock quantities and says which comparisons are
legitimate. This file answers a different question: when the emergent clock
moves, **how much of that is EpiHiper and how much is ours?** It matters
because only one of the two is a knob we hold, and because a reviewer asking
"isn't your mutation rate just an artefact of the contact network?" deserves a
number rather than a hand-wave.

Short answer: the ABM owns exactly one term, it is now measured rather than
assumed, and it was wrong by 1.7x in every projection made before 2026-10-12.

## The identity

    mu_truth  =  (mutation_rate_per_cycle x L x E[cycles])  /  GI  x  s_path
                 \_____________ lambda ________________/      \__/   \____/
                        placed subs per transmission           ABM  emergent

| term | value | who owns it | how we know |
|---|---|---|---|
| `mutation_rate_per_cycle` | 3.40e-6 | us (painter constant) | set in `rate_limited.py` |
| `L` | 29,903 | the reference genome | fixed |
| `E[cycles]` | 5.72 | us, via `initial_viral_load`, `peak_viral_load`, `rt_early_population_threshold`, burst draw | computed over the burst distribution |
| **`GI`** | **6.872 d** | **EpiHiper, entirely** | **measured, `benchmark sequence`** |
| `s_path` | 0.851 | emergent | measured, `mu_truth` / implied |

`lambda` is the only term the painter controls, `GI` the only term the ABM
controls, and `s_path` is a consequence of both.

## It closes, at three independent parameter settings

This is the check that makes the decomposition more than bookkeeping. Nothing
below is fitted; every term is measured or computed from first principles.

**lambda against one-step divergence.** `lambda` is a count of mutations the
painter *places*; parent->child Hamming distance is a count of differences
*observable* after one transmission. At this divergence the within-transmission
collision rate is ~1.1e-5 (see `clock_modes.md` §1c), so the two should agree:

| parameters | E[cycles] | lambda predicted | one-step observed | agreement |
|---|---|---|---|---|
| uniform 10-1000, init 10, peak 1e6 | 2.02 | 0.2054 | 0.2054 | **0.0%** |
| exponential(100), init 2, peak 1e11 | 5.72 | 0.5816 | 0.5834 | **0.3%** |
| lognormal(med 1e3, s 0.7), init 2, peak 1e11 | 3.37 | 0.3430 | 0.3428 | **0.06%** |

The third row is the strongest of the three: `lambda` was printed by
`burst_summary()` from the config *before* the run, and the painter then placed
that many mutations independently. Two code paths reading the same config and
agreeing to 0.06% is a different claim from a post-hoc fit.

**The full identity against `mu_truth`.** Same two settings, `GI` measured,
`s_path` measured:

| | lambda | / GI | x s_path | predicted | measured `mu_truth` | agreement |
|---|---|---|---|---|---|---|
| uniform | 0.2054 | 6.872 | 0.9080 | 9.90 | 9.90 | **0.0%** |
| exponential | 0.5816 | 6.872 | 0.8510 | 26.29 | 26.37 | **-0.3%** |
| lognormal | 0.3428 | 6.872 | 0.8853 | 16.12 | 16.12 | **0.0%** |

(subs/genome/year. Divide by `L` for per-site: 26.37 -> 8.82e-4.)

So the emergent clock is fully accounted for. There is no residual to attribute
to the contact network, the mixing structure, or anything else about the ABM
beyond the generation interval.

### The interval did not move when lambda did

The lognormal run is a controlled experiment that came for free. It repainted
the *same* ABM output with `lambda` 1.7x lower, and the measured interval came
back **byte-identical in every state**:

| state | GI, exponential run | GI, lognormal run |
|---|---|---|
| ga | 6.93 | 6.93 |
| ma | 6.88 | 6.88 |
| mn | 6.76 | 6.76 |
| va | 6.88 | 6.88 |
| wa | 6.91 | 6.91 |

Not "close" -- every digit. That is the sharpest available statement of the
claim this document makes: the interval is a property of the ABM output alone,
and nothing the painter does can touch it. A 1.7x change in the one term the
painter owns moved the other term not at all.

`s_path` also behaved as predicted. Lower `lambda` means fewer repeat hits per
site, so less should be lost along root-to-tip paths: predicted ~0.89 before
the run, measured **0.8853**, a 0.6% miss.

## What EpiHiper actually determines

**The generation interval, and only that.** Infector exposure to infectee
exposure. Measured over 23.7M transmission edges across five twins:

| state | edges | GI (d) | subs/transmission | implied subs/genome/yr |
|---|---|---|---|---|
| ga | 6,126,993 | 6.93 | 0.583 | 30.71 |
| ma | 3,861,499 | 6.88 | 0.582 | 30.87 |
| mn | 4,287,236 | 6.76 | 0.582 | 31.40 |
| va | 5,292,666 | 6.88 | 0.584 | 30.99 |
| wa | 4,166,779 | 6.91 | 0.586 | 30.98 |

Three things to read off that table:

1. **The interval is stable: 6.76-6.93 d, a 2.5% spread.** EpiHiper's
   progression model behaves identically across five different states,
   populations and seeding schedules.
2. **Epidemic size is irrelevant to it.** Edge counts vary 59% (ma to ga) while
   both per-transmission quantities barely move. Whatever differs between these
   twins, it is not changing how fast transmission chains turn over.
3. **The clock's cross-state variation is the interval's, not the painter's.**
   The painter contributes 0.7% spread, the ABM 2.5%, and the implied clock's
   2.2% tracks the latter -- mn has the shortest interval and the fastest clock,
   and 0.582/6.76*365 = 31.42 recovers its value from the denominator alone.

So: if the emergent clock ever shifts between twins or between runs, **look at
the generation interval first.** It is the only route by which the ABM can move
it, and a change there means the disease-progression parameterisation changed.

What EpiHiper does *not* determine: anything about which substitutions appear
or how many. Attack rate, epidemic size, county mixing and network homophily
are invisible to this identity.

### Where in the config the interval actually comes from

There is no generation-interval parameter. It is emergent from six `dwellTime`
fields and the per-state `infectivity` weights in `cfg/exp1/disease.json`,
identical in all ten age x strain combinations:

| transition | dwell | role | infectivity |
|---|---|---|---|
| `E->P` (p=0.65) | fixed 3 | latent, symptomatic branch | 0 |
| `E->A` (p=0.35) | normal(5, 1) | latent, asymptomatic branch | 0 |
| `P->I` | fixed 2 | presymptomatic infectious | 1.6 |
| `I->rM` (p=0.958) | discrete, mean 4.375 | symptomatic infectious | 1.6 |
| `A->R` | normal(5, 1) | asymptomatic infectious | 0.8 |

Mean latent period = 0.65*3 + 0.35*5 = **3.70 d**, which is 54% of the
interval. The infectious window is 6.24 d at infectivity 1.6 (symptomatic) or
5 d at 0.8 (asymptomatic). Computing the interval as latent plus mean time into
the infectious run, branches weighted by expected transmissions (infectivity x
duration), gives **6.31 d** against the 6.872 d measured. The ~0.5 d gap is
what integer ticks predict (`startTick 0`, `endTick 700`): with transmission
only on day boundaries the mean time into the window is (n+1)/2 rather than
n/2, putting the estimate at 6.81.

So the measurement is consistent with the disease model, and **the fields that
would move the clock are the latent block** -- the `fixed: 3`, the
`normal(5,1)`, and the 0.65/0.35 split.

Two consequences:

* **Strain 1 and strain 2 have identical timing** (6.311 vs 6.309 d). Delta's
  1.6x is pure hazard, not tempo, so variant composition cannot move the clock.
* **The one per-state calibrated knob cannot touch timing.** `calibration.csv`
  tunes `transmissibility` and nothing else -- a hazard scalar. Every twin
  shares the same `disease.json`, so the modelled timing is byte-identical
  across all five.

### How much a recalibration could move the clock

`transmissibility` varies 64% across the twins while the interval varies 2.5%,
and the two are inversely ordered (Spearman rho = **-0.975**):

| state | transmissibility | edges | GI (d) |
|---|---|---|---|
| ga | 0.1512 | 6,126,993 | 6.93 |
| wa | 0.1691 | 4,166,779 | 6.91 |
| va | 0.1898 | 5,292,666 | 6.88 |
| ma | 0.1977 | 3,861,499 | 6.88 |
| mn | 0.2772 | 4,287,236 | 6.76 |

The plausible mechanism is that a higher hazard lands the first successful
transmission earlier in the infectious window, compressing the *realised*
interval while the modelled timing is unchanged. Elasticity is **-0.029**: an
83% change in transmissibility moves the interval 2.5%, and since
`mu_truth` ~ 1/GI, the clock by 2.5% the other way. Real and directional, but
very inelastic.

**This correlation is not established as causal, and a rival explanation is
large enough to produce it on its own.** The interval varies by age group,
because the `I->hM` and `I->dM` branches shorten the infectious run for older
people:

| age group | GI_est (strain 2) |
|---|---|
| school-aged (5-17) | 6.389 |
| preschool / adult | 6.363 |
| older adult (50-64) | 6.300 |
| golden-aged (65+) | 6.127 |

That is a 0.26 d span from age mix alone, against a 0.17 d observed cross-state
span -- 1.5x larger than the effect being explained. A state whose infections
skew older gets a shorter interval for reasons unrelated to transmissibility,
and the five twins differ in population, network, age structure, seeding and
geography all at once. `Open` has how to separate them.

## Which measures are ABM-resilient

Ranked by how much an ABM change would move them:

| measure | ABM-resilient? | why |
|---|---|---|
| one-step divergence (subs/transmission) | **yes** | per-transmission; the interval does not enter. Equals `lambda`, a painter quantity. |
| `lambda`, E[cycles] | **yes, fully** | computed from painter config alone |
| `mu_truth` | **no** | has `GI` in its denominator; 1/GI scaling, directly |
| `mu_sim` | **no** | tracks `mu_truth` (0.3-6.5% below it) and inherits the same dependence |
| `s_path` | partly | depends on how many mutations accumulate per unit time, so on both |

The practical consequence: **report subs/transmission when you want a statement
about the mutational model, and `mu_truth` when you want one comparable to real
surveillance.** The literature's ~8e-4 subs/site/year is a root-to-tip estimate,
so `mu_truth` is the only one of these comparable to it.

## Where the twin currently sits

| arm | lambda | GI | `mu_truth` | vs real (~8.0e-4/site/yr) |
|---|---|---|---|---|
| uniform painting | 0.2054 | 6.872 | 3.31e-4 | 0.41x |
| exponential painting | 0.5816 | 6.872 | 8.86e-4 | 1.11x |
| **lognormal painting (current)** | **0.3428** | **6.872** | **5.42e-4** | **0.67x** |
| lognormal + exp2 ABM (not yet run) | 0.3428 | ~4.69 | ~7.9e-4 | ~0.99x |

The current 0.67x is a deliberate transit state, not a regression: `lambda` is
now literature-grounded (Bar-On burst size at the median) and the interval is
still 1.5x too long. The exp2 disease model shortens the interval; see
`design/exp2.md`. Fixing the interval first instead would have pushed the clock
to 1.6x and left per-transmission divergence untouched, so this ordering is the
only one that converges.

`mu_sim` (augur refine, no clock prior) sat 0.3-6.5% below `mu_truth` in all
five states at the exponential setting -- all low, never high.

### The cross-state spread is the estimator, not the twins

`mu_truth` spreads 6.5-8.4% across the five states, far more than any of its
inputs. The lognormal run settles where that comes from, because it changed
`lambda` by 1.7x and left the spread untouched:

| | exponential | lognormal |
|---|---|---|
| `lambda` spread | 0.7% | 0.9% |
| GI spread | 2.5% | 2.5% |
| implied clock (lambda/GI) spread | 2.3% | 3.1% |
| **`s_path` spread** | **8.8%** | **8.8%** |
| `mu_truth` spread | 8.4% | 6.5% |
| mean r2 of the `mu_truth` fits | 0.46 | 0.31 |

**`s_path`'s spread is identical to the decimal across a 1.7x change in
`lambda`.** A physical property of the twins had every opportunity to move and
did not. Meanwhile r2 *fell* to 0.31 (ga 0.425 down to wa 0.216), because
fewer substitutions per genome means less signal per tip. The spread tracks fit
quality, not saturation.

The practical consequence: **no change to the disease model or the painter will
tighten this.** They move `lambda` and GI, whose spreads are under 3%. The
6.5-8.8% lives in the pooled root-to-tip estimator, which
`scripts/clock_within_component.py` already diagnoses -- each chain carries its
seed's own divergence as a constant, and continued importation attenuates the
pooled slope. The within-chain decomposition is the only lever on it.

One caution for comparing arms: at lower `lambda`, `mu_truth` is a *weaker*
measurement (r2 0.31 vs 0.46). Per-state values deserve error bars rather than
being read as points, and exp1-vs-exp2 comparisons of `mu_truth` are comparing
estimates of different precision.

## A correction worth recording

The projection table in `rate_limited.py` predicted ~28 subs/genome/yr for
exponential(100). The realised value is 26.37, so the projection looks good.
**It was right for the wrong reasons.** Reverse-engineering it:

    projected:  lambda / 4 d       x 365 x 0.54   = 28.7
    actual:     lambda / 6.872 d   x 365 x 0.851  = 26.3

It assumed a 4-day generation interval (**short by 1.72x**) and an
observable/raw fraction of 0.54 (**low by 1.58x**). The two errors run in
opposite directions and cancel to 1.09x, which is why the number landed.

That cancellation is luck and will not survive the next parameter change. The
comment in `rate_limited.py` still carries both superseded assumptions and its
projection column should be recomputed with GI = 6.872 and s_path = 0.851
before anyone uses it to choose a parameter again.

## Open

* **`s_path` is not a tight constant.** Per-state range 10.4-17.9% loss against
  a 14.9% mean, and the loose term is `mu_truth` itself, whose root-to-tip
  regressions carry r2 of only 0.36-0.53. Treat ~15% as central and the spread
  as fit noise, not as real variation between twins. `ma` is the outlier on both
  its `mu_truth` (9.30e-4, highest) and its ratio (0.896, worst r2 at 0.37).
* **Serial interval is not measured.** The interval here is generation
  (exposure to exposure). Serial (onset to onset) would need the symptom
  columns and is what epidemiological estimates usually report, so the two
  should not be compared directly.
* **Is transmissibility -> interval causal?** The cross-state rho of -0.975 is
  confounded six ways, and the age-mix rival above is larger than the effect.
  Two ways to test it, cheapest first.
  *Free:* the all-events tables carry `age_group`. Compute each state's
  infection age mix, weight the per-age table above by it, and see whether that
  alone predicts the 6.76 -> 6.93 ordering. If it does, the transmissibility
  explanation is dead without running anything.
  *One extra run:* take va twice -- once at its calibrated 0.1898, once at mn's
  0.2772 -- holding population, network, seeds, `disease.json` and the RNG seed
  fixed. Two runs differing in one variable beat five differing in six. If the
  interval drops ~2.5%, any future transmissibility recalibration silently
  shifts the clock by a known elasticity; if it does not move, recalibration is
  clock-neutral.
  The interval needs only the ABM output -- no painting, trees or Nextstrain --
  so this costs one EpiHiper run plus seconds of arithmetic. Note the code gap:
  `generation_interval()` does not depend on the painter, but it is only
  reachable through `benchmark sequence`, which demands a painted FASTA because
  it also computes divergence. Running this experiment wants a subcommand that
  takes an all-events file alone.
* **The 1.108x overshoot is not explained.** The identity closes, so it is
  not a bookkeeping error -- `lambda` would need to come down ~10%, which
  means `mean_burst_size` up (E[cycles] down). Whether to chase that or
  accept 1.1x is a judgement call, and 100 is the biologically motivated
  value rather than a fitted one.

## Reproducing any of this

    phylogas benchmark sequence --painted <painted.fasta.xz> \
        --infections <linelist_allevents.csv.xz> --out sequence_divergence.csv

writes `gi_mean_days`, `mean_substitutions` and
`implied_subs_per_genome_per_year` together, from one pass over one file, so
numerator and denominator always describe the same run. `mu_truth` comes from
`05_benchmarks/clock_estimates.csv` -- read it with a csv-aware parser, since
the `method` column contains a comma and `awk -F,` shifts every later field.
