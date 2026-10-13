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

## It closes, at two independent parameter settings

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

**The full identity against `mu_truth`.** Same two settings, `GI` measured,
`s_path` measured:

| | lambda | / GI | x s_path | predicted | measured `mu_truth` | agreement |
|---|---|---|---|---|---|---|
| old | 0.2054 | 6.872 | 0.9080 | 9.90 | 9.90 | **0.0%** |
| new | 0.5816 | 6.872 | 0.8510 | 26.29 | 26.37 | **-0.3%** |

(subs/genome/year. Divide by `L` for per-site: 26.37 -> 8.82e-4.)

So the emergent clock is fully accounted for. There is no residual to attribute
to the contact network, the mixing structure, or anything else about the ABM
beyond the generation interval.

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

| | value | vs real SARS-CoV-2 (~8.0e-4/site/yr) |
|---|---|---|
| `mu_truth` mean | 8.864e-4 /site/yr (26.37 subs/genome/yr) | **1.108x** |
| per-state range | 8.56e-4 - 9.30e-4 | 1.07x - 1.16x |
| before recalibration | 3.31e-4 (9.90 subs/genome/yr) | 0.41x |

`mu_sim` (augur refine, no clock prior) sits 0.3-6.5% below `mu_truth` in all
five states -- all low, never high.

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
