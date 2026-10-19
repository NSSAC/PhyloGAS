# EpiHiper Simulation Design of exp2

Adapted from [exp1](exp1.md). Three files differ from `cfg/exp1/` --
`disease.json`, `runParameters.json` and `config.json`; every other file is
byte-identical.

Three changes: two to transition dwell times in `disease.json`, and the run
length. Nothing else about the design moves: the cross-infection cell
structure, the 60% transmissibility advantage of the new variant, vaccination,
seeding, networks and replicates are all inherited from exp1 unchanged.

## 1. Generation interval, now variant-specific

**Why.** The realised generation interval of the exp1 twins was measured at
**6.87 days** (mean over 23.7M transmission edges, five states, spread 2.5%);
see `docs/abm_contribution.md`. Published estimates put the generation interval
near **4.7 d for Delta** and **5.5 d for Alpha** (Hart et al., Lancet Infect
Dis 2022), so exp1 ran roughly 1.3-1.5x slow. Because the emergent molecular
clock is `mutations-per-transmission / generation interval`, a long interval
drags the clock down and forces compensating error into the mutational model.

exp1 also gave the two variants **identical timing** (6.311 vs 6.309 d): the
1.6x advantage of the new variant is pure hazard, not tempo. Real Alpha and
Delta differ in tempo as well, so exp2 separates them.

| transition | exp1 (both variants) | exp2 variant 1 (Alpha) | exp2 variant 2 (Delta) |
|---|---|---|---|
| `E->P` (p=0.65) | fixed 3 | fixed 2 | fixed 2 |
| `E->A` (p=0.35) | normal(5, 1) | normal(3.5, 1) | normal(2.5, 1) |
| `P->I` | fixed 2 | fixed 2 | fixed 2 |
| `I->rM` (p=0.958) | discrete, mean 4.375 | discrete, mean 3.607 | discrete, mean 2.459 |
| `A->R` | normal(5, 1) | normal(5.0, 1) | normal(3.5, 1) |
| **generation interval** | **6.81 d** | **5.46 d** | **4.69 d** |

Notes on the parameterisation:

* **`fixed` dwells are kept integer.** Ticks are whole days -- `exposure_tick`
  and `exposure_date` agree to the day with zero residual across 5.3M rows --
  so a fractional `fixed` value would be ambiguous. Fractional means are
  carried by the `normal` and `discrete` forms, which are sampled and rounded.
* **`I->rM` keeps its declining shape**, refitted as a truncated geometric over
  days 1-10 rather than replaced by a point value, so the "most clear early, a
  few linger" structure of exp1 survives.
* **Infectivity values are untouched** (1.0/1.0/0.5 for variant 1, 1.6/1.6/0.8
  for variant 2), so the 60% transmissibility advantage is exactly as in exp1.
* **The asymptomatic latent period no longer exceeds the symptomatic one.** In
  exp1 `E->A` (5 d) ran 2 d longer than `E->P` (3 d), which has no clear
  biological reading -- asymptomatic cases are not infected later, they simply
  never develop symptoms -- and contributed ~0.6 d to the interval.
* **The latent block alone cannot reach the target.** It is 54% of the interval
  in exp1, but hitting 4.7 d through latency alone would need a mean latent
  period near 1.3 d, well under the 2-3 d reported for SARS-CoV-2. The
  infectious window has to shorten too, which is why `I->rM` and `A->R` move.

Predicted generation intervals come from the analytic model in
`docs/abm_contribution.md` (latent plus mean time into the infectious run,
branches weighted by infectivity x duration, plus 0.5 d for integer ticks).
That model predicted **6.86 d** for exp1 against **6.87 d** measured, so it is
trustworthy at this scale. **The exp2 figures above are predictions and must be
re-measured** with `phylogas benchmark sequence` once the runs exist.

## 2. Immunity waning, reshaped at constant mean

**This is a deliberate departure from exp1's documented design**, which
specifies "W is sampled from an exponential distribution of mean 6 months".

**Why.** The exponential puts its mode at zero, so a large share of recoveries
wane almost immediately: 15.1% within 30 days, 28.0% within 60. In the exp1
twins, 25.3% of people with a variant-2 infection between ticks 128-428 had
more than one in that window, with observed gaps as short as 39 days.
Reinfection at three weeks is not supported by SARS-CoV-2 observation;
reinfection at 4-6 months is.

`R2 -> W` is the **only** route by which a variant-2 recovery returns to the
susceptible pool -- `R2` has exactly one exit and the vaccination intervention
does not target it -- so this distribution alone gates every repeat variant-2
infection.

    exp1:  "dwellTime": { "gamma": { "alpha": 1.0, "beta": 183.0 } }
    exp2:  "dwellTime": { "gamma": { "alpha": 4.0, "beta": 45.75 } }

The mean is held at 183 days (`beta = 183 / alpha`); only the shape moves.
Applies to all **15** such transitions: `R1_*->W_*`, `R2_*->W_*`, `V_*->W_*`
across the five age strata.

| waned within | exp1 (alpha=1) | exp2 (alpha=4) |
|---|---|---|
| 30 d | 15.1% | 0.5% |
| 60 d | 28.0% | 4.4% |
| 90 d | 38.8% | 13.7% |
| 150 d | 55.9% | 41.5% |
| still immune > 1 yr | 13.6% | 4.3% |

**What this does and does not buy.** It fixes the *timing* of repeat
infections, not the eventual count: over the 300-tick window the aggregate
repeat fraction only falls from ~25.3% to roughly 19%. And beyond ~250 days the
sign reverses -- exp2 produces *more* repeats than exp1, because the
exponential's heavy right tail (13.6% still immune at a year) is removed. Both
tails of the exp1 distribution are hard to defend; exp2 fixes the one that
matters for a one-year wave.

**Side effect to expect.** Less early waning means a smaller susceptible pool
during the wave, which should shorten the realised generation interval slightly
and raise the emergent clock. Direction is predictable, magnitude is not.

## 3. Run length cut from 700 to 500 ticks

`cfg/exp2/runParameters.json` `endTick` 700 -> **500**, and
`cfg/exp2/config.json` `projectionEndDate` / `simulationEndDate`
2022-10-31 -> **2022-04-14** to match (`tickZero` 2020-11-30 + 500 d). The
date window and the tick count are coupled because the `%tick%` templates
(`template-seeding-*.json`, `template-vax.json`) are expanded by the pipeline
over that window.

**Why 500 is safe.** Nothing in the analysis reads past tick ~430. The painted
window is set by `genetic_painter.start_tick: "auto"` -- the painted variant's
first importation, read from `abm.seed_schedule` -- plus
`genetic_painter.num_ticks: 300`. Measured first variant-2 seed per state:

| state | v2 first seed | painted window | headroom to 500 |
|---|---|---|---|
| ga | 130 | 130-430 | 70 |
| ma | 103 | 103-403 | 97 |
| mn | 124 | 124-424 | 76 |
| va | 128 | 128-428 | 72 |
| wa | 124 | 124-424 | 76 |

The tree, the clock (`mu_truth`, `mu_lineage`), the line lists and the sampling
calendar all live inside that window. The root-to-tip lever arm is set by
`num_ticks: 300`, not by `endTick`, so the cut does not touch `Sxx` and the
cross-state `mu_truth` spread is unaffected.

**What exp1 was discarding.** Binning exposures in
`data/example_data/0.5_crossinfection/output.csv.gz` (ga, 89.5M rows, ticks
0-699): alpha is extinct by tick ~300, and Delta does not die out -- it
plateaus at ~1.02M exposures per 50 ticks from tick 400 through 699. Ticks
500-699 hold 4.08M of 11.58M Delta exposures, 35% of the total, none of which
exp1 ever painted. The last 270 ticks were simulated and thrown away.

**Truncated schedules, and why it is harmless.** `wa-2.csv` seeds to tick 693
and `va-2.csv` to 569; vaccination events likewise extend past 500. All fire
after the painted window closes, so they cannot causally affect anything inside
it.

**Why not 450.** The binding state is ga at 430. 500 leaves 70 ticks of slack,
which is what absorbs a later variant-2 importation in the exp2 calibration or
any future increase in `num_ticks`. At 450 a shifted window would truncate
silently and surface as a short alignment rather than an error.

**Replicates are not a `cfg/` lever.** No file in `cfg/exp1` or `cfg/exp2`
mentions replicates, and nothing in the repo expands the `%tick%` templates
either. The count (exp1: 60 per cell) is set in the EpiHiper pipeline batch
layer on Rivanna, so a reduction to 30 happens there, not here.

## What exp2 deliberately does NOT change

* **The cross-infection sweep.** `R1_*` susceptibility stays at **1.00**, which
  is the `alpha = 1.00` cell of exp1's four-level sweep
  (alpha in {0.25, 0.50, 0.75, 1.00}, one cell per level per state). Pinning it
  to a single intermediate value would collapse a swept dimension into a
  constant. exp2 inherits the cell structure and should be run across the same
  four levels.
* Transmissibility calibration, vaccination, seeding, networks, replicates.

**Recalibration.** Shortening the generation interval raises `R_0` for a given
transmissibility, so the exp1 calibration (current variant against `R_0 = 4.0`)
does not carry over unchanged. Transmissibility should be re-fitted before the
exp2 twins are treated as comparable to exp1.

## Resolved: the exp1 duration discrepancy

`design/exp1.md` previously stated a duration of 365 days while
`cfg/exp1/runParameters.json` set `endTick: 700`. The output file settles it --
exp1 ran ticks 0-699. `design/exp1.md` has been corrected; 365 was never what
ran.
