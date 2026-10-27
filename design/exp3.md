# EpiHiper Simulation Design of exp3

Adapted from [exp2](exp2.md) and [exp1](exp1.md). Only **one file** differs from
`cfg/exp2/`: `calibration.csv`, which reverts to the `cfg/exp1/calibration.csv`
per-state `transmissibility` values. Every other file in `cfg/exp3/`
(`disease.json`, `runParameters.json`, `config.json`, `trait.json`,
`initialization-base.json`, `intervention-base-noise.json`, and the seeding/vax
templates) is byte-identical to `cfg/exp2/`.

Simulation path on Rivanna: `/project/bii_nssac/epihiper-simulations/pipeline-jc/run/20261026_1/output_root/proj/20261026_1/batch_1/`. The output files are `output*.csv.gz`.

## 1. What `exp2` Fixed (and `exp3` Preserves Unchanged)

Evaluated across all five states (`GA, MA, MN, VA, WA`, cell `0.25`) in the
completed `exp2` runs (`/project/bii_nssac/epihiper-simulations/pipeline-jc/run/20261022_2/output_root/proj/20261022_2/batch_1/`),
Changes 1–3 of [exp2](exp2.md) worked as intended:

1. **Variant-specific Generation Intervals (`disease.json` dwell times):**
   * **Variant 1 (`E1`, Alpha):** Realised generation interval measured at
     **4.92–4.96 d** across states (mean **4.94 d**, vs. 5.46 d analytic
     prediction and 6.87 d in `exp1`).
   * **Variant 2 (`E2`, Delta):** Realised generation interval measured at
     **4.56–4.59 d** across states (mean **4.58 d**, 0.7% cross-state spread,
     vs. 4.69 d analytic prediction and 4.7 d in Hart et al. 2022).
   * At $\lambda = 0.3428\text{ subs/tx}$, $\text{GI} = 4.58\text{ d}$ yields a
     raw substitution clock of $27.3\text{ subs/genome/yr}$
     ($9.14 \times 10^{-4}\text{ subs/site/yr}$) and $\mu_{\text{truth}} \approx
     8.0 \times 10^{-4}\text{ subs/site/yr}$ (**1.00× real SARS-CoV-2**).
2. **Reshaped Gamma Immunity Waning (`alpha = 4.0, beta = 45.75`, mean 183 d, mode 137 d):**
   * Repeat-infection intervals $\le 30\text{ d}$ fell from **15.1% (`exp1`) to
     0.02%–0.03% (`exp2`)**; intervals $\le 60\text{ d}$ fell from **28.0% to
     0.65%–0.75%**.
   * Maximum `E2` infections in a single host fell from **10 (`exp1`) to 5 (`exp2`)**,
     and windowed `E2` reinfection fractions dropped to **22.1%–27.8%**.
3. **500-tick horizon (`runParameters.json`, `config.json`):**
   * Retains full coverage of the 300-tick painted window plus 70–97 ticks of
     headroom.

---

## 2. What Broke in `exp2`: Calendar Phase Shift Against Empirical Importation Schedules

In `exp2`, `calibration.csv` was re-fitted upward against a target naive-population
$R_0(\text{v1}) = 4.0$ (`commit 92dd1b6`). While the first-order infectious
duration reduction of Variant 1 in `exp2` required only a **1.107×** (+10.7%)
increase in `transmissibility` to hold per-case transmissions constant, the
empirical `exp2` calibration raised `transmissibility` by **1.33× to 2.11×**:

| State | `exp1` `transmissibility` (`exp3`) | `exp2` `transmissibility` | `exp2 / exp1` Multiplier |
| :--- | :---: | :---: | :---: |
| **ga** | `0.15117165154352094` | `0.317152439012276` | **2.10×** |
| **ma** | `0.19768989966250955` | `0.364169845830287` | **1.84×** |
| **mn** | `0.2771846159090046`  | `0.368448550109758` | **1.33×** |
| **va** | `0.1898064352024446`  | `0.330233099449603` | **1.74×** |
| **wa** | `0.1691184588837505`  | `0.35715540839982`  | **2.11×** |

Combining a **28%–33% shorter Generation Interval** with a **1.33×–2.11× higher
`transmissibility`** nearly doubled the real-time exponential growth rate
$r \approx \frac{\ln R_e}{\text{GI}}$ (`~0.14 d⁻¹` in `exp1` $\to$ `~0.24 d⁻¹`
in `exp2`), shifting both the `E1` (Alpha) and `E2` (Delta) waves **30 to 66 days
earlier** than the real 2021 calendar (`tickZero = 2020-11-30`):

| State | Variant | `exp1` `Re(50k)` | `exp2` `Re(50k)` | `exp1` Peak Tick (Calendar) | `exp2` Peak Tick (Calendar) | Shift |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **GA** | `E1` | 2.43 | 2.98 | tick 125 (`2021-04-04`) | tick 59 (`2021-01-28`) | **-66 d** |
| **MA** | `E1` | 2.60 | 3.32 | tick 109 (`2021-03-19`) | tick 88 (`2021-02-26`) | **-21 d** |
| **MN** | `E1` | 2.67 | 3.11 | tick 115 (`2021-03-25`) | tick 94 (`2021-03-04`) | **-21 d** |
| **VA** | `E1` | 2.78 | 3.06 | tick 122 (`2021-04-01`) | tick 102 (`2021-03-12`) | **-20 d** |
| **WA** | `E1` | 2.59 | 3.00 | tick 122 (`2021-04-01`) | tick 95 (`2021-03-05`) | **-27 d** |
| **GA** | `E2` | 1.88 | 2.12 | tick 254 (`2021-08-11`) | tick 207 (`2021-06-25`) | **-47 d** |
| **MA** | `E2` | 2.07 | 1.94 | tick 230 (`2021-07-18`) | tick 190 (`2021-06-08`) | **-40 d** |
| **MN** | `E2` | 1.61 | 1.43 | tick 262 (`2021-08-19`) | tick 216 (`2021-07-04`) | **-46 d** |
| **VA** | `E2` | 1.83 | 2.11 | tick 247 (`2021-08-04`) | tick 217 (`2021-07-05`) | **-30 d** |
| **WA** | `E2` | 1.70 | 1.81 | tick 239 (`2021-07-27`) | tick 200 (`2021-06-18`) | **-39 d** |

### Consequences of the Calendar Shift in `exp2`

Because the importation schedules (`data/seed_schedule/<state>-{1,2}.csv`) and
vaccination schedules (`data/vax_schedule/`) are anchored to **fixed real-world
2021 calendar dates**, accelerating the epidemic clock created two structural
distortions:

1. **Loss of `E1`–`E2` Co-Circulation Overlap:**
   In the real US record (and in `exp1`), Alpha (`E1`) peaked in early April
   2021 (`ticks 110–130`) and declined through May–June (`ticks 150–200`),
   directly overlapping the March–June take-off of Delta (`E2`). In `exp2`,
   `E1` burns through 65%–92% of the population by March and reaches complete
   extinction (`0 cases/day`) before `E2` begins its main rise.
2. **Severe Founder Bottleneck on `E2` Importations:**
   In `data/seed_schedule/<state>-2.csv`, 90% of empirical Delta importations
   arrive after `tick 234–247` (late July to September 2021, median `tick 286–309`).
   In `exp2`, the primary `E2` wave peaks in mid-June to early July (`ticks 190–217`),
   when **only 0.7% to 2.4% of scheduled `E2` seeds (`36 to 89 seeds` per state)**
   have been introduced:

| State (`E2`) | Total `E2` Seeds | Schedule `p10 / Median / p75` | Seeds Arriving $\le$ **`exp1` Peak** ($R_e > 1$) | Seeds Arriving $\le$ **`exp2` Peak** ($R_e > 1$) |
| :--- | :---: | :---: | :---: | :---: |
| **GA** | 3,778 | `241 / 286 / 330` | **898 (23.8%)** $\le$ tick 254 | **45 (1.2%)** $\le$ tick 207 |
| **MA** | 4,977 | `247 / 309 / 351` | **159 (3.2%)** $\le$ tick 230 | **36 (0.7%)** $\le$ tick 190 |
| **MN** | 4,473 | `246 / 301 / 341` | **971 (21.7%)** $\le$ tick 262 | **81 (1.8%)** $\le$ tick 216 |
| **VA** | 3,322 | `245 / 309 / 353` | **415 (12.5%)** $\le$ tick 247 | **64 (1.9%)** $\le$ tick 217 |
| **WA** | 3,739 | `234 / 297 / 345` | **480 (12.8%)** $\le$ tick 239 | **89 (2.4%)** $\le$ tick 200 |

In `exp2`, the entire primary `E2` wave of 3M–6M infections is founded by those
first 36–89 March–May seeds, while the main July–August importation surge
(`ticks 240–330`) arrives *after* the first `E2` wave has already crashed into
its post-peak herd-immunity trough ($R_e < 1$), fizzling out until Gamma waning
(`mode = 137 d`) triggers a secondary autumn rebound (`ticks 350–450`).

---

## 3. Justification for `exp3`: Holding Real-Time Epidemic Growth Rate ($r$) Invariant Under Implicit NPIs

### Why Unmitigated $R_0 = 4.0$ Cannot Simultaneously Match $\text{GI}$ and Calendar Wave Timing

In any epidemic model, the generation interval ($\text{GI}$), the per-generation
reproduction number ($R$), and the real-time exponential growth rate ($r$, the
Malthusian parameter in $\text{days}^{-1}$) are coupled by:

$$r \approx \frac{\ln R}{\text{GI}}$$

As noted in [exp1](exp1.md), the ABM does not explicitly model time-varying
non-pharmaceutical interventions (NPIs), masking, remote work, or seasonal
contact reductions on the synthetic contact network. In the real US during
Spring–Summer 2021, those mitigations damped effective reproduction numbers so
that Alpha and Delta waves grew over 3–4 months ($r \approx 0.10\text{--}0.14\text{ d}^{-1}$),
peaking in April and August 2021.

* **`exp1` matched the 2021 calendar by compensating errors:** An artificially
  long generation interval ($\text{GI} = 6.87\text{ d}$) divided $\ln R_e$ by a
  1.4× larger denominator, yielding a realistic calendar growth speed
  $r_{\text{v1}} \approx \frac{\ln(2.61)}{6.87\text{ d}} = 0.140\text{ d}^{-1}$
  and $r_{\text{v2}} \approx \frac{\ln(2.11)}{6.87\text{ d}} = 0.109\text{ d}^{-1}$—at
  the cost of a 1.4× too-slow molecular clock.
* **`exp2` fixed $\text{GI}$ (`4.58–4.94 d`) and raised $R$, doubling $r$:**
  Enforcing unmitigated $R_0 \approx 4.0\text{--}5.0$ at $\text{GI} = 4.6\text{--}4.9\text{ d}$
  without an NPI layer forced the simulated waves to grow at unmitigated speed
  ($r \approx 0.24\text{ d}^{-1}$, doubling every ~2.9 days), outrunning the
  calendar seed and vaccination schedules.

### Standard ABM Convention (`exp3`)

When an ABM omits explicit NPI interventions while driving seeding and
vaccination with empirical calendar schedules, global `transmissibility`
absorbs the constant background mitigation factor $(1 - c_{\text{NPI}})$:

$$\beta_{\text{eff}} = (1 - c_{\text{NPI}}) \beta_0$$

Crucially, because `exp2`'s `disease.json` shortened infectious durations
(`I -> rM` and `A -> R`) alongside the generation interval, **retaining `exp1`'s
`calibration.csv` values in `exp3` automatically scales down per-generation $R_e$
in proportion to the shorter $\text{GI}$, keeping the real-time growth rate $r$
almost invariant:**

* **Variant 2 (`E2`, Delta):** In `exp2`/`exp3`, Variant 2's infectivity-weighted
  infectious duration is ~22%–29% shorter than in `exp1` (reducing the
  `v2/v1` duration-weighted ratio from `1.60` to `1.25`), while its Generation
  Interval is **33.3% shorter** (`4.58 d` vs. `6.87 d`). At `exp1`'s
  `transmissibility`, a ~22% reduction in per-generation $R_e$ (from `2.11` at
  10k cases in `exp1` to `~1.65` in `exp3`) paired with $\text{GI} = 4.58\text{ d}$
  gives:
  $$r_{\text{v2, exp3}} \approx \frac{\ln(1.65)}{4.58\text{ d}} = \mathbf{0.1093\text{ d}^{-1}} \quad \text{vs.} \quad r_{\text{v2, exp1}} = \frac{\ln(2.11)}{6.87\text{ d}} = \mathbf{0.1087\text{ d}^{-1}}$$
  matching `exp1`'s calendar wave speed to within **0.5%** while preserving the
  **4.58 d** Generation Interval.
* **Variant 1 (`E1`, Alpha):** Variant 1's infectious duration is **9.7% shorter**
  (`1 / 1.107`) than in `exp1` and its Generation Interval is **4.94 d**,
  eliminating the `1.33×–2.11×` `transmissibility` boost that caused `E1` to
  peak in January–February in `exp2`.
* **Intrinsic Variant Ratio Preserved:** Because `calibration.csv` is a single
  global scalar per state applied equally to both variants, the `exp2` commitment
  to an intrinsic naive ratio of **$R_0(\text{v2}) / R_0(\text{v1}) \approx 1.25$**
  (and the `1.6×` per-contact hazard ratio) is preserved identically in `exp3`.
