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
     **4.92–4.96 d** in `exp2` (mean **4.94 d**) and **5.09–5.30 d** in `exp3`
     (mean **5.22 d**, vs. 5.46 d analytic prediction and 5.5 d in Hart et al. 2022).
   * **Variant 2 (`E2`, Delta):** Realised generation interval measured at
     **4.56–4.59 d** in `exp2` (mean **4.58 d**) and **4.62–4.74 d** in `exp3`
     (mean **4.70 d**, vs. 4.69 d analytic prediction and 4.7 d in Hart et al. 2022).
   * At $\lambda = 0.3430\text{ subs/tx}$, $\text{GI} = 4.70\text{ d}$ yields an
     implied substitution clock of $26.6\text{ subs/genome/yr}$
     ($8.91 \times 10^{-4}\text{ subs/site/yr}$) and, at $s_{\text{path}} \approx 0.885$,
     $\mu_{\text{truth}} \approx 23.6\text{ subs/genome/yr}$
     ($7.89 \times 10^{-4}\text{ subs/site/yr}$, **0.99× real SARS-CoV-2**).
2. **Reshaped Gamma Immunity Waning (`alpha = 4.0, beta = 45.75`, mean 183 d, mode 137 d):**
   * Repeat-infection intervals $\le 30\text{ d}$ fell from **15.1% (`exp1`) to
     0.02%–0.03% (`exp2` and `exp3`)**; intervals $\le 60\text{ d}$ fell from
     **28.0% (`exp1`) to 0.65%–0.75% (`exp2`) and 0.44%–0.63% (`exp3`)**.
   * Maximum `E2` infections in a single host fell from **10 (`exp1`) to 5 (`exp2`)
     and 4 (`exp3`, 5 in MN)**, and windowed `E2` reinfection fractions dropped
     from **25.3% (`exp1`) to 22.1%–27.8% (`exp2`) and 11.4%–16.9% (`exp3`)**.
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

| State | `exp1` & `exp3` `transmissibility` | `exp2` `transmissibility` | `exp2 / exp1` Multiplier |
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
earlier** than the real 2021 calendar (`tickZero = 2020-11-30`).

### Consequences of the Calendar Shift in `exp2`

Because the importation schedules (`data/seed_schedule/<state>-{1,2}.csv`) and
vaccination schedules (`data/vax_schedule/`) are anchored to **fixed real-world
2021 calendar dates**, accelerating the epidemic clock created two structural
distortions:

1. **Loss of `E1`–`E2` Co-Circulation Overlap:**
   In the real US record (and in `exp1`), Alpha (`E1`) peaked in early April
   2021 (`ticks 110–130`) and declined through May–June (`ticks 150–200`),
   directly overlapping the March–June take-off of Delta (`E2`). In `exp2`,
   `E1` burned through 65%–92% of the population by March and reached complete
   extinction (`0 cases/day`) before `E2` began its main rise.
2. **Severe Founder Bottleneck on `E2` Importations:**
   In `data/seed_schedule/<state>-2.csv`, 90% of empirical Delta importations
   arrive after `tick 234–247` (late July to September 2021, median `tick 286–309`).
   In `exp2`, the primary `E2` wave peaked in mid-June to early July (`ticks 190–217`),
   when **only 0.7% to 2.4% of scheduled `E2` seeds (`36 to 89 seeds` per state)**
   had been introduced—stranding the main July–August importation surge inside
   the post-peak herd-immunity trough ($R_e < 1$).

---

## 3. Justification for `exp3`: Holding Real-Time Epidemic Growth Rate ($r$) Invariant Under Implicit NPIs

In any epidemic model, the generation interval ($\text{GI}$), the per-generation
reproduction number ($R$), and the real-time exponential growth rate ($r$, the
Malthusian parameter in $\text{days}^{-1}$) are coupled by:

$$r \approx \frac{\ln R}{\text{GI}}$$

As noted in [exp1](exp1.md), the ABM does not explicitly model time-varying
non-pharmaceutical interventions (NPIs), masking, remote work, or seasonal
contact reductions on the synthetic contact network. When an ABM omits explicit
NPI interventions while driving seeding and vaccination with empirical calendar
schedules, global `transmissibility` absorbs the constant background mitigation
factor $(1 - c_{\text{NPI}})$:

$$\beta_{\text{eff}} = (1 - c_{\text{NPI}}) \beta_0$$

Crucially, because `exp2`'s `disease.json` shortened infectious durations
(`I -> rM` and `A -> R`) alongside the generation interval, **retaining `exp1`'s
`calibration.csv` values in `exp3` automatically scales down per-generation $R_e$
in proportion to the shorter $\text{GI}$, keeping the real-time growth rate $r$
almost invariant:**

* **Variant 2 (`E2`, Delta):** In `exp2`/`exp3`, Variant 2's infectivity-weighted
  infectious duration is ~22%–29% shorter than in `exp1` (reducing the
  `v2/v1` duration-weighted ratio from `1.60` to `1.25`), while its Generation
  Interval is **31.6% shorter** (`4.70 d` vs. `6.87 d`). At `exp1`'s
  `transmissibility`, a ~24% reduction in 10k-cohort $R_e$ (5-state mean `2.11`
  in `exp1` $\to$ `1.60` in `exp3`) paired with $\text{GI} = 4.70\text{ d}$ gives:
  $$r_{\text{v2, exp3}} \approx \frac{\ln(1.60)}{4.70\text{ d}} = \mathbf{0.100\text{ d}^{-1}} \quad \text{vs.} \quad r_{\text{v2, exp1}} = \frac{\ln(2.11)}{6.87\text{ d}} = \mathbf{0.109\text{ d}^{-1}}$$
  preserving `exp1`'s August calendar wave peak while delivering the **4.70 d**
  Generation Interval.
* **Intrinsic Variant Ratio Preserved:** Because `calibration.csv` is a single
  global scalar per state applied equally to both variants, the `exp2` commitment
  to an intrinsic naive ratio of **$R_0(\text{v2}) / R_0(\text{v1}) \approx 1.25$**
  (and the `1.6×` per-contact hazard ratio) is preserved identically in `exp3`.

---

## 4. Measured `exp3` Verification Results (`exp1` vs. `exp2` vs. `exp3`, cell `0.25`)

### A. Generation Intervals, Reinfection Fractions, and Waning Gaps (`exp3`, replicate 0)

| State | `E1` Total | `E1` GI (d) | `E2` Total | `E2` GI (d) | `E2` Painted-Window Edges | `E2` Unique Hosts | `reinf_all` | `reinf_win` (300d) | `max_inf` | `gaps <= 30d` | `gaps <= 60d` |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **GA** | 4,408,224 | **5.30 d** | 5,342,451 | **4.74 d** | 3,917,729 | 4,156,497 | 27.6% | **16.0%** | 4 | **0.03%** | **0.44%** |
| **MA** | 3,254,688 | **5.23 d** | 3,615,712 | **4.71 d** | 2,103,538 | 2,868,302 | 24.9% | **11.4%** | 4 | **0.03%** | **0.63%** |
| **MN** | 3,513,248 | **5.09 d** | 4,108,415 | **4.62 d** | 3,031,622 | 3,194,733 | 27.4% | **16.9%** | 5 | **0.03%** | **0.62%** |
| **VA** | 3,808,552 | **5.20 d** | 4,617,163 | **4.70 d** | 3,213,682 | 3,639,635 | 26.0% | **15.1%** | 4 | **0.03%** | **0.56%** |
| **WA** | 3,100,280 | **5.29 d** | 3,548,756 | **4.74 d** | 2,362,479 | 2,841,641 | 24.1% | **12.5%** | 4 | **0.03%** | **0.55%** |
| **Mean** | 3,616,998 | **5.22 d** | 4,246,499 | **4.70 d** | 2,925,810 | 3,340,162 | 26.0% | **14.4%** | 4.2 | **0.03%** | **0.56%** |

* **Generation Interval:** The 5-state mean `E2` (Delta) GI in `exp3` is **4.70 d**
  (range `4.62–4.74 d`, 2.5% spread), matching the analytic prediction (**4.69 d**)
  to 0.01 d and Hart et al. (2022)'s published estimate (**4.7 d**). The 5-state
  mean `E1` (Alpha) GI is **5.22 d** (vs. `4.94 d` in `exp2`, where higher
  `transmissibility` compressed household transmission into day 1).

### B. Cohort $R_0 / R_e$, Wave Peak Timing, and Pre-Peak Seed Alignment (`scripts/calc_re.py`)

#### Variant 1 (`E1`, Alpha — `tickZero = 2020-11-30`)

| State | Exp | Total `E1` | GI (d) | `R0(1k)` | `R0(10k)` | `Re(50k)` | `Peak7dRe` | Peak Tick (Calendar) | Peak / day |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **GA** | `exp1` | 4,962,374 | 6.87 | 3.34 | 2.21 | 2.43 | 3.58 | tick 125 (`2021-04-04`) | 94,303 |
| **GA** | `exp2` | 7,040,048 | 4.96 | 4.17 | 2.80 | 2.98 | 7.89 | tick  59 (`2021-01-28`) | 226,541 |
| **GA** | **`exp3`** | **4,408,224** | **5.30** | **3.41** | **2.01** | **2.22** | **5.23** | **tick  94 (`2021-03-04`)** | **109,739** |
| **MA** | `exp1` | 3,746,129 | 6.87 | 4.07 | 2.14 | 2.60 | 6.61 | tick 109 (`2021-03-19`) | 103,001 |
| **MA** | `exp2` | 4,819,306 | 4.94 | 3.55 | 3.18 | 3.32 | 6.46 | tick  88 (`2021-02-26`) | 205,902 |
| **MA** | **`exp3`** | **3,254,688** | **5.23** | **2.94** | **2.39** | **2.57** | **3.87** | **tick 104 (`2021-03-14`)** | **116,827** |
| **MN** | `exp1` | 3,944,451 | 6.87 | 4.14 | 2.56 | 2.67 | 6.99 | tick 115 (`2021-03-25`) | 97,433 |
| **MN** | `exp2` | 4,037,878 | 4.93 | 4.18 | 3.20 | 3.11 | 4.92 | tick  94 (`2021-03-04`) | 168,779 |
| **MN** | **`exp3`** | **3,513,248** | **5.09** | **3.06** | **2.67** | **2.72** | **4.26** | **tick  97 (`2021-03-07`)** | **116,934** |
| **VA** | `exp1` | 4,491,567 | 6.87 | 3.52 | 3.52 | 2.78 | 5.98 | tick 122 (`2021-04-01`) | 81,108 |
| **VA** | `exp2` | 5,435,026 | 4.93 | 5.03 | 3.60 | 3.06 | 6.36 | tick 102 (`2021-03-12`) | 201,577 |
| **VA** | **`exp3`** | **3,808,552** | **5.20** | **3.45** | **2.34** | **2.22** | **4.84** | **tick 116 (`2021-03-26`)** | **110,774** |
| **WA** | `exp1` | 3,370,603 | 6.87 | 3.10 | 3.08 | 2.59 | 3.52 | tick 122 (`2021-04-01`) | 60,378 |
| **WA** | `exp2` | 5,046,956 | 4.92 | 4.74 | 2.91 | 3.00 | 6.60 | tick  95 (`2021-03-05`) | 166,489 |
| **WA** | **`exp3`** | **3,100,280** | **5.29** | **3.88** | **1.90** | **2.18** | **5.45** | **tick 105 (`2021-03-15`)** | **78,092** |

#### Variant 2 (`E2`, Delta — Painted Wave)

| State | Exp | Total `E2` | GI (d) | `R0(1k)` | `R0(10k)` | `Re(50k)` | `Peak7dRe` | Peak Tick (Calendar) | Peak / day | Scheduled `E2` Seeds $\le$ Peak |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **GA** | `exp1` | 12,200,297 | 6.87 | 2.27 | 1.85 | 1.88 | 3.04 | tick 254 (`2021-08-11`) | 62,115 | 898 (23.8%) |
| **GA** | `exp2` |  9,030,134 | 4.56 | 3.88 | 2.10 | 2.12 | 4.89 | tick 207 (`2021-06-25`) | 125,494 | 45 (1.2%) |
| **GA** | **`exp3`** | **5,342,451** | **4.74** | **2.18** | **1.71** | **1.49** | **2.31** | **tick 239 (`2021-07-27`)** | **44,795** | **307 (8.1%)** |
| **MA** | `exp1` |  8,751,161 | 6.87 | 2.39 | 2.40 | 2.07 | 2.63 | tick 230 (`2021-07-18`) | 46,723 | 159 (3.2%) |
| **MA** | `exp2` |  5,786,679 | 4.58 | 1.98 | 1.76 | 1.94 | 2.87 | tick 190 (`2021-06-08`) | 58,089 | 36 (0.7%) |
| **MA** | **`exp3`** | **3,615,712** | **4.71** | **1.92** | **1.52** | **1.50** | **2.19** | **tick 243 (`2021-07-31`)** | **29,200** | **392 (7.9%)** |
| **MN** | `exp1` |  8,795,126 | 6.87 | 2.84 | 2.24 | 1.61 | 3.34 | tick 262 (`2021-08-19`) | 63,230 | 971 (21.7%) |
| **MN** | `exp2` |  4,938,020 | 4.57 | 2.37 | 1.66 | 1.43 | 3.73 | tick 216 (`2021-07-04`) | 54,479 | 81 (1.8%) |
| **MN** | **`exp3`** | **4,108,415** | **4.62** | **2.12** | **1.55** | **1.46** | **2.98** | **tick 243 (`2021-07-31`)** | **46,196** | **348 (7.8%)** |
| **VA** | `exp1` | 10,913,015 | 6.87 | 2.23 | 1.97 | 1.83 | 2.31 | tick 247 (`2021-08-04`) | 51,860 | 415 (12.5%) |
| **VA** | `exp2` |  6,850,610 | 4.59 | 2.42 | 2.15 | 2.11 | 2.70 | tick 217 (`2021-07-05`) | 68,709 | 64 (1.9%) |
| **VA** | **`exp3`** | **4,617,163** | **4.70** | **1.70** | **1.65** | **1.73** | **1.82** | **tick 266 (`2021-08-23`)** | **36,253** | **749 (22.5%)** |
| **WA** | `exp1` |  8,820,893 | 6.87 | 2.96 | 2.08 | 1.70 | 3.32 | tick 239 (`2021-07-27`) | 38,953 | 480 (12.8%) |
| **WA** | `exp2` |  6,235,721 | 4.58 | 2.75 | 1.75 | 1.81 | 3.50 | tick 200 (`2021-06-18`) | 52,841 | 89 (2.4%) |
| **WA** | **`exp3`** | **3,548,756** | **4.74** | **2.12** | **1.59** | **1.36** | **2.25** | **tick 248 (`2021-08-05`)** | **24,119** | **702 (18.8%)** |

* **5-State Mean `E2` Peak Tick:** `exp3` peaks at **tick 247.8 (`2021-08-05`)** vs.
  **tick 246.4 (`2021-08-03`)** in `exp1` (a 1.4-day difference) and **tick 206.0
  (`2021-06-24`)** in `exp2` (+41.8 days later than `exp2`), restoring `4.3×–11.7×`
  more pre-peak `E2` importations (`307–749 seeds` per state).
