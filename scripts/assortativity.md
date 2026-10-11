# assortativity.py — what it measures, and how to read it

Standalone, runnable from the PhyloGAS repo root, writes a small CSV.
Does not modify pipeline outputs.

## Does a sampling strategy escape the epidemic's homophily?

### Context

The poster claims that transmission networks are strongly homophilic, that
Simple Uniform Random Sampling (SURS) gets trapped inside densely-ascertained
demographically-uniform clusters, and that LASSO-Greedy uses demographic
targeting as a topological compass to cross those boundaries.

Graph coverage cannot test that claim: both strategies can touch many nodes.
What separates them is the *character of the transmission edges they retain*.

### What it does

1. Builds the true transmission graph from the ABM all-events file — one
   undirected edge per infection event, joining `alias_contact` (the
   infector's infection id) to `alias_pid` (the infectee's). Root importations
   have no infector and contribute no edge.
2. Joins demographics onto every infection via `sim_pid`.
3. Computes Newman's categorical assortativity `r` on the **full** graph. This
   is the epidemic's natural homophily and is the only meaningful reference.
4. For each sampling strategy, induces the subgraph on sampled infections —
   keeping an edge only when **both** endpoints were sampled — and recomputes.

### Running it

```
python scripts/assortativity.py \
    --infections results/va_delta_wave/05_benchmarks/allevents_variants.csv.xz \
    --samples 'results/va_delta_wave/03_sampled_datasets/*_samples.csv.xz' \
    --state va \
    --out poster_project/data/assortativity_va.csv
```

Or all five states at once, from a compute node:

```
parallel -j 5 --tag --line-buffer --joblog slurm/logs/assort.joblog \
  'python scripts/assortativity.py \
     --infections results/{}_delta_wave/05_benchmarks/allevents_variants.csv.xz \
     --samples "results/{}_delta_wave/03_sampled_datasets/*_samples.csv.xz" \
     --state {} --out poster_project/data/assortativity_{}.csv' \
  ::: va ga ma mn wa
```

`--population` is optional and usually unnecessary: the all-events file
already carries `age_group, smh_race, race, county, county_fips, rucc_code,
ses_category, sex`. Pass it only for an attribute that file lacks.

Attributes are auto-detected; override with `--attributes age_group,smh_race,county`.
Redundant pairs are dropped automatically (`age` when `age_group` is present,
`race` when `smh_race` is, fips variants when `county` is).

### How to interpret — read this before plotting anything

Every row is one sampling strategy. The row with `strategy=__full_graph__` is
the unsampled epidemic: **the baseline every other row must be compared
against.** An absolute assortativity value means nothing on its own.

| column | meaning |
|---|---|
| `n_nodes_sampled` | infections in this sample |
| `n_edges_retained` | true transmission edges with **both** endpoints sampled |
| `edge_retention` | that count over all true edges |
| `r_<attr>` | Newman's assortativity over retained edges, −1…+1 |
| `within_<attr>_frac` | share of retained edges joining two same-category infections |
| `cross_county_edge_frac` | share of retained edges whose endpoints differ in county |
| `n_counties_touched` | distinct counties appearing in retained edges |

**Direction of merit.** There is no universally "good" direction here; the
comparison against baseline is the finding:

- `r` **above** baseline → the sample is *more* homophilic than the epidemic.
  It settled inside demographic islands. This is the predicted SURS signature.
- `r` **below** baseline → the sample spans demographic boundaries more than
  the disease itself does. It broke homophilic gravity. This is the predicted
  LASSO-Greedy signature, and the poster's mechanistic claim.

**`edge_retention` is not optional context.** A strategy that scatters samples
so thinly that few true edges retain both endpoints will report a low `r` for
a reason that is a *failure*: it shattered the transmission tree instead of
traversing it. This is the expected behaviour of LASSO-Greedy with very many
micro-strata (~1500 buckets), and it is why that arm belongs in the panel —
without it, low assortativity reads as unqualified success. Never plot `r`
without `edge_retention` beside it. Below `--min-edges` (default 200) `r` is
written empty rather than reported as noise.

**Why both `r_<attr>` and `within_<attr>_frac`.** `r` is normalised against the
category marginals, which makes it comparable across attributes but also
saturates near 1 on a strongly clustered graph — differences then land in the
third decimal and plot poorly. `within_<attr>_frac` is the raw mixing-matrix
diagonal: less principled, more interpretable, and it keeps moving at the
ceiling. Use `r` for the claim and the fraction for the axis if `r` compresses.

**The spatial columns often discriminate hardest.** In validation, strategies
separated ~10× more on `cross_county_edge_frac` and `n_counties_touched` than
on `r`. These are also the columns comparable to the earlier cross-replicate
sweep's "boundary node" (34.4% vs 29.0%) and "counties touched" (109 vs 127)
figures — but note that sweep was pre-seed-QC, so do not place its numbers and
these in the same panel unlabelled.

### Validation

`newman_r` was checked against analytic cases: all-within-category edges give
r = +1; a fully crossing two-category graph gives r = −1; a single category
returns `None` (undefined, not 0); uniform random mixing over 200k edges gives
r = −0.0007; edges with a missing attribute are dropped rather than pooled
into a category.

End-to-end, the script was run on a synthetic islands-and-bridges graph: eight
demographic islands of 500 infections, sparse cross-island bridges, and three
sampling styles at equal budget. It reproduced the expected signature — a
concentrated sample landing above baseline while touching 2 of 8 counties with
0.1% cross-county edges, and spread samples landing below baseline while
touching all 8 with ~10× the cross-county fraction. Rebuild that fixture from
the commands in this repo's session notes if you want to re-verify.

### Known limitations

- Assortativity is computed on the *true* transmission edges restricted to
  sampled nodes. It describes what the sample could in principle reveal, not
  what the phylogenetic reconstruction actually recovered. The reconstruction
  side is `Mugration_Metrics.csv` (`cosine_similarity`, `topological_f1`).
- Categorical assortativity ignores ordering, so for age it treats adjacent
  bands as no more similar than distant ones. For an ordinal reading, a
  numeric assortativity on midpoint age would be a separate measure.
- One replicate per state unless you pass samples from several.

---

## Note: `stride_variant_prevalence_error` is degenerate in this regime

In the five-state run of 2026-10-10 this metric returned a single identical
value across all eight sampling strategies in four of five states (2.000; MA
had two distinct values). It reads the correct column -- `variant_benchmark`,
which is populated -- so this is not a plumbing fault. The metric was designed
when the benchmark carried ten distinct variant labels; with every infection
now a Delta sublineage there is effectively one variant circulating, the
estimation task is trivial, and the error saturates at its ceiling regardless
of which cases were sampled. Do not report it without first re-deriving it
over sublineage prevalence.
