# Salvage audit: `synthetic_biosurveillance`

Performed 2026-09-28, in response to the concern in `vision.txt`:

> Some residuals remain in both the old and probably need to be cleaned up
> but i'm concerned that we not leave behind useful functionality.

Scope: `../synthetic_biosurveillance`. Every `.py` with no same-named
counterpart in PhyloGAS was inspected.

---

## Verdict summary

| Item | Verdict | Action |
|---|---|---|
| `src/genetic_painter.py` | **Stale duplicate** | Delete from sb |
| `src/freq_v02.py`, `src/mutational_models/`, `src/cosine_entro.py`, `src/sampling.py` | **Stale duplicates** | Delete from sb |
| `data/merge_persontrait_with_person.py` | **Valuable — rescued** | Copied to `phylogas/popprep/merge_persontrait.py` |
| `data/get_state_county_lat_longs.py` | **Useful — rescued** | Copied to `phylogas/popprep/county_centroids.py` |
| `hmm/*/msaprocessing_*.py` | **Research dead end** | Leave in sb; do not port |
| `AlternativeMutationalModel/*/freq.py` | **Superseded ancestor** | Leave in sb as provenance |
| `network/network.py` | **Throwaway** | Leave; do not port |
| `junk/pairwise_comparisons.py` | **Third-party** | Leave; not ours |

---

## Detail

### 1. `src/genetic_painter.py` — stale duplicate, safe to delete

984 lines vs 1227 in PhyloGAS. Predates the mutational-model registry
refactor, the reinfection/`infection_id` tracking, the re-importation fix, and
all of the 2026-09 performance work. `src/freq_v02.py` is its symlink twin.

**Nothing to rescue.** PhyloGAS is strictly ahead. Deleting these prevents
someone running the old copy and getting silently wrong lineages (the
re-importation bug alone mis-parented ~392 transmissions per 5.35M-record run).

Recommended:
```bash
cd ../synthetic_biosurveillance
git rm -r src/genetic_painter.py src/freq_v02.py src/mutational_models
git commit -m "genetic painter now lives in PhyloGAS; remove stale copy"
```
Leave a pointer in that repo's README.

### 2. `data/merge_persontrait_with_person.py` — rescued

**This was the most important find.** It builds the demographics CSV that the
painter's `--add_metadata` flag consumes, and nothing in PhyloGAS reproduced it.
Without it, the provenance of `va_2_4_0_demographics.csv` was undocumented.

It joins the EpiHiper persontrait file with the synthetic person file and derives:

- `smh_race` — the Scenario Modeling Hub race categorisation (Latino takes
  precedence over race code)
- `latino` — boolean from `hispanic_code > 1`
- `race` — human-readable labels from the 9-way census race code
- `age_group` — `p/s/a/o/g` expanded to "Preschool (0-4)" etc.
- county name via a FIPS lookup join

Those are exactly the columns in the production `--add_metadata` string:
`gender,county,home_latitude,home_longitude,latino,race,smh_race,age_group`.

Now at `src/phylogas/popprep/merge_persontrait.py`.

### 3. `data/get_state_county_lat_longs.py` — rescued

24 lines. Filters a national county-centroid table down to one state and emits
the FIPS lookup that the script above joins against. Trivial but a required
input, so it travels with its consumer.

Now at `src/phylogas/popprep/county_centroids.py`.

### 4. `hmm/` — research dead end, do not port

Four `msaprocessing_*.py` scripts building HMM libraries per variant
(Global / Virginia / Virginia_delta / Virginia_omicron). That repo's own README
says the global build never ran:

> Although msaprocessing_global.py should work in theory, the global multiple
> sequence alignment is too large to run even on a large memory partition.

Described as an alternative mutational model that "could eventually be used" —
it was not. The entropy-threshold approach won and is what PhyloGAS implements.
Historically interesting, not live functionality. Leave it where it is.

### 5. `AlternativeMutationalModel/` — the painter's ancestor

`va_delta/freq.py` and `va_omicron/freq.py` are earlier generations of the
genetic painter (`freq.py` → `freq_v02.py` → `genetic_painter.py`). The
directory's README documents the feature history: compression, `persontrait_file`,
`add_metadata`, `proportional`.

**All of it is in PhyloGAS already.** Valuable only as provenance for the
lineage of the mutation engine. Leave in sb; do not port.

### 6. `network/network.py` — throwaway

38 lines, hardcoded `output_abridges.csv`, hardcoded node id `476724`, plots a
small contact network with matplotlib. A one-off figure script. PhyloGAS's
`painter/analysis/tree_analysis.ipynb` covers this ground more generally.

### 7. `junk/pairwise_comparisons.py` — third-party

605 lines, header says *"@author: Dr Charles Foster, Virology Research Lab"*.
Someone else's SNP-distance tool, vendored in. Not ours to move; licensing
unclear. It lives in a directory named `junk`.

---

## Duplication between TwinSampler and BeyondBaseline

Not strictly in scope, but the higher risk. Both repos still contain
`scripts/linelist_generation/` **and** `scripts/scenarios_simulation/`, and the
shared files have diverged:

| File | Divergence |
|---|---|
| `scripts/scenarios_simulation/run_all_scenarios.py` | 1526 diff lines |
| `scripts/linelist_generation/simulate_linelist.py` | 314 diff lines |
| `scripts/scenarios_simulation/sampling_algorithms.py` | 58 diff lines |
| `utils/sampling.py` | identical |

Per the new READMEs the intended split is:

- **TwinSampler** owns ascertainment and linelist generation
  (`linelist_generation/`)
- **BeyondBaseline** owns sampling algorithms and evaluation
  (`scenarios_simulation/`)

So each repo should delete the half it does not own. This matters more than the
sb cleanup: two diverging copies of `run_all_scenarios.py` means a result can no
longer be attributed to a specific implementation.

**Recommendation:** before deleting either copy, diff them and confirm the
surviving one is the newer. The 1526-line divergence is large enough that real
functionality may exist on only one side.
