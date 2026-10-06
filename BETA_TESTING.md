# Beta testing PhyloGAS

Thanks for trying this. This page is the shortest path from a clean machine to
a painted genome set, plus an honest account of what does not work yet.

Please report problems with the output of `phylogas status` and
`pytest tests/ -v` attached.

---

## 1. Install

```bash
# clone and create the environment
git clone https://github.com/NSSAC/PhyloGAS.git
cd PhyloGAS
mamba env create -f environment.yml
conda activate phylogas_env
pip install -e .

# confirm the install
pytest tests/ -v
phylogas --version
```

`environment.yml` also installs the two satellite repositories
(`TwinSampler`, `BeyondBaseline`) straight from GitHub. If pip fails on either,
that is the bug -- tell us, do not work around it.

**Expected:** 67 passing tests, roughly 10 seconds. They need no downloaded
data. If any fail, stop and send us the output.

---

## 2. Find out what to do next

```bash
cp config.template.yaml config.yaml
phylogas status
```

`status` walks the pipeline stages and prints the single next command. Run it
whenever you are unsure where you are:

```
PhyloGAS pipeline status
====================================================
  [okay] Configuration            config.yaml
  [MISS] Synthetic population     data/va/va_2_4_0_demographics.csv
  [MISS] EpiHiper output          data/va/va_replicate_0_output.csv.gz
  ...
Next step:
  phylogas fetch-data --states va --with-simulations
```

---

## 3. Get the data

```bash
phylogas fetch-data --states va --with-simulations --dry-run   # sizes first
phylogas fetch-data --states va --with-simulations
```

Roughly 0.62 GB for Virginia: four synthetic-population files from UVA
Dataverse plus one EpiHiper replicate from Zenodo
(`doi:10.5281/zenodo.23067670`). Everything is MD5-verified; re-running
re-verifies rather than blindly trusting what is on disk.

Available states: `va`, `ga`, `ma`, `mn`, `wa`. No California replicate has
been published.

The command also builds `<state>_2_4_0_demographics.csv`, which is a derived
join -- **not** the raw EpiHiper persontrait file. See section 7.

---

## 4. Run the pipeline

```bash
phylogas run --config config.yaml --cores all
phylogas run --config config.yaml --cores all --until paint_only   # stop after painting
phylogas run --config config.yaml --profile slurm                  # on a cluster
```

Or drive one stage at a time:

```bash
phylogas train --config config.yaml   # entropy map from an MSA
phylogas paint --config config.yaml   # paint genomes onto the network
```

Add `--dry-run` to any command to see what would execute.

**Rough expectations** (measured on a laptop, not a cluster):

| | |
|---|---|
| paint throughput | ~6,000-8,000 records/s |
| 390k records | ~1 minute |
| 5.3M records | ~20 minutes |
| peak memory | ~17 GB for a 5.3M-record run |

If painting is dramatically slower, check `compression_level` in your config.
It defaults to `1`; the lzma default of `6` is roughly 37x slower here.

The full chain — painting, ascertainment, sampling, one Nextstrain build and
the clock benchmark — has been run on a SLURM cluster for the Virginia Delta
replicate (5.35M infections). The clock benchmark over all 5.35M painted
sequences takes a few minutes.

### Turning on Nextstrain

Off by default. Set `nextstrain.enabled: true` and `nextstrain.dir` to an
`ncov` checkout (README, "Optional: Nextstrain"), then pick the trees:

```yaml
sampling:
  nextstrain_recipes: ["4S__surs"]    # or several, or "all"
```

Each entry yields `results/<project>/04_nextstrain_builds/<recipe>/auspice.json`.
`scenarios-recipes` lists the valid ids.

---

## 5. Sample and benchmark

`phylogas run` does the sampling for you (rule `sample_scenarios`). Each
sample file is named after its recipe — `<scenario>__<algorithm>_samples.csv.xz`,
e.g. `4S-4_LL-P__lasso_greedy_samples.csv.xz` or `4S__surs_samples.csv.xz`.

By hand:

```bash
# selection: strategies over the linelist (no ABM truth needed)
beyond-baseline-sweep --linelist ... --population ... --save-samples --outdir runs/

# scoring against ABM truth (needs the simulation)
phylogas benchmark mugration --truth <abm.json> \
    --samples 'runs/*_samples.csv.xz' --infections <allevents.csv.xz>

# both at once, ranked
phylogas compare-strategies --truth <abm.json> --infections <allevents.csv.xz> \
    --linelist ... --population ... --algorithms surs lasso_greedy
```

Sanity check on the painted genomes themselves:

```bash
phylogas benchmark sequence --painted <fasta.xz> --infections <allevents.csv.xz>
```

Healthy output is a low mean (~0.2 substitutions per transmission) and
`pairs_over_8 == 0`. A heavy tail means children are inheriting the wrong
parent genome -- please report it.

### The clock benchmark

With Nextstrain enabled, `rule all` also writes
`05_benchmarks/clock_estimates.csv`. Standalone:

```bash
phylogas benchmark clock --config config.yaml --mode truth
```

How to read it (full detail in [`docs/clock_modes.md`](docs/clock_modes.md)):

| row | what it is |
|---|---|
| `mu_truth` | substitutions/site/year the painter actually produced, pooled root-to-tip over every painted sequence, no tree |
| `mu_lineage` | the same, fitted within transmission chains only — the per-lineage rate |
| `mu_between_chains` | the part of `mu_truth` that comes from differences between imported seeds; `mu_truth = w·mu_lineage + (1−w)·mu_between_chains` |
| `mu_seed_trend*` | divergence of the imported seeds against their import date |
| `mu_operational` | the rate the Nextstrain build assumed (0.0008) |
| `mu_sim` | TreeTime's unconstrained estimate on the build — needs augur, see Known gaps |

These numbers characterise **the simulation**, and how far its rate is from
the one Nextstrain was told to use. They are not, on their own, a validation
of Nextstrain (that needs `mu_sim`) or of the simulation against reality
(that needs the same measurement on real sequences, `--kind real`).

---

## 6. Known gaps

Please do not spend time on these; they are known.

| Area | Status |
|---|---|
| UCSC cluster TSV | Their TLS certificate expired 2025-07-02. Use `--insecure-download`, or `--input_file` with the copy in `data/importations/`. |
| Nextstrain clock | Builds use a fixed 0.0008 subs/site/yr. The painter produces ~2.8–2.9e-4, so dated trees are stretched to the assumed rate. Deliberate for now (it is what a health department would run); the gap is what `clock_estimates.csv` reports. |
| `mu_sim` | `benchmark clock --mode inferred` needs `augur` on `PATH`. Augur lives in the Nextstrain environment, not PhyloGAS's, so in the pipeline this row is skipped with a note. Run it inside `nextstrain shell`. |
| Skyline coalescent | Fails on these trees (no gradient). `cfg/nextstrain/ncov/base.yaml` sets `coalescent: opt`. |
| Sampling replicates | `sampling.replicates` is read but not wired: one sample set per recipe. |
| Reporting delay | Not modelled. Line-list dates are state onset, not collection or report date. |
| Seed dating | Imported-seed divergence trends *down* with import date in the VA run, which suggests the importation schedule does not pair each seed with its real collection date. Under investigation. |
| Metadata ownership | Chain/component ids come from TwinSampler while sequences come from the painter; who should own the joined metadata is open — [`docs/design/metadata_ownership.md`](docs/design/metadata_ownership.md). |
| `benchmark compare` | Exercised only with two augur files. The simulated-parsimony side needs `--save-matrices` output from a real run. |
| UVA Dataverse | Unstable as of 2026-09-30: website 502, metadata API 500, file API ~50% 503. Downloads usually succeed on retry. If you have the files on a cluster, see [`docs/staging_local_data.md`](docs/staging_local_data.md). |
| Seed coverage | The bundled seed FASTA is Virginia Delta (3,322 seqs). Other states under-cover -- WA has 3,721 E2 importations. Use `--with-seeds`. |
| Scale | One state, one replicate (VA, 5.35M infections) verified end to end on a cluster. Other states and multi-replicate runs are not. |
| First Nextstrain build | Needs internet: ncov downloads a nextclade dataset into the checkout. Run it from a login node if compute nodes are offline. |
| Docker/Apptainer | Not built. |

---

## 7. Things that will confuse you

**`demographics_file` is not the persontrait file.** The raw EpiHiper
persontrait opens with a JSON schema line, carries `county_fips` rather than
`county`, `hispanic` rather than `latino`, and in v2.4.0 has no
`home_latitude`/`home_longitude` at all -- those live in
`residence_locations`, reached via `household.rlid`. `phylogas fetch-data`
builds the joined table for you; `phylogas build-demographics` does it
standalone. Hand the painter a raw persontrait file and it will now refuse,
with the command you need.

**Replicates carry two variants.** Each EpiHiper run has both `E1` and `E2`
waves. `genetic_painter.painted_prefix` selects which one. The default `E2` is
the second (Delta-like) wave.

**Sublineages are included by default.** `--pango B.1.617.2` means
`B.1.617.2*`, i.e. the whole Delta clade -- 21,034 Virginia samples versus 836
for the exact lineage. Use `--no-include-sublineages` to narrow it.

**Recipes, not algorithms, name a tree.** "SURS" alone does not identify a
sample: `1S__surs` and `4S__surs` are different sets. A recipe id is
`<scenario>__<algorithm slug>`. SURS ignores the target distribution, so its
ids carry the stride only; `4S-4_LL-P__surs` is accepted as an alias for
`4S__surs`. `phylogas nextstrain-config` takes `--recipe`, not `--algo`.

**Deleting an intermediate does not make Snakemake rebuild it** when the
final targets are already newer. Use `--forcerun <rule>`, which also reruns
everything downstream:
`snakemake --configfile config.yaml --cores all --forcerun simulate_linelist -n`.

**`allevents` is one row per infection.** `linelist_allevents.csv.xz` used to
carry one row per clinical state (~1.7 per infection), and reinfections after a
person's first detection were silently never ascertained. Both are fixed in
TwinSampler; line lists, prevalence and any count from older outputs will
change when regenerated. Regenerate rather than mixing old and new.

**`alias_pid` is the infection, `sim_pid` is the person.** `alias_pid`
(`"{pid}.{exposure_tick}"`) is 1:1 with the sequence `strain`. A person
reinfected has two. Reinfection in this disease model follows exponential
waning of immunity (mean ~183 days) — see `cfg/exp1/disease.json`.

**`date` is onset, not collection.** TwinSampler records the onset of the
clinical state that was ascertained. `exposure_date` is when the infection
was acquired. Nextstrain treats `date` as the sampling date.

**Older SURS samples will not reproduce bit-for-bit.** SURS is now seeded
from `--seed`; before, it ignored the seed entirely.

**`compression_level` defaults to 1, not 6.** Much faster, files a few times
larger. Set `6` with `compression_threads: 8` if archive size matters more
than wall time.

---

## 8. What we most want to hear about

1. **Install failures** -- especially the satellite repos.
2. **Anything `phylogas status` gets wrong** about where you are.
3. **Error messages that do not tell you what to do next.** Most now suggest a
   fix; the ones that do not are bugs.
4. **Biologically implausible output** -- `benchmark sequence` reporting a
   heavy divergence tail, or mutation rates that look wrong against real
   surveillance.
5. **Behaviour at cluster scale** beyond the one verified state/replicate —
   other states, several replicates, `nextstrain_recipes: "all"`.
6. **Clock numbers that look off** — paste `clock_estimates.csv`.

Useful context for any report:

```bash
phylogas --version
phylogas status --config config.yaml
pytest tests/ -v
conda list | grep -E "phylogas|beyond-baseline|twin-sampler|pandas|numpy|biopython"
```
