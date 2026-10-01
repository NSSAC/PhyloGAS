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

**Expected:** 33 passing tests, roughly 12 seconds. They need no downloaded
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

---

## 5. Sample and benchmark

```bash
# selection: strategies over the linelist (no ABM truth needed)
beyond-baseline-sweep --linelist ... --population ... --save-samples --outdir runs/

# scoring against ABM truth (needs the simulation)
phylogas benchmark mugration --truth <abm.json> \
    --samples 'runs/*_samples.csv.xz' --infections <allevents.csv.xz>

# both at once, ranked
phylogas compare-strategies --truth <abm.json> --infections <allevents.csv.xz> \
    --linelist ... --population ... --algorithms surs LASSO-Greedy
```

Sanity check on the painted genomes themselves:

```bash
phylogas benchmark sequence --painted <fasta.xz> --infections <allevents.csv.xz>
```

Healthy output is a low mean (~0.2 substitutions per transmission) and
`pairs_over_8 == 0`. A heavy tail means children are inheriting the wrong
parent genome -- please report it.

---

## 6. Known gaps

Please do not spend time on these; they are known.

| Area | Status |
|---|---|
| UCSC cluster TSV | Their TLS certificate expired 2025-07-02. Use `--insecure-download`, or `--input_file` with the copy in `data/importations/`. |
| Nextstrain rule | Written but untested; needs a local `ncov` checkout. Off by default. |
| `benchmark compare` | Exercised only with two augur files. The simulated-parsimony side needs `--save-matrices` output from a real run. |
| RUCC file | Not on Dataverse. Fetch from USDA ERS yourself. |
| Seed coverage | The bundled seed FASTA is Virginia Delta (3,322 seqs). Other states under-cover -- WA has 3,721 E2 importations. Use `--with-seeds`. |
| Scale | Largest local test was 5.3M records. Cluster-scale runs are unverified. |
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
5. **Behaviour at cluster scale**, which we have not been able to test.

Useful context for any report:

```bash
phylogas --version
phylogas status --config config.yaml
pytest tests/ -v
conda list | grep -E "phylogas|beyond-baseline|twin-sampler|pandas|numpy|biopython"
```
