# Staging data from a cluster (Dataverse workaround)

As of 2026-09-30 the UVA Dataverse is unstable. If you already have the digital
twin files on an HPC filesystem, symlink them into the repo instead of
downloading. This page gives the exact layout.

---

## Is Dataverse actually down?

Partly. Measured 2026-09-30:

| endpoint | result |
|---|---|
| website (`/`) | **HTTP 502** |
| metadata API (`/api/datasets/...`) | **HTTP 500** after ~30 s |
| **file access API** (`/api/access/datafile/<id>`) | **303 or 503, roughly 50/50** |

So the download endpoint `fetch-data` uses is *intermittent*, not dead. A full
20.2 MB download succeeded on retry and passed `xz -t`:

```
downloading va_household.csv.xz (id=120566) ...
  HTTP 503; retrying in 5s (attempt 2/4)
done: va_household.csv.xz (20.2 MB)
```

**`phylogas fetch-data` is often still worth trying.** It retries with
exponential backoff, and since 2026-09-30 it retries the *whole transfer*, so a
drop part way through a 500 MB body restarts rather than aborting.

Zenodo (the EpiHiper replicates) is unaffected.

---

## Layout

Place or symlink files here. These are the names `fetch-data` would have
written, and what the default config expects.

```
PhyloGAS/data/
├── county_fips.csv                       # committed to the repo
├── reference/reference.fasta             # committed to the repo
├── Ruralurbancontinuumcodes2023.csv      # USDA ERS, not Dataverse
└── va/                                   # one directory per state code
    ├── va_2_4_0_demographics.csv         # <-- the important one
    ├── va_persontrait_epihiper.txt[.xz]
    ├── va_person.csv[.xz]
    ├── va_household.csv[.xz]
    ├── va_residence_locations.csv[.xz]
    └── va_replicate_0_output.csv.gz      # Zenodo; let fetch-data get this
```

`[.xz]` means **either form works**. The checkers accept a configured
`foo.csv.xz` when the disk has `foo.csv`, and vice versa, because every
consumer reads both through pandas. Cluster copies are usually decompressed.

---

## The shortcut: one symlink

The cluster already holds the **derived** demographics table — the output of
the persontrait → person → household → residence_locations join:

```
/project/biocomplexity/vdh_genomics/synthetic_biosurveillance/
    SARS-Cov2-Biosurveillance-Simulation/data/merged_population_files/va_2_4_0_demographics.csv
```

Symlink that and the four raw files become unnecessary:

```bash
cd PhyloGAS
mkdir -p data/va
ln -s /project/.../merged_population_files/va_2_4_0_demographics.csv data/va/

# EpiHiper replicate still comes from Zenodo, which is up
phylogas fetch-data --states va --with-simulations --no-build-demographics

phylogas status          # should now show only genuinely missing inputs
```

`--no-build-demographics` stops `fetch-data` trying to rebuild a table you
already have.

---

## The long way: raw files

If you need to rebuild the demographics table (a different state, or a changed
join), link all four raw inputs and run the builder:

```bash
cd PhyloGAS
mkdir -p data/va
CLUSTER=/project/bii_nssac/biocomplexity/c4gc_asw3xp/LineList   # adjust

for f in va_persontrait_epihiper.txt va_person.csv va_household.csv va_residence_locations.csv; do
    [ -e "$CLUSTER/$f" ] && ln -sf "$CLUSTER/$f" "data/va/$f" || echo "MISSING: $f"
done

phylogas build-demographics \
    --persontrait data/va/va_persontrait_epihiper.txt \
    --person      data/va/va_person.csv \
    --household   data/va/va_household.csv \
    --residence   data/va/va_residence_locations.csv \
    --fips        data/county_fips.csv \
    --out         data/va/va_2_4_0_demographics.csv
```

**`residence_locations` is not optional** under the v2.4.0 schema: it is the
only source of `home_latitude` / `home_longitude`, which the painter's
`--add_metadata` requires. Without it those columns come out empty and the
builder warns. If it is not on the cluster, use the shortcut above instead.

---

## Multiple states

```bash
CLUSTER=/project/...          # adjust
for st in va ga ma mn wa; do
    mkdir -p "data/$st"
    ln -sf "$CLUSTER/${st}_2_4_0_demographics.csv" "data/$st/" 2>/dev/null \
        || echo "no demographics for $st"
done
phylogas fetch-data --states va ga ma mn wa --with-simulations --no-build-demographics
```

Then point `population.state` in your config at whichever state you are
running. No California replicate is published.

---

## Checking your work

```bash
phylogas validate-config --config config.yaml
```

When a file is found under a different extension than configured, it is
reported and the config name shown beneath:

```
  [okay] demographics (derived)           data/va/va_2_4_0_demographics.csv
  [okay] persontrait                      data/va/va_persontrait_epihiper.txt
         (config says data/va/va_persontrait_epihiper.txt.xz)
```

`phylogas status` then prints the single next command.

---

## Notes

* **Symlinks are fine.** `Path.exists()`, pandas and the placeholder expansion
  all follow them; covered by `tests/test_smoke.py`.
* **Do not copy into `data/example_data/`.** That path is gitignored and
  reserved for the bundled sample slice.
* **`.gitignore` already excludes** `data/<state>/`-style artifacts, so linked
  cluster files will not be accidentally committed. Check with
  `git status --short` if unsure.
* **Prefer symlinks to copies** for the 400-600 MB replicates, unless the
  cluster filesystem is slow from your compute node.
