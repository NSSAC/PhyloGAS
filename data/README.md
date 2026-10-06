# Bundled data

Small reference files committed to the repository. Everything large is fetched
from Zenodo by `phylogas fetch-data`; see `docs/data_acquisition.md`.

| file | size | source |
|---|---|---|
| `county_fips.csv` | 57 KB | derived from the US Census county table |
| `reference/reference.fasta` | 30 KB | SARS-CoV-2 reference, Wuhan/Hu-1/2019 (id matches ncov's refine.root) |

---

## Rural-urban continuum codes (moved)

`Ruralurbancontinuumcodes2023.csv` now ships inside TwinSampler, its only
consumer (`linelist_generation/data/`, with provenance and rights in the
README there). `simulate-linelist` uses that copy by default. To override it,
set `population.rucc_file`; `phylogas fetch-data --with-rucc` downloads the
current USDA edition into `data/` for that purpose.

---

## county_fips.csv

Maps 5-digit county FIPS codes to names. Used by
`phylogas build-demographics` to turn the persontrait file's `county_fips`
into the human-readable `county` column the painter attaches to sequence
metadata.

---

## reference/reference.fasta

The SARS-CoV-2 reference genome (Wuhan/Hu-1/2019, 29,903 bp), optionally prepended
to painted output via `--reference` so downstream phylogenetics has a rooting
sequence.
