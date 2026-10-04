# Bundled data

Small reference files committed to the repository. Everything large is fetched
from Zenodo by `phylogas fetch-data`; see `docs/data_acquisition.md`.

| file | size | source |
|---|---|---|
| `Ruralurbancontinuumcodes2023.csv` | 629 KB | USDA Economic Research Service |
| `county_fips.csv` | 57 KB | derived from the US Census county table |
| `reference/reference.fasta` | 30 KB | SARS-CoV-2 reference, Wuhan/Hu-1/2019 (id matches ncov's refine.root) |

---

## Ruralurbancontinuumcodes2023.csv

USDA Rural-Urban Continuum Codes, 2023 edition. Classifies every US county on
a 1-9 scale from "metro, 1 million+" to "completely rural, not adjacent to a
metro area".

Used by TwinSampler's ascertainment model as the geographic-access modifier:
rural counties get a lower detection probability, reflecting reduced testing
access.

**Source**
<https://www.ers.usda.gov/data-products/rural-urban-continuum-codes>

Direct download (what this copy came from, retrieved 2026-10-01):
```
https://www.ers.usda.gov/media/5768/2023-rural-urban-continuum-codes.csv?v=27068
```

**Rights.** Produced by the USDA Economic Research Service, an agency of the
United States federal government. Works of the US government are not subject
to domestic copyright protection (17 U.S.C. 105) and are in the public domain.
Redistributed here unmodified for reproducibility; please credit USDA ERS.

**Shape.** Long format, one row per county per attribute:

```
FIPS,State,County_Name,Attribute,Value
01001,AL,Autauga County,Population_2020,58805
01001,AL,Autauga County,RUCC_2023,2
```

9,703 rows covering 3,235 counties. `rucc_utils.load_and_pivot_rucc()` pivots
this to wide and renames `RUCC_2023` to `rucc_code`.

**Refreshing.** The copy here is pinned so that results stay reproducible. To
pull the current version instead:

```bash
phylogas fetch-data --with-rucc
```

Note that RUCC is revised roughly every decade (2003, 2013, 2023), and county
FIPS codes do change — Connecticut replaced its counties with planning regions
in 2022, for instance. Re-running an old analysis against a new RUCC file may
not reproduce.

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
