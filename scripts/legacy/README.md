# Legacy scripts

Kept for provenance and for the record of where the Dataverse file IDs came
from. **Superseded — prefer the `phylogas` CLI.**

## `download_pgcoe_dataverse-v2.sh`

Replaced by `phylogas fetch-data`, which was built from this script's file-ID
table. The table was verified against the live Dataverse API on 2026-09-28 and
is now maintained in `src/phylogas/dataverse.py`.

Do not use this script for new work. Three concrete problems:

1. **It cannot produce a usable demographics file.** It is missing
   `<state>_residence_locations.csv.xz`, and under the v2.4.0 schema that file
   is the only source of `home_latitude` / `home_longitude`. Those columns are
   required by the painter's `--add_metadata`.
2. **It downloads ~7.9 GB of contact networks nobody reads.** Confirmed by
   searching all three repos: PhyloGAS consumes EpiHiper *output*, not the
   contact-network *input*.
3. **`unxz -k` keeps both forms**, roughly doubling an already large footprint
   (~40 GB for all six states). pandas reads `.xz` natively, so decompression
   is usually unnecessary.

Equivalent modern invocation:

```bash
phylogas fetch-data --states va ga            # + builds the demographics table
phylogas fetch-data --states va --dry-run     # preview sizes first
```
