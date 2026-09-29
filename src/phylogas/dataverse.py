"""UVA Dataverse acquisition for the PGCOE synthetic populations.

Each PGCOE state is a **separate** Dataverse deposit; there is no single DOI
covering all six. The file-ID table below was resolved from the datasets API
and verified against the served ``Content-Disposition`` filenames on
2026-09-28.

The IDs are hardcoded deliberately. Discovering them by querying each DOI adds
six network round-trips and an extra failure mode, and the IDs are stable
identifiers by design.
"""

from __future__ import annotations

import lzma
import shutil
import sys
from pathlib import Path

API = "https://dataverse.lib.virginia.edu/api/access/datafile"
DATASET_API = "https://dataverse.lib.virginia.edu/api/datasets/:persistentId/?persistentId="

# state -> Dataverse DOI
STATE_DOIS = {
    "va": "doi:10.18130/V3/5LSDCY",
    "ca": "doi:10.18130/V3/A7DQWM",
    "ga": "doi:10.18130/V3/4DWNRH",
    "ma": "doi:10.18130/V3/ZB0SGL",
    "mn": "doi:10.18130/V3/SB2PWT",
    "wa": "doi:10.18130/V3/PNGMRJ",
}

# The four files needed to build the demographics table used by the painter
# and by TwinSampler. Sizes are the compressed (.xz) sizes in MB.
CORE_FILES = {
    "va": {
        "va_persontrait_epihiper.txt.xz": (121039, 25.5),
        "va_person.csv.xz": (120565, 25.4),
        "va_household.csv.xz": (120566, 20.2),
        "va_residence_locations.csv.xz": (120568, 32.9),
    },
    "ca": {
        "ca_persontrait_epihiper.txt.xz": (121117, 106.5),
        "ca_person.csv.xz": (120625, 107.3),
        "ca_household.csv.xz": (120623, 124.9),
        "ca_residence_locations.csv.xz": (120621, 116.7),
    },
    "ga": {
        "ga_persontrait_epihiper.txt.xz": (121106, 29.9),
        "ga_person.csv.xz": (120598, 29.8),
        "ga_household.csv.xz": (120595, 23.8),
        "ga_residence_locations.csv.xz": (120596, 39.8),
    },
    "ma": {
        "ma_persontrait_epihiper.txt.xz": (121110, 20.0),
        "ma_person.csv.xz": (120589, 20.2),
        "ma_household.csv.xz": (120591, 16.5),
        "ma_residence_locations.csv.xz": (120590, 22.4),
    },
    "mn": {
        "mn_persontrait_epihiper.txt.xz": (121112, 16.8),
        "mn_person.csv.xz": (120577, 17.0),
        "mn_household.csv.xz": (120581, 13.9),
        "mn_residence_locations.csv.xz": (120580, 27.2),
    },
    "wa": {
        "wa_persontrait_epihiper.txt.xz": (121104, 21.4),
        "wa_person.csv.xz": (120604, 21.5),
        "wa_household.csv.xz": (120607, 17.9),
        "wa_residence_locations.csv.xz": (120606, 30.2),
    },
}

# EpiHiper *inputs*. PhyloGAS consumes EpiHiper OUTPUT (output.csv.gz), so
# these are not needed unless you intend to run the ABM yourself. The contact
# network alone is 7.1 GB across the six states - 88% of the total download -
# which is why it is opt-in.
EPIHIPER_INPUT_FILES = {
    "va": {"va_contact_network_epihiper.txt.xz": (121037, 907.2)},
    "ca": {"ca_contact_network_epihiper.txt.xz": (121116, 4025.4)},
    "ga": {"ga_contact_network_epihiper.txt.xz": (121107, 1094.8)},
    "ma": {"ma_contact_network_epihiper.txt.xz": (121109, 719.9)},
    "mn": {"mn_contact_network_epihiper.txt.xz": (121114, 589.4)},
    "wa": {"wa_contact_network_epihiper.txt.xz": (121105, 788.5)},
}


def file_plan(states, with_epihiper_inputs: bool = False):
    """Return [(state, filename, file_id, size_mb), ...] for the requested states."""
    plan = []
    for st in states:
        st = st.lower()
        if st not in CORE_FILES:
            raise KeyError(st)
        for name, (fid, mb) in CORE_FILES[st].items():
            plan.append((st, name, fid, mb))
        if with_epihiper_inputs:
            for name, (fid, mb) in EPIHIPER_INPUT_FILES.get(st, {}).items():
                plan.append((st, name, fid, mb))
    return plan


def download_file(file_id: int, dest: Path, timeout: int = 1800) -> bool:
    """Stream one Dataverse file to ``dest``. Returns True if newly downloaded."""
    if dest.exists() and dest.stat().st_size > 0:
        print(f"    exists, skipping: {dest.name}")
        return False

    try:
        import requests
    except ImportError:
        raise SystemExit("ERROR: `requests` is required. pip install requests")

    tmp = dest.with_suffix(dest.suffix + ".part")
    url = f"{API}/{file_id}"
    print(f"    downloading {dest.name} (id={file_id}) ...", flush=True)
    with requests.get(url, stream=True, timeout=timeout, allow_redirects=True) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        done = 0
        step = max(total // 20, 8 << 20) if total else (8 << 20)
        nxt = step
        with tmp.open("wb") as fh:
            for chunk in r.iter_content(chunk_size=1 << 20):
                fh.write(chunk)
                done += len(chunk)
                if done >= nxt:
                    pct = f"{100 * done / total:5.1f}%" if total else "     "
                    print(f"      {pct}  {done / 1048576:8.1f} MB", flush=True)
                    nxt += step
    tmp.rename(dest)
    print(f"    done: {dest.name} ({dest.stat().st_size / 1048576:.1f} MB)")
    return True


def decompress_xz(path: Path, keep: bool = False) -> Path:
    """Decompress a .xz file. Returns the decompressed path.

    Note: pandas reads .xz natively, so decompression is only needed for tools
    that cannot. It is off by default in the CLI for that reason - the v2.4.0
    population files expand roughly 6x.
    """
    if path.suffix != ".xz":
        return path
    out = path.with_suffix("")
    if out.exists():
        print(f"    already decompressed: {out.name}")
        return out
    print(f"    decompressing {path.name} ...", flush=True)
    with lzma.open(path, "rb") as fin, out.open("wb") as fout:
        shutil.copyfileobj(fin, fout, length=16 << 20)
    if not keep:
        path.unlink()
    return out
