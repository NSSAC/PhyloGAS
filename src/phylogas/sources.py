"""Remote data acquisition for PhyloGAS.

Two upstream sources, each holding a different half of the inputs:

* **UVA Dataverse** - the synthetic populations (persontrait, person,
  household, residence_locations). One deposit per PGCOE state.
* **Zenodo** - the EpiHiper agent-based simulation outputs that PhyloGAS
  paints onto. A single record holds one replicate per state.

Original module docstring follows.

UVA Dataverse acquisition for the PGCOE synthetic populations.

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
    with _get_with_retry(url, timeout) as r:
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


# ==========================================================================
# Zenodo: EpiHiper simulation outputs
# ==========================================================================
# "Phylogeographic Analysis Similars - Agent Based Simulations"
# Verified open-access 2026-09-30; MD5s confirmed against the upload log.
ZENODO_RECORD = "23067670"
ZENODO_DOI = "10.5281/zenodo.23067670"
ZENODO_CONCEPT_DOI = "10.5281/zenodo.23067669"   # always resolves to the latest version
ZENODO_FILES_URL = "https://zenodo.org/records/{record}/files/{key}?download=1"

# state -> (filename, size_mb, md5)
SIMULATION_FILES = {
    "va": ("va_replicate_0_output.csv.gz", 530.3, "2971f450fde914fd12b420096b90848b"),
    "ga": ("ga_replicate_0_output.csv.gz", 583.7, "508eab5b6a04ec44982ceba62306ae28"),
    "ma": ("ma_replicate_0_output.csv.gz", 431.8, "bfdcc576f006b10ba114165033ac5088"),
    "mn": ("mn_replicate_0_output.csv.gz", 428.1, "ac6f078664fedade69058d9bf0b8fb9e"),
    "wa": ("wa_replicate_0_output.csv.gz", 425.2, "9be12c1f872654e5af64cb2352fcb263"),
    # NOTE: no California replicate in this record.
}


def simulation_url(state: str) -> str:
    key, _, _ = SIMULATION_FILES[state.lower()]
    return ZENODO_FILES_URL.format(record=ZENODO_RECORD, key=key)


def simulation_plan(states):
    """Return [(state, filename, url, size_mb, md5), ...] for available states."""
    plan = []
    for st in states:
        st = st.lower()
        if st not in SIMULATION_FILES:
            continue
        key, mb, md5 = SIMULATION_FILES[st]
        plan.append((st, key, simulation_url(st), mb, md5))
    return plan


def _get_with_retry(url: str, timeout: int, attempts: int = 4):
    """GET with exponential backoff.

    Dataverse in particular returns transient 503s under load; a single
    failure part way through a multi-state fetch is not worth aborting on.
    """
    import time

    import requests

    last = None
    for i in range(attempts):
        try:
            r = requests.get(url, stream=True, timeout=timeout, allow_redirects=True)
            if r.status_code in (500, 502, 503, 504):
                last = f"HTTP {r.status_code}"
                r.close()
            else:
                r.raise_for_status()
                return r
        except Exception as exc:  # network error, timeout, etc.
            last = str(exc)
        if i < attempts - 1:
            wait = 5 * (2 ** i)
            print(f"      {last}; retrying in {wait}s "
                  f"(attempt {i + 2}/{attempts})", flush=True)
            time.sleep(wait)
    raise RuntimeError(f"failed after {attempts} attempts: {last}")


def md5sum(path: Path, chunk: int = 1 << 20) -> str:
    import hashlib

    h = hashlib.md5()
    with Path(path).open("rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def download_url(url: str, dest: Path, expect_md5: str | None = None,
                 timeout: int = 3600) -> bool:
    """Stream a URL to ``dest``, optionally verifying its MD5.

    Returns True if newly downloaded. A file that already exists is verified
    (when a checksum is known) rather than blindly trusted, so a truncated
    earlier attempt is caught instead of silently reused.
    """
    if dest.exists() and dest.stat().st_size > 0:
        if expect_md5:
            print(f"    verifying existing {dest.name} ...", flush=True)
            got = md5sum(dest)
            if got == expect_md5:
                print(f"    exists and checksum matches: {dest.name}")
                return False
            print(f"    checksum MISMATCH on existing file (got {got[:12]}..., "
                  f"expected {expect_md5[:12]}...); re-downloading", flush=True)
            dest.unlink()
        else:
            print(f"    exists, skipping: {dest.name}")
            return False

    try:
        import requests
    except ImportError:
        raise SystemExit("ERROR: `requests` is required. pip install requests")

    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"    downloading {dest.name} ...", flush=True)
    with _get_with_retry(url, timeout) as r:
        total = int(r.headers.get("content-length", 0))
        done = 0
        step = max(total // 20, 32 << 20) if total else (32 << 20)
        nxt = step
        with tmp.open("wb") as fh:
            for block in r.iter_content(chunk_size=1 << 20):
                fh.write(block)
                done += len(block)
                if done >= nxt:
                    pct = f"{100 * done / total:5.1f}%" if total else "     "
                    print(f"      {pct}  {done / 1048576:8.1f} MB", flush=True)
                    nxt += step
    tmp.rename(dest)

    if expect_md5:
        print(f"    verifying checksum ...", flush=True)
        got = md5sum(dest)
        if got != expect_md5:
            dest.unlink()
            raise SystemExit(
                f"ERROR: checksum mismatch for {dest.name}\n"
                f"       expected {expect_md5}\n"
                f"       got      {got}\n"
                "       The download was removed; re-run to try again."
            )
        print(f"    checksum OK")
    print(f"    done: {dest.name} ({dest.stat().st_size / 1048576:.1f} MB)")
    return True
