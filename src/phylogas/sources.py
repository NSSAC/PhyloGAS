"""Remote data acquisition for PhyloGAS.

Everything is fetched from **Zenodo**, in two groups of records:

* **Synthetic populations** - one record per PGCOE state, mirroring the UVA
  Dataverse deposits. Supplies persontrait, person, household and
  residence_locations, which are joined into the demographics table.
* **EpiHiper simulation outputs** - a single record holding one transmission
  replicate per state. This is what ``phylogas paint`` consumes.

Zenodo is the primary source because UVA Dataverse proved unreliable: measured
2026-09-30 its website returned 502, its metadata API 500, and its file API
roughly 50% 503. The Dataverse identifiers are retained below as provenance
and as a manual fallback, but nothing fetches from them by default.

Record IDs, filenames, sizes and MD5 checksums in the tables below were
harvested directly from the Zenodo API on 2026-10-01, not transcribed by hand.
"""

from __future__ import annotations

import lzma
import shutil
import sys
from pathlib import Path

# --------------------------------------------------------------------------
# Zenodo: synthetic populations (one record per state)
# --------------------------------------------------------------------------
ZENODO_FILES_URL = "https://zenodo.org/records/{record}/files/{key}?download=1"

ZENODO_POPULATION_RECORDS = {
    "va": ("23088534", "10.5281/zenodo.23088534"),
    "ga": ("23088629", "10.5281/zenodo.23088629"),
    "ma": ("23088723", "10.5281/zenodo.23088723"),
    "mn": ("23088773", "10.5281/zenodo.23088773"),
    "wa": ("23088830", "10.5281/zenodo.23088830"),
    "ca": ("23088872", "10.5281/zenodo.23088872"),
}


# Original UVA Dataverse deposits. Kept for provenance and citation; not used
# for fetching (see the module docstring).
DATAVERSE_DOIS = {
    "va": "doi:10.18130/V3/5LSDCY",
    "ca": "doi:10.18130/V3/A7DQWM",
    "ga": "doi:10.18130/V3/4DWNRH",
    "ma": "doi:10.18130/V3/ZB0SGL",
    "mn": "doi:10.18130/V3/SB2PWT",
    "wa": "doi:10.18130/V3/PNGMRJ",
}
# Backwards-compatible alias.
STATE_DOIS = {k: f"doi:{v[1]}" for k, v in ZENODO_POPULATION_RECORDS.items()}


def population_record(state: str):
    """(record_id, doi) for a state's synthetic population on Zenodo."""
    return ZENODO_POPULATION_RECORDS.get(state.lower())


def population_mirror_url(state: str) -> "str | None":
    rec = population_record(state)
    return f"https://doi.org/{rec[1]}" if rec else None


def population_file_url(state: str, key: str) -> str:
    rid, _ = ZENODO_POPULATION_RECORDS[state.lower()]
    return ZENODO_FILES_URL.format(record=rid, key=key)


# The four files needed to build the demographics table used by the painter
# and by TwinSampler. Sizes are MB; checksums are MD5 from the Zenodo record.
CORE_FILES = {
    "va": {
        "va_persontrait_epihiper.txt.xz": (25.5, "13bcf3726b4b2d5ec23b5de5593ba5d2"),
        "va_person.csv.xz": (25.4, "9eea2e5dfc15711b7d5426d73048dd21"),
        "va_household.csv.xz": (20.2, "dc7aa659ca1512c2afc0c81eed2df103"),
        "va_residence_locations.csv.xz": (32.9, "cab93db949662cfa6789e39f58d583fb"),
    },
    "ca": {
        "ca_persontrait_epihiper.txt.xz": (106.5, "50e8826ee13b10bb0330a4cd1a7ff82a"),
        "ca_person.csv.xz": (107.3, "3b490678c54cb124c0cfe981cc860951"),
        "ca_household.csv.xz": (124.9, "6e2c0d5b681d7ee4be84bdebe49a318a"),
        "ca_residence_locations.csv.xz": (116.7, "75f4af815779cbc37e48bf17394c2dba"),
    },
    "ga": {
        "ga_persontrait_epihiper.txt.xz": (29.9, "eb34c2407c3c7e558994949dd8df8f5f"),
        "ga_person.csv.xz": (29.8, "86f6d0d78356e0fac5014548091b30b3"),
        "ga_household.csv.xz": (23.8, "467a26c9a1a53a387f1a7db45a145147"),
        "ga_residence_locations.csv.xz": (39.8, "f863e3a1d5b1c14a13878071af9aad01"),
    },
    "ma": {
        "ma_persontrait_epihiper.txt.xz": (20.0, "7c383b10311da52a0b00578ce6178d1d"),
        "ma_person.csv.xz": (20.2, "21e95d083c71f67e98f91398c07fed0f"),
        "ma_household.csv.xz": (16.5, "385ee72ec2e67e0fecdef966c68ff003"),
        "ma_residence_locations.csv.xz": (22.4, "377cf7e5e5cd2dd30a3c3b3684135131"),
    },
    "mn": {
        "mn_persontrait_epihiper.txt.xz": (16.8, "0c4d0ef351a7497460e7928a45426389"),
        "mn_person.csv.xz": (17.0, "c4cb647b074cbec026633b614752f459"),
        "mn_household.csv.xz": (13.9, "796b94c10b07039f7b89558807ccd662"),
        "mn_residence_locations.csv.xz": (27.2, "ec8d19e8df34da872e1fce05cb2232af"),
    },
    "wa": {
        "wa_persontrait_epihiper.txt.xz": (21.4, "5bb66bd43ea1f06512201c59a4be022c"),
        "wa_person.csv.xz": (21.5, "47b40bf915e13095c5d443f763835543"),
        "wa_household.csv.xz": (17.9, "dc7bbcdec87652926a123ec20ac85af9"),
        "wa_residence_locations.csv.xz": (30.2, "0cbf4215894cbac24c54421f19588ddf"),
    },
}

# EpiHiper *inputs*. PhyloGAS consumes EpiHiper OUTPUT (the replicates), so
# these are not needed unless you intend to run the ABM yourself. The contact
# networks alone are several GB, which is why they are opt-in.
EPIHIPER_INPUT_FILES = {
    "va": {
        "va_contact_network_epihiper.txt.xz": (907.2, "975f856e6f44f4e37345ebeaf470c253"),
        "va_persontrait_epihiper_db.tar.gz": (344.9, "b8cb3b2e8b455f8657706dfcf007e31c"),
    },
    "ca": {
        "ca_contact_network_epihiper.txt.xz": (4025.4, "97a64208da87450328034a6842b67644"),
        "ca_persontrait_epihiper_db.tar.gz": (834.0, "17fcb13d07df3ea9057f7c95c71ad5ff"),
    },
    "ga": {
        "ga_persontrait_epihiper_db.tar.gz": (385.3, "85173505a94141dcb5226fcd1d88ae30"),
        "ga_contact_network_epihiper.txt.xz": (1094.8, "6a79bac251c1212632db326a8c437872"),
    },
    "ma": {
        "ma_contact_network_epihiper.txt.xz": (719.9, "45d61c4fe699c3c520ab2f4c5183cbf2"),
        "ma_persontrait_epihiper_db.tar.gz": (334.7, "eef0f82a6d54eef04a388198e7301ff3"),
    },
    "mn": {
        "mn_persontrait_epihiper_db.tar.gz": (301.8, "426931d5c3e1bebfc6984e09aec2b674"),
        "mn_contact_network_epihiper.txt.xz": (589.4, "b65d3da62f983de54d5ecbbf50a1792c"),
    },
    "wa": {
        "wa_contact_network_epihiper.txt.xz": (788.5, "f8f4bd2d8c5180c2a16ef107b9de9537"),
        "wa_persontrait_epihiper_db.tar.gz": (329.8, "e0a786be5910c9dc1c0e4f808ea0a902"),
    },
}


def file_plan(states, with_epihiper_inputs: bool = False):
    """Return [(state, filename, url, size_mb, md5), ...] for the requested states."""
    plan = []
    for st in states:
        st = st.lower()
        if st not in CORE_FILES:
            raise KeyError(st)
        for name, (mb, md5) in CORE_FILES[st].items():
            plan.append((st, name, population_file_url(st, name), mb, md5))
        if with_epihiper_inputs:
            for name, (mb, md5) in EPIHIPER_INPUT_FILES.get(st, {}).items():
                plan.append((st, name, population_file_url(st, name), mb, md5))
    return plan


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


def _stream_to_file(url: str, dest: Path, timeout: int, label: str,
                    attempts: int = 4) -> None:
    """Stream a URL to ``dest``, retrying the WHOLE transfer on failure.

    _get_with_retry only covers establishing the connection. A drop part way
    through a 500 MB body would otherwise discard the partial file and abort,
    which matters because Dataverse is currently returning intermittent 503s
    (roughly half of requests as of 2026-09-30).
    """
    import time

    tmp = dest.with_suffix(dest.suffix + ".part")
    last = None
    for attempt in range(attempts):
        try:
            with _get_with_retry(url, timeout, attempts=2) as r:
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
            if total and done < total:
                raise IOError(f"truncated: got {done} of {total} bytes")
            tmp.rename(dest)
            return
        except Exception as exc:
            last = exc
            tmp.unlink(missing_ok=True)
            if attempt < attempts - 1:
                wait = 10 * (2 ** attempt)
                print(f"      transfer failed ({exc}); restarting in {wait}s "
                      f"(attempt {attempt + 2}/{attempts})", flush=True)
                time.sleep(wait)
    raise RuntimeError(f"{label}: failed after {attempts} attempts: {last}")


# --------------------------------------------------------------------------
# USDA Economic Research Service: rural-urban continuum codes
# --------------------------------------------------------------------------
# A pinned copy ships inside TwinSampler (linelist_generation/data/), which
# uses it by default; this URL is only for fetching a newer edition. USDA ERS works
# are US federal government publications and therefore public domain
# (17 U.S.C. 105). No checksum is pinned because the agency revises the file
# in place without versioning the URL.
RUCC_URL = ("https://www.ers.usda.gov/media/5768/"
            "2023-rural-urban-continuum-codes.csv?v=27068")
RUCC_FILENAME = "Ruralurbancontinuumcodes2023.csv"


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

    print(f"    downloading {dest.name} ...", flush=True)
    _stream_to_file(url, dest, timeout, dest.name)

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
