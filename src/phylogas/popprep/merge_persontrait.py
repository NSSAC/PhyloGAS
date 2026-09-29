#!/usr/bin/env python3
"""Build the demographics table that the Genetic Painter joins into metadata.

Rescued from `synthetic_biosurveillance/data/merge_persontrait_with_person.py`
and rewritten for the v2.4.0 synthetic population schema published on UVA
Dataverse. The original targeted v1.9.0 and no longer matches the deposited
files.

Why this exists
---------------
The file the pipeline configures as ``demographics_file`` is NOT the EpiHiper
persontrait file, despite frequently being called one. It is a derived join.
The painter's ``--add_metadata`` expects::

    gender, county, home_latitude, home_longitude, latino, race, smh_race, age_group

but the v2.4.0 persontrait file provides only ``gender``, ``race``,
``smh_race`` and ``age_group``. The rest have to be recovered:

* ``county``        - persontrait has ``county_fips``; mapped via a FIPS table.
* ``latino``        - derived from the person file's ``hispanic`` code.
* ``home_latitude`` - not present in v2.4.0 at all. Reached by joining
  ``home_longitude``   persontrait -> household (on ``hid``)
                       -> residence_locations (on ``rlid``).

v1.9.0 carried ``home_latitude``/``home_longitude`` directly on the persontrait
row, which is why the original script never needed the household join.

The persontrait file also opens with a one-line EpiHiper JSON schema header,
which must be skipped or pandas treats it as the column row.

Usage
-----
    python -m phylogas.popprep.merge_persontrait \\
        --persontrait va_persontrait_epihiper.txt.xz \\
        --person      va_person.csv.xz \\
        --household   va_household.csv.xz \\
        --residence   va_residence_locations.csv.xz \\
        --fips        county_fips.csv \\
        --out         va_2_4_0_demographics.csv

``--household`` and ``--residence`` may be omitted, in which case the
coordinate columns are emitted empty and a warning is printed.
All inputs may be plain or ``.xz``/``.gz`` compressed; pandas infers from the
extension.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Census race codes used by the synthetic population.
RACE_LOOKUP = {
    1: "White alone",
    2: "Black alone",
    3: "Native American alone",
    4: "Alaska Native alone",
    5: "Native American and Alaska Native, tribe specified",
    6: "Asian alone",
    7: "Native Hawaiian and Other Pacific Islander alone",
    8: "Some other race alone",
    9: "Two or more races",
}

AGE_GROUP_LOOKUP = {
    "p": "Preschool (0-4)",
    "s": "Student (5-17)",
    "a": "Adult (18-49)",
    "o": "Older adult (50-64)",
    "g": "Senior (65+)",
}

# Columns the painter's --add_metadata expects to find.
REQUIRED_OUTPUT_COLUMNS = [
    "pid", "gender", "county", "home_latitude", "home_longitude",
    "latino", "race", "smh_race", "age_group",
]


def calculate_smh_race(row) -> str:
    """Scenario Modeling Hub race categorisation.

    Latino takes precedence over race code, matching the original
    synthetic_biosurveillance implementation.
    """
    try:
        if int(row["hispanic_code"]) > 1:
            return "Latino"
        rc = int(row["race_code"])
    except (ValueError, TypeError, KeyError):
        return "Other"
    if rc == 1:
        return "White"
    if rc == 2:
        return "Black"
    if rc == 6:
        return "Asian"
    return "Other"


def _read(path: str | Path, skiprows: int = 0, usecols=None) -> pd.DataFrame:
    """Read a CSV that may be plain, .xz or .gz (pandas infers by extension)."""
    p = Path(path)
    if not p.exists():
        raise SystemExit(f"ERROR: input not found: {p}")
    try:
        return pd.read_csv(p, skiprows=skiprows, usecols=usecols)
    except ValueError as exc:
        raise SystemExit(f"ERROR: could not parse {p}: {exc}") from exc


def _detect_schema_header(path: str | Path) -> int:
    """Return 1 if the file opens with the EpiHiper JSON schema line, else 0.

    Detected rather than assumed: some locally-prepared persontrait files have
    already had the line stripped, and skipping unconditionally would silently
    discard the real header.
    """
    import gzip
    import lzma

    p = str(path)
    opener = lzma.open if p.endswith(".xz") else gzip.open if p.endswith(".gz") else open
    try:
        with opener(p, "rt") as fh:
            first = fh.readline().lstrip()
    except OSError as exc:
        raise SystemExit(f"ERROR: could not open {p}: {exc}") from exc
    return 1 if first.startswith("{") else 0


def build_demographics(
    persontrait: str,
    person: str,
    fips: str,
    out: str,
    household: str | None = None,
    residence: str | None = None,
) -> pd.DataFrame:
    """Join the population files into the painter's demographics table."""
    skip = _detect_schema_header(persontrait)
    if skip:
        print("  Detected EpiHiper JSON schema header; skipping first line.")
    pt = _read(persontrait, skiprows=skip)
    print(f"  persontrait: {len(pt):,} rows, columns={list(pt.columns)}")

    if "pid" not in pt.columns:
        raise SystemExit(
            f"ERROR: no 'pid' column in {persontrait} after skipping {skip} line(s).\n"
            f"       Found: {list(pt.columns)[:8]}"
        )

    # --- person file: race + hispanic --------------------------------------
    pf = _read(person, usecols=lambda c: c in {"pid", "race", "hispanic"})
    pf = pf.rename(columns={"race": "race_code", "hispanic": "hispanic_code"})
    pt = pt.drop(columns=[c for c in ("race", "hispanic", "smh_race") if c in pt.columns])
    pt = pt.merge(pf, on="pid", how="inner")
    print(f"  + person      -> {len(pt):,} rows")

    # --- derived demographic columns ---------------------------------------
    pt["smh_race"] = pt.apply(calculate_smh_race, axis=1)
    pt["latino"] = np.where(pd.to_numeric(pt["hispanic_code"], errors="coerce") > 1, True, False)
    pt["race"] = pd.to_numeric(pt["race_code"], errors="coerce").map(RACE_LOOKUP)

    if "age_group" in pt.columns:
        pt = pt.rename(columns={"age_group": "age_group_code"})
        pt["age_group"] = pt["age_group_code"].map(AGE_GROUP_LOOKUP)

    # --- county name from FIPS ---------------------------------------------
    if "county_fips" in pt.columns:
        fips_df = _read(fips, usecols=lambda c: c.lower() in {"fips", "county"})
        fips_df.columns = [c.lower() for c in fips_df.columns]
        if {"fips", "county"} <= set(fips_df.columns):
            fips_df["fips"] = fips_df["fips"].astype(str).str.zfill(5)
            pt["_fips5"] = pt["county_fips"].astype(str).str.zfill(5)
            pt = pt.merge(
                fips_df.rename(columns={"fips": "_fips5"}), on="_fips5", how="left"
            )
            missing = pt["county"].isna().sum()
            if missing:
                print(f"  WARNING: {missing:,} rows had no county name for their FIPS code")
            pt = pt.drop(columns=["_fips5"])
            print(f"  + county names from FIPS")
        else:
            print(f"  WARNING: {fips} lacks 'FIPS'/'county' columns; county left empty")
            pt["county"] = pd.NA
    else:
        pt["county"] = pd.NA

    # --- coordinates: persontrait -> household -> residence_locations -------
    # v2.4.0 dropped home_latitude/home_longitude from the persontrait file.
    if "home_latitude" in pt.columns and "home_longitude" in pt.columns:
        print("  Coordinates already present on persontrait (v1.9.0-style schema)")
    elif household and residence:
        hh = _read(household, usecols=lambda c: c in {"hid", "rlid"})
        rl = _read(residence, usecols=lambda c: c in {"rlid", "latitude", "longitude"})
        if "hid" not in pt.columns:
            print("  WARNING: persontrait has no 'hid'; cannot join coordinates")
            pt["home_latitude"] = pd.NA
            pt["home_longitude"] = pd.NA
        else:
            before = len(pt)
            pt = pt.merge(hh, on="hid", how="left").merge(rl, on="rlid", how="left")
            pt = pt.rename(columns={"latitude": "home_latitude", "longitude": "home_longitude"})
            got = pt["home_latitude"].notna().sum()
            print(f"  + household -> residence_locations: {got:,}/{before:,} rows got coordinates")
            if "rlid" in pt.columns:
                pt = pt.drop(columns=["rlid"])
    else:
        print(
            "  WARNING: --household/--residence not supplied. v2.4.0 persontrait\n"
            "           has no coordinates, so home_latitude/home_longitude will be empty."
        )
        pt["home_latitude"] = pd.NA
        pt["home_longitude"] = pd.NA

    # --- tidy ----------------------------------------------------------------
    pt = pt.drop(
        columns=[c for c in ("race_code", "hispanic_code", "age_group_code", "county_fips")
                 if c in pt.columns]
    )

    for col in REQUIRED_OUTPUT_COLUMNS:
        if col not in pt.columns:
            print(f"  WARNING: expected output column '{col}' could not be produced")
            pt[col] = pd.NA

    ordered = REQUIRED_OUTPUT_COLUMNS + [c for c in pt.columns if c not in REQUIRED_OUTPUT_COLUMNS]
    pt = pt[ordered]

    Path(out).parent.mkdir(parents=True, exist_ok=True)
    pt.to_csv(out, index=False)
    print(f"  Wrote {len(pt):,} rows x {len(pt.columns)} cols -> {out}")
    return pt


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        prog="phylogas build-demographics",
        description="Join the synthetic population files into the painter's demographics table.",
    )
    ap.add_argument("--persontrait", required=True, help="<state>_persontrait_epihiper.txt[.xz]")
    ap.add_argument("--person", required=True, help="<state>_person.csv[.xz]")
    ap.add_argument("--fips", required=True, help="county FIPS -> name lookup CSV")
    ap.add_argument("--out", required=True, help="output demographics CSV")
    ap.add_argument("--household", default=None, help="<state>_household.csv[.xz] (for coordinates)")
    ap.add_argument("--residence", default=None, help="<state>_residence_locations.csv[.xz]")
    args = ap.parse_args(argv)

    build_demographics(
        persontrait=args.persontrait,
        person=args.person,
        fips=args.fips,
        out=args.out,
        household=args.household,
        residence=args.residence,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
