"""Canonical construction of the identifiers that link sequences to metadata.

Three strings tie the pipeline together, and every one of them is built from
the same two facts -- a person id and the tick they were exposed on:

``alias_pid``
    ``"{pid}.{tick}"``. The join key of record. Emitted by the painter in its
    metadata and by TwinSampler in the line list, so downstream joins are a
    key match rather than a parse of the strain id.

``strain``
    ``"{country}/{divisionAbbr}-EHip-{pid}.{tick}/{year}"``. Nextstrain's
    required identifier: it is the FASTA header the painter writes *and* the
    ``strain`` column ncov matches metadata on. If these two disagree by one
    character, ncov silently drops the sequence.

``year``
    Taken from the exposure date, which is ``start_date + (tick -
    start_tick)`` days. Both the painter and TwinSampler derive it that way.

This module is the canonical implementation for PhyloGAS. TwinSampler keeps
its own copy on purpose -- it must stay runnable without PhyloGAS installed --
but both are pinned to the same fixture table (``tests/data/strain_ids.json``
here, and the matching test in TwinSampler), so a change to either formula
fails a test instead of silently emptying a tree.

Why this exists: before it, ``alias_pid`` was built in three places with three
different type treatments (no cast, ``astype(str)``, ``int()``), and the
geography in ``strain`` came from config on one side and a hardcoded default
on the other.
"""

from __future__ import annotations

import datetime as _dt


def _as_int(value, what: str) -> int:
    """Coerce a pid or tick to int, refusing anything lossy.

    The whole point of this module is that ``123.0`` and ``123`` must not
    produce different identifiers, so a non-integral value is an error rather
    than something to round.
    """
    try:
        out = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{what} is not an integer: {value!r}") from exc
    if isinstance(value, float) and not float(value).is_integer():
        raise ValueError(f"{what} is not integral: {value!r}")
    if isinstance(value, str) and str(value).strip() != str(out):
        raise ValueError(f"{what} is not a plain integer string: {value!r}")
    return out


def alias_pid(pid, tick) -> str:
    """Return the ``"{pid}.{tick}"`` join key.

    >>> alias_pid(123, 45)
    '123.45'
    >>> alias_pid("123", "45")
    '123.45'
    >>> alias_pid(123.0, 45.0)
    '123.45'

    Float ticks used to produce ``'123.45.0'`` from one code path and
    ``'123.45'`` from another:

    >>> alias_pid(123, 45.5)
    Traceback (most recent call last):
        ...
    ValueError: tick is not integral: 45.5
    """
    return f"{_as_int(pid, 'pid')}.{_as_int(tick, 'tick')}"


def exposure_date(start_date, start_tick, tick) -> _dt.date:
    """The calendar date a tick corresponds to.

    >>> exposure_date("2021-04-07", 128, 128)
    datetime.date(2021, 4, 7)
    >>> exposure_date("2021-04-07", 128, 428)
    datetime.date(2022, 2, 1)
    """
    base = (start_date if isinstance(start_date, _dt.date)
            else _dt.date.fromisoformat(str(start_date)[:10]))
    return base + _dt.timedelta(days=_as_int(tick, "tick") - _as_int(start_tick, "start_tick"))


def exposure_year(start_date, start_tick, tick) -> int:
    """The year component of a strain id.

    >>> exposure_year("2021-04-07", 128, 128)
    2021
    >>> exposure_year("2021-04-07", 128, 428)
    2022
    """
    return exposure_date(start_date, start_tick, tick).year


def strain_id(country, division_abbr, pid, tick, year) -> str:
    """Return the Nextstrain ``strain`` identifier.

    ``year`` is passed in rather than derived so callers that already have a
    tick-to-year table do not rebuild it per record.

    >>> strain_id("USA", "VA", 123, 45, 2021)
    'USA/VA-EHip-123.45/2021'
    """
    return f"{country}/{division_abbr}-EHip-{alias_pid(pid, tick)}/{_as_int(year, 'year')}"


def strain_id_from_alias(country, division_abbr, alias, year) -> str:
    """``strain_id`` for a caller that already holds the ``alias_pid``.

    The painter builds the alias once per record and reuses it, so this avoids
    splitting it back apart.

    >>> strain_id_from_alias("USA", "VA", "123.45", 2021)
    'USA/VA-EHip-123.45/2021'
    """
    return f"{country}/{division_abbr}-EHip-{alias}/{_as_int(year, 'year')}"


def strain_id_from_tick(country, division_abbr, pid, tick, start_date, start_tick) -> str:
    """``strain_id`` with the year derived from the tick.

    >>> strain_id_from_tick("USA", "VA", 123, 428, "2021-04-07", 128)
    'USA/VA-EHip-123.428/2022'
    """
    return strain_id(country, division_abbr, pid, tick,
                     exposure_year(start_date, start_tick, tick))
