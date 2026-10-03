"""The strain id is a cross-repo contract; pin it.

The painter writes `strain` as a FASTA header. TwinSampler writes the same
string into the line list's `strain` column, from its own implementation.
ncov joins sequences to metadata on that string, so a one-character
disagreement empties the tree -- sometimes without an error.

TwinSampler has a mirror of this test reading the same fixture. Change one
formula and at least one of the two fails.
"""

import json
from pathlib import Path

import pytest

from phylogas import ids

FIXTURE = Path(__file__).parent / "data" / "strain_ids.json"


def _fixture():
    d = json.loads(FIXTURE.read_text())
    return d, d["cases"]


def test_fixture_is_readable():
    d, cases = _fixture()
    assert cases, "fixture has no cases"
    assert d["start_date"] and d["start_tick"] is not None


@pytest.mark.parametrize("case", _fixture()[1], ids=lambda c: c["strain"])
def test_alias_pid_matches_fixture(case):
    assert ids.alias_pid(case["pid"], case["tick"]) == case["alias_pid"]


@pytest.mark.parametrize("case", _fixture()[1], ids=lambda c: c["strain"])
def test_strain_id_matches_fixture(case):
    assert ids.strain_id(
        case["country"], case["divisionAbbr"], case["pid"], case["tick"], case["year"]
    ) == case["strain"]


@pytest.mark.parametrize("case", _fixture()[1], ids=lambda c: c["strain"])
def test_year_derives_from_the_tick(case):
    d, _ = _fixture()
    assert ids.exposure_year(d["start_date"], d["start_tick"], case["tick"]) == case["year"]
    assert ids.strain_id_from_tick(
        case["country"], case["divisionAbbr"], case["pid"], case["tick"],
        d["start_date"], d["start_tick"],
    ) == case["strain"]


def test_integral_floats_and_strings_agree_with_ints():
    """The three old implementations differed exactly here.

    One used no cast, one `astype(str)`, one `int()`. A float tick gave
    '123.45.0' from one and '123.45' from another.
    """
    assert ids.alias_pid(123, 45) == ids.alias_pid(123.0, 45.0) == ids.alias_pid("123", "45")


@pytest.mark.parametrize("pid,tick", [(123, 45.5), (1.5, 10), ("12x", 3)])
def test_non_integral_input_is_refused(pid, tick):
    with pytest.raises(ValueError):
        ids.alias_pid(pid, tick)
