"""Seed screening: does it reject a genome that would poison a whole chain?

The case these guard is Virginia's second importation,
`USA/VA-CDC-LC0041349/2021` -- a mixed Alpha/Delta assembly that founded the
giant transmission component, after which the ncov build's diagnostics
excluded 71% of that state's sampled genomes.

    pytest tests/test_seed_qc.py -v
"""

from __future__ import annotations

import pytest

qc = pytest.importorskip("phylogas.seqprep.seed_qc")


REF = "ACGT" * 7475 + "ACG"          # 29,903 bases, like an aligned genome


def _with(bases: dict, length: int = 29903, fill: str = "A") -> str:
    seq = list(fill * length)
    for pos, base in bases.items():
        seq[pos - 1] = base
    return "".join(seq)


# --------------------------------------------------------------------------
# tool-free screens
# --------------------------------------------------------------------------
def test_ambiguity_counts_non_acgt():
    assert qc.ambiguity("ACGT") == 0.0
    assert qc.ambiguity("ACGN") == pytest.approx(0.25)
    assert qc.ambiguity("NNNN") == 1.0
    assert qc.ambiguity("acgt") == 0.0            # case-insensitive
    assert qc.ambiguity("") == 1.0


def test_terminal_missing_finds_a_lost_end():
    """A genome can pass on overall ambiguity while an entire end is gone --
    the Virginia seed was 47% ambiguous across its 3' 2kb."""
    good = "A" * 29903
    assert qc.terminal_missing(good) == 0.0
    # Passes overall ambiguity, yet most of the 3' window is gone -- the shape
    # of the Virginia seed, which was 47% ambiguous across its 3' 2kb.
    lost_3prime = "A" * 28503 + "N" * 1400
    assert qc.ambiguity(lost_3prime) < qc.MAX_AMBIGUOUS
    assert qc.terminal_missing(lost_3prime) == pytest.approx(0.70)
    lost_5prime = "N" * 1400 + "A" * 28503
    assert qc.terminal_missing(lost_5prime) == pytest.approx(0.70)


def test_basic_screen_rejects_and_reports_a_reason():
    recs = {"clean": "A" * 29903,
            "ambiguous": "N" * 3000 + "A" * 26903,
            # under the 5% ambiguity bar, but 70% of its 3' window is missing
            "truncated": "A" * 28503 + "N" * 1400}
    fails = qc.basic_failures(recs)
    assert "clean" not in fails
    assert "ambiguous" in fails and "ambiguous bases" in fails["ambiguous"]
    assert "truncated" in fails and "terminal" in fails["truncated"]


def test_consensus_outliers_need_a_population_and_equal_lengths():
    assert qc.consensus_outliers({"a": "ACGT", "b": "ACGT"}) == {}      # too few
    recs = {f"s{i}": "A" * 500 for i in range(60)}
    recs["ragged"] = "A" * 499
    assert qc.consensus_outliers(recs) == {}                            # unaligned


def test_consensus_outlier_flags_a_divergent_genome():
    recs = {f"s{i}": _with({}, 500) for i in range(80)}
    recs["odd"] = _with({p: "C" for p in range(10, 200, 4)}, 500)       # ~48 mismatches
    out = qc.consensus_outliers(recs)
    assert "odd" in out and out["odd"] > 10
    assert not any(k.startswith("s") for k in out)


# --------------------------------------------------------------------------
# lineage membership (pango_aliasor, not a table of clades)
# --------------------------------------------------------------------------
def fake_base_map(bases, recombinant=False):
    """Stand-in for pango_aliasor: AY.* and XBB.* roll up to their ancestor."""
    m = {b: b for b in bases}
    if "B.1.617.2" in bases:
        m.update({f"AY.{i}": "B.1.617.2" for i in range(1, 200)})
        m["B.1.617.2.1"] = "B.1.617.2"
    if "B.1.1.7" in bases:
        m["Q.1"] = "B.1.1.7"
    if "XBB" in bases:
        m.update({"XBB.1.5": "XBB", "XBB.1.16": "XBB"})
    return m


@pytest.mark.parametrize("called,expected", [
    ("B.1.617.2", True),          # the lineage itself
    ("AY.103", True),             # a sublineage, rolled up by the aliasor
    ("AY.127", True),
    ("B.1.1.7", False),           # Alpha: the Virginia seed's backbone
    ("20I", False),               # a clade name is not a lineage call
    ("", None),                   # unassigned: never a rejection
    ("?", None),
    ("recombinant", None),
])
def test_lineage_matcher_rolls_calls_up(called, expected):
    inside = qc.lineage_matcher(["B.1.617.2"], base_map_fn=fake_base_map)
    assert inside(called) is expected


def test_lineage_matcher_handles_any_lineage_without_a_table():
    """The case a hardcoded clade table could not cover: recombinants."""
    inside = qc.lineage_matcher(["XBB"], base_map_fn=fake_base_map)
    assert inside("XBB.1.5") is True
    assert inside("B.1.617.2") is False


def test_lineage_matcher_strips_the_query_star():
    inside = qc.lineage_matcher(["B.1.617.2*"], base_map_fn=fake_base_map)
    assert inside("AY.4") is True


def test_lineage_matcher_declines_to_judge_without_an_aliasor():
    """No pango_aliasor: unknown calls must not be rejected on a guess."""
    inside = qc.lineage_matcher(["B.1.617.2"], base_map_fn=None)
    assert inside("B.1.617.2") is True        # exact match still decidable
    assert inside("AY.103") is None
    assert inside("B.1.1.7") is None


def test_lineage_matcher_survives_a_broken_aliasor():
    def boom(bases, recombinant=False):
        raise RuntimeError("pango_aliasor exploded")
    inside = qc.lineage_matcher(["B.1.617.2"], base_map_fn=boom)
    assert inside("AY.103") is None


def test_lineage_matcher_accepts_an_aliasor_without_the_recombinant_flag():
    inside = qc.lineage_matcher(["B.1.617.2"],
                                base_map_fn=lambda bases: {"AY.4": "B.1.617.2"})
    assert inside("AY.4") is True


def test_no_lineages_requested_means_no_check():
    assert qc.lineage_matcher([], base_map_fn=fake_base_map) is None


# --------------------------------------------------------------------------
# dispatch
# --------------------------------------------------------------------------
def test_off_mode_screens_nothing():
    recs = {"ambiguous": "N" * 29903}
    fails, status = qc.screen(recs, ["B.1.617.2"], mode="off")
    assert fails == {} and "disabled" in status


def test_nextclade_mode_falls_back_when_the_tool_is_missing(tmp_path, monkeypatch):
    import os
    monkeypatch.setattr(qc.shutil, "which", lambda _n: None)
    monkeypatch.setattr(os, "environ", {"NEXTSTRAIN_HOME": str(tmp_path / "absent")})
    recs = {"clean": "A" * 29903, "ambiguous": "N" * 3000 + "A" * 26903}
    fails, status = qc.screen(recs, ["B.1.617.2"], mode="nextclade", ncov_dir=tmp_path)
    assert status.startswith("WARNING")
    assert "neither on PATH nor in a Nextstrain conda runtime" in status
    assert "cannot detect mixed-lineage" in status       # the fallback's real limit
    assert "ambiguous" in fails and "clean" not in fails


def test_nextclade_needs_its_dataset(tmp_path, monkeypatch):
    monkeypatch.setattr(qc.shutil, "which", lambda n: "/usr/bin/" + n)
    fails, status = qc.nextclade_failures({"a": "ACGT"}, ["B.1.617.2"], ncov_dir=tmp_path)
    assert fails == {} and "dataset" in status


def test_nextclade_is_found_on_path_or_in_the_runtime(tmp_path, monkeypatch):
    """`--exec` is not usable: it exists only in newer CLI releases, and an
    older one forwards it to Snakemake as --executor. So the tool is located
    by path first; `test_runtime_pipes_*` covers the route taken when there
    is no path to find."""
    import os, stat
    monkeypatch.setattr(qc.shutil, "which", lambda n: "/usr/bin/nextclade"
                        if n == "nextclade" else None)
    assert qc.nextclade_command() == ["/usr/bin/nextclade"]

    monkeypatch.setattr(qc.shutil, "which", lambda _n: None)
    binp = tmp_path / "runtimes" / "conda" / "env" / "bin"
    binp.mkdir(parents=True)
    tool = binp / "nextclade"
    tool.write_text("#!/bin/sh\n")
    tool.chmod(tool.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setattr(os, "environ", {"NEXTSTRAIN_HOME": str(tmp_path)})
    assert qc.nextclade_command() == [str(tool)]

    monkeypatch.setattr(os, "environ", {"NEXTSTRAIN_HOME": str(tmp_path / "absent")})
    assert qc.nextclade_command() is None


def test_nextclade_parses_its_tsv_and_applies_ncov_criteria(tmp_path, monkeypatch):
    """Criteria mirror ncov scripts/diagnostic.py, so a seed that passes here
    cannot be excluded later for a reason it already carried."""
    dataset = tmp_path / "data" / "sars-cov-2-nextclade-defaults.zip"
    dataset.parent.mkdir(parents=True)
    dataset.write_text("")
    out = tmp_path / "nextclade.tsv"

    rows = [
        # seqName, Nextclade_pango, reversions, labeled, mixed-sites status
        ("good", "AY.103", "1", "1", "good"),          # a Delta sublineage
        ("contaminated", "AY.103", "7", "3", "good"),  # 10 > 5
        ("alpha", "B.1.1.7", "0", "0", "good"),        # clean QC, wrong lineage
        ("mixed", "AY.103", "0", "0", "bad"),
        ("uncalled", "?", "0", "0", "good"),           # never rejected on doubt
    ]

    def fake_run(cmd, **kw):
        header = ("seqName\tNextclade_pango\t"
                  "privateNucMutations.reversionSubstitutions\t"
                  "privateNucMutations.labeledSubstitutions\tqc.mixedSites.status\n")
        out.write_text(header + "".join("\t".join(r) + "\n" for r in rows))
        return type("P", (), {"returncode": 0, "stdout": "", "stderr": ""})()

    monkeypatch.setattr(qc.shutil, "which", lambda n: "/usr/bin/" + n)
    monkeypatch.setattr(qc.subprocess, "run", fake_run)
    recs = {n: "A" * 29903 for n, *_ in rows}
    fails, status = qc.nextclade_failures(recs, ["B.1.617.2"], ncov_dir=tmp_path,
                                          workdir=tmp_path, base_map_fn=fake_base_map)
    assert "nextclade screened 5" in status
    assert "good" not in fails
    assert "uncalled" not in fails                  # doubt is not a rejection
    assert "reversions+contaminants" in fails["contaminated"]
    assert "is outside B.1.617.2" in fails["alpha"]
    assert "mixed sites" in fails["mixed"]


def test_report_is_written_with_one_row_per_candidate(tmp_path):
    path = tmp_path / "VA_seed_qc.csv"
    qc.write_report(path, [
        dict(cluster_id="c1", attempt=1, strain="bad", verdict="reject", reason="why"),
        dict(cluster_id="c1", attempt=2, strain="ok", verdict="accept", reason=""),
    ])
    text = path.read_text().splitlines()
    assert text[0].startswith("cluster_id,attempt,strain,verdict,reason")
    assert len(text) == 3
    assert "reject" in text[1] and "accept" in text[2]


def test_early_importation_sort_does_not_compare_dicts():
    """Two rejects at the same position must not fall through to comparing
    the dicts themselves -- that raised TypeError mid-run."""
    order = {"c1": 0, "c2": 0}
    rejected = [dict(cluster_id="c1", strain="a"), dict(cluster_id="c2", strain="b")]
    out = sorted(((order.get(r["cluster_id"], 10 ** 9), r) for r in rejected),
                 key=lambda t: t[0])
    assert [r["strain"] for _, r in out] == ["a", "b"]


# --------------------------------------------------------------------------
# reaching the runtime when there is no binary to point at
# --------------------------------------------------------------------------
@pytest.fixture
def runtime():
    return pytest.importorskip("phylogas.nextstrain_runtime")


def test_runtime_prefers_a_direct_binary(runtime, monkeypatch):
    seen = {}

    def fake_run(cmd, **kw):
        seen["cmd"], seen["kw"] = cmd, kw
        return type("P", (), {"returncode": 0, "stdout": "", "stderr": ""})()

    monkeypatch.setattr(runtime.shutil, "which", lambda n: "/usr/bin/" + n)
    monkeypatch.setattr(qc.subprocess, "run", fake_run)
    proc, source = runtime.run("augur", ["refine", "--tree", "t.nwk"],
                               ncov_dir="/some/ncov")
    assert proc.returncode == 0 and source == "PATH"
    assert seen["cmd"] == ["/usr/bin/augur", "refine", "--tree", "t.nwk"]
    assert "nextstrain" not in " ".join(seen["cmd"])


def test_runtime_pipes_into_nextstrain_shell_when_no_binary_is_reachable(
        runtime, tmp_path, monkeypatch):
    """The docker and singularity runtimes have no bin directory on the host,
    so the only route is the CLI. `nextstrain shell` is not interactive-only:
    it runs what arrives on stdin."""
    import os
    seen = {}

    def fake_run(cmd, **kw):
        seen["cmd"], seen["kw"] = cmd, kw
        return type("P", (), {"returncode": 0, "stdout": "augur 34.1.4", "stderr": ""})()

    monkeypatch.setattr(runtime.shutil, "which",
                        lambda n: "/usr/bin/nextstrain" if n == "nextstrain" else None)
    monkeypatch.setattr(os, "environ", {"NEXTSTRAIN_HOME": str(tmp_path / "absent")})
    monkeypatch.setattr(qc.subprocess, "run", fake_run)

    proc, source = runtime.run("augur", ["refine", "--tree", "a b.nwk"],
                               ncov_dir=tmp_path)
    assert proc.returncode == 0 and source == "nextstrain shell"
    assert seen["cmd"] == ["nextstrain", "shell", str(tmp_path)]
    # The arguments become a shell line, so a path with a space must survive.
    assert seen["kw"]["input"].strip() == "augur refine --tree 'a b.nwk'"


def test_runtime_reports_why_it_could_not_run_anything(runtime, tmp_path, monkeypatch):
    import os
    monkeypatch.setattr(runtime.shutil, "which", lambda _n: None)
    monkeypatch.setattr(os, "environ", {"NEXTSTRAIN_HOME": str(tmp_path / "absent")})
    proc, reason = runtime.run("augur", ["--version"], ncov_dir=tmp_path)
    assert proc is None
    assert "neither on PATH nor in a Nextstrain conda runtime" in reason
    assert "Nextstrain CLI is not available" in reason


def test_runtime_will_not_pipe_without_a_checkout_to_enter(runtime, tmp_path,
                                                           monkeypatch):
    """`nextstrain shell` needs a build directory; with none to name, say so
    rather than guessing one."""
    import os
    monkeypatch.setattr(runtime.shutil, "which",
                        lambda n: "/usr/bin/nextstrain" if n == "nextstrain" else None)
    monkeypatch.setattr(os, "environ", {"NEXTSTRAIN_HOME": str(tmp_path / "absent")})
    proc, reason = runtime.run("augur", ["--version"])
    assert proc is None and "neither on PATH" in reason
