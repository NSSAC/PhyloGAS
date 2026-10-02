"""Smoke tests: can a fresh install actually run?

Deliberately narrow. These check wiring -- imports resolve, entry points exist,
config expands, the mutation engine produces sane output -- not science. They
should pass on any machine with the package installed and no data downloaded.

    pytest tests/ -v
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest


# --------------------------------------------------------------------------
# imports
# --------------------------------------------------------------------------
def test_package_imports():
    import phylogas
    assert phylogas.__version__


@pytest.mark.parametrize("mod", [
    "phylogas.cli",
    "phylogas.config",
    "phylogas.sources",
    "phylogas.benchmark.mugration",
    "phylogas.benchmark.scoring",
    "phylogas.benchmark.runner",
    "phylogas.popprep.merge_persontrait",
])
def test_submodules_import(mod):
    __import__(mod)


def test_mutational_models_registered():
    from phylogas.painter.mutational_models import registry
    for name in ("rate_limited", "simple", "poor"):
        assert registry[name] is not None


# --------------------------------------------------------------------------
# CLI surface
# --------------------------------------------------------------------------
def test_phylogas_on_path():
    assert shutil.which("phylogas"), "phylogas not on PATH; was the package installed?"


@pytest.mark.parametrize("cmd", [
    "status", "train", "paint", "prep-seeds", "subset-fasta", "run",
    "benchmark", "compare-strategies", "validate-config", "fetch-data",
    "build-demographics",
])
def test_subcommand_help(cmd):
    r = subprocess.run(["phylogas", cmd, "--help"], capture_output=True, text=True)
    assert r.returncode == 0, f"`phylogas {cmd} --help` failed:\n{r.stdout}\n{r.stderr}"


@pytest.mark.parametrize("what", ["mugration", "sequence", "compare"])
def test_benchmark_subcommand_help(what):
    r = subprocess.run(["phylogas", "benchmark", what, "--help"],
                       capture_output=True, text=True)
    assert r.returncode == 0


# --------------------------------------------------------------------------
# config
# --------------------------------------------------------------------------
def test_config_template_loads_and_expands(tmp_path):
    from pathlib import Path

    from phylogas.config import Config

    template = Path(__file__).resolve().parents[1] / "config.template.yaml"
    assert template.is_file(), "config.template.yaml is missing from the repo root"

    cfg = Config.load(template)
    assert cfg.get("project_name")

    # {data_dir} and {population.state} must be resolved, not left literal.
    out = str(cfg.get("genetic_painter.output_prefix"))
    assert "{" not in out, f"unexpanded placeholder in output_prefix: {out}"
    epihiper = str(cfg.get("epihiper.output_csv"))
    assert "{" not in epihiper, f"unexpanded placeholder in output_csv: {epihiper}"


def test_config_missing_key_raises():
    from phylogas.config import Config, ConfigError

    cfg = Config({"a": {"b": 1}})
    assert cfg.get("a.b") == 1
    assert cfg.get("a.nope", default="fallback") == "fallback"
    with pytest.raises(ConfigError):
        cfg.get("a.nope")


# --------------------------------------------------------------------------
# data source tables
# --------------------------------------------------------------------------
def test_source_tables_consistent():
    from phylogas import sources

    for state, files in sources.CORE_FILES.items():
        assert state in sources.ZENODO_POPULATION_RECORDS, f"{state} has files but no record"
        # Every state needs all four files to build a demographics table.
        kinds = {"persontrait", "person.csv", "household", "residence_locations"}
        for kind in kinds:
            assert any(kind in f for f in files), f"{state} missing {kind}"

    for state in sources.SIMULATION_FILES:
        assert state in sources.CORE_FILES, f"{state} has a replicate but no population"


def test_population_urls_well_formed():
    """Every core file must resolve to a Zenodo download URL with a checksum."""
    from phylogas import sources

    for state in sources.CORE_FILES:
        rid, doi = sources.population_record(state)
        assert doi.startswith("10.5281/zenodo."), f"{state}: {doi}"
        for name, (mb, md5) in sources.CORE_FILES[state].items():
            url = sources.population_file_url(state, name)
            assert url.startswith(f"https://zenodo.org/records/{rid}/files/")
            assert name in url
            assert len(md5) == 32, f"{state}/{name}: md5 looks wrong ({md5!r})"
            assert mb > 0


def test_file_plan_shape():
    """file_plan must yield (state, name, url, size_mb, md5)."""
    from phylogas import sources

    plan = sources.file_plan(["va"])
    assert len(plan) == 4
    for st, name, url, mb, md5 in plan:
        assert st == "va" and name.startswith("va_")
        assert url.startswith("https://zenodo.org/")
        assert len(md5) == 32


def test_simulation_urls_well_formed():
    from phylogas import sources

    for state in sources.SIMULATION_FILES:
        url = sources.simulation_url(state)
        assert url.startswith("https://zenodo.org/records/")
        assert state in url


# --------------------------------------------------------------------------
# mutation engine
# --------------------------------------------------------------------------
def test_rate_limited_model_produces_valid_sequences():
    import numpy as np

    from phylogas.painter.mutational_models import registry

    rng = np.random.default_rng(0)
    n = 2000
    # Mostly conserved, with a mutable tail -- mirrors a real threshold file.
    thresholds = np.full(n, 100.0)
    thresholds[:200] = 70.0

    letters = np.array(list("ACGTNRKSYMWBHDV"))
    prob = np.zeros((len(letters), n))
    prob[:4, :] = 0.25

    model = registry["rate_limited"](10.0, thresholds, prob, letters)
    seq = rng.choice(list("ACGT"), n).astype("S1")

    mutated = 0
    for _ in range(300):
        out = model.mutate(seq)
        assert len(out) == n
        diff = np.where(out != seq)[0]
        mutated += len(diff) > 0
        for i in diff:
            assert out[i] in (b"A", b"C", b"G", b"T"), "non-canonical base produced"
            assert out[i] != seq[i], "a 'mutation' left the base unchanged"
            assert thresholds[i] < 100.0, "mutation landed on a fully conserved site"

    assert 0 < mutated < 300, f"implausible mutation frequency: {mutated}/300"


def test_scoring_identical_matrices_scores_perfectly():
    import numpy as np

    from phylogas.benchmark.scoring import score_offdiagonals

    v = np.array([0.0, 0.2, 0.5, 0.0, 0.3])
    s = score_offdiagonals(v, v)
    assert s["cosine_similarity"] == pytest.approx(1.0)
    assert s["topological_f1"] == pytest.approx(1.0)
    assert s["masked_mae"] == pytest.approx(0.0)


def test_scoring_rejects_mismatched_shapes():
    import numpy as np

    from phylogas.benchmark.scoring import score_offdiagonals

    with pytest.raises(ValueError):
        score_offdiagonals(np.zeros(4), np.zeros(5))


# --------------------------------------------------------------------------
# optional integrations
# --------------------------------------------------------------------------
def test_beyondbaseline_available():
    """PhyloGAS hard-depends on BeyondBaseline for sample selection."""
    pytest.importorskip("scenarios_simulation")
    assert shutil.which("beyond-baseline-sweep"), \
        "BeyondBaseline is importable but its CLI is not on PATH"


def test_snakefile_present_and_parses():
    from pathlib import Path

    sf = Path(__file__).resolve().parents[1] / "Snakefile"
    assert sf.is_file()
    text = sf.read_text()
    for rule in ("rule all", "rule train_entropy", "rule paint_network"):
        assert rule in text, f"{rule} missing from the Snakefile"


# --------------------------------------------------------------------------
# compression-extension tolerance
# --------------------------------------------------------------------------
def test_resolve_variant_handles_extension_drift(tmp_path):
    """Inputs copied from a cluster are often decompressed, while the config
    names the .xz form (or vice versa). Both must be found."""
    from phylogas.config import Config

    plain = tmp_path / "va_person.csv"
    plain.write_text("pid\n1\n")
    # config says .xz, disk has plain
    assert Config.resolve_variant(str(plain) + ".xz") == plain
    # exact match
    assert Config.resolve_variant(str(plain)) == plain

    comp = tmp_path / "va_household.csv.xz"
    comp.write_bytes(b"\xfd7zXZ\x00")
    # config says plain, disk has .xz
    assert Config.resolve_variant(str(tmp_path / "va_household.csv")) == comp

    # genuinely absent stays absent
    assert Config.resolve_variant(str(tmp_path / "nope.csv")) is None
    assert Config.resolve_variant(None) is None


def test_resolve_variant_follows_symlinks(tmp_path):
    """Cluster inputs are normally symlinked, not copied."""
    from phylogas.config import Config

    real = tmp_path / "real.csv"
    real.write_text("pid\n1\n")
    link_dir = tmp_path / "data"
    link_dir.mkdir()
    link = link_dir / "va_person.csv"
    link.symlink_to(real)

    assert Config.resolve_variant(str(link)) == link
    assert Config.resolve_variant(str(link) + ".xz") == link


# --------------------------------------------------------------------------
# painter output filtering
# --------------------------------------------------------------------------
def test_filter_key_detection(tmp_path):
    """--linelist-filter must recognise the identifier column in any of the
    file shapes it is handed: painter metadata, TwinSampler linelist, or a
    BeyondBaseline samples file."""
    from phylogas.painter.genetic_painter import _read_filter_keys

    # BeyondBaseline samples: alias_pid
    p = tmp_path / "surs_samples.csv"
    p.write_text("alias_pid,county\n123.45,Fairfax\n678.90,Arlington\n")
    ids, pids, col = _read_filter_keys(str(p))
    assert col == "alias_pid" and ids == {"123.45", "678.90"} and not pids

    # painter metadata: strain
    p = tmp_path / "meta.csv"
    p.write_text("strain,date\nUSA/VA-EHip-123.45/2021,2021-06-01\n")
    ids, pids, col = _read_filter_keys(str(p))
    assert col == "strain" and ids == {"123.45"}

    # bare pids
    p = tmp_path / "pids.csv"
    p.write_text("pid,x\n123,1\n456,2\n")
    ids, pids, col = _read_filter_keys(str(p))
    assert col == "pid" and pids == {"123", "456"} and not ids


def test_output_paths_match_compression():
    from phylogas.painter.genetic_painter import _output_paths

    assert _output_paths("p", "xz") == ("p.sequences.fasta.xz", "p.metadata.tsv.xz")
    assert _output_paths("p", "bgzf") == ("p.sequences.fasta.gz", "p.metadata.tsv.gz")
    assert _output_paths("p", None) == ("p.sequences.fasta", "p.metadata.tsv")


def test_rucc_file_present_and_parses():
    """The pinned USDA RUCC copy must be committed and loadable."""
    import pandas as pd

    rucc = Path(__file__).resolve().parents[1] / "data" / "Ruralurbancontinuumcodes2023.csv"
    assert rucc.is_file(), "data/Ruralurbancontinuumcodes2023.csv is missing"
    df = pd.read_csv(rucc, encoding="latin1")
    assert {"FIPS", "State", "County_Name", "Attribute", "Value"} <= set(df.columns)
    assert (df["Attribute"] == "RUCC_2023").any(), "no RUCC_2023 rows"


# --------------------------------------------------------------------------
# benchmark modules moved from BeyondBaseline / TwinSampler
# --------------------------------------------------------------------------
def test_benchmark_modules_import():
    import phylogas.benchmark.variants  # noqa: F401
    import phylogas.benchmark.truth_metrics  # noqa: F401
    import phylogas.benchmark.truth_runner  # noqa: F401


def test_variant_matchers_present_and_deterministic():
    """Both matchers moved from TwinSampler must be callable and stable."""
    import pandas as pd

    from phylogas.benchmark.variants import (
        mode1_temporal_match, mode2_bipartite_match, unroll_schedule,
        summarize_components, VARIANT_COLUMN,
    )

    comps = pd.DataFrame({
        "component_id": [1, 2, 3],
        "first_tick": [0, 20, 60],
        "component_size": [10, 5, 30],
    })
    sched = pd.DataFrame({
        "tick": [0, 20], "date": ["2021-06-01", "2021-06-21"],
        "variant": ["B.1.617.2", "AY.44"], "clusters": [1, 1], "sample_count": [10, 5],
    })
    real = unroll_schedule(sched)
    assert len(real) == 2

    for fn in (mode1_temporal_match, mode2_bipartite_match):
        a = fn(comps.copy(), real.copy())
        b = fn(comps.copy(), real.copy())
        assert a == b, f"{fn.__name__} is not deterministic"
        assert set(a.values()) <= {"B.1.617.2", "AY.44"}

    assert VARIANT_COLUMN == "variant_benchmark"


def test_summarize_components_rebuilds_from_events():
    """PhyloGAS must reconstruct component_summary without re-running
    component detection -- TwinSampler already did that."""
    import pandas as pd

    from phylogas.benchmark.variants import summarize_components

    events = pd.DataFrame({
        "component_id": [1, 1, 2],
        "exposure_tick": [5, 9, 30],
        "pid": ["a", "b", "c"],
    })
    s = summarize_components(events).set_index("component_id")
    assert s.loc[1, "first_tick"] == 5 and s.loc[1, "component_size"] == 2
    assert s.loc[2, "first_tick"] == 30 and s.loc[2, "component_size"] == 1


def test_variant_column_fallback():
    """Readers must accept the legacy column name."""
    import pandas as pd

    from phylogas.benchmark.truth_metrics import variant_column

    assert variant_column(pd.DataFrame(columns=["variant_benchmark"])) == "variant_benchmark"
    assert variant_column(pd.DataFrame(columns=["variant_label"])) == "variant_label"
    assert variant_column(pd.DataFrame(columns=["x"])) is None
