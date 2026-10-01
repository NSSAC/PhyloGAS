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
