# ===========================================================================
# PhyloGAS pipeline
#
#   phylogas run --config config.yaml --cores all
#   phylogas run --config config.yaml --profile slurm
#
# or directly:
#   snakemake --configfile config.yaml --cores all
#
# Stage numbering matches the results/ directory layout in the README.
# ===========================================================================

import functools as _functools
import json as _json
import os
import sys
from pathlib import Path

try:                                    # the import path is stable, but do not
    from snakemake.exceptions import WorkflowError   # die at parse time over it
except ImportError:                     # pragma: no cover
    WorkflowError = ValueError

# The default config, only when none is given on the command line. As an
# unconditional directive it was loaded underneath every --configfile, so keys
# present only in config.yaml leaked into other runs.
if not config:
    configfile: "config.yaml"



# --------------------------------------------------------------------------
# Config access helpers
# --------------------------------------------------------------------------
def cfg(dotted, default=None):
    """Look up a dotted config key, with {placeholder} expansion."""
    node = config
    for part in dotted.split("."):
        if not isinstance(node, dict) or part not in node:
            return default
        node = node[part]
    return _expand(node) if isinstance(node, str) else node


def _expand(s, _depth=0):
    """Resolve {dotted.key} references the same way phylogas.config does."""
    import re
    if _depth > 10:
        return s
    def sub(m):
        node = config
        for part in m.group(1).split("."):
            if not isinstance(node, dict) or part not in node:
                return m.group(0)
            node = node[part]
        return str(node) if not isinstance(node, (dict, list)) else m.group(0)
    out = re.sub(r"\{([a-zA-Z0-9_.]+)\}", sub, s)
    return out if out == s else _expand(out, _depth + 1)


import os
import re


def _ascertainment_base():
    out = cfg("ascertainment.output", "")
    return re.sub(r"\.csv(\.gz|\.xz)?$", "", str(out)) if out else ""


def resolve_benchmark(key, kind):
    """'auto' | 'none' | <path> -> a path to depend on, or None.

    'auto' is the file this pipeline's line-list step writes, whether or not
    it exists yet. An explicit path is used if it exists (or is that same
    output). 'none' disables the benchmarks that need it.
    """
    raw = cfg(key, "auto")
    raw = "auto" if raw in (None, "") else str(raw)
    if raw.lower() in ("none", "skip", "off", "false"):
        return None
    produced = {"allevents": ALLEVENTS_OUT, "mugration": MUGRATION_OUT}[kind]
    if raw.lower() == "auto":
        # The pipeline's own line-list step writes both files, so 'auto' names
        # them as outputs-to-be rather than files that must already exist.
        # Snakemake then schedules each benchmark after ascertainment in the
        # same run. (Requiring existence made benchmarks vanish from a first
        # run and appear only on a second.)
        return produced
    cand = raw
    if os.path.normpath(cand) == os.path.normpath(produced):
        return produced
    alt = cand[:-3] if cand.endswith(".xz") else cand + ".xz"
    for p in (cand, alt):
        if os.path.exists(p):
            return p
    return None


RESULTS       = cfg("results_dir", "results")
# Derived settings -- location, start_tick, start_date -- filled by the same
# rule `phylogas` applies when it loads the config, so the painter and the
# line list cannot be handed different calendars or geography.
try:
    from phylogas.config import resolve_derived as _resolve_derived, ConfigError as _ConfigError
except ImportError:                     # pragma: no cover
    _resolve_derived = None
if _resolve_derived is not None:
    try:
        _resolve_derived(config, Path(workflow.configfiles[0]).parent
                         if workflow.configfiles else None, expand=_expand)
    except _ConfigError as _exc:
        raise WorkflowError(str(_exc))

PROJECT       = cfg("project_name", "phylogas")

THRESHOLD     = cfg("genetic_painter.entropy_thresholds")
PROBMATRIX    = cfg("genetic_painter.probability_matrix")
def _existing_variant(path):
    """The configured file, or its compressed / decompressed sibling on disk.

    fetch-data compresses the training alignment to .fasta.xz, while configs
    commonly name the plain .fasta; train reads either.
    """
    if not path or os.path.exists(path):
        return path
    for ext in (".xz", ".gz"):
        if os.path.exists(path + ext):
            return path + ext
        if path.endswith(ext) and os.path.exists(path[: -len(ext)]):
            return path[: -len(ext)]
    return path


ALIGN_FASTA   = _existing_variant(cfg("genetic_painter.align_fasta"))
PAINT_PREFIX  = cfg("genetic_painter.output_prefix", f"{RESULTS}/01_synthetic_genomes/{PROJECT}")

# genetic_painter.location is a mapping in YAML but the painter takes it as a
# JSON string, so it is read here once and shared with simulate_linelist.
_loc = cfg("genetic_painter.location", {}) or {}
if isinstance(_loc, str):
    import json as _json
    _loc = _json.loads(_loc)
LOCATION = _loc

_COMP         = cfg("genetic_painter.compression", "xz")
_EXT          = {"xz": ".xz", "bgzf": ".gz"}.get(_COMP, "")
PAINTED_FASTA = f"{PAINT_PREFIX}.sequences.fasta{_EXT}"
PAINTED_META  = f"{PAINT_PREFIX}.metadata.tsv{_EXT}"

LINELIST      = cfg("ascertainment.output", f"{RESULTS}/02_simulated_linelists/linelist.csv")
SAMPLE_DIR    = cfg("sampling.outdir", f"{RESULTS}/03_sampled_datasets")
BENCH_DIR     = cfg("benchmark.outdir", f"{RESULTS}/05_benchmarks")
ALGORITHMS    = cfg("sampling.algorithms", ["surs"])

# Resolved once so every rule and the printed summary agree.
# Two different questions, so two names.
#
# ALLEVENTS_OUT is where simulate_linelist WILL write the all-events table:
# derived from ascertainment.output, no existence check. Declared as that
# rule's output and depended on by the all_infections arm, so the arm works on
# a clean project.
#
# ALLEVENTS is whether a file is there to SCORE, which is what the optional
# benchmark rules ask -- `benchmark.allevents: none` opts out, and `auto` on a
# project that has not run yet simply skips them. Resolving the arm's
# dependency that way meant it silently vanished on a first run and only
# appeared on a second, because the file it needs is produced during the run
# that would have used it.
# Derived from LINELIST, not from the raw config key: TwinSampler builds this
# sibling from whatever --out it was given, and --out is LINELIST. Deriving it
# the same way means the two cannot disagree, and it inherits LINELIST's
# fallback instead of collapsing to "" when the key is absent.
ALLEVENTS_OUT = re.sub(r"\.csv(\.gz|\.xz)?$", "", str(LINELIST)) + "_allevents.csv.xz"
MUGRATION_OUT = re.sub(r"\.csv(\.gz|\.xz)?$", "", str(LINELIST)) + "_mugration.json"
ALLEVENTS = resolve_benchmark("benchmark.allevents", "allevents")
MUGRATION = resolve_benchmark("benchmark.truth_mugration", "mugration")


# --------------------------------------------------------------------------
# Targets
# --------------------------------------------------------------------------
# Clock mode -> which `phylogas benchmark clock` invocations to make.
# docs/clock_modes.md explains why these are not interchangeable. Defined here
# rather than inside the nextstrain block because `rule all` needs it too, and
# the two must agree: a target nothing produces is a DAG error.
_CLOCK_MODES = {
    "operational": ["operational"],
    "generative":  ["truth", "inferred"],
    "both":        ["all"],
}


def _clock_invocations():
    mode = str(cfg("nextstrain.clock.mode", "both"))
    if mode not in _CLOCK_MODES:
        raise ValueError(
            f"nextstrain.clock.mode is '{mode}'; expected one of "
            f"{', '.join(sorted(_CLOCK_MODES))}. See docs/clock_modes.md.")
    return _CLOCK_MODES[mode]


def _nextstrain_targets():
    """Auspice JSONs and the clock CSV, when nextstrain is enabled.

    A function because the arm list and the clock mode are both read from
    config inside the `if cfg("nextstrain.enabled")` block below, which has
    not run yet when `rule all` is defined.
    """
    if not cfg("nextstrain.enabled", False):
        return []
    out = cfg("nextstrain.outdir", f"{RESULTS}/04_nextstrain_builds")
    targets = [f"{out}/{a}/auspice.json" for a in _ns_arms()]
    # Only ask for the CSV when a rule will exist to produce it.
    if _clock_invocations():
        targets.append(f"{BENCH_DIR}/clock_estimates.csv")
    return targets


def _target_variant() -> str:
    """Which variant's states to keep.

    Defaults to the variant token in genetic_painter.painted_prefix ("E2" ->
    "2") so the painter and the line list cannot end up describing different
    epidemics, and so an existing config needs no new key. Empty for a
    single-variant model whose states carry no variant token.
    """
    explicit = str(cfg("ascertainment.target_variant", "") or "").strip()
    if explicit:
        return explicit
    painted = str(cfg("genetic_painter.painted_prefix", "") or "")
    m = re.search(r"(\d+)$", painted)
    return m.group(1) if m else ""



def _recipe_scenario(rid: str) -> str:
    """Scenario half of a recipe id."""
    return rid.split("__", 1)[0] if "__" in rid else rid


def _sample_recipes() -> list:
    """Every recipe BeyondBaseline writes a sample file for, this run.

    The registry decides, not a scenario x algorithm cross product: a
    target-blind algorithm contributes one recipe per stride rather than one
    per scenario, so the product would declare outputs that are never
    written. Restricted to sampling.algorithms, which is what the runner is
    told to run.
    """
    want = {str(a).strip().lower() for a in (ALGORITHMS or []) if str(a).strip()}
    try:
        from scenarios_simulation.recipes import all_recipes
    except ImportError:
        # No registry: fall back to the recipes the tree needs, which is the
        # minimum that keeps the DAG buildable with PhyloGAS alone.
        return [r for r in _ns_recipes() if _recipe_algo(r) in want]

    return [r for r in all_recipes() if _recipe_algo(r) in want]


def _sample_files() -> list:
    """Paths of those sample files."""
    return [f"{SAMPLE_DIR}/{r}_samples.csv.xz" for r in _sample_recipes()]


@_functools.lru_cache(maxsize=None)
def _ns_recipes() -> tuple:
    """Recipe ids carried into the phylodynamic stage, one tree each.

    "all" means every id BeyondBaseline's registry knows. An id that the
    registry has collapsed -- a target-blind algorithm names only the axes it
    reads -- is translated to its canonical form here, so an older config
    keeps working and the DAG still builds exactly one tree per sample.
    """
    want = cfg("sampling.nextstrain_recipes", []) or []
    if isinstance(want, str):
        want = [want]
    want = [str(w).strip() for w in want if str(w).strip()]
    if not want:
        raise WorkflowError(
            'sampling.nextstrain_recipes is empty. Set it to one or more '
            'recipe ids from `scenarios-recipes`, e.g. ["4S__surs"], or '
            '"all".')

    try:
        from scenarios_simulation.recipes import all_recipes, aliases
    except ImportError:
        # PhyloGAS-only checkout: take the config at its word. "all" cannot be
        # expanded without the registry, so say which import is missing.
        if any(w.lower() == "all" for w in want):
            raise WorkflowError(
                'sampling.nextstrain_recipes is "all", which needs '
                "BeyondBaseline's recipe registry to expand. Install "
                "BeyondBaseline, or list the recipe ids explicitly.")
        return tuple(dict.fromkeys(want))

    known, alias = all_recipes(), aliases()
    if any(w.lower() == "all" for w in want):
        if len(want) > 1:
            raise WorkflowError(
                'sampling.nextstrain_recipes mixes "all" with specific ids; '
                "use one or the other.")
        # Only the algorithms this run samples: a recipe for an algorithm
        # absent from sampling.algorithms has no sample file to build from.
        algos = {str(a).strip().lower() for a in (ALGORITHMS or [])}
        every = [r for r in known if _recipe_algo(r) in algos]
        if not every:
            raise WorkflowError(
                f'sampling.nextstrain_recipes is "all", but no recipe matches '
                f"sampling.algorithms ({sorted(algos)}). Known algorithms: "
                f"{sorted({_recipe_algo(r) for r in known})}.")
        return tuple(every)

    out = []
    for w in want:
        rid = w if w in known else alias.get(w)
        if rid is None:
            raise WorkflowError(
                f"sampling.nextstrain_recipes has {w!r}, which is not a "
                f"recipe id. `scenarios-recipes` lists all {len(known)}.")
        if rid != w:
            # Collapsed: several grid cells share one sample, so two entries
            # could name the same tree. Deduplicated below.
            print(f"note: nextstrain_recipes {w!r} -> {rid!r} "
                  f"(that algorithm ignores the target distribution)",
                  file=sys.stderr)
        out.append(rid)
    return tuple(dict.fromkeys(out))


def _recipe_algo(rid: str) -> str:
    """Algorithm half of a recipe id, which is what BeyondBaseline samples by."""
    return rid.split("__", 1)[1] if "__" in rid else rid


def _ns_all_infections() -> bool:
    """Whether the control arm is requested.

    nextstrain.builds no longer carries an algorithm list -- which recipes get
    a tree is sampling.nextstrain_recipes -- so the only thing left to read
    from it is whether the all_infections arm is on.
    """
    for b in (cfg("nextstrain.builds", []) or []):
        if isinstance(b, dict) and b.get("type") == "all_infections":
            return True
    return False


def _ns_arms() -> list:
    """Every arm that gets a tree: one per recipe, plus the control if on.

    An arm name is a recipe id, which is also the sample filename stem and the
    ncov build suffix, so one string identifies the tree end to end.
    """
    if not cfg("nextstrain.enabled", False):
        return []
    return list(_ns_recipes()) + (["all_infections"]
                                  if _ns_all_infections() else [])


def _pkg_version(dist: str) -> str:
    """Installed version of a sibling package, for use as a rule `params`.

    Snakemake invalidates an output when a rule's params change, but it tracks
    files -- not the code that produced them. TwinSampler and BeyondBaseline
    are installed separately, so reinstalling one left its outputs looking
    current: a line list built with a swapped age map, or samples keyed on an
    all-NaN denominator, both survived a rerun untouched and had to be deleted
    by hand. Both now derive their version from git, so every commit changes
    this string and the affected stage rebuilds on its own.
    """
    from importlib.metadata import PackageNotFoundError, version

    try:
        return version(dist)
    except PackageNotFoundError:
        # The rule's own command will fail with something clearer.
        return "absent"


def _rucc_override(wildcards) -> dict:
    """population.rucc_file as a named input, or nothing.

    TwinSampler ships the USDA 2023 table and uses it when --rucc is omitted,
    so this key is only for running against a different edition.
    """
    value = str(cfg("population.rucc_file", "") or "").strip()
    return {"rucc": value} if value else {}


def _required_path(dotted: str) -> str:
    """A config path that must be set, for use in an `input:` block.

    An empty value reaching `input:` becomes a nameless missing file, so fail
    at parse time naming the key instead.
    """
    value = str(cfg(dotted, "") or "").strip()
    if not value:
        raise WorkflowError(f"{dotted} is not set, and simulate_linelist requires it.")
    return value


rule all:
    """Default target: painted genomes plus every requested sampled subset."""
    input:
        PAINTED_FASTA,
        PAINTED_META,
        # One subset per recipe that gets a tree. The subset exists to feed
        # the tree, so when nextstrain is off there is nothing to subset.
        expand(f"{SAMPLE_DIR}/{PROJECT}.{{recipe}}.fasta{_EXT}",
               recipe=[a for a in _ns_arms() if a != "all_infections"]),
        _nextstrain_targets(),


rule paint_only:
    """Stop after the genetic painter (stages 1-2)."""
    input:
        PAINTED_FASTA,
        PAINTED_META,


# --------------------------------------------------------------------------
# Stage 0: seed acquisition
# --------------------------------------------------------------------------
rule prep_seeds:
    """Fetch seed sequences for the variant of interest and build a schedule."""
    output:
        directory(cfg("seeds.output_folder", "data/importations/sequences")),
    params:
        state=cfg("population.state_name", "Virginia"),
        pango=cfg("variant.pango", "B.1.617.2"),
        outlier=cfg("seeds.outlier_method", "none"),
        config=workflow.configfiles[0] if workflow.configfiles else "config.yaml",
        # Absolute ticks come from the ABM's own config, so the seed schedule
        # and the EpiHiper output share one calendar.
        abm_arg=(f" --abm-config {cfg('abm.config')}" if cfg("abm.config", "") else ""),
    shell:
        "phylogas prep-seeds --config {params.config} "
        "--state {params.state} --pango {params.pango} "
        "--outlier-method {params.outlier} --seed-mode "
        "--output-folder {output}{params.abm_arg}"


# --------------------------------------------------------------------------
# Stage 1: entropy training ("the Map")
# --------------------------------------------------------------------------
rule train_entropy:
    """Derive per-site entropy thresholds and the substitution matrix from an MSA.

    Deterministic given the same MSA, so Snakemake will skip it on reruns.
    """
    input:
        msa=ALIGN_FASTA,
    output:
        thresholds=THRESHOLD,
        probmatrix=PROBMATRIX,
    params:
        config=workflow.configfiles[0] if workflow.configfiles else "config.yaml",
    shell:
        "phylogas train --config {params.config} --align-fasta {input.msa}"


# --------------------------------------------------------------------------
# Stage 2: the Genetic Painter
# --------------------------------------------------------------------------
rule paint_network:
    """Overlay evolving viral genomes onto the EpiHiper transmission network.

    Single-threaded in the painting loop; `threads` is sized for the external
    multi-threaded xz writer, which is otherwise the throughput bottleneck.
    """
    input:
        thresholds=THRESHOLD,
        probmatrix=PROBMATRIX,
        graph=cfg("epihiper.output_csv"),
        seeds=cfg("genetic_painter.seed_fasta"),
    output:
        fasta=PAINTED_FASTA,
        metadata=PAINTED_META,
    threads: cfg("genetic_painter.compression_threads", 8) or 1
    resources:
        mem_mb=cfg("resources.paint_mem_mb", 64000),
        runtime=cfg("resources.paint_runtime_min", 240),
    params:
        config=workflow.configfiles[0] if workflow.configfiles else "config.yaml",
    shell:
        "phylogas paint --config {params.config} "
        "--analysis-type generate_sequence "
        "--compression-threads {threads}"


# --------------------------------------------------------------------------
# Stage 3: ascertainment  (TwinSampler)
# --------------------------------------------------------------------------
rule simulate_linelist:
    """Apply the 'lens of ascertainment' to produce a realistic, biased linelist.

    Provided by the TwinSampler package (installed via environment.yml).
    """
    input:
        # Optional override of the RUCC table TwinSampler bundles. Present
        # as an input only when configured, so a refreshed file still
        # triggers a rerun.
        unpack(_rucc_override),
        # Required by simulate_linelist, so inputs rather than params:
        # Snakemake checks them when it builds the DAG.
        graph=cfg("epihiper.output_csv"),
        people=cfg("population.demographics_file", cfg("population.persontrait_file", "")),
        households=_required_path("population.household_file"),
        ascertain=_required_path("ascertainment.parameters"),
    output:
        linelist=LINELIST,
        # Written because the shell always passes --output_all_events.
        # Declared so Snakemake tracks and cleans it, and so the
        # all_infections arm can depend on it rather than on its existence.
        allevents=ALLEVENTS_OUT,
        # Also always written (TwinSampler derives it from --out). Declared so
        # benchmark_mugration can depend on it within the same run.
        mugration=MUGRATION_OUT,
    params:
        start_date=cfg("genetic_painter.start_date"),
        start_tick=cfg("genetic_painter.start_tick"),
        stop_tick=(cfg("genetic_painter.start_tick", 0) + cfg("genetic_painter.num_ticks", 0)),
        seed=cfg("random_seed", 42),
        schedule=cfg("ascertainment.schedule_input", ""),
        # TwinSampler changed its own default to just_components when
        # benchmark variant assignment moved to `phylogas assign-variants`.
        # This fallback used to say variant_bipartite, which would have
        # silently re-enabled the moved step for a config without the key.
        variant_mode=cfg("ascertainment.variant_mode", "just_components"),
        # Geography must match what the painter puts in its FASTA headers:
        # both build the Nextstrain `strain` id from these, and ncov joins
        # sequences to metadata on that string. Single source of truth.
        country=LOCATION.get("country", "USA"),
        region=LOCATION.get("region", "North America"),
        division=LOCATION.get("division", "Virginia"),
        division_abbr=LOCATION.get("divisionAbbr", "VA"),
        # One variant only, applied before the time filter.
        target_variant_arg=(f" --target_variant {_target_variant()}"
                            if _target_variant() else ""),
        prefix_override=_json.dumps(
            cfg("ascertainment.prefix_override",
                ["A", "P", "I", "dM", "hM"])),
        # Not passed to the command; present so a TwinSampler reinstall
        # invalidates this line list.
        twin_sampler_version=_pkg_version("twin-sampler"),
        rucc_arg=lambda wc, input: (f" --rucc {input.rucc}"
                                    if "rucc" in input.keys() else ""),
    shell:
        "simulate_linelist --epihiper {input.graph} --people {input.people} "
        "--households {input.households}{params.rucc_arg} "
        "--ascertain {input.ascertain} --start_date {params.start_date} "
        "--start_tick {params.start_tick} --stop_tick {params.stop_tick} "
        "--out {output.linelist} --seed {params.seed} --output_all_events "
        "--schedule_input {params.schedule} --variant_mode {params.variant_mode} "
        "--country {params.country:q} --region {params.region:q} "
        "--division {params.division:q} --division_abbr {params.division_abbr:q} "
        "--prefix_override {params.prefix_override:q}"
        "{params.target_variant_arg}"


# --------------------------------------------------------------------------
# Stage 4: adaptive sampling  (BeyondBaseline)
# --------------------------------------------------------------------------
rule sample_scenarios:
    """Select which cases get sequenced, under a fixed budget.

    Provided by the BeyondBaseline package (installed via environment.yml).
    """
    input:
        linelist=LINELIST,
        population=cfg("population.demographics_file", cfg("population.persontrait_file", "")),
    output:
        # Every sample file the run produces, from BeyondBaseline's registry.
        # One rule rather than one per algorithm: the runner loops scenarios
        # x algorithms internally, so invoking it once per algorithm would
        # re-walk the whole scenario sweep each time. Registry-derived rather
        # than a cross product because a target-blind algorithm writes one
        # file per stride, not one per scenario.
        samples=_sample_files(),
    params:
        outdir=SAMPLE_DIR,
        batch=cfg("sampling.batch_size", 400),
        seed=cfg("random_seed", 42),
        norepl="--no-replacement" if cfg("sampling.no_replacement", True) else "",
        algorithms=" ".join(str(a) for a in ALGORITHMS),
        # Not passed to the command; present so a BeyondBaseline reinstall
        # invalidates these samples.
        beyond_baseline_version=_pkg_version("beyond-baseline"),
    shell:
        "scenarios-runner --linelist {input.linelist} "
        "--population {input.population} --outdir {params.outdir} "
        "--batch-size {params.batch} --seed {params.seed} {params.norepl} "
        "--save-samples --algorithms {params.algorithms}"


# --------------------------------------------------------------------------
# Stage 5: FASTA subsetting
# --------------------------------------------------------------------------
rule subset_fasta:
    """Pull the selected strains out of the full ground-truth FASTA."""
    input:
        fasta=PAINTED_FASTA,
        samples=f"{SAMPLE_DIR}/{{recipe}}_samples.csv.xz",
    output:
        # Named by recipe, not algorithm: two recipes can share an algorithm
        # (1S__surs and 4S__surs), and an algorithm-named output would have
        # them overwrite each other.
        fasta=f"{SAMPLE_DIR}/{PROJECT}.{{recipe}}.fasta{_EXT}",
    shell:
        "phylogas subset-fasta -m {input.samples} -f {input.fasta} -o {output.fasta}"


# --------------------------------------------------------------------------
# Stage 6: phylodynamic reconstruction  (Nextstrain)
# --------------------------------------------------------------------------
if cfg("nextstrain.enabled", False):

    NS_DIR    = cfg("nextstrain.dir", "")
    NS_STAGE  = cfg("nextstrain.stage_subdir", "data/phylogas")
    NS_RUNNER = cfg("nextstrain.runner", "nextstrain")
    # Only meaningful for runner: snakemake. Snakemake builds ncov's own
    # workflow/envs/nextstrain.yaml per working directory unless told to share
    # one prefix, and that env is a multi-GB solve we do not want repeated per
    # arm. Empty means "let snakemake decide" (.snakemake/conda under the ncov
    # checkout).
    NS_CONDA_PREFIX = cfg("nextstrain.conda_prefix", "")
    NS_OUT    = cfg("nextstrain.outdir", f"{RESULTS}/04_nextstrain_builds")

    # nextstrain.builds decides which arms exist. A `strategy` entry expands to
    # one arm per algorithm (defaulting to sampling.algorithms); an
    # `all_infections` entry is a single arm with its own inputs and its own
    # subsampling scheme.
    NS_STRATEGY_ARMS = _ns_recipes()
    NS_ALL_INFECTIONS = _ns_all_infections()
    NS_ARMS = _ns_arms()

    NS_CLOCK = _clock_invocations()

    rule nextstrain_config_strategy:
        """Render, stage and validate the ncov config for one sampling arm.

        Validation is the point of making this its own rule: it refuses a
        build whose sequences do not join to its metadata, or whose traits
        cannot be scored, before any hours are spent.
        """
        wildcard_constraints:
            recipe=("|".join(re.escape(r) for r in NS_STRATEGY_ARMS)
                    if NS_STRATEGY_ARMS else "$^")
        input:
            fasta=f"{SAMPLE_DIR}/{PROJECT}.{{recipe}}.fasta{_EXT}",
            metadata=f"{SAMPLE_DIR}/{{recipe}}_samples.csv.xz",
        output:
            configfile=f"{NS_DIR}/{NS_STAGE}/{{recipe}}/config.yaml",
        params:
            config=workflow.configfiles[0] if workflow.configfiles else "config.yaml",
        shell:
            "phylogas nextstrain-config --config {params.config} "
            "--build-type strategy --recipe {wildcards.recipe}"

    if NS_ALL_INFECTIONS:

        rule nextstrain_config_all_infections:
            """The specialised arm: every painted infection plus the all-events
            line list, capped by the country_17k scheme. Separate from the
            strategy rule because its inputs differ, not just its subsampling.
            """
            input:
                fasta=PAINTED_FASTA,
                metadata=ALLEVENTS_OUT,
            output:
                configfile=f"{NS_DIR}/{NS_STAGE}/all_infections/config.yaml",
            params:
                config=workflow.configfiles[0] if workflow.configfiles else "config.yaml",
            shell:
                "phylogas nextstrain-config --config {params.config} "
                "--build-type all_infections"

    rule nextstrain_build:
        """Run the ncov workflow through the Nextstrain CLI.

        `nextstrain build` resolves the runtime (docker / conda / singularity /
        ambient) from `nextstrain check-setup --set-default`, so the same rule
        works on a laptop and on a cluster. It also mounts only the build
        directory, which is why nextstrain_config stages inputs inside the
        checkout and writes relative paths.

        Set nextstrain.runner: snakemake to bypass the CLI and call snakemake
        directly, which needs the full ncov environment already active.
        """
        wildcard_constraints:
            arm=("|".join(re.escape(a) for a in NS_ARMS)
                 if NS_ARMS else "$^")
        input:
            configfile=f"{NS_DIR}/{NS_STAGE}/{{arm}}/config.yaml",
        output:
            auspice=f"{NS_OUT}/{{arm}}/auspice.json",
        params:
            ns_dir=NS_DIR,
            rel_config=f"{NS_STAGE}/{{arm}}/config.yaml",
            runner=NS_RUNNER,
            conda_prefix=(f"--conda-prefix {NS_CONDA_PREFIX}"
                          if NS_CONDA_PREFIX else ""),
            project=PROJECT,
            json_prefix=cfg("nextstrain.auspice_json_prefix", "ncov"),
        threads: 8
        shell:
            r"""
            if [ "{params.runner}" = "nextstrain" ]; then
                nextstrain build {params.ns_dir} \
                    --configfile {params.rel_config} \
                    --cores {threads} --rerun-incomplete
            else
                snakemake --snakefile {params.ns_dir}/Snakefile \
                    --directory {params.ns_dir} \
                    --configfile {params.ns_dir}/{params.rel_config} \
                    --use-conda {params.conda_prefix} \
                    --cores {threads} --rerun-incomplete
            fi
            mkdir -p $(dirname {output.auspice})
            # Spelled out rather than globbed on the arm suffix: a glob also
            # matches the _tip-frequencies and _root-sequence siblings, which
            # would cp three files onto one destination path.
            cp {params.ns_dir}/auspice/{params.json_prefix}_{params.project}_{wildcards.arm}.json \
               {output.auspice}
            """

    if NS_CLOCK:

        rule clock_estimates:
            """Measure, infer and record the molecular clock rate per arm.

            One rule rather than one per arm, because the CLI accumulates into
            a single CSV -- replacing any existing row for the same
            (arm, quantity) -- so concatenating ragged per-arm files is
            unnecessary.

            The builds are inputs only when the mode needs them: `truth` is
            tree-free and reads the painted sequences directly, while
            `inferred` and `operational` need a finished ncov build.
            """
            input:
                fasta=PAINTED_FASTA,
                metadata=PAINTED_META,
                auspice=([f"{NS_OUT}/{a}/auspice.json" for a in NS_ARMS]
                         if any(m in ("all", "inferred", "operational")
                                for m in NS_CLOCK) else []),
            output:
                csv=f"{BENCH_DIR}/clock_estimates.csv",
            params:
                config=workflow.configfiles[0] if workflow.configfiles else "config.yaml",
                arms=" ".join(NS_ARMS),
                modes=" ".join(NS_CLOCK),
                ns_dir=NS_DIR,
                project=PROJECT,
            shell:
                r"""
                for arm in {params.arms}; do
                  for mode in {params.modes}; do
                    phylogas benchmark clock --config {params.config}                         --arm "$arm" --mode "$mode"                         --build-dir {params.ns_dir}/results/{params.project}_"$arm"                         --out {output.csv} || true
                  done
                done
                # The loop tolerates per-arm failures -- a missing augur, say --
                # but the target must exist for the DAG, and an empty file is
                # more honest than a stale one.
                test -f {output.csv} || : > {output.csv}
                """


# --------------------------------------------------------------------------
# Stage 7: benchmarking
# --------------------------------------------------------------------------
# Compares the inferred tree against the painter's ground truth
# (topological F1, mugration cosine similarity).
#
if ALLEVENTS:

 rule assign_variants:
     """Attach benchmark variant labels to transmission components.

     Moved here from TwinSampler: these labels are ground truth for prevalence
     estimation, matched against a real importation schedule, and deliberately
     independent of the lineage a genome implies.
     """
     input:
         allevents=ALLEVENTS or "",
         # benchmark.variant_schedule: the prep-seeds sublineage schedule by
         # default, or an overlay. ascertainment.schedule_input is the older
         # key and is still honoured when the new one is absent.
         schedule=(cfg("benchmark.variant_schedule", "")
                   or cfg("ascertainment.schedule_input", "")),
     output:
         csv=f"{BENCH_DIR}/allevents_variants.csv.xz",
     params:
         mode=cfg("benchmark.variant_mode", "bipartite"),
     shell:
         "phylogas assign-variants --allevents {input.allevents} "
         "--schedule {input.schedule} --mode {params.mode} --out {output.csv}"


 rule benchmark_truth:
      """Score sampled sets against ABM ground truth."""
      input:
          samples=_sample_files(),
          infections=f"{BENCH_DIR}/allevents_variants.csv.xz",
      output:
          csv=f"{BENCH_DIR}/AUC_truth_rankings.csv",
      params:
          glob=f"{SAMPLE_DIR}/*_samples.csv.xz",
      shell:
          "phylogas benchmark truth --samples '{params.glob}' "
          "--infections {input.infections} --out {output.csv}"


if MUGRATION and ALLEVENTS:

    rule benchmark_mugration:
        """Score every sampled strategy against the ABM mugration truth.

        The transmission graph is built once and reused across all strategies,
        so this costs about the same as the old inline version in
        BeyondBaseline's run_all_scenarios.py.
        """
        input:
            truth=MUGRATION or "",
            infections=ALLEVENTS or "",
            samples=_sample_files(),
        output:
            csv=f"{BENCH_DIR}/Mugration_Metrics.csv",
        params:
            glob=f"{SAMPLE_DIR}/*_samples.csv.xz",
        shell:
            "phylogas benchmark mugration --truth {input.truth} "
            "--samples '{params.glob}' --infections {input.infections} "
            "--out {output.csv}"


if ALLEVENTS:

 rule benchmark_sequence:
     """Parent->child divergence in the painted genomes (sanity + rate check)."""
     input:
         painted=PAINTED_FASTA,
         infections=ALLEVENTS or "",
     output:
         csv=f"{BENCH_DIR}/sequence_divergence.csv",
     shell:
         "phylogas benchmark sequence --painted {input.painted} "
         "--infections {input.infections} --out {output.csv}"


# --------------------------------------------------------------------------
# End-to-end targets
#
# `all` (the default) stops at painted genomes, sampled subsets and -- with
# nextstrain enabled -- the trees and clock_estimates.csv. The ground-truth
# benchmarks depend on the line-list step's all-events and mugration outputs,
# which Snakemake schedules first, so one invocation runs everything:
#
#   phylogas run --config config.yaml --cores all benchmarks
#   phylogas run --config config.yaml --cores all full       # all + benchmarks
# --------------------------------------------------------------------------
def _benchmark_targets():
    targets = []
    if ALLEVENTS:
        targets += [f"{BENCH_DIR}/AUC_truth_rankings.csv",
                    f"{BENCH_DIR}/sequence_divergence.csv"]
    if MUGRATION and ALLEVENTS:
        targets.append(f"{BENCH_DIR}/Mugration_Metrics.csv")
    return targets


rule benchmarks:
    """Every ground-truth benchmark the configured inputs allow."""
    input:
        _benchmark_targets(),


rule full:
    """The default target plus the ground-truth benchmarks."""
    input:
        rules.all.input,
        _benchmark_targets(),
