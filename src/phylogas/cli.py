"""Unified `phylogas` command line interface.

Wraps the previously-disparate scripts (`seq_prep.py`, `genetic_painter.py`,
`subset_fasta_streaming.py`, ...) behind one verb-based entry point, so the
suite behaves like a cohesive tool rather than a set of academic scripts.

Every subcommand accepts ``--config config.yaml``. Where a stage previously
took a long list of flags, those flags remain available and override the
config, so existing invocations can be migrated incrementally. See
``command_map.md`` for the old-to-new mapping.

Implemented with argparse rather than Typer/Click deliberately: it is in the
standard library, so `phylogas --help` works in a bare Python environment
before any third-party dependency is resolved.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

from . import __version__

_PKG = Path(__file__).resolve().parent


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def _load_config(args) -> "object":
    from .config import Config, ConfigError

    try:
        return Config.load(args.config)
    except ConfigError as exc:
        sys.exit(f"ERROR: {exc}")


def _run(cmd: list[str], dry_run: bool = False) -> int:
    """Echo and execute a subprocess, returning its exit status."""
    printable = " ".join(str(c) for c in cmd)
    print(f"+ {printable}", flush=True)
    if dry_run:
        return 0
    return subprocess.call([str(c) for c in cmd])


def _cfg_or_flag(cfg, dotted: str, flag_value, default=None):
    """Command-line flag wins; otherwise fall back to config, then default."""
    if flag_value is not None:
        return flag_value
    if cfg is None:
        return default
    return cfg.get(dotted, default=default)


# --------------------------------------------------------------------------
# phylogas paint
# --------------------------------------------------------------------------
def cmd_paint(args) -> int:
    """Stage 2: overlay viral genomes onto the transmission network."""
    cfg = _load_config(args) if args.config else None
    painter = _PKG / "painter" / "genetic_painter.py"

    def pick(dotted, val, default=None):
        return _cfg_or_flag(cfg, dotted, val, default)

    gp = "genetic_painter."
    cmd = [
        sys.executable, "-u", str(painter),
        "--analysis_type", args.analysis_type,
        "--threshold_file", pick(gp + "entropy_thresholds", args.threshold_file),
        "--base_threshold_df", pick(gp + "probability_matrix", args.base_threshold_df),
        "--random_number_seed", str(pick("random_seed", args.seed, 43)),
    ]

    optional = [
        ("--align_fasta", pick(gp + "align_fasta", args.align_fasta)),
        ("--seed_fasta", pick(gp + "seed_fasta", args.seed_fasta)),
        ("--input_graph_csv", pick("epihiper.output_csv", args.input_graph_csv)),
        ("--output_prefix", pick(gp + "output_prefix", args.output_prefix)),
        ("--start_date", pick(gp + "start_date", args.start_date)),
        ("--start_tick", pick(gp + "start_tick", args.start_tick)),
        ("--num_ticks", pick(gp + "num_ticks", args.num_ticks)),
        ("--reference", pick(gp + "reference_fasta", args.reference)),
        ("--compression", pick(gp + "compression", args.compression)),
        ("--compression_level", pick(gp + "compression_level", args.compression_level)),
        ("--compression_threads", pick(gp + "compression_threads", args.compression_threads)),
        ("--persontrait_file", pick("population.demographics_file", args.persontrait_file)
                               or pick("population.persontrait_file", None)),
        ("--add_metadata", pick(gp + "add_metadata", args.add_metadata)),
        ("--input_graph_painted_prefix", pick(gp + "painted_prefix", args.painted_prefix)),
        ("--initial_viral_load", pick(gp + "initial_viral_load", args.initial_viral_load)),
    ]
    for flag, value in optional:
        if value is not None:
            cmd += [flag, str(value)]

    for spec in (args.linelist_filter or []):
        cmd += ["--linelist_filter", str(spec)]
    for spec in (_cfg_or_flag(cfg, gp + "linelist_filters", None, []) or []):
        if spec not in (args.linelist_filter or []):
            cmd += ["--linelist_filter", str(spec)]

    model = pick(gp + "mutation_model", args.mutation_model, "rate_limit")
    if model in ("rate_limit", "rate_limited"):
        cmd.append("--rate_limit")
    elif model == "poor":
        cmd.append("--poor")

    if pick(gp + "proportional", args.proportional, True):
        cmd.append("--proportional")
    else:
        cmd.append("--neutral")

    return _run(cmd, args.dry_run)


# --------------------------------------------------------------------------
# phylogas train
# --------------------------------------------------------------------------
def cmd_train(args) -> int:
    """Stage 1: derive the entropy 'fitness map' from a real-world MSA."""
    args.analysis_type = "entropy_analysis"
    for attr in (
        "input_graph_csv", "output_prefix", "start_date", "start_tick", "num_ticks",
        "reference", "compression", "compression_level", "compression_threads",
        "persontrait_file", "add_metadata", "painted_prefix", "initial_viral_load",
        "seed_fasta", "linelist_filter",
    ):
        setattr(args, attr, None)
    args.mutation_model = None
    return cmd_paint(args)


# --------------------------------------------------------------------------
# phylogas prep-seeds
# --------------------------------------------------------------------------
def cmd_prep_seeds(args) -> int:
    """Stage 0: fetch seed sequences and build an importation schedule."""
    cfg = _load_config(args) if args.config else None
    script = _PKG / "seqprep" / "seq_prep.py"
    cmd = [sys.executable, "-u", str(script)]
    pairs = [
        ("--state", _cfg_or_flag(cfg, "population.state_name", args.state)),
        ("--pango", _cfg_or_flag(cfg, "variant.pango", args.pango)),
        ("--output_folder", _cfg_or_flag(cfg, "seeds.output_folder", args.output_folder)),
        ("--outlier_method", _cfg_or_flag(cfg, "seeds.outlier_method", args.outlier_method)),
    ]
    for flag, value in pairs:
        if value is not None:
            cmd += [flag, str(value)]
    if args.seed_mode:
        cmd.append("--seed_mode")
    if _cfg_or_flag(cfg, "seeds.insecure_download", args.insecure_download, False):
        cmd.append("--insecure_download")

    # Sublineage scope. Default ON (query "B.1.617.2*" and roll descendants up).
    include_sub = _cfg_or_flag(cfg, "variant.include_sublineages",
                               args.include_sublineages, True)
    cmd.append("--include_sublineages" if include_sub else "--no_include_sublineages")
    return _run(cmd, args.dry_run)


# --------------------------------------------------------------------------
# phylogas subset-fasta
# --------------------------------------------------------------------------
def cmd_subset_fasta(args) -> int:
    """Stage 5: pull the sampled strains out of the full ground-truth FASTA."""
    script = _PKG / "seqprep" / "subset_fasta_streaming.py"
    cmd = [
        sys.executable, "-u", str(script),
        "-m", args.metadata,
        "-f", args.fasta,
        "-o", args.output,
    ]
    return _run(cmd, args.dry_run)


# --------------------------------------------------------------------------
# phylogas run  (Snakemake driver)
# --------------------------------------------------------------------------
def cmd_run(args) -> int:
    """Run the whole pipeline via Snakemake."""
    if shutil.which("snakemake") is None:
        sys.exit(
            "ERROR: snakemake is not on PATH.\n"
            "  mamba env create -f environment.yml && conda activate phylogas_env"
        )
    snakefile = Path(args.snakefile) if args.snakefile else Path("Snakefile")
    if not snakefile.is_file():
        sys.exit(f"ERROR: Snakefile not found: {snakefile}")

    cmd = ["snakemake", "--snakefile", str(snakefile), "--configfile", str(args.config)]
    cmd += ["--cores", str(args.cores)]
    if args.profile:
        cmd += ["--profile", args.profile]
    if args.dry_run:
        cmd.append("--dry-run")
    if args.until:
        cmd += ["--until", args.until]
    cmd += args.extra
    return _run(cmd, dry_run=False)


# --------------------------------------------------------------------------
# phylogas validate-config
# --------------------------------------------------------------------------
def cmd_validate_config(args) -> int:
    """Load the config, expand placeholders, and report which inputs exist."""
    cfg = _load_config(args)
    print(f"Config OK: {cfg}")

    checks = [
        ("genetic_painter.entropy_thresholds", "entropy thresholds"),
        ("genetic_painter.probability_matrix", "probability matrix"),
        ("genetic_painter.reference_fasta", "reference FASTA"),
        ("genetic_painter.seed_fasta", "seed FASTA"),
        ("epihiper.output_csv", "EpiHiper transmission log"),
        ("population.demographics_file", "demographics (derived)"),
        ("population.fips_file", "county FIPS lookup"),
    ]
    missing = 0
    print("\nInput files:")
    for key, label in checks:
        raw = cfg.get(key, default=None)
        if raw is None:
            print(f"  [ -- ] {label:32s} (not configured)")
            continue
        found = cfg.resolve_variant(raw)
        if found is None:
            print(f"  [MISS] {label:32s} {raw}")
            missing += 1
        elif str(found) != str(Path(str(raw)).expanduser()):
            # Same file, different compression extension than configured.
            print(f"  [okay] {label:32s} {found}")
            print(f"         (config says {raw})")
        else:
            print(f"  [okay] {label:32s} {raw}")

    if missing:
        print(f"\n{missing} configured input(s) are missing.")
        print("Run `phylogas fetch-data --with-simulations` to download them from Zenodo.")
        return 1
    print("\nAll configured inputs are present.")
    return 0


# --------------------------------------------------------------------------
# phylogas fetch-data
# --------------------------------------------------------------------------
def cmd_build_demographics(args) -> int:
    """Join the population files into the painter's demographics table.

    The file the config calls `demographics_file` is NOT the EpiHiper
    persontrait file: it is this derived join. See the module docstring of
    phylogas.popprep.merge_persontrait for why.
    """
    from .popprep.merge_persontrait import build_demographics

    build_demographics(
        persontrait=args.persontrait,
        person=args.person,
        fips=args.fips,
        out=args.out,
        household=args.household,
        residence=args.residence,
    )
    return 0


_STATE_NAMES = {"va": "Virginia", "ca": "California", "ga": "Georgia",
                "ma": "Massachusetts", "mn": "Minnesota", "wa": "Washington"}


def cmd_fetch_data(args) -> int:
    """Download the synthetic population files from the UVA Dataverse."""
    from . import sources as dv

    cfg = _load_config(args) if args.config else None
    dest_root = Path(_cfg_or_flag(cfg, "data_dir", args.dest, "data"))
    states = args.states or ([cfg.get("population.state", default="va")] if cfg else ["va"])
    states = [s.lower() for s in states]

    unknown = [s for s in states if s not in dv.CORE_FILES]
    if unknown:
        sys.exit(
            f"ERROR: no Dataverse deposit registered for: {', '.join(unknown)}\n"
            f"       Known states: {', '.join(sorted(dv.CORE_FILES))}"
        )

    plan = dv.file_plan(states, with_epihiper_inputs=args.with_epihiper_inputs)
    total_mb = sum(mb for *_, mb, _ in plan)

    print(f"Destination : {dest_root.resolve()}")
    print(f"States      : {', '.join(states)}")
    for st in states:
        rec = dv.population_record(st)
        print(f"  {st}: https://doi.org/{rec[1]}")
    print(f"\nFiles to fetch ({len(plan)}, {total_mb / 1024:.2f} GB compressed):")
    for st, name, url, mb, md5 in plan:
        print(f"  {mb:9.1f} MB  {name}")
    if not args.with_epihiper_inputs:
        skipped = sum(
            mb for st in states for mb, _ in dv.EPIHIPER_INPUT_FILES.get(st, {}).values()
        )
        print(f"\n  (skipping {skipped / 1024:.2f} GB of EpiHiper contact networks;")
        print( "   nothing in PhyloGAS reads them. Use --with-epihiper-inputs if you")
        print( "   intend to run the ABM yourself.)")

    sim_plan = dv.simulation_plan(states) if args.with_simulations else []
    if args.with_simulations:
        missing_sim = [s for s in states if s not in dv.SIMULATION_FILES]
        print(f"\nEpiHiper simulation replicates (Zenodo {dv.ZENODO_DOI}):")
        for _st, key, _url, mb, _md5 in sim_plan:
            print(f"  {mb:9.1f} MB  {key}")
        if missing_sim:
            print(f"  (no replicate published for: {', '.join(missing_sim)})")
        total_mb += sum(mb for *_, mb, _ in sim_plan)
    else:
        avail = [s for s in states if s in dv.SIMULATION_FILES]
        if avail:
            sim_mb = sum(dv.SIMULATION_FILES[s][1] for s in avail)
            print(f"\n  (skipping {sim_mb / 1024:.2f} GB of EpiHiper simulation output;")
            print( "   add --with-simulations to fetch the transmission networks that")
            print( "   `phylogas paint` consumes.)")

    if args.dry_run:
        print("\n--dry-run: nothing downloaded.")
        return 0

    for st, name, url, _mb, md5 in plan:
        state_dir = dest_root / st
        state_dir.mkdir(parents=True, exist_ok=True)
        target = state_dir / name
        plain = target.with_suffix("")
        if plain.exists():
            print(f"    exists (decompressed), skipping: {plain.name}")
            continue
        try:
            dv.download_url(url, target, expect_md5=None if args.no_verify else md5)
        except SystemExit:
            raise
        except Exception as exc:
            print(f"    ERROR downloading {name}: {exc}", file=sys.stderr)
            return 1
        if args.decompress:
            dv.decompress_xz(target, keep=args.keep_compressed)

    # EpiHiper simulation replicates (Zenodo).
    for st, key, url, _mb, md5 in sim_plan:
        state_dir = dest_root / st
        state_dir.mkdir(parents=True, exist_ok=True)
        try:
            dv.download_url(url, state_dir / key, expect_md5=None if args.no_verify else md5)
        except SystemExit:
            raise
        except Exception as exc:
            print(f"    ERROR downloading {key}: {exc}", file=sys.stderr)
            return 1

    # USDA rural-urban continuum codes. A pinned copy ships in data/; this
    # refreshes it from the agency.
    if args.with_rucc:
        dest_root.mkdir(parents=True, exist_ok=True)
        target = dest_root / dv.RUCC_FILENAME
        print(f"\nFetching RUCC codes from USDA ERS ...")
        try:
            dv.download_url(dv.RUCC_URL, target, expect_md5=None)
            print(f"    NOTE: USDA revises this file in place, so no checksum is")
            print(f"          pinned. Re-running an old analysis against a refreshed")
            print(f"          file may not reproduce; see data/README.md.")
        except Exception as exc:
            print(f"    ERROR fetching RUCC: {exc}", file=sys.stderr)

    # Seed sequences (Cov-Spectrum). Public data, so no redistribution
    # constraint; fetched on demand rather than committed.
    if args.with_seeds:
        for st in states:
            print(f"\nFetching seed sequences for {st} ...")
            seed_dir = dest_root / st / "seeds"
            seed_dir.mkdir(parents=True, exist_ok=True)
            sub = _cfg_or_flag(cfg, "variant.include_sublineages", args.include_sublineages, True)
            seed_cmd = [
                sys.executable, "-u", str(_PKG / "seqprep" / "seq_prep.py"),
                "--state", _STATE_NAMES.get(st, st),
                "--pango", str(_cfg_or_flag(cfg, "variant.pango", args.pango, "B.1.617.2")),
                "--output_folder", str(seed_dir),
                "--outlier_method", str(_cfg_or_flag(cfg, "seeds.outlier_method", None, "chaining")),
                "--seed_mode",
                "--include_sublineages" if sub else "--no_include_sublineages",
            ]
            if _cfg_or_flag(cfg, "seeds.insecure_download", args.insecure_download, False):
                seed_cmd.append("--insecure_download")
            if _run(seed_cmd) != 0:
                print(f"    WARNING: seed acquisition failed for {st}", file=sys.stderr)

    # Build the demographics table the pipeline actually reads.
    if args.build_demographics:
        fips = args.fips or _find_fips(dest_root)
        if fips is None:
            print(
                "\nWARNING: no county FIPS lookup found; skipping demographics build.\n"
                "         Pass --fips <county_fips.csv>, or place it at data/county_fips.csv",
                file=sys.stderr,
            )
        else:
            for st in states:
                print(f"\nBuilding demographics for {st} ...")
                d = dest_root / st
                out = d / f"{st}_2_4_0_demographics.csv"
                if out.exists():
                    print(f"  exists, skipping: {out.name}")
                    continue
                try:
                    from .popprep.merge_persontrait import build_demographics
                    build_demographics(
                        persontrait=str(_pick(d, f"{st}_persontrait_epihiper.txt")),
                        person=str(_pick(d, f"{st}_person.csv")),
                        household=str(_pick(d, f"{st}_household.csv")),
                        residence=str(_pick(d, f"{st}_residence_locations.csv")),
                        fips=str(fips),
                        out=str(out),
                    )
                except SystemExit as exc:
                    print(f"  demographics build failed: {exc}", file=sys.stderr)

    print("\nStill required from elsewhere (see docs/data_acquisition.md):")
    print("  - Ruralurbancontinuumcodes2023.csv  (USDA ERS)")
    if not args.with_simulations:
        print("  - EpiHiper simulation replicates    (--with-simulations)")
    print(f"\nVerify with:  phylogas validate-config --config {args.config or 'config.yaml'}")
    return 0


def _pick(d: Path, stem: str) -> Path:
    """Return the decompressed file if present, else the .xz (pandas reads both)."""
    plain = d / stem
    return plain if plain.exists() else d / (stem + ".xz")


def _find_fips(root: Path):
    for cand in (root / "county_fips.csv", root / "Data" / "county_fips.csv",
                 Path("data/county_fips.csv")):
        if cand.exists():
            return cand
    return None


# --------------------------------------------------------------------------
# parser
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# phylogas benchmark / compare-strategies
# --------------------------------------------------------------------------
def cmd_benchmark(args) -> int:
    """Score inferred phylodynamics against the simulation's ground truth."""
    from .benchmark import runner

    if args.what == "mugration":
        runner.benchmark_mugration(
            truth_json=args.truth, samples=args.samples, infections=args.infections,
            out_csv=args.out, state_col=args.state_col,
            duration_years=args.duration_years, save_matrices=args.save_matrices,
        )
    elif args.what == "sequence":
        runner.benchmark_sequence(
            painted_fasta=args.painted, infections=args.infections,
            out_csv=args.out, max_pairs=args.max_pairs,
        )
    elif args.what == "compare":
        runner.benchmark_compare(
            truth_json=args.truth, simulated_json=args.simulated,
            augur_json=args.augur, out_csv=args.out, state_col=args.state_col,
        )
    return 0


def cmd_compare_strategies(args) -> int:
    """Run a BeyondBaseline sweep, then score every strategy against ABM truth.

    Shells out to `beyond-baseline-sweep` (consistent with how `phylogas paint`
    invokes the painter), then scores all saved sample sets in one pass so the
    transmission graph is built once rather than per strategy.
    """
    from .benchmark import runner

    if shutil.which("beyond-baseline-sweep") is None:
        sys.exit(
            "ERROR: `beyond-baseline-sweep` not on PATH.\n"
            "  PhyloGAS depends on BeyondBaseline for sample selection:\n"
            "    pip install git+https://github.com/NSSAC/BeyondBaseline.git"
        )

    cfg = _load_config(args) if args.config else None
    outdir = Path(_cfg_or_flag(cfg, "sampling.outdir", args.outdir, "results/03_sampled_datasets"))
    outdir.mkdir(parents=True, exist_ok=True)
    algorithms = args.algorithms or _cfg_or_flag(cfg, "sampling.algorithms", None, ["surs"])
    batch = _cfg_or_flag(cfg, "sampling.batch_size", args.batch_size, 400)
    linelist = _cfg_or_flag(cfg, "ascertainment.output", args.linelist, None)
    population = _cfg_or_flag(cfg, "population.demographics_file", args.population, None)

    for name, val in (("linelist", linelist), ("population", population)):
        if not val:
            sys.exit(f"ERROR: --{name} not supplied and not resolvable from the config.")

    # ---- 1. selection (BeyondBaseline) ----------------------------------
    sweep = [
        "beyond-baseline-sweep",
        "--linelist", str(linelist),
        "--population", str(population),
        "--outdir", str(outdir),
        "--batch-size", str(batch),
        "--seed", str(_cfg_or_flag(cfg, "random_seed", args.seed, 42)),
        "--save-samples",
        "--algorithms", *[str(a) for a in algorithms],
    ]
    if args.no_replacement:
        sweep.append("--no-replacement")
    if args.stratifiers:
        sweep += ["--stratifiers", *args.stratifiers]
    if args.sweep_args:
        sweep += args.sweep_args

    if not args.skip_sweep:
        rc = _run(sweep, args.dry_run)
        if rc != 0:
            return rc
    else:
        print("--skip-sweep: reusing existing sample sets in the output directory")

    if args.dry_run:
        print(f"+ phylogas benchmark mugration --truth {args.truth} "
              f"--samples '{outdir}/*_samples.csv.xz' --infections {args.infections}")
        return 0

    # ---- 2. scoring (PhyloGAS), one graph build for all strategies -------
    print("\n" + "=" * 60)
    print("Scoring every strategy against ABM truth")
    print("=" * 60)
    df = runner.benchmark_mugration(
        truth_json=args.truth,
        samples=[str(outdir / "*_samples.csv.xz"), str(outdir / "*_samples.csv")],
        infections=args.infections,
        out_csv=args.out or str(outdir / "strategy_comparison.csv"),
        state_col=args.state_col,
    )
    print("\nRanking (best cosine similarity first):")
    print(df[["label", "cosine_similarity", "topological_f1", "n_samples"]].to_string(index=False))
    return 0


# --------------------------------------------------------------------------
# phylogas status  -- the "where am I / what next" entry point
# --------------------------------------------------------------------------
# Stage table: (key, label, config key holding the artefact, next command).
# `config key` may be None for stages whose output is not a single file.
_STAGES = [
    ("config",  "Configuration",        None,
     "cp config.template.yaml config.yaml"),
    ("pop",     "Synthetic population", "population.demographics_file",
     "phylogas fetch-data --states {state} --with-simulations"),
    ("map",     "Entropy map",          "genetic_painter.probability_matrix",
     "phylogas train --config {config}"),
    ("seeds",   "Seed sequences",       "genetic_painter.seed_fasta",
     "phylogas prep-seeds --config {config}"),
    ("abm",     "EpiHiper output",      "epihiper.output_csv",
     "phylogas fetch-data --states {state} --with-simulations"),
    ("painted", "Painted genomes",      None,
     "phylogas paint --config {config}"),
]


def cmd_status(args) -> int:
    """Report pipeline readiness and print the single next command to run."""
    cfg_path = Path(args.config)
    print("PhyloGAS pipeline status\n" + "=" * 52)

    if not cfg_path.is_file():
        print(f"  [ -- ] Configuration            {cfg_path} not found")
        print("\nNext step:")
        print(f"  cp config.template.yaml {cfg_path}")
        print(f"  $EDITOR {cfg_path}")
        return 1

    cfg = _load_config(args)
    state = cfg.get("population.state", default="va")
    print(f"  [okay] Configuration            {cfg_path}")

    first_missing = None
    for _key, label, cfgkey, nxt in _STAGES[1:]:
        if cfgkey is None:
            # Derived stage: check the painter's declared output prefix.
            prefix = cfg.get("genetic_painter.output_prefix", default=None)
            if prefix is None:
                continue
            comp = cfg.get("genetic_painter.compression", default="xz")
            ext = {"xz": ".xz", "bgzf": ".gz"}.get(comp, "")
            path = cfg.resolve_variant(f"{prefix}.sequences.fasta{ext}")
            shown = path or Path(f"{prefix}.sequences.fasta{ext}")
        else:
            raw = cfg.get(cfgkey, default=None)
            if raw is None:
                print(f"  [ -- ] {label:28s} (not configured)")
                continue
            path = cfg.resolve_variant(raw)
            shown = path or Path(str(raw)).expanduser()

        if path is not None:
            size = path.stat().st_size / 1048576
            print(f"  [okay] {label:28s} {size:8.1f} MB  {path}")
        else:
            print(f"  [MISS] {label:28s} {shown}")
            if first_missing is None:
                first_missing = nxt.format(state=state, config=cfg_path)

    print()
    if first_missing is None:
        print("Everything is in place. Run the full pipeline with:")
        print(f"  phylogas run --config {cfg_path} --cores all")
        return 0

    print("Next step:")
    print(f"  {first_missing}")
    print(f"\nFull input report:  phylogas validate-config --config {cfg_path}")
    return 1


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="phylogas",
        description="PhyloGAS: digital-twin benchmarking for genomic surveillance.",
        epilog="See command_map.md for the mapping from the older standalone scripts.",
    )
    p.add_argument("--version", action="version", version=f"phylogas {__version__}")
    sub = p.add_subparsers(dest="command", metavar="<command>")

    # -- status (start here) ----------------------------------------------
    sp = sub.add_parser("status",
                        help="Show pipeline readiness and the next command to run")
    sp.add_argument("--config", "-c", default="config.yaml")
    sp.set_defaults(func=cmd_status)

    def add_common(sp):
        sp.add_argument("--config", "-c", default=None, metavar="YAML",
                        help="Pipeline configuration file (default: config.yaml)")
        sp.add_argument("--dry-run", "-n", action="store_true",
                        help="Print the command that would run, then stop")

    # -- train ------------------------------------------------------------
    sp = sub.add_parser("train", help="Stage 1: build the entropy map from an MSA")
    add_common(sp)
    sp.add_argument("--align-fasta", dest="align_fasta")
    sp.add_argument("--threshold-file", dest="threshold_file")
    sp.add_argument("--base-threshold-df", dest="base_threshold_df")
    sp.add_argument("--seed", type=int)
    sp.add_argument("--proportional", action="store_true", default=None)
    sp.set_defaults(func=cmd_train)

    # -- paint ------------------------------------------------------------
    sp = sub.add_parser("paint", help="Stage 2: paint genomes onto the network")
    add_common(sp)
    sp.add_argument("--analysis-type", dest="analysis_type", default="generate_sequence",
                    choices=["entropy_analysis", "generate_sequence", "both"])
    sp.add_argument("--align-fasta", dest="align_fasta")
    sp.add_argument("--seed-fasta", dest="seed_fasta")
    sp.add_argument("--threshold-file", dest="threshold_file")
    sp.add_argument("--base-threshold-df", dest="base_threshold_df")
    sp.add_argument("--input-graph-csv", dest="input_graph_csv")
    sp.add_argument("--output-prefix", dest="output_prefix")
    sp.add_argument("--start-date", dest="start_date")
    sp.add_argument("--start-tick", dest="start_tick", type=int)
    sp.add_argument("--num-ticks", dest="num_ticks", type=int)
    sp.add_argument("--reference")
    sp.add_argument("--compression", choices=["None", "xz", "bgzf"])
    sp.add_argument("--compression-level", dest="compression_level", type=int)
    sp.add_argument("--compression-threads", dest="compression_threads", type=int)
    sp.add_argument("--persontrait-file", dest="persontrait_file")
    sp.add_argument("--add-metadata", dest="add_metadata")
    sp.add_argument("--painted-prefix", dest="painted_prefix")
    sp.add_argument("--linelist-filter", dest="linelist_filter", action="append",
                    default=None, metavar="[LABEL=]FILE",
                    help="Also write a FASTA/metadata pair restricted to the infections "
                         "named in FILE. Repeatable: one output per filter, plus the "
                         "unrestricted set.")
    sp.add_argument("--mutation-model", dest="mutation_model",
                    choices=["rate_limit", "simple", "poor"])
    sp.add_argument("--initial-viral-load", dest="initial_viral_load", type=float)
    sp.add_argument("--seed", type=int)
    sp.add_argument("--proportional", action="store_true", default=None)
    sp.set_defaults(func=cmd_paint)

    # -- prep-seeds -------------------------------------------------------
    sp = sub.add_parser("prep-seeds", help="Stage 0: acquire seed sequences")
    add_common(sp)
    sp.add_argument("--state")
    sp.add_argument("--pango")
    sp.add_argument("--output-folder", dest="output_folder")
    sp.add_argument("--outlier-method", dest="outlier_method",
                    choices=["iqr", "zscore", "chaining"])
    sp.add_argument("--seed-mode", dest="seed_mode", action="store_true")
    sp.add_argument("--include-sublineages", dest="include_sublineages",
                    action="store_true", default=None,
                    help="Include descendant lineages (query 'B.1.617.2*'). Default: ON.")
    sp.add_argument("--no-include-sublineages", dest="include_sublineages",
                    action="store_false",
                    help="Match the named lineage exactly, excluding descendants.")
    sp.add_argument("--insecure-download", dest="insecure_download",
                    action="store_true", default=None,
                    help="Skip TLS verification for the cluster TSV. Workaround for the "
                         "expired clustertracker.gi.ucsc.edu certificate.")
    sp.set_defaults(func=cmd_prep_seeds)

    # -- subset-fasta -----------------------------------------------------
    sp = sub.add_parser("subset-fasta", help="Stage 5: subset FASTA by sampled metadata")
    add_common(sp)
    sp.add_argument("-m", "--metadata", required=True)
    sp.add_argument("-f", "--fasta", required=True)
    sp.add_argument("-o", "--output", required=True)
    sp.set_defaults(func=cmd_subset_fasta)

    # -- run --------------------------------------------------------------
    sp = sub.add_parser("run", help="Run the full pipeline via Snakemake")
    sp.add_argument("--config", "-c", default="config.yaml")
    sp.add_argument("--snakefile", default=None)
    sp.add_argument("--cores", default="all")
    sp.add_argument("--profile", default=None, help="Snakemake profile, e.g. slurm")
    sp.add_argument("--until", default=None, help="Stop after this rule")
    sp.add_argument("--dry-run", "-n", action="store_true")
    sp.add_argument("extra", nargs=argparse.REMAINDER,
                    help="Additional arguments passed through to snakemake")
    sp.set_defaults(func=cmd_run)

    # -- benchmark ---------------------------------------------------------
    sp = sub.add_parser("benchmark",
                        help="Score inferred phylodynamics against ABM ground truth")
    bsub = sp.add_subparsers(dest="what", metavar="<what>", required=True)

    b = bsub.add_parser("mugration",
                        help="Score sampled sets against the ABM mugration truth")
    b.add_argument("--truth", required=True, help="ABM truth traits JSON (linelist_mugration.json)")
    b.add_argument("--samples", required=True, nargs="+",
                   help="Sample CSVs from BeyondBaseline --save-samples (globs allowed)")
    b.add_argument("--infections", required=True, help="ABM all-events transmission file")
    b.add_argument("--out", default=None, help="Output CSV")
    b.add_argument("--state-col", dest="state_col", default="county")
    b.add_argument("--duration-years", dest="duration_years", type=float, default=1.0)
    b.add_argument("--save-matrices", dest="save_matrices", default=None,
                   help="Directory to write each inferred transition matrix as JSON")
    b.set_defaults(func=cmd_benchmark)

    b = bsub.add_parser("sequence",
                        help="Measure parent->child divergence in the painted genomes")
    b.add_argument("--painted", required=True, help="Painted FASTA (.xz/.gz/plain)")
    b.add_argument("--infections", required=True, help="ABM all-events transmission file")
    b.add_argument("--out", default=None)
    b.add_argument("--max-pairs", dest="max_pairs", type=int, default=200000)
    b.set_defaults(func=cmd_benchmark)

    b = bsub.add_parser("compare",
                        help="Simulated parsimony vs augur inference, against the same truth")
    b.add_argument("--truth", required=True)
    b.add_argument("--simulated", required=True, help="Traits JSON from simulated parsimony")
    b.add_argument("--augur", required=True, help="Traits JSON from `augur traits`")
    b.add_argument("--out", default=None)
    b.add_argument("--state-col", dest="state_col", default="county")
    b.set_defaults(func=cmd_benchmark)

    # -- compare-strategies ------------------------------------------------
    sp = sub.add_parser("compare-strategies",
                        help="Run a BeyondBaseline sweep, then rank every strategy")
    sp.add_argument("--config", "-c", default=None)
    sp.add_argument("--truth", required=True, help="ABM truth traits JSON")
    sp.add_argument("--infections", required=True, help="ABM all-events transmission file")
    sp.add_argument("--linelist", default=None)
    sp.add_argument("--population", default=None)
    sp.add_argument("--algorithms", nargs="+", default=None)
    sp.add_argument("--stratifiers", nargs="+", default=None)
    sp.add_argument("--batch-size", dest="batch_size", type=int, default=None)
    sp.add_argument("--no-replacement", dest="no_replacement", action="store_true")
    sp.add_argument("--outdir", default=None)
    sp.add_argument("--out", default=None, help="Comparison CSV")
    sp.add_argument("--state-col", dest="state_col", default="county")
    sp.add_argument("--seed", type=int, default=None)
    sp.add_argument("--skip-sweep", dest="skip_sweep", action="store_true",
                    help="Score existing sample sets without re-running selection")
    sp.add_argument("--dry-run", "-n", dest="dry_run", action="store_true")
    sp.add_argument("--sweep-args", dest="sweep_args", nargs=argparse.REMAINDER, default=None,
                    help="Extra arguments passed through to beyond-baseline-sweep")
    sp.set_defaults(func=cmd_compare_strategies)

    # -- validate-config --------------------------------------------------
    sp = sub.add_parser("validate-config", help="Check the config and its inputs")
    sp.add_argument("--config", "-c", default="config.yaml")
    sp.set_defaults(func=cmd_validate_config)

    # -- fetch-data -------------------------------------------------------
    sp = sub.add_parser("fetch-data",
                        help="Download synthetic population files from the UVA Dataverse")
    sp.add_argument("--config", "-c", default=None)
    sp.add_argument("--states", nargs="+", metavar="ST",
                    help="State codes to fetch (va ca ga ma mn wa). Default: config population.state")
    sp.add_argument("--dest", default=None, help="Destination root (default: config data_dir)")
    sp.add_argument("--with-epihiper-inputs", dest="with_epihiper_inputs", action="store_true",
                    help="Also fetch the EpiHiper contact networks and persontrait "
                         "databases. Nothing in PhyloGAS reads these; only needed to "
                         "run the ABM yourself.")
    sp.add_argument("--decompress", action="store_true",
                    help="Decompress .xz after download. Off by default: pandas reads "
                         ".xz natively and these files expand ~6x.")
    sp.add_argument("--keep-compressed", dest="keep_compressed", action="store_true",
                    help="With --decompress, also keep the .xz original.")
    sp.add_argument("--no-build-demographics", dest="build_demographics",
                    action="store_false", default=True,
                    help="Skip building <state>_2_4_0_demographics.csv after download.")
    sp.add_argument("--fips", default=None, help="county FIPS -> name lookup CSV")
    sp.add_argument("--with-simulations", dest="with_simulations", action="store_true",
                    help="Also fetch the EpiHiper simulation replicates from Zenodo "
                         "(~0.4-0.6 GB per state). These are the transmission networks "
                         "`phylogas paint` consumes.")
    sp.add_argument("--with-rucc", dest="with_rucc", action="store_true",
                    help="Refresh Ruralurbancontinuumcodes2023.csv from USDA ERS. "
                         "A pinned copy already ships in data/.")
    sp.add_argument("--with-seeds", dest="with_seeds", action="store_true",
                    help="Also fetch per-state seed sequences from Cov-Spectrum. "
                         "Needed for states other than VA: the bundled seed FASTA is "
                         "Virginia Delta and will under-cover other states.")
    sp.add_argument("--pango", default=None, help="Lineage for --with-seeds (default: config)")
    sp.add_argument("--include-sublineages", dest="include_sublineages",
                    action="store_true", default=None,
                    help="With --with-seeds: include descendant lineages. Default: ON.")
    sp.add_argument("--no-include-sublineages", dest="include_sublineages",
                    action="store_false",
                    help="With --with-seeds: match the named lineage exactly.")
    sp.add_argument("--insecure-download", dest="insecure_download",
                    action="store_true", default=None,
                    help="With --with-seeds: skip TLS verification for the cluster TSV "
                         "(expired UCSC certificate workaround).")
    sp.add_argument("--no-verify", dest="no_verify", action="store_true",
                    help="Skip MD5 verification of downloaded files.")
    sp.add_argument("--dry-run", "-n", dest="dry_run", action="store_true",
                    help="List what would be downloaded, then stop.")
    sp.set_defaults(func=cmd_fetch_data)

    # -- build-demographics -----------------------------------------------
    sp = sub.add_parser("build-demographics",
                        help="Join population files into the painter's demographics table")
    sp.add_argument("--persontrait", required=True)
    sp.add_argument("--person", required=True)
    sp.add_argument("--fips", required=True)
    sp.add_argument("--out", required=True)
    sp.add_argument("--household", default=None, help="needed for home_latitude/longitude")
    sp.add_argument("--residence", default=None, help="needed for home_latitude/longitude")
    sp.set_defaults(func=cmd_build_demographics)

    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, "command", None):
        parser.print_help()
        print("\nNot sure where to begin?  phylogas status")
        return 1
    return args.func(args) or 0


if __name__ == "__main__":
    raise SystemExit(main())
