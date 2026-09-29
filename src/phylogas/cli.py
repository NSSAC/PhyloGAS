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
        "seed_fasta",
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
        exists = Path(str(raw)).expanduser().exists()
        print(f"  [{'okay' if exists else 'MISS'}] {label:32s} {raw}")
        missing += (not exists)

    if missing:
        print(f"\n{missing} configured input(s) are missing.")
        print("Run `phylogas fetch-data` for Dataverse-hosted inputs.")
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


def cmd_fetch_data(args) -> int:
    """Download the synthetic population files from the UVA Dataverse."""
    from . import dataverse as dv

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
    total_mb = sum(mb for _, _, _, mb in plan)

    print(f"Destination : {dest_root.resolve()}")
    print(f"States      : {', '.join(states)}")
    for st in states:
        print(f"  {st}: {dv.STATE_DOIS[st]}")
    print(f"\nFiles to fetch ({len(plan)}, {total_mb / 1024:.2f} GB compressed):")
    for st, name, fid, mb in plan:
        print(f"  {mb:9.1f} MB  {name}")
    if not args.with_epihiper_inputs:
        skipped = sum(
            mb for st in states for _, mb in dv.EPIHIPER_INPUT_FILES.get(st, {}).values()
        )
        print(f"\n  (skipping {skipped / 1024:.2f} GB of EpiHiper contact networks;")
        print( "   nothing in PhyloGAS reads them. Use --with-epihiper-inputs if you")
        print( "   intend to run the ABM yourself.)")

    if args.dry_run:
        print("\n--dry-run: nothing downloaded.")
        return 0

    for st, name, fid, _mb in plan:
        state_dir = dest_root / st
        state_dir.mkdir(parents=True, exist_ok=True)
        target = state_dir / name
        plain = target.with_suffix("")
        if plain.exists():
            print(f"    exists (decompressed), skipping: {plain.name}")
            continue
        try:
            dv.download_file(fid, target)
        except Exception as exc:
            print(f"    ERROR downloading {name}: {exc}", file=sys.stderr)
            return 1
        if args.decompress:
            dv.decompress_xz(target, keep=args.keep_compressed)

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
    print("  - EpiHiper simulation replicates    (deposit pending)")
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
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="phylogas",
        description="PhyloGAS: digital-twin benchmarking for genomic surveillance.",
        epilog="See command_map.md for the mapping from the older standalone scripts.",
    )
    p.add_argument("--version", action="version", version=f"phylogas {__version__}")
    sub = p.add_subparsers(dest="command", metavar="<command>")

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
                    help="Also fetch the EpiHiper contact networks (7.1 GB for all six "
                         "states). Nothing in PhyloGAS reads these; only needed to run the ABM.")
    sp.add_argument("--decompress", action="store_true",
                    help="Decompress .xz after download. Off by default: pandas reads "
                         ".xz natively and these files expand ~6x.")
    sp.add_argument("--keep-compressed", dest="keep_compressed", action="store_true",
                    help="With --decompress, also keep the .xz original.")
    sp.add_argument("--no-build-demographics", dest="build_demographics",
                    action="store_false", default=True,
                    help="Skip building <state>_2_4_0_demographics.csv after download.")
    sp.add_argument("--fips", default=None, help="county FIPS -> name lookup CSV")
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
        return 1
    return args.func(args) or 0


if __name__ == "__main__":
    raise SystemExit(main())
