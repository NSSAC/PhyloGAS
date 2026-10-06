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
# <repo>/src/phylogas -> <repo>, so bundled data/ can be found when the
# package is installed with `pip install -e .` and run from elsewhere.
_REPO_ROOT = _PKG.parent.parent


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

    # genetic_painter.location / .reference_location are mappings in YAML but
    # the painter takes them as JSON strings. These were never passed, so the
    # painter silently used its own argparse defaults (USA / Virginia / VA)
    # and the config key was dead -- which mattered once TwinSampler started
    # honouring the same values to build matching strain ids.
    def _location_json(dotted, flag_value):
        raw = flag_value if flag_value is not None else (
            cfg.get(dotted, default=None) if cfg is not None else None)
        if raw is None:
            return None
        if isinstance(raw, str):
            return raw            # already JSON, or passed through verbatim
        import json
        return json.dumps(raw, separators=(",", ":"), sort_keys=True)

    location = _location_json(gp + "location", getattr(args, "location", None))
    reference_location = _location_json(gp + "reference_location",
                                        getattr(args, "reference_location", None))

    # The seed FASTA and the training alignment are produced by prep-seeds /
    # fetch-data under their own keys, so resolve them rather than trusting
    # genetic_painter.{seed,align}_fasta to have been edited by hand.
    align_fasta = pick(gp + "align_fasta", args.align_fasta)
    seed_fasta = args.seed_fasta or pick(gp + "seed_fasta", None)
    if cfg is not None:
        if not align_fasta:
            found, _ = _training_fasta_path(cfg, args)
            align_fasta = str(found) if found else None
        found, where = _seed_fasta_path(cfg, args)
        if found is not None and where != "genetic_painter.seed_fasta":
            seed_fasta = str(found)
            print(f"  seed FASTA: using {found}  (via {where})")
        elif found is not None:
            seed_fasta = str(found)

    optional = [
        ("--align_fasta", align_fasta),
        ("--seed_fasta", seed_fasta),
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
        ("--location", location),
        ("--reference_location", reference_location),
    ]
    # Create the directories the painter writes into. It does not make them
    # itself, so `train` failed outright when the configured entropy-map
    # folder did not exist yet.
    written = [pick(gp + "entropy_thresholds", args.threshold_file),
               pick(gp + "probability_matrix", args.base_threshold_df)]
    if args.analysis_type != "entropy_analysis":
        written.append(pick(gp + "output_prefix", args.output_prefix))
    for raw in written:
        if raw and not args.dry_run:
            Path(str(raw)).expanduser().parent.mkdir(parents=True, exist_ok=True)

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

    if pick(gp + "emit_reference", getattr(args, "emit_reference", None), False):
        cmd.append("--emit_reference")

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
# phylogas check-join
# --------------------------------------------------------------------------
def _parse_date(value):
    """ISO date from a string or a date; None when absent or unparseable."""
    import datetime as _dt

    if value is None or value == "":
        return None
    if isinstance(value, _dt.datetime):
        return value.date()
    if isinstance(value, _dt.date):
        return value
    try:
        return _dt.date.fromisoformat(str(value).strip()[:10])
    except ValueError:
        sys.exit(f"ERROR: could not read '{value}' as an ISO date (YYYY-MM-DD).")


def _context_ids(explicit=None):
    """Sequence ids that legitimately have no line-list row.

    The reference and the clade anchor are analysis *context*, not simulated
    infections, so they never appear in the TwinSampler line list. With
    `genetic_painter.emit_reference` on, the painted FASTA carries the
    reference -- and a join check that counted it as unmatched would block a
    build over a record that is supposed to be there.

    Read from cfg/nextstrain/reference_id.txt, which is also what ncov's
    subsampling `include` uses, so the two cannot drift apart.
    """
    if explicit:
        return {i.strip() for i in explicit if i.strip()}
    out = set()
    f = _REPO_ROOT / "cfg" / "nextstrain" / "reference_id.txt"
    if f.is_file():
        for line in f.read_text().splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                out.add(line)
    return out


def _open_text(path):
    """Open a possibly-compressed text file. These stages mix .xz, .gz and plain."""
    import gzip
    import lzma

    p = str(path)
    if p.endswith(".xz"):
        return lzma.open(p, "rt")
    if p.endswith(".gz"):
        return gzip.open(p, "rt")
    return open(p, "r")


def _fasta_ids(path):
    """Sequence ids, in file order. The id is the first whitespace-delimited
    token, because the painter appends metadata after a space."""
    out = []
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                out.append(line[1:].strip().split(None, 1)[0])
    return out


def _metadata_ids(path, column):
    """Values of one metadata column, in file order, plus the column list.

    The delimiter is sniffed from the header rather than taken from the
    extension, since these files are .tsv in one stage and .csv in another.
    """
    with _open_text(path) as fh:
        first = fh.readline().rstrip("\n").rstrip("\r")
        delim = "\t" if first.count("\t") >= first.count(",") else ","
        cols = first.split(delim)
        if column not in cols:
            return None, cols
        idx = cols.index(column)
        vals = []
        for line in fh:
            parts = line.rstrip("\n").rstrip("\r").split(delim)
            if len(parts) > idx:
                v = parts[idx].strip()
                if v:
                    vals.append(v)
    return vals, cols


def cmd_check_join(args) -> int:
    """Verify that every FASTA header has a metadata row.

    Nextstrain joins sequences to metadata on the `strain` string, and the two
    sides are built by different code in different repositories -- the painter
    writes the FASTA headers, TwinSampler writes the line list. A one-character
    disagreement makes ncov drop sequences, and depending on the step it may do
    so without an error, leaving a tree that is quietly missing data.

    So this is a precondition for a build, not a diagnostic.
    """
    fasta = Path(args.fasta)
    meta = Path(args.metadata)
    for f in (fasta, meta):
        if not f.exists():
            sys.exit(f"ERROR: not found: {f}")

    headers = _fasta_ids(fasta)
    header_set = set(headers)

    meta_vals, cols = _metadata_ids(meta, args.column)
    if meta_vals is None:
        sys.exit(f"ERROR: {meta} has no '{args.column}' column.\n"
                 f"       columns: {', '.join(cols[:12])}"
                 f"{' ...' if len(cols) > 12 else ''}")
    meta_ids = set(meta_vals)

    context = _context_ids(args.ignore_ids)
    exempt = header_set & context
    missing = (header_set - meta_ids) - context
    extra = meta_ids - header_set
    matched = len(header_set) - len(missing) - len(exempt)
    denom = len(header_set) - len(exempt)
    pct = 100.0 * matched / denom if denom else 0.0

    print(f"FASTA    : {fasta}")
    print(f"           {len(headers):,} records, {len(header_set):,} distinct ids")
    if len(headers) != len(header_set):
        print(f"           WARNING: {len(headers) - len(header_set):,} duplicate headers")
    print(f"Metadata : {meta}  (column '{args.column}', {len(meta_ids):,} distinct ids)")
    if len(meta_vals) != len(meta_ids):
        # ncov expects one metadata row per strain. An all-events line list has
        # several rows per person (E, P, I, ...), so it needs de-duplicating
        # before it can serve as metadata.
        print(f"           WARNING: {len(meta_vals) - len(meta_ids):,} duplicate "
              f"'{args.column}' values ({len(meta_vals):,} rows). ncov expects "
              f"one row per sequence.")
    print(f"\nMatched  : {matched:,} / {denom:,}  ({pct:.2f}%)")
    if exempt:
        print(f"Context  : {len(exempt):,} exempt (reference/clade anchor, no "
              f"line-list row by design): {', '.join(sorted(exempt))}")
    print(f"Unmatched: {len(missing):,} sequences with no metadata row")
    print(f"Spare    : {len(extra):,} metadata rows with no sequence  (harmless; ncov filters them)")

    if missing:
        print("\nSequences ncov would drop:")
        for sid in sorted(missing)[: args.max_report]:
            print(f"  {sid}")
        if len(missing) > args.max_report:
            print(f"  ... and {len(missing) - args.max_report:,} more")
        # A total mismatch is nearly always the geography or the date anchor,
        # since those shift every id at once.
        if matched == 0 and extra:
            print("\nNothing matched at all. The usual causes, in order:")
            print("  1. country/divisionAbbr differ between the painter's")
            print("     genetic_painter.location and what simulate_linelist was given")
            print("     (compare the prefixes below).")
            print("  2. start_date / start_tick differ between the two runs, which")
            print("     shifts the year in every id.")
            print(f"  FASTA    e.g.  {sorted(header_set)[0]}")
            print(f"  metadata e.g.  {sorted(extra)[0]}")
        print(f"\nFAIL: {len(missing):,} of {denom:,} sequences would be dropped.")
        return 1

    print("\nOK: every sequence has a metadata row.")
    return 0


# --------------------------------------------------------------------------
# phylogas nextstrain-config
# --------------------------------------------------------------------------
_NS_SHARED = ("custom_lat_longs.tsv", "auspice_augmented.json",
              "reference_id.txt", "references_metadata.tsv",
              "references_sequences.fasta")


def _deep_merge(base: dict, over: dict) -> dict:
    """Recursive dict update, the same shape of merge ncov does on its own
    defaults. `over` wins."""
    out = dict(base)
    for k, v in over.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = v
    return out


def _stage(src: Path, dest: Path) -> str:
    """Put a file inside the ncov checkout, hardlinking when possible.

    `nextstrain build <dir>` mounts only that directory, so under the docker
    and singularity runtimes nothing outside it is visible. Pointing at an
    absolute path elsewhere works under `ambient`/`conda` and fails under the
    runtime most likely on a cluster, which is why these are staged rather
    than referenced. Hardlinks because the FASTAs are large.
    """
    import os

    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        dest.unlink()
    try:
        os.link(src, dest)
        how = "hardlink"
    except OSError:
        import shutil
        shutil.copy2(src, dest)
        how = "copy"
    return how


def _resolve_recipe(rid: str) -> str:
    """Canonical form of a recipe id, via BeyondBaseline's registry.

    An id the registry has collapsed -- a target-blind algorithm is named by
    the axes it reads, so "4S-4_LL-P__surs" became "4S__surs" -- is
    translated, and an id it does not know at all is an error naming the
    index. Without the registry installed the id is taken at its word, which
    keeps this usable in a PhyloGAS-only checkout.
    """
    try:
        from scenarios_simulation.recipes import aliases, all_recipes
    except ImportError:
        return rid
    known = all_recipes()
    if rid in known:
        return rid
    alias = aliases().get(rid)
    if alias is not None:
        print(f"  note: {rid} -> {alias} "
              f"(that algorithm ignores the target distribution)")
        return alias
    sys.exit(f"ERROR: {rid!r} is not a recipe id.\n"
             f"       `scenarios-recipes` lists all {len(known)}.")


def cmd_nextstrain_config(args) -> int:
    """Render one ncov config for one arm, stage its inputs, and validate.

    The output is base.yaml plus an `inputs`/`builds` overlay, written inside
    the ncov checkout with paths relative to it. Validation runs first: a build
    that cannot be scored, or whose sequences do not join to its metadata, is
    refused here rather than discovered hours later.
    """
    import yaml

    cfg = _load_config(args)
    ns_dir = Path(str(_cfg_or_flag(cfg, "nextstrain.dir", args.dir, ""))).expanduser()
    if not ns_dir.is_dir():
        sys.exit("ERROR: set nextstrain.dir (or --dir) to your ncov checkout.\n"
                 f"       got: {ns_dir or '(unset)'}")
    if not (ns_dir / "Snakefile").is_file():
        sys.exit(f"ERROR: {ns_dir} has no Snakefile; is it an ncov checkout?")

    base_path = Path(args.base) if args.base else (
        _REPO_ROOT / "cfg" / "nextstrain" / "ncov" / "base.yaml")
    if not base_path.is_file():
        sys.exit(f"ERROR: base config not found: {base_path}")
    base = yaml.safe_load(base_path.read_text()) or {}

    # --- what this arm is -------------------------------------------------
    project = str(cfg.get("project_name", default="phylogas"))
    if args.build_type == "strategy":
        if getattr(args, "algo", None) and not args.recipe:
            sys.exit("ERROR: --algo was replaced by --recipe, which names a "
                     "scenario and an algorithm together.\n"
                     "       An algorithm alone no longer identifies a build: "
                     "1S__surs and 4S__surs are different trees.\n"
                     "       `scenarios-recipes` lists the ids.")
        if not args.recipe:
            sys.exit("ERROR: --recipe is required for --build-type strategy "
                     "(it names which sampled set to build).")
        arm = _resolve_recipe(str(args.recipe))
    else:
        arm = "all_infections"
    build_name = f"{project}_{arm}"

    loc = cfg.get("genetic_painter.location", default={}) or {}
    if isinstance(loc, str):
        import json
        loc = json.loads(loc)
    division = loc.get("division", "Virginia")

    # --- the two files this arm is built from ------------------------------
    if args.build_type == "strategy":
        sample_dir = Path(str(cfg.get("sampling.outdir",
                                      default=f"{cfg.get('results_dir', default='results')}"
                                              "/03_sampled_datasets")))
        comp = cfg.get("genetic_painter.compression", default="xz")
        ext = {"xz": ".xz", "bgzf": ".gz"}.get(comp, "")
        # Both named by the recipe: subset_fasta writes
        # <project>.<recipe>.fasta and BeyondBaseline writes
        # <recipe>_samples.csv.xz, so one id locates the pair.
        aligned = (Path(args.aligned) if args.aligned
                   else sample_dir / f"{project}.{arm}.fasta{ext}")
        metadata = (Path(args.metadata) if args.metadata
                    else sample_dir / f"{arm}_samples.csv.xz")
        scheme = args.scheme or "strategy_focal_context"
    else:
        prefix = cfg.get("genetic_painter.output_prefix", default=None)
        comp = cfg.get("genetic_painter.compression", default="xz")
        ext = {"xz": ".xz", "bgzf": ".gz"}.get(comp, "")
        aligned = Path(args.aligned) if args.aligned else Path(f"{prefix}.sequences.fasta{ext}")
        # Every infection, not only the ascertained ones, so the all-events
        # table rather than the line list. It has several rows per person, so
        # it must be de-duplicated on `strain`; check-join warns if it is not.
        derived = cfg.ascertainment_outputs().get("allevents")
        metadata = Path(args.metadata) if args.metadata else Path(str(derived or ""))
        scheme = args.scheme or "country_17k"

    for label, f in (("sequences", aligned), ("metadata", metadata)):
        if not f or not Path(f).exists():
            sys.exit(f"ERROR: {label} not found for arm '{arm}': {f}\n"
                     f"       pass --aligned/--metadata to override.")

    # --- validate before writing anything ---------------------------------
    cols = base.get("traits", {}).get("default", {}).get("columns", [])
    if "county" not in cols:
        sys.exit("ERROR: base.yaml traits.default.columns must include 'county'.\n"
                 "       phylogas benchmark mugration reads the traits JSON as\n"
                 "       models['county'], so a build without it cannot be scored.")

    print(f"Arm      : {arm}  (build '{build_name}', scheme '{scheme}')")
    print(f"Sequences: {aligned}")
    print(f"Metadata : {metadata}")

    seq_ids = _fasta_ids(aligned)
    meta_vals, meta_cols = _metadata_ids(metadata, args.column)
    if meta_vals is None:
        sys.exit(f"ERROR: {metadata} has no '{args.column}' column.\n"
                 f"       columns: {', '.join(meta_cols[:12])}"
                 f"{' ...' if len(meta_cols) > 12 else ''}")
    context = _context_ids(args.ignore_ids)
    exempt = set(seq_ids) & context
    missing = (set(seq_ids) - set(meta_vals)) - context
    denom = len(set(seq_ids)) - len(exempt)
    print(f"Join     : {denom - len(missing):,} / {denom:,} sequences have a "
          f"metadata row"
          + (f"  ({len(exempt)} context record(s) exempt)" if exempt else ""))
    if exempt and args.build_type != "all_infections":
        # Both inputs then carry the root. data/reference/reference.fasta is
        # byte-identical to ncov's copy, so whichever ncov keeps is the same
        # sequence -- harmless, just redundant.
        print(f"  note: this arm already contains {', '.join(sorted(exempt))}, "
              f"which reference_data also supplies (identical sequence, so "
              f"harmless).")
    if missing:
        print("\nSequences ncov would drop:", file=sys.stderr)
        for sid in sorted(missing)[:5]:
            print(f"  {sid}", file=sys.stderr)
        if len(missing) > 5:
            print(f"  ... and {len(missing) - 5:,} more", file=sys.stderr)
        sys.exit(f"ERROR: {len(missing):,} sequences have no metadata row. "
                 f"Run:\n       phylogas check-join -f {aligned} -m {metadata}")
    if len(meta_vals) != len(set(meta_vals)):
        print(f"  note: metadata has {len(meta_vals) - len(set(meta_vals)):,} "
              f"duplicate '{args.column}' rows; ncov wants one per sequence")

    # --- stage ------------------------------------------------------------
    stage_root = Path(str(cfg.get("nextstrain.stage_subdir", default="data/phylogas")))
    shared_rel = stage_root / "_shared"
    arm_rel = stage_root / arm

    for name in _NS_SHARED:
        src = _REPO_ROOT / "cfg" / "nextstrain" / name
        if src.is_file():
            _stage(src, ns_dir / shared_rel / name)
        else:
            print(f"  WARNING: {src.name} not in cfg/nextstrain/; "
                  f"base.yaml may reference it", file=sys.stderr)

    aligned_rel = arm_rel / Path(aligned).name
    metadata_rel = arm_rel / Path(metadata).name
    how_a = _stage(Path(aligned), ns_dir / aligned_rel)
    how_m = _stage(Path(metadata), ns_dir / metadata_rel)
    print(f"Staged   : {ns_dir / arm_rel}  ({how_a}/{how_m})")

    # --- render -----------------------------------------------------------
    overlay = {
        "inputs": [
            {"name": "reference_data",
             "metadata": str(shared_rel / "references_metadata.tsv"),
             "aligned": str(shared_rel / "references_sequences.fasta")},
            {"name": build_name,
             "metadata": str(metadata_rel),
             "aligned": str(aligned_rel)},
        ],
        "builds": {
            build_name: {
                "subsampling_scheme": scheme,
                # Substituted into the scheme's {division} queries.
                "division": division,
            },
        },
    }
    merged = _deep_merge(base, overlay)

    out = ns_dir / arm_rel / "config.yaml"
    header = (f"# Generated by `phylogas nextstrain-config` -- do not edit.\n"
              f"# arm: {arm}   build: {build_name}   scheme: {scheme}\n"
              f"# base: {base_path}\n"
              f"# ncov deep merges its defaults/parameters.yaml under this file.\n")
    out.write_text(header + yaml.safe_dump(merged, sort_keys=False, default_flow_style=False))
    print(f"Wrote    : {out}")

    runner = str(cfg.get("nextstrain.runner", default="nextstrain"))
    if runner == "nextstrain":
        cmd = f"nextstrain build {ns_dir} --configfile {arm_rel / 'config.yaml'}"
    else:
        # --directory because ncov writes auspice/ and results/ relative to the
        # working directory, not to the Snakefile. --use-conda because every
        # ncov rule carries `conda: config["conda_environment"]`, which is
        # inert without it -- and with it, ncov's own pinned augur/nextclade
        # are used instead of whatever happens to be on PATH.
        prefix = str(cfg.get("nextstrain.conda_prefix", default="") or "")
        cmd = (f"snakemake --snakefile {ns_dir / 'Snakefile'} "
               f"--directory {ns_dir} "
               f"--configfile {ns_dir / arm_rel / 'config.yaml'} "
               f"--use-conda "
               + (f"--conda-prefix {prefix} " if prefix else "")
               + "--cores all")
    print(f"\nRun with:\n  {cmd}")
    return 0


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
        ("epihiper.output_csv", "EpiHiper transmission log"),
        ("population.demographics_file", "demographics (derived)"),
        ("population.fips_file", "county FIPS lookup"),
    ]
    missing = 0
    print("\nInput files:")

    # Resolved rather than looked up directly: these are produced under one
    # config key and read under another (see _seed_fasta_path).
    seed_found, seed_where = _seed_fasta_path(cfg)
    train_found, train_target = _training_fasta_path(cfg)
    for label, found, shown, note in (
        ("training sequences", train_found, train_target, None),
        ("seed FASTA", seed_found,
         seed_found or cfg.get("genetic_painter.seed_fasta", default=None),
         None if seed_where == "genetic_painter.seed_fasta" else seed_where),
    ):
        if found is not None:
            print(f"  [okay] {label:32s} {found}")
            if note:
                print(f"         (found via {note})")
        elif shown is None:
            print(f"  [ -- ] {label:32s} (not configured)")
        else:
            print(f"  [MISS] {label:32s} {shown}")
            missing += 1
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

    # --- benchmark inputs, with auto-resolution ---------------------------
    print("\nBenchmarking inputs:")
    for key, kind, label in (
        ("benchmark.allevents", "allevents", "all-events (ABM truth)"),
        ("benchmark.truth_mugration", "mugration", "mugration truth"),
    ):
        path, why = cfg.resolve_benchmark_input(key, kind)
        if path is None:
            print(f"  [ -- ] {label:32s} skipped -- {why}")
        else:
            import datetime as _dt
            age = _dt.datetime.now() - _dt.datetime.fromtimestamp(path.stat().st_mtime)
            hrs = age.total_seconds() / 3600
            when = f"{hrs:.0f}h ago" if hrs < 48 else f"{hrs/24:.0f}d ago"
            print(f"  [okay] {label:32s} {path}")
            print(f"         {path.stat().st_size/1048576:.1f} MB, modified {when}  [{why}]")

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


# --------------------------------------------------------------------------
# Derived input locations
#
# Several stages write a file under one config key and are read back under
# another: `prep-seeds` writes into `seeds.output_folder`, `fetch-data
# --with-seeds` writes into <data_dir>/<state>/seeds, and the painter reads
# `genetic_painter.seed_fasta`. Rather than have a producer rewrite the user's
# config, the readers resolve the configured key first and then fall back to
# the places the producers actually write to.
# --------------------------------------------------------------------------
def _sanitize_pango(pango) -> str:
    """Match seq_prep.py's filename sanitisation of a lineage name."""
    return str(pango).split(",")[0].strip().replace(".", "_").replace("/", "_")


def _state_name(cfg, args=None) -> str:
    """Full state name as seq_prep.py expects it (e.g. 'Virginia')."""
    name = _cfg_or_flag(cfg, "population.state_name", getattr(args, "state", None), None)
    if name:
        return str(name)
    code = str(_cfg_or_flag(cfg, "population.state", None, "va")).lower()
    return _STATE_NAMES.get(code, code)


def _seed_fasta_path(cfg, args=None):
    """Locate the seed FASTA. Returns ``(path_or_None, note)``.

    ``note`` says where it came from, so `status` can report a file that was
    found somewhere other than the configured key.
    """
    raw = cfg.get("genetic_painter.seed_fasta", default=None)
    if raw is not None:
        found = cfg.resolve_variant(raw)
        if found is not None:
            return found, "genetic_painter.seed_fasta"

    state = _state_name(cfg, args)
    pango = _sanitize_pango(_cfg_or_flag(cfg, "variant.pango", None, "B.1.617.2"))
    stem = f"{state.replace(' ', '_')}_{pango}_seed_sequences.fasta"

    data_dir = str(cfg.get("data_dir", default="data"))
    code = str(cfg.get("population.state", default="va")).lower()
    for folder, note in (
        (cfg.get("seeds.output_folder", default=None), "seeds.output_folder"),
        (f"{data_dir}/{code}/seeds", "fetch-data --with-seeds output"),
        (f"{data_dir}/importations/sequences", "bundled importations folder"),
    ):
        if not folder:
            continue
        found = cfg.resolve_variant(Path(str(folder)) / stem)
        if found is not None:
            return found, note

    # Nothing on disk: report the configured target if there is one, so the
    # caller can print a meaningful path.
    return None, "genetic_painter.seed_fasta" if raw is not None else "not configured"


def _training_window(cfg):
    """``(date_from, date_to)`` for the bulk training-sequence download.

    Defaults to the window the painter simulates, which is what the existing
    runs used. A mutational model is often better trained on a window that
    starts earlier than the simulation; set ``training.date_from`` /
    ``training.date_to`` explicitly to do that.
    """
    import datetime as _dt

    d_from = cfg.get("training.date_from", default=None)
    d_to = cfg.get("training.date_to", default=None)
    if d_from and d_to:
        return str(d_from), str(d_to)

    start = cfg.get("genetic_painter.start_date", default=None)
    ticks = cfg.get("genetic_painter.num_ticks", default=None)
    if not start:
        return (str(d_from) if d_from else None, str(d_to) if d_to else None)
    s = _dt.date.fromisoformat(str(start))
    e = s + _dt.timedelta(days=int(ticks)) if ticks else None
    return (str(d_from) if d_from else s.isoformat(),
            str(d_to) if d_to else (e.isoformat() if e else None))


def _bulk_output_name(state, pango, d_from, d_to, include_sublineages=True) -> str:
    """The filename seq_prep.py's bulk mode writes, so callers can find it."""
    stem = f"{str(state).replace(' ', '_')}_{_sanitize_pango(pango)}"
    span = f"{d_from}_{d_to}" if d_from and d_to else "all-dates"
    name = f"{stem}_{span}.fasta"
    return name if include_sublineages else name.replace(".fasta", "_exact.fasta")


def _training_fasta_path(cfg, args=None):
    """Locate the bulk training alignment. Returns ``(found_or_None, target)``.

    ``target`` is where `fetch-data --with-training-sequences` would put it:
    the configured ``genetic_painter.align_fasta`` if set, otherwise a name
    derived the way seq_prep.py's bulk mode names its output.
    """
    raw = cfg.get("genetic_painter.align_fasta", default=None)
    if raw is not None:
        return cfg.resolve_variant(raw), Path(str(raw)).expanduser()

    state = _state_name(cfg, args).replace(" ", "_")
    pango = _sanitize_pango(_cfg_or_flag(cfg, "variant.pango", None, "B.1.617.2"))
    d_from, d_to = _training_window(cfg)
    span = f"{d_from}_{d_to}" if d_from and d_to else "all-dates"
    data_dir = str(cfg.get("data_dir", default="data"))
    target = Path(f"{data_dir}/training_sequences/{state}_{pango}_{span}.fasta")
    return cfg.resolve_variant(target), target


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

    # Training sequences for the mutational model (Cov-Spectrum bulk mode).
    # seq_prep.py's bulk path is independent of --seed_mode: it queries
    # Cov-Spectrum over a date range and writes one FASTA, which `train` reads
    # as genetic_painter.align_fasta.
    if args.with_training_sequences:
        if cfg is None:
            print("\nERROR: --with-training-sequences needs --config to know the "
                  "date range and the target path.", file=sys.stderr)
            return 1
        found, target = _training_fasta_path(cfg, args)
        d_from, d_to = _training_window(cfg)
        print(f"\nTraining sequences (Cov-Spectrum bulk mode):")
        print(f"  target: {target}")
        if d_from and d_to:
            src = ("training.date_from/date_to" if cfg.get("training.date_from", default=None)
                   else "genetic_painter.start_date + num_ticks; set "
                        "training.date_from/date_to to widen")
            print(f"  window: {d_from} .. {d_to}   ({src})")
        else:
            print("  window: all dates (no training window or start_date configured)")
        if found is not None:
            print(f"  already present: {found}")

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

    # USDA rural-urban continuum codes. TwinSampler bundles the 2023 table;
    # this fetches the agency's current file for use via population.rucc_file.
    if args.with_rucc:
        dest_root.mkdir(parents=True, exist_ok=True)
        target = dest_root / dv.RUCC_FILENAME
        print(f"\nFetching RUCC codes from USDA ERS ...")
        try:
            dv.download_url(dv.RUCC_URL, target, expect_md5=None)
            print(f"    To use it, set population.rucc_file: {target}")
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

    # Training sequences: the download itself (plan printed further up).
    if args.with_training_sequences:
        found, target = _training_fasta_path(cfg, args)
        d_from, d_to = _training_window(cfg)
        if found is None:   # already reported as present in the plan above
            print(f"\nFetching training sequences for {_state_name(cfg, args)} ...")
            target.parent.mkdir(parents=True, exist_ok=True)
            state = _state_name(cfg, args)
            pango = str(_cfg_or_flag(cfg, "variant.pango", args.pango, "B.1.617.2"))
            sub = _cfg_or_flag(cfg, "variant.include_sublineages",
                               args.include_sublineages, True)
            bulk_cmd = [
                sys.executable, "-u", str(_PKG / "seqprep" / "seq_prep.py"),
                "--state", state,
                "--pango", pango,
                "--output_folder", str(target.parent),
                "--include_sublineages" if sub else "--no_include_sublineages",
            ]
            if d_from and d_to:
                bulk_cmd += ["--date_from", d_from, "--date_to", d_to]
            if _run(bulk_cmd) != 0:
                print("    WARNING: training-sequence download failed",
                      file=sys.stderr)
            else:
                # Bulk mode names its own output; move it onto the configured
                # align_fasta path so `train` finds it without config edits.
                produced = target.parent / _bulk_output_name(state, pango, d_from, d_to, sub)
                if produced.exists() and produced.resolve() == target.resolve():
                    print(f"    saved: {target}")
                elif produced.exists():
                    if target.suffix in (".xz", ".gz", ".bz2", ".zst"):
                        # Bulk mode writes plain FASTA; renaming it onto a
                        # compressed name would misreport the format.
                        print(f"    wrote {produced}")
                        print(f"    NOTE: genetic_painter.align_fasta names a "
                              f"{target.suffix} file. Compress it, or point the "
                              f"key at {produced}.")
                    else:
                        produced.replace(target)
                        # seq_prep.py has just announced its own filename, so
                        # be explicit that the file moved and why.
                        print(f"    renamed {produced.name} -> {target}")
                        print(f"            (genetic_painter.align_fasta, which "
                              f"`phylogas train` reads)")
                if cfg.get("genetic_painter.align_fasta", default=None) is None:
                    print(f"    NOTE: set genetic_painter.align_fasta to {target}")

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

    # Only name what is actually absent. RUCC is not listed: TwinSampler
    # bundles it, and an explicitly configured override is checked by
    # Snakemake as an input.
    outstanding = []
    if not args.with_simulations:
        outstanding.append("- EpiHiper simulation replicates    (--with-simulations)")
    if cfg is not None:
        found, target = _training_fasta_path(cfg, args)
        if found is None and not args.with_training_sequences:
            outstanding.append(
                "- training sequences for the mutational model "
                "(--with-training-sequences)")
    if outstanding:
        print("\nStill required from elsewhere (see docs/data_acquisition.md):")
        for line in outstanding:
            print(f"  {line}")
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
    # Imported per branch, not up front: `benchmark clock` needs neither
    # networkx nor the mugration machinery, and requiring them would make a
    # clock estimate impossible in an environment that only has augur.
    if args.what == "truth":
        from .benchmark.truth_runner import score_samples
        score_samples(
            samples_globs=args.samples, infections=args.infections,
            linelist=args.linelist, date_field=args.date_field,
            out_csv=args.out, stride_weeks=args.stride_weeks,
        )
    elif args.what == "mugration":
        from .benchmark import runner
        runner.benchmark_mugration(
            truth_json=args.truth, samples=args.samples, infections=args.infections,
            out_csv=args.out, state_col=args.state_col,
            duration_years=args.duration_years, save_matrices=args.save_matrices,
        )
    elif args.what == "sequence":
        from .benchmark import runner
        runner.benchmark_sequence(
            painted_fasta=args.painted, infections=args.infections,
            out_csv=args.out, max_pairs=args.max_pairs,
        )
    elif args.what == "compare":
        from .benchmark import runner
        runner.benchmark_compare(
            truth_json=args.truth, simulated_json=args.simulated,
            augur_json=args.augur, out_csv=args.out, state_col=args.state_col,
        )
    elif args.what == "clock":
        return _benchmark_clock(args)
    return 0


def _benchmark_clock(args) -> int:
    """Measure, infer and record the molecular clock rate.

    Three quantities, kept apart on purpose. docs/clock_modes.md says which
    answers which question; the short version is that comparing TreeTime's
    simulated estimate directly to its real-data estimate conflates the
    simulation's physics with TreeTime's inference error, and mu_truth is what
    separates them.

    The arm matters as much as the mode: root-to-tip regression depends on
    which tips are present, so a rate fitted on a sampled arm is partly a
    property of the sampler, not of the simulation.
    """
    import csv

    from .benchmark import clock as clk

    # --config now defaults to config.yaml. An explicitly named file that is
    # missing is still an error; the default quietly degrades to no config so
    # the command keeps working outside a project directory.
    if args.config and (Path(args.config).exists() or args.config != "config.yaml"):
        cfg = _load_config(args)
    else:
        cfg = None
    arm = args.arm or "all_infections"
    want = (["truth", "inferred", "operational"] if args.mode == "all"
            else [args.mode])

    def _cfgget(dotted, default=None):
        return cfg.get(dotted, default=default) if cfg is not None else default

    comp = _cfgget("genetic_painter.compression", "xz")
    ext = {"xz": ".xz", "bgzf": ".gz"}.get(comp, "")
    prefix = _cfgget("genetic_painter.output_prefix")

    results = []

    # --- mu_truth: no tree, no inference --------------------------------
    if "truth" in want:
        fasta = args.fasta or (f"{prefix}.sequences.fasta{ext}" if prefix else None)
        meta = args.metadata or (f"{prefix}.metadata.tsv{ext}" if prefix else None)
        ref = args.reference or _cfgget("genetic_painter.reference_fasta")
        missing = [n for n, v in (("--fasta", fasta), ("--metadata", meta),
                                  ("--reference", ref)) if not v or not Path(v).exists()]
        if missing:
            print(f"  truth: skipped, need {', '.join(missing)}", file=sys.stderr)
        else:
            # The reference is never a data point: its divergence from itself
            # is structurally 0, and it carries Wuhan-Hu-1's real 2019-12-26
            # collection date, ~500 days before any simulated infection. Left
            # in, it reports a 537-day span for a 70-day study.
            drop = set() if args.keep_context else _context_ids(args.exclude_id)

            # Bound the window to the simulated period unless told otherwise.
            # Only for simulated data -- clipping real sequences to the
            # simulation's start date would be meaningless.
            dmin, dmax = _parse_date(args.date_min), _parse_date(args.date_max)
            if dmin is None and args.kind == "simulated":
                dmin = _parse_date(_cfgget("genetic_painter.start_date"))
                if dmin is not None:
                    print(f"         window: >= {dmin} "
                          f"(genetic_painter.start_date; --date-min to override)")

            print(f"  truth: {fasta}")
            res = clk.root_to_tip(
                fasta, ref, meta,
                date_col=args.date_field, id_col=args.column,
                max_records=(args.max_records
                             if args.max_records is not None
                             else _cfgget("nextstrain.clock.max_records")),
                reservoir_seed=int(_cfgget("random_seed", 0) or 0),
                date_min=dmin, date_max=dmax, exclude_ids=drop,
                # Same estimator either way; only the label differs. On real
                # sequences this is the cheapest and most comparable number
                # available -- no TreeTime on either side -- and it needs only
                # a `strain` and a `date` column.
                quantity="mu_truth" if args.kind == "simulated" else "mu_real_observed")
            results.append(clk.format_row(arm, res))

            # --- decompose the pooled rate ---------------------------
            # mu_truth mixes the painter's per-lineage rate with the
            # importation regime: the painter starts chains from real seed
            # genomes that already carry ~42 substitutions of their own, and
            # accumulates along time since each chain's import rather than
            # along calendar date. Pooled is still the right number to
            # compare against mu_sim, which fits the same shape on the same
            # Wuhan root -- but on its own it cannot say what rate the
            # painter produced. These rows split it.
            #
            # Simulated data only: real sequences have no component table and
            # no real_strain, so there is nothing to join.
            if args.decompose and args.kind == "simulated":
                from .benchmark import lineage as lin

                comp_path = args.components
                if comp_path is None and cfg is not None:
                    try:
                        comp_path = cfg.ascertainment_outputs().get("allevents")
                    except Exception:
                        comp_path = None
                comps = None
                if comp_path and Path(comp_path).exists():
                    print(f"  chains: {comp_path}")
                    try:
                        comps = lin.read_components(comp_path)
                        print(f"          {len(comps):,} infections mapped "
                              f"to chains")
                    except ValueError as exc:
                        print(f"  WARNING: {exc}", file=sys.stderr)
                else:
                    print("  chains: not found, so mu_lineage is skipped "
                          "(pass --components <allevents csv>)",
                          file=sys.stderr)
                try:
                    for row in lin.decompose(
                            fasta, ref, meta, components=comps,
                            date_col=args.date_field, id_col=args.column,
                            date_min=dmin, date_max=dmax, exclude_ids=drop,
                            takeoff_min=args.takeoff_min,
                            max_records=(args.max_records
                                         if args.max_records is not None
                                         else _cfgget(
                                             "nextstrain.clock.max_records"))):
                        results.append(clk.format_row(arm, row))
                except Exception as exc:          # noqa: BLE001
                    # The pooled number is the one the DAG depends on; a
                    # failure here must not take it down with it.
                    print(f"  WARNING: decomposition failed ({type(exc).__name__}: "
                          f"{exc}); mu_truth is unaffected", file=sys.stderr)

    # --- mu_sim / mu_real: unconstrained augur refine --------------------
    if "inferred" in want:
        build = args.build_dir
        if not build:
            ns_dir = _cfgget("nextstrain.dir", "")
            project = _cfgget("project_name", "phylogas")
            if ns_dir:
                build = str(Path(ns_dir) / "results" / f"{project}_{arm}")
        if not build or not Path(build).is_dir():
            print(f"  inferred: skipped, no ncov build directory "
                  f"(pass --build-dir; expected results/<project>_<arm>)",
                  file=sys.stderr)
        else:
            res = _refine_unconstrained(Path(build), args, clk)
            if res:
                results.append(clk.format_row(arm, res))

    # --- mu_operational: what ncov assumed -------------------------------
    if "operational" in want:
        bl = args.branch_lengths
        if not bl and args.build_dir:
            bl = str(Path(args.build_dir) / "branch_lengths.json")
        if not bl or not Path(bl).exists():
            print("  operational: skipped, no branch_lengths.json "
                  "(pass --branch-lengths)", file=sys.stderr)
        else:
            results.append(clk.format_row(arm, clk.operational_rate(bl)))

    if not results:
        print("\nNothing computed. Every requested quantity was missing its "
              "inputs; see the messages above.", file=sys.stderr)
        return 1

    # --- report ----------------------------------------------------------
    print()
    for r in results:
        if r.get("error"):
            print(f"  {r.get('quantity', '?'):<16} ERROR: {r['error']}")
            continue
        rate = r.get("rate_subs_per_site_per_year")
        r2 = r.get("r_squared")
        line = f"  {r['quantity']:<24} {rate:.4e} subs/site/year"
        if r2 is not None and not isinstance(r2, str):
            line += f"   R2={r2:.3f}"
        if r.get("tips"):
            line += f"   n={r['tips']:,}"
        if r.get("date_span_days") is not None:
            line += f"  span={r['date_span_days']}d"
        if r.get("chains_fitted"):
            line += f"  chains={r['chains_fitted']:,}"
        if r.get("imports") and not r.get("chains_fitted"):
            line += f"  imports={r['imports']:,}"
        print(line)

    # Only measurements go in the results file; failures were reported above.
    measured = [r for r in results if not r.get("error")]
    if not measured:
        print("\nNo quantity could be computed.", file=sys.stderr)
        return 1

    out = args.out or str(Path(str(_cfgget("benchmark.outdir", "results/05_benchmarks")))
                          / "clock_estimates.csv")
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    # Accumulate: one row per (arm, quantity), so arms can be added one at a
    # time and the sampling-bias comparison reads straight off the file.
    existing, fieldnames = [], []
    if Path(out).exists():
        with open(out) as fh:
            rd = csv.DictReader(fh)
            fieldnames = list(rd.fieldnames or [])
            keys = {(r["arm"], r.get("quantity")) for r in measured
                    if r.get("quantity")}
            existing = [r for r in rd if (r.get("arm"), r.get("quantity")) not in keys]
    for r in measured:
        for k in r:
            if k not in fieldnames:
                fieldnames.append(k)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in existing + measured:
            w.writerow(r)
    print(f"\nWrote: {out}  ({len(existing) + len(measured)} rows)")

    by_q = {r.get("quantity"): r for r in measured}
    if "mu_lineage" in by_q and "mu_between_chains" in by_q:
        w = by_q["mu_lineage"].get("within_date_variance_share")
        pooled = by_q.get("mu_truth", {}).get("rate_subs_per_site_per_year")
        lin_r = by_q["mu_lineage"]["rate_subs_per_site_per_year"]
        btw_r = by_q["mu_between_chains"]["rate_subs_per_site_per_year"]
        print(f"\nDecomposition: mu_truth = w*mu_lineage + (1-w)*mu_between_chains")
        if isinstance(w, float):
            print(f"      w = {w:.4f} of the date variance lies within chains")
            print(f"      {w:.4f}*{lin_r:.4e} + {1 - w:.4f}*{btw_r:.4e} "
                  f"= {w * lin_r + (1 - w) * btw_r:.4e}")
        if pooled:
            print(f"      pooled mu_truth = {pooled:.4e}")
        print("      mu_lineage is the painter's rate; pooled mu_truth is the\n"
              "      wave's, and is what mu_sim is comparable to. See\n"
              "      docs/clock_modes.md.")

    if any(r.get("quantity") == "mu_operational" for r in measured):
        print("\nNote: mu_operational was assumed by ncov, not measured. It "
              "cannot support\n      a claim about the simulation's rate -- "
              "see docs/clock_modes.md.")
    return 0


def _refine_unconstrained(build: Path, args, clk):
    """Run augur refine with no --clock-rate over ncov's own intermediates.

    ncov cannot be configured to do this: its refine rule reads
    config["refine"]["clock_rate"] with [] and passes --clock-rate
    unconditionally, so omitting the key leaves its 0.0008 default in force.
    Hence a separate invocation over the files that rule consumes.
    """
    import shutil
    import subprocess
    import tempfile

    tree = build / "tree_raw.nwk"
    aln = build / "filtered.fasta"
    meta = next((build / n for n in ("metadata_adjusted.tsv.xz",
                                     "metadata_adjusted.tsv") if (build / n).exists()),
                None)
    for label, f in (("tree_raw.nwk", tree), ("filtered.fasta", aln)):
        if not f.exists():
            print(f"  inferred: skipped, {build}/{label} not found "
                  f"(has the ncov build run?)", file=sys.stderr)
            return None
    if meta is None:
        print(f"  inferred: skipped, no metadata_adjusted.tsv[.xz] in {build}",
              file=sys.stderr)
        return None
    if shutil.which("augur") is None:
        print("  inferred: skipped, `augur` not on PATH. Inside the Nextstrain "
              "runtime:\n             nextstrain shell <ncov-dir>",
              file=sys.stderr)
        return None

    tmp = Path(tempfile.mkdtemp(prefix="phylogas-clock-"))
    node_data = tmp / "branch_lengths.json"
    cmd = [
        "augur", "refine",
        "--tree", str(tree),
        "--alignment", str(aln),
        "--metadata", str(meta),
        "--output-tree", str(tmp / "tree.nwk"),
        "--output-node-data", str(node_data),
        "--timetree",
        "--coalescent", str(args.coalescent),
        "--date-inference", "marginal",
        "--divergence-unit", "mutations",
        "--date-confidence",
        "--no-covariance",
    ]
    # clock_filter_iqd deviation from the fitted line prunes valid branches
    # when the rate is being fitted rather than assumed, and prunes different
    # amounts per arm. Off by default here; --clock-filter-iqd to re-enable.
    if args.clock_filter_iqd:
        cmd += ["--clock-filter-iqd", str(args.clock_filter_iqd)]

    print(f"  inferred: augur refine (unconstrained) on {build.name}")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-4:]
        print("  inferred: augur refine failed:\n    " + "\n    ".join(tail),
              file=sys.stderr)
        return None

    res = clk.operational_rate(node_data)
    if res.get("error"):
        return res
    res["quantity"] = "mu_sim" if args.kind == "simulated" else "mu_real"
    res["method"] = "augur refine, no --clock-rate (root-to-tip fit by TreeTime)"
    res["node_data"] = str(node_data)
    return res


def cmd_assign_variants(args) -> int:
    """Attach benchmark variant labels to transmission components.

    Moved here from TwinSampler: these labels are ground truth for prevalence
    estimation, matched against a real importation schedule. TwinSampler still
    finds the components; PhyloGAS decides what variant each one is.
    """
    import pandas as pd

    from .benchmark.variants import assign_variants, VARIANT_COLUMN

    print(f"Events   : {args.allevents}")
    events = pd.read_csv(args.allevents, dtype={"alias_pid": str, "alias_contact": str,
                                                "sim_pid": str, "pid": str})
    events.columns = [c.strip() for c in events.columns]
    print(f"  {len(events):,} rows")

    print(f"Schedule : {args.schedule}")
    schedule = pd.read_csv(args.schedule)

    out = assign_variants(events, schedule, mode=args.mode, column=args.column)

    dest = Path(args.out)
    dest.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(dest, index=False,
               compression="xz" if str(dest).endswith(".xz") else None)
    print(f"\nWrote {len(out):,} rows with '{args.column}' -> {dest}")
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
# Pipeline stages in dependency order: each one's inputs are produced by the
# stages above it. `status` walks this list, so the order here is what decides
# the order of the report and which command is named as the next step.
#
#   key      -> internal id
#   label    -> what the report calls it
#   cfgkey   -> config key holding the path, or None for a derived location
#               (see _resolve_stage)
#   nxt      -> the command that produces it
_STAGES = [
    ("config",   "Configuration",         None,
     "cp config.template.yaml config.yaml"),
    ("pop",      "Synthetic population",  "population.demographics_file",
     "phylogas fetch-data --config {config} --states {state} --with-simulations"),
    ("training", "Training sequences",    None,
     "phylogas fetch-data --config {config} --states {state} --with-training-sequences"),
    ("map",      "Entropy map",           "genetic_painter.probability_matrix",
     "phylogas train --config {config}"),
    ("seeds",    "Seed sequences",        None,
     "phylogas prep-seeds --config {config} --seed-mode"),
    ("abm",      "EpiHiper output",       "epihiper.output_csv",
     "phylogas fetch-data --config {config} --states {state} --with-simulations"),
    ("painted",  "Painted genomes",       None,
     "phylogas paint --config {config}"),
]


def _resolve_stage(key, cfgkey, cfg):
    """Locate one stage's artefact. Returns ``(found_or_None, shown, note)``.

    ``shown`` is the path to display when it is missing; ``note`` is an
    optional remark, used when a file turns up somewhere other than the
    configured key.
    """
    if cfgkey is not None:
        raw = cfg.get(cfgkey, default=None)
        if raw is None:
            return None, None, None          # not configured
        return (cfg.resolve_variant(raw),
                Path(str(raw)).expanduser(), None)

    if key == "training":
        found, target = _training_fasta_path(cfg)
        note = None
        if found is None and cfg.get("genetic_painter.align_fasta", default=None) is None:
            note = "align_fasta not set; would be written here"
        return found, target, note

    if key == "seeds":
        found, where = _seed_fasta_path(cfg)
        raw = cfg.get("genetic_painter.seed_fasta", default=None)
        if found is None:
            return None, (Path(str(raw)).expanduser() if raw else None), None
        note = None
        if where != "genetic_painter.seed_fasta":
            note = f"found via {where}; config points elsewhere"
        return found, found, note

    if key == "painted":
        prefix = cfg.get("genetic_painter.output_prefix", default=None)
        if prefix is None:
            return None, None, None
        comp = cfg.get("genetic_painter.compression", default="xz")
        ext = {"xz": ".xz", "bgzf": ".gz"}.get(comp, "")
        raw = f"{prefix}.sequences.fasta{ext}"
        return cfg.resolve_variant(raw), Path(raw), None

    return None, None, None


# Flags the Snakefile passes to simulate_linelist that older TwinSampler
# builds do not accept. argparse rejects an unknown flag with exit status 2
# and "unrecognized arguments", which names the flag but not the cause, so
# `status` probes for it instead of letting a sweep discover it per job.
_TS_REQUIRED_FLAGS = ("--division_abbr",)


def _report_twinsampler() -> None:
    """Report the TwinSampler CLI, and whether it is new enough.

    PhyloGAS and TwinSampler are separate repos installed separately, so an
    install can lag a PhyloGAS pull. Being on PATH is not enough: the
    geography flags arrived in TwinSampler bf8e68b, and without them
    `simulate_linelist` exits 2 before doing any work.
    """
    import shutil
    import subprocess

    label = "TwinSampler"
    exe = shutil.which("simulate_linelist")
    if exe is None:
        print(f"  [MISS] {label:28s} simulate_linelist not on PATH")
        print("         pip install -e /path/to/TwinSampler")
        return
    try:
        helptext = subprocess.run([exe, "--help"], capture_output=True,
                                  text=True, timeout=60).stdout
    except Exception as exc:
        print(f"  [MISS] {label:28s} simulate_linelist not runnable: {exc}")
        return

    stale = [f for f in _TS_REQUIRED_FLAGS if f not in helptext]
    if stale:
        print(f"  [MISS] {label:28s} too old: no {', '.join(stale)}")
        print("         cd /path/to/TwinSampler && git pull && pip install -e .")
    else:
        print(f"  [okay] {label:28s} {exe}")


def _report_nextstrain(cfg) -> None:
    """Report the external Nextstrain dependency: the checkout and the CLI.

    Diagnosed rather than assumed. "It is on PATH" is fine as an interface
    decision, but a silent assumption turns into a confusing failure hours
    into a run, so `status` says which of the three states you are in.
    """
    import shutil
    import subprocess

    ns_dir = str(cfg.get("nextstrain.dir", default="") or "")
    label = "ncov checkout"
    if not ns_dir:
        print(f"  [ -- ] {label:28s} (nextstrain.dir not set)")
    elif not (Path(ns_dir).expanduser() / "Snakefile").is_file():
        print(f"  [MISS] {label:28s} {ns_dir}  (no Snakefile)")
    else:
        print(f"  [okay] {label:28s} {ns_dir}")

    runner = str(cfg.get("nextstrain.runner", default="nextstrain"))
    if runner == "snakemake":
        # No CLI needed: snakemake --use-conda builds ncov's own environment,
        # so what has to exist is snakemake and something to build it with.
        label = "runner: snakemake"
        missing = [t for t in ("snakemake",) if shutil.which(t) is None]
        if not (shutil.which("conda") or shutil.which("mamba")):
            missing.append("conda or mamba")
        if missing:
            print(f"  [MISS] {label:28s} not on PATH: {', '.join(missing)}")
        else:
            print(f"  [okay] {label:28s} snakemake + conda found")
        return

    label = "nextstrain CLI"
    if shutil.which("nextstrain") is None:
        print(f"  [MISS] {label:28s} not on PATH")
        plat = "mac" if sys.platform == "darwin" else "linux"
        print(f"         curl -fsSL --proto '=https' "
              f"https://nextstrain.org/cli/installer/{plat} | bash")
        return
    try:
        ver = subprocess.run(["nextstrain", "--version"], capture_output=True,
                             text=True, timeout=30).stdout.strip()
    except Exception as exc:
        print(f"  [MISS] {label:28s} present but not runnable: {exc}")
        return

    # Installed is not the same as usable: a runtime has to be set up too.
    try:
        chk = subprocess.run(["nextstrain", "check-setup"], capture_output=True,
                             text=True, timeout=180)
        ok = chk.returncode == 0
    except Exception:
        ok = False
    if ok:
        print(f"  [okay] {label:28s} {ver}")
    else:
        print(f"  [MISS] {label:28s} {ver}, but no runtime is set up")
        print("         nextstrain setup --set-default conda   # or docker")


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

    missing = []   # (label, command) in dependency order
    for key, label, cfgkey, nxt in _STAGES[1:]:
        path, shown, note = _resolve_stage(key, cfgkey, cfg)

        if path is None and shown is None:
            print(f"  [ -- ] {label:28s} (not configured)")
            continue

        if path is not None:
            size = path.stat().st_size / 1048576
            print(f"  [okay] {label:28s} {size:8.1f} MB  {path}")
            if note:
                print(f"         {note}")
        else:
            print(f"  [MISS] {label:28s} {shown}")
            if note:
                print(f"         {note}")
            missing.append((label, nxt.format(state=state, config=cfg_path)))

    _report_twinsampler()

    if cfg.get("nextstrain.enabled", default=False):
        _report_nextstrain(cfg)

    print()
    if not missing:
        print("Everything is in place. Run the full pipeline with:")
        print(f"  phylogas run --config {cfg_path} --cores all")
        return 0

    # Report every gap, not just the first. A single "next step" sent people
    # to `train` while the seed FASTA was also absent, so they had to work the
    # rest out by backtracking.
    print("Next step:")
    print(f"  {missing[0][1]}")
    if len(missing) > 1:
        print(f"\nThen, for the other {len(missing) - 1} missing input"
              f"{'s' if len(missing) > 2 else ''}:")
        seen = {missing[0][1]}
        for label, cmd in missing[1:]:
            shown = cmd if cmd not in seen else "(same command as above)"
            print(f"  {label:22s}  {shown}")
            seen.add(cmd)
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

    # -- check-join -------------------------------------------------------
    sp = sub.add_parser("check-join",
                        help="Verify every FASTA header has a metadata row "
                             "(ncov joins on the strain id)")
    sp.add_argument("-f", "--fasta", required=True,
                    help="painted or subset FASTA (.fasta, .xz, .gz)")
    sp.add_argument("-m", "--metadata", required=True,
                    help="metadata table; TwinSampler line list or painter metadata")
    sp.add_argument("--column", default="strain",
                    help="metadata column holding the id (default: strain)")
    sp.add_argument("--max-report", dest="max_report", type=int, default=10,
                    help="how many unmatched ids to print (default: 10)")
    sp.add_argument("--ignore-ids", dest="ignore_ids", nargs="+", default=None,
                    help="ids allowed to have no metadata row. Defaults to "
                         "cfg/nextstrain/reference_id.txt (the reference and "
                         "clade anchor, which are context rather than infections).")
    sp.set_defaults(func=cmd_check_join)

    # -- nextstrain-config -------------------------------------------------
    sp = sub.add_parser("nextstrain-config",
                        help="Render, stage and validate an ncov config for one arm")
    sp.add_argument("--config", "-c", default="config.yaml")
    sp.add_argument("--build-type", dest="build_type", required=True,
                    choices=["strategy", "all_infections"],
                    help="strategy: one BeyondBaseline recipe (needs --recipe). "
                         "all_infections: the full painted set with the 17k cap.")
    sp.add_argument("--recipe", default=None,
                    help="recipe id, for --build-type strategy, as printed by "
                         "`scenarios-recipes` (e.g. 4S__surs)")
    sp.add_argument("--algo", default=None, help=argparse.SUPPRESS)
    sp.add_argument("--dir", default=None, help="ncov checkout (default: nextstrain.dir)")
    sp.add_argument("--base", default=None,
                    help="base config (default: cfg/nextstrain/ncov/base.yaml)")
    sp.add_argument("--aligned", default=None, help="override the sequences file")
    sp.add_argument("--metadata", default=None, help="override the metadata file")
    sp.add_argument("--scheme", default=None, help="override the subsampling scheme")
    sp.add_argument("--column", default="strain",
                    help="metadata column holding the sequence id (default: strain)")
    sp.add_argument("--ignore-ids", dest="ignore_ids", nargs="+", default=None,
                    help="ids allowed to have no metadata row (default: "
                         "cfg/nextstrain/reference_id.txt)")
    sp.set_defaults(func=cmd_nextstrain_config)

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
    sp.add_argument("--emit-reference", dest="emit_reference",
                    action="store_true", default=None,
                    help="Write the reference as the first output record. Off by "
                         "default; Nextstrain supplies its own root.")
    sp.add_argument("--location", default=None,
                    help="JSON location block; overrides genetic_painter.location")
    sp.add_argument("--reference-location", dest="reference_location", default=None,
                    help="JSON location block for the reference record")
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

    b = bsub.add_parser("truth",
                        help="Score sampled sets against ABM ground truth "
                             "(prevalence error, component coverage)")
    b.add_argument("--samples", required=True, nargs="+",
                   help="Sample CSVs from BeyondBaseline --save-samples (globs allowed)")
    b.add_argument("--infections", required=True,
                   help="ABM all-events file, after `phylogas assign-variants`")
    b.add_argument("--linelist", default=None,
                   help="Linelist, if its edges should define coverage instead")
    b.add_argument("--date-field", dest="date_field", default="date")
    b.add_argument("--stride-weeks", dest="stride_weeks", type=int, default=4)
    b.add_argument("--out", default=None, help="Output CSV (AUC_truth_rankings.csv)")
    b.set_defaults(func=cmd_benchmark)

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

    b = bsub.add_parser("clock",
                        help="Molecular clock rate: measured, inferred and assumed")
    # Default to config.yaml as `status` does. Without it every lookup fell
    # through to a hardcoded default, and the CSV landed in
    # results/05_benchmarks rather than the project's results/<name>/.
    b.add_argument("--config", "-c", default="config.yaml")
    b.add_argument("--mode", default="all",
                   choices=["truth", "inferred", "operational", "all"],
                   help="truth: tree-free root-to-tip from the painted sequences. "
                        "inferred: unconstrained augur refine on an ncov build. "
                        "operational: the rate ncov assumed. Default: all.")
    b.add_argument("--arm", default=None,
                   help="which arm these numbers describe (default: all_infections). "
                        "Fitting on a sampled arm measures the sampler too -- see "
                        "docs/clock_modes.md")
    b.add_argument("--kind", default="simulated", choices=["simulated", "real"],
                   help="labels the inferred rate mu_sim or mu_real")
    b.add_argument("--fasta", default=None, help="sequences for --mode truth")
    b.add_argument("--metadata", default=None, help="metadata for --mode truth")
    b.add_argument("--reference", default=None,
                   help="reference FASTA (default: genetic_painter.reference_fasta)")
    b.add_argument("--build-dir", dest="build_dir", default=None,
                   help="ncov results/<build> directory, for --mode inferred")
    b.add_argument("--branch-lengths", dest="branch_lengths", default=None,
                   help="branch_lengths.json, for --mode operational")
    b.add_argument("--date-field", dest="date_field", default="date")
    b.add_argument("--column", default="strain", help="metadata id column")
    b.add_argument("--coalescent", default="skyline")
    b.add_argument("--clock-filter-iqd", dest="clock_filter_iqd", type=float, default=0,
                   help="0 (default) disables it. With the rate being fitted, the "
                        "filter prunes valid branches, and by differing amounts per arm.")
    b.add_argument("--max-records", dest="max_records", type=int, default=None,
                   # Sampled across the file, not the first N -- see
                   # root_to_tip. Falls back to nextstrain.clock.max_records.
                   help="fit on at most N sequences (for a quick look at a "
                        "huge FASTA). Sampled uniformly across the file, so "
                        "the date range is preserved")
    b.add_argument("--components", default=None,
                   help="all-events table carrying component_id, joined on "
                        "alias_pid to split mu_truth into mu_lineage and "
                        "mu_between_chains (default: the configured "
                        "ascertainment allevents output)")
    b.add_argument("--takeoff-min", dest="takeoff_min", type=int, default=10,
                   help="tips a chain needs before it counts as having taken "
                        "off, for the seed-trend rows (default: 10)")
    b.add_argument("--no-decompose", dest="decompose", action="store_false",
                   default=True,
                   help="skip the mu_lineage / seed-trend rows; report the "
                        "pooled mu_truth alone")
    b.add_argument("--date-min", dest="date_min", default=None,
                   help="earliest date to fit (default for --kind simulated: "
                        "genetic_painter.start_date)")
    b.add_argument("--date-max", dest="date_max", default=None,
                   help="latest date to fit (default: unbounded)")
    b.add_argument("--exclude-id", dest="exclude_id", action="append", default=None,
                   help="drop this sequence id; repeatable. Defaults to the ids in "
                        "cfg/nextstrain/reference_id.txt")
    b.add_argument("--keep-context", dest="keep_context", action="store_true",
                   help="keep the reference and other context sequences in the fit")
    b.add_argument("--out", default=None, help="CSV (default: <benchmark.outdir>/clock_estimates.csv)")
    b.set_defaults(func=cmd_benchmark)

    # -- assign-variants ---------------------------------------------------
    sp = sub.add_parser("assign-variants",
                        help="Attach benchmark variant labels to transmission components")
    sp.add_argument("--allevents", required=True,
                    help="TwinSampler all-events file (must carry component_id)")
    sp.add_argument("--schedule", required=True,
                    help="Importation schedule: tick,date,variant,clusters,sample_count")
    sp.add_argument("--mode", default="bipartite", choices=["bipartite", "temporal"],
                    help="Matching strategy (default: bipartite)")
    sp.add_argument("--column", default="variant_benchmark",
                    help="Output column name (default: variant_benchmark)")
    sp.add_argument("--out", required=True, help="Output CSV (.xz supported)")
    sp.set_defaults(func=cmd_assign_variants)

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
                    help="Download the current Ruralurbancontinuumcodes2023.csv from "
                         "USDA ERS into data/. Optional: TwinSampler bundles the 2023 "
                         "table; point population.rucc_file here to use this copy.")
    sp.add_argument("--with-seeds", dest="with_seeds", action="store_true",
                    help="Also fetch per-state seed sequences from Cov-Spectrum. "
                         "Needed for states other than VA: the bundled seed FASTA is "
                         "Virginia Delta and will under-cover other states.")
    sp.add_argument("--with-training-sequences", "--training-sequences",
                    dest="with_training_sequences", action="store_true",
                    help="Also fetch the bulk Cov-Spectrum alignment that `phylogas "
                         "train` reads as genetic_painter.align_fasta. Uses "
                         "seq_prep.py's bulk mode over training.date_from/date_to, "
                         "defaulting to the painter's own window. Skipped if the "
                         "file is already present.")
    sp.add_argument("--pango", default=None,
                    help="Lineage for --with-seeds / --with-training-sequences "
                         "(default: config)")
    sp.add_argument("--include-sublineages", dest="include_sublineages",
                    action="store_true", default=None,
                    help="With --with-seeds / --with-training-sequences: include descendant lineages. Default: ON.")
    sp.add_argument("--no-include-sublineages", dest="include_sublineages",
                    action="store_false",
                    help="With --with-seeds / --with-training-sequences: match the named lineage exactly.")
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
