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

import os
from pathlib import Path

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


RESULTS       = cfg("results_dir", "results")
PROJECT       = cfg("project_name", "phylogas")

THRESHOLD     = cfg("genetic_painter.entropy_thresholds")
PROBMATRIX    = cfg("genetic_painter.probability_matrix")
ALIGN_FASTA   = cfg("genetic_painter.align_fasta")
PAINT_PREFIX  = cfg("genetic_painter.output_prefix", f"{RESULTS}/01_synthetic_genomes/{PROJECT}")

_COMP         = cfg("genetic_painter.compression", "xz")
_EXT          = {"xz": ".xz", "bgzf": ".gz"}.get(_COMP, "")
PAINTED_FASTA = f"{PAINT_PREFIX}.sequences.fasta{_EXT}"
PAINTED_META  = f"{PAINT_PREFIX}.metadata.tsv{_EXT}"

LINELIST      = cfg("ascertainment.output", f"{RESULTS}/02_simulated_linelists/linelist.csv")
SAMPLE_DIR    = cfg("sampling.outdir", f"{RESULTS}/03_sampled_datasets")
BENCH_DIR     = cfg("benchmark.outdir", f"{RESULTS}/05_benchmarks")
ALGORITHMS    = cfg("sampling.algorithms", ["surs"])


# --------------------------------------------------------------------------
# Targets
# --------------------------------------------------------------------------
rule all:
    """Default target: painted genomes plus every requested sampled subset."""
    input:
        PAINTED_FASTA,
        PAINTED_META,
        expand(f"{SAMPLE_DIR}/{PROJECT}.{{algo}}.fasta{_EXT}", algo=ALGORITHMS),


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
        outlier=cfg("seeds.outlier_method", "chaining"),
        config=workflow.configfiles[0] if workflow.configfiles else "config.yaml",
    shell:
        "phylogas prep-seeds --config {params.config} "
        "--state {params.state} --pango {params.pango} "
        "--outlier-method {params.outlier} --seed-mode "
        "--output-folder {output}"


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
        graph=cfg("epihiper.output_csv"),
        people=cfg("population.persontrait_file"),
    output:
        linelist=LINELIST,
    params:
        households=cfg("population.household_file", ""),
        rucc=cfg("population.rucc_file", ""),
        ascertain=cfg("ascertainment.parameters", ""),
        start_date=cfg("genetic_painter.start_date"),
        start_tick=cfg("genetic_painter.start_tick"),
        stop_tick=(cfg("genetic_painter.start_tick", 0) + cfg("genetic_painter.num_ticks", 0)),
        seed=cfg("random_seed", 42),
        schedule=cfg("ascertainment.schedule_input", ""),
        variant_mode=cfg("ascertainment.variant_mode", 2),
    shell:
        "simulate_linelist --epihiper {input.graph} --people {input.people} "
        "--households {params.households} --rucc {params.rucc} "
        "--ascertain {params.ascertain} --start_date {params.start_date} "
        "--start_tick {params.start_tick} --stop_tick {params.stop_tick} "
        "--out {output.linelist} --seed {params.seed} --output_all_events "
        "--schedule_input {params.schedule} --variant_mode {params.variant_mode}"


# --------------------------------------------------------------------------
# Stage 4: adaptive sampling  (BeyondBaseline)
# --------------------------------------------------------------------------
rule sample_scenarios:
    """Select which cases get sequenced, under a fixed budget.

    Provided by the BeyondBaseline package (installed via environment.yml).
    """
    input:
        linelist=LINELIST,
        population=cfg("population.persontrait_file"),
    output:
        samples=f"{SAMPLE_DIR}/{{algo}}_samples.csv.xz",
    params:
        outdir=SAMPLE_DIR,
        batch=cfg("sampling.batch_size", 400),
        seed=cfg("random_seed", 42),
        norepl="--no-replacement" if cfg("sampling.no_replacement", True) else "",
    shell:
        "scenarios-runner --linelist {input.linelist} "
        "--population {input.population} --outdir {params.outdir} "
        "--batch-size {params.batch} --seed {params.seed} {params.norepl} "
        "--save-samples --algorithms {wildcards.algo}"


# --------------------------------------------------------------------------
# Stage 5: FASTA subsetting
# --------------------------------------------------------------------------
rule subset_fasta:
    """Pull the selected strains out of the full ground-truth FASTA."""
    input:
        fasta=PAINTED_FASTA,
        samples=f"{SAMPLE_DIR}/{{algo}}_samples.csv.xz",
    output:
        fasta=f"{SAMPLE_DIR}/{PROJECT}.{{algo}}.fasta{_EXT}",
    shell:
        "phylogas subset-fasta -m {input.samples} -f {input.fasta} -o {output.fasta}"


# --------------------------------------------------------------------------
# Stage 6: phylodynamic reconstruction  (Nextstrain)
# --------------------------------------------------------------------------
if cfg("nextstrain.enabled", False):

    rule nextstrain_build:
        """Hand the sampled FASTA to an external Nextstrain/ncov workflow."""
        input:
            fasta=f"{SAMPLE_DIR}/{PROJECT}.{{algo}}.fasta{_EXT}",
        output:
            auspice=f"{cfg('nextstrain.outdir')}/{{algo}}/auspice.json",
        params:
            snakefile=cfg("nextstrain.snakefile"),
            configfile=cfg("nextstrain.configfile"),
        threads: 8
        shell:
            "snakemake --snakefile {params.snakefile} "
            "--configfile {params.configfile} --cores {threads} "
            "--rerun-incomplete"


# --------------------------------------------------------------------------
# Stage 7: benchmarking
# --------------------------------------------------------------------------
# Compares the inferred tree against the painter's ground truth
# (topological F1, mugration cosine similarity).
#
# NOT YET IMPLEMENTED: the metric code does not exist in this repository. The
# mugration half lives in BeyondBaseline/scripts/scenarios_simulation/
# mugration_station.py and should be promoted into a `phylogas benchmark`
# subcommand rather than reimplemented here. Left as an explicit stub so the
# gap is visible instead of silently missing.
