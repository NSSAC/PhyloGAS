#!/bin/bash
#SBATCH -p bii                   # Partition name
#SBATCH -A bii_nssac             # Account
#SBATCH -c 8                     # Cores for augur tree (IQ-TREE 2) & streaming
#SBATCH -t 02:30:00              # Wall-clock limit (~10-20 min expected)
#SBATCH -J phylogas_phylo        # Job name
#SBATCH --array=0-4              # 5 states: 0=va, 1=ga, 2=ma, 3=mn, 4=wa
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --mem=48G

# Stage 1 (5-state array):
#   1) Score existing full-state Nextstrain Auspice builds (`04_nextstrain_builds/`)
#      against the pruned ABM ground-truth genealogy -> `05_benchmarks/Phylogeny_Metrics.csv`
#      (for VA, also scores any finished Bryan `nova3_jun_aug` JUNIPER runs if present).
#   2) Export a curated `06_microbetrace/` showcase bundle (`N = 600`, `1S__surs`) per state.
#   3) Extract 6 focal genetically-isolated county outbreaks per state (`3 counties x pi={0.25,0.50}`),
#      run `Nextstrain` (`IQ-TREE 2 + TreeTime`) on each, export `dated_tree_years.nwk` for
#      `TransPhylo`, and write initial `MicrobeTrace_1SNP` + `Nextstrain` scores.
#
# Usage:
#   JID=$(sbatch --parsable slurm/run_phylo_bench.sbatch)
#   sbatch --dependency=afterok:${JID} slurm/run_outbreak_bench.sbatch

cd "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
[ -f slurm/prepare_inputs.sh ] || cd ..
[ -f slurm/prepare_inputs.sh ] || {
    echo "submit from the PhyloGAS checkout (or its slurm/ folder); job ran in $(pwd)" \
        > "${SLURM_SUBMIT_DIR:-.}/PHYLOGAS_PHYLO_${SLURM_JOB_ID:-local}_FAILED.txt"
    exit 1
}

STATES=(va ga ma mn wa)
if [ -z "${STATE:-}" ]; then
    IDX="${SLURM_ARRAY_TASK_ID:-0}"
    STATE="${STATES[$IDX]}"
fi

JOB_TAG="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-local}}"
mkdir -p slurm/logs
exec >"slurm/logs/phylo_bench_${STATE}_${JOB_TAG}.out" 2>"slurm/logs/phylo_bench_${STATE}_${JOB_TAG}.err"

set -eo pipefail

source /project/biocomplexity/asw3xp/miniconda3/bin/activate
conda activate phylogas_env
unset PYTHONPATH

set -u

trap 'echo "=== seff ${SLURM_JOB_ID:-} ==="; seff "${SLURM_JOB_ID:-}" 2>/dev/null || true' EXIT

CONFIG="configs/${STATE}.yaml"
RDIR="results/${STATE}_delta_wave"
CPUS="${SLURM_CPUS_PER_TASK:-8}"
echo "state: ${STATE}   config: ${CONFIG}   cpus: ${CPUS}   host: $(hostname)"

# Optional hook for Virginia: if Bryan's nova3_jun_aug JUNIPER runs have finished
# `results.RData` or `nodes.csv`, stage/link them for automatic comparison!
JUN_ARG=()
BRYAN_NOVA="/project/bii_nssac/people/bl4zc/pgcoe_genepi_analytics/juniper/runs/nova3_jun_aug"
if [ "${STATE}" = "va" ] && [ -d "${BRYAN_NOVA}" ]; then
    STAGE_JUN="${RDIR}/05_benchmarks/_bryan_nova3_stage"
    mkdir -p "${STAGE_JUN}"
    for bdir in "${BRYAN_NOVA}"/seq*; do
        [ -d "${bdir}" ] || continue
        bname="$(basename "${bdir}")"
        if [ -f "${bdir}/nodes.csv" ] && [ -f "${bdir}/input_data/sample_lookup.csv" ]; then
            mkdir -p "${STAGE_JUN}/${bname}/input_data"
            cp -f "${bdir}/nodes.csv" "${STAGE_JUN}/${bname}/nodes.csv"
            cp -f "${bdir}/input_data/sample_lookup.csv" "${STAGE_JUN}/${bname}/input_data/sample_lookup.csv"
        fi
    done
    if compgen -G "${STAGE_JUN}/*/nodes.csv" >/dev/null; then
        JUN_ARG=(--juniper-runs-dir "${STAGE_JUN}")
    fi
fi

# 1. Score statewide Nextstrain builds (and any staged JUNIPER sub-studies)
echo "=== [1/3] Scoring statewide Nextstrain Auspice trees for ${STATE} ==="
phylogas benchmark phylogeny --config "${CONFIG}" "${JUN_ARG[@]}"

# 2. Export verified MicrobeTrace showcase bundle (N=600, cluster-enriched on 1S__surs)
echo "=== [2/3] Exporting MicrobeTrace showcase bundle for ${STATE} ==="
PAINTED_FA="${RDIR}/01_synthetic_genomes/${STATE}_delta_wave.sequences.fasta.xz"
SURS_SAMP="${RDIR}/03_sampled_datasets/1S__surs_samples.csv.xz"
ALLEV="${RDIR}/05_benchmarks/allevents_variants.csv.xz"
if [ -f "${PAINTED_FA}" ] && [ -f "${SURS_SAMP}" ] && [ -f "${ALLEV}" ]; then
    phylogas export-microbetrace \
        --fasta "${PAINTED_FA}" \
        --samples "${SURS_SAMP}" \
        --infections "${ALLEV}" \
        --outdir "${RDIR}/06_microbetrace" \
        --prefix "${STATE}_delta_wave_1S__surs_microbetrace" \
        --selection-mode cluster-enriched \
        --max-sequences 600
fi

# 3. Prepare 6 focal outbreak datasets (3 counties x pi={0.25, 0.50}) & run Nextstrain/TreeTime
echo "=== [3/3] Preparing focal outbreak benchmarks & running Nextstrain/TreeTime for ${STATE} ==="
phylogas benchmark outbreaks --config "${CONFIG}" --step prep --n-clusters 3 --sampling-rates 0.25 0.50

for ob_dir in "${RDIR}/07_outbreak_bench/${STATE}_c"*; do
    [ -d "${ob_dir}" ] || continue
    echo "--- Running augur tree + refine on $(basename "${ob_dir}") ---"
    phylogas benchmark outbreaks --config "${CONFIG}" --step nextstrain --workdir "${ob_dir}" --threads "${CPUS}"
    phylogas benchmark outbreaks --config "${CONFIG}" --step score --workdir "${ob_dir}"
done

phylogas benchmark outbreaks --config "${CONFIG}" --step aggregate
echo "=== Completed Stage 1 phylogeny & outbreak staging for ${STATE} ==="
