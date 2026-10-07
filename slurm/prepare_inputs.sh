#!/bin/bash
# Stage and check every input one state needs. Called at the top of each
# slurm/run_<state>.sbatch; can also be run by hand:
#
#   bash slurm/prepare_inputs.sh ma
#
# Safe to re-run: downloads already on disk are verified rather than fetched
# again, and the seeds are rebuilt only when missing or wrong -- rewriting them
# would make Snakemake repaint the state.
#
# Needs outbound internet (Dataverse, Zenodo, Cov-Spectrum).
set -eo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

st="${1:?usage: prepare_inputs.sh <state: ga|ma|mn|va|wa>}"
declare -A NAME=([ga]=Georgia [ma]=Massachusetts [mn]=Minnesota [va]=Virginia [wa]=Washington)
name="${NAME[$st]:?unknown state: $st}"
cfg="configs/${st}.yaml"
[ -f "$cfg" ] || { echo "missing $cfg -- run scripts/make_state_configs.py"; exit 1; }

echo "=== [${st}] inputs: population, EpiHiper replicate, training alignment ==="
# The training alignment's window is the painted window (first Delta import +
# num_ticks). An alignment already on disk is kept as is.
phylogas fetch-data --config "$cfg" --states "$st" \
    --with-simulations --with-training-sequences

# Alignments downloaded before fetch-data compressed them are still plain on
# disk; compress in place. xz keeps the modification time, so this does not
# trigger retraining, and `train` reads .fasta.xz directly.
for f in data/training_sequences/*.fasta; do
    [ -e "$f" ] || continue
    [ -e "$f.xz" ] && { echo "both $f and $f.xz exist; leaving them for you to resolve"; continue; }
    echo "compressing $f"
    xz -T0 "$f"
done

echo "=== [${st}] seeds ==="
seqdir="data/importations/sequences"
fasta="${seqdir}/${name}_B_1_617_2_seed_sequences.fasta"
manifest="${seqdir}/${name}_B_1_617_2_seed_manifest.csv"
ticks="${seqdir}/${name}_B_1_617_2_ticks.csv"
ran="data/seed_schedule/${st}-2.csv"
imports="${seqdir}/${name}_importation_schedule.csv"
if [ -s "$fasta" ] && [ -s "$manifest" ] && [ -s "$imports" ] \
        && diff -q "$ticks" "$ran" >/dev/null 2>&1; then
    echo "seeds present and matching ${ran}; leaving them untouched"
else
    phylogas prep-seeds --config "$cfg" --seed-mode
    # The seeds must cover exactly the importations EpiHiper made, day for
    # day, or the painter pairs importations with the wrong genomes.
    if ! diff -q "$ticks" "$ran" >/dev/null; then
        echo "seed schedule DIFFERS from ${ran}: some clusters had no retrievable"
        echo "sequence (see the prep-seeds banner above). Not painting."
        exit 1
    fi
    echo "seed schedule matches ${ran}"
fi

echo "=== [${st}] config check ==="
phylogas validate-config --config "$cfg"
