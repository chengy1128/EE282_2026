#!/usr/bin/env bash
set -euo pipefail

# Homework 4 Part 3: Genome assembly with hifiasm (run on HPC3)
# Usage: srun -A class_ee282 --cpus-per-task=16 --mem=32G --time=2:00:00 --pty bash hw4_assembly.sh

WORKDIR="$HOME/ee282/homework4"
DATADIR="${WORKDIR}/data"
ASMDIR="${WORKDIR}/assembly"
mkdir -p "$DATADIR" "$ASMDIR"

# Step 1: Use HiFi reads directly from shared class directory
READS="${DATADIR}/ISO_HiFi_Shukla2025.fasta.gz"
#"/pub/jje/ee282/ISO_HiFi_Shukla2025.fasta.gz"

# Step 2: Run hifiasm (homozygous/inbred mode with -l0)
THREADS="${1:-$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 16)}"
echo "=== Running hifiasm with ${THREADS} threads ==="
hifiasm -o "${ASMDIR}/dmel_hifi" -t "$THREADS" -l0 "$READS"

# Step 3: Convert primary contigs GFA to FASTA
echo "=== Converting GFA to FASTA ==="
GFA="${ASMDIR}/dmel_hifi.bp.p_ctg.gfa"
FASTA="${ASMDIR}/dmel_hifi.bp.p_ctg.fa"
awk '/^S/{print ">"$2; print $3}' "$GFA" > "$FASTA"

echo "=== Assembly complete ==="
echo "Primary contigs FASTA: $FASTA"
echo ""
echo "Quick stats:"
bioawk -c fastx '{sum+=length($seq); n++} END{print "Sequences:", n; print "Total bp:", sum}' "$FASTA"
