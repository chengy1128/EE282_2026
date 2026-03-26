#!/usr/bin/env bash
set -euo pipefail

# Homework 4 Part 1: Summarize partitions of Drosophila melanogaster genome
# Partition: sequences <= 100kb and sequences > 100kb

FASTA="homework3/data/fasta/dmel-all-chromosome-r6.66.fasta.gz"
OUTDIR="homework4/results"
mkdir -p "$OUTDIR"

THRESHOLD=100000

echo "=== Partitioning genome at ${THRESHOLD} bp threshold ==="

# Use bioawk to compute per-sequence: name, length, Ns, and partition label
# Then aggregate with awk
bioawk -c fastx '{
    len = length($seq);
    ns  = gsub(/[Nn]/, "", $seq);
    if (len <= 100000)
        print "le100kb", len, ns;
    else
        print "gt100kb", len, ns;
}' <(gunzip -c "$FASTA") \
| awk '
BEGIN { OFS="\t" }
{
    partition = $1; len = $2; ns = $3;
    total_nt[partition] += len;
    total_ns[partition] += ns;
    total_seq[partition] += 1;
}
END {
    print "partition", "total_nucleotides", "total_Ns", "total_sequences";
    for (p in total_nt)
        print p, total_nt[p], total_ns[p], total_seq[p];
}' > "${OUTDIR}/genome_partition_summary.tsv"

echo "Results written to ${OUTDIR}/genome_partition_summary.tsv"
cat "${OUTDIR}/genome_partition_summary.tsv" | column -t
