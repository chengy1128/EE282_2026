#!/usr/bin/env bash
set -euo pipefail

mkdir -p homework3/results
FASTA="homework3/data/fasta/dmel-all-chromosome-r6.66.fasta.gz"

# seqkit replaces faSize: num_seqs and sum_len come from seqkit stats.
seqkit stats -T "$FASTA" > homework3/results/genome_stats.tsv

NUM_SEQS=$(awk -F '\t' 'NR==2 {print $4}' homework3/results/genome_stats.tsv)
TOTAL_NT=$(awk -F '\t' 'NR==2 {print $5}' homework3/results/genome_stats.tsv)
TOTAL_NS=$(seqkit seq -s "$FASTA" | awk '{n+=gsub(/[Nn]/,"&")} END{print n+0}')

cat > homework3/results/genome_summary.tsv <<EOT
metric	value
Total number of nucleotides	${TOTAL_NT}
Total number of Ns	${TOTAL_NS}
Total number of sequences	${NUM_SEQS}
EOT
