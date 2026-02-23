#!/usr/bin/env bash
set -euo pipefail

mkdir -p homework3/results

# Verify genome FASTA md5 using just the target line from FlyBase md5sum.txt.
grep ' dmel-all-chromosome-r6.66.fasta.gz$' homework3/data/fasta/md5sum.txt > homework3/results/genome.md5
(
  cd homework3/data/fasta
  md5sum -c ../../results/genome.md5
) | tee homework3/results/genome_integrity.txt

# Verify annotation GTF md5 from FlyBase md5sum.txt.
(
  cd homework3/data/gtf
  md5sum -c md5sum.txt
) | tee homework3/results/annotation_integrity.txt
