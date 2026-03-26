#!/usr/bin/env bash
set -euo pipefail

# FlyBase release used for this homework (current at time of completion).
RELEASE="dmel_r6.66_FB2025_05"
BASE_URL="https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/${RELEASE}"

mkdir -p homework3/data/fasta homework3/data/gtf

# Genome assembly FASTA (all chromosomes) + checksum file.
curl -L --fail --retry 3 -o homework3/data/fasta/dmel-all-chromosome-r6.66.fasta.gz \
  "${BASE_URL}/fasta/dmel-all-chromosome-r6.66.fasta.gz"
curl -L --fail --retry 3 -o homework3/data/fasta/md5sum.txt \
  "${BASE_URL}/fasta/md5sum.txt"

# Annotation GTF + checksum file.
curl -L --fail --retry 3 -o homework3/data/gtf/dmel-all-r6.66.gtf.gz \
  "${BASE_URL}/gtf/dmel-all-r6.66.gtf.gz"
curl -L --fail --retry 3 -o homework3/data/gtf/md5sum.txt \
  "${BASE_URL}/gtf/md5sum.txt"
