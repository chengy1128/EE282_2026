#!/usr/bin/env bash
set -euo pipefail

mkdir -p homework3/results
GTF="homework3/data/gtf/dmel-all-r6.66.gtf.gz"

# Feature counts sorted from most common to least common.
gzip -dc "$GTF" \
  | awk -F '\t' '!/^#/ {print $3}' \
  | sort \
  | uniq -c \
  | sort -nr \
  | awk 'BEGIN{print "count\tfeature"} {print $1"\t"$2}' \
  > homework3/results/annotation_feature_counts.tsv

# Gene counts on chromosome arms requested in the homework.
gzip -dc "$GTF" \
  | awk -F '\t' '
    BEGIN {print "arm\tgenes"}
    !/^#/ && $3=="gene" {if ($1=="X"||$1=="Y"||$1=="2L"||$1=="2R"||$1=="3L"||$1=="3R"||$1=="4") c[$1]++}
    END {
      print "X\t"  c["X"]+0
      print "Y\t"  c["Y"]+0
      print "2L\t" c["2L"]+0
      print "2R\t" c["2R"]+0
      print "3L\t" c["3L"]+0
      print "3R\t" c["3R"]+0
      print "4\t"  c["4"]+0
    }
  ' \
  > homework3/results/annotation_genes_per_arm.tsv
