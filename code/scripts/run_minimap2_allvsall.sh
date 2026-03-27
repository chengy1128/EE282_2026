#!/bin/bash
# All-vs-all alignment of oriented subtelomere sequences using minimap2
# -X: skip self and dual mappings (avoid aligning a sequence to itself)
# -c: output CIGAR in PAF
# -N: max number of secondary alignments per query (set high for all-vs-all)

INPUT="/data/T2T_CHM13v2.0/subtelomere_comparison/arm_sequences_oriented_no_rdna.fa"
OUTPUT_DIR="/data/T2T_CHM13v2.0/subtelomere_comparison"

# Run minimap2 all-vs-all alignment
echo "Running minimap2 all-vs-all alignment with 64 threads..."
minimap2 -x asm20 -X -c -N 10000 -t 64 \
    ${INPUT} ${INPUT} \
    > ${OUTPUT_DIR}/subtelomere_allvsall.paf \
    2> ${OUTPUT_DIR}/minimap2.log

echo "Alignment complete: ${OUTPUT_DIR}/subtelomere_allvsall.paf"

# Generate summary statistics
echo ""
echo "=== Alignment Summary ==="
echo "Total alignments: $(wc -l < ${OUTPUT_DIR}/subtelomere_allvsall.paf)"
echo ""
echo "Top 20 alignments by alignment block length:"
sort -k11,11nr ${OUTPUT_DIR}/subtelomere_allvsall.paf | head -20 | \
    awk 'BEGIN{OFS="\t"; print "Query", "Target", "QueryLen", "TargetLen", "AlignLen", "MapQ"}
         {print $1, $6, $2, $7, $11, $12}'
