#!/bin/bash
# Extract first 500kb (p) and last 500kb (q) from each chromosome

FASTA="/data/T2T_CHM13v2.0/chm13v2.0.fa"
FAI="/data/T2T_CHM13v2.0/chm13v2.0.fa.fai"
OUTDIR="/data/T2T_CHM13v2.0/chr_terminal_500kb"
SIZE=500000

# Create output directory
mkdir -p "$OUTDIR"

echo "Extracting terminal 500kb regions..."
echo "Output directory: $OUTDIR"
echo ""

# Read each chromosome from .fai file
while read -r chr length rest; do
    echo "Processing $chr (length: $length bp)"

    # Extract first 500kb (p arm)
    if [ $length -ge $SIZE ]; then
        p_end=$SIZE
    else
        p_end=$length
    fi
    samtools faidx "$FASTA" "${chr}:1-${p_end}" | gzip > "${OUTDIR}/${chr}_p.fa.gz"
    echo "  Saved ${chr}_p.fa.gz (1-${p_end})"

    # Extract last 500kb (q arm)
    if [ $length -ge $SIZE ]; then
        q_start=$((length - SIZE + 1))
    else
        q_start=1
    fi
    samtools faidx "$FASTA" "${chr}:${q_start}-${length}" | gzip > "${OUTDIR}/${chr}_q.fa.gz"
    echo "  Saved ${chr}_q.fa.gz (${q_start}-${length})"

done < "$FAI"

echo ""
echo "Extraction complete!"
echo "Files saved in: $OUTDIR"
