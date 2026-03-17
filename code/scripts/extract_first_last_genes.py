#!/usr/bin/env python3
"""
Extract the first and last protein coding gene from each chromosome.
Only include gene features (no transcripts, exons, CDS).
Calculate distance from chromosome start and end using .fai file.
"""
import re
from collections import defaultdict

input_file = '/data/T2T_CHM13v2.0/chm13v2.0_RefSeq_Liftoff_v5.2.gff3'
fai_file = '/data/T2T_CHM13v2.0/chm13v2.0.fa.fai'
output_file = '/data/T2T_CHM13v2.0/chm13v2.0_RefSeq_Liftoff_v5.2_protein_coding_first_last.gff3'
summary_file = '/data/T2T_CHM13v2.0/chm13v2.0_RefSeq_Liftoff_v5.2_protein_coding_first_last_summary.txt'

print("Reading chromosome lengths from .fai file...")

# Read chromosome lengths from .fai file
chr_lengths = {}
with open(fai_file, 'r') as f:
    for line in f:
        fields = line.strip().split('\t')
        if len(fields) >= 2:
            chr_name = fields[0]
            chr_length = int(fields[1])
            chr_lengths[chr_name] = chr_length

print(f"Found {len(chr_lengths)} chromosomes in .fai file")

print("\nReading protein coding genes...")

# Dictionary to store genes by chromosome
# Structure: {chr: [(start, end, gene_id, gene_name, line), ...]}
genes_by_chr = defaultdict(list)

# Read all protein coding genes
header_lines = []
with open(input_file, 'r') as f:
    for line in f:
        if line.startswith('#'):
            header_lines.append(line)
            continue

        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue

        chrom = fields[0]
        feature_type = fields[2]
        start = int(fields[3])
        end = int(fields[4])
        attributes = fields[8]

        # Only process gene features with protein_coding biotype
        if feature_type != 'gene':
            continue
        if 'gene_biotype=protein_coding' not in attributes:
            continue

        # Extract gene ID and name
        id_match = re.search(r'ID=([^;]+)', attributes)
        name_match = re.search(r'gene_name=([^;]+)', attributes)

        if not id_match:
            continue

        gene_id = id_match.group(1)
        gene_name = name_match.group(1) if name_match else gene_id

        # Store gene information
        genes_by_chr[chrom].append((start, end, gene_id, gene_name, line))

print(f"Found protein coding genes on {len(genes_by_chr)} chromosomes")

# Find first and last gene for each chromosome
selected_gene_lines = []
summary_data = []

for chrom in sorted(genes_by_chr.keys(), key=lambda x: (x.replace('chr', '').replace('X', '23').replace('Y', '24').replace('M', '25').zfill(2))):
    genes = genes_by_chr[chrom]
    if not genes:
        continue

    # Get chromosome length from .fai file
    chr_length = chr_lengths.get(chrom, 0)
    if chr_length == 0:
        print(f"Warning: {chrom} not found in .fai file, skipping...")
        continue

    # Sort by start position
    genes.sort(key=lambda x: x[0])

    # First gene (smallest start position)
    first_gene = genes[0]
    first_start, first_end, first_id, first_name, first_line = first_gene
    dist_from_start = first_start - 1  # Distance from position 1

    selected_gene_lines.append(first_line)

    # Last gene (largest end position)
    genes_sorted_by_end = sorted(genes, key=lambda x: x[1], reverse=True)
    last_gene = genes_sorted_by_end[0]
    last_start, last_end, last_id, last_name, last_line = last_gene
    dist_from_end = chr_length - last_end

    # Make sure first and last are different genes
    if first_id != last_id:
        selected_gene_lines.append(last_line)

    print(f"\n{chrom}: Length={chr_length:,} bp")
    print(f"  First: {first_name} ({first_id})")
    print(f"         Position: {first_start:,}-{first_end:,}")
    print(f"         Distance from chr start: {dist_from_start:,} bp")
    print(f"  Last:  {last_name} ({last_id})")
    print(f"         Position: {last_start:,}-{last_end:,}")
    print(f"         Distance from chr end: {dist_from_end:,} bp")

    # Store summary data - separate rows for first and last genes
    summary_data.append({
        'chr': chrom,
        'chr_length': chr_length,
        'position': 'first',
        'gene_name': first_name,
        'gene_id': first_id,
        'start': first_start,
        'end': first_end,
        'distance_from_edge': dist_from_start
    })

    summary_data.append({
        'chr': chrom,
        'chr_length': chr_length,
        'position': 'last',
        'gene_name': last_name,
        'gene_id': last_id,
        'start': last_start,
        'end': last_end,
        'distance_from_edge': dist_from_end
    })

print(f"\n\nTotal selected genes: {len(selected_gene_lines)}")

# Write summary file
print(f"\nWriting summary to: {summary_file}")
with open(summary_file, 'w') as f:
    # Write header
    f.write("Chromosome\tChr_Length\tPosition\tGene_Name\tGene_ID\tStart\tEnd\tDistance_from_Edge\n")

    # Write data - each gene on its own row
    for row in summary_data:
        f.write(f"{row['chr']}\t{row['chr_length']}\t{row['position']}\t{row['gene_name']}\t")
        f.write(f"{row['gene_id']}\t{row['start']}\t{row['end']}\t{row['distance_from_edge']}\n")

# Write GFF3 output file
print(f"Writing GFF3 output to: {output_file}")
with open(output_file, 'w') as f_out:
    # Write header
    for header in header_lines:
        f_out.write(header)

    # Write selected gene lines
    for line in selected_gene_lines:
        f_out.write(line)

print(f"\nExtraction complete!")
print(f"GFF3 file: {output_file}")
print(f"Summary file: {summary_file}")
