#!/usr/bin/env python3
"""
Extract sequences from chromosome p and q arms based on first protein coding gene list.
Excludes telomere sequences:
- p arm: extract from (telomere_length + 1) to (telomere_length + distance_from_edge)
- q arm: extract from (chrom_len - telomere_length - distance_from_edge + 1) to (chrom_len - telomere_length)
"""

import argparse
from pysam import FastaFile

def parse_gene_list(gene_list_file):
    """Parse the first protein coding gene list file with telomere length."""
    entries = []
    with open(gene_list_file, 'r') as f:
        header = f.readline()  # Skip header
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            # New format: Chromosome, Chr_Length, arm, Gene_Name, Start, End, telomere_length, Distance_from_Edge
            if len(parts) >= 8:
                chrom = parts[0]
                chrom_len = int(parts[1])
                arm = parts[2]
                gene_name = parts[3]
                telomere_length = int(parts[6])
                distance_from_edge = int(parts[7])
                entries.append({
                    'chrom': chrom,
                    'chrom_len': chrom_len,
                    'arm': arm,
                    'gene_name': gene_name,
                    'telomere_length': telomere_length,
                    'distance_from_edge': distance_from_edge
                })
    return entries

def extract_sequences(fasta_file, gene_list_file, output_file):
    """Extract p and q arm sequences based on gene list, excluding telomeres."""

    # Parse gene list
    entries = parse_gene_list(gene_list_file)

    # Open FASTA file
    fasta = FastaFile(fasta_file)

    # Get chromosome lengths from FASTA for validation
    fasta_chrom_lengths = dict(zip(fasta.references, fasta.lengths))

    with open(output_file, 'w') as out:
        for entry in entries:
            chrom = entry['chrom']
            arm = entry['arm']
            gene_name = entry['gene_name']
            distance = entry['distance_from_edge']
            telomere_len = entry['telomere_length']
            chrom_len = entry['chrom_len']

            if chrom not in fasta_chrom_lengths:
                print(f"Warning: {chrom} not found in FASTA, skipping")
                continue

            if arm == 'p':
                # p arm: exclude telomere at beginning
                # Extract from (telomere_length) to (telomere_length + distance_from_edge)
                # pysam uses 0-based coordinates
                start = telomere_len  # 0-based start (after telomere)
                end = telomere_len + distance  # 0-based exclusive end

                if end > chrom_len:
                    print(f"Warning: {chrom} p arm end ({end}) exceeds chromosome length ({chrom_len}), using chromosome length")
                    end = chrom_len

                # Extract sequence
                seq = fasta.fetch(chrom, start, end)

                # 1-based coordinates for header (exact coordinates from reference)
                header = f">{chrom}_p_arm|{gene_name}|{start + 1}-{end}"

            elif arm == 'q':
                # q arm: exclude telomere at end
                # Extract from (chrom_len - telomere_len - distance) to (chrom_len - telomere_len)
                # pysam uses 0-based coordinates
                start = chrom_len - telomere_len - distance  # 0-based start
                end = chrom_len - telomere_len  # 0-based exclusive end (before telomere)

                if start < 0:
                    print(f"Warning: {chrom} q arm start ({start}) is negative, using position 0")
                    start = 0

                # Extract sequence
                seq = fasta.fetch(chrom, start, end)

                # 1-based coordinates for header (exact coordinates from reference)
                header = f">{chrom}_q_arm|{gene_name}|{start + 1}-{end}"

            else:
                print(f"Warning: Unknown arm type '{arm}' for {chrom}, skipping")
                continue

            # Write to output
            out.write(header + '\n')
            # Write sequence in lines of 60 characters
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + '\n')

    fasta.close()
    print(f"Sequences extracted to {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description='Extract p and q arm sequences based on first protein coding gene list'
    )
    parser.add_argument(
        '-f', '--fasta',
        default='/data/T2T_CHM13v2.0/chm13v2.0.fa',
        help='Reference FASTA file (default: /data/T2T_CHM13v2.0/chm13v2.0.fa)'
    )
    parser.add_argument(
        '-g', '--gene-list',
        default='/data/T2T_CHM13v2.0/first_protein_coding_gene_list_copy.txt',
        help='First protein coding gene list file'
    )
    parser.add_argument(
        '-o', '--output',
        default='/data/T2T_CHM13v2.0/arm_sequences.fa',
        help='Output FASTA file (default: /data/T2T_CHM13v2.0/arm_sequences.fa)'
    )

    args = parser.parse_args()

    extract_sequences(args.fasta, args.gene_list, args.output)

if __name__ == '__main__':
    main()
