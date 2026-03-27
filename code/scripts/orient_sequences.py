#!/usr/bin/env python3
"""
Orient all subtelomere sequences in telomere→centromere direction.
- p arm: keep as is (already telomere→centromere)
- q arm: reverse (convert from centromere→telomere to telomere→centromere)

Note: Using reverse only (not reverse complement) to compare structural
features on the same strand, aligned by distance from telomere.
"""

import argparse

def reverse_sequence(seq):
    """Return reversed sequence (same strand, flipped order)."""
    return seq[::-1]

def parse_fasta(fasta_file):
    """Parse FASTA file and yield (header, sequence) tuples."""
    header = None
    seq_lines = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq_lines)
                header = line[1:]  # Remove '>'
                seq_lines = []
            else:
                seq_lines.append(line)

        if header is not None:
            yield header, ''.join(seq_lines)

def orient_sequences(input_fasta, output_fasta):
    """Orient all sequences to telomere→centromere direction."""

    with open(output_fasta, 'w') as out:
        for header, seq in parse_fasta(input_fasta):
            # Parse header to determine arm type
            # Format: chrN_p_arm|GeneName|range or chrN_q_arm|GeneName|range
            parts = header.split('|')
            seq_id = parts[0]

            if '_q_arm' in seq_id:
                # Reverse q arm sequences (same strand, flipped order)
                seq = reverse_sequence(seq)
                # Update header to indicate orientation
                new_header = f">{header}|reversed|telomere_to_centromere"
            else:
                # p arm - keep as is
                new_header = f">{header}|forward|telomere_to_centromere"

            # Write to output
            out.write(new_header + '\n')
            # Write sequence in lines of 60 characters
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + '\n')

    print(f"Oriented sequences saved to {output_fasta}")

def main():
    parser = argparse.ArgumentParser(
        description='Orient subtelomere sequences to telomere→centromere direction'
    )
    parser.add_argument(
        '-i', '--input',
        default='/data/T2T_CHM13v2.0/arm_sequences.fa',
        help='Input FASTA file with p and q arm sequences'
    )
    parser.add_argument(
        '-o', '--output',
        default='/data/T2T_CHM13v2.0/subtelomere_comparison/arm_sequences_oriented.fa',
        help='Output FASTA file with oriented sequences'
    )

    args = parser.parse_args()
    orient_sequences(args.input, args.output)

if __name__ == '__main__':
    main()
