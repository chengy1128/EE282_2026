#!/usr/bin/env python3
"""
Generate a Circos-style plot showing individual pairwise alignments.
Each alignment is drawn as a separate link between the specific regions.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import re

try:
    from pycirclize import Circos
    HAS_PYCIRCLIZE = True
except ImportError:
    HAS_PYCIRCLIZE = False
    print("Warning: pycirclize not found. Install with: pip install pycirclize")

def parse_paf(paf_file, min_align_len=5000):
    """Parse PAF file and extract individual alignment information."""
    alignments = []
    with open(paf_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 12:
                query = fields[0]
                query_len = int(fields[1])
                query_start = int(fields[2])
                query_end = int(fields[3])
                target = fields[5]
                target_len = int(fields[6])
                target_start = int(fields[7])
                target_end = int(fields[8])
                matches = int(fields[9])
                align_len = int(fields[10])

                # Filter by minimum alignment length
                if align_len < min_align_len:
                    continue

                # Skip self-alignments (same sequence)
                q_simple = simplify_name(query)
                t_simple = simplify_name(target)
                if q_simple == t_simple:
                    continue

                alignments.append({
                    'query': query,
                    'query_simple': q_simple,
                    'query_len': query_len,
                    'query_start': query_start,
                    'query_end': query_end,
                    'target': target,
                    'target_simple': t_simple,
                    'target_len': target_len,
                    'target_start': target_start,
                    'target_end': target_end,
                    'matches': matches,
                    'align_len': align_len,
                    'identity': matches / align_len if align_len > 0 else 0
                })
    return alignments

def simplify_name(full_name):
    """Extract simple chromosome arm name from full header."""
    match = re.match(r'(chr[0-9XY]+_[pq]_arm)', full_name)
    if match:
        return match.group(1)
    return full_name

def get_seq_info(alignments):
    """Get sequence names and lengths."""
    seq_lengths = {}
    for aln in alignments:
        seq_lengths[aln['query_simple']] = aln['query_len']
        seq_lengths[aln['target_simple']] = aln['target_len']

    # Sort sequences by chromosome number then arm
    sequences = sorted(seq_lengths.keys(), key=lambda x: (
        int(re.search(r'chr(\d+)', x).group(1)) if re.search(r'chr(\d+)', x) else
        (23 if 'chrX' in x else 24 if 'chrY' in x else 25),
        0 if '_p_' in x else 1
    ))

    return sequences, seq_lengths

def get_identity_color(identity, alpha=0.6):
    """Get color based on sequence identity."""
    if identity >= 0.95:
        return (0.8, 0.2, 0.2, alpha)  # Red
    elif identity >= 0.85:
        return (0.8, 0.5, 0.2, alpha)  # Orange
    elif identity >= 0.75:
        return (0.2, 0.6, 0.8, alpha)  # Blue
    else:
        return (0.5, 0.5, 0.5, alpha)  # Gray

def plot_circos_individual(sequences, seq_lengths, alignments, output_file):
    """Generate Circos plot with individual alignments using pycirclize."""

    # Prepare sector data
    sectors = {seq: seq_lengths[seq] for seq in sequences}

    # Create Circos instance
    circos = Circos(sectors, space=3)

    # Color palette for chromosomes
    colors = plt.cm.tab20(np.linspace(0, 1, 24))
    chr_colors = {}
    for seq in sequences:
        chr_num = re.search(r'chr(\d+|X|Y)', seq).group(1)
        if chr_num == 'X':
            idx = 22
        elif chr_num == 'Y':
            idx = 23
        else:
            idx = int(chr_num) - 1
        chr_colors[seq] = colors[idx % 24]

    # Add tracks to each sector
    for sector in circos.sectors:
        # Outer track with sector name
        track = sector.add_track((95, 100))
        track.axis(fc=chr_colors[sector.name], ec="black", lw=0.5)

        # Add sector label
        label = sector.name.replace('chr', '').replace('_arm', '').replace('_', '')
        track.text(label, size=30, r=115)

        # Inner track showing scale/position
        track2 = sector.add_track((88, 94))
        track2.axis(fc="lightgray", ec="none")

        # Add tick marks every 50kb
        major_ticks = list(range(0, seq_lengths[sector.name], 50000))
        track2.xticks(major_ticks, label_size=0)

    # Sort alignments by length (draw smaller ones last so they're on top)
    alignments_sorted = sorted(alignments, key=lambda x: x['align_len'], reverse=True)

    # Draw individual alignment links
    for aln in alignments_sorted:
        q_name = aln['query_simple']
        t_name = aln['target_simple']

        # Get color based on identity
        color = get_identity_color(aln['identity'], alpha=0.5)

        # Line width based on alignment length (log scale)
        lw = max(0.3, np.log10(aln['align_len']) - 2.5)

        # Add link between specific regions
        circos.link(
            (q_name, aln['query_start'], aln['query_end']),
            (t_name, aln['target_start'], aln['target_end']),
            color=color,
            lw=lw
        )

    # Plot and save
    fig = circos.plotfig(figsize=(30, 30))

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=(0.8, 0.2, 0.2, 0.7), label='Identity ≥95%'),
        Patch(facecolor=(0.8, 0.5, 0.2, 0.6), label='Identity 85-95%'),
        Patch(facecolor=(0.2, 0.6, 0.8, 0.5), label='Identity 75-85%'),
        Patch(facecolor=(0.5, 0.5, 0.5, 0.4), label='Identity <75%'),
    ]
    fig.legend(handles=legend_elements, loc='lower right', fontsize=20)

    plt.title('Subtelomere Individual Pairwise Alignments', fontsize=40, y=1.1)
    fig.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Circos plot saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description='Generate Circos-style plot showing individual pairwise alignments'
    )
    parser.add_argument(
        '-i', '--input',
        default='/data/T2T_CHM13v2.0/subtelomere_comparison/subtelomere_allvsall.paf',
        help='Input PAF file from minimap2'
    )
    parser.add_argument(
        '-o', '--output',
        default='/data/T2T_CHM13v2.0/subtelomere_comparison/subtelomere_circos_individual.png',
        help='Output PNG file'
    )
    parser.add_argument(
        '-m', '--min-align',
        type=int,
        default=1000,
        help='Minimum alignment length to display (default: 5000 bp)'
    )

    args = parser.parse_args()

    if not HAS_PYCIRCLIZE:
        print("Error: pycirclize is required for this script.")
        print("Install with: pip install pycirclize")
        return

    print("Parsing PAF file...")
    alignments = parse_paf(args.input, min_align_len=args.min_align)
    print(f"Loaded {len(alignments)} alignments (≥{args.min_align} bp, excluding self)")

    print("Getting sequence information...")
    sequences, seq_lengths = get_seq_info(alignments)
    print(f"Found {len(sequences)} sequences")

    print("Generating Circos plot with individual alignments...")
    plot_circos_individual(sequences, seq_lengths, alignments, args.output)

if __name__ == '__main__':
    main()
