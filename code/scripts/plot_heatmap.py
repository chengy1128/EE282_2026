#!/usr/bin/env python3
"""
Generate a heatmap of pairwise subtelomere similarity from minimap2 PAF output.
Clusters sequences by similarity to reveal subtelomere families.
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from collections import defaultdict
import re

def parse_paf(paf_file):
    """Parse PAF file and extract alignment information."""
    alignments = []
    with open(paf_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 12:
                query = fields[0]
                query_len = int(fields[1])
                target = fields[5]
                target_len = int(fields[6])
                matches = int(fields[9])
                align_len = int(fields[10])

                alignments.append({
                    'query': query,
                    'query_len': query_len,
                    'target': target,
                    'target_len': target_len,
                    'matches': matches,
                    'align_len': align_len,
                    'identity': matches / align_len if align_len > 0 else 0
                })
    return alignments

def simplify_name(full_name):
    """Extract simple chromosome arm name from full header."""
    # Extract chrN_p_arm or chrN_q_arm from the full name
    match = re.match(r'(chr[0-9XY]+_[pq]_arm)', full_name)
    if match:
        return match.group(1)
    return full_name

def build_similarity_matrix(alignments, metric='coverage'):
    """
    Build pairwise similarity matrix.

    metric options:
    - 'coverage': total aligned bases / min(seq1_len, seq2_len)
    - 'define_coverage': total matches / min(seq1_len, seq2_len)
    - 'coverage_long': total matches / max(seq1_len, seq2_len)
    - 'identity': weighted average identity
    - 'align_len': total alignment length
    """
    # Get unique sequences and their lengths
    seq_lengths = {}
    for aln in alignments:
        q_simple = simplify_name(aln['query'])
        t_simple = simplify_name(aln['target'])
        seq_lengths[q_simple] = aln['query_len']
        seq_lengths[t_simple] = aln['target_len']

    sequences = sorted(seq_lengths.keys(), key=lambda x: (
        # Sort by chromosome number, then by arm
        int(re.search(r'chr(\d+)', x).group(1)) if re.search(r'chr(\d+)', x) else
        (23 if 'chrX' in x else 24 if 'chrY' in x else 25),
        x.split('_')[1]  # p before q
    ))

    # Initialize matrix
    n = len(sequences)
    matrix = np.zeros((n, n))
    seq_to_idx = {seq: i for i, seq in enumerate(sequences)}

    # Aggregate alignments
    pair_data = defaultdict(lambda: {'matches': 0, 'align_len': 0})

    for aln in alignments:
        q_simple = simplify_name(aln['query'])
        t_simple = simplify_name(aln['target'])

        # Skip self-alignments for off-diagonal
        if q_simple == t_simple:
            continue

        key = tuple(sorted([q_simple, t_simple]))
        pair_data[key]['matches'] += aln['matches']
        pair_data[key]['align_len'] += aln['align_len']

    # Fill matrix
    for (seq1, seq2), data in pair_data.items():
        i, j = seq_to_idx[seq1], seq_to_idx[seq2]

        if metric == 'coverage':
            min_len = min(seq_lengths[seq1], seq_lengths[seq2])
            value = data['align_len'] / min_len if min_len > 0 else 0
        elif metric == 'define_coverage':
            min_len = min(seq_lengths[seq1], seq_lengths[seq2])
            value = data['matches'] / min_len if min_len > 0 else 0
        elif metric == 'coverage_long':
            max_len = max(seq_lengths[seq1], seq_lengths[seq2])
            value = data['matches'] / max_len if max_len > 0 else 0
        elif metric == 'identity':
            value = data['matches'] / data['align_len'] if data['align_len'] > 0 else 0
        else:  # align_len
            value = data['align_len']

        matrix[i, j] = value
        matrix[j, i] = value

    # Fill diagonal with 1 for coverage/identity, or seq_length for align_len
    for i, seq in enumerate(sequences):
        if metric in ['coverage', 'identity', 'define_coverage', 'coverage_long']:
            matrix[i, i] = 1.0
        else:
            matrix[i, i] = seq_lengths[seq]

    return matrix, sequences

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    """Truncate a colormap to use a subset of the range."""
    new_cmap = mcolors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def plot_heatmap(matrix, sequences, output_file, metric='coverage', cluster=True):
    """Generate clustered heatmap."""

    # Create DataFrame
    df = pd.DataFrame(matrix, index=sequences, columns=sequences)

    # Simplify labels for display
    labels = [s.replace('_arm', '').replace('chr', '') for s in sequences]
    df.index = labels
    df.columns = labels

    # Set up figure
    plt.figure(figsize=(30, 30))

    # Create clustermap
    if metric == 'coverage':
        vmax = min(1.0, np.percentile(matrix[matrix < 1], 99)) if np.any(matrix < 1) else 1.0
        vmin = 0.0
        # Use YlOrRd but stop at 0.8 (remove darkest red)
        cmap = truncate_colormap(plt.get_cmap('YlOrRd'), 0.0, 0.75)
        cbar_label = 'Alignment Coverage'
    elif metric == 'define_coverage':
        # Calculate stats for diagnostics
        flat_data = matrix.flatten()
        print(f"Metrics Stats - Min: {flat_data.min():.4f}, Max: {flat_data.max():.4f}")
        
        # Determine logical Vmax
        relevant_data = matrix[matrix < 1.0]
        if relevant_data.size > 0:
             calc_vmax = np.percentile(relevant_data, 99)
             vmax = min(1.0, calc_vmax)
             # Ensure a minimum range so we don't collapse to 0 if data is all 0
             vmax = max(vmax, 0.01)
        else:
             # Default to 1.0 if no small values found (only 0 or >1)
             vmax = 1.0
             
        vmin = 0.0
        print(f"Color Map Range: 0.0 (White) -> {vmax:.4f} (Dark Red). Values > {vmax:.4f} are saturated.")
        
        # Create custom colormap starting with white
        # Use a list with from_list to ensure white occupies the start
        # We sample YlOrRd from 0.0 (lightest) to 0.75 to get a smooth transition
        base_colors = plt.get_cmap('YlOrRd')(np.linspace(0.0, 0.75, 256))
        
        # Force the very first color to be pure white (1,1,1,1)
        # YlOrRd(0.0) is usually off-white/pale-yellow, so we override it
        base_colors[0] = [1.0, 1.0, 1.0, 1.0]
        
        cmap = mcolors.LinearSegmentedColormap.from_list('white_ylorrd', base_colors)
        
        cbar_label = f'Defined Coverage (Matches / ShorterSeq Len)\n[0.0=White, >{vmax:.2f}=Saturated]'
    elif metric == 'coverage_long':
        # Calculate stats for diagnostics
        flat_data = matrix.flatten()
        print(f"Metrics Stats - Min: {flat_data.min():.4f}, Max: {flat_data.max():.4f}")
        
        # Determine logical Vmax
        relevant_data = matrix[matrix < 1.0]
        if relevant_data.size > 0:
             calc_vmax = np.percentile(relevant_data, 99)
             vmax = min(1.0, calc_vmax)
             vmax = max(vmax, 0.01)
        else:
             vmax = 1.0
             
        vmin = 0.0
        print(f"Color Map Range: 0.0 (White) -> {vmax:.4f} (Dark Red). Values > {vmax:.4f} are saturated.")
        
        # Create custom colormap starting with white
        base_colors = plt.get_cmap('YlOrRd')(np.linspace(0.0, 0.75, 256))
        base_colors[0] = [1.0, 1.0, 1.0, 1.0]
        cmap = mcolors.LinearSegmentedColormap.from_list('white_ylorrd', base_colors)
        
        cbar_label = f'Defined Coverage (Matches / LongerSeq Len)\n[0.0=White, >{vmax:.2f}=Saturated]'
    elif metric == 'identity':
        vmax = 1.0
        vmin = 0.6
        cmap = mcolors.LinearSegmentedColormap.from_list("white_pink_red", ["white", "pink", "red"])
        cbar_label = 'Sequence Identity'
    else:
        vmax = None
        vmin = 0.0
        cmap = 'viridis'
        cbar_label = 'Alignment Length (bp)'

    # Create clustered heatmap
    g = sns.clustermap(
        df,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        figsize=(30, 30),
        dendrogram_ratio=(0.1, 0.1),
        cbar_pos=(-0.05, 0.8, 0.03, 0.15),
        xticklabels=True,
        yticklabels=True,
        linewidths=0.5,
        linecolor='white',
        row_cluster=cluster,
        col_cluster=cluster
    )

    # Adjust labels
    g.ax_heatmap.set_xlabel('Subtelomere', fontsize=28)
    g.ax_heatmap.set_ylabel('Subtelomere', fontsize=28)
    g.ax_heatmap.tick_params(axis='both', labelsize=28)

    # Rotate x labels
    plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)

    # Add title
    g.fig.suptitle(f'Subtelomere Pairwise Similarity ({cbar_label})',
                   fontsize=50, y=1.02)

    # Save
    g.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Heatmap saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description='Generate heatmap of subtelomere pairwise similarity'
    )
    parser.add_argument(
        '-i', '--input',
        default='/data/T2T_CHM13v2.0/subtelomere_comparison/subtelomere_allvsall.paf',
        help='Input PAF file from minimap2'
    )
    parser.add_argument(
        '-o', '--output',
        default='/data/T2T_CHM13v2.0/subtelomere_comparison/subtelomere_heatmap.png',
        help='Output PNG file'
    )
    parser.add_argument(
        '-m', '--metric',
        choices=['coverage', 'identity', 'align_len', 'define_coverage', 'coverage_long'],
        default='coverage',
        help='Similarity metric (default: coverage)'
    )
    parser.add_argument(
        '--no-cluster',
        action='store_false',
        dest='cluster',
        help='Disable clustering to sort by chromosome number'
    )
    parser.set_defaults(cluster=True)

    parser.add_argument(
        '--save-matrix',
        help='Path to save the similarity matrix as CSV'
    )

    args = parser.parse_args()

    print("Parsing PAF file...")
    alignments = parse_paf(args.input)
    print(f"Loaded {len(alignments)} alignments")

    print(f"Building similarity matrix (metric: {args.metric})...")
    matrix, sequences = build_similarity_matrix(alignments, args.metric)
    print(f"Matrix size: {len(sequences)}x{len(sequences)}")

    if args.save_matrix:
        print(f"Saving matrix to {args.save_matrix}...")
        df_out = pd.DataFrame(matrix, index=sequences, columns=sequences)
        df_out.to_csv(args.save_matrix)

    print(f"Generating heatmap (clustering: {args.cluster})...")
    plot_heatmap(matrix, sequences, args.output, args.metric, args.cluster)

if __name__ == '__main__':
    main()
