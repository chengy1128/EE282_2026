#!/usr/bin/env python3
"""
Plot distance from chromosome edge with broken Y-axis (3 sections)
"""

import os
os.environ['MPLBACKEND'] = 'Agg'

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the data
df = pd.read_csv('/Users/yunsawa/Desktop/T2T/CHM13/first_protein_coding_gene_list.txt', sep='\t')

# Clean the distance column (remove commas and convert to int)
df['Distance_from_Edge(no telomere)'] = df['Distance_from_Edge(no telomere)'].astype(str).str.replace(',', '').astype(int)

# Sort by chromosome number and arm
def chr_sort_key(chr_name):
    chr_part = chr_name.replace('chr', '')
    if chr_part == 'X':
        return 23
    elif chr_part == 'Y':
        return 24
    else:
        return int(chr_part)

df['chr_order'] = df['Chromosome'].apply(chr_sort_key)
df = df.sort_values(['chr_order', 'arm'])

# Get unique chromosomes in order
chromosomes = df['Chromosome'].drop_duplicates().tolist()

# Prepare data for grouped bar chart (p and q arms side by side)
p_values = []
q_values = []
x_labels = []

for chrom in chromosomes:
    chrom_data = df[df['Chromosome'] == chrom]
    p_data = chrom_data[chrom_data['arm'] == 'p']['Distance_from_Edge(no telomere)']
    q_data = chrom_data[chrom_data['arm'] == 'q']['Distance_from_Edge(no telomere)']

    # Convert to kb, use 0 if no data
    p_values.append(p_data.values[0] / 1000 if len(p_data) > 0 else 0)
    q_values.append(q_data.values[0] / 1000 if len(q_data) > 0 else 0)
    x_labels.append(chrom)

x_pos = np.arange(len(chromosomes))
bar_width = 0.35

# Define break points (in kb)
# Bottom section: 0-50 kb
# Middle section: 50-1850 kb
# Top section: 5000-6000 kb
bottom_upper = 50
middle_lower = 50
middle_upper = 1850
top_lower = 5000

# Create figure with three subplots (broken axis)
fig, (ax_top, ax_mid, ax_bottom) = plt.subplots(3, 1, sharex=True, figsize=(36, 24),
                                                 gridspec_kw={'height_ratios': [1, 3, 1.5], 'hspace': 0.05})

# Plot grouped bars on all three axes (p arm in blue, q arm in orange)
# p arm bars (left side of each group)
ax_top.bar(x_pos - bar_width/2, p_values, bar_width, color='#1f77b4', edgecolor='black', linewidth=0.5, label='p arm')
ax_mid.bar(x_pos - bar_width/2, p_values, bar_width, color='#1f77b4', edgecolor='black', linewidth=0.5)
ax_bottom.bar(x_pos - bar_width/2, p_values, bar_width, color='#1f77b4', edgecolor='black', linewidth=0.5)

# q arm bars (right side of each group)
ax_top.bar(x_pos + bar_width/2, q_values, bar_width, color='#ff7f0e', edgecolor='black', linewidth=0.5, label='q arm')
ax_mid.bar(x_pos + bar_width/2, q_values, bar_width, color='#ff7f0e', edgecolor='black', linewidth=0.5)
ax_bottom.bar(x_pos + bar_width/2, q_values, bar_width, color='#ff7f0e', edgecolor='black', linewidth=0.5)

# Combined y_values for setting axis limits
all_values = p_values + q_values

# Set axis limits (with extra room for labels)
ax_top.set_ylim(top_lower, max(all_values) * 1.15)
ax_mid.set_ylim(middle_lower, middle_upper + 300)
ax_bottom.set_ylim(0, bottom_upper + 15)

# Hide the spines between axes
ax_top.spines['bottom'].set_visible(False)
ax_mid.spines['top'].set_visible(False)
ax_mid.spines['bottom'].set_visible(False)
ax_bottom.spines['top'].set_visible(False)

ax_top.xaxis.tick_top()
ax_top.tick_params(labeltop=False)
ax_mid.tick_params(labeltop=False)

# Add break marks (using matplotlib's recommended approach)
d = 0.5  # proportion of vertical to horizontal extent of the slanted line
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=20,
              linestyle="none", color='k', mec='k', mew=2, clip_on=False)

# Break between top and middle
ax_top.plot([0, 1], [0, 0], transform=ax_top.transAxes, **kwargs)
ax_mid.plot([0, 1], [1, 1], transform=ax_mid.transAxes, **kwargs)

# Break between middle and bottom
ax_mid.plot([0, 1], [0, 0], transform=ax_mid.transAxes, **kwargs)
ax_bottom.plot([0, 1], [1, 1], transform=ax_bottom.transAxes, **kwargs)

# Add distance labels on top of each bar
def add_bar_label(ax_top, ax_mid, ax_bottom, x, y, top_lower, middle_lower, middle_upper):
    if y == 0:
        return
    label = f'{int(y):,}'
    # Determine which axis to add the label to based on y value
    if y > top_lower:
        # Values in top section
        ax_top.text(x, y + 50, label, ha='center', va='bottom', fontsize=28, rotation=90, color='black')
    elif y > middle_upper:
        # Values in the gap between middle and top - place label at top of middle axis
        ax_mid.text(x, middle_upper + 250, label, ha='center', va='bottom', fontsize=28, rotation=90, color='black')
    elif y > middle_lower:
        ax_mid.text(x, y + 30, label, ha='center', va='bottom', fontsize=28, rotation=90, color='black')
    else:
        ax_bottom.text(x, y + 1, label, ha='center', va='bottom', fontsize=28, rotation=90, color='black')

# Add labels for p arm bars
for i, (x, y) in enumerate(zip(x_pos - bar_width/2, p_values)):
    add_bar_label(ax_top, ax_mid, ax_bottom, x, y, top_lower, middle_lower, middle_upper)

# Add labels for q arm bars
for i, (x, y) in enumerate(zip(x_pos + bar_width/2, q_values)):
    add_bar_label(ax_top, ax_mid, ax_bottom, x, y, top_lower, middle_lower, middle_upper)

# Set x-axis labels
ax_bottom.set_xticks(x_pos)
ax_bottom.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=32)

# Set tick label font size for y-axis
ax_top.tick_params(axis='y', labelsize=32)
ax_mid.tick_params(axis='y', labelsize=32)
ax_bottom.tick_params(axis='y', labelsize=32)

# Labels
fig.text(0.04, 0.5, 'Distance from Chromosome Edge (kb)', va='center', rotation='vertical', fontsize=32)
ax_bottom.set_xlabel('Chromosome', fontsize=32)
ax_top.set_title('Distance of First Protein-Coding Gene from subtelomere-telomere boundary', fontsize=72, fontweight='bold')

# Add legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='#1f77b4', edgecolor='black', label='p arm'),
                   Patch(facecolor='#ff7f0e', edgecolor='black', label='q arm')]
ax_top.legend(handles=legend_elements, loc='upper right', fontsize=32)

# Format y-axis with commas for readability
ax_top.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: format(int(x), ',')))
ax_mid.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: format(int(x), ',')))
ax_bottom.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: format(int(x), ',')))

plt.subplots_adjust(left=0.12, right=0.95, bottom=0.15, top=0.90)
plt.savefig('/Users/yunsawa/Desktop/T2T/CHM13/first_protein_coding_gene_distance_broken_axis.png', dpi=150, bbox_inches='tight')
print("Plot saved to: /Users/yunsawa/Desktop/T2T/CHM13/first_protein_coding_gene_distance_broken_axis.png")
