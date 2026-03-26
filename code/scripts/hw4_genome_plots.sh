#!/usr/bin/env bash
set -euo pipefail

# Homework 4 Part 2: Plots for genome partitions
# For each partition (<=100kb, >100kb): length histogram, GC% histogram, cumulative size plot

FASTA="homework3/data/fasta/dmel-all-chromosome-r6.66.fasta.gz"
OUTDIR="homework4/plots"
DATADIR="homework4/results"
mkdir -p "$OUTDIR" "$DATADIR"

THRESHOLD=100000

echo "=== Extracting per-sequence length and GC% with bioawk ==="

# bioawk gc() returns GC fraction; compute length, gc%, and partition
bioawk -c fastx '{
    len = length($seq);
    gcpct = gc($seq) * 100;
    if (len <= 100000)
        print "le100kb", len, gcpct;
    else
        print "gt100kb", len, gcpct;
}' <(gunzip -c "$FASTA") > "${DATADIR}/seq_length_gc.tsv"

echo "=== Generating plots with R ==="

Rscript - <<'RSCRIPT'
library(ggplot2)

data <- read.table("homework4/results/seq_length_gc.tsv",
                    header = FALSE, sep = "\t",
                    col.names = c("partition", "length", "gc_pct"))

partitions <- list(
    list(label = "le100kb", title = "Sequences <= 100kb"),
    list(label = "gt100kb", title = "Sequences > 100kb")
)

for (p in partitions) {
    sub <- data[data$partition == p$label, ]

    # 1. Sequence length distribution (log scale)
    png(paste0("homework4/plots/length_hist_", p$label, ".png"),
        width = 800, height = 600)
    print(
        ggplot(sub, aes(x = length)) +
            geom_histogram(bins = 50, fill = "steelblue", color = "black") +
            scale_x_log10() +
            labs(title = paste("Sequence Length Distribution -", p$title),
                 x = "Sequence Length (bp, log scale)",
                 y = "Count") +
            theme_minimal()
    )
    dev.off()

    # 2. GC% distribution
    png(paste0("homework4/plots/gc_hist_", p$label, ".png"),
        width = 800, height = 600)
    print(
        ggplot(sub, aes(x = gc_pct)) +
            geom_histogram(bins = 50, fill = "darkorange", color = "black") +
            labs(title = paste("GC% Distribution -", p$title),
                 x = "GC%",
                 y = "Count") +
            theme_minimal()
    )
    dev.off()

    # 3. Cumulative sequence size (sorted largest to smallest)
    sorted_lens <- sort(sub$length, decreasing = TRUE)
    cumsum_lens <- cumsum(as.numeric(sorted_lens))
    cdf_data <- data.frame(index = seq_along(cumsum_lens),
                           cumsize = cumsum_lens)

    png(paste0("homework4/plots/cumsize_", p$label, ".png"),
        width = 800, height = 600)
    print(
        ggplot(cdf_data, aes(x = index, y = cumsize)) +
            geom_line(color = "darkgreen", linewidth = 1.2) +
            labs(title = paste("Cumulative Sequence Size -", p$title),
                 x = "Sequence Index (sorted largest to smallest)",
                 y = "Cumulative Size (bp)") +
            theme_minimal()
    )
    dev.off()
}

cat("All 6 plots generated.\n")
RSCRIPT

echo "=== Done. Plots saved to ${OUTDIR}/ ==="
