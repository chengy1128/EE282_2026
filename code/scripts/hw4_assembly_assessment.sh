#!/usr/bin/env bash
set -euo pipefail

# Homework 4 Part 4: Assembly assessment (run on HPC3)
# Prerequisites: hw4_assembly.sh has been run successfully
# Usage: srun -A class_ee282 --cpus-per-task=8 --mem=16G --time=2:00:00 --pty bash hw4_assembly_assessment.sh

WORKDIR="$HOME/ee282/homework4"
ASMDIR="${WORKDIR}/assembly"
RESULTSDIR="${WORKDIR}/results"
PLOTSDIR="${WORKDIR}/plots"
mkdir -p "$RESULTSDIR" "$PLOTSDIR"

MY_ASM="${ASMDIR}/dmel_hifi.bp.p_ctg.fa"
FLYBASE_FASTA="${WORKDIR}/data/dmel-all-chromosome-r6.66.fasta.gz"

# --- Download FlyBase reference if not present ---
if [ ! -f "$FLYBASE_FASTA" ]; then
    echo "=== Downloading FlyBase genome ==="
    mkdir -p "${WORKDIR}/data"
    wget -O "$FLYBASE_FASTA" \
        "https://ftp.flybase.net/releases/current/dmel_r6.66/fasta/dmel-all-chromosome-r6.66.fasta.gz"
fi

# ===================================================================
# Part 4a: Calculate N50 of your assembly
# ===================================================================
echo "=== Calculating N50 ==="

bioawk -c fastx '{print length($seq)}' "$MY_ASM" \
| sort -rn \
| awk '{
    sum += $1; lengths[NR] = $1
}
END {
    half = sum / 2;
    cumsum = 0;
    for (i = 1; i <= NR; i++) {
        cumsum += lengths[i];
        if (cumsum >= half) {
            print "Assembly N50:", lengths[i];
            print "Total assembly size:", sum;
            print "Number of contigs:", NR;
            break;
        }
    }
}' | tee "${RESULTSDIR}/n50.txt"

echo ""
echo "FlyBase Drosophila melanogaster reference contig N50 can be found at:"
echo "https://flybase.org/reports/FBgn0000001 -> genome stats page"
echo "(Check FlyBase genome page for current contig N50 for comparison)"

# ===================================================================
# Part 4b: Contiguity plot comparison
# ===================================================================
echo ""
echo "=== Generating contiguity plots ==="

# My assembly contig sizes
bioawk -c fastx '{print length($seq)}' "$MY_ASM" | sort -rn \
    > "${RESULTSDIR}/my_asm_sizes.txt"

# FlyBase scaffold sizes
bioawk -c fastx '{print length($seq)}' <(gunzip -c "$FLYBASE_FASTA") | sort -rn \
    > "${RESULTSDIR}/flybase_scaffold_sizes.txt"

# FlyBase contig sizes (split scaffolds at runs of Ns)
bioawk -c fastx '{
    n = split($seq, contigs, /[Nn]+/);
    for (i = 1; i <= n; i++) {
        len = length(contigs[i]);
        if (len > 0) print len;
    }
}' <(gunzip -c "$FLYBASE_FASTA") | sort -rn \
    > "${RESULTSDIR}/flybase_contig_sizes.txt"

# Generate contiguity plot with R
Rscript - "${RESULTSDIR}" "${PLOTSDIR}" <<'RSCRIPT'
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
resultsdir <- args[1]
plotsdir <- args[2]

my_sizes    <- scan(file.path(resultsdir, "my_asm_sizes.txt"), quiet = TRUE)
fb_scaffold <- scan(file.path(resultsdir, "flybase_scaffold_sizes.txt"), quiet = TRUE)
fb_contig   <- scan(file.path(resultsdir, "flybase_contig_sizes.txt"), quiet = TRUE)

make_cdf <- function(sizes) {
    sorted <- sort(sizes, decreasing = TRUE)
    data.frame(size = sorted, cumsum = cumsum(as.numeric(sorted)))
}

cdf_my <- make_cdf(my_sizes);         cdf_my$label <- "My HiFi Assembly"
cdf_fb_scaf <- make_cdf(fb_scaffold); cdf_fb_scaf$label <- "FlyBase Scaffold"
cdf_fb_ctg <- make_cdf(fb_contig);    cdf_fb_ctg$label <- "FlyBase Contig"

all_data <- rbind(cdf_my, cdf_fb_scaf, cdf_fb_ctg)

png(file.path(plotsdir, "contiguity_plot.png"), width = 900, height = 600)
ggplot(all_data, aes(x = cumsum, y = size, color = label)) +
    geom_line(linewidth = 1.2) +
    scale_x_continuous(labels = function(x) paste0(x / 1e6, "M")) +
    scale_y_log10(labels = function(x) paste0(x / 1e6, "M")) +
    labs(title = "Contiguity Plot: Assembly Comparison",
         x = "Cumulative Assembly Size (Mbp)",
         y = "Contig/Scaffold Size (bp, log scale)",
         color = "Assembly") +
    theme_minimal()
dev.off()
cat("Contiguity plot saved.\n")
RSCRIPT

# ===================================================================
# Part 4c: BUSCO assessment
# ===================================================================
echo ""
echo "=== Running BUSCO ==="
echo "Prerequisite: conda create -n busco -c conda-forge -c bioconda busco"

# Activate BUSCO environment
# Uncomment and adjust if needed:
# conda activate busco

# BUSCO for my assembly
echo "--- BUSCO on my assembly ---"
busco -i "$MY_ASM" \
      -o busco_my_asm \
      -m genome \
      -l diptera_odb10 \
      --out_path "${RESULTSDIR}" \
      -c 8 \
      -f

# Prepare FlyBase contig-level FASTA for fair comparison
FLYBASE_CONTIGS="${WORKDIR}/data/flybase_contigs.fa"
if [ ! -f "$FLYBASE_CONTIGS" ]; then
    echo "--- Splitting FlyBase scaffolds into contigs ---"
    bioawk -c fastx '{
        n = split($seq, contigs, /[Nn]+/);
        for (i = 1; i <= n; i++) {
            if (length(contigs[i]) > 0) {
                print ">" $name "_contig" i;
                print contigs[i];
            }
        }
    }' <(gunzip -c "$FLYBASE_FASTA") > "$FLYBASE_CONTIGS"
fi

# BUSCO for FlyBase reference
echo "--- BUSCO on FlyBase reference ---"
busco -i "$FLYBASE_CONTIGS" \
      -o busco_flybase \
      -m genome \
      -l diptera_odb10 \
      --out_path "${RESULTSDIR}" \
      -c 8 \
      -f

echo ""
echo "=== BUSCO results ==="
echo "My assembly:"
cat "${RESULTSDIR}/busco_my_asm/short_summary"*.txt 2>/dev/null \
    || echo "Check ${RESULTSDIR}/busco_my_asm/ for results"
echo ""
echo "FlyBase reference:"
cat "${RESULTSDIR}/busco_flybase/short_summary"*.txt 2>/dev/null \
    || echo "Check ${RESULTSDIR}/busco_flybase/ for results"

echo ""
echo "=== Assessment complete ==="
echo "Copy results back to your repo: homework4/results/ and homework4/plots/"
