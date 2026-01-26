# Comparative Analysis of Subtelomeric Structural Variation between CHM13 and HG002 T2T Assemblies

## Introduction  - short

The completion of the first "Telomere-to-Telomere" (T2T) human genome assembly, CHM13, marked a paradigm shift in genomics by resolving previously inaccessible "dark matter" regions, including centromeres and subtelomeres (Nurk et al. 2022). Subtelomeres—the regions immediately proximal to the hexameric  telomeric repeats—are characterized by extreme structural variation, segmental duplications, and rapid evolutionary flux. These regions often harbor gene families critical for immunity and environmental adaptation, such as the olfactory receptor (OR) and zinc finger (ZNF) families. Despite their biological importance, these regions remained unresolved in the GRCh38 reference due to the limitations of short-read sequencing.

With the recent availability of the diploid T2T assembly from the HG002 cell line, we now have the unprecedented opportunity to perform a high-resolution comparative analysis between two complete human genomes. **The primary aim of this project is to quantify the conservation and variation of subtelomeric sequences between the CHM13 and HG002 assemblies. I hypothesize that subtelomeric regions will exhibit chromosome-specific patterns of divergence, with acrocentric chromosomes showing higher levels of structural variation due to frequent inter-chromosomal exchanges.**

## Data and Methodology

The primary datasets for this research include the **CHM13v2.0** (complete hydatidiform mole) and **HG002v1.1** (Ashkenazi Jewish male) T2T assemblies. These datasets represent the current gold standard for genomic completeness and are provided in FASTA format, accompanied by GFF3 annotation files.

### Unix-Based Data Processing

To manage these large-scale genomic files, I will utilize the **Unix/Bash** shell techniques emphasized in this course. Processing T2T assemblies requires efficient stream editing and file parsing to avoid memory overflows. I will use `grep` and `awk` to filter the GFF3 files for telomere-associated features and `sed` to standardize chromosome nomenclature across both datasets. By leveraging the command line, I will extract the terminal 500 kb of each chromosome arm to create a focused subset of the genome for downstream analysis. This "data munging" phase is critical for transforming raw genomic data into a format suitable for statistical analysis in R.

### Analysis and Visualization in R

Following data extraction, I will employ **R** for Exploratory Data Analysis (EDA). The central theme of the analysis will be the use of visualization to identify biological patterns. I will use the `ggplot2` package to generate:

1. **Distribution Plots:** To visualize the frequency and size distribution of structural variants (SVs) across the subtelomeres.
2. **Identity Heatmaps:** To compare sequence conservation percentages between CHM13 and HG002 across all 46 chromosomes.
These visualizations are not merely illustrative; they are essential EDA tools that allow for the detection of outliers and assembly artifacts, directly aligning with the course’s focus on using R to summarize and interpret complex data.

### Version Control and Reproducibility

In accordance with professional bioinformatics standards, all scripts, configuration files, and documentation will be managed via **Git and GitHub**. I will use a branch-based workflow to track changes and ensure that the final analysis is fully reproducible. The project will be hosted in a public GitHub repository, utilizing Markdown for all reporting and documentation.

## Research Questions

The proposed analysis will address the following specific questions:

1. **Conservation vs. Divergence:** Which specific human chromosome ends exhibit the highest degree of sequence conservation, and which show the most significant structural divergence between these two benchmark individuals?
2. **Structural Variant (SV) Profiling:** What is the frequency and size distribution of large-scale insertions or deletions within the subtelomeric regions, and do these SVs overlap with known coding sequences?
3. **Visualization as an Analytical Tool:** Can R-based visualizations effectively identify "hotspots" of variation that would be missed by simple summary statistics?

## Feasibility and Anticipated Challenges

This project is highly feasible within the timeframe of this course. The required T2T assemblies are publicly available through the Human Pangenome Reference Consortium (HPRC) and are already "finished" and polished, which eliminates the high computational cost of de novo assembly.

However, I anticipate certain technical challenges. Managing the high disk-usage requirements of T2T FASTA files (approximately 3GB per haploid genome) requires careful file management within the Unix environment to avoid exceeding storage quotas. Additionally, the highly repetitive nature of subtelomeric segmental duplications may complicate initial alignments. To mitigate this, I will focus on high-confidence "anchored" sequences within the subtelomeres and utilize the version control methods learned in class to document any filtering steps I take to ensure reproducibility.

## Conclusion

Comparing the subtelomeres of CHM13 and HG002 provides a unique opportunity to apply foundational bioinformatics skills to a cutting-edge genomic problem. By integrating **Unix-based file processing**, **R-based statistical visualization**, and **Git-based project management**, this project fulfills the core learning outcomes of the course while addressing a significant gap in our understanding of human genomic architecture. Identifying patterns of variation at the ends of human chromosomes will contribute to the broader goal of characterizing the "complete" human pangenome, demonstrating the power of T2T genomics in personal medicine and evolutionary biology.

## References

* Nurk, S., Koren, S., Rhie, A., Rautiainen, M., Bzikadze, A. V., Mikheenko, A., ... & Phillippy, A. M. (2022). The complete sequence of a human genome. *Science* **376**: 44-53.
* Jarvis, E. D., Formenti, G., Rhie, A., Guarracino, A., Yang, C., Wood, J., ... & Phillippy, A. M. (2022). Semi-automated assembly of complete diploid human genomes. *Nature* **602**: 110-120.
* Wang, T., Antonacci-Fulton, L., Howe, K., Lawson, H. A., Lucas, J. K., Phillippy, A. M., ... & Miga, K. H. (2022). The Human Pangenome Project: a resilient reference for the modern genome era. *Nature* **604**: 437-442.
