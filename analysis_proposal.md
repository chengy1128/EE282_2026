# Analysis Proposal: Comparative Subtelomeric Variation in CHM13 and HG002 T2T Assemblies

The completion of telomere-to-telomere (T2T) genome assemblies has transformed the study of the most complex and repetitive regions of the human genome (Nurk et al., 2022). Subtelomeres—defined here as the region from the telomeric repeats inward to the first protein-coding gene on each chromosome arm—are especially rich in segmental duplications, rapidly evolving gene families, and structural variants that can influence adaptation and disease risk (Stong et al., 2014). My project asks a focused comparative question: how conserved are subtelomeric regions between the CHM13 (haploid) and HG002 (diploid) T2T assemblies, and where do they differ most strongly? The core claim of this proposal is that a targeted, data-driven comparison of subtelomeric sequences across all chromosomes will reveal chromosome-specific patterns of divergence, particularly on acrocentric chromosomes where inter-chromosomal exchange is frequent. This analysis is both biologically meaningful and methodologically feasible with the Unix/R workflows emphasized in this course.

The analysis will be grounded in two specific, publicly available datasets: the CHM13v2.0 T2T assembly (Nurk et al., 2022) and the HG002v1.1 T2T assembly (Rhie et al., 2023), both distributed in FASTA format with accompanying GFF3 annotations. These resources provide complete sequence coverage of each chromosome arm and high-quality annotations that can be used to contextualize variation. I will acquire the assemblies from the Human Pangenome Reference Consortium (HPRC; Liao et al., 2023) or official T2T repository mirrors, then store them within the project’s data directory with clear versioned filenames. The emphasis on complete, finished assemblies is essential because subtelomeric regions were poorly represented in older references such as GRCh38 (Schneider et al., 2017); the T2T datasets make this study possible at all.

My analysis strategy is structured in three stages: (1) focused data extraction and harmonization, (2) comparative alignment and feature quantification, and (3) visualization-driven interpretation.

First, I will extract the subtelomeric segments for all chromosome arms. Using Unix utilities (`grep`, `awk`, `cut`, and `sed`), I will parse FASTA headers and GFF3 annotations to standardize chromosome naming, identify telomeric repeat boundaries, and locate the first protein-coding gene on each arm. I will then create per-chromosome-arm subtelomeric FASTA files that span from the telomeric repeats to that first protein-coding gene. This definition grounds the analysis in functional genome architecture while keeping the extraction procedure transparent and reproducible. I will record each extraction step in a shell script or Markdown code block to support reproducibility.

Second, I will conduct a comparative sequence analysis between CHM13 and HG002 subtelomeric regions. My primary comparative metric will be sequence identity and structural variation. For identity, I will align each CHM13 subtelomeric segment to its corresponding HG002 segment using a whole-genome aligner that is efficient for long sequences (e.g., `minimap2` (Li, 2018) or `nucmer` (Marçais et al., 2018), depending on availability). I will parse the alignment output into tabular form, summarizing percent identity and coverage for each chromosome arm. For structural variation, I will compute the size and frequency of insertion/deletion events derived from alignment gaps, and I will cross-reference these events with subtelomeric gene annotations from the GFF3 files. This allows me to test whether high-variation subtelomeres overlap regions enriched for gene families such as olfactory receptors (OR) or zinc finger (ZNF) clusters.

Third, I will use R for exploratory data analysis and visualization, treating visualization as a primary analytic tool. The following figures will be produced:

1. Distribution plots of subtelomeric alignment identity across chromosomes, highlighting outliers and overall dispersion.
2. Heatmaps of percent identity and coverage by chromosome arm, enabling rapid identification of hotspots of divergence.
3. Bar charts or density plots of insertion/deletion size distributions, grouped by chromosome class (e.g., acrocentric vs. non-acrocentric).
4. Overlay plots indicating whether structural variants intersect annotated gene families, which directly connect sequence variation to functional genomic content.

These plots will be generated with `ggplot2` (Wickham, 2016) and will be accompanied by concise quantitative summaries (means, medians, and interquartile ranges) to anchor interpretation. The visualizations are designed to make the data persuasive: if my hypothesis about acrocentric divergence is correct, it should appear as a clear pattern in the heatmaps and as a skewed distribution in acrocentric vs. non-acrocentric comparisons.

This project is feasible within the constraints of the course. The data sources are well-defined, and the processing steps are rooted in standard Unix text and file manipulation. The subtelomere definition (terminal 500 kb) keeps data sizes tractable while preserving the biological signal of interest. Additionally, alignment and parsing can be done at the level of chromosome arms rather than full genomes, which reduces computational requirements without sacrificing interpretability. The only anticipated technical challenge is the repetitive nature of subtelomeric segments, which can lead to multi-mapping or ambiguous alignments. I will mitigate this by prioritizing the longest, highest-confidence alignments and explicitly documenting any filtering thresholds used (e.g., minimum alignment length, minimum identity). This approach provides methodological transparency and keeps the analysis reproducible.

The proposed analyses directly serve the research goals of the project. By quantifying conservation and divergence in a targeted region of the genome, I will answer which chromosome ends are most stable and which are most variable. By measuring the size and frequency of insertions and deletions, I will describe the structural landscape of subtelomeric variation. By integrating annotation overlap, I will connect structural divergence to functional gene families. Finally, by emphasizing visualization, I will provide an interpretable and persuasive argument for how subtelomeric variation differs between two canonical T2T assemblies.

In conclusion, this analysis proposal presents a concrete, feasible strategy to interrogate subtelomeric variation at unprecedented resolution. The approach uses the strengths of the course—Unix-based data handling, R-based visualization, and Git/GitHub reproducibility—to address a high-impact biological question. Comparing CHM13 and HG002 at the subtelomere level will produce both descriptive insights and a compelling set of figures, demonstrating how complete genome assemblies can be leveraged to understand the most dynamic regions of the human genome. The resulting analysis will be concise, methodologically rigorous, and persuasive in its scientific framing.

## References

Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, 34(18), 3094–3100. https://doi.org/10.1093/bioinformatics/bty191

Liao, W.-W., Asri, M., Ebler, J., et al. (2023). A draft human pangenome reference. *Nature*, 617(7960), 312–324. https://doi.org/10.1038/s41586-023-05896-x

Marçais, G., Delcher, A. L., Phillippy, A. M., Coston, R., Salzberg, S. L., & Zimin, A. (2018). MUMmer4: A fast and versatile genome alignment system. *PLOS Computational Biology*, 14(1), e1005944. https://doi.org/10.1371/journal.pcbi.1005944

Nurk, S., Koren, S., Rhie, A., et al. (2022). The complete sequence of a human genome. *Science*, 376(6588), 44–53. https://doi.org/10.1126/science.abj6987

Rhie, A., Nurk, S., Cechova, M., et al. (2023). The complete sequence of a human Y chromosome. *Nature*, 621(7978), 344–354. https://doi.org/10.1038/s41586-023-06457-y

Schneider, V. A., Graves-Lindsay, T., Howe, K., et al. (2017). Evaluation of GRCh38 and de novo haploid genome assemblies demonstrates the enduring quality of the reference assembly. *Genome Research*, 27(5), 849–864. https://doi.org/10.1101/gr.213611.116

Stong, N., Deng, Z., Gupta, R., et al. (2014). Subtelomeric CTCF and cohesin binding site organization using improved subtelomere assemblies and a novel annotation pipeline. *Genome Research*, 24(6), 1039–1050. https://doi.org/10.1101/gr.166983.113

Wickham, H. (2016). *ggplot2: Elegant Graphics for Data Analysis*. Springer-Verlag New York. https://doi.org/10.1007/978-3-319-24277-4
