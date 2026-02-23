# Homework 3: Pipelines (Drosophila melanogaster)

## Notes on Tooling
`faSize` was not available in my environment, so I used `seqkit` as a replacement for genome size/sequence summaries.

- `seqkit stats` was used for:
  - total number of nucleotides (`sum_len`)
  - total number of sequences (`num_seqs`)
- `seqkit seq -s | awk` was used for:
  - total number of `N`/`n` bases

`bioawk` was also not available, so equivalent `awk` commands were used to parse GTF columns.

## Data Sources (FlyBase current release used)
- Release: `dmel_r6.66_FB2025_05`
- Genome (all chromosomes): `dmel-all-chromosome-r6.66.fasta.gz`
- Annotation: `dmel-all-r6.66.gtf.gz`

## Scripts
All scripts are in `homework3/scripts`:

1. `homework3/scripts/01_download_data.sh` - download genome/annotation + checksum files
2. `homework3/scripts/02_verify_integrity.sh` - checksum verification
3. `homework3/scripts/03_summarize_genome.sh` - genome summaries using `seqkit`
4. `homework3/scripts/04_summarize_annotation.sh` - annotation summaries using `awk`
5. `homework3/scripts/05_run_all.sh` - run full workflow end-to-end

## File Integrity
### Genome FASTA
- Expected MD5: `ccb86e94117eb4eeaaf70efb6be1b6b9` for `dmel-all-chromosome-r6.66.fasta.gz`
- Verification result: `OK` (`homework3/results/genome_integrity.txt`)

### Annotation GTF
- Expected MD5: `ea600dbb86f1779463f69082131753cd` for `dmel-all-r6.66.gtf.gz`
- Verification result: `OK` (`homework3/results/annotation_integrity.txt`)

## Genome Assembly Summary
From `homework3/results/genome_summary.tsv`:

- Total number of nucleotides: **143726002**
- Total number of Ns: **1152978**
- Total number of sequences: **1870**

## Annotation Summary
### Total number of features by type (most to least common)
From `homework3/results/annotation_feature_counts.tsv`:

| count | feature |
|---:|---|
| 190176 | exon |
| 163377 | CDS |
| 46856 | 5UTR |
| 33778 | 3UTR |
| 30922 | start_codon |
| 30862 | stop_codon |
| 30836 | mRNA |
| 17872 | gene |
| 3059 | ncRNA |
| 485 | miRNA |
| 365 | pseudogene |
| 312 | tRNA |
| 270 | snoRNA |
| 262 | pre_miRNA |
| 115 | rRNA |
| 32 | snRNA |

### Total number of genes per chromosome arm (X, Y, 2L, 2R, 3L, 3R, 4)
From `homework3/results/annotation_genes_per_arm.tsv`:

| chromosome arm | genes |
|---|---:|
| X | 2704 |
| Y | 113 |
| 2L | 3508 |
| 2R | 3649 |
| 3L | 3481 |
| 3R | 4226 |
| 4 | 114 |

## Reproducibility
Run the full pipeline from repo root:

```bash
bash homework3/scripts/05_run_all.sh
```
