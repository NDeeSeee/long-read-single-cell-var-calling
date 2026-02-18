# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository implements an RNA-seq variant calling pipeline for detecting somatic mutations in hematologic malignancies, leveraging **AlphaGenome** for genome-wide analysis and variant interpretation. The pipeline processes single-cell and bulk RNA-seq data to identify clinically relevant variants in genes commonly mutated in myeloid disorders (BCOR, RUNX1, ASXL1, SRSF2, CBL, MEF2B, PHF6, SETBP1) and provides variant annotation with allele frequency estimation.

### Key Design Considerations

- **AlphaGenome Integration**: Variant predictions and genomic context interpretation should leverage AlphaGenome models/annotations where applicable for improved accuracy and clinical utility
- **Multi-sample Scale**: Pipeline designed to handle large batch processing from the Grimes scRNA-Seq KINNEX dataset
- **Hematologic Oncology Focus**: Emphasis on mutations relevant to AML, other myeloid malignancies, and leukemias
- **VAF Quantification**: Accurate variant allele frequency (VAF) estimation critical for clonal heterogeneity assessment
- **sc/bulk RNA Duality**: Support both single-cell and bulk RNA-seq input formats

## Test Data

Test BAM files are located at: `/data/salomonis-archive/BAMs/Grimes/scRNA-Seq/KINNEX/5801-diagnosis`

Validation variants are in `test_variants.txt` - includes 8 recurrent mutations across key cancer-associated genes with expected VAF ranges (8-37%).

## Development Environment

- **Conda Environment**: `bio-cli` (check activation in local setup)
- **Language**: Python 3.x with bioinformatics libraries (BioPython, pysam, HTSeq, etc.)
- **Supporting Tools**: samtools, bcftools, potentially GATK for variant refinement
- **Pipeline Framework**: (To be determined - Snakemake, Nextflow, or custom orchestration)

## Pipeline Architecture (High-Level)

1. **Input**: BAM files (indexed), reference genome, GTF/GFF annotation
2. **Variant Detection**: SNV/indel calling from RNA-seq alignments
3. **AlphaGenome Interpretation**: Functional annotation and prediction scoring
4. **Annotation**: VEP or SnpEff integration with AlphaGenome scores
5. **Filtering**: Clinical relevance filtering (known hotspots, VAF thresholds)
6. **Output**: Annotated VCF with VAF, AlphaGenome scores, gene context

## Common Commands (To Be Implemented)

Once pipeline is established:
- `python pipeline.py --bam <file> --output <dir>` - Run full variant calling
- `python extract_variants.py --bam <file> --genes BCOR,RUNX1` - Gene-specific extraction
- `python validate_variants.py --vcf <file> --truth test_variants.txt` - Benchmarking
- Testing: `pytest tests/` or similar

## Important Notes

- Variant coordinates in test_variants.txt are 1-based genomic positions
- VAF "UNK" indicates variants where frequency was not measured in original data
- Clinical variants include frameshift mutations (fs), nonsense mutations (*), and missense changes
- Integration with AlphaGenome should prioritize genome-wide context, not just gene-level predictions
