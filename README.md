# RNA-seq Variant Caller Comparison Pipeline

Comprehensive comparison of four variant calling approaches for detecting somatic mutations in hematologic malignancies using PacBio Iso-Seq single-cell RNA-seq data.

## Quick Start

```bash
# Run all variant callers
cd /data/salomonis-archive/FASTQs/NCI-R01/rna_seq_varcall
./scripts/run_callers.sh all

# Run individual callers
./scripts/run_callers.sh lab_supervised      # 1-2 hours
./scripts/run_callers.sh lab_unsupervised    # 30 min (chr21)
./scripts/run_callers.sh gatk                # 30-60 min
./scripts/run_callers.sh deepvariant         # 1-2 hours

# Evaluate results
python scripts/evaluate_callers.py

# Generate report
python scripts/generate_report.py
```

---

## Overview

### What This Pipeline Does

Compares four variant calling methods on a truth set of 8 known cancer-associated variants:

| Caller | Type | Focus | Language | Runtime |
|--------|------|-------|----------|---------|
| **Lab Scripts (Supervised)** | Targeted | Known mutations + cell-level genotyping | Python | 1-2h |
| **Lab Scripts (Unsupervised)** | Discovery | De novo SNV discovery | Python | 30min-8h |
| **GATK HaplotypeCaller** | Standard | Industry-standard RNA-seq variant calling | Java | 30-60min |
| **DeepVariant** | ML-based | Deep learning variant predictions | Singularity | 1-2h |

### Truth Set (8 Variants)

Located in: `truth_set/truth_set.vcf.gz`

| Gene | Position | Type | VAF |
|------|----------|------|-----|
| RUNX1 | chr21:34799432 | Nonsense | 35% |
| ASXL1 | chr20:32434638 | Frameshift | 8% |
| SETBP1 | chr18:44951948 | Missense | 28% |
| SRSF2 | chr17:76736877 | Missense | 37% |
| CBL | chr11:119242488 | Frameshift | UNK |
| MEF2B | chr19:19145882 | Missense | UNK |
| BCOR | chrX:40057210 | Nonsense | UNK |
| PHF6 | chrX:134415107 | Missense | UNK |

### Input Data

- **BAM**: `/data/salomonis-archive/BAMs/Grimes/scRNA-Seq/KINNEX/5801-diagnosis/scisoseq.mapped.bam` (11GB, indexed)
- **Reference**: `/data/salomonis-archive/genomes/hg38/genome.fa` (3GB)
- **Target regions**: `truth_set/truth_regions.bed` (±1kb around variants)

### Expected Outputs

```
results/
├── supervised_extraction/          # Lab supervised results
│   ├── 5801-diagnosis_complete_analysis.tsv
│   ├── 5801-diagnosis_mutation_matrix.csv
│   └── lab_supervised.vcf.gz
├── global_snv/                     # Lab unsupervised results
│   ├── chr21_output.txt
│   └── lab_unsupervised.vcf.gz
├── gatk/                           # GATK results
│   ├── variants_raw.vcf.gz
│   ├── variants_filtered.vcf.gz
│   └── variants_normalized.vcf.gz
├── deepvariant/                    # DeepVariant results
│   ├── variants.vcf.gz
│   ├── variants.g.vcf.gz
│   └── variants_normalized.vcf.gz
└── comparison/                     # Evaluation results
    ├── summary.json
    ├── comparison.csv
    └── COMPARISON_REPORT.md
```

---

## Detailed Usage

### Phase 1: Environment Setup ✓

**Status**: COMPLETE

Required tools:
- ✓ samtools (installed)
- ✓ GATK 4.x (installed)
- ⏳ bcftools 1.23 (installing)
- ✓ Reference genome (3GB, ready)
- ✓ Test BAM (11GB, indexed, ready)

### Phase 2: Truth Set Preparation ✓

**Status**: COMPLETE

Truth set files created:
- `truth_set/truth_set.vcf` - Uncompressed VCF with 8 variants
- `truth_set/truth_set.vcf.gz` - Compressed VCF
- `truth_set/truth_regions.bed` - Target regions (±1kb)
- `truth_set/test_variants_formatted.txt` - Lab script format

### Phase 3: Run Variant Callers

#### 3.1 Lab Supervised (variant_extraction.py)

```bash
./scripts/run_callers.sh lab_supervised
```

**What it does**:
- Searches for specified mutations in BAM file
- Extracts read-level details per cell barcode
- Creates mutation × cell matrix
- Outputs: complete_analysis.tsv, mutation_matrix.csv

**Runtime**: 1-2 hours
**Output**: `results/supervised_extraction/`

**Manual run** (if needed):
```bash
python /data/salomonis-archive/LabFiles/Nathan/Revio/altanalyze3/altanalyze3/components/bam/variant_extraction.py \
  --sample 5801-diagnosis \
  --bam /data/salomonis-archive/BAMs/Grimes/scRNA-Seq/KINNEX/5801-diagnosis/scisoseq.mapped.bam \
  --mutations truth_set/test_variants_formatted.txt \
  --reference /data/salomonis-archive/genomes/hg38/genome.fa \
  --output-dir results/supervised_extraction
```

#### 3.2 Lab Unsupervised (global_snv.py)

```bash
./scripts/run_callers.sh lab_unsupervised
```

**What it does**:
- De novo SNV discovery (no prior knowledge)
- Filters paralogs automatically
- Reports frequency distributions
- **Note**: SNVs only, misses frameshifts

**Runtime**:
- Chr21 test: 10-30 min
- Full genome: 4-8 hours

**Output**: `results/global_snv/`

**To run full genome** (after chr21 test):
```bash
python /data/salomonis-archive/LabFiles/Nathan/Revio/altanalyze3/altanalyze3/components/bam/global_snv.py \
  /data/salomonis-archive/BAMs/Grimes/scRNA-Seq/KINNEX/5801-diagnosis/scisoseq.mapped.bam \
  /data/salomonis-archive/genomes/hg38/genome.fa \
  results/global_snv/genome_output.txt \
  --min_reads 50 \
  --min_percent 8.0
```

#### 3.3 GATK HaplotypeCaller

```bash
./scripts/run_callers.sh gatk
```

**What it does**:
- Standard GATK RNA-seq variant calling
- Calls both SNVs and indels
- Applies quality filters

**Runtime**: 30-60 minutes
**Output**: `results/gatk/variants_filtered.vcf.gz`

**Manual run** (if needed):
```bash
# Step 1: Call variants
gatk HaplotypeCaller \
  -R /data/salomonis-archive/genomes/hg38/genome.fa \
  -I /data/salomonis-archive/BAMs/Grimes/scRNA-Seq/KINNEX/5801-diagnosis/scisoseq.mapped.bam \
  -O results/gatk/variants_raw.vcf.gz \
  --dont-use-soft-clipped-bases \
  --standard-min-confidence-threshold-for-calling 20.0 \
  --min-base-quality-score 10 \
  -L truth_set/truth_regions.bed \
  --native-pair-hmm-threads 8

# Step 2: Filter
gatk VariantFiltration \
  -R /data/salomonis-archive/genomes/hg38/genome.fa \
  -V results/gatk/variants_raw.vcf.gz \
  -O results/gatk/variants_filtered.vcf.gz \
  --filter-name "LowQual" --filter-expression "QUAL < 30.0" \
  --filter-name "LowDepth" --filter-expression "DP < 10"
```

#### 3.4 DeepVariant

```bash
./scripts/run_callers.sh deepvariant
```

**What it does**:
- Deep learning-based variant calling
- Uses WGS model (no RNA-seq specific model)
- Generates both VCF and gVCF

**Runtime**: 1-2 hours
**Output**: `results/deepvariant/variants.vcf.gz`

**First-time setup** (downloads container):
- Container auto-downloads on first run
- Downloads to: `containers/deepvariant_1.6.1.sif`
- Size: ~2GB

**Manual run** (if needed):
```bash
singularity exec \
  --bind /data:/data \
  containers/deepvariant_1.6.1.sif \
  /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref=/data/salomonis-archive/genomes/hg38/genome.fa \
    --reads=/data/salomonis-archive/BAMs/Grimes/scRNA-Seq/KINNEX/5801-diagnosis/scisoseq.mapped.bam \
    --regions=truth_set/truth_regions.bed \
    --output_vcf=results/deepvariant/variants.vcf.gz \
    --output_gvcf=results/deepvariant/variants.g.vcf.gz \
    --num_shards=8 \
    --make_examples_extra_args="min_mapping_quality=1"
```

### Phase 4: Evaluate Results

```bash
python scripts/evaluate_callers.py
```

**What it does**:
- Loads truth VCF and all caller VCFs
- Matches variants (exact for SNVs, fuzzy for indels)
- Calculates sensitivity, precision, F1
- Measures VAF accuracy (correlation, MAE)
- Outputs JSON and CSV summaries

**Output Files**:
- `results/comparison/summary.json` - Full metrics in JSON
- `results/comparison/comparison.csv` - Comparison table

**Output Example**:
```
Caller              TP FP FN Sensitivity Precision F1    VAF_Corr VAF_MAE
Lab_Supervised     8  0  0   1.000       1.000    1.000 0.988    0.025
Lab_Unsupervised   6  0  2   0.750       1.000    0.857 0.992    0.012
GATK               7  1  1   0.875       0.875    0.875 0.876    0.035
DeepVariant        6  2  2   0.750       0.750    0.750 0.765    0.048
```

### Phase 5: Generate Report

```bash
python scripts/generate_report.py
```

**Output**: `results/comparison/COMPARISON_REPORT.md`

**Report includes**:
- Executive summary
- Caller strengths/weaknesses
- Detailed performance metrics
- Technical considerations
- Production recommendations
- Appendix with software versions

---

## Advanced Usage

### Run Only Specific Chromosomes

**Lab Supervised** (by modifying mutation file):
```bash
# Edit truth_set/test_variants_formatted.txt to include only chr21
python /data/salomonis-archive/LabFiles/Nathan/Revio/altanalyze3/altanalyze3/components/bam/variant_extraction.py ...
```

**Lab Unsupervised** (already has --test_chromosome):
```bash
python /path/to/global_snv.py \
  scisoseq.mapped.bam genome.fa output.txt \
  --test_chromosome chr21
```

### Custom Filtering Parameters

**Lab Unsupervised**:
```bash
# More stringent
--min_reads 100 --min_percent 15.0

# More sensitive
--min_reads 5 --min_percent 2.0
```

**GATK**:
```bash
# Higher quality threshold
--standard-min-confidence-threshold-for-calling 30.0

# More base pairs (removes low quality)
--min-base-quality-score 20
```

### Manual VCF Normalization

If bcftools is available:
```bash
# Sort
bcftools sort -Oz variants.vcf.gz > sorted.vcf.gz

# Normalize (left-align indels)
bcftools norm -f genome.fa normalized.vcf.gz > final.vcf.gz

# Index
bcftools index -t final.vcf.gz
```

### Inspect Specific Variant

```bash
# Check if variant was found by each caller
for vcf in results/*/variants*.vcf.gz; do
  echo "=== $vcf ==="
  zcat $vcf | grep "chr21" | grep "34799432"
done

# Manual inspection with samtools
samtools view -b scisoseq.mapped.bam chr21:34799430-34799435 | samtools tview - genome.fa
```

---

## Troubleshooting

### bcftools Installation Hangs

**Solution**: Install in separate environment
```bash
conda create -n bcftools-env -c bioconda bcftools=1.23 -y
conda activate bcftools-env
# Then use bcftools commands
```

### GATK Out of Memory

**Solution**: Increase memory allocation
```bash
export _JAVA_OPTIONS="-Xmx16g"
gatk HaplotypeCaller ...
```

### DeepVariant Container Download Fails

**Solution**: Manual download
```bash
cd containers/
singularity pull docker://google/deepvariant:1.6.1
# Wait for ~5-10 minutes, 2GB download
```

### Lab Script Import Errors

**Issue**: Missing dependencies in altanalyze3

**Solution**: Check Python environment
```bash
python -c "import pysam; import HTSeq; print('OK')"
```

### Low Variant Detection

**Possible causes**:
1. MAPQ filter too stringent (43% of reads are MAPQ=0)
2. VAF threshold too high
3. Depth too low in regions

**Solution**: Lower thresholds
- `--min-mapping-quality=1` for DeepVariant
- `--min-base-quality-score 10` for GATK
- `--min_percent 5.0` for lab_unsupervised

---

## Timing Guide

| Phase | Component | Runtime | Total |
|-------|-----------|---------|-------|
| 1 | Setup | instant | instant |
| 2 | Truth prep | ~5 min | 5 min |
| 3 | Lab supervised | 1-2 hours | 1-2 hours |
| 3 | Lab unsupervised (chr21) | 30 min | 30 min |
| 3 | Lab unsupervised (full) | 4-8 hours | 4-8 hours |
| 3 | GATK | 30-60 min | 30-60 min |
| 3 | DeepVariant | 1-2 hours | 1-2 hours |
| 4 | Evaluation | ~10 min | 10 min |
| 5 | Report | ~5 min | 5 min |
| **Total (targeted)** | **All except full unsupervised** | | **~4 hours** |
| **Total (comprehensive)** | **Include full unsupervised** | | **~10 hours** |

---

## Key Files

```
/data/salomonis-archive/FASTQs/NCI-R01/rna_seq_varcall/
├── README.md (this file)
├── CLAUDE.md (project instructions)
├── IMPLEMENTATION_STATUS.md (detailed status)
├── test_variants.txt (original truth variants)
├── scripts/
│   ├── run_callers.sh (main runner)
│   ├── prepare_truth_vcf.py (created truth VCF)
│   ├── tsv_to_vcf.py (convert lab outputs)
│   ├── evaluate_callers.py (evaluation)
│   ├── generate_report.py (reporting)
│   └── compress_vcf.py (VCF compression)
├── truth_set/
│   ├── truth_set.vcf (8 variants)
│   ├── truth_set.vcf.gz (compressed)
│   ├── truth_regions.bed (target regions)
│   └── test_variants_formatted.txt (lab format)
└── results/
    ├── supervised_extraction/ (lab supervised)
    ├── global_snv/ (lab unsupervised)
    ├── gatk/ (GATK results)
    ├── deepvariant/ (DeepVariant results)
    └── comparison/ (evaluation results)
```

---

## Dependencies

**Required**:
- Python 3.x with pysam, scipy
- samtools
- GATK 4.x
- Reference genome (3GB)
- Test BAM (11GB, indexed)

**Optional**:
- bcftools (for VCF operations)
- Singularity (for DeepVariant)
- IGV (for manual inspection)

---

## Output Metrics Explained

### Sensitivity (Recall)
- **Formula**: TP / (TP + FN)
- **Meaning**: Of the 8 true variants, how many did the caller find?
- **Range**: 0-1 (higher is better)
- **Example**: 0.875 = found 7 of 8 variants

### Precision
- **Formula**: TP / (TP + FP)
- **Meaning**: Of the variants called, how many are true?
- **Range**: 0-1 (higher is better)
- **Example**: 0.875 = 7 correct, 1 false positive

### F1 Score
- **Formula**: 2 × (Precision × Sensitivity) / (Precision + Sensitivity)
- **Meaning**: Harmonic mean of precision and sensitivity
- **Range**: 0-1 (higher is better)
- **Use**: Balanced metric when both precision and recall matter

### VAF Correlation
- **Formula**: Pearson correlation coefficient
- **Meaning**: How well does called VAF match expected VAF?
- **Range**: -1 to 1 (higher is better)
- **Example**: 0.95 = excellent correlation

### VAF MAE (Mean Absolute Error)
- **Formula**: Average absolute difference between expected and called VAF
- **Meaning**: Average deviation in VAF quantification
- **Range**: 0 to 1 (lower is better)
- **Example**: 0.025 = average error of 2.5% VAF

---

## Contact & Support

For issues with:
- **Lab scripts**: Contact Nathan (Lab Files owner)
- **GATK**: See Broad Institute documentation
- **DeepVariant**: See Google DeepVariant documentation
- **This pipeline**: Check IMPLEMENTATION_STATUS.md

---

**Last Updated**: 2026-02-13
**Status**: Ready for execution
