# Variant Caller Comparison - Implementation Status

## Summary
Implementation of comprehensive variant caller comparison pipeline for RNA-seq data from patient 5801-diagnosis.

**Goal**: Identify best variant caller (GATK, DeepVariant, or Lab Scripts) for production use.

---

## Completed Phases

### Phase 1: Environment Setup ✓
- **Tools Available**:
  - ✓ samtools (installed in bio-cli)
  - ✓ gatk (installed in bio-cli)
  - ⏳ bcftools (installing via conda)
  - ✓ Reference genome (3GB, `/data/salomonis-archive/genomes/hg38/genome.fa`)
  - ✓ Test BAM (11GB, indexed, at `/data/salomonis-archive/BAMs/Grimes/scRNA-Seq/KINNEX/5801-diagnosis/scisoseq.mapped.bam`)

### Phase 2: Truth Set Preparation ✓
- **Created Scripts**:
  - `scripts/prepare_truth_vcf.py` - Converts test_variants.txt to VCF format
  - `scripts/compress_vcf.py` - Compresses VCF files
  - `truth_set/test_variants_formatted.txt` - Formatted variants for lab scripts

- **Created Files**:
  - `truth_set/truth_set.vcf` - Uncompressed truth VCF (8 variants)
  - `truth_set/truth_set.vcf.gz` - Compressed truth VCF
  - `truth_set/truth_regions.bed` - Target regions for callers (±1kb around variants)

**Truth Set Details** (8 known variants):
| Gene | Position | Protein Change | VAF | Type |
|------|----------|-----------------|-----|------|
| RUNX1 | chr21:34799432 | p.W279* | 35% | Nonsense |
| ASXL1 | chr20:32434638 | p.G643fs | 8% | Frameshift |
| SETBP1 | chr18:44951948 | p.G870S | 28% | Missense |
| SRSF2 | chr17:76736877 | p.P95R | 37% | Missense |
| CBL | chr11:119242488 | p.S606fs | UNK | Frameshift |
| MEF2B | chr19:19145882 | p.L341V | UNK | Missense |
| BCOR | chrX:40057210 | p.R1480* | UNK | Nonsense |
| PHF6 | chrX:134415107 | p.R274P | UNK | Missense |

---

## In-Progress Phases

### Phase 3: Run Variant Callers (⏳ Next)

#### 3.1 Lab Script 1: variant_extraction.py (Supervised)

**Location**: `/data/salomonis-archive/LabFiles/Nathan/Revio/altanalyze3/altanalyze3/components/bam/variant_extraction.py`

**Command**:
```bash
python /data/salomonis-archive/LabFiles/Nathan/Revio/altanalyze3/altanalyze3/components/bam/variant_extraction.py \
  --sample 5801-diagnosis \
  --bam /data/salomonis-archive/BAMs/Grimes/scRNA-Seq/KINNEX/5801-diagnosis/scisoseq.mapped.bam \
  --mutations truth_set/test_variants_formatted.txt \
  --reference /data/salomonis-archive/genomes/hg38/genome.fa \
  --output-dir results/supervised_extraction
```

**Expected Output**:
- `results/supervised_extraction/5801-diagnosis_complete_analysis.tsv` (read-level details)
- `results/supervised_extraction/5801-diagnosis_mutation_matrix.csv` (cell × variant matrix)

**Runtime**: 1-2 hours

#### 3.2 Lab Script 2: global_snv.py (Unsupervised Discovery)

**Location**: `/data/salomonis-archive/LabFiles/Nathan/Revio/altanalyze3/altanalyze3/components/bam/global_snv.py`

**Test run (Chr21 only)**:
```bash
python /data/salomonis-archive/LabFiles/Nathan/Revio/altanalyze3/altanalyze3/components/bam/global_snv.py \
  /data/salomonis-archive/BAMs/Grimes/scRNA-Seq/KINNEX/5801-diagnosis/scisoseq.mapped.bam \
  /data/salomonis-archive/genomes/hg38/genome.fa \
  results/global_snv/chr21_output.txt \
  --test_chromosome chr21 \
  --min_reads 10 \
  --min_percent 5.0
```

**Runtime**: 10-30 min (chr21), 4-8 hours (full genome)

#### 3.3 GATK HaplotypeCaller (RNA-seq mode)

**Command**:
```bash
gatk HaplotypeCaller \
  -R /data/salomonis-archive/genomes/hg38/genome.fa \
  -I /data/salomonis-archive/BAMs/Grimes/scRNA-Seq/KINNEX/5801-diagnosis/scisoseq.mapped.bam \
  -O results/gatk/variants_raw.vcf.gz \
  --dont-use-soft-clipped-bases \
  --standard-min-confidence-threshold-for-calling 20.0 \
  --min-base-quality-score 10 \
  -L truth_set/truth_regions.bed \
  --native-pair-hmm-threads 8
```

**Filtering**:
```bash
gatk VariantFiltration \
  -R /data/salomonis-archive/genomes/hg38/genome.fa \
  -V results/gatk/variants_raw.vcf.gz \
  -O results/gatk/variants_filtered.vcf.gz \
  --filter-name "LowQual" --filter-expression "QUAL < 30.0" \
  --filter-name "LowDepth" --filter-expression "DP < 10"
```

**Runtime**: 30-60 min

#### 3.4 DeepVariant (via Singularity)

**Setup** (if needed):
```bash
cd containers/
singularity pull docker://google/deepvariant:1.6.1
```

**Command**:
```bash
singularity exec \
  --bind /data:/data \
  containers/deepvariant_1.6.1.sif \
  /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref=/data/salomonis-archive/genomes/hg38/genome.fa \
    --reads=/data/salomonis-archive/BAMs/Grimes/scRNA-Seq/KINNEX/5801-diagnosis/scisoseq.mapped.bam \
    --regions=/data/salomonis-archive/FASTQs/NCI-R01/rna_seq_varcall/truth_set/truth_regions.bed \
    --output_vcf=results/deepvariant/variants.vcf.gz \
    --output_gvcf=results/deepvariant/variants.g.vcf.gz \
    --num_shards=8 \
    --make_examples_extra_args="min_mapping_quality=1"
```

**Runtime**: 1-2 hours

---

## Planned Phases

### Phase 4: Convert Outputs to VCF (1 hour)
- **Script**: `scripts/tsv_to_vcf.py`
- Convert lab script TSV outputs to standard VCF format
- Normalize indels and coordinate representations

### Phase 5: Evaluation & Comparison (2-4 hours)
- **Script**: `scripts/evaluate_callers.py`
- Calculate: Sensitivity, Precision, F1, VAF accuracy
- Generate comparison matrices and metrics

### Phase 6: Generate Report (1 hour)
- Create markdown comparison report
- Production recommendations
- Performance analysis

---

## Project Structure

```
/data/salomonis-archive/FASTQs/NCI-R01/rna_seq_varcall/
├── CLAUDE.md                           # Project instructions
├── IMPLEMENTATION_STATUS.md            # This file
├── test_variants.txt                   # Original test data (8 variants)
├── scripts/
│   ├── prepare_truth_vcf.py           # ✓ Create truth set VCF
│   ├── compress_vcf.py                 # ✓ Compress VCF files
│   ├── tsv_to_vcf.py                  # Convert lab outputs
│   ├── evaluate_callers.py            # Compare callers
│   └── generate_report.py             # Generate report
├── reference/
│   └── [symlink or copy of genome.fa needed]
├── truth_set/
│   ├── truth_set.vcf                  # ✓ VCF with 8 variants
│   ├── truth_set.vcf.gz               # ✓ Compressed VCF
│   ├── truth_regions.bed              # ✓ Target regions (±1kb)
│   └── test_variants_formatted.txt    # ✓ Formatted for lab scripts
├── results/
│   ├── supervised_extraction/         # Lab supervised output
│   ├── global_snv/                    # Lab unsupervised output
│   ├── gatk/                          # GATK output
│   ├── deepvariant/                   # DeepVariant output
│   └── comparison/                    # Evaluation results
└── containers/
    └── [DeepVariant singularity image if needed]
```

---

## Key Data & Files

**Input Data**:
- BAM: `/data/salomonis-archive/BAMs/Grimes/scRNA-Seq/KINNEX/5801-diagnosis/scisoseq.mapped.bam` (11GB, indexed)
- Genome: `/data/salomonis-archive/genomes/hg38/genome.fa` (3GB, indexed)
- Truth variants: `test_variants.txt` → `truth_set/truth_set.vcf`

**Lab Scripts**:
- Supervised: `/data/salomonis-archive/LabFiles/Nathan/Revio/altanalyze3/altanalyze3/components/bam/variant_extraction.py`
- Unsupervised: `/data/salomonis-archive/LabFiles/Nathan/Revio/altanalyze3/altanalyze3/components/bam/global_snv.py`

---

## Next Steps

1. **Wait for bcftools installation** to complete
2. **Run variant callers in sequence**:
   - Lab supervised (1-2h)
   - Lab unsupervised chr21 test (30min)
   - GATK (30-60min)
   - DeepVariant (1-2h)
3. **Convert outputs** to normalized VCF
4. **Evaluate & compare** using evaluate_callers.py
5. **Generate final report** with recommendations

---

## Expected Outcomes

### Minimum Success Criteria:
- All 4 callers execute without errors
- Detect ≥6/8 truth variants with at least one caller
- Clear comparison metrics available

### Ideal Success Criteria:
- Detect 8/8 variants with concordance across callers
- VAF correlation r > 0.9 for known-VAF variants
- Clear production recommendation

---

## Notes

- **Mapping Quality**: 43% of reads have MAPQ=0 (multi-mappers) - will affect SNV detection
- **Read Length**: ~802bp average (long reads advantage)
- **Cell Count**: ~8-10K cells per sample
- **VAF Known**: 4 of 8 variants (8%, 28%, 35%, 37%)
- **Frameshifts**: 2 variants (ASXL1, CBL) - SNV-only tools will miss these

---

**Last Updated**: 2026-02-13
**Status**: Phase 2 Complete, Phase 3 Ready to Begin
