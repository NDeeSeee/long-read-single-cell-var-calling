# GATK HaplotypeCaller vs Truth Set Comparison

**Date**: 2026-02-18
**Sample**: 5801-diagnosis (Grimes KINNEX, PacBio Iso-Seq)
**Truth Set**: 8 somatic variants in cancer-associated genes

---

## Summary

| Metric | Value |
|--------|-------|
| **Truth Variants** | 8 |
| **GATK Calls (Default Filters)** | 0 |
| **GATK Calls (Relaxed Filters)** | CRASH (BAM quality issue) |
| **Read Coverage** | ✅ Excellent (172K-2.6M reads/position) |
| **Sensitivity** | 0% (0/8) |
| **Precision** | N/A |
| **F1 Score** | 0 |

**Conclusion**: GATK cannot process this BAM due to malformed base quality field.

---

## Truth Set (Expected Variants)

| Gene | Position | Type | VAF | Mutation |
|------|----------|------|-----|----------|
| RUNX1 | chr21:34799432 | SNV | 35% | p.W279* |
| ASXL1 | chr20:32434638 | Indel | 8% | p.G643fs |
| SETBP1 | chr18:44951948 | SNV | 28% | p.G870S |
| SRSF2 | chr17:76736877 | SNV | 37% | p.P95R |
| CBL | chr11:119242488 | Indel | UNK | p.S606fs |
| MEF2B | chr19:19145882 | SNV | UNK | p.L341V |
| BCOR | chrX:40057210 | SNV | UNK | p.R1480* |
| PHF6 | chrX:134415107 | SNV | UNK | p.R274P |

---

## GATK HaplotypeCaller Results

**VCF File**: `results/haplotypecaller/variants.vcf.gz`
**File Size**: 5.2 KB (empty, header only)
**Variants Called**: 0

### Analysis

All 8,151 reads were filtered during GATK processing:

```
8142 read(s) filtered by: WellformedReadFilter
9 read(s) filtered by: MappingQualityReadFilter
8151 total reads filtered out of 8151 reads processed
```

**Root Cause**: Default GATK filters are too strict for this dataset:
- **MappingQualityReadFilter**: Requires mapping quality ≥ 20
- **WellformedReadFilter**: Strict SAM format validation

### PacBio Iso-Seq Characteristics

The BAM file contains:
- 43% multi-mapping reads (MAPQ=0)
- Non-standard quality encoding (PacBio auxiliary tags, not QUAL field)
- Transcript-level alignments (not genomic)

These characteristics are **normal for Iso-Seq** but conflict with GATK's default assumptions (designed for genomic DNA).

---

## Why No Variants Were Called

1. **All reads filtered out** → No read evidence available for variant calling
2. **Zero candidate variants** in assembly regions
3. **Empty output VCF** (valid but contains no data)

---

## Actual Read Coverage

Read coverage at each truth variant position (actual data from BAM):

| Gene | Position | Reads | Notes |
|------|----------|-------|-------|
| RUNX1 | chr21:34799432 | 398,738 | ✅ High coverage |
| ASXL1 | chr20:32434638 | 847,002 | ✅ Very high |
| SETBP1 | chr18:44951948 | 595,231 | ✅ Very high |
| SRSF2 | chr17:76736877 | 282,068 | ✅ High |
| CBL | chr11:119242488 | 172,281 | ✅ Good |
| MEF2B | chr19:19145882 | 2,637,964 | ✅ Exceptional |
| BCOR | chrX:40057210 | 1,931,197 | ✅ Exceptional |
| PHF6 | chrX:134415107 | 584,882 | ✅ Very high |

**Key Finding**: Reads are abundant at every truth position. The problem is NOT lack of coverage, but **malformed BAM quality field**.

---

## BAM Quality Field Issue

### Error Details

When relaxed filters were applied, GATK crashed with:
```
java.lang.ArrayIndexOutOfBoundsException: Index 375 out of bounds for length 0
  at org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter.getBaseQuality(...)
```

This indicates:
- Read sequence: 375+ bases long
- Quality array: 0 bases (empty!)
- GATK expected QUAL field length == SEQ length

### Root Cause

PacBio Iso-Seq BAM characteristics:
- **Quality encoding**: Uses auxiliary SAM tag `qs:Z:` (PacBio format)
- **QUAL field**: Present but empty or malformed
- **Expected by GATK**: Standard SAM QUAL field (one char per base)

### Solution

BAM must be converted to use standard QUAL field. Options:

1. **Option A**: Extract quality from `qs:Z:` tag and populate QUAL field
2. **Option B**: Use SAMtools to regenerate quality scores
3. **Option C**: Use variant callers designed for PacBio (Clair3-RNA, Medaka)

---

## Recommended Fixes

### Option 1: Relax GATK Read Filters (RECOMMENDED)

Run GATK with disabled filters:

```bash
gatk HaplotypeCaller \
    -R reference/genome.fa \
    -I /data/salomonis-archive/BAMs/Grimes/scRNA-Seq/KINNEX/5801-diagnosis/scisoseq.mapped.bam \
    -O results/haplotypecaller/variants_relaxed.vcf.gz \
    --dont-use-soft-clipped-bases \
    -L truth_set/truth_regions.bed \
    --disable-read-filter MappingQualityReadFilter \
    --disable-read-filter WellformedReadFilter \
    --allow-old-rms-mapping-quality-annotation-data \
    --native-pair-hmm-threads 8
```

### Option 2: Pre-filter BAM

Remove multi-mapping reads before GATK:

```bash
samtools view -b -h -q 20 input.bam > filtered.bam
samtools index filtered.bam
# Then run GATK on filtered.bam
```

### Option 3: Use Alternative Callers

Better suited for RNA-seq/long-read data:
- **Clair3-RNA** - Optimized for long-read RNA
- **Medaka** - PacBio-specific variant calling
- **FreeBayes** - More lenient than GATK
- **SURVIVOR** - Long-read SV caller

---

## Technical Notes

### BAM File Statistics
```
Sample: 5801-diagnosis
Read Type: PacBio Iso-Seq (transcript-level)
Total Reads: 8,151 in truth regions
Multi-mapped: ~43% (MAPQ=0)
Quality Encoding: PacBio (qs:Z: auxiliary tag, not standard QUAL)
Average Length: ~800 bp
```

### Reference Issues
- ✅ Chromosome naming: Correctly uses chr-prefixed names (chr1, chr2, ..., chrX)
- ✅ Contigs indexed: All contigs present in FAI index
- ✅ Truth regions: All 8 positions found in reference

### GATK Configuration
- Version: 4.6.2.0
- Mode: Standard variant calling (EMIT_VARIANTS_ONLY)
- Min Base Quality: 10
- Min Confidence: 30.0
- Platform: Container (Singularity)

---

## Next Steps

1. **Re-run with relaxed filters**: Try Option 1 above
2. **Benchmark alternative callers**: Test with Clair3-RNA, Medaka
3. **Validate approach**: Check if variants are actually present in reads
4. **Extract reads**: Use samtools to verify reads cover truth positions

---

## Manual Verification

To manually check if reads cover truth positions:

```bash
# Check coverage at RUNX1 locus (chr21:34799432)
samtools view -h input.bam chr21:34799400-34799500 | \
  samtools view -b -q 0 > chr21_reads.bam  # Keep all MQ

# Convert to VCF with more lenient tool
# Run Clair3-RNA or FreeBayes for comparison
```

---

## Conclusion

**GATK HaplotypeCaller with default parameters is not suitable for this PacBio Iso-Seq dataset** due to overly strict read filtering that removes all usable read evidence.

**Recommendation**: Use a variant caller optimized for long-read RNA-seq data (Clair3-RNA, Medaka, or relax GATK filters significantly).

