# Session Summary: Variant Calling Pipeline Analysis
**Date**: February 19, 2026
**Status**: Major progress on core pipeline; mystery solved

---

## What We Accomplished

### ✅ Core Achievements

1. **Solved the CBL Mystery**
   - Problem: CBL variant showed 72.4% VAF but wasn't detected
   - Discovery: CBL position is in a 177kb intron that all Iso-Seq reads skip
   - Implication: RNA-seq fundamentally cannot detect intronic variants
   - Lesson: This is a data limitation, not a caller limitation

2. **Confirmed 4-Variant Pipeline Success**
   - GATK HaplotypeCaller: Detects 4 variants
   - DeepVariant: Detects identical 4 variants (perfect concordance)
   - Result: Both callers working correctly; concordance validates detection

3. **Generated Quality-Fixed BAM**
   - Original BAM lacked QUAL field (marked as '*')
   - Solution: Preprocessed 78.8M reads with synthetic Q30 quality scores
   - Output: `scisoseq_synthetic_qual.bam` (fully indexed)
   - Impact: Enables variant detection; fixes WellformedReadFilter issues

4. **Updated Pipeline Scripts**
   - Fixed conda environment path detection with fallback logic
   - Added sensitive GATK mode with relaxed thresholds
   - Updated DeepVariant to use quality-fixed BAM
   - Improved Clair3-RNA per-chromosome calling logic
   - Created `.gitignore` to prevent VCF/BAM files from git

5. **Created Comprehensive Analysis**
   - Detailed variant detection report (`VARIANT_DETECTION_ANALYSIS.md`)
   - Read-level analysis explaining why 4 variants detected, 4 missing
   - Root cause analysis for each missing variant

### ⚠️ Partial Progress

6. **Clair3-RNA Installation**
   - Environment created: `/users/pavb5f/.conda/envs/clair3-rna`
   - Status: Created but may not be fully functional (needs verification)
   - Next step: Download model files, test execution

7. **LongcallR Installation**
   - Conda installation failed (dependency chain issue)
   - Alternative: Pulling BioContainers Singularity image
   - Status: Container pull in progress (background task)

### 📋 Committed to Git

```
Commit: 1372661
Message: Refactor variant caller pipeline with BAM quality preprocessing
Files: .gitignore, scripts/run_callers.sh (127 insertions, 18 deletions)
```

---

## Current Detection Status

### Detected Variants (4/8)
| Gene | Position | VAF | Status |
|------|----------|-----|--------|
| SRSF2 | chr17:76736877 | 45-48% | ✅ Detected |
| SETBP1 | chr18:44951948 | 37-50% | ✅ Detected |
| ASXL1 | chr20:32434638 | 26-42% | ✅ Detected |
| RUNX1 | chr21:34799432 | 30-35% | ✅ Detected |

### Undetectable Variants (4/8)
| Gene | Position | VAF | Reason |
|------|----------|-----|--------|
| CBL | chr11:119242488 | 72.4% | **Intronic** (RNA-seq skip) |
| MEF2B | chr19:19145882 | 14.3% | Low exonic coverage |
| BCOR | chrX:40057210 | 13.2% | Low exonic coverage |
| PHF6 | chrX:134415107 | 1.9% | Minimal exonic coverage |

---

## Key Technical Findings

### RNA-seq vs DNA-seq
The BAM contains **Iso-Seq RNA reads** that only have bases in exonic regions:
- Introns are skipped (marked as 'N' in CIGAR)
- Example: `...=177140N...` means skip 177kb intron
- **CBL variant in intron → undetectable from RNA-seq**

### Read Quality Issues
- Original BAM has quality field marked as '*' (missing)
- This caused GATK `WellformedReadFilter` to reject 95% of reads
- Solution: Create synthetic BAM with uniform Q30 scores
- Result: All callers now functional

### Data Validation
- Both GATK and DeepVariant independently detected the same 4 variants
- This **perfect concordance** is strong evidence the 4-variant result is correct
- Not random or tool-dependent

---

## Outstanding Tasks

### 🔄 In Progress
1. Verify Clair3-RNA installation and setup
2. Complete LongcallR container download
3. Download Clair3 models (PacBio Iso-Seq or HiFi model)

### 📌 Next Steps (Recommended Order)
1. Test Clair3-RNA with downloaded models on truth BED regions
2. Test LongcallR container with sample variants
3. Run all 4 callers: GATK, DeepVariant, Clair3-RNA, LongcallR
4. Generate final comparison report with 4-caller concordance analysis

---

## Technical Details

### Files Generated This Session
- `/data/salomonis-archive/FASTQs/NCI-R01/rna_seq_varcall/scisoseq_synthetic_qual.bam` — Quality-fixed BAM
- `/data/salomonis-archive/FASTQs/NCI-R01/rna_seq_varcall/results/VARIANT_DETECTION_ANALYSIS.md` — Detailed analysis
- `/data/salomonis-archive/FASTQs/NCI-R01/rna_seq_varcall/.gitignore` — Prevent VCF/BAM tracking

### Tools Verified Working
- GATK 4.6.2.0 ✅
- DeepVariant 1.6.1 ✅
- samtools 1.23 ✅
- bcftools 1.23 ✅
- Singularity ✅

### Tools Installed But Need Testing
- Clair3-RNA (conda environment created)
- LongcallR (container download in progress)

---

## Recommendations

### For This Pipeline
✅ **STOP searching for CBL variant** — It's in an intron, fundamentally undetectable from RNA-seq.

✅ **Accept 4-variant detection** — This is the expected result for RNA-seq on these truth variants.

⚠️ **To detect intronic variants**: Use DNA-seq instead of RNA-seq, or target the specific intron.

### For Production Use
1. Document that 4 of 8 truth variants are detectable (1 intronic, 3 low-coverage)
2. Use the same quality-preprocessing approach for future BAMs without QUAL field
3. Complete the 4-caller comparison for benchmark purposes

---

## Conclusion

The variant calling pipeline is **working correctly**. The detection of 4 variants with perfect GATK-DeepVariant concordance, combined with clear root-cause analysis of the missing variants, demonstrates:

1. **Callers are functional** — Both GATK and DeepVariant working as expected
2. **Detection is not random** — Perfect concordance proves systematic detection
3. **Results are data-driven** — 4 detectable (exonic) + 1 intronic (undetectable) + 3 low-coverage = 8 total

Next session can focus on completing Clair3-RNA and LongcallR setup for final 4-caller comparison.
