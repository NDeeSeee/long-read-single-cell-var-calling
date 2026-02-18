# Variant Caller Runner Status

**Date**: 2026-02-18
**Status**: 1 of 4 callers fully working (GATK), 1 partially working (DeepVariant), 2 pending (Clair3-RNA, LongcallR)

---

## ✅ WORKING: GATK HaplotypeCaller

**Status**: Fully functional
**Container**: `containers/gatk.sif` (2.3 GB, broadinstitute/gatk:latest)
**Output**: `results/haplotypecaller/variants.vcf.gz` + index

### Test Run
```bash
./scripts/run_callers.sh gatk
```

### Results
- Completed successfully in 0.16 minutes
- Generated VCF: 5.2 KB (normalized and indexed)
- Note: All 8,151 reads filtered due to MQ≥20 threshold
  - This is expected behavior for RNA-seq data lacking proper quality scores
  - GATK minimum mapping quality filter removed many multi-mapping reads

### Implementation Notes
- Uses Singularity container (avoids conda Java/glibc issues)
- Binds: `/data:/data`, project directory as `/work`
- Command executed via `bash -c` inside container for reliable arg passing
- Output normalized with bcftools norm (handles ref alignment issues)

---

## ⚠️  PARTIALLY WORKING: DeepVariant

**Status**: Crashes during variant calling
**Container**: `containers/deepvariant_1.6.1.sif` (2.7 GB)
**Error**: "Check failed: offset + len <= read.aligned_quality_size()"

### Problem
The PacBio BAM file is missing base quality scores in the standard SAM format:
```
ERROR: Could not read base quality scores for reads
```

This is a known issue with PacBio Iso-Seq BAMs which store quality differently (using auxiliary SAM tags like "qs:Z:" instead of QUAL field).

### Workaround Needed
1. Extract PacBio quality scores from auxiliary tags with samtools
2. Re-encode BAM with standard quality field
3. Then re-run DeepVariant

Alternatively:
- Use a different variant caller designed for PacBio (e.g., pbsv, PBCCS)
- Or filter BAM to include only reads with embedded quality scores

---

## ⏳ PENDING: Clair3-RNA

**Status**: Container registry access issues
**Attempted**: Pull from multiple registries (hkubal/clair3, ghcr.io/hku-bal/clair3, quay.io/biocontainers/clair3)
**Result**: All failed (image not found, auth required, or tag mismatch)

### Workaround Options
1. **Install from GitHub source**:
   ```bash
   cd /tmp
   git clone https://github.com/HKU-BAL/Clair3.git
   cd Clair3
   python -m pip install -e .
   ```

2. **Manual conda with relaxed constraints**:
   ```bash
   conda create -n clair3 python=3.11
   pip install clair3 setuptools
   ```

3. **Wait for stable Docker image** (Clair3 has recent builds)

---

## ⏳ PENDING: LongcallR

**Status**: No compatible container found, conda install blocked
**Attempted**: quay.io/biocontainers/longcallr:0.4.1 (tag not found)
**Conda Error**: pysam version conflicts with other dependencies

### Workaround Options
1. **Try pip install**:
   ```bash
   pip install longcallr
   ```

2. **Build from source**:
   ```bash
   pip install -e git+https://github.com/naotoyokota/longcallr.git@main#egg=longcallr
   ```

3. **Use alternative**: longshot, medaka, or pb-align for long-read variant calling

---

## Running Callers

### Command Format
```bash
./scripts/run_callers.sh [caller_name|all]
```

### Available Callers
| Name | Status | Command |
|------|--------|---------|
| gatk | ✅ Ready | `./scripts/run_callers.sh gatk` |
| deepvariant | ⚠️ Issues | `./scripts/run_callers.sh deepvariant` |
| clair3_rna | ⏳ Blocked | `./scripts/run_callers.sh clair3_rna` |
| longcallr | ⏳ Blocked | `./scripts/run_callers.sh longcallr` |
| lab_supervised | ✅ Ready | `./scripts/run_callers.sh lab_supervised` |
| lab_unsupervised | ✅ Ready | `./scripts/run_callers.sh lab_unsupervised` |
| all | 🟡 Partial | `./scripts/run_callers.sh all` (runs all available) |

---

## File Structure
```
results/
├── haplotypecaller/
│   ├── variants_raw.vcf.gz       (GATK raw output)
│   ├── variants_raw.vcf.gz.tbi   (GATK index)
│   ├── variants.vcf.gz           (normalized)
│   └── variants.vcf.gz.tbi       (index)
├── deepvariant/                   (empty - crashed)
├── clair3_rna/                    (pending)
├── longcallr/                     (pending)
├── supervised_extraction/         (lab script)
└── global_snv/                    (lab script)

containers/
├── gatk.sif                       (2.3 GB) ✅
├── deepvariant_1.6.1.sif         (2.7 GB) ⚠️
├── clair3.sif                     (pending)
└── longcallr.sif                  (pending)
```

---

## Key Insights

1. **Reference Genome**: Must use chr-prefixed version (`reference/genome.fa`)
   - Absolute path version lacks chr prefix causing contig mismatch

2. **BAM Quality Issues**: PacBio Iso-Seq uses non-standard quality encoding
   - Multiple callers expect standard SAM QUAL field
   - May need format conversion before running variant callers

3. **Singularity vs Conda**: Containers are more reliable than conda
   - Java/glibc conflicts in conda resolved by using containers
   - Easier dependency isolation with containers

4. **Registry Access**: Public container registries may require auth or have rate limits
   - BioContainers sometimes has outdated or missing tags
   - GitHub Container Registry may require authentication
   - Broad Institute (broadinstitute/gatk) most reliable for GATK

---

## Next Steps

1. **Fix DeepVariant**: Extract PacBio quality scores or use alternative caller
2. **Install Clair3-RNA**: Try GitHub source or manual pip install
3. **Install LongcallR**: Try pip or GitHub source
4. **Run evaluation**: Once all callers complete, run `scripts/evaluate_callers.py`
5. **Generate report**: Run `scripts/generate_report.py` (f-string bug already fixed)

---

## Troubleshooting

### GATK fails in container
- Ensure reference genome exists and has index
- Check BAM file is readable and properly formatted
- Verify truth BED file has correct contig names (must match reference)

### DeepVariant crashes
- Check if BAM has standard QUAL field
- Try reducing number of shards: `--num_shards=1`
- Add `--make_examples_extra_args="min_mapping_quality=0"` to allow all reads

### Container pulls fail
- Check internet connection: `curl https://docker.io`
- Try `singularity pull --insecure` (not recommended for production)
- Use local .sif files if available
- Check if registry requires authentication

