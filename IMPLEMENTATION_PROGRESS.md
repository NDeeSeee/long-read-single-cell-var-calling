# Implementation Progress - Variant Caller Setup

**Date**: 2026-02-18
**Status**: Partial - Infrastructure ready, Tools pending installation/fix

## Completed ✅

### 1. DeepVariant Container Setup
- **Status**: ✅ COMPLETE
- **Action**: Copied `deepvariant_1.6.1.sif` from project root to `containers/` directory
- **Path**: `/data/salomonis-archive/FASTQs/NCI-R01/rna_seq_varcall/containers/deepvariant_1.6.1.sif`
- **Model**: Changed from WGS to PACBIO model in run_callers.sh (appropriate for PacBio long reads)
- **Ready to Run**: ✓ YES

### 2. Code Updates
- **f-string Bug Fix**: ✅ FIXED
  - File: `scripts/generate_report.py`
  - Issue: Line 279 had `{samtools_version}` but variable was undefined
  - Solution: Added `get_tool_version()` function and defined `samtools_version` variable
  - Now extracts samtools version dynamically from system

- **run_callers.sh Updates**: ✅ COMPLETE
  - Added `run_gatk()` with isolated environment support
  - Updated `run_deepvariant()` with correct container path and model type
  - Added `run_clair3_rna()` function stub
  - Added `run_longcallr()` function stub
  - Updated main case statement to include all 4 callers
  - Changed GATK output directory: `results/gatk/` → `results/haplotypecaller/`

### 3. Directory Structure
- **Created**: `containers/` with DeepVariant container
- **Ready**: `results/` directory exists for all output

---

## Pending ⏳

### 1. GATK HaplotypeCaller
- **Status**: ⚠️ BLOCKED - Java/libc library mismatch
- **Problem**: Both `gatk-env` and `gatk-env-v2` conda environments fail with:
  ```
  java: symbol lookup error: java: undefined symbol: JLI_StringDup
  ```
- **Root Cause**: Java binary in conda environment has glibc incompatibility with system
- **Environment Created**:
  - `gatk-env`: OpenJDK 17 + GATK 4.6.2.0
  - `gatk-env-v2`: OpenJDK 21 + GATK 4.6.2.0
  - Both fail with the same symbol error
- **Workarounds Needed**:
  - Option A: Use GATK Docker container instead of conda
  - Option B: Fix glibc/libc compatibility (requires system libs investigation)
  - Option C: Use alternative Java distribution (e.g., Eclipse Temurin)
  - Option D: Use haplotype-caller specific container from BioContainers

### 2. Clair3-RNA
- **Status**: ❌ BLOCKED - Conda dependency conflicts
- **Problem**: Cannot resolve dependencies due to:
  - `libcurl` version conflicts (8.13+/8.14.1+ needed)
  - `libnghttp2` version chain conflicts
  - `c-ares` version incompatibilities
  - Incompatible `whatshap` dependencies for older Clair3 versions
- **Environment**: `clair3-rna` created with Python 3.10, but cannot install Clair3 package
- **Workarounds Needed**:
  - Option A: Use Clair3 Singularity container from BioContainers
  - Option B: Install from GitHub source instead of conda
  - Option C: Use micromamba/conda-lock for environment resolution

### 3. LongcallR
- **Status**: ❌ BLOCKED - Conda dependency conflicts
- **Problem**: Cannot resolve dependencies:
  - `pysam >=0.23.0,<0.24.0a0` required but conflicts with other packages
  - `libdeflate` version chains incompatible
- **Environment**: `longcallr` not created due to solver conflicts
- **Workarounds Needed**:
  - Option A: Use LongcallR Singularity container from BioContainers
  - Option B: Install from PyPI or GitHub source
  - Option C: Create minimal isolated environment with specific versions

---

## Next Steps Recommended

### Priority 1: Fix GATK
1. Try Eclipse Temurin Java: `conda create -n gatk-temurin -c conda-forge temurin-jdk=17`
2. Or pull GATK container: `singularity pull gatk.sif docker://broadinstitute/gatk:latest`
3. Test with: `containers/deepvariant_1.6.1.sif` as model for containerized approach

### Priority 2: Install Clair3-RNA & LongcallR via Containers
1. Search BioContainers for available versions:
   - `singularity pull clair3.sif docker://quay.io/biocontainers/clair3:latest`
   - `singularity pull longcallr.sif docker://quay.io/biocontainers/longcallr:latest`
2. Create wrapper scripts in `scripts/` directory to invoke containers
3. Update `run_callers.sh` to use container paths

### Priority 3: Test Pipeline
1. Verify reference genome symlink: `reference/genome.fa`
2. Verify truth set: `truth_set/truth_regions.bed`, `truth_set/test_variants.txt`
3. Verify BAM file: `/data/salomonis-archive/BAMs/Grimes/scRNA-Seq/KINNEX/5801-diagnosis/scisoseq.mapped.bam`
4. Run DeepVariant (ready to go):
   ```bash
   ./scripts/run_callers.sh deepvariant
   ```

---

## System Information
- **Shell Issue**: zsh in bio-cli has warnings about libtinfow.so.6 (non-fatal)
- **Singularity**: ✅ Available at `/users/pavb5f/.conda/envs/bio-cli/bin/singularity`
- **Conda**: ✅ Functional for environment creation
- **samtools 1.23**: ✅ Working in bio-cli environment
- **bcftools**: ✅ Available in bio-cli environment
- **bgzip**: ✅ Available in bio-cli environment

---

## File Changes Summary

| File | Change | Status |
|------|--------|--------|
| `scripts/generate_report.py` | Fix f-string, add samtools_version | ✅ Complete |
| `scripts/run_callers.sh` | Add Clair3/LongcallR, fix paths | ✅ Complete |
| `containers/deepvariant_1.6.1.sif` | Copied from root | ✅ Complete |
| GATK environment | Two conda envs created but broken | ⚠️ Pending fix |
| Clair3-RNA environment | Python 3.10 env created, install blocked | ⚠️ Pending fix |
| LongcallR environment | Not created due to conflicts | ⚠️ Pending creation |

---

## Testing Checklist

- [ ] DeepVariant runs successfully with PACBIO model
- [ ] GATK HaplotypeCaller works (pending Java fix)
- [ ] Clair3-RNA runs successfully (pending install)
- [ ] LongcallR runs successfully (pending install)
- [ ] All 4 VCFs generated in results/ directories
- [ ] evaluate_callers.py completes successfully
- [ ] generate_report.py produces markdown report

