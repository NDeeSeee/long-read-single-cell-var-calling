# Running RNA-seq Variant Calling Pipeline on MacBook

## Quick Setup (5 minutes)

```bash
# 1. Clone the repository
git clone /path/to/rna_seq_varcall
cd rna_seq_varcall

# 2. Create conda environment
conda create -n variant-callers -c bioconda \
  samtools bcftools gatk4=4.6.2.0 \
  -y

conda activate variant-callers

# 3. Download test data and reference (if not already available)
# Copy from server:
# - BAM: /data/salomonis-archive/BAMs/Grimes/scRNA-Seq/KINNEX/5801-diagnosis/scisoseq.mapped.bam
# - Reference: /data/salomonis-archive/genomes/hg38/genome.fa

# 4. Run all callers
python scripts/evaluate_callers.py
```

## Detailed Setup

### Prerequisites

**On MacBook with Homebrew**:
```bash
# Install conda (if not already installed)
brew install miniforge

# Restart terminal and verify
conda --version
```

### Step 1: Clone Repository

```bash
# Clone from git
git clone <repo-url> rna_seq_varcall
cd rna_seq_varcall

# Or copy from server
scp -r <user>@<server>:/data/salomonis-archive/FASTQs/NCI-R01/rna_seq_varcall .
```

### Step 2: Create Environment

```bash
# Create fresh environment
conda create -n variant-callers python=3.10 -y
conda activate variant-callers

# Install core tools
conda install -c bioconda -c conda-forge \
  samtools=1.21 \
  bcftools=1.23 \
  gatk4=4.6.2.0 \
  -y

# Verify installations
samtools --version
bcftools --version
gatk --version
```

### Step 3: Prepare Data

Copy or download required files:

```bash
# BAM file (11GB)
# Option A: Copy from server
scp <user>@<server>:/data/salomonis-archive/BAMs/Grimes/scRNA-Seq/KINNEX/5801-diagnosis/scisoseq.mapped.bam \
  ~/data/scisoseq.mapped.bam

# Option B: Download if available online
# (contact lab for credentials/link)

# Reference genome (3GB)
scp <user>@<server>:/data/salomonis-archive/genomes/hg38/genome.fa \
  ~/data/genome.fa

# Index reference
samtools faidx ~/data/genome.fa

# Directory structure
mkdir -p ~/variant-calling-work
cd ~/variant-calling-work
cp -r /path/to/cloned/rna_seq_varcall/* .
```

### Step 4: Update Paths

Edit scripts to use local paths:

**In `scripts/run_callers.sh`**:
```bash
# Change these paths:
BAM_PATH="/Users/yourname/data/scisoseq.mapped.bam"
REFERENCE="/Users/yourname/data/genome.fa"
OUTPUT_DIR="./results"
```

**In `CLAUDE.md`**:
```bash
# Update BAM and reference paths to local locations
```

### Step 5: Run Variant Callers

#### Option A: Run Individual Callers

```bash
# Lab Supervised (1-2 hours)
python /path/to/altanalyze3/variant_extraction.py \
  --sample 5801-diagnosis \
  --bam ~/data/scisoseq.mapped.bam \
  --mutations truth_set/test_variants_formatted.txt \
  --reference ~/data/genome.fa \
  --output-dir results/supervised_extraction

# GATK HaplotypeCaller (30-60 min)
gatk HaplotypeCaller \
  -R ~/data/genome.fa \
  -I ~/data/scisoseq.mapped.bam \
  -O results/gatk/variants.vcf.gz \
  --dont-use-soft-clipped-bases \
  -L truth_set/truth_regions.bed \
  --native-pair-hmm-threads 8

# DeepVariant (requires Singularity/Docker)
# See DeepVariant section below
```

#### Option B: Run Custom Script

```bash
# Copy test data
cp ~/data/scisoseq.mapped.bam ./data/
cp ~/data/genome.fa ./data/

# Update paths in scripts
sed -i 's|/data/salomonis-archive|./data|g' scripts/*.py

# Run evaluation
python scripts/evaluate_callers.py
```

## Tool-Specific Setup

### GATK (Java)

**Problem**: Old Conda Java often has issues on MacBook
**Solution**: Use system Java

```bash
# Check system Java
java -version

# If needed, install via Homebrew
brew install openjdk@17

# Set JAVA_HOME
export JAVA_HOME=$(/usr/libexec/java_home -v 17)

# Verify GATK works
gatk --version
```

### DeepVariant (Requires Docker/Singularity)

**Option 1: Docker (Easier on MacBook)**

```bash
# Install Docker Desktop (if not already)
# Download from: https://www.docker.com/products/docker-desktop

# Pull DeepVariant image
docker pull google/deepvariant:1.6.1

# Run variant calling
docker run -it \
  -v ~/data:/data \
  -v $(pwd)/results:/results \
  google/deepvariant:1.6.1 \
  /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref=/data/genome.fa \
    --reads=/data/scisoseq.mapped.bam \
    --output_vcf=/results/deepvariant/variants.vcf.gz \
    --num_shards=4
```

**Option 2: Singularity (If available)**

```bash
# Install Singularity (requires Linux VM on Mac)
# Or use Lima: brew install lima

# Pull and run
singularity pull docker://google/deepvariant:1.6.1
singularity exec deepvariant_1.6.1.sif /opt/deepvariant/bin/run_deepvariant ...
```

### Clair3-RNA & LongcallR

**Note**: These have conda dependency issues. Try pip if conda fails:

```bash
# Create fresh environment
conda create -n clair3 python=3.10 -y
conda activate clair3

# Try conda first (fastest)
conda install -c bioconda clair3 -y

# If conda fails, try pip
pip install clair3

# Run
run_clair3.py \
  --bam_fn ~/data/scisoseq.mapped.bam \
  --ref_fn ~/data/genome.fa \
  --output_dir results/clair3 \
  --model_path ~/.clair3
```

## Running Tests

### Verify Installation

```bash
# Test all tools
python scripts/verify_setup.sh
```

### Test with Subset

Run on just chr21 (contains RUNX1 variant):

```bash
# Edit truth_set/truth_regions_nochr.bed to only include:
# chr21	34799431	34799433

# Run GATK on reduced regions
gatk HaplotypeCaller \
  -R ~/data/genome.fa \
  -I ~/data/scisoseq.mapped.bam \
  -O results/gatk/chr21_only.vcf.gz \
  -L truth_set/truth_regions_nochr.bed:chr21
```

## Evaluation & Results

```bash
# After running callers, evaluate
python scripts/evaluate_callers.py

# View results
cat results/comparison/summary.json
open results/comparison/COMPARISON_REPORT.md
```

## Troubleshooting MacBook-Specific Issues

### Issue: "Command not found" after conda install

**Solution**:
```bash
# Activate environment
conda activate variant-callers

# Verify PATH
which samtools
which gatk
```

### Issue: "Java not found" with GATK

**Solution**:
```bash
# Set JAVA_HOME
export JAVA_HOME=$(/usr/libexec/java_home)
export PATH="$JAVA_HOME/bin:$PATH"

# Add to ~/.zshrc or ~/.bash_profile
echo 'export JAVA_HOME=$(/usr/libexec/java_home)' >> ~/.zshrc
source ~/.zshrc
```

### Issue: "Permission denied" when reading BAM

**Solution**:
```bash
# Check permissions
ls -l ~/data/scisoseq.mapped.bam

# Fix if needed
chmod 644 ~/data/scisoseq.mapped.bam
```

### Issue: "Out of memory" during variant calling

**Solution**: Increase Java heap size
```bash
export _JAVA_OPTIONS="-Xmx16g"
gatk HaplotypeCaller ...
```

### Issue: Docker/Singularity not available

**Solution**: Skip DeepVariant or use lab's variant_extraction.py instead

## Recommended Workflow

For fastest results on MacBook:

1. **Run Lab Supervised** (1-2 hours)
   - Highest quality results
   - Designed for this data

2. **Run GATK** (30-60 min)
   - For comparison
   - Industry standard

3. **Evaluate** (10 min)
   ```bash
   python scripts/evaluate_callers.py
   ```

4. **Skip DeepVariant** (optional, slow)
   - Requires Docker setup
   - Not necessary for comparison

**Total time**: ~2.5 hours

## Output

Results will be in `./results/`:
- `supervised_extraction/` - Lab results
- `gatk/` - GATK VCFs
- `comparison/` - Performance metrics + report

View report:
```bash
open results/comparison/COMPARISON_REPORT.md
```

## Support

If you encounter issues:
1. Check `IMPLEMENTATION_STATUS.md` for known issues
2. Review `README.md` troubleshooting section
3. Check conda environment with `conda list`
4. View conda logs with `conda log`

---

**Last Updated**: 2026-02-18
**Status**: Ready for MacBook execution
