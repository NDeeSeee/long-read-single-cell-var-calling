#!/bin/bash
#
# Run all variant callers for comparison
# Usage: ./run_callers.sh [caller_name|all]
# Example: ./run_callers.sh all
#          ./run_callers.sh haplotypecaller
#          ./run_callers.sh deepvariant
#          ./run_callers.sh clair3_rna
#          ./run_callers.sh longcallr
#

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Working directory: all data, results, models live here
WORKDIR="/data/salomonis-archive/FASTQs/NCI-R01/rna_seq_varcall"

BAM="/data/salomonis-archive/BAMs/Grimes/scRNA-Seq/KINNEX/5801-diagnosis/scisoseq.mapped.bam"
REFERENCE="$WORKDIR/reference/genome.fa"
TRUTH_BED="$WORKDIR/truth_set/truth_regions.bed"

LOG_DIR="$WORKDIR/logs"
mkdir -p "$LOG_DIR"

CALLER="${1:-all}"

echo "[INFO] Starting variant caller comparison pipeline"
echo "[INFO] BAM: $BAM"
echo "[INFO] Reference: $REFERENCE"
echo "[INFO] Target caller: $CALLER"

# ── Tool paths ─────────────────────────────────────────────────────────────
# Detect conda envs directory dynamically (works on Linux server and macOS)
_CONDA_BIN="${CONDA_EXE:-$(command -v conda 2>/dev/null)}"
if [ -n "$_CONDA_BIN" ]; then
    CONDA_BASE="$("$_CONDA_BIN" info --base 2>/dev/null)/envs"
else
    # Fallback: common server path
    CONDA_BASE="/users/pavb5f/.conda/envs"
fi

GATK_BIN="$CONDA_BASE/gatk-env/bin/gatk"
GATK_ENV_BIN="$CONDA_BASE/gatk-env/bin"  # GATK wrapper needs python on PATH

# Clair3: prefer Singularity on Linux; fall back to Docker on macOS
CLAIR3_DOCKER_IMAGE="hkubal/clair3:latest"
CLAIR3_SIF="$WORKDIR/containers/clair3_latest.sif"

# LongcallR: original C++ tool is not on bioconda; use BioContainers image
# (longcallR_nn in conda env is a different deep-learning tool requiring feature extraction)
LONGCALLR_DOCKER_IMAGE="quay.io/biocontainers/longcallr:1.12.0--py312h67e1f27_0"
LONGCALLR_SIF="$WORKDIR/containers/longcallr_1.12.0.sif"

BGZIP_BIN="$CONDA_BASE/bio-cli/bin/bgzip"
BCFTOOLS_BIN="$CONDA_BASE/bio-cli/bin/bcftools"
SINGULARITY_BIN="$CONDA_BASE/bio-cli/bin/singularity"
DV_SIF="$WORKDIR/containers/deepvariant_1.6.1.sif"
CLAIR3_MODEL="$WORKDIR/models/clair3_rna"

# ============================================================================
# GATK HaplotypeCaller
# ============================================================================
run_haplotypecaller() {
    echo ""
    echo "========================================================================"
    echo "[CALLER] GATK HaplotypeCaller"
    echo "========================================================================"

    OUTPUT_DIR="$WORKDIR/results/haplotypecaller"
    LOG_FILE="$LOG_DIR/haplotypecaller.log"

    mkdir -p "$OUTPUT_DIR"

    if [ ! -f "$GATK_BIN" ]; then
        echo "[ERROR] GATK not found at $GATK_BIN"
        echo "[INFO] Create env: conda create -n gatk-env -c conda-forge -c bioconda openjdk=17 gatk4=4.6.2.0 -y"
        exit 1
    fi

    echo "[INFO] Running GATK HaplotypeCaller..."
    echo "[INFO] Output: $OUTPUT_DIR"
    echo "[INFO] Log: $LOG_FILE"

    # Prepend gatk-env bin so the gatk wrapper finds 'python'
    PATH="$GATK_ENV_BIN:$PATH" "$GATK_BIN" HaplotypeCaller \
        -R "$REFERENCE" \
        -I "$BAM" \
        -O "$OUTPUT_DIR/variants_raw.vcf.gz" \
        -L "$TRUTH_BED" \
        --dont-use-soft-clipped-bases \
        --min-base-quality-score 10 \
        --native-pair-hmm-threads 8 \
        --sample-name 5801-diagnosis \
        2>&1 | tee "$LOG_FILE"

    echo "[STEP] Normalizing..." | tee -a "$LOG_FILE"
    "$BCFTOOLS_BIN" norm -f "$REFERENCE" -m -any \
        "$OUTPUT_DIR/variants_raw.vcf.gz" | \
        "$BCFTOOLS_BIN" sort -Oz -o "$OUTPUT_DIR/variants.vcf.gz"
    "$BCFTOOLS_BIN" index -t "$OUTPUT_DIR/variants.vcf.gz"

    echo "[SUCCESS] HaplotypeCaller analysis complete"
    echo "[INFO] Output: $OUTPUT_DIR/variants.vcf.gz"
}

# ============================================================================
# DeepVariant (PacBio model via Singularity)
# ============================================================================
run_deepvariant() {
    echo ""
    echo "========================================================================"
    echo "[CALLER] DeepVariant (PACBIO model via Singularity)"
    echo "========================================================================"

    OUTPUT_DIR="$WORKDIR/results/deepvariant"
    LOG_FILE="$LOG_DIR/deepvariant.log"

    mkdir -p "$OUTPUT_DIR"

    # If container is still at project root, move it into containers/
    if [ ! -f "$DV_SIF" ]; then
        ROOT_SIF="$PROJECT_DIR/deepvariant_1.6.1.sif"
        if [ -f "$ROOT_SIF" ]; then
            echo "[INFO] Moving container from project root to containers/ directory..."
            mkdir -p "$(dirname "$DV_SIF")"
            mv "$ROOT_SIF" "$DV_SIF"
        else
            echo "[INFO] DeepVariant container not found. Pulling from Docker Hub..."
            mkdir -p "$(dirname "$DV_SIF")"
            "$SINGULARITY_BIN" pull --dir "$(dirname "$DV_SIF")" \
                docker://google/deepvariant:1.6.1
        fi
    fi

    echo "[INFO] Running DeepVariant..."
    echo "[INFO] Container: $DV_SIF"
    echo "[INFO] Output: $OUTPUT_DIR"
    echo "[INFO] Log: $LOG_FILE"

    "$SINGULARITY_BIN" exec \
        --bind /data:/data \
        --bind "$WORKDIR:$WORKDIR" \
        "$DV_SIF" \
        /opt/deepvariant/bin/run_deepvariant \
            --model_type=PACBIO \
            --ref="$REFERENCE" \
            --reads="$BAM" \
            --regions="$TRUTH_BED" \
            --output_vcf="$OUTPUT_DIR/variants.vcf.gz" \
            --num_shards=8 \
            2>&1 | tee "$LOG_FILE"

    "$BCFTOOLS_BIN" index -t "$OUTPUT_DIR/variants.vcf.gz"

    echo "[SUCCESS] DeepVariant analysis complete"
    echo "[INFO] Output: $OUTPUT_DIR/variants.vcf.gz"
}

# ============================================================================
# Clair3-RNA (PacBio HiFi model)
# Singularity on Linux servers; Docker on macOS
# ============================================================================
run_clair3_rna() {
    echo ""
    echo "========================================================================"
    echo "[CALLER] Clair3 (PacBio HiFi / hifi platform)"
    echo "========================================================================"

    OUTPUT_DIR="$WORKDIR/results/clair3_rna"
    LOG_FILE="$LOG_DIR/clair3_rna.log"

    mkdir -p "$OUTPUT_DIR"

    if [ ! -d "$CLAIR3_MODEL" ]; then
        echo "[ERROR] Clair3 model directory not found at $CLAIR3_MODEL"
        echo "[INFO] Run: bash scripts/setup_envs.sh clair3"
        echo "[INFO] Or download from: https://github.com/HKU-BAL/Clair3#pre-trained-models"
        exit 1
    fi

    # Determine runtime: Singularity (Linux) or Docker (macOS)
    if [ -f "$CLAIR3_SIF" ] && command -v singularity &>/dev/null; then
        echo "[INFO] Using Singularity container: $CLAIR3_SIF"
        _CLAIR3_RUN() {
            "$SINGULARITY_BIN" exec \
                --bind /data:/data \
                --bind "$WORKDIR:$WORKDIR" \
                "$CLAIR3_SIF" \
                /usr/local/bin/run_clair3.sh "$@"
        }
    elif command -v docker &>/dev/null; then
        echo "[INFO] Using Docker image: $CLAIR3_DOCKER_IMAGE"
        _CLAIR3_RUN() {
            docker run --rm --platform linux/amd64 \
                -v /data:/data \
                -v "$WORKDIR:$WORKDIR" \
                "$CLAIR3_DOCKER_IMAGE" \
                run_clair3.sh "$@"
        }
    else
        echo "[ERROR] Neither Singularity nor Docker found. Install one to run Clair3."
        echo "[INFO] macOS: brew install --cask docker"
        echo "[INFO] Linux: singularity pull $CLAIR3_SIF docker://$CLAIR3_DOCKER_IMAGE"
        exit 1
    fi

    echo "[INFO] Output: $OUTPUT_DIR"
    echo "[INFO] Log: $LOG_FILE"

    _CLAIR3_RUN \
        --bam_fn="$BAM" \
        --ref_fn="$REFERENCE" \
        --threads=16 \
        --platform="hifi" \
        --model_path="$CLAIR3_MODEL" \
        --output="$OUTPUT_DIR" \
        --bed_fn="$TRUTH_BED" \
        --include_all_ctgs \
        --no_phasing_for_fa \
        --haploid_sensitive \
        2>&1 | tee "$LOG_FILE"

    echo "[SUCCESS] Clair3-RNA analysis complete"
    echo "[INFO] Output: $OUTPUT_DIR/merge_output.vcf.gz"
}

# ============================================================================
# LongcallR (BioContainers: Singularity on Linux, Docker on macOS)
# ============================================================================
run_longcallr() {
    echo ""
    echo "========================================================================"
    echo "[CALLER] LongcallR"
    echo "========================================================================"

    OUTPUT_DIR="$WORKDIR/results/longcallr"
    LOG_FILE="$LOG_DIR/longcallr.log"

    mkdir -p "$OUTPUT_DIR"

    # Determine runtime: Singularity (Linux) or Docker (macOS)
    if [ -f "$LONGCALLR_SIF" ] && command -v singularity &>/dev/null; then
        echo "[INFO] Using Singularity container: $LONGCALLR_SIF"
        _LONGCALLR_RUN() {
            "$SINGULARITY_BIN" exec \
                --bind /data:/data \
                --bind "$WORKDIR:$WORKDIR" \
                "$LONGCALLR_SIF" \
                longcallr "$@"
        }
    elif command -v docker &>/dev/null; then
        echo "[INFO] Using Docker image: $LONGCALLR_DOCKER_IMAGE"
        if ! docker image inspect "$LONGCALLR_DOCKER_IMAGE" &>/dev/null; then
            echo "[INFO] Pulling Docker image..."
            docker pull --platform linux/amd64 "$LONGCALLR_DOCKER_IMAGE"
        fi
        _LONGCALLR_RUN() {
            docker run --rm --platform linux/amd64 \
                -v /data:/data \
                -v "$WORKDIR:$WORKDIR" \
                "$LONGCALLR_DOCKER_IMAGE" \
                longcallR "$@"
        }
    else
        echo "[ERROR] Neither Singularity SIF nor Docker found."
        echo "[INFO] Linux: singularity pull $LONGCALLR_SIF docker://$LONGCALLR_DOCKER_IMAGE"
        echo "[INFO] macOS: start Docker Desktop, then re-run."
        exit 1
    fi

    echo "[INFO] Output: $OUTPUT_DIR"
    echo "[INFO] Log: $LOG_FILE"

    # Chromosomes containing truth variants
    CONTIGS="chr11 chr17 chr18 chr19 chr20 chr21 chrX"

    _LONGCALLR_RUN \
        --bam-path "$BAM" \
        --ref-path "$REFERENCE" \
        --output "$OUTPUT_DIR/variants" \
        --preset hifi-isoseq \
        --contigs $CONTIGS \
        --threads 16 \
        --min-allele-freq 0.05 \
        --low-allele-frac-cutoff 0.03 \
        --min-mapq 1 \
        --no-bam-output \
        2>&1 | tee "$LOG_FILE"

    # longcallR writes output as prefix.vcf — compress and index
    if [ -f "$OUTPUT_DIR/variants.vcf" ]; then
        "$BGZIP_BIN" "$OUTPUT_DIR/variants.vcf"
    fi
    "$BCFTOOLS_BIN" index -t "$OUTPUT_DIR/variants.vcf.gz"

    echo "[SUCCESS] LongcallR analysis complete"
    echo "[INFO] Output: $OUTPUT_DIR/variants.vcf.gz"
}

# ============================================================================
# Main entry point
# ============================================================================
case "$CALLER" in
    haplotypecaller|gatk)
        run_haplotypecaller
        ;;
    deepvariant)
        run_deepvariant
        ;;
    clair3_rna)
        run_clair3_rna
        ;;
    longcallr)
        run_longcallr
        ;;
    all)
        run_haplotypecaller
        run_deepvariant
        run_clair3_rna
        run_longcallr
        ;;
    *)
        echo "[ERROR] Unknown caller: $CALLER"
        echo "Usage: $0 [haplotypecaller|deepvariant|clair3_rna|longcallr|all]"
        exit 1
        ;;
esac

echo ""
echo "========================================================================"
echo "[COMPLETE] Caller execution finished: $CALLER"
echo "========================================================================"
echo "[INFO] Results available in: $WORKDIR/results/"
echo ""
