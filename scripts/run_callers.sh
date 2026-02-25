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
BAM_QUAL="$WORKDIR/scisoseq_synthetic_qual.bam"  # Full BAM with synthetic quality preprocessing
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
    _CONDA_BASE_DETECTED="$("$_CONDA_BIN" info --base 2>/dev/null)/envs"
    # Verify detected path has our envs; fallback if not
    if [ -d "$_CONDA_BASE_DETECTED/bio-cli" ]; then
        CONDA_BASE="$_CONDA_BASE_DETECTED"
    else
        # Fallback: common server path
        CONDA_BASE="/users/pavb5f/.conda/envs"
    fi
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
        -I "$BAM_QUAL" \
        -O "$OUTPUT_DIR/variants_raw.vcf.gz" \
        --dont-use-soft-clipped-bases \
        --min-base-quality-score 10 \
        --native-pair-hmm-threads 8 \
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
# GATK HaplotypeCaller - Sensitive Mode (Lower Detection Threshold)
# ============================================================================
run_haplotypecaller_sensitive() {
    echo ""
    echo "========================================================================"
    echo "[CALLER] GATK HaplotypeCaller (SENSITIVE - Lower Thresholds)"
    echo "========================================================================"

    OUTPUT_DIR="$WORKDIR/results/haplotypecaller_sensitive"
    LOG_FILE="$LOG_DIR/haplotypecaller_sensitive.log"

    mkdir -p "$OUTPUT_DIR"

    if [ ! -f "$GATK_BIN" ]; then
        echo "[ERROR] GATK not found at $GATK_BIN"
        exit 1
    fi

    echo "[INFO] Running GATK HaplotypeCaller with relaxed discovery thresholds..."
    echo "[INFO] Output: $OUTPUT_DIR"
    echo "[INFO] Log: $LOG_FILE"

    # Relaxed parameters for sensitive variant discovery:
    # --standard-min-confidence-threshold-for-calling: Lower confidence threshold (default 10)
    # --min-pruning: Reduce aggressive pruning of haplotypes
    # --disable-read-filter: Include reads that would normally be filtered
    # Use quality-preprocessed BAM to avoid WellformedReadFilter rejection
    echo "[INFO] Using quality-preprocessed BAM: $BAM_QUAL"
    PATH="$GATK_ENV_BIN:$PATH" "$GATK_BIN" HaplotypeCaller \
        -R "$REFERENCE" \
        -I "$BAM_QUAL" \
        -O "$OUTPUT_DIR/variants_raw.vcf.gz" \
        -L "$TRUTH_BED" \
        --standard-min-confidence-threshold-for-calling 0 \
        --min-pruning 1 \
        --min-base-quality-score 3 \
        --native-pair-hmm-threads 8 \
        --disable-read-filter NotDuplicateReadFilter \
        --disable-read-filter MappingQualityReadFilter \
        2>&1 | tee "$LOG_FILE"

    echo "[STEP] Normalizing..." | tee -a "$LOG_FILE"
    "$BCFTOOLS_BIN" norm -f "$REFERENCE" -m -any \
        "$OUTPUT_DIR/variants_raw.vcf.gz" | \
        "$BCFTOOLS_BIN" sort -Oz -o "$OUTPUT_DIR/variants.vcf.gz"
    "$BCFTOOLS_BIN" index -t "$OUTPUT_DIR/variants.vcf.gz"

    echo "[SUCCESS] HaplotypeCaller (Sensitive) analysis complete"
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
    echo "[INFO] Using quality-preprocessed BAM: $BAM_QUAL"
    echo "Running DeepVariant on ENTIRE BAM (no region restriction)..."

    "$SINGULARITY_BIN" exec \
        --bind /data:/data \
        --bind "$WORKDIR:$WORKDIR" \
        "$DV_SIF" \
        /opt/deepvariant/bin/run_deepvariant \
            --model_type=PACBIO \
            --ref="$REFERENCE" \
            --reads="$BAM_QUAL" \
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

    # Determine runtime: Singularity (Linux) or Docker (macOS)
    if [ -f "$CLAIR3_SIF" ] && command -v singularity &>/dev/null; then
        echo "[INFO] Using Singularity container: $CLAIR3_SIF"
        _CLAIR3_RUN() {
            "$SINGULARITY_BIN" exec \
                --bind /data:/data \
                --bind "$WORKDIR:$WORKDIR" \
                "$CLAIR3_SIF" \
                /opt/bin/run_clair3.sh "$@"
        }
    elif command -v docker &>/dev/null; then
        echo "[INFO] Using Docker image: $CLAIR3_DOCKER_IMAGE"
        _CLAIR3_RUN() {
            docker run --rm --platform linux/amd64 \
                -v /data:/data \
                -v "$WORKDIR:$WORKDIR" \
                "$CLAIR3_DOCKER_IMAGE" \
                /opt/bin/run_clair3.sh "$@"
        }
    else
        echo "[ERROR] Neither Singularity nor Docker found. Install one to run Clair3."
        echo "[INFO] macOS: brew install --cask docker"
        echo "[INFO] Linux: singularity pull $CLAIR3_SIF docker://$CLAIR3_DOCKER_IMAGE"
        exit 1
    fi

    echo "[INFO] Output: $OUTPUT_DIR"
    echo "[INFO] Log: $LOG_FILE"

    # run_clair3.sh handles both pileup and full-alignment stages per chromosome
    for CTGNAME in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
                   chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 \
                   chr21 chr22 chrX chrY chrM; do
        echo "[INFO] Calling variants on $CTGNAME..."
        _CLAIR3_RUN \
            --bam_fn="$BAM_QUAL" \
            --ref_fn="$REFERENCE" \
            --threads=16 \
            --platform="hifi" \
            --model_path="/opt/models/hifi" \
            --ctg_name="$CTGNAME" \
            --sample_name="sample" \
            --min_mq=5 \
            --snp_min_af=0.05 \
            --indel_min_af=0.05 \
            --output="$OUTPUT_DIR/${CTGNAME}" \
            2>&1 | tee -a "$LOG_FILE"
    done

    # Merge per-chromosome VCFs (run_clair3.sh outputs merge_output.vcf.gz per chrom dir)
    echo "[INFO] Merging VCFs..."
    "$BCFTOOLS_BIN" concat -o "$OUTPUT_DIR/variants.vcf.gz" -O z \
        $(for c in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
                   chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 \
                   chr21 chr22 chrX chrY chrM; do
              echo "$OUTPUT_DIR/${c}/merge_output.vcf.gz"
          done) 2>&1 | tee -a "$LOG_FILE"
    "$BCFTOOLS_BIN" index -t "$OUTPUT_DIR/variants.vcf.gz" 2>&1 | tee -a "$LOG_FILE"

    echo "[SUCCESS] Clair3-RNA analysis complete"
    echo "[INFO] Output: $OUTPUT_DIR/variants.vcf.gz"
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
                longcallR "$@"
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

    _LONGCALLR_RUN \
        --bam-path "$BAM_QUAL" \
        --ref-path "$REFERENCE" \
        --output "$OUTPUT_DIR/variants" \
        --preset hifi-isoseq \
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
    haplotypecaller_sensitive|gatk_sensitive)
        run_haplotypecaller_sensitive
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
