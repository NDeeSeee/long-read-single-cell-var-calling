#!/bin/bash
#
# Run all variant callers for comparison
# Usage: ./run_callers.sh [caller_name|all]
# Example: ./run_callers.sh all
#          ./run_callers.sh lab_supervised
#

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_DIR"

BAM="/data/salomonis-archive/BAMs/Grimes/scRNA-Seq/KINNEX/5801-diagnosis/scisoseq.mapped.bam"
REFERENCE="/data/salomonis-archive/genomes/hg38/genome.fa"
TRUTH_BED="$PROJECT_DIR/truth_set/truth_regions.bed"
TRUTH_VARIANTS="$PROJECT_DIR/truth_set/test_variants_formatted.txt"

LOG_DIR="$PROJECT_DIR/logs"
mkdir -p "$LOG_DIR"

CALLER="${1:-all}"

echo "[INFO] Starting variant caller comparison pipeline"
echo "[INFO] BAM: $BAM"
echo "[INFO] Reference: $REFERENCE"
echo "[INFO] Target caller: $CALLER"

# ============================================================================
# LAB SCRIPT 1: Supervised variant extraction
# ============================================================================
run_lab_supervised() {
    echo ""
    echo "========================================================================"
    echo "[CALLER] Lab Script - Supervised (variant_extraction.py)"
    echo "========================================================================"

    CALLER_SCRIPT="/data/salomonis-archive/LabFiles/Nathan/Revio/altanalyze3/altanalyze3/components/bam/variant_extraction.py"
    OUTPUT_DIR="$PROJECT_DIR/results/supervised_extraction"
    LOG_FILE="$LOG_DIR/lab_supervised.log"

    mkdir -p "$OUTPUT_DIR"

    echo "[INFO] Running variant_extraction.py..."
    echo "[INFO] Output: $OUTPUT_DIR"
    echo "[INFO] Log: $LOG_FILE"

    python "$CALLER_SCRIPT" \
        --sample 5801-diagnosis \
        --bam "$BAM" \
        --mutations "$TRUTH_VARIANTS" \
        --reference "$REFERENCE" \
        --output-dir "$OUTPUT_DIR" \
        2>&1 | tee "$LOG_FILE"

    echo "[SUCCESS] Lab supervised extraction complete"
    echo "[INFO] Converting TSV to VCF..."
    python "$SCRIPT_DIR/tsv_to_vcf.py" \
        "$OUTPUT_DIR/5801-diagnosis_complete_analysis.tsv" \
        "$OUTPUT_DIR/lab_supervised.vcf" 2>&1 | tee -a "$LOG_FILE"

    echo "[INFO] Compressing VCF..."
    python "$SCRIPT_DIR/compress_vcf.py" "$OUTPUT_DIR/lab_supervised.vcf"
}

# ============================================================================
# LAB SCRIPT 2: Unsupervised SNV discovery
# ============================================================================
run_lab_unsupervised() {
    echo ""
    echo "========================================================================"
    echo "[CALLER] Lab Script - Unsupervised (global_snv.py)"
    echo "========================================================================"

    CALLER_SCRIPT="/data/salomonis-archive/LabFiles/Nathan/Revio/altanalyze3/altanalyze3/components/bam/global_snv.py"
    OUTPUT_DIR="$PROJECT_DIR/results/global_snv"
    LOG_FILE="$LOG_DIR/lab_unsupervised.log"

    mkdir -p "$OUTPUT_DIR"

    echo "[INFO] Running global_snv.py (chr21 test)..."
    echo "[INFO] Output: $OUTPUT_DIR"
    echo "[INFO] Log: $LOG_FILE"

    python "$CALLER_SCRIPT" \
        "$BAM" \
        "$REFERENCE" \
        "$OUTPUT_DIR/chr21_output.txt" \
        --test_chromosome chr21 \
        --min_reads 10 \
        --min_percent 5.0 \
        2>&1 | tee "$LOG_FILE"

    echo "[SUCCESS] Lab unsupervised extraction complete"
    echo "[INFO] Chr21 test complete. To run full genome:"
    echo "[INFO]   python $CALLER_SCRIPT \\"
    echo "[INFO]     $BAM \\"
    echo "[INFO]     $REFERENCE \\"
    echo "[INFO]     $OUTPUT_DIR/genome_output.txt \\"
    echo "[INFO]     --min_reads 50 --min_percent 8.0"
}

# ============================================================================
# GATK HaplotypeCaller
# ============================================================================
run_gatk() {
    echo ""
    echo "========================================================================"
    echo "[CALLER] GATK HaplotypeCaller (RNA-seq mode)"
    echo "========================================================================"

    OUTPUT_DIR="$PROJECT_DIR/results/haplotypecaller"
    LOG_FILE="$LOG_DIR/gatk.log"

    mkdir -p "$OUTPUT_DIR"

    echo "[INFO] Running GATK HaplotypeCaller..."
    echo "[INFO] Output: $OUTPUT_DIR"
    echo "[INFO] Log: $LOG_FILE"

    # Try to use GATK from isolated environment
    GATK_ENV="/users/pavb5f/.conda/envs/gatk-env-v2/bin/gatk"

    if [ ! -x "$GATK_ENV" ]; then
        echo "[ERROR] GATK not found at $GATK_ENV"
        echo "[ERROR] Please create gatk-env-v2 conda environment first"
        return 1
    fi

    # Call variants
    echo "[STEP] HaplotypeCaller..." | tee "$LOG_FILE"
    $GATK_ENV HaplotypeCaller \
        -R "$REFERENCE" \
        -I "$BAM" \
        -O "$OUTPUT_DIR/variants_raw.vcf.gz" \
        --dont-use-soft-clipped-bases \
        --min-base-quality-score 10 \
        -L "$TRUTH_BED" \
        --native-pair-hmm-threads 8 \
        --sample-name 5801-diagnosis \
        2>&1 | tee -a "$LOG_FILE"

    # Normalize and index with bcftools
    echo "[STEP] Normalizing VCF..." | tee -a "$LOG_FILE"
    /users/pavb5f/.conda/envs/bio-cli/bin/bcftools norm -f "$REFERENCE" -m -any \
        "$OUTPUT_DIR/variants_raw.vcf.gz" | \
        /users/pavb5f/.conda/envs/bio-cli/bin/bcftools sort -Oz -o "$OUTPUT_DIR/variants.vcf.gz" \
        2>&1 | tee -a "$LOG_FILE"

    /users/pavb5f/.conda/envs/bio-cli/bin/bcftools index -t "$OUTPUT_DIR/variants.vcf.gz"

    echo "[SUCCESS] GATK analysis complete"
    echo "[INFO] Output: $OUTPUT_DIR/variants.vcf.gz"
}

# ============================================================================
# DeepVariant
# ============================================================================
run_deepvariant() {
    echo ""
    echo "========================================================================"
    echo "[CALLER] DeepVariant (PacBio model via Singularity)"
    echo "========================================================================"

    OUTPUT_DIR="$PROJECT_DIR/results/deepvariant"
    LOG_FILE="$LOG_DIR/deepvariant.log"
    CONTAINER="$PROJECT_DIR/containers/deepvariant_1.6.1.sif"
    SINGULARITY="/users/pavb5f/.conda/envs/bio-cli/bin/singularity"

    mkdir -p "$OUTPUT_DIR"

    # Check if container exists
    if [ ! -f "$CONTAINER" ]; then
        echo "[ERROR] DeepVariant container not found at $CONTAINER"
        echo "[ERROR] Please copy deepvariant_1.6.1.sif to containers/ directory"
        return 1
    fi

    echo "[INFO] Running DeepVariant..."
    echo "[INFO] Container: $CONTAINER"
    echo "[INFO] Output: $OUTPUT_DIR"
    echo "[INFO] Log: $LOG_FILE"

    $SINGULARITY exec \
        --bind /data:/data \
        --bind "$PROJECT_DIR:$PROJECT_DIR" \
        "$CONTAINER" \
        /opt/deepvariant/bin/run_deepvariant \
            --model_type=PACBIO \
            --ref="$REFERENCE" \
            --reads="$BAM" \
            --regions="$TRUTH_BED" \
            --output_vcf="$OUTPUT_DIR/variants.vcf.gz" \
            --num_shards=8 \
            2>&1 | tee "$LOG_FILE"

    # Index output VCF
    /users/pavb5f/.conda/envs/bio-cli/bin/bcftools index -t "$OUTPUT_DIR/variants.vcf.gz" 2>&1 | tee -a "$LOG_FILE"

    echo "[SUCCESS] DeepVariant analysis complete"
    echo "[INFO] Output: $OUTPUT_DIR/variants.vcf.gz"
}

# ============================================================================
# Clair3-RNA
# ============================================================================
run_clair3_rna() {
    echo ""
    echo "========================================================================"
    echo "[CALLER] Clair3-RNA (Long-read RNA variant calling)"
    echo "========================================================================"

    OUTPUT_DIR="$PROJECT_DIR/results/clair3_rna"
    LOG_FILE="$LOG_DIR/clair3_rna.log"

    mkdir -p "$OUTPUT_DIR"

    echo "[INFO] Running Clair3-RNA..."
    echo "[INFO] Output: $OUTPUT_DIR"
    echo "[INFO] Log: $LOG_FILE"

    CLAIR3_ENV="/users/pavb5f/.conda/envs/clair3-rna/bin/run_clair3.sh"

    if [ ! -x "$CLAIR3_ENV" ]; then
        echo "[ERROR] Clair3-RNA not found at $CLAIR3_ENV"
        echo "[ERROR] Please create clair3-rna conda environment first"
        return 1
    fi

    # Run Clair3-RNA
    $CLAIR3_ENV \
        --bam_fn="$BAM" \
        --ref_fn="$REFERENCE" \
        --threads=16 \
        --platform="hifi" \
        --output="$OUTPUT_DIR" \
        --bed_fn="$TRUTH_BED" \
        2>&1 | tee "$LOG_FILE"

    # Index output VCF
    if [ -f "$OUTPUT_DIR/merge_output.vcf.gz" ]; then
        /users/pavb5f/.conda/envs/bio-cli/bin/bcftools index -t "$OUTPUT_DIR/merge_output.vcf.gz" 2>&1 | tee -a "$LOG_FILE"
        echo "[SUCCESS] Clair3-RNA analysis complete"
        echo "[INFO] Output: $OUTPUT_DIR/merge_output.vcf.gz"
    else
        echo "[ERROR] Clair3-RNA did not produce output VCF"
    fi
}

# ============================================================================
# LongcallR
# ============================================================================
run_longcallr() {
    echo ""
    echo "========================================================================"
    echo "[CALLER] LongcallR (Long-read variant calling)"
    echo "========================================================================"

    OUTPUT_DIR="$PROJECT_DIR/results/longcallr"
    LOG_FILE="$LOG_DIR/longcallr.log"

    mkdir -p "$OUTPUT_DIR"

    echo "[INFO] Running LongcallR..."
    echo "[INFO] Output: $OUTPUT_DIR"
    echo "[INFO] Log: $LOG_FILE"

    LONGCALLR_BIN="/users/pavb5f/.conda/envs/longcallr/bin/longcallr"

    if [ ! -x "$LONGCALLR_BIN" ]; then
        echo "[ERROR] LongcallR not found at $LONGCALLR_BIN"
        echo "[ERROR] Please create longcallr conda environment first"
        return 1
    fi

    # Run LongcallR
    $LONGCALLR_BIN variants \
        -f "$REFERENCE" \
        -b "$BAM" \
        -o "$OUTPUT_DIR/variants.vcf" \
        -r "$TRUTH_BED" \
        --threads 16 \
        2>&1 | tee "$LOG_FILE"

    # Compress and index
    if [ -f "$OUTPUT_DIR/variants.vcf" ]; then
        /users/pavb5f/.conda/envs/bio-cli/bin/bgzip -f "$OUTPUT_DIR/variants.vcf"
        /users/pavb5f/.conda/envs/bio-cli/bin/bcftools index -t "$OUTPUT_DIR/variants.vcf.gz" 2>&1 | tee -a "$LOG_FILE"
        echo "[SUCCESS] LongcallR analysis complete"
        echo "[INFO] Output: $OUTPUT_DIR/variants.vcf.gz"
    else
        echo "[ERROR] LongcallR did not produce output VCF"
    fi
}

# ============================================================================
# Main entry point
# ============================================================================
case "$CALLER" in
    lab_supervised)
        run_lab_supervised
        ;;
    lab_unsupervised)
        run_lab_unsupervised
        ;;
    gatk)
        run_gatk
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
        run_lab_supervised
        run_lab_unsupervised
        run_gatk
        run_deepvariant
        run_clair3_rna
        run_longcallr
        ;;
    *)
        echo "[ERROR] Unknown caller: $CALLER"
        echo "Usage: $0 [lab_supervised|lab_unsupervised|gatk|deepvariant|clair3_rna|longcallr|all]"
        exit 1
        ;;
esac

echo ""
echo "========================================================================"
echo "[COMPLETE] Caller execution finished: $CALLER"
echo "========================================================================"
echo "[INFO] Results available in: $PROJECT_DIR/results/"
echo ""
