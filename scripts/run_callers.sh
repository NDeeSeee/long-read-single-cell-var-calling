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

    OUTPUT_DIR="$PROJECT_DIR/results/gatk"
    LOG_FILE="$LOG_DIR/gatk.log"

    mkdir -p "$OUTPUT_DIR"

    echo "[INFO] Running GATK HaplotypeCaller..."
    echo "[INFO] Output: $OUTPUT_DIR"
    echo "[INFO] Log: $LOG_FILE"

    # Call variants
    echo "[STEP] HaplotypeCaller..." | tee "$LOG_FILE"
    gatk HaplotypeCaller \
        -R "$REFERENCE" \
        -I "$BAM" \
        -O "$OUTPUT_DIR/variants_raw.vcf.gz" \
        --dont-use-soft-clipped-bases \
        --standard-min-confidence-threshold-for-calling 20.0 \
        --min-base-quality-score 10 \
        -L "$TRUTH_BED" \
        --native-pair-hmm-threads 8 \
        2>&1 | tee -a "$LOG_FILE"

    # Filter variants
    echo "[STEP] VariantFiltration..." | tee -a "$LOG_FILE"
    gatk VariantFiltration \
        -R "$REFERENCE" \
        -V "$OUTPUT_DIR/variants_raw.vcf.gz" \
        -O "$OUTPUT_DIR/variants_filtered.vcf.gz" \
        --filter-name "LowQual" --filter-expression "QUAL < 30.0" \
        --filter-name "LowDepth" --filter-expression "DP < 10" \
        2>&1 | tee -a "$LOG_FILE"

    echo "[SUCCESS] GATK analysis complete"
    echo "[INFO] Output: $OUTPUT_DIR/variants_filtered.vcf.gz"
}

# ============================================================================
# DeepVariant
# ============================================================================
run_deepvariant() {
    echo ""
    echo "========================================================================"
    echo "[CALLER] DeepVariant (WGS model via Singularity)"
    echo "========================================================================"

    OUTPUT_DIR="$PROJECT_DIR/results/deepvariant"
    LOG_FILE="$LOG_DIR/deepvariant.log"
    CONTAINER="$PROJECT_DIR/containers/deepvariant_1.6.1.sif"

    mkdir -p "$OUTPUT_DIR"

    # Check if container exists
    if [ ! -f "$CONTAINER" ]; then
        echo "[INFO] DeepVariant container not found. Downloading..."
        mkdir -p "$PROJECT_DIR/containers"
        cd "$PROJECT_DIR/containers"
        singularity pull docker://google/deepvariant:1.6.1
        cd "$PROJECT_DIR"
    fi

    echo "[INFO] Running DeepVariant..."
    echo "[INFO] Container: $CONTAINER"
    echo "[INFO] Output: $OUTPUT_DIR"
    echo "[INFO] Log: $LOG_FILE"

    singularity exec \
        --bind /data:/data \
        "$CONTAINER" \
        /opt/deepvariant/bin/run_deepvariant \
            --model_type=WGS \
            --ref="$REFERENCE" \
            --reads="$BAM" \
            --regions="$TRUTH_BED" \
            --output_vcf="$OUTPUT_DIR/variants.vcf.gz" \
            --output_gvcf="$OUTPUT_DIR/variants.g.vcf.gz" \
            --num_shards=8 \
            --make_examples_extra_args="min_mapping_quality=1" \
            2>&1 | tee "$LOG_FILE"

    echo "[SUCCESS] DeepVariant analysis complete"
    echo "[INFO] Output: $OUTPUT_DIR/variants.vcf.gz"
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
    all)
        run_lab_supervised
        run_lab_unsupervised
        run_gatk
        run_deepvariant
        ;;
    *)
        echo "[ERROR] Unknown caller: $CALLER"
        echo "Usage: $0 [lab_supervised|lab_unsupervised|gatk|deepvariant|all]"
        exit 1
        ;;
esac

echo ""
echo "========================================================================"
echo "[COMPLETE] Caller execution finished: $CALLER"
echo "========================================================================"
echo "[INFO] Results available in: $PROJECT_DIR/results/"
echo ""
