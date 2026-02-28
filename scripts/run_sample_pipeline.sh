#!/usr/bin/env bash
# =============================================================================
# run_sample_pipeline.sh
#
# Full end-to-end variant calling pipeline for a single PacBio Iso-Seq BAM.
# Checkpoint-aware: each step is skipped if its primary output already exists.
#
# Usage:
#   bash scripts/run_sample_pipeline.sh <bam_path>
#
# Example:
#   bash scripts/run_sample_pipeline.sh bams/MDS5801-1229IR_CD34.bam
#
# Pipeline steps:
#   1. Synthetic quality preprocessing (QUAL=* → Q30, required by GATK)
#   2. GATK HaplotypeCaller
#   3. DeepVariant (PACBIO model)
#   4. Clair3-RNA (hifi model, per-chromosome)
#   5. LongcallR (hifi-isoseq preset)
#   6. VEP annotation + somatic classification (per caller)
#   7. Multi-caller concordance merge
#   8. mosdepth BAM_DP > 1500 filter
#   9. GENE + BAM allele count annotation
#  10. TSV export (human-readable table)
#
# Output layout:
#   results/{SAMPLE}/haplotypecaller/variants_raw.vcf
#   results/{SAMPLE}/deepvariant/variants.vcf.gz
#   results/{SAMPLE}/clair3_rna/variants.vcf.gz
#   results/{SAMPLE}/longcallr/variants.vcf.gz
#   results/{SAMPLE}/merged/high_confidence_final.{vcf.gz,tsv}
#   results/{SAMPLE}/merged/medium_confidence_final.{vcf.gz,tsv}
#   results/{SAMPLE}/merged/low_confidence_final.{vcf.gz,tsv}
#   annotation/{SAMPLE}/{gatk,deepvariant,clair3,longcallr}/annotated.vcf.gz
#   annotation/{SAMPLE}/{gatk,deepvariant,clair3,longcallr}/somatic_known.vcf.gz ...
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Arguments and sample name
# ---------------------------------------------------------------------------
INPUT_BAM="${1:?Usage: $0 <bam_path>}"
INPUT_BAM="$(realpath "$INPUT_BAM")"
SAMPLE="$(basename "$INPUT_BAM" .bam)"

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
WORKDIR="/data/salomonis-archive/FASTQs/NCI-R01/rna_seq_varcall"
REFERENCE="$WORKDIR/reference/genome.fa"

# Per-sample directories
RESULTS="$WORKDIR/results/$SAMPLE"
ANNOT="$WORKDIR/annotation/$SAMPLE"
LOG_DIR="$WORKDIR/logs/$SAMPLE"
mkdir -p "$RESULTS" "$ANNOT" "$LOG_DIR"

# Preprocessed BAM (synthetic Q30 for GATK)
BAM_QUAL="$RESULTS/${SAMPLE}_qual.bam"

# Tool paths
PYTHON="/users/pavb5f/.conda/envs/bio-cli/bin/python"
SAMTOOLS="/users/pavb5f/.conda/envs/bio-cli/bin/samtools"
BGZIP="/users/pavb5f/.conda/envs/bio-cli/bin/bgzip"
BCFTOOLS="/users/pavb5f/.conda/envs/bio-cli/bin/bcftools"
SINGULARITY="/users/pavb5f/.conda/envs/bio-cli/bin/singularity"
MOSDEPTH="/usr/local/anaconda3-2020/envs/mosdepth/bin/mosdepth"
GATK="/users/pavb5f/.conda/envs/gatk-env/bin/gatk"

# Container paths
DV_SIF="$WORKDIR/containers/deepvariant_1.6.1.sif"
CLAIR3_SIF="$WORKDIR/containers/clair3_latest.sif"
LONGCALLR_SIF="$WORKDIR/containers/longcallr_1.12.0.sif"
VEP_SIF="$WORKDIR/containers/ensembl-vep_release_112.0.sif"
CLAIR3_MODEL="/opt/models/hifi"   # bundled inside clair3_latest.sif; same as original run_callers.sh
VEP_CACHE="$WORKDIR/vep_cache"
VEP_PLUGINS="$WORKDIR/vep_plugins"

THREADS=16

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
log()  { echo "[$(date '+%H:%M:%S')] [$SAMPLE] $*"; }
skip() { echo "[$(date '+%H:%M:%S')] [$SAMPLE] SKIP (exists): $1"; }
die()  { echo "[ERROR] $*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Step 1: Synthetic quality preprocessing
# ---------------------------------------------------------------------------
log "=== STEP 1: Synthetic quality ==="
if [ -f "$BAM_QUAL" ] && [ -f "${BAM_QUAL}.bai" ]; then
    skip "$BAM_QUAL"
else
    log "Adding synthetic Q30 scores to $INPUT_BAM ..."
    "$PYTHON" "$WORKDIR/scripts/add_synthetic_quality.py" \
        "$INPUT_BAM" "$BAM_QUAL" \
        2>&1 | tee "$LOG_DIR/synthetic_qual.log"
    "$SAMTOOLS" index "$BAM_QUAL"
    log "Synthetic quality done."
fi

# ---------------------------------------------------------------------------
# Step 2: GATK HaplotypeCaller
# ---------------------------------------------------------------------------
log "=== STEP 2: GATK HaplotypeCaller ==="
HC_DIR="$RESULTS/haplotypecaller"
HC_VCF="$HC_DIR/variants_raw.vcf"
mkdir -p "$HC_DIR"

HC_VCF_NORM="$HC_DIR/variants.vcf.gz"

if [ -f "${HC_VCF}.gz" ] || [ -f "$HC_VCF" ]; then
    skip "$HC_VCF"
else
    log "Running GATK HaplotypeCaller (canonical chromosomes only)..."
    PATH="/users/pavb5f/.conda/envs/gatk-env/bin:$PATH" \
    "$GATK" HaplotypeCaller \
        -R "$REFERENCE" \
        -I "$BAM_QUAL" \
        -O "$HC_VCF" \
        -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 \
        -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 \
        -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 \
        -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 \
        -L chr21 -L chr22 -L chrX -L chrY -L chrM \
        --dont-use-soft-clipped-bases \
        --min-base-quality-score 10 \
        --native-pair-hmm-threads "$THREADS" \
        2>&1 | tee "$LOG_DIR/haplotypecaller.log"
    log "GATK done."
fi

# Normalize GATK output (matching run_callers.sh behaviour)
if [ ! -f "$HC_VCF_NORM" ]; then
    _hc_input="${HC_VCF}.gz"; [ ! -f "$_hc_input" ] && _hc_input="$HC_VCF"
    log "Normalizing GATK VCF..."
    "$BCFTOOLS" norm -f "$REFERENCE" -m -any "$_hc_input" | \
        "$BCFTOOLS" sort -Oz -o "$HC_VCF_NORM"
    "$BCFTOOLS" index -t "$HC_VCF_NORM"
fi

# ---------------------------------------------------------------------------
# Step 3: DeepVariant
# ---------------------------------------------------------------------------
log "=== STEP 3: DeepVariant ==="
DV_DIR="$RESULTS/deepvariant"
DV_VCF="$DV_DIR/variants.vcf.gz"
mkdir -p "$DV_DIR"

if [ -f "$DV_VCF" ]; then
    skip "$DV_VCF"
else
    log "Running DeepVariant (PACBIO model, full genome)..."
    "$SINGULARITY" exec \
        --bind /data:/data \
        --bind "$WORKDIR:$WORKDIR" \
        "$DV_SIF" \
        /opt/deepvariant/bin/run_deepvariant \
            --model_type=PACBIO \
            --ref="$REFERENCE" \
            --reads="$BAM_QUAL" \
            --output_vcf="$DV_VCF" \
            --num_shards="$THREADS" \
            --regions "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM" \
            --intermediate_results_dir="$DV_DIR/tmp" \
        2>&1 | tee "$LOG_DIR/deepvariant.log"
    log "DeepVariant done."
fi

# ---------------------------------------------------------------------------
# Step 4: Clair3-RNA
# ---------------------------------------------------------------------------
log "=== STEP 4: Clair3-RNA ==="
C3_DIR="$RESULTS/clair3_rna"
C3_VCF="$C3_DIR/variants.vcf.gz"
mkdir -p "$C3_DIR"

if [ -f "$C3_VCF" ]; then
    skip "$C3_VCF"
else
    log "Running Clair3-RNA (hifi model, full genome)..."
    "$SINGULARITY" exec \
        --bind /data:/data \
        --bind "$WORKDIR:$WORKDIR" \
        "$CLAIR3_SIF" \
        /opt/bin/run_clair3.sh \
            --bam_fn="$BAM_QUAL" \
            --ref_fn="$REFERENCE" \
            --threads="$THREADS" \
            --platform="hifi" \
            --model_path="$CLAIR3_MODEL" \
            --output="$C3_DIR" \
            --sample_name="$SAMPLE" \
            --min_mq=5 \
            --snp_min_af=0.05 \
            --indel_min_af=0.05 \
        2>&1 | tee "$LOG_DIR/clair3.log"

    # Clair3 outputs merge_output.vcf.gz — symlink to standard name
    if [ -f "$C3_DIR/merge_output.vcf.gz" ] && [ ! -f "$C3_VCF" ]; then
        cp "$C3_DIR/merge_output.vcf.gz" "$C3_VCF"
        cp "$C3_DIR/merge_output.vcf.gz.tbi" "${C3_VCF}.tbi" 2>/dev/null || \
            "$BCFTOOLS" index -t "$C3_VCF"
    fi
    log "Clair3-RNA done."
fi

# ---------------------------------------------------------------------------
# Step 5: LongcallR
# ---------------------------------------------------------------------------
log "=== STEP 5: LongcallR ==="
LCR_DIR="$RESULTS/longcallr"
LCR_VCF="$LCR_DIR/variants.vcf.gz"
mkdir -p "$LCR_DIR"

if [ -f "$LCR_VCF" ]; then
    skip "$LCR_VCF"
else
    log "Running LongcallR (hifi-isoseq preset, full genome)..."
    "$SINGULARITY" exec \
        --bind /data:/data \
        --bind "$WORKDIR:$WORKDIR" \
        "$LONGCALLR_SIF" \
        longcallR \
            --bam-path "$BAM_QUAL" \
            --ref-path "$REFERENCE" \
            --output "$LCR_DIR/variants" \
            --preset hifi-isoseq \
            --threads "$THREADS" \
            --min-allele-freq 0.05 \
            --low-allele-frac-cutoff 0.03 \
            --min-mapq 1 \
            --no-bam-output \
        2>&1 | tee "$LOG_DIR/longcallr.log"

    # LongcallR outputs an unsorted VCF — sort before bgzip+index
    if [ -f "$LCR_DIR/variants.vcf" ]; then
        "$BCFTOOLS" sort -Oz -o "$LCR_VCF" "$LCR_DIR/variants.vcf"
        rm -f "$LCR_DIR/variants.vcf"
    else
        # Already bgzipped but possibly unsorted — re-sort in place
        "$BCFTOOLS" sort -Oz -o "${LCR_VCF}.sorted" "$LCR_VCF" && mv "${LCR_VCF}.sorted" "$LCR_VCF"
    fi
    "$BCFTOOLS" index -t "$LCR_VCF"
    log "LongcallR done."
fi

# ---------------------------------------------------------------------------
# Step 6: VEP annotation + somatic classification (per caller)
# ---------------------------------------------------------------------------
log "=== STEP 6: VEP annotation + classification ==="

declare -A CALLER_VCFS=(
    [gatk]="$HC_VCF_NORM"
    [deepvariant]="$DV_VCF"
    [clair3]="$C3_VCF"
    [longcallr]="$LCR_VCF"
)

for caller in gatk deepvariant clair3 longcallr; do
    CALLER_ANNOT_DIR="$ANNOT/$caller"
    mkdir -p "$CALLER_ANNOT_DIR"

    if [ -f "$CALLER_ANNOT_DIR/somatic_known.vcf.gz" ]; then
        skip "$CALLER_ANNOT_DIR/annotated.vcf.gz"
        continue
    fi

    INPUT_VCF="${CALLER_VCFS[$caller]}"
    # GATK outputs plain VCF — check both forms
    [ "$caller" = "gatk" ] && [ -f "${INPUT_VCF}.gz" ] && INPUT_VCF="${INPUT_VCF}.gz"
    [ ! -f "$INPUT_VCF" ] && { log "WARN: $caller VCF not found, skipping annotation."; continue; }

    log "Annotating $caller with VEP..."
    "$SINGULARITY" exec \
        --bind /data:/data \
        --bind "$WORKDIR:$WORKDIR" \
        "$VEP_SIF" \
        vep \
            --input_file "$INPUT_VCF" \
            --output_file "$CALLER_ANNOT_DIR/annotated.vcf" \
            --format vcf \
            --vcf \
            --cache \
            --dir_cache "$VEP_CACHE" \
            --dir_plugins "$VEP_PLUGINS" \
            --assembly GRCh38 \
            --species homo_sapiens \
            --fasta "$REFERENCE" \
            --fork "$THREADS" \
            --everything \
            --af_gnomADe \
            --af_gnomADg \
            --check_existing \
            --no_escape \
            --canonical \
            --symbol \
            --hgvs \
            --offline \
        2>&1 | tee "$CALLER_ANNOT_DIR/vep.log"

    "$BCFTOOLS" sort "$CALLER_ANNOT_DIR/annotated.vcf" \
        -O z -o "$CALLER_ANNOT_DIR/annotated.vcf.gz"
    "$BCFTOOLS" index -t "$CALLER_ANNOT_DIR/annotated.vcf.gz"
    rm "$CALLER_ANNOT_DIR/annotated.vcf"

    log "Classifying $caller variants..."
    "$PYTHON" "$WORKDIR/scripts/classify_variants.py" \
        --input  "$CALLER_ANNOT_DIR/annotated.vcf.gz" \
        --output "$CALLER_ANNOT_DIR" \
        --sample "${SAMPLE}_${caller}" \
        2>&1 | tee -a "$CALLER_ANNOT_DIR/vep.log"

    log "$caller annotation + classification done."
done

# ---------------------------------------------------------------------------
# Step 7: Multi-caller concordance merge
# ---------------------------------------------------------------------------
log "=== STEP 7: Multi-caller merge ==="
MERGED_DIR="$RESULTS/merged"
mkdir -p "$MERGED_DIR"

if [ -f "$MERGED_DIR/all_candidates.vcf.gz" ]; then
    skip "$MERGED_DIR/all_candidates.vcf.gz"
else
    log "Running merge_callers.py..."
    "$PYTHON" "$WORKDIR/scripts/merge_callers.py" \
        --annotation_dir "$ANNOT" \
        --output "$MERGED_DIR" \
        2>&1 | tee "$LOG_DIR/merge.log"
    log "Merge done."
fi

# ---------------------------------------------------------------------------
# Step 8: mosdepth BAM_DP > 1500 filter
# ---------------------------------------------------------------------------
log "=== STEP 8: mosdepth BAM_DP > 1500 filter ==="

_apply_bamdp_filter() {
    local label="$1"
    local input_vcf="$2"
    local output_vcf="$3"
    local tmp_prefix="$RESULTS/tmp_mosdepth_${label}"

    [ -f "$output_vcf" ] && { skip "$output_vcf"; return; }

    log "  Extracting positions for $label..."
    mkdir -p "$(dirname "$tmp_prefix")"
    "$BCFTOOLS" query -f '%CHROM\t%POS0\t%POS\n' "$input_vcf" \
        | sort -k1,1 -k2,2n > "${tmp_prefix}.bed"

    log "  Running mosdepth for $label..."
    "$MOSDEPTH" \
        --by "${tmp_prefix}.bed" \
        --no-per-base \
        --threads 8 \
        "${tmp_prefix}_md" \
        "$BAM_QUAL"

    log "  Building depth TSV..."
    zcat "${tmp_prefix}_md.regions.bed.gz" \
        | awk 'BEGIN{OFS="\t"} {print $1, $3, int($4)}' \
        | "$BGZIP" -c > "${tmp_prefix}_depth.tsv.gz"
    /users/pavb5f/.conda/envs/bio-cli/bin/tabix \
        -s1 -b2 -e2 "${tmp_prefix}_depth.tsv.gz"

    echo '##INFO=<ID=BAM_DP,Number=1,Type=Integer,Description="Raw BAM depth from mosdepth">' \
        > "${tmp_prefix}_hdr.txt"

    log "  Annotating + filtering $label (BAM_DP > 1500)..."
    "$BCFTOOLS" annotate \
        -a "${tmp_prefix}_depth.tsv.gz" \
        -c CHROM,POS,INFO/BAM_DP \
        -h "${tmp_prefix}_hdr.txt" \
        "$input_vcf" \
        -Oz -o "${tmp_prefix}_annotated.vcf.gz"
    "$BCFTOOLS" index -t "${tmp_prefix}_annotated.vcf.gz"

    "$BCFTOOLS" filter \
        -i 'INFO/BAM_DP > 1500' \
        "${tmp_prefix}_annotated.vcf.gz" \
        -Oz -o "$output_vcf"
    "$BCFTOOLS" index -t "$output_vcf"

    # Clean up temp files
    rm -f "${tmp_prefix}".bed "${tmp_prefix}_md".* "${tmp_prefix}_depth".* \
          "${tmp_prefix}_hdr.txt" "${tmp_prefix}_annotated".vcf.gz*

    N=$("$BCFTOOLS" view -H "$output_vcf" | wc -l)
    log "  $label: $N variants with BAM_DP > 1500"
}

_run_gatk_nongermline_filter() {
    local out="$ANNOT/gatk/non_germline_bamdp1500.vcf.gz"
    [ -f "$out" ] && { skip "$out"; return; }

    log "  Concatenating GATK non-germline buckets..."
    TMP_CONCAT="$RESULTS/tmp_gatk_nongermline.vcf.gz"
    "$BCFTOOLS" concat -a \
        "$ANNOT/gatk/somatic_known.vcf.gz" \
        "$ANNOT/gatk/novel_candidates.vcf.gz" \
        "$ANNOT/gatk/rare_candidates.vcf.gz" \
        "$ANNOT/gatk/uncertain.vcf.gz" \
        -O z -o "$TMP_CONCAT"
    "$BCFTOOLS" index -t "$TMP_CONCAT"

    local tmp_prefix="$RESULTS/tmp_mosdepth_gatk_ng"
    "$BCFTOOLS" query -f '%CHROM\t%POS0\t%POS\n' "$TMP_CONCAT" \
        | sort -k1,1 -k2,2n > "${tmp_prefix}.bed"
    "$MOSDEPTH" --by "${tmp_prefix}.bed" --no-per-base --threads 8 \
        "${tmp_prefix}_md" "$BAM_QUAL"
    zcat "${tmp_prefix}_md.regions.bed.gz" \
        | awk 'BEGIN{OFS="\t"} {print $1, $3, int($4)}' \
        | "$BGZIP" -c > "${tmp_prefix}_depth.tsv.gz"
    /users/pavb5f/.conda/envs/bio-cli/bin/tabix -s1 -b2 -e2 "${tmp_prefix}_depth.tsv.gz"
    echo '##INFO=<ID=BAM_DP,Number=1,Type=Integer,Description="Raw BAM depth from mosdepth">' \
        > "${tmp_prefix}_hdr.txt"
    "$BCFTOOLS" annotate \
        -a "${tmp_prefix}_depth.tsv.gz" -c CHROM,POS,INFO/BAM_DP \
        -h "${tmp_prefix}_hdr.txt" "$TMP_CONCAT" -Oz -o "${tmp_prefix}_ann.vcf.gz"
    "$BCFTOOLS" index -t "${tmp_prefix}_ann.vcf.gz"
    "$BCFTOOLS" filter -i 'INFO/BAM_DP > 1500' "${tmp_prefix}_ann.vcf.gz" \
        -Oz -o "$out"
    "$BCFTOOLS" index -t "$out"
    rm -f "$TMP_CONCAT" "${TMP_CONCAT}.tbi" "${tmp_prefix}".* "${tmp_prefix}_md".* \
          "${tmp_prefix}_depth".* "${tmp_prefix}_hdr.txt" "${tmp_prefix}_ann".vcf.gz*
    N=$("$BCFTOOLS" view -H "$out" | wc -l)
    log "  GATK non-germline BAM_DP>1500: $N variants"
}
_run_gatk_nongermline_filter

# Merged confidence tiers
for tier in high_confidence medium_confidence low_confidence; do
    _apply_bamdp_filter "$tier" \
        "$MERGED_DIR/${tier}.vcf.gz" \
        "$MERGED_DIR/${tier}_bamdp1500.vcf.gz"
done

# ---------------------------------------------------------------------------
# Step 9: GENE + BAM allele count annotation
# ---------------------------------------------------------------------------
log "=== STEP 9: BAM allele annotations ==="

for tier in non_germline high_confidence medium_confidence low_confidence; do
    if [ "$tier" = "non_germline" ]; then
        input="$ANNOT/gatk/non_germline_bamdp1500.vcf.gz"
        output="$ANNOT/gatk/non_germline_final.vcf.gz"
    else
        input="$MERGED_DIR/${tier}_bamdp1500.vcf.gz"
        output="$MERGED_DIR/${tier}_final.vcf.gz"
    fi

    [ -f "$output" ] && { skip "$output"; continue; }
    [ ! -f "$input" ] && { log "WARN: $input not found, skipping."; continue; }

    N=$("$BCFTOOLS" view -H "$input" | wc -l)
    log "  Annotating $tier ($N variants)..."
    "$PYTHON" "$WORKDIR/scripts/add_bam_annotations.py" \
        --input  "$input" \
        --output "$output" \
        --bam    "$BAM_QUAL" \
        --ref    "$REFERENCE" \
        2>&1 | tee "$LOG_DIR/bam_annot_${tier}.log"
done

# ---------------------------------------------------------------------------
# Step 10: TSV export
# ---------------------------------------------------------------------------
log "=== STEP 10: TSV export ==="

for tier in high_confidence medium_confidence low_confidence; do
    vcf="$MERGED_DIR/${tier}_final.vcf.gz"
    tsv="$MERGED_DIR/${tier}_final.tsv"
    [ -f "$tsv" ] && { skip "$tsv"; continue; }
    [ ! -f "$vcf" ] && continue
    "$PYTHON" "$WORKDIR/scripts/vcf_to_table.py" \
        --input "$vcf" --output "$tsv" \
        2>&1 | tee -a "$LOG_DIR/tsv_export.log"
done

vcf="$ANNOT/gatk/non_germline_final.vcf.gz"
tsv="$ANNOT/gatk/non_germline_final.tsv"
if [ -f "$vcf" ] && [ ! -f "$tsv" ]; then
    "$PYTHON" "$WORKDIR/scripts/vcf_to_table.py" \
        --input "$vcf" --output "$tsv" \
        2>&1 | tee -a "$LOG_DIR/tsv_export.log"
fi

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------
log "========================================================"
log "PIPELINE COMPLETE for sample: $SAMPLE"
log "Primary outputs:"
log "  $MERGED_DIR/high_confidence_final.tsv"
log "  $MERGED_DIR/high_confidence_final.vcf.gz"
log "  $ANNOT/gatk/non_germline_final.tsv"
log "========================================================"
