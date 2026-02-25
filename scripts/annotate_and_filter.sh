#!/usr/bin/env bash
# =============================================================================
# annotate_and_filter.sh
# Annotate VCFs with VEP (gnomAD, dbSNP, ClinVar, CancerHotspots, COSMIC)
# then classify variants into SOMATIC_KNOWN / NOVEL_CANDIDATE / GERMLINE
#
# Usage:
#   bash scripts/annotate_and_filter.sh <input.vcf.gz> <sample_name>
#   bash scripts/annotate_and_filter.sh results/haplotypecaller/variants_raw.vcf gatk
# =============================================================================
set -euo pipefail

WORKDIR="/data/salomonis-archive/FASTQs/NCI-R01/rna_seq_varcall"
VEP_SIF="$WORKDIR/containers/ensembl-vep_release_112.0.sif"
VEP_CACHE="$WORKDIR/vep_cache"
VEP_PLUGINS="$WORKDIR/vep_plugins"
REFERENCE="$WORKDIR/reference/genome.fa"
SINGULARITY_BIN="/users/pavb5f/.conda/envs/bio-cli/bin/singularity"
BCFTOOLS_BIN="/users/pavb5f/.conda/envs/bio-cli/bin/bcftools"

INPUT_VCF="${1:?Usage: $0 <input.vcf.gz> <sample_name>}"
SAMPLE="${2:?Usage: $0 <input.vcf.gz> <sample_name>}"
OUT_DIR="$WORKDIR/annotation/${SAMPLE}"
mkdir -p "$OUT_DIR"

echo "========================================================================"
echo "[VEP] Annotating: $INPUT_VCF"
echo "[VEP] Sample:     $SAMPLE"
echo "[VEP] Output:     $OUT_DIR"
echo "========================================================================"

# -----------------------------------------------------------------------------
# Step 1: VEP annotation
# Enables: dbSNP IDs, gnomAD AF (exome + genome), ClinVar, canonical transcript,
#          consequence, gene symbols, SIFT/PolyPhen, existing variant lookup
# -----------------------------------------------------------------------------
"$SINGULARITY_BIN" exec \
    --bind /data:/data \
    --bind "$WORKDIR:$WORKDIR" \
    "$VEP_SIF" \
    vep \
    --input_file "$INPUT_VCF" \
    --output_file "$OUT_DIR/annotated.vcf" \
    --format vcf \
    --vcf \
    --cache \
    --dir_cache "$VEP_CACHE" \
    --dir_plugins "$VEP_PLUGINS" \
    --assembly GRCh38 \
    --species homo_sapiens \
    --fasta "$REFERENCE" \
    --fork 16 \
    --everything \
    --af_gnomADe \
    --af_gnomADg \
    --check_existing \
    --no_escape \
    --canonical \
    --symbol \
    --hgvs \
    --offline \
    2>&1 | tee "$OUT_DIR/vep.log"

# Compress and index
"$BCFTOOLS_BIN" sort "$OUT_DIR/annotated.vcf" \
    -O z -o "$OUT_DIR/annotated.vcf.gz"
"$BCFTOOLS_BIN" index -t "$OUT_DIR/annotated.vcf.gz"
rm "$OUT_DIR/annotated.vcf"

echo "[VEP] Annotation complete: $OUT_DIR/annotated.vcf.gz"

# -----------------------------------------------------------------------------
# Step 2: Classification via Python decision logic
# -----------------------------------------------------------------------------
python3 "$WORKDIR/scripts/classify_variants.py" \
    --input    "$OUT_DIR/annotated.vcf.gz" \
    --output   "$OUT_DIR" \
    --sample   "$SAMPLE"

echo "========================================================================"
echo "[DONE] Results in $OUT_DIR/"
echo "  somatic_known.vcf.gz      — in COSMIC/CancerHotspots/ClinVar Pathogenic"
echo "  novel_candidates.vcf.gz   — absent from gnomAD, high quality"
echo "  rare_candidates.vcf.gz    — gnomAD AF < 0.0001, not confirmed germline"
echo "  germline.vcf.gz           — gnomAD AF > 0.001, excluded"
echo "========================================================================"
