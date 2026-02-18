#!/bin/bash
#
# Verify that the variant caller comparison pipeline is properly set up
#

set +e  # Don't exit on errors - we want to report all issues

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

echo "========================================================================"
echo "Variant Caller Comparison Pipeline - Setup Verification"
echo "========================================================================"
echo ""

ERRORS=0
WARNINGS=0

# Function to print success
success() {
    echo "  ✓ $1"
}

# Function to print warning
warn() {
    echo "  ⚠ $1"
    ((WARNINGS++))
}

# Function to print error
error() {
    echo "  ✗ $1"
    ((ERRORS++))
}

# ============================================================================
# Check Project Structure
# ============================================================================
echo "1. Project Structure"

[ -f "$PROJECT_DIR/README.md" ] && success "README.md" || error "README.md not found"
[ -f "$PROJECT_DIR/CLAUDE.md" ] && success "CLAUDE.md" || error "CLAUDE.md not found"
[ -f "$PROJECT_DIR/IMPLEMENTATION_STATUS.md" ] && success "IMPLEMENTATION_STATUS.md" || error "IMPLEMENTATION_STATUS.md not found"
[ -f "$PROJECT_DIR/test_variants.txt" ] && success "test_variants.txt" || error "test_variants.txt not found"

echo ""
echo "2. Scripts"

for script in prepare_truth_vcf.py compress_vcf.py tsv_to_vcf.py evaluate_callers.py generate_report.py run_callers.sh; do
    [ -f "$PROJECT_DIR/scripts/$script" ] && success "$script" || error "$script not found"
done

echo ""
echo "3. Truth Set Files"

[ -f "$PROJECT_DIR/truth_set/truth_set.vcf" ] && success "truth_set.vcf" || error "truth_set.vcf not found"
[ -f "$PROJECT_DIR/truth_set/truth_set.vcf.gz" ] && success "truth_set.vcf.gz" || error "truth_set.vcf.gz not found"
[ -f "$PROJECT_DIR/truth_set/truth_regions.bed" ] && success "truth_regions.bed" || error "truth_regions.bed not found"
[ -f "$PROJECT_DIR/truth_set/test_variants_formatted.txt" ] && success "test_variants_formatted.txt" || error "test_variants_formatted.txt not found"

echo ""
echo "4. Result Directories"

for dir in supervised_extraction global_snv gatk deepvariant comparison; do
    [ -d "$PROJECT_DIR/results/$dir" ] && success "results/$dir/" || error "results/$dir/ not found"
done

# ============================================================================
# Check External Tools & Data
# ============================================================================
echo ""
echo "5. External Tools"

which samtools > /dev/null 2>&1 && success "samtools" || error "samtools not found"
which gatk > /dev/null 2>&1 && success "gatk" || error "gatk not found"
which bcftools > /dev/null 2>&1 && success "bcftools" || warn "bcftools not found (will install)"
which python3 > /dev/null 2>&1 && success "python3" || error "python3 not found"
which singularity > /dev/null 2>&1 && success "singularity" || warn "singularity not found (needed for DeepVariant)"

echo ""
echo "6. Python Modules"

python3 -c "import pysam" > /dev/null 2>&1 && success "pysam" || error "pysam not installed"
python3 -c "import scipy" > /dev/null 2>&1 && success "scipy" || error "scipy not installed"

echo ""
echo "7. Reference Data"

[ -f "/data/salomonis-archive/genomes/hg38/genome.fa" ] && success "Reference genome (3GB)" || error "Reference genome not found"
[ -f "/data/salomonis-archive/genomes/hg38/genome.fa.fai" ] && success "Reference index" || error "Reference index not found"
[ -f "/data/salomonis-archive/BAMs/Grimes/scRNA-Seq/KINNEX/5801-diagnosis/scisoseq.mapped.bam" ] && success "Test BAM (11GB)" || error "Test BAM not found"
[ -f "/data/salomonis-archive/BAMs/Grimes/scRNA-Seq/KINNEX/5801-diagnosis/scisoseq.mapped.bam.bai" ] && success "BAM index" || error "BAM index not found"

echo ""
echo "8. Lab Script Availability"

[ -f "/data/salomonis-archive/LabFiles/Nathan/Revio/altanalyze3/altanalyze3/components/bam/variant_extraction.py" ] && \
    success "variant_extraction.py (supervised)" || \
    error "variant_extraction.py not found"

[ -f "/data/salomonis-archive/LabFiles/Nathan/Revio/altanalyze3/altanalyze3/components/bam/global_snv.py" ] && \
    success "global_snv.py (unsupervised)" || \
    error "global_snv.py not found"

# ============================================================================
# Verification Summary
# ============================================================================
echo ""
echo "========================================================================"
echo "Verification Summary"
echo "========================================================================"
echo ""

if [ $ERRORS -eq 0 ] && [ $WARNINGS -eq 0 ]; then
    echo "✓ All checks passed! Pipeline is ready to run."
    echo ""
    echo "Next steps:"
    echo "  1. cd $PROJECT_DIR"
    echo "  2. ./scripts/run_callers.sh all"
    echo "  3. python scripts/evaluate_callers.py"
    echo "  4. python scripts/generate_report.py"
    exit 0
elif [ $ERRORS -eq 0 ]; then
    echo "✓ All critical checks passed ($WARNINGS warnings)"
    echo ""
    echo "Warnings:"
    echo "  - bcftools not installed (will install automatically)"
    echo "  - singularity not found (needed for DeepVariant)"
    echo ""
    echo "You can still run the pipeline. To install missing tools:"
    echo "  conda install -c bioconda bcftools=1.23"
    exit 0
else
    echo "✗ $ERRORS critical issue(s) found, $WARNINGS warning(s)"
    echo ""
    echo "Please address the errors above before running the pipeline."
    exit 1
fi
