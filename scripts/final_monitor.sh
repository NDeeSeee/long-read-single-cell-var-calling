#!/bin/bash
# Final monitoring and comparison script

echo ""
echo "════════════════════════════════════════════════════════════════"
echo "          SOMATIC VARIANT CALLING - FINAL PROGRESS"
echo "════════════════════════════════════════════════════════════════"
echo ""
echo "Updated: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# Truth set
echo "TARGET: 8 known somatic variants"
echo "  - RUNX1 p.W279* (35% VAF)"
echo "  - ASXL1 p.G643fs (8% VAF)"
echo "  - SETBP1 p.G870S (28% VAF)"
echo "  - SRSF2 p.P95R (37% VAF)"
echo "  - BCOR, CBL, MEF2B, PHF6 (VAF unknown)"
echo ""

# DeepVariant
echo "──────────────────────────────────────────────────────────────"
echo "DeepVariant (Deep Learning Caller)"
echo "──────────────────────────────────────────────────────────────"
if [ -f "results/deepvariant/variants.vcf.gz" ]; then
  variant_count=$(bcftools view -H results/deepvariant/variants.vcf.gz 2>/dev/null | wc -l)
  echo "✓ COMPLETE: $variant_count variants called"
  echo ""
  echo "Sample variants (first 5):"
  bcftools view -H results/deepvariant/variants.vcf.gz 2>/dev/null | head -5 | while read line; do
    echo "  $line" | awk '{printf "    %s:%s %s→%s (AF=", $1, $2, $4, $5; if($8 ~ /AF=/) {match($8, /AF=([^;]*)/,a); printf "%s", a[1]} else {printf "?"} print ")"}'
  done
else
  echo "⏳ Running... (converting SIF, preparing examples)"
fi
echo ""

# Freebayes
echo "──────────────────────────────────────────────────────────────"
echo "Freebayes (Bayesian Caller - Mutect2 Alternative)"
echo "──────────────────────────────────────────────────────────────"
if [ -f "results/freebayes/variants.vcf.gz" ]; then
  variant_count=$(bcftools view -H results/freebayes/variants.vcf.gz 2>/dev/null | wc -l)
  echo "✓ COMPLETE: $variant_count variants called"
  echo ""
  echo "Sample variants (first 5):"
  bcftools view -H results/freebayes/variants.vcf.gz 2>/dev/null | head -5 | while read line; do
    echo "  $line" | awk '{printf "    %s:%s %s→%s (AF=", $1, $2, $4, $5; if($8 ~ /AF=/) {match($8, /AF=([^;]*)/,a); printf "%s", a[1]} else {printf "?"} print ")"}'
  done
else
  echo "⏳ Running..."
fi
echo ""

# Mutect2
echo "──────────────────────────────────────────────────────────────"
echo "Mutect2 (GATK Somatic Caller)"
echo "──────────────────────────────────────────────────────────────"
if [ -f "results/mutect2/variants_raw.vcf.gz" ]; then
  echo "✓ COMPLETE"
else
  echo "✗ BLOCKED: Java library incompatibility (JLI_StringDup)"
  echo "  → Using Freebayes as alternative"
fi
echo ""

echo "════════════════════════════════════════════════════════════════"
echo "Next Steps:"
echo "  1. Wait for DeepVariant and Freebayes to complete"
echo "  2. Run comparison: bash scripts/compare_final.sh"
echo "  3. View results in results/comparison/"
echo "════════════════════════════════════════════════════════════════"
echo ""
