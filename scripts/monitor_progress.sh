#!/bin/bash
# Monitor variant caller progress and results

echo "========================================"
echo "SOMATIC VARIANT CALLING PROGRESS"
echo "========================================"
echo ""
echo "Timestamp: $(date)"
echo ""

# Check Mutect2 status
echo "Mutect2 (Somatic Caller):"
if [ -f "results/mutect2/variants_raw.vcf.gz" ]; then
  vcf_lines=$(bcftools view results/mutect2/variants_raw.vcf.gz 2>/dev/null | grep -v "^#" | wc -l)
  echo "  ✓ VCF Created: $vcf_lines variants called"
else
  echo "  ⏳ Running... (checking for output file)"
fi
echo ""

# Check DeepVariant status  
echo "DeepVariant (Comparison Caller):"
if [ -f "results/deepvariant/variants.vcf.gz" ]; then
  vcf_lines=$(bcftools view results/deepvariant/variants.vcf.gz 2>/dev/null | grep -v "^#" | wc -l)
  echo "  ✓ VCF Created: $vcf_lines variants called"
else
  echo "  ⏳ Running..."
fi
echo ""

# Show truth set for comparison
echo "Truth Set (8 expected variants):"
bcftools view truth_set/truth_set.vcf.gz 2>/dev/null | grep -v "^##" | awk '{print "  " $0}'
echo ""

echo "========================================"
echo "Run this script again to check progress"
echo "========================================"
