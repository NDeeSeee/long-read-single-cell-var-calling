#!/usr/bin/env python3
"""
Compare VCF outputs from GATK and DeepVariant against truth set.
"""

import subprocess
import os

def get_variants_from_vcf(vcf_file):
    """Extract variants (chr:pos:ref:alt) from VCF"""
    variants = set()
    try:
        result = subprocess.run(
            f'bcftools query -f "%CHROM:%POS:%REF:%ALT\\n" {vcf_file}',
            shell=True,
            capture_output=True,
            text=True,
            timeout=60
        )
        for line in result.stdout.strip().split('\n'):
            if line:
                variants.add(line)
    except Exception as e:
        print(f"Error reading {vcf_file}: {e}")
    return variants

def get_vaf_from_vcf(vcf_file):
    """Extract VAF values from VCF (from AF field or computed from AD)"""
    vafs = {}
    try:
        result = subprocess.run(
            f'bcftools query -f "%CHROM:%POS:%REF:%ALT\t%AF\\n" {vcf_file}',
            shell=True,
            capture_output=True,
            text=True,
            timeout=60
        )
        for line in result.stdout.strip().split('\n'):
            if line:
                parts = line.split('\t')
                if len(parts) == 2:
                    key = parts[0]
                    try:
                        vaf = float(parts[1])
                        vafs[key] = vaf
                    except:
                        pass
    except Exception as e:
        print(f"Error reading VAF from {vcf_file}: {e}")
    return vafs

def main():
    truth_vcf = 'truth_set/truth_set.vcf.gz'
    gatk_vcf = 'results/gatk/variants.vcf.gz'
    deepvariant_vcf = 'results/deepvariant/variants.vcf.gz'

    # Check which files exist
    vcfs_exist = {
        'truth': os.path.exists(truth_vcf),
        'gatk': os.path.exists(gatk_vcf),
        'deepvariant': os.path.exists(deepvariant_vcf)
    }

    print("=" * 70)
    print("VARIANT CALLER COMPARISON REPORT")
    print("=" * 70)
    print()

    print("VCF Files Status:")
    for name, exists in vcfs_exist.items():
        status = "✓ EXISTS" if exists else "✗ MISSING"
        print(f"  {name:15} {status}")
    print()

    if not vcfs_exist['truth']:
        print("ERROR: Truth set VCF not found!")
        return

    # Get truth set variants
    truth_vars = get_variants_from_vcf(truth_vcf)
    print(f"Truth Set: {len(truth_vars)} variants")
    for v in sorted(truth_vars):
        print(f"  {v}")
    print()

    # Compare each caller
    results = {}
    for caller_name, vcf_path in [('GATK', gatk_vcf), ('DeepVariant', deepvariant_vcf)]:
        if not vcfs_exist[caller_name.lower()]:
            results[caller_name] = {
                'file_status': 'MISSING',
                'detected': 0,
                'total': len(truth_vars)
            }
            continue

        caller_vars = get_variants_from_vcf(vcf_path)
        caller_vafs = get_vaf_from_vcf(vcf_path)
        
        detected = truth_vars & caller_vars  # Intersection
        
        results[caller_name] = {
            'file_status': 'EXISTS',
            'detected': len(detected),
            'total': len(truth_vars),
            'sensitivity': len(detected) / len(truth_vars) if truth_vars else 0,
            'variants': caller_vars,
            'vafs': caller_vafs
        }

    # Print results
    print("=" * 70)
    print("COMPARISON RESULTS")
    print("=" * 70)
    print()

    for caller_name in ['GATK', 'DeepVariant']:
        result = results[caller_name]
        print(f"{caller_name}:")
        
        if result['file_status'] == 'MISSING':
            print(f"  Status: VCF FILE NOT YET AVAILABLE")
        else:
            detected = result['detected']
            total = result['total']
            sensitivity = result['sensitivity']
            
            print(f"  Detected: {detected}/{total} variants")
            print(f"  Sensitivity: {sensitivity:.1%}")
            print(f"  VCF: {['gatk', 'deepvariant'][caller_name == 'DeepVariant']}/variants.vcf.gz")
            
            if result['vafs']:
                print(f"  Sample VAF values:")
                for var, vaf in list(result['vafs'].items())[:3]:
                    print(f"    {var}: {vaf:.3f}")
        print()

    print("=" * 70)
    print("NOTE: Run this script again once both callers complete to see full results")
    print("=" * 70)

if __name__ == '__main__':
    main()
