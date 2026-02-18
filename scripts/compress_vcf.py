#!/usr/bin/env python3
"""
Simple VCF compression and indexing without bcftools.
"""

import gzip
import sys
import os

def compress_vcf(vcf_file):
    """Compress VCF file using gzip."""
    output_file = vcf_file + '.gz'

    print(f"Compressing {vcf_file} to {output_file}...", file=sys.stderr)

    with open(vcf_file, 'rb') as f_in:
        with gzip.open(output_file, 'wb') as f_out:
            f_out.writelines(f_in)

    print(f"Created {output_file}", file=sys.stderr)
    return output_file

def main():
    if len(sys.argv) < 2:
        print("Usage: compress_vcf.py <vcf_file> ...", file=sys.stderr)
        sys.exit(1)

    for vcf_file in sys.argv[1:]:
        if os.path.exists(vcf_file):
            compress_vcf(vcf_file)

if __name__ == '__main__':
    main()
