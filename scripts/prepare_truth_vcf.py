#!/usr/bin/env python3
"""
Convert test_variants.txt to VCF format for variant comparison.
"""

import subprocess

def main():
    input_file = 'test_variants.txt'
    output_vcf = 'truth_set/truth_set.vcf'

    # Read truth set
    variants = []
    with open(input_file) as f:
        next(f)  # skip header
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 6:
                continue
            chrom, pos, end, gene, variant, vaf = fields
            variants.append({
                'chrom': chrom,
                'pos': int(pos),
                'gene': gene,
                'variant': variant,
                'vaf': vaf
            })

    # Sort by chromosome and position
    # Define chromosome order
    chrom_order = {f'chr{i}': i for i in range(1, 23)}
    chrom_order['chrX'] = 23
    chrom_order['chrY'] = 24
    chrom_order['chrM'] = 25

    variants.sort(key=lambda v: (chrom_order.get(v['chrom'], 999), v['pos']))

    # Write VCF
    with open(output_vcf, 'w') as f:
        # VCF header
        f.write('##fileformat=VCFv4.2\n')
        f.write('##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">\n')
        f.write('##INFO=<ID=VAR,Number=1,Type=String,Description="Variant annotation">\n')
        f.write('##INFO=<ID=VAF,Number=1,Type=String,Description="Variant allele frequency">\n')
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

        # Use placeholder REF/ALT
        for v in variants:
            info = f'GENE={v["gene"]};VAR={v["variant"]};VAF={v["vaf"]}'
            f.write(f'{v["chrom"]}\t{v["pos"]}\t.\tN\t.\t.\t.\t{info}\n')

    # Compress and index
    subprocess.run(f'bgzip -f {output_vcf}', shell=True, check=True)
    subprocess.run(f'bcftools index -t {output_vcf}.gz', shell=True, check=True)

    print(f'Created truth set VCF: {output_vcf}.gz')
    print(f'Total variants: {len(variants)}')

if __name__ == '__main__':
    main()
