#!/usr/bin/env python3
"""
Convert variant_extraction.py TSV output to standard VCF format.
Aggregates by (chr, pos, ref, alt) and calculates metrics.
"""

import sys
import os
from collections import defaultdict

def parse_extraction_tsv(tsv_file):
    """
    Parse TSV from variant_extraction.py.
    Expected columns: Sample, Gene, Variant, Chromosome, Position, RefBase, AltBase,
                     CellBarcode, ReadID, MutReadCount, WTReadCount, ...
    Returns dict of (chr, pos, ref, alt) -> list of records
    """
    variants = defaultdict(list)

    if not os.path.exists(tsv_file):
        print(f"Warning: {tsv_file} not found", file=sys.stderr)
        return variants

    with open(tsv_file, 'r') as f:
        header = f.readline().strip().split('\t')

        # Find column indices
        try:
            chr_idx = header.index('Chromosome')
            pos_idx = header.index('Position')
            ref_idx = header.index('RefBase')
            alt_idx = header.index('AltBase')
            gene_idx = header.index('Gene')
            var_idx = header.index('Variant')
            mut_count_idx = header.index('MutReadCount')
            wt_count_idx = header.index('WTReadCount')
        except ValueError as e:
            print(f"Error: Required column not found in TSV: {e}", file=sys.stderr)
            return variants

        for line in f:
            if not line.strip():
                continue
            fields = line.strip().split('\t')

            try:
                chrom = fields[chr_idx]
                pos = int(fields[pos_idx])
                ref = fields[ref_idx]
                alt = fields[alt_idx]
                gene = fields[gene_idx]
                variant = fields[var_idx]
                mut_count = int(fields[mut_count_idx])
                wt_count = int(fields[wt_count_idx])

                key = (chrom, pos, ref, alt)
                variants[key].append({
                    'gene': gene,
                    'variant': variant,
                    'mut_reads': mut_count,
                    'wt_reads': wt_count
                })
            except (ValueError, IndexError) as e:
                print(f"Warning: Could not parse line: {e}", file=sys.stderr)
                continue

    return variants

def aggregate_variants(variants_dict):
    """
    Aggregate variant reads and calculate metrics.
    Returns list of VCF records.
    """
    vcf_records = []

    for (chrom, pos, ref, alt), records in variants_dict.items():
        # Aggregate reads
        total_mut = sum(r['mut_reads'] for r in records)
        total_wt = sum(r['wt_reads'] for r in records)
        total_depth = total_mut + total_wt

        if total_depth == 0:
            continue

        vaf = total_mut / total_depth

        # Get gene and variant info (use first record)
        gene = records[0]['gene']
        variant = records[0]['variant']
        ncells = len(set((chrom, pos, ref, alt) for r in records))

        # Build VCF record
        var_id = f"{gene}_{variant}"
        qual = min(100, max(20, int(vaf * 100)))  # Quality proportional to VAF

        info = f"DP={total_depth};MUT={total_mut};WT={total_wt};Gene={gene};Variant={variant};NCells={len(records)}"

        # FORMAT: GT:DP:AD:VAF
        ad = f"{total_wt},{total_mut}"
        format_str = "GT:DP:AD:VAF"
        sample_str = f"0/1:{total_depth}:{ad}:{vaf:.3f}"

        vcf_line = {
            'chrom': chrom,
            'pos': pos,
            'id': var_id,
            'ref': ref,
            'alt': alt,
            'qual': qual,
            'filter': '.',
            'info': info,
            'format': format_str,
            'sample': sample_str
        }

        vcf_records.append(vcf_line)

    return vcf_records

def write_vcf(records, output_file, reference_fasta):
    """Write VCF file with header."""
    header = f"""##fileformat=VCFv4.2
##source=tsv_to_vcf.py
##reference={reference_fasta}
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##INFO=<ID=MUT,Number=1,Type=Integer,Description="Mutant reads">
##INFO=<ID=WT,Number=1,Type=Integer,Description="Wild-type reads">
##INFO=<ID=Gene,Number=1,Type=String,Description="Gene name">
##INFO=<ID=Variant,Number=1,Type=String,Description="Variant annotation">
##INFO=<ID=NCells,Number=1,Type=Integer,Description="Number of cells with variant">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">
##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t5801-diagnosis
"""

    # Sort records
    records = sorted(records, key=lambda r: (r['chrom'], r['pos']))

    with open(output_file, 'w') as f:
        f.write(header)
        for record in records:
            vcf_line = f"{record['chrom']}\t{record['pos']}\t{record['id']}\t{record['ref']}\t{record['alt']}\t{record['qual']}\t{record['filter']}\t{record['info']}\t{record['format']}\t{record['sample']}\n"
            f.write(vcf_line)

    print(f"Wrote {len(records)} variants to {output_file}", file=sys.stderr)

def main():
    if len(sys.argv) < 2:
        print("Usage: tsv_to_vcf.py <tsv_file> [output_vcf]", file=sys.stderr)
        sys.exit(1)

    tsv_file = sys.argv[1]
    output_vcf = sys.argv[2] if len(sys.argv) > 2 else tsv_file.replace('.tsv', '.vcf')
    reference_fa = '/data/salomonis-archive/genomes/hg38/genome.fa'

    print(f"Parsing {tsv_file}...", file=sys.stderr)
    variants = parse_extraction_tsv(tsv_file)

    print(f"Aggregating {len(variants)} unique variants...", file=sys.stderr)
    vcf_records = aggregate_variants(variants)

    print(f"Writing to {output_vcf}...", file=sys.stderr)
    write_vcf(vcf_records, output_vcf, reference_fa)

if __name__ == '__main__':
    main()
