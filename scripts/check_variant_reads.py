#!/usr/bin/env python3
"""
Check how many reads in the BAM actually carry each truth variant allele.
"""

import pysam
import sys
from collections import defaultdict

def check_variant_reads(bam_file, truth_positions):
    """
    Check variant alleles in reads at truth positions.

    Args:
        bam_file: Path to BAM file
        truth_positions: List of (chrom, pos, ref, alt, gene) tuples (0-indexed)
    """

    print(f"Analyzing {len(truth_positions)} truth positions in {bam_file}")
    print()

    results = {}

    for chrom, pos_1based, ref, alt, gene in truth_positions:
        pos = pos_1based - 1  # Convert to 0-indexed

        reads_total = 0
        reads_with_ref = 0
        reads_with_alt = 0
        reads_with_other = 0
        base_counts = defaultdict(int)

        try:
            # Fetch reads at this position
            for read in pysam.AlignmentFile(bam_file, "rb").fetch(chrom, pos, pos + 1):
                reads_total += 1

                # Get the base at this position in the read
                if read.is_unmapped or read.query_alignment_start is None:
                    continue

                # Get position within read sequence
                read_pos = read.get_reference_positions(full_length=True)

                try:
                    # Find which position in the read corresponds to our genomic position
                    query_pos = None
                    for i, ref_pos in enumerate(read_pos):
                        if ref_pos == pos:
                            query_pos = i
                            break

                    if query_pos is None:
                        reads_with_other += 1
                        continue

                    # Get the base at this position
                    base = read.seq[query_pos]
                    base_counts[base] += 1

                    if base == ref:
                        reads_with_ref += 1
                    elif base == alt:
                        reads_with_alt += 1
                    else:
                        reads_with_other += 1

                except (IndexError, TypeError):
                    reads_with_other += 1

        except Exception as e:
            print(f"Error processing {chrom}:{pos_1based}: {e}")
            continue

        # Calculate percentages
        pct_ref = (reads_with_ref / reads_total * 100) if reads_total > 0 else 0
        pct_alt = (reads_with_alt / reads_total * 100) if reads_total > 0 else 0

        results[gene] = {
            'position': f"{chrom}:{pos_1based}",
            'total_reads': reads_total,
            'ref_reads': reads_with_ref,
            'alt_reads': reads_with_alt,
            'other_reads': reads_with_other,
            'pct_ref': pct_ref,
            'pct_alt': pct_alt,
            'ref_allele': ref,
            'alt_allele': alt,
            'base_distribution': dict(base_counts)
        }

    return results

if __name__ == '__main__':
    bam_file = '/data/salomonis-archive/FASTQs/NCI-R01/rna_seq_varcall/results/scisoseq_truth_regions.bam'

    # Truth positions (1-based, as in the truth set)
    truth_positions = [
        ('chr21', 34799432, 'C', 'T', 'RUNX1'),
        ('chr20', 32434638, 'A', 'G', 'ASXL1'),
        ('chr18', 44951948, 'G', 'A', 'SETBP1'),
        ('chr17', 76736877, 'G', 'C', 'SRSF2'),
        ('chr11', 119242488, 'N', 'N', 'CBL'),  # Frameshift - skip
        ('chr19', 19145882, 'N', 'N', 'MEF2B'),  # Unknown - skip
        ('chrX', 40057210, 'N', 'N', 'BCOR'),  # Unknown - skip
        ('chrX', 134415107, 'N', 'N', 'PHF6'),  # Unknown - skip
    ]

    # Filter to known positions only
    known_positions = [p for p in truth_positions if p[2] != 'N']

    results = check_variant_reads(bam_file, known_positions)

    # Print results
    print("=" * 100)
    print(f"{'Gene':<12} {'Position':<20} {'Total':<8} {'REF':<8} {'ALT':<8} {'Other':<8} {'VAF%':<8} {'Alleles':<30}")
    print("=" * 100)

    total_all = 0
    alt_all = 0

    for gene in sorted(results.keys(), key=lambda x: results[x]['position']):
        r = results[gene]
        total_all += r['total_reads']
        alt_all += r['alt_reads']

        vaf = (r['alt_reads'] / r['total_reads'] * 100) if r['total_reads'] > 0 else 0

        # Show base distribution
        bases = ', '.join([f"{b}:{c}" for b, c in sorted(r['base_distribution'].items())])

        print(f"{gene:<12} {r['position']:<20} {r['total_reads']:<8} {r['ref_reads']:<8} {r['alt_reads']:<8} {r['other_reads']:<8} {vaf:<8.1f} {bases:<30}")

    print("=" * 100)
    print(f"{'TOTAL':<12} {'':<20} {total_all:<8} {'':<8} {alt_all:<8} {'':<8} {(alt_all/total_all*100 if total_all > 0 else 0):<8.1f}")
    print()
    print(f"Summary: {alt_all:,} reads with alt allele out of {total_all:,} total reads ({alt_all/total_all*100:.1f}%)")
