#!/usr/bin/env python3
"""
Extract PacBio quality scores from qs:Z: tag and populate standard QUAL field.

PacBio Iso-Seq BAMs store quality as ASCII-encoded string in qs:Z: auxiliary tag
rather than standard SAM QUAL field. This script converts them for compatibility
with tools like GATK that expect standard QUAL format.
"""

import sys
import pysam
from pathlib import Path

def extract_pacbio_quality(input_bam, output_bam, region=None):
    """
    Extract quality scores from qs:Z: tag and populate QUAL field.

    Args:
        input_bam: Path to input BAM with qs:Z: tag
        output_bam: Path to output BAM with QUAL field populated
        region: Optional region string (e.g., 'chr21:1-1000000')
    """

    print(f"Opening input BAM: {input_bam}")
    infile = pysam.AlignmentFile(input_bam, "rb")

    print(f"Creating output BAM: {output_bam}")
    outfile = pysam.AlignmentFile(output_bam, "wb", template=infile)

    processed = 0
    with_qs_tag = 0
    updated = 0
    errors = 0

    try:
        # Iterate over reads
        iterator = infile.fetch(region=region) if region else infile

        for read in iterator:
            processed += 1

            # Check if read has qs:Z: tag
            if read.has_tag('qs'):
                with_qs_tag += 1
                try:
                    # Extract quality string
                    qs = read.get_tag('qs')

                    # Verify length matches sequence
                    if len(qs) != len(read.seq):
                        print(f"WARNING: Read {read.query_name} has mismatched lengths:")
                        print(f"  Sequence: {len(read.seq)} bp")
                        print(f"  qs tag: {len(qs)} chars")
                        errors += 1
                        # Still try to use what we have
                        qs = qs[:len(read.seq)]  # Truncate if needed

                    # Set QUAL field from qs:Z: tag
                    read.query_qualities = pysam.qualitystring_to_array(qs)
                    updated += 1

                except Exception as e:
                    print(f"ERROR processing read {read.query_name}: {e}")
                    errors += 1

            # Write read to output
            outfile.write(read)

            if processed % 100000 == 0:
                print(f"  Processed {processed:,} reads ({with_qs_tag:,} with qs tag, {updated:,} updated)")

        print(f"\nCompleted!")
        print(f"  Total reads: {processed:,}")
        print(f"  Reads with qs:Z: tag: {with_qs_tag:,}")
        print(f"  Reads updated: {updated:,}")
        print(f"  Errors: {errors:,}")

    finally:
        infile.close()
        outfile.close()

    # Index output BAM
    print(f"Indexing output BAM...")
    pysam.index(output_bam)
    print(f"Done! Output indexed at {output_bam}.bai")

    return processed, with_qs_tag, updated, errors

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python extract_pacbio_quality.py <input.bam> <output.bam> [region]")
        print("\nExample:")
        print("  python extract_pacbio_quality.py input.bam output_qual_fixed.bam")
        print("  python extract_pacbio_quality.py input.bam output_qual_fixed.bam chr21:1-50000000")
        sys.exit(1)

    input_bam = sys.argv[1]
    output_bam = sys.argv[2]
    region = sys.argv[3] if len(sys.argv) > 3 else None

    # Verify input exists
    if not Path(input_bam).exists():
        print(f"ERROR: Input BAM not found: {input_bam}")
        sys.exit(1)

    # Run extraction
    processed, with_qs, updated, errors = extract_pacbio_quality(input_bam, output_bam, region)

    # Check results
    if updated == 0:
        print("\nWARNING: No reads were updated. Check if BAM has qs:Z: tags.")
        sys.exit(1)

    if errors > 0 and errors > processed * 0.1:  # More than 10% errors
        print("\nWARNING: High error rate detected. Output may be incomplete.")
        sys.exit(1)

    print(f"\nSUCCESS: Output saved to {output_bam}")
