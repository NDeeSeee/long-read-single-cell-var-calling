#!/usr/bin/env python3
"""
Add synthetic quality scores to BAM that has empty QUAL field.

Creates reasonable quality scores (Q30 = high confidence) based on read length.
This allows tools like GATK to process the BAM.
"""

import sys
import pysam
from pathlib import Path

def add_synthetic_quality(input_bam, output_bam, default_quality=30, region=None):
    """
    Add synthetic quality scores to BAM with empty QUAL field.

    Args:
        input_bam: Path to input BAM
        output_bam: Path to output BAM with synthetic quality
        default_quality: Default Phred quality score to assign (0-60)
        region: Optional region string
    """

    print(f"Opening input BAM: {input_bam}")
    infile = pysam.AlignmentFile(input_bam, "rb")

    print(f"Creating output BAM: {output_bam}")
    outfile = pysam.AlignmentFile(output_bam, "wb", template=infile)

    processed = 0
    updated = 0

    try:
        iterator = infile.fetch(region=region) if region else infile

        for read in iterator:
            processed += 1

            # Check if QUAL is empty or missing
            if read.query_qualities is None or len(read.query_qualities) == 0:
                # Create synthetic quality array
                # Use ASCII character for quality score (e.g., '>' = Q30)
                qual_char = chr(33 + default_quality)  # ASCII 33 = '!', +30 = '>'
                synthetic_qual = qual_char * len(read.seq)

                # Set the quality
                read.query_qualities = pysam.qualitystring_to_array(synthetic_qual)
                updated += 1

            # Write read
            outfile.write(read)

            if processed % 100000 == 0:
                print(f"  Processed {processed:,} reads ({updated:,} updated)")

        print(f"\nCompleted!")
        print(f"  Total reads: {processed:,}")
        print(f"  Reads updated: {updated:,}")
        print(f"  Default quality: Q{default_quality}")

    finally:
        infile.close()
        outfile.close()

    # Index output
    print(f"Indexing output BAM...")
    pysam.index(output_bam)
    print(f"Done! Output indexed at {output_bam}.bai")

    return processed, updated

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python add_synthetic_quality.py <input.bam> <output.bam> [quality]")
        print("\nExample:")
        print("  python add_synthetic_quality.py input.bam output_qual.bam 30")
        print("  python add_synthetic_quality.py input.bam output_qual.bam 30 chr21:1-50000000")
        sys.exit(1)

    input_bam = sys.argv[1]
    output_bam = sys.argv[2]
    default_quality = int(sys.argv[3]) if len(sys.argv) > 3 else 30
    region = sys.argv[4] if len(sys.argv) > 4 else None

    # Verify input
    if not Path(input_bam).exists():
        print(f"ERROR: Input BAM not found: {input_bam}")
        sys.exit(1)

    if not (0 <= default_quality <= 60):
        print(f"ERROR: Quality must be 0-60, got {default_quality}")
        sys.exit(1)

    # Add synthetic quality
    processed, updated = add_synthetic_quality(input_bam, output_bam, default_quality, region)

    if updated == 0:
        print("\nWARNING: No reads were updated.")
        sys.exit(1)

    print(f"\nSUCCESS: Output saved to {output_bam}")
