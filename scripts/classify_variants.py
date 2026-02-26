#!/usr/bin/env python3
"""
classify_variants.py

Reads a VEP-annotated VCF and classifies each variant into:
  - SOMATIC_KNOWN     : in COSMIC / CancerHotspots / ClinVar Pathogenic
  - NOVEL_CANDIDATE   : absent from gnomAD entirely, QUAL >= 30
  - RARE_CANDIDATE    : gnomAD AF < 0.0001 (1 in 10,000), not confirmed germline
  - GERMLINE          : gnomAD AF >= 0.001 in any population

Decision priority:  SOMATIC_KNOWN > NOVEL_CANDIDATE > RARE_CANDIDATE > GERMLINE
"""

import argparse
import re
import gzip
import subprocess
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# AF thresholds
# ---------------------------------------------------------------------------
GERMLINE_AF    = 0.001    # gnomAD AF >= 0.1%  → exclude as germline
RARE_AF        = 0.0001   # gnomAD AF < 0.01%  → rare candidate
MIN_QUAL       = 30       # minimum QUAL for novel candidates
MIN_ALT_DEPTH  = 5        # minimum ALT supporting reads (FORMAT AD field)

# ClinVar significance terms that indicate pathogenicity
CLINVAR_PATHOGENIC = {
    "pathogenic", "likely_pathogenic",
    "pathogenic/likely_pathogenic",
}

# Known somatic databases annotated in VEP CSQ field
COSMIC_PATTERN       = re.compile(r'COSV\d+|COSM\d+')
CANCERHOTSPOT_PATTERN = re.compile(r'cancer_hotspot', re.IGNORECASE)


def parse_csq_header(header_lines):
    """Extract CSQ field order from VEP header."""
    for line in header_lines:
        if 'ID=CSQ' in line:
            m = re.search(r'Format: ([^"]+)"', line)
            if m:
                return m.group(1).split('|')
    return []


def parse_csq(csq_string, fields):
    """Parse a single CSQ entry into a dict."""
    vals = csq_string.split('|')
    return dict(zip(fields, vals + [''] * (len(fields) - len(vals))))


def get_max_gnomad_af(csq_entries, fields):
    """Return the maximum gnomAD AF across all transcripts and populations."""
    af_fields = [f for f in fields if 'gnomAD' in f and 'AF' in f]
    max_af = 0.0
    for entry in csq_entries:
        d = parse_csq(entry, fields)
        for f in af_fields:
            val = d.get(f, '')
            if val and val != '.':
                try:
                    max_af = max(max_af, float(val))
                except ValueError:
                    pass
    return max_af


def is_somatic_known(csq_entries, fields, existing_var):
    """
    Check for somatic evidence:
    - Existing variant ID matches COSMIC pattern
    - ClinVar significance is pathogenic
    """
    # Check existing variant IDs (dbSNP, COSMIC, etc.)
    if existing_var and COSMIC_PATTERN.search(existing_var):
        return True, "COSMIC"

    for entry in csq_entries:
        d = parse_csq(entry, fields)

        # ClinVar
        clin_sig = d.get('CLIN_SIG', '').lower()
        if any(term in clin_sig for term in CLINVAR_PATHOGENIC):
            return True, f"ClinVar:{d.get('CLIN_SIG','')}"

        # COSMIC via existing variant annotation
        existing = d.get('Existing_variation', '')
        if existing and COSMIC_PATTERN.search(existing):
            return True, "COSMIC"

    return False, ""


def classify_variant(qual, gnomad_af, somatic, somatic_source, alt_depth=0):
    """Apply tiered decision logic."""
    # Low-support calls are unreliable regardless of annotation.
    # alt_depth == -1 means no depth field present — skip filter.
    if alt_depth != -1 and alt_depth < MIN_ALT_DEPTH:
        return "GERMLINE", f"low_depth_AD={alt_depth}"
    if somatic:
        # COSMIC always overrides gnomAD — somatic mutations don't appear in population DBs
        # ClinVar pathogenic alone defers to gnomAD: ClinVar annotates germline Mendelian
        # variants (BRCA1, CFTR…) that are common in the population and NOT somatic.
        is_cosmic = 'COSMIC' in somatic_source
        if is_cosmic or gnomad_af < GERMLINE_AF:
            return "SOMATIC_KNOWN", somatic_source
    if gnomad_af >= GERMLINE_AF:
        return "GERMLINE", f"gnomAD_AF={gnomad_af:.4f}"
    if gnomad_af == 0.0 and qual >= MIN_QUAL:
        return "NOVEL_CANDIDATE", "not_in_gnomAD"
    if gnomad_af < RARE_AF:
        return "RARE_CANDIDATE", f"gnomAD_AF={gnomad_af:.6f}"
    return "UNCERTAIN", f"gnomAD_AF={gnomad_af:.4f}"


def get_alt_depth(cols):
    """
    Parse ALT allele depth from FORMAT fields.
    Priority:
      1. AD field (REF,ALT counts) — GATK, DeepVariant, Clair3
      2. DP * AF                   — LongcallR (GT:GQ:PS:DP:AF:PQ)
    Returns -1 if no depth info is available (depth filter skipped).
    """
    if len(cols) < 10:
        return -1
    fmt    = cols[8].split(':')
    sample = cols[9].split(':')

    def field(name):
        if name not in fmt:
            return None
        idx = fmt.index(name)
        return sample[idx] if idx < len(sample) else None

    # Option 1: AD field
    ad_raw = field('AD')
    if ad_raw and ad_raw not in ('.', ''):
        ad_vals = ad_raw.split(',')
        if len(ad_vals) >= 2:
            try:
                return int(ad_vals[1])
            except ValueError:
                pass

    # Option 2: DP * AF (LongcallR)
    dp_raw = field('DP')
    af_raw = field('AF')
    if dp_raw and af_raw and dp_raw not in ('.', '') and af_raw not in ('.', ''):
        try:
            return int(float(dp_raw) * float(af_raw))
        except ValueError:
            pass

    return -1   # no depth info — skip depth filter


def open_vcf(path):
    p = str(path)
    if p.endswith('.gz'):
        return gzip.open(p, 'rt')
    return open(p, 'r')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',  required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--sample', required=True)
    args = parser.parse_args()

    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Output file handles (uncompressed, will bgzip+index after)
    buckets = {
        'SOMATIC_KNOWN':   open(out_dir / 'somatic_known.vcf',   'w'),
        'NOVEL_CANDIDATE': open(out_dir / 'novel_candidates.vcf', 'w'),
        'RARE_CANDIDATE':  open(out_dir / 'rare_candidates.vcf',  'w'),
        'GERMLINE':        open(out_dir / 'germline.vcf',         'w'),
        'UNCERTAIN':       open(out_dir / 'uncertain.vcf',        'w'),
    }

    counts = {k: 0 for k in buckets}
    header_lines = []
    csq_fields = []

    SOMATIC_INFO_HEADERS = (
        '##INFO=<ID=SOMATIC_CLASS,Number=1,Type=String,'
        'Description="Somatic classification tier: SOMATIC_KNOWN, NOVEL_CANDIDATE, '
        'RARE_CANDIDATE, GERMLINE, or UNCERTAIN">\n'
        '##INFO=<ID=SOMATIC_REASON,Number=1,Type=String,'
        'Description="Reason for SOMATIC_CLASS assignment">\n'
    )

    with open_vcf(args.input) as fh:
        for line in fh:
            if line.startswith('#'):
                header_lines.append(line)
                if not csq_fields:
                    csq_fields = parse_csq_header([line])
                # Insert SOMATIC INFO headers just before the #CHROM line
                if line.startswith('#CHROM'):
                    for fout in buckets.values():
                        fout.write(SOMATIC_INFO_HEADERS)
                for fout in buckets.values():
                    fout.write(line)
                continue

            # Parse variant line
            cols = line.rstrip('\n').split('\t')
            chrom, pos, vid, ref, alt, qual_str = cols[:6]
            info = cols[7] if len(cols) > 7 else ''

            try:
                qual = float(qual_str) if qual_str not in ('.', '') else 0.0
            except ValueError:
                qual = 0.0

            # Extract CSQ entries
            csq_entries = []
            for field in info.split(';'):
                if field.startswith('CSQ='):
                    csq_entries = field[4:].split(',')
                    break

            # Existing variant IDs (dbSNP rs, COSMIC, etc.)
            existing_var = vid  # VID column often has rs/COSV IDs

            gnomad_af  = get_max_gnomad_af(csq_entries, csq_fields)
            alt_depth  = get_alt_depth(cols)
            somatic, src = is_somatic_known(csq_entries, csq_fields, existing_var)
            category, reason = classify_variant(qual, gnomad_af, somatic, src, alt_depth)

            # Append classification to INFO
            cols[7] = info + f';SOMATIC_CLASS={category};SOMATIC_REASON={reason}'
            buckets[category].write('\t'.join(cols) + '\n')
            counts[category] += 1

    for fout in buckets.values():
        fout.close()

    # Compress and index each output
    # Input is already sorted (VEP preserves order), so bgzip directly.
    # bcftools sort crashes on header-only (empty) VCFs, so we skip it.
    bgzip   = '/users/pavb5f/.conda/envs/bio-cli/bin/bgzip'
    bcftools = '/users/pavb5f/.conda/envs/bio-cli/bin/bcftools'
    for category, vcf in [
        ('SOMATIC_KNOWN',   'somatic_known.vcf'),
        ('NOVEL_CANDIDATE', 'novel_candidates.vcf'),
        ('RARE_CANDIDATE',  'rare_candidates.vcf'),
        ('GERMLINE',        'germline.vcf'),
        ('UNCERTAIN',       'uncertain.vcf'),
    ]:
        vcf_path = out_dir / vcf
        gz_path  = out_dir / (vcf + '.gz')
        # bgzip in-place (-f overwrite if exists)
        subprocess.run([bgzip, '-f', str(vcf_path)],
                       check=True, capture_output=True)
        # Only index if file has variants (non-trivial size)
        if gz_path.stat().st_size > 1024:
            subprocess.run([bcftools, 'index', '-t', str(gz_path)],
                           check=True, capture_output=True)

    # Summary
    print("\n=== Classification Summary ===")
    total = sum(counts.values())
    for cat, n in sorted(counts.items(), key=lambda x: -x[1]):
        pct = 100 * n / total if total else 0
        print(f"  {cat:<20} {n:>8,}  ({pct:.1f}%)")
    print(f"  {'TOTAL':<20} {total:>8,}")
    print(f"\nOutputs written to: {out_dir}/")


if __name__ == '__main__':
    main()
