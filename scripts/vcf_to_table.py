#!/usr/bin/env python3
"""
vcf_to_table.py

Convert an annotated _final.vcf.gz to a human-readable TSV.

Extracts per-variant (canonical transcript preferred):
  CHROM, POS, REF, ALT, GENE, CONSEQUENCE, HGVSp, HGVSc,
  COSMIC, CLINVAR, gnomAD_AF_max,
  SOMATIC_CLASS, SOMATIC_REASON,
  BAM_REF, BAM_ALT, BAM_VAF, BAM_DP,
  CALLER_COUNT, CALLERS   (merged files only)

Usage:
    python scripts/vcf_to_table.py \\
        --input  annotation/gatk/non_germline_final.vcf.gz \\
        --output annotation/gatk/non_germline_final.tsv

    python scripts/vcf_to_table.py \\
        --input  results/merged/high_confidence_final.vcf.gz \\
        --output results/merged/high_confidence_final.tsv
"""

import argparse
import gzip
import re
import sys
from pathlib import Path

IMPACT_RANK = {'HIGH': 3, 'MODERATE': 2, 'LOW': 1, 'MODIFIER': 0}
COSMIC_RE   = re.compile(r'COSV\d+|COSM\d+')

# Within MODIFIER, rank by how close/inside the gene the variant is.
CONSEQUENCE_SUBRANK = {
    'non_coding_transcript_exon_variant': 5,
    'non_coding_transcript_variant':      4,
    '3_prime_UTR_variant':                3,
    '5_prime_UTR_variant':                3,
    'intron_variant':                     2,
    'upstream_gene_variant':              1,
    'downstream_gene_variant':            1,
    'regulatory_region_variant':          1,
    'intergenic_variant':                 0,
}


def _csq_rank(impact_str, consequence_str):
    imp = IMPACT_RANK.get(impact_str, 0)
    sub = CONSEQUENCE_SUBRANK.get(consequence_str.split('&')[0], 0) if imp == 0 else 0
    return (imp, sub)


def open_vcf(path):
    return gzip.open(str(path), 'rt') if str(path).endswith('.gz') else open(path)


def parse_csq_header(header_lines):
    for line in header_lines:
        if 'ID=CSQ' in line:
            m = re.search(r'Format: ([^"]+)"', line)
            if m:
                return m.group(1).split('|')
    return []


def best_csq(csq_string, fields):
    """
    Return the single best CSQ entry as a dict.
    Prefers canonical + highest (IMPACT, consequence_subrank);
    falls back to any highest-rank entry.
    """
    if not csq_string or not fields:
        return {}

    can_idx  = fields.index('CANONICAL')   if 'CANONICAL'   in fields else -1
    imp_idx  = fields.index('IMPACT')      if 'IMPACT'      in fields else -1
    cons_idx = fields.index('Consequence') if 'Consequence' in fields else -1

    best_canonical = (None, (-1, -1))
    best_any       = (None, (-1, -1))

    for entry in csq_string.split(','):
        vals = entry.split('|')
        d    = dict(zip(fields, vals + [''] * max(0, len(fields) - len(vals))))
        imp_str  = d.get('IMPACT', '')
        cons_str = vals[cons_idx] if cons_idx != -1 and cons_idx < len(vals) else ''
        rank = _csq_rank(imp_str, cons_str)
        is_can = can_idx != -1 and can_idx < len(vals) and vals[can_idx] == 'YES'
        if rank > best_any[1]:
            best_any = (d, rank)
        if is_can and rank > best_canonical[1]:
            best_canonical = (d, rank)

    chosen = best_canonical[0] if best_canonical[0] else best_any[0]
    return chosen or {}


def get_info(info_str, key):
    """Extract a single INFO key value."""
    m = re.search(r'(?:^|;)' + re.escape(key) + r'=([^;]+)', info_str)
    return m.group(1) if m else ''


def cosmic_ids(existing_variation):
    """Pull COSV/COSM IDs from Existing_variation field."""
    return '&'.join(COSMIC_RE.findall(existing_variation)) or ''


def max_gnomad_af(csq_string, fields):
    """Max gnomAD AF across all transcripts and populations."""
    af_fields = [f for f in fields if 'gnomAD' in f and f.endswith('AF')]
    max_af = 0.0
    for entry in csq_string.split(','):
        vals = entry.split('|')
        d    = dict(zip(fields, vals + [''] * max(0, len(fields) - len(vals))))
        for f in af_fields:
            v = d.get(f, '')
            if v and v != '.':
                try:
                    max_af = max(max_af, float(v))
                except ValueError:
                    pass
    return f'{max_af:.6g}' if max_af > 0 else '0'


COLUMNS = [
    'CHROM', 'POS', 'REF', 'ALT',
    'GENE', 'CONSEQUENCE', 'HGVSp', 'HGVSc',
    'COSMIC', 'CLINVAR', 'gnomAD_AF_max',
    'SOMATIC_CLASS', 'SOMATIC_REASON',
    'BAM_REF', 'BAM_ALT', 'BAM_VAF', 'BAM_DP',
    'CALLER_COUNT', 'CALLERS',
]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--input',  required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    header_lines = []
    csq_fields   = []
    rows = []

    with open_vcf(args.input) as fh:
        for line in fh:
            if line.startswith('#'):
                header_lines.append(line)
                if not csq_fields:
                    csq_fields = parse_csq_header([line])
                continue

            cols = line.rstrip('\n').split('\t')
            if len(cols) < 8:
                continue

            chrom, pos, _, ref, alt = cols[:5]
            info = cols[7]

            # CSQ
            csq_raw = ''
            for field in info.split(';'):
                if field.startswith('CSQ='):
                    csq_raw = field[4:]
                    break

            d = best_csq(csq_raw, csq_fields)

            consequence = d.get('Consequence', '').split('&')[0]   # first effect only
            hgvsp       = d.get('HGVSp',  '')
            hgvsc       = d.get('HGVSc',  '')
            existing    = d.get('Existing_variation', '')
            clin_sig    = d.get('CLIN_SIG', '')

            row = {
                'CHROM':         chrom,
                'POS':           pos,
                'REF':           ref,
                'ALT':           alt.split(',')[0],
                'GENE':          get_info(info, 'GENE'),
                'CONSEQUENCE':   consequence,
                'HGVSp':         hgvsp,
                'HGVSc':         hgvsc,
                'COSMIC':        cosmic_ids(existing),
                'CLINVAR':       clin_sig,
                'gnomAD_AF_max': max_gnomad_af(csq_raw, csq_fields) if csq_raw else '0',
                'SOMATIC_CLASS': get_info(info, 'SOMATIC_CLASS'),
                'SOMATIC_REASON':get_info(info, 'SOMATIC_REASON'),
                'BAM_REF':       get_info(info, 'BAM_REF'),
                'BAM_ALT':       get_info(info, 'BAM_ALT'),
                'BAM_VAF':       get_info(info, 'BAM_VAF'),
                'BAM_DP':        get_info(info, 'BAM_DP'),
                'CALLER_COUNT':  get_info(info, 'CALLER_COUNT'),
                'CALLERS':       get_info(info, 'CALLERS'),
            }
            rows.append(row)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, 'w') as fh:
        fh.write('\t'.join(COLUMNS) + '\n')
        for row in rows:
            fh.write('\t'.join(row.get(c, '') for c in COLUMNS) + '\n')

    print(f'Wrote {len(rows):,} variants → {out}')


if __name__ == '__main__':
    main()
