#!/usr/bin/env python3
"""
add_bam_annotations.py

Adds three annotation layers to any classified/annotated VCF:

  1. INFO/GENE       — Hugo gene symbol (from VEP CSQ, canonical transcript first)
  2. INFO/BAM_REF    — raw REF read count directly from BAM (samtools mpileup)
  3. INFO/BAM_ALT    — raw ALT read count directly from BAM
  4. INFO/BAM_VAF    — raw VAF = BAM_ALT / (BAM_REF + BAM_ALT)

BAM counts use samtools mpileup with no quality or mapping-quality filters
(-B -Q 0 -q 0 -A), so MQ=255 PacBio reads excluded by GATK are included.
This gives true per-allele depth consistent with mosdepth BAM_DP.

Usage:
    python scripts/add_bam_annotations.py \\
        --input  annotation/gatk/non_germline_bamdp1500.vcf.gz \\
        --output annotation/gatk/non_germline_final.vcf.gz \\
        --bam    scisoseq_synthetic_qual.bam \\
        [--ref   reference/genome.fa]
"""

import argparse
import gzip
import re
import subprocess
import sys
import tempfile
from pathlib import Path

WORKDIR  = Path('/data/salomonis-archive/FASTQs/NCI-R01/rna_seq_varcall')
SAMTOOLS = '/users/pavb5f/.conda/envs/bio-cli/bin/samtools'
BGZIP    = '/users/pavb5f/.conda/envs/bio-cli/bin/bgzip'
BCFTOOLS = '/users/pavb5f/.conda/envs/bio-cli/bin/bcftools'


# ---------------------------------------------------------------------------
# GENE extraction from VEP CSQ
# ---------------------------------------------------------------------------

def parse_csq_fields(header_lines):
    for line in header_lines:
        if 'ID=CSQ' in line:
            m = re.search(r'Format: ([^"]+)"', line)
            if m:
                return m.group(1).split('|')
    return []


IMPACT_RANK = {'HIGH': 3, 'MODERATE': 2, 'LOW': 1, 'MODIFIER': 0}

# Within MODIFIER, rank by how close/inside the gene the variant is.
# non_coding_transcript_exon_variant (inside gene body) > upstream/downstream.
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


def _rank(impact_str, consequence_str):
    """Return (impact_rank, consequence_subrank) tuple for comparison."""
    imp = IMPACT_RANK.get(impact_str, 0)
    # Sub-rank only acts as tiebreaker within the same impact tier
    sub = CONSEQUENCE_SUBRANK.get(consequence_str.split('&')[0], 0) if imp == 0 else 0
    return (imp, sub)


def extract_gene(csq_string, fields):
    """
    Return Hugo gene symbol from VEP CSQ field.
    Among canonical transcripts, prefer highest IMPACT then consequence specificity.
    Falls back to highest-rank non-canonical entry, then first non-empty SYMBOL.
    """
    if not fields or not csq_string:
        return '.'
    sym_idx  = fields.index('SYMBOL')      if 'SYMBOL'      in fields else -1
    can_idx  = fields.index('CANONICAL')   if 'CANONICAL'   in fields else -1
    imp_idx  = fields.index('IMPACT')      if 'IMPACT'      in fields else -1
    cons_idx = fields.index('Consequence') if 'Consequence' in fields else -1
    if sym_idx == -1:
        return '.'

    best_canonical = ('', (-1, -1))
    best_any       = ('', (-1, -1))
    first_gene     = ''

    for entry in csq_string.split(','):
        vals = entry.split('|')
        if sym_idx >= len(vals):
            continue
        symbol = vals[sym_idx]
        if not symbol:
            continue
        if not first_gene:
            first_gene = symbol
        imp_str  = vals[imp_idx]  if imp_idx  != -1 and imp_idx  < len(vals) else ''
        cons_str = vals[cons_idx] if cons_idx != -1 and cons_idx < len(vals) else ''
        r = _rank(imp_str, cons_str)
        if r > best_any[1]:
            best_any = (symbol, r)
        is_canonical = can_idx != -1 and can_idx < len(vals) and vals[can_idx] == 'YES'
        if is_canonical and r > best_canonical[1]:
            best_canonical = (symbol, r)

    if best_canonical[0]:
        return best_canonical[0]
    if best_any[0]:
        return best_any[0]
    return first_gene or '.'


# ---------------------------------------------------------------------------
# Pileup string parsers
# ---------------------------------------------------------------------------

def _skip_indel(s, i):
    """Advance i past a +N[seq] or -N[seq] annotation. Returns new i."""
    i += 1  # skip +/-
    num = ''
    while i < len(s) and s[i].isdigit():
        num += s[i]
        i += 1
    return i + (int(num) if num else 0)


def count_snv(pileup, alt):
    """Count REF (., ,) and ALT base in pileup string for an SNV."""
    s = re.sub(r'\^.', '', pileup).replace('$', '')
    ref_n = alt_n = 0
    i = 0
    while i < len(s):
        c = s[i]
        if c in '+-':
            i = _skip_indel(s, i)
        elif c in '.,':
            ref_n += 1
            i += 1
        elif c == '*':
            i += 1
        elif c.upper() == alt.upper():
            alt_n += 1
            i += 1
        else:
            i += 1
    return ref_n, alt_n


def count_insertion(pileup, ins_seq):
    """Count reads supporting REF vs insertion ALT."""
    s = re.sub(r'\^.', '', pileup).replace('$', '')
    ref_n = alt_n = 0
    i = 0
    while i < len(s):
        c = s[i]
        i += 1
        if c in '.,ACGTNacgtn':
            if i < len(s) and s[i] in '+-':
                op = s[i]
                j  = i + 1
                num = ''
                while j < len(s) and s[j].isdigit():
                    num += s[j]
                    j += 1
                n = int(num) if num else 0
                seq = s[j:j+n].upper()
                if op == '+' and seq == ins_seq.upper():
                    alt_n += 1
                elif c in '.,':
                    ref_n += 1
                i = j + n
            elif c in '.,':
                ref_n += 1
        elif c == '-':
            num = ''
            while i < len(s) and s[i].isdigit():
                num += s[i]
                i += 1
            i += int(num) if num else 0
        # * = deleted placeholder, skip
    return ref_n, alt_n


def count_deletion(pileup, del_seq):
    """Count reads supporting REF vs deletion ALT."""
    s = re.sub(r'\^.', '', pileup).replace('$', '')
    ref_n = alt_n = 0
    i = 0
    while i < len(s):
        c = s[i]
        i += 1
        if c in '.,ACGTNacgtn':
            if i < len(s) and s[i] == '-':
                j  = i + 1
                num = ''
                while j < len(s) and s[j].isdigit():
                    num += s[j]
                    j += 1
                n   = int(num) if num else 0
                seq = s[j:j+n].upper()
                if seq == del_seq.upper():
                    alt_n += 1
                else:
                    if c in '.,':
                        ref_n += 1
                i = j + n
            elif c in '.,':
                ref_n += 1
        elif c == '+':
            num = ''
            while i < len(s) and s[i].isdigit():
                num += s[i]
                i += 1
            i += int(num) if num else 0
        elif c == '*':
            alt_n += 1  # deleted base marker = deletion support
    return ref_n, alt_n


# ---------------------------------------------------------------------------
# Batch BAM pileup query
# ---------------------------------------------------------------------------

def batch_bam_counts(bam_path, ref_fa, variants):
    """
    variants : list of (chrom, pos[1-based], ref, alt)
    Returns  : dict {(chrom, pos, ref, alt): (bam_ref, bam_alt)}
    """
    if not variants:
        return {}

    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as fh:
        for chrom, pos, ref, alt in variants:
            fh.write(f'{chrom}\t{pos-1}\t{pos}\n')
        bed = fh.name

    cmd = [SAMTOOLS, 'mpileup',
           '-l', bed,
           '-f', str(ref_fa),
           '-B',        # disable BAQ
           '-Q', '0',   # no min base quality
           '-q', '0',   # no min mapping quality (include MQ=255 PacBio reads)
           '-A',        # count orphan reads
           str(bam_path)]

    result = subprocess.run(cmd, capture_output=True, text=True)
    Path(bed).unlink()

    if result.returncode != 0:
        print(f'[WARN] samtools mpileup error: {result.stderr[:200]}', file=sys.stderr)
        return {}

    # Parse pileup output → dict keyed by (chrom, pos)
    pileup_map = {}
    for line in result.stdout.splitlines():
        parts = line.split('\t')
        if len(parts) >= 5:
            pileup_map[(parts[0], int(parts[1]))] = parts[4]  # pileup string

    counts = {}
    for chrom, pos, ref, alt in variants:
        pileup = pileup_map.get((chrom, pos), '')
        if not pileup:
            counts[(chrom, pos, ref, alt)] = (0, 0)
            continue

        # First ALT allele only (handles multi-allelics)
        a = alt.split(',')[0]

        if len(ref) == 1 and len(a) == 1:
            r, a_cnt = count_snv(pileup, a)
        elif len(a) > len(ref):
            r, a_cnt = count_insertion(pileup, a[len(ref):])
        else:
            r, a_cnt = count_deletion(pileup, ref[len(a):])

        counts[(chrom, pos, ref, alt)] = (r, a_cnt)

    return counts


# ---------------------------------------------------------------------------
# New INFO header definitions
# ---------------------------------------------------------------------------

NEW_HEADERS = (
    '##INFO=<ID=GENE,Number=1,Type=String,'
    'Description="Hugo gene symbol from VEP CSQ (canonical transcript preferred)">\n'
    '##INFO=<ID=BAM_REF,Number=1,Type=Integer,'
    'Description="Raw REF read count from BAM via samtools mpileup (no quality filters)">\n'
    '##INFO=<ID=BAM_ALT,Number=1,Type=Integer,'
    'Description="Raw ALT read count from BAM via samtools mpileup (no quality filters)">\n'
    '##INFO=<ID=BAM_VAF,Number=1,Type=Float,'
    'Description="Raw VAF = BAM_ALT / (BAM_REF + BAM_ALT) from unfiltered BAM counts">\n'
)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def open_vcf(path):
    return gzip.open(str(path), 'rt') if str(path).endswith('.gz') else open(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input',  required=True,  help='Input VCF (plain or .gz)')
    parser.add_argument('--output', required=True,  help='Output VCF.gz')
    parser.add_argument('--bam',    required=True,  help='BAM file')
    parser.add_argument('--ref',    default=str(WORKDIR / 'reference/genome.fa'),
                        help='Reference FASTA')
    args = parser.parse_args()

    out_path = Path(args.output)
    tmp_vcf  = str(out_path).replace('.vcf.gz', '.tmp.vcf')

    print(f'Input : {args.input}')
    print(f'BAM   : {args.bam}')
    print(f'Output: {args.output}')

    # --- Pass 1: read header + data lines ---
    header_lines = []
    data_lines   = []
    csq_fields   = []

    with open_vcf(args.input) as fh:
        for line in fh:
            if line.startswith('#'):
                header_lines.append(line)
                if not csq_fields:
                    csq_fields = parse_csq_fields([line])
            else:
                data_lines.append(line.rstrip('\n'))

    # --- Pass 2: collect variants for batch pileup ---
    variants = []
    for line in data_lines:
        c = line.split('\t')
        variants.append((c[0], int(c[1]), c[3], c[4]))

    print(f'Querying BAM for {len(variants):,} variants ...')
    counts = batch_bam_counts(args.bam, args.ref, variants)
    print('BAM query complete.')

    # --- Pass 3: write annotated VCF ---
    with open(tmp_vcf, 'w') as out:
        for line in header_lines:
            if line.startswith('#CHROM'):
                out.write(NEW_HEADERS)
            out.write(line)

        for i, line in enumerate(data_lines):
            cols  = line.split('\t')
            chrom = cols[0]
            pos   = int(cols[1])
            ref   = cols[3]
            alt   = cols[4]
            info  = cols[7] if len(cols) > 7 else '.'

            # 1. GENE
            gene = '.'
            for field in info.split(';'):
                if field.startswith('CSQ='):
                    gene = extract_gene(field[4:], csq_fields)
                    break

            # 2. BAM allele counts
            bam_ref, bam_alt = counts.get((chrom, pos, ref, alt), (0, 0))
            total   = bam_ref + bam_alt
            bam_vaf = round(bam_alt / total, 4) if total > 0 else 0.0

            # Append to INFO
            cols[7] = (info.rstrip(';')
                       + f';GENE={gene}'
                       + f';BAM_REF={bam_ref}'
                       + f';BAM_ALT={bam_alt}'
                       + f';BAM_VAF={bam_vaf}')

            out.write('\t'.join(cols) + '\n')

    # --- Compress + index ---
    subprocess.run([BGZIP, '-f', tmp_vcf], check=True)
    # bgzip creates tmp_vcf + '.gz'; rename to final output path
    tmp_gz = Path(tmp_vcf + '.gz')
    if tmp_gz != out_path:
        tmp_gz.rename(out_path)
    subprocess.run([BCFTOOLS, 'index', '-t', str(out_path)], check=True)
    print(f'\nWrote: {out_path}')

    # --- Sanity check on SRSF2 ---
    hit = subprocess.run(
        [BCFTOOLS, 'view', '-H', str(out_path), 'chr17:76736877-76736877'],
        capture_output=True, text=True
    ).stdout.strip()
    if hit:
        gene_m  = re.search(r'GENE=([^;]+)',    hit)
        bref_m  = re.search(r'BAM_REF=([^;]+)', hit)
        balt_m  = re.search(r'BAM_ALT=([^;]+)', hit)
        bvaf_m  = re.search(r'BAM_VAF=([^;]+)', hit)
        print(f'\nSRSF2 check → GENE={gene_m.group(1) if gene_m else "?"}'
              f'  BAM_REF={bref_m.group(1) if bref_m else "?"}'
              f'  BAM_ALT={balt_m.group(1) if balt_m else "?"}'
              f'  BAM_VAF={bvaf_m.group(1) if bvaf_m else "?"}')


if __name__ == '__main__':
    main()
