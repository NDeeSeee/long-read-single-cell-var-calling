#!/usr/bin/env python3
"""
merge_callers.py

Merges somatic candidate VCFs from multiple callers (post-VEP annotation
outputs), counts how many callers agree on each variant site, and writes
confidence-tiered output VCFs.

Usage:
    python scripts/merge_callers.py \
        [--callers gatk deepvariant clair3 longcallr] \
        [--annotation_dir annotation/] \
        [--output results/merged/]

Output tiers:
    high_confidence.vcf.gz    ≥3 callers agree  — highest reliability
    medium_confidence.vcf.gz   2 callers agree  — solid candidates
    low_confidence.vcf.gz      1 caller only    — likely noise / filtered
    all_candidates.vcf.gz      full set tagged with CALLER_COUNT and CALLERS
"""

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

WORKDIR  = Path('/data/salomonis-archive/FASTQs/NCI-R01/rna_seq_varcall')
BCFTOOLS = '/users/pavb5f/.conda/envs/bio-cli/bin/bcftools'
BGZIP    = '/users/pavb5f/.conda/envs/bio-cli/bin/bgzip'

# Somatic buckets to pull from each caller's annotation output.
# GERMLINE is intentionally excluded.
SOMATIC_BUCKETS = ['somatic_known', 'novel_candidates', 'rare_candidates', 'uncertain']


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run(cmd, check=True):
    result = subprocess.run(
        [str(c) for c in cmd], check=False,
        capture_output=True, text=True
    )
    if check and result.returncode != 0:
        print(f"\nERROR running:\n  {' '.join(str(c) for c in cmd)}", file=sys.stderr)
        print(result.stderr[-2000:], file=sys.stderr)
        sys.exit(result.returncode)
    return result


def is_alt_call(gt):
    """True if the genotype contains at least one ALT allele (not ref, not missing)."""
    if not gt or gt in ('.', './.', '.|.'):
        return False
    alleles = gt.replace('|', '/').split('/')
    return any(a not in ('0', '.') for a in alleles)


# ---------------------------------------------------------------------------
# Per-caller preparation
# ---------------------------------------------------------------------------

def combine_caller_vcfs(caller, ann_dir, tmp_dir):
    """
    Concatenate and sort all somatic bucket VCFs for one caller.
    Returns path to sorted VCF.gz, or None if no data found.
    """
    caller_dir = ann_dir / caller
    inputs = []
    for bucket in SOMATIC_BUCKETS:
        vcf = caller_dir / f'{bucket}.vcf.gz'
        # Skip empty/missing files (header-only files are ~few hundred bytes)
        if vcf.exists() and vcf.stat().st_size > 2048:
            inputs.append(str(vcf))

    if not inputs:
        print(f"  [SKIP] {caller}: no somatic VCFs found in {caller_dir}/",
              file=sys.stderr)
        return None

    combined = tmp_dir / f'{caller}_combined.vcf.gz'
    sorted_  = tmp_dir / f'{caller}_sorted.vcf.gz'

    run([BCFTOOLS, 'concat', '-a', '-D', '-O', 'z', '-o', combined] + inputs)
    # Use workdir as sort temp space (avoids small /tmp)
    sort_tmp = tmp_dir / f'{caller}_sorttmp'
    sort_tmp.mkdir(exist_ok=True)
    run([BCFTOOLS, 'sort', '-T', sort_tmp, '-O', 'z', '-o', sorted_, combined])
    run([BCFTOOLS, 'index', '-t', sorted_])

    n = int(run([BCFTOOLS, 'view', '-H', sorted_]).stdout.count('\n'))
    print(f"  {caller:<15} {n:>9,} candidates  ({len(inputs)} bucket(s))")
    return sorted_


# ---------------------------------------------------------------------------
# Merge and score
# ---------------------------------------------------------------------------

def score_concordance(merged_vcf, active_callers, out_dir):
    """
    Read the multi-sample merged VCF, count ALT calls per site,
    and write tiered output VCFs.
    """
    fouts = {
        'high':   open(out_dir / 'high_confidence.vcf',   'w'),
        'medium': open(out_dir / 'medium_confidence.vcf',  'w'),
        'low':    open(out_dir / 'low_confidence.vcf',     'w'),
        'all':    open(out_dir / 'all_candidates.vcf',     'w'),
    }
    counts = {k: 0 for k in fouts}
    n_callers = len(active_callers)

    INFO_HEADERS = (
        '##INFO=<ID=CALLER_COUNT,Number=1,Type=Integer,'
        'Description="Number of callers with an ALT call at this site">\n'
        '##INFO=<ID=CALLERS,Number=1,Type=String,'
        'Description="Comma-separated list of supporting callers">\n'
    )

    with open(merged_vcf) as fh:
        for line in fh:
            if line.startswith('##'):
                for f in fouts.values():
                    f.write(line)
                continue

            if line.startswith('#CHROM'):
                for f in fouts.values():
                    f.write(INFO_HEADERS)
                    f.write(line)
                continue

            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue

            fmt    = cols[8].split(':')
            gt_idx = fmt.index('GT') if 'GT' in fmt else 0

            supporting = [
                active_callers[i]
                for i in range(n_callers)
                if (9 + i) < len(cols) and
                   is_alt_call(cols[9 + i].split(':')[gt_idx])
            ]
            n = len(supporting)

            info = cols[7] if cols[7] != '.' else ''
            cols[7] = (info + f';CALLER_COUNT={n};CALLERS={",".join(supporting)}'
                       ).lstrip(';')
            out_line = '\t'.join(cols) + '\n'

            fouts['all'].write(out_line)
            counts['all'] += 1

            if n >= 3:
                fouts['high'].write(out_line)
                counts['high'] += 1
            elif n == 2:
                fouts['medium'].write(out_line)
                counts['medium'] += 1
            else:
                fouts['low'].write(out_line)
                counts['low'] += 1

    for f in fouts.values():
        f.close()
    return counts


def compress_outputs(out_dir):
    for name in ['high_confidence', 'medium_confidence', 'low_confidence', 'all_candidates']:
        vcf = out_dir / f'{name}.vcf'
        gz  = out_dir / f'{name}.vcf.gz'
        if not vcf.exists():
            continue
        run([BGZIP, '-f', vcf])
        if gz.exists() and gz.stat().st_size > 2048:
            run([BCFTOOLS, 'index', '-t', gz])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--callers', nargs='+',
                        default=['gatk', 'deepvariant', 'clair3', 'longcallr'],
                        help='Caller names (must match annotation/ subdirectory names)')
    parser.add_argument('--annotation_dir', default='annotation/',
                        help='Base annotation directory  (default: annotation/)')
    parser.add_argument('--output', default='results/merged/',
                        help='Output directory  (default: results/merged/)')
    args = parser.parse_args()

    ann_dir = Path(args.annotation_dir)
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=== Multi-Caller Concordance Merge ===")
    print(f"Callers        : {', '.join(args.callers)}")
    print(f"Annotation dir : {ann_dir}/")
    print(f"Output dir     : {out_dir}/\n")

    with tempfile.TemporaryDirectory(dir=WORKDIR / 'tmp') as tmp:
        tmp_dir = Path(tmp)
        (WORKDIR / 'tmp').mkdir(exist_ok=True)

        # --- Step 1: combine buckets per caller ---
        print("[1] Combining somatic buckets per caller ...")
        combined_vcfs  = []
        active_callers = []
        for caller in args.callers:
            vcf = combine_caller_vcfs(caller, ann_dir, tmp_dir)
            if vcf:
                combined_vcfs.append(str(vcf))
                active_callers.append(caller)

        if len(combined_vcfs) < 2:
            sys.exit("ERROR: at least 2 callers with annotated results are required.")

        # --- Step 2: multi-sample merge ---
        print(f"\n[2] Merging {len(combined_vcfs)} caller VCFs ...")
        merged_gz  = tmp_dir / 'merged.vcf.gz'
        merged_vcf = tmp_dir / 'merged.vcf'

        run([BCFTOOLS, 'merge',
             '--merge', 'none',     # keep multiallelic records separate
             '--force-samples',     # allow duplicate sample names
             '-O', 'z', '-o', merged_gz] + combined_vcfs)
        run([BCFTOOLS, 'view', '-o', merged_vcf, merged_gz])

        n_sites = run([BCFTOOLS, 'view', '-H', str(merged_gz)]).stdout.count('\n')
        print(f"  Merged site count: {n_sites:,}")

        # --- Step 3: score concordance and write tiers ---
        print(f"\n[3] Scoring concordance (callers: {active_callers}) ...")
        counts = score_concordance(merged_vcf, active_callers, out_dir)

        # --- Step 4: compress and index ---
        print("\n[4] Compressing outputs ...")
        compress_outputs(out_dir)

    # --- Summary ---
    print("\n=== Concordance Summary ===")
    print(f"  Active callers        : {', '.join(active_callers)}")
    print(f"  High   (≥3 callers)   : {counts['high']:>9,}")
    print(f"  Medium  (2 callers)   : {counts['medium']:>9,}")
    print(f"  Low     (1 caller)    : {counts['low']:>9,}")
    print(f"  All candidates        : {counts['all']:>9,}")
    print(f"\nOutputs in {out_dir}/")
    print("  high_confidence.vcf.gz    ≥3 callers — highest reliability")
    print("  medium_confidence.vcf.gz   2 callers — good candidates")
    print("  low_confidence.vcf.gz      1 caller  — likely noise")
    print("  all_candidates.vcf.gz      full set with CALLER_COUNT tag")


if __name__ == '__main__':
    main()
