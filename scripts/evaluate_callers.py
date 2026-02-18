#!/usr/bin/env python3
"""
Evaluate and compare variant callers.
Calculates sensitivity, precision, F1 score, and VAF accuracy.
"""

import os
import sys
import json
import csv
from collections import defaultdict
from scipy.stats import pearsonr
import warnings
warnings.filterwarnings('ignore')

def parse_vcf(vcf_file):
    """
    Parse VCF file and return variants.
    Returns dict of (chrom, pos, ref, alt) -> {vaf, qual, ...}
    """
    variants = {}

    if not os.path.exists(vcf_file):
        print(f"Warning: {vcf_file} not found", file=sys.stderr)
        return variants

    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')

            if len(fields) < 8:
                continue

            try:
                chrom = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]
                qual = float(fields[5]) if fields[5] != '.' else 0

                # Parse sample columns if present
                vaf = None
                if len(fields) >= 10:
                    # Try to extract VAF from FORMAT/sample columns
                    format_fields = fields[8].split(':')
                    sample_fields = fields[9].split(':')

                    if 'VAF' in format_fields:
                        vaf_idx = format_fields.index('VAF')
                        try:
                            vaf = float(sample_fields[vaf_idx])
                        except (ValueError, IndexError):
                            pass

                # Extract from INFO if available
                if vaf is None and 'VAF=' in fields[7]:
                    info = fields[7]
                    for field in info.split(';'):
                        if field.startswith('VAF='):
                            try:
                                vaf = float(field.split('=')[1])
                            except ValueError:
                                pass

                key = (chrom, pos, ref, alt)
                variants[key] = {
                    'qual': qual,
                    'vaf': vaf,
                    'line': line.strip()
                }

            except (ValueError, IndexError) as e:
                print(f"Warning: Could not parse VCF line: {e}", file=sys.stderr)
                continue

    return variants

def fuzzy_match_indel(key1, key2, max_offset=5):
    """
    Check if two variants are the same allowing for indel representation differences.
    Returns True if they match within max_offset bases.
    """
    chrom1, pos1, ref1, alt1 = key1
    chrom2, pos2, ref2, alt2 = key2

    if chrom1 != chrom2:
        return False

    if abs(pos1 - pos2) > max_offset:
        return False

    # If both are SNVs and close, check
    if len(ref1) == len(ref2) == 1 and len(alt1) == len(alt2) == 1:
        return pos1 == pos2 and ref1 == ref2 and alt1 == alt2

    # For indels, allow fuzzy matching
    if pos1 == pos2:
        return True

    return False

def match_variants(truth_variants, caller_variants):
    """
    Match truth set variants with caller variants.
    Returns (tp, fp, fn) and mapping of truth to caller variants.
    """
    truth_keys = set(truth_variants.keys())
    caller_keys = set(caller_variants.keys())

    tp = 0
    fp = 0
    fn = 0
    matched_truth = set()
    matched_caller = set()
    truth_to_caller = {}

    # First pass: exact matches
    exact_matches = truth_keys & caller_keys
    for key in exact_matches:
        tp += 1
        matched_truth.add(key)
        matched_caller.add(key)
        truth_to_caller[key] = key

    # Second pass: fuzzy matches for unmatched indels
    unmatched_truth = truth_keys - matched_truth
    unmatched_caller = caller_keys - matched_caller

    for t_key in unmatched_truth:
        for c_key in unmatched_caller:
            if fuzzy_match_indel(t_key, c_key):
                tp += 1
                matched_truth.add(t_key)
                matched_caller.add(c_key)
                truth_to_caller[t_key] = c_key
                break

    # Count remaining
    fn = len(truth_keys - matched_truth)
    fp = len(caller_keys - matched_caller)

    return tp, fp, fn, truth_to_caller

def calculate_metrics(tp, fp, fn):
    """Calculate sensitivity, precision, F1 score."""
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    f1 = 2 * (precision * sensitivity) / (precision + sensitivity) if (precision + sensitivity) > 0 else 0

    return {
        'sensitivity': sensitivity,
        'precision': precision,
        'f1': f1,
        'tp': tp,
        'fp': fp,
        'fn': fn
    }

def calculate_vaf_accuracy(truth_variants, caller_variants, truth_to_caller):
    """
    Calculate VAF accuracy for variants with known VAF values.
    """
    # Collect VAF pairs for variants with known truth VAF
    truth_vafs = []
    caller_vafs = []

    for truth_key, caller_key in truth_to_caller.items():
        truth_vaf = truth_variants[truth_key]['vaf']
        caller_vaf = caller_variants[caller_key]['vaf']

        # Only include if both have VAF values
        if truth_vaf is not None and caller_vaf is not None:
            truth_vafs.append(truth_vaf)
            caller_vafs.append(caller_vaf)

    if len(truth_vafs) < 2:
        return {
            'correlation': None,
            'pvalue': None,
            'mae': None,
            'n_variants': len(truth_vafs)
        }

    # Calculate correlation
    try:
        r, pval = pearsonr(truth_vafs, caller_vafs)
    except:
        r, pval = None, None

    # Calculate MAE
    mae = sum(abs(t - c) for t, c in zip(truth_vafs, caller_vafs)) / len(truth_vafs)

    return {
        'correlation': r,
        'pvalue': pval,
        'mae': mae,
        'n_variants': len(truth_vafs),
        'truth_vafs': truth_vafs,
        'caller_vafs': caller_vafs
    }

def main():
    # Define caller VCF files
    WORKDIR = '/data/salomonis-archive/FASTQs/NCI-R01/rna_seq_varcall'
    truth_vcf = f'{WORKDIR}/truth_set/truth_set.vcf'
    caller_vcfs = {
        'HaplotypeCaller': f'{WORKDIR}/results/haplotypecaller/variants.vcf.gz',
        'DeepVariant':     f'{WORKDIR}/results/deepvariant/variants.vcf.gz',
        'Clair3_RNA':      f'{WORKDIR}/results/clair3_rna/merge_output.vcf.gz',
        'LongcallR':       f'{WORKDIR}/results/longcallr/variants.vcf.gz',
    }

    output_dir = f'{WORKDIR}/results/comparison'
    os.makedirs(output_dir, exist_ok=True)

    print(f"Loading truth set from {truth_vcf}...", file=sys.stderr)
    truth_variants = parse_vcf(truth_vcf)
    print(f"Loaded {len(truth_variants)} truth variants", file=sys.stderr)

    results = {}

    print("\nEvaluating callers...", file=sys.stderr)
    for caller_name, caller_vcf in caller_vcfs.items():
        print(f"  {caller_name}: {caller_vcf}", file=sys.stderr)

        caller_variants = parse_vcf(caller_vcf)
        print(f"    Loaded {len(caller_variants)} variants", file=sys.stderr)

        # Match variants
        tp, fp, fn, truth_to_caller = match_variants(truth_variants, caller_variants)

        # Calculate metrics
        metrics = calculate_metrics(tp, fp, fn)

        # Calculate VAF accuracy
        vaf_metrics = calculate_vaf_accuracy(truth_variants, caller_variants, truth_to_caller)

        results[caller_name] = {
            'metrics': metrics,
            'vaf_metrics': vaf_metrics,
            'truth_to_caller': {str(k): str(v) for k, v in truth_to_caller.items()}
        }

        print(f"    Sensitivity: {metrics['sensitivity']:.2%}", file=sys.stderr)
        print(f"    Precision: {metrics['precision']:.2%}", file=sys.stderr)
        print(f"    F1: {metrics['f1']:.3f}", file=sys.stderr)
        if vaf_metrics['correlation'] is not None:
            print(f"    VAF correlation: {vaf_metrics['correlation']:.3f}", file=sys.stderr)
            print(f"    VAF MAE: {vaf_metrics['mae']:.4f}", file=sys.stderr)

    # Write summary JSON
    summary_file = os.path.join(output_dir, 'summary.json')
    print(f"\nWriting summary to {summary_file}...", file=sys.stderr)

    # Convert Nones to strings for JSON serialization
    json_results = {}
    for caller_name, caller_results in results.items():
        json_results[caller_name] = {
            'metrics': caller_results['metrics'],
            'vaf_metrics': {k: v for k, v in caller_results['vaf_metrics'].items()
                           if k not in ['truth_vafs', 'caller_vafs']}
        }

    with open(summary_file, 'w') as f:
        json.dump(json_results, f, indent=2)

    # Write comparison table CSV
    table_file = os.path.join(output_dir, 'comparison.csv')
    print(f"Writing comparison table to {table_file}...", file=sys.stderr)

    with open(table_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Caller', 'TP', 'FP', 'FN', 'Sensitivity', 'Precision', 'F1', 'VAF_Corr', 'VAF_MAE'])

        for caller_name in sorted(results.keys()):
            res = results[caller_name]
            metrics = res['metrics']
            vaf_metrics = res['vaf_metrics']

            writer.writerow([
                caller_name,
                metrics['tp'],
                metrics['fp'],
                metrics['fn'],
                f"{metrics['sensitivity']:.3f}",
                f"{metrics['precision']:.3f}",
                f"{metrics['f1']:.3f}",
                f"{vaf_metrics['correlation']:.3f}" if vaf_metrics['correlation'] is not None else 'N/A',
                f"{vaf_metrics['mae']:.4f}" if vaf_metrics['mae'] is not None else 'N/A'
            ])

    print(f"\nEvaluation complete. Results in {output_dir}/", file=sys.stderr)

if __name__ == '__main__':
    main()
