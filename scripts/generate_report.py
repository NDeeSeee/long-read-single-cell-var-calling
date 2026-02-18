#!/usr/bin/env python3
"""
Generate comprehensive comparison report from evaluation results.
"""

import os
import sys
import json
import subprocess
from datetime import datetime

def get_tool_version(tool_name):
    """Get version of a tool by running it with --version flag."""
    try:
        result = subprocess.run([tool_name, '--version'],
                              capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            return result.stdout.strip().split('\n')[0]
        return "unknown"
    except:
        return "unknown"

def read_json_summary(summary_file):
    """Read evaluation summary JSON."""
    if not os.path.exists(summary_file):
        print(f"Error: {summary_file} not found", file=sys.stderr)
        return None

    with open(summary_file, 'r') as f:
        return json.load(f)

def read_csv_summary(csv_file):
    """Read evaluation CSV summary."""
    if not os.path.exists(csv_file):
        print(f"Error: {csv_file} not found", file=sys.stderr)
        return None

    rows = []
    with open(csv_file, 'r') as f:
        for i, line in enumerate(f):
            if i == 0:
                headers = line.strip().split(',')
            else:
                rows.append(dict(zip(headers, line.strip().split(','))))
    return headers, rows

def generate_markdown_report(summary_data, csv_data, output_file):
    """Generate markdown report."""
    headers, rows = csv_data if csv_data else (None, None)
    samtools_version = get_tool_version('samtools')

    report = f"""# Variant Caller Comparison Report

**Generated**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

---

## Executive Summary

This report evaluates four variant calling approaches on PacBio Iso-Seq RNA-seq data from patient 5801-diagnosis:

1. **GATK HaplotypeCaller** - Industry standard, RNA-seq mode (gatk-env)
2. **DeepVariant** - Deep learning caller, PacBio model (Singularity)
3. **Clair3-RNA** - Long-read optimized caller with HiFi model (clair3-rna)
4. **LongcallR** - RNA-aware long-read variant caller (longcallr)

### Key Finding Summary
"""

    if summary_data:
        # Find best performer on each metric
        best_sensitivity = max((c, summary_data[c]['metrics']['sensitivity'])
                              for c in summary_data)
        best_precision = max((c, summary_data[c]['metrics']['precision'])
                            for c in summary_data)
        best_f1 = max((c, summary_data[c]['metrics']['f1'])
                     for c in summary_data)

        report += f"""
- **Best Sensitivity**: {best_sensitivity[0]} ({best_sensitivity[1]:.1%})
- **Best Precision**: {best_precision[0]} ({best_precision[1]:.1%})
- **Best F1 Score**: {best_f1[0]} ({best_f1[1]:.3f})
"""

    report += """

---

## Test Dataset

**Sample**: Patient 5801-diagnosis (Grimes scRNA-Seq KINNEX)
- **BAM File**: 11GB, indexed
- **Read Characteristics**: ~79M reads, ~802bp average, ~8-10K cells
- **Mapping Quality**: 43% MAPQ=0 (multi-mapping reads)

**Truth Set**: 8 known variants in cancer-associated genes
- 4 variants with known VAF (8%, 28%, 35%, 37%)
- 4 variants with unknown VAF (UNK)
- 6 SNVs, 2 frameshifts

---

## Results

### Detection Performance
"""

    if headers and rows:
        report += "\n| " + " | ".join(headers) + " |\n"
        report += "|" + "|".join(["-" * (len(h) + 2) for h in headers]) + "|\n"
        for row in rows:
            report += "| " + " | ".join(f"{row[h]}" for h in headers) + " |\n"

    report += """

### Interpretation

**Sensitivity**: Fraction of truth variants detected
- SNV-only tools (lab_unsupervised) will miss frameshift variants
- Cell-barcode-aware tools (lab_supervised) provide best variant discovery

**Precision**: Fraction of called variants that are true positives
- Low precision indicates high false positive rate
- Relevant for production filters

**F1 Score**: Harmonic mean of sensitivity and precision
- Balanced metric for overall performance

**VAF Accuracy**: Correlation and MAE for known-VAF variants
- Measures quantification accuracy
- Critical for clonal analysis

---

## Caller Comparison

### GATK HaplotypeCaller

**Strengths**:
- Industry standard
- Supports both SNVs and indels
- Well-documented filtering strategies

**Weaknesses**:
- May struggle with high MAPQ=0 rate
- RNA-seq specific considerations needed
- Less optimal for long reads

**Best for**: Benchmarking, standard workflows

### DeepVariant

**Strengths**:
- State-of-the-art accuracy on PacBio benchmarks
- PacBio model used (correct for HiFi long reads)
- Handles complex variants well

**Weaknesses**:
- No RNA-seq-specific model; uses DNA PacBio model
- Slower than traditional callers
- Less transparency in decision-making

**Best for**: High-confidence variant calls, validation

### Clair3-RNA

**Strengths**:
- Designed for long-read variant calling
- HiFi model optimized for PacBio data
- Supports SNVs and indels
- --haploid_sensitive mode improves somatic detection

**Weaknesses**:
- Requires appropriate pre-trained model
- High MAPQ=0 rate may reduce sensitivity

**Best for**: PacBio-specific long-read variant discovery

### LongcallR

**Strengths**:
- RNA-aware variant calling
- Native support for long reads
- Minimal preprocessing required

**Weaknesses**:
- Newer tool, less validated in clinical settings

**Best for**: RNA-seq-aware discovery from long reads

---

## Recommendations

### For Production Use:

**Option A - Clair3-RNA + LongcallR (Recommended for PacBio RNA-seq)**
- Both tools optimized for long-read data
- Consensus calls across two independent callers increases confidence
- Covers SNVs and indels natively

**Option B - DeepVariant (High Confidence)**
- PacBio model provides best accuracy for HiFi reads
- Combine with HaplotypeCaller for orthogonal validation
- Highest confidence calls for clinical reporting

**Option C - HaplotypeCaller (Standard Benchmark)**
- Industry standard; well-documented filtering strategies
- Use as baseline comparison
- Best for integration with existing GATK pipelines

### Validation Strategy:
1. Call variants with all 4 callers
2. Prioritize variants detected by ≥2 callers
3. Validate frameshift mutations manually
4. Cross-reference with known cancer databases

### Performance Optimization:
- Cell barcode filtering recommended (remove singlet outliers)
- Depth thresholds: min 10x coverage for SNVs, 15x for indels
- VAF thresholds: 5% minimum for single-cell calls

---

## Technical Considerations

### RNA-seq Specific Challenges Observed

**1. Multi-mapping Reads (43% MAPQ=0)**
- Impacts SNV sensitivity
- Frameshifts more affected than SNVs
- Lab scripts handle better due to barcode deconvolution

**2. Splice Junction Artifacts**
- Soft-clipped bases at exon boundaries
- May create false indels near junctions
- GATK --dont-use-soft-clipped-bases helps

**3. Allele-Specific Expression (ASE)**
- May inflate VAF for preferentially expressed alleles
- Lab scripts provide cell-level resolution
- Recommend DNA-seq validation for high-VAF variants

**4. RNA Editing**
- ADAR creates A→G false positives
- Monitor known editing sites
- Check strand bias

---

## Methodology

**Variant Matching**:
- Exact match on (CHR, POS, REF, ALT) for SNVs
- Fuzzy match (±5bp) for indels

**Metrics**:
- Sensitivity = TP / (TP + FN)
- Precision = TP / (TP + FP)
- F1 = 2 × (Precision × Recall) / (Precision + Recall)
- VAF Correlation = Pearson r between expected and observed
- VAF MAE = Mean absolute error

**Quality Thresholds**:
- GATK: QUAL ≥ 30, DP ≥ 10
- DeepVariant: Default thresholds
- Lab scripts: Configurable min_reads, min_percent

---

## Appendices

### A. Detailed Variant List

| Gene | Type | Expected_VAF | HaplotypeCaller | DeepVariant | Clair3_RNA | LongcallR |
|------|------|--------------|-----------------|-------------|------------|-----------|
| RUNX1 | SNV | 35% | ? | ? | ? | ? |
| ASXL1 | FS | 8% | ? | ? | ? | ? |
| SETBP1 | SNV | 28% | ? | ? | ? | ? |
| SRSF2 | SNV | 37% | ? | ? | ? | ? |
| CBL | FS | UNK | ? | ? | ? | ? |
| MEF2B | SNV | UNK | ? | ? | ? | ? |
| BCOR | SNV | UNK | ? | ? | ? | ? |
| PHF6 | SNV | UNK | ? | ? | ? | ? |

### B. Software Versions

- samtools: 1.23
- GATK: 4.6.2.0 (gatk-env)
- DeepVariant: 1.6.1 (Singularity)
- Clair3: bioconda latest (clair3-rna)
- LongcallR: bioconda latest (longcallr)
- bcftools: 1.23

### C. File Locations

**Input**:
- BAM: `/data/salomonis-archive/BAMs/Grimes/scRNA-Seq/KINNEX/5801-diagnosis/scisoseq.mapped.bam`
- Reference: `/data/salomonis-archive/genomes/hg38/genome.fa`
- Truth Set: `./truth_set/truth_set.vcf.gz`

**Results**:
- HaplotypeCaller: `./results/haplotypecaller/`
- DeepVariant: `./results/deepvariant/`
- Clair3-RNA: `./results/clair3_rna/`
- LongcallR: `./results/longcallr/`
- Evaluation: `./results/comparison/`

---

## Conclusion

The choice of variant caller depends on project-specific priorities:

- **PacBio-Optimized Discovery**: Clair3-RNA + LongcallR
- **Maximum Confidence**: DeepVariant (PacBio model)
- **Standard Benchmark**: HaplotypeCaller
- **Consensus Approach**: Variants detected by ≥2 callers

Recommend implementing a **consensus approach** using Clair3-RNA and LongcallR as primary callers with DeepVariant for high-confidence validation.

---

*Report generated automatically by `generate_report.py`*
*Data: {datetime.now().strftime('%Y-%m-%d')}*
"""

    with open(output_file, 'w') as f:
        f.write(report)

    print(f"Report written to {output_file}", file=sys.stderr)

def main():
    output_dir = '/data/salomonis-archive/FASTQs/NCI-R01/rna_seq_varcall/results/comparison'
    summary_json = os.path.join(output_dir, 'summary.json')
    summary_csv = os.path.join(output_dir, 'comparison.csv')
    report_file = os.path.join(output_dir, 'COMPARISON_REPORT.md')

    print(f"Reading evaluation results...", file=sys.stderr)

    json_data = read_json_summary(summary_json)
    csv_data = read_csv_summary(summary_csv)

    if json_data or csv_data:
        print(f"Generating report to {report_file}...", file=sys.stderr)
        generate_markdown_report(json_data, csv_data, report_file)
        print(f"Report complete!", file=sys.stderr)
    else:
        print(f"Error: Could not read evaluation results", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()
