# RNA-seq Somatic Variant Calling Pipeline

Detection of somatic mutations in hematologic malignancies using PacBio Iso-Seq (KINNEX) single-cell RNA-seq data. The pipeline runs four independent variant callers, annotates against population databases via VEP, and applies a multi-caller concordance filter to isolate high-confidence somatic candidates.

---

## Sample

| Field | Value |
|---|---|
| Sample | 5801-diagnosis |
| Data type | PacBio Iso-Seq / KINNEX scRNA-seq |
| BAM | `/data/salomonis-archive/BAMs/Grimes/scRNA-Seq/KINNEX/5801-diagnosis/scisoseq.mapped.bam` |
| Quality-fixed BAM | `scisoseq_synthetic_qual.bam` (synthetic Q30, required by GATK) |
| Reference | `reference/genome.fa` → GRCh38 from 10x Genomics Space Ranger |

---

## Truth Set (8 Known Somatic Variants)

| Gene | Position (GRCh38) | Mutation | VAF |
|---|---|---|---|
| SRSF2 | chr17:76736877 | p.Pro95Arg (missense) | 37% |
| SETBP1 | chr18:44951948 | missense | 28% |
| RUNX1 | chr21:34799432 | nonsense | 35% |
| ASXL1 | chr20:32434638 | frameshift | 8% |
| CBL | chr11:119242488 | frameshift | UNK |
| MEF2B | chr19:19145882 | missense | UNK |
| BCOR | chrX:40057210 | nonsense | UNK |
| PHF6 | chrX:134415107 | missense | UNK |

> **Note — CBL**: This variant is in a 177 kb intron. Iso-Seq reads skip intronic sequence (CIGAR `N` operations), making it fundamentally undetectable by RNA-seq. This is a data limitation, not a caller bug.

---

## Pipeline Overview

```
BAM (scisoseq_synthetic_qual.bam)
        │
        ├─► GATK HaplotypeCaller  ─┐
        ├─► DeepVariant (PACBIO)  ─┤
        ├─► Clair3-RNA (hifi)     ─┤── VEP annotation ──► classify_variants.py
        └─► LongcallR (hifi-isoseq)┘         │
                                              ▼
                                    somatic_known / novel_candidates /
                                    rare_candidates / uncertain / germline
                                              │
                                    merge_callers.py (≥2 caller concordance)
                                              │
                                    mosdepth BAM depth filter (BAM_DP > 1500)
                                              │
                                    high_confidence_bamdp1500.vcf.gz  ◄── final output
```

---

## Quick Start

```bash
cd /data/salomonis-archive/FASTQs/NCI-R01/rna_seq_varcall

# Step 1: Run all 4 callers (full genome)
bash scripts/run_callers.sh all

# Step 2: Annotate each caller's VCF with VEP + gnomAD + COSMIC
bash scripts/annotate_and_filter.sh results/haplotypecaller/variants_raw.vcf     gatk
bash scripts/annotate_and_filter.sh results/deepvariant/variants.vcf.gz          deepvariant
bash scripts/annotate_and_filter.sh results/clair3_rna/variants.vcf.gz           clair3
bash scripts/annotate_and_filter.sh results/longcallr/variants.vcf.gz            longcallr

# Step 3: Classify variants (runs automatically inside Step 2, but can be re-run)
python scripts/classify_variants.py \
    --input annotation/gatk/annotated.vcf.gz \
    --output annotation/gatk --sample gatk

# Step 4: Multi-caller concordance merge
python scripts/merge_callers.py

# Step 5: Apply raw BAM depth filter (BAM_DP > 1500) using mosdepth
#         (see Filtering section below for the full script)
```

---

## Variant Callers

| Caller | Container / Env | Model | SNPs | Indels | Full-genome output |
|---|---|---|---|---|---|
| GATK HaplotypeCaller | `gatk-env` conda | RNA-seq | ✓ | ✓ | `results/haplotypecaller/variants_raw.vcf` |
| DeepVariant | `containers/deepvariant_1.6.1.sif` | PACBIO | ✓ | ✓ | `results/deepvariant/variants.vcf.gz` |
| Clair3-RNA | `containers/clair3_latest.sif` | hifi (per-chrom) | ✓ | ✓ | `results/clair3_rna/variants.vcf.gz` |
| LongcallR | `containers/longcallr_latest.sif` | hifi-isoseq | ✓ | ✗ | `results/longcallr/variants.vcf.gz` |

### Full-genome variant counts

| Caller | Total | SNPs | Indels | Ts/Tv |
|---|---|---|---|---|
| GATK | 1,619,585 | 1,100,978 | 482,378 | 3.05 |
| DeepVariant | 2,449,434 | ~2.2M | ~250K | 3.53 |
| Clair3-RNA | 3,611,894 | 3,334,002 | 279,377 | — |
| LongcallR | 441,710 | 441,710 | 0 | 3.92 |

> **LongcallR** is SNP-only by design. It cannot call indels.

---

## VEP Annotation

VEP (Ensembl v112, GRCh38) is run via Singularity:

```
containers/ensembl-vep_release_112.0.sif
vep_cache/homo_sapiens/112_GRCh38/   (1.1 GB)
```

**Databases annotated per variant:**
- gnomAD exomes (r2.1.1) and genomes (v3.1.2) — 14 population AFs
- dbSNP 156 — rsIDs
- COSMIC 98 — somatic mutation IDs (COSV/COSM)
- ClinVar (Oct 2023) — clinical significance
- SIFT / PolyPhen — protein function predictions
- HGVS notation (HGVSc, HGVSp)
- Consequence, gene symbol, canonical transcript, MANE Select

Output per caller: `annotation/{caller}/annotated.vcf.gz`

---

## Variant Classification

`scripts/classify_variants.py` reads each annotated VCF and assigns one of five tiers:

| Class | Criteria | Action |
|---|---|---|
| **GERMLINE** | gnomAD AF ≥ 0.001 **or** ALT reads < 5 | Excluded from somatic analysis |
| **SOMATIC_KNOWN** | COSMIC match (any gnomAD AF) **or** ClinVar pathogenic + gnomAD AF < 0.001 | High-priority candidate |
| **NOVEL_CANDIDATE** | gnomAD AF = 0, QUAL ≥ 30, ALT reads ≥ 5 | Candidate — absent from population |
| **RARE_CANDIDATE** | gnomAD AF < 0.0001, ALT reads ≥ 5 | Candidate — very rare in population |
| **UNCERTAIN** | gnomAD AF 0.0001–0.001, ALT reads ≥ 5 | Grey zone — rare but present in population |

**Decision priority**: SOMATIC_KNOWN > gnomAD filter > NOVEL > RARE > UNCERTAIN

**Key thresholds:**
- `GERMLINE_AF = 0.001` — variants above this frequency in any gnomAD population are excluded
- `MIN_ALT_DEPTH = 5` — minimum ALT-supporting reads (uses FORMAT/AD; falls back to DP×AF for LongcallR)
- `MIN_QUAL = 30` — minimum QUAL score for NOVEL_CANDIDATE

Classification outputs per caller:
```
annotation/{caller}/
├── annotated.vcf.gz          # VEP-annotated input
├── somatic_known.vcf.gz
├── novel_candidates.vcf.gz
├── rare_candidates.vcf.gz
├── uncertain.vcf.gz
└── germline.vcf.gz           # excluded
```

### Classification summary (post-filter)

| Class | GATK | DeepVariant | Clair3 | LongcallR |
|---|---|---|---|---|
| GERMLINE | 85.4% | 89.9% | 91.6% | 73.1% |
| NOVEL_CANDIDATE | 7.6% | 0.1% | 0.2% | 16.8% |
| SOMATIC_KNOWN | 3.1% | 1.9% | 1.6% | 3.9% |
| UNCERTAIN | 2.4% | 0.9% | 0.6% | 1.5% |
| RARE_CANDIDATE | 1.5% | 7.2% | 6.0% | 4.7% |

---

## Multi-Caller Concordance Merge

`scripts/merge_callers.py` merges all four callers' somatic buckets (excluding germline) and counts how many independent callers support each variant site:

```bash
python scripts/merge_callers.py \
    [--callers gatk deepvariant clair3 longcallr] \
    [--annotation_dir annotation/] \
    [--output results/merged/]
```

**Output tiers:**

| File | Variants | Criterion |
|---|---|---|
| `high_confidence.vcf.gz` | 66,280 | ≥ 3 callers agree |
| `medium_confidence.vcf.gz` | 118,114 | 2 callers agree |
| `low_confidence.vcf.gz` | 332,620 | 1 caller only (likely noise) |
| `all_candidates.vcf.gz` | 517,014 | All, tagged with `CALLER_COUNT` and `CALLERS` |

All output VCFs carry the full VEP `CSQ` annotation, `SOMATIC_CLASS`, `SOMATIC_REASON`, `CALLER_COUNT`, and `CALLERS` INFO fields.

---

## Depth Filtering — Critical Note for PacBio RNA-seq

### Problem: GATK FORMAT/DP severely underestimates true depth

GATK excludes reads with MQ=255 (PacBio's "mapping quality not computed" flag for splice-junction-spanning reads) and soft-clipped bases (`--dont-use-soft-clipped-bases`). This means GATK's `FORMAT/DP` typically reflects only **2–5% of actual BAM reads**.

**Example — SRSF2 (chr17:76736877):**

| Source | Depth |
|---|---|
| Raw BAM (mosdepth) | 5,548 |
| GATK FORMAT/DP | 162 (2.9%) |
| GATK AD (REF,ALT) | 88, 74 |

Filtering on `FORMAT/DP > 1500` would incorrectly discard SRSF2.

### Solution: Use mosdepth for true BAM depth

`mosdepth` counts all reads in the BAM without GATK's internal filters. The pipeline annotates each variant with `INFO/BAM_DP` from mosdepth and filters on that:

```bash
MOSDEPTH=/usr/local/anaconda3-2020/envs/mosdepth/bin/mosdepth

# Run mosdepth over variant positions
$MOSDEPTH --by positions.bed --no-per-base --threads 8 \
    tmp/mosdepth_out scisoseq_synthetic_qual.bam

# Annotate VCF with BAM_DP and filter
bcftools annotate -a bam_depth.tsv.gz -c CHROM,POS,INFO/BAM_DP ...
bcftools filter -i 'INFO/BAM_DP > 1500' ...
```

### Final filtered output files (BAM_DP > 1500)

| File | Variants | SNPs | Indels | Description |
|---|---|---|---|---|
| `annotation/gatk/non_germline_bamdp1500.vcf.gz` | 1,351 | 1,056 | 295 | All GATK non-germline, raw depth > 1500 |
| `results/merged/high_confidence_bamdp1500.vcf.gz` | **983** | 968 | 15 | ≥3 callers, raw depth > 1500 — **primary output** |
| `results/merged/medium_confidence_bamdp1500.vcf.gz` | 289 | 246 | 43 | 2 callers, raw depth > 1500 |
| `results/merged/low_confidence_bamdp1500.vcf.gz` | 609 | 224 | 387 | 1 caller, raw depth > 1500 |

---

## Key Scripts

| Script | Purpose |
|---|---|
| `scripts/run_callers.sh` | Run all 4 variant callers (full genome) |
| `scripts/annotate_and_filter.sh` | VEP annotation → somatic classification |
| `scripts/classify_variants.py` | Tiered somatic classification with gnomAD + COSMIC |
| `scripts/merge_callers.py` | Multi-caller concordance merge → confidence tiers |

---

## Key Biological Caveats

1. **No matched normal**: Without paired germline DNA, rare germline variants are indistinguishable from somatic mutations by population frequency alone. COSMIC and ClinVar annotations are the primary evidence for somatic origin.

2. **RNA editing (A→I)**: ADAR-mediated RNA editing produces A>G variants that appear as real mutations. A future filter using REDIportal (~4.5M known editing sites) would reduce this noise.

3. **Intronic variants undetectable**: Iso-Seq reads skip introns entirely (CIGAR `N`). Any mutation in an intron (e.g., CBL chr11:119242488 in a 177 kb intron) cannot be detected from RNA-seq data regardless of caller or depth.

4. **UNCERTAIN bucket**: The 116 GATK variants (9 in high-confidence) with gnomAD AF 0.01–0.1% are included but represent a grey zone between rare germline and somatic.

---

## Software Versions

| Tool | Version | Location |
|---|---|---|
| samtools | 1.23 | `bio-cli` conda env |
| bcftools | 1.23 | `bio-cli` conda env |
| GATK | 4.6.2.0 | `gatk-env` conda env |
| DeepVariant | 1.6.1 | `containers/deepvariant_1.6.1.sif` |
| Clair3 | latest | `containers/clair3_latest.sif` |
| LongcallR | latest | `containers/longcallr_latest.sif` |
| VEP | 112.0 | `containers/ensembl-vep_release_112.0.sif` |
| mosdepth | 0.3.3 | `mosdepth` conda env |
| Singularity | — | `bio-cli` conda env |

---

**Last updated**: 2026-02-26
