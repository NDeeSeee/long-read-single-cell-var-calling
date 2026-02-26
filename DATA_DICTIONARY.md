# Data Dictionary — RNA-seq Somatic Variant Calling Pipeline

**Sample**: 5801-diagnosis (PacBio Iso-Seq / KINNEX scRNA-seq, Grimes dataset)
**Reference**: GRCh38 (chr-prefixed, 10x Genomics Space Ranger build)
**Last updated**: 2026-02-26

---

## Output File Catalog

### Primary Outputs (use these)

| File | Variants | Format | Description |
|---|---|---|---|
| `results/merged/high_confidence_final.vcf.gz` | 983 | VCF.gz | **Primary output.** ≥3 callers concordant + BAM_DP>1500. Includes GENE, BAM allele counts, VAF. |
| `results/merged/high_confidence_final.tsv` | 983 | TSV | Human-readable version of above. One row per variant, key fields only. |
| `annotation/gatk/non_germline_final.vcf.gz` | 1,351 | VCF.gz | GATK-only non-germline variants with BAM_DP>1500, GENE, and BAM counts. |
| `annotation/gatk/non_germline_final.tsv` | 1,351 | TSV | Human-readable version of above. |
| `annotation/gatk/germline_bamdp1500.vcf.gz` | 1,552 | VCF.gz | GATK germline variants (gnomAD AF≥0.001) with BAM_DP>1500. High-expression common variants. |

### Secondary Outputs

| File | Variants | Description |
|---|---|---|
| `results/merged/medium_confidence_final.vcf.gz` | 289 | 2-caller concordance, BAM_DP>1500, with GENE + BAM counts |
| `results/merged/medium_confidence_final.tsv` | 289 | Human-readable version |
| `results/merged/low_confidence_final.vcf.gz` | 609 | Single-caller only, BAM_DP>1500, with GENE + BAM counts |
| `annotation/gatk/non_germline_bamdp1500.vcf.gz` | 1,351 | Pre-annotation version of non_germline_final (no GENE/BAM counts) |
| `results/merged/high_confidence_bamdp1500.vcf.gz` | 983 | Pre-annotation version of high_confidence_final |

### Per-Caller Annotated Files

Each caller produces a full VEP-annotated file and five somatic classification buckets:

```
annotation/{gatk,deepvariant,clair3,longcallr}/
├── annotated.vcf.gz          — VEP-annotated, all variants
├── somatic_known.vcf.gz      — COSMIC / ClinVar pathogenic
├── novel_candidates.vcf.gz   — gnomAD AF = 0, QUAL ≥ 30, AD ≥ 5
├── rare_candidates.vcf.gz    — gnomAD AF < 0.0001
├── uncertain.vcf.gz          — gnomAD AF 0.0001–0.001
└── germline.vcf.gz           — gnomAD AF ≥ 0.001 or AD < 5 (excluded)
```

---

## TSV Column Definitions

These columns appear in all `*_final.tsv` files:

| Column | Type | Description |
|---|---|---|
| `CHROM` | string | Chromosome (chr-prefixed GRCh38) |
| `POS` | int | 1-based genomic position |
| `REF` | string | Reference allele |
| `ALT` | string | Alternate allele (first ALT for multi-allelic sites) |
| `GENE` | string | Hugo gene symbol from VEP CSQ. Among all VEP transcripts, picks the canonical transcript with highest functional impact (HIGH > MODERATE > LOW > MODIFIER). |
| `CONSEQUENCE` | string | Sequence ontology consequence term from highest-impact VEP transcript (e.g. `missense_variant`, `stop_gained`, `frameshift_variant`) |
| `HGVSp` | string | Protein-level HGVS notation from canonical transcript (e.g. `p.Pro95Arg`) |
| `HGVSc` | string | cDNA-level HGVS notation from canonical transcript |
| `COSMIC` | string | COSMIC mutation IDs (COSV/COSM) from VEP Existing_variation field; `&`-separated if multiple |
| `CLINVAR` | string | ClinVar clinical significance from VEP annotation |
| `gnomAD_AF_max` | float | Maximum gnomAD allele frequency across all populations and transcripts annotated by VEP. Zero if absent from gnomAD. |
| `SOMATIC_CLASS` | string | Somatic classification tier (see table below) |
| `SOMATIC_REASON` | string | Evidence string supporting SOMATIC_CLASS (e.g. `COSMIC`, `ClinVar:likely_pathogenic`, `not_in_gnomAD`, `gnomAD_AF=0.0003`) |
| `BAM_REF` | int | REF allele read count directly from BAM via `samtools mpileup -B -Q 0 -q 0 -A`. Includes MQ=255 PacBio reads excluded by GATK. |
| `BAM_ALT` | int | ALT allele read count directly from BAM (same flags as above) |
| `BAM_VAF` | float | BAM_ALT / (BAM_REF + BAM_ALT). Raw VAF from unfiltered BAM counts. |
| `BAM_DP` | int | Total BAM depth from mosdepth (all reads, no quality filters). Used for depth filtering. |
| `CALLER_COUNT` | int | Number of independent callers that called this variant (merged files only; blank in GATK-only files) |
| `CALLERS` | string | Comma-separated caller names that called this variant (merged files only) |

---

## Somatic Classification Tiers

Applied by `scripts/classify_variants.py` to each VEP-annotated VCF:

| Class | Criteria | Notes |
|---|---|---|
| `SOMATIC_KNOWN` | COSMIC match (any gnomAD AF) **or** ClinVar pathogenic + gnomAD AF < 0.001 | COSMIC unconditionally overrides gnomAD — somatic mutations are absent from population databases by definition. ClinVar alone defers to gnomAD (BRCA1/CFTR are pathogenic germline variants). |
| `NOVEL_CANDIDATE` | gnomAD AF = 0, QUAL ≥ 30, AD ≥ 5 | Absent from all gnomAD populations. Could be somatic or ultra-rare germline. |
| `RARE_CANDIDATE` | gnomAD AF < 0.0001, AD ≥ 5 | Very rare but present in population. |
| `UNCERTAIN` | gnomAD AF 0.0001–0.001, AD ≥ 5 | Grey zone — rare but measurable in population. |
| `GERMLINE` | gnomAD AF ≥ 0.001 **or** AD < 5 | Excluded from somatic analysis. |

**Decision priority**: SOMATIC_KNOWN > gnomAD filter > NOVEL > RARE > UNCERTAIN

---

## Filtering Pipeline

```
All variants per caller
        │
        ▼
VEP annotation (gnomAD, dbSNP, COSMIC, ClinVar, SIFT/PolyPhen)
        │
        ▼
classify_variants.py
  — gnomAD AF ≥ 0.001 OR AD < 5  →  GERMLINE bucket (excluded)
  — COSMIC / ClinVar pathogenic   →  SOMATIC_KNOWN bucket
  — gnomAD AF = 0, QUAL ≥ 30     →  NOVEL_CANDIDATE bucket
  — gnomAD AF < 0.0001            →  RARE_CANDIDATE bucket
  — gnomAD AF 0.0001–0.001        →  UNCERTAIN bucket
        │
        ▼
merge_callers.py
  — Merges SOMATIC_KNOWN + NOVEL + RARE + UNCERTAIN from all 4 callers
  — CALLER_COUNT ≥ 3  →  high_confidence
  — CALLER_COUNT = 2  →  medium_confidence
  — CALLER_COUNT = 1  →  low_confidence
        │
        ▼
mosdepth BAM_DP filter (BAM_DP > 1500)
  — Removes low-coverage positions
  — Uses raw BAM depth (all reads), NOT GATK FORMAT/DP
  — Critical: GATK FORMAT/DP = ~3% of true depth for PacBio data
        │
        ▼
add_bam_annotations.py
  — Adds GENE (Hugo symbol, IMPACT-priority canonical transcript)
  — Adds BAM_REF, BAM_ALT, BAM_VAF via samtools mpileup (no quality filters)
        │
        ▼
vcf_to_table.py  →  *_final.tsv  (human-readable)
```

---

## Why BAM_DP ≠ GATK FORMAT/DP

GATK HaplotypeCaller excludes:
1. Reads with MQ=255 (PacBio's code for "mapping quality not computed" on splice-junction reads)
2. Soft-clipped bases (`--dont-use-soft-clipped-bases`)

For PacBio Iso-Seq data where most reads span splice junctions, this discards **~97% of reads**.

**Example — SRSF2 (chr17:76736877):**

| Source | Depth |
|---|---|
| mosdepth (all BAM reads) | 5,548 |
| GATK FORMAT/DP | 162 (2.9%) |
| BAM_REF (mpileup, no filters) | 2,859 |
| BAM_ALT (mpileup, no filters) | 2,688 |
| BAM_VAF | 0.485 (48.5%) |

Filtering on `FORMAT/DP > 1500` would incorrectly discard SRSF2. The pipeline uses `BAM_DP > 1500` from mosdepth instead.

---

## Truth Set Detection Status

8 known somatic variants from clinical report (sample 5801-diagnosis):

| Gene | Position (GRCh38) | Mutation | Expected VAF | Detected | Callers | Reason if missing |
|---|---|---|---|---|---|---|
| **SRSF2** | chr17:76736877 | p.Pro95Arg (missense) | 37% | ✅ In high_confidence_final | GATK, DV, Clair3, LongcallR | — |
| **SETBP1** | chr18:44951948 | p.Gly870Ser (missense) | 28% | ⚠️ Called, filtered | GATK, DV, Clair3, LongcallR | BAM_DP ≤ 1500 at this position (low RNA coverage) |
| **RUNX1** | chr21:34799432 | p.Trp279Ter (nonsense) | 35% | ⚠️ Called, filtered | GATK, DV, Clair3, LongcallR | BAM_DP ≤ 1500 (GATK AD=14,6; Clair3 AD=468,263) |
| **ASXL1** | chr20:32434638 | Frameshift | 8% | ⚠️ Called, filtered | GATK, DV, Clair3 | BAM_DP ≤ 1500 + VEP annotates as synonymous (coordinate may differ from frameshift) |
| **CBL** | chr11:119242488 | Frameshift | UNK | ❌ Not callable | — | Variant is in a 177 kb intron. Iso-Seq reads skip introns (CIGAR `N`). Fundamentally undetectable by RNA-seq. |
| **MEF2B** | chr19:19145882 | Missense | UNK | ❌ Not detected | — | Not called by any caller; likely insufficient RNA expression or coverage at this position |
| **BCOR** | chrX:40057210 | Nonsense | UNK | ❌ Not detected | — | Not called by any caller |
| **PHF6** | chrX:134415107 | Missense | UNK | ❌ Not detected | — | Not called by any caller |

**Summary**: 4/8 variants callable from RNA-seq (CBL is intronic; 3 others have insufficient coverage after BAM_DP>1500 filter). SRSF2 fully validated with all 4 callers and 48.5% BAM VAF.

To recover SETBP1/RUNX1/ASXL1: lower the BAM_DP threshold (e.g. >200) or examine `medium_confidence_final` files without the depth filter.

---

## Key Scripts

| Script | Input → Output | Description |
|---|---|---|
| `scripts/run_callers.sh` | BAM → raw VCFs | Run all 4 callers on `scisoseq_synthetic_qual.bam` |
| `scripts/annotate_and_filter.sh` | raw VCF → annotated + classified | VEP annotation + classify_variants.py |
| `scripts/classify_variants.py` | annotated VCF → 5 buckets | Tiered somatic classification |
| `scripts/merge_callers.py` | 4× classified VCFs → merged tiers | Multi-caller concordance |
| `scripts/add_bam_annotations.py` | VCF + BAM → VCF with GENE/BAM counts | Adds GENE, BAM_REF, BAM_ALT, BAM_VAF |
| `scripts/vcf_to_table.py` | `*_final.vcf.gz` → TSV | Human-readable flat table |

---

## Known Limitations

1. **No matched normal**: Rare germline variants (gnomAD AF < 0.001) are indistinguishable from somatic without paired DNA. COSMIC and ClinVar are the primary somatic evidence.

2. **RNA editing (A→I)**: ADAR editing produces A>G variants mimicking real mutations. A future filter with REDIportal (~4.5M known sites) would reduce this noise. Estimate: some fraction of NOVEL_CANDIDATE variants are editing artefacts.

3. **Intronic variants**: RNA-seq cannot detect mutations in introns. Iso-Seq reads skip intronic sequence entirely (CIGAR `N` operations). CBL is the affected truth variant here.

4. **BAM_DP > 1500 filter**: Conservative threshold chosen to match GATK's high-confidence calls. Reduces noise but excludes low-expressed genes (SETBP1, RUNX1, ASXL1 at these specific positions).

5. **GATK FORMAT/DP vs BAM_DP**: All depth filtering in this pipeline uses mosdepth `BAM_DP`, not GATK's `FORMAT/DP`. Never filter on `FORMAT/DP` for PacBio data.

6. **LongcallR SNP-only**: LongcallR cannot call indels. ASXL1 frameshift and other indels are not supported by LongcallR's CALLER_COUNT.
