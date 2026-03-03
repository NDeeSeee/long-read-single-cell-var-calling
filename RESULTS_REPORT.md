# RNA-seq Variant Calling Pipeline — Results Report

**Date**: 2026-03-03
**Cohort**: Grimes scRNA-Seq KINNEX dataset — MDS/AML patients
**Reference**: GRCh38 (chr-prefixed)
**Pipeline**: `scripts/run_sample_pipeline.sh` / `scripts/run_all_samples.sh`

---

## 1. Overview

This pipeline performs somatic variant calling on PacBio Iso-Seq (KINNEX) long-read RNA-seq BAMs. Four independent variant callers are run per sample, followed by VEP annotation, somatic classification, germline/somatic separation, BAM-level allele counting, and multi-caller consensus filtering.

### Cohort — 8 Samples

| # | Sample | Status |
|---|--------|--------|
| 1 | MDS5801-1229IR_CD34 | ✅ Complete |
| 2 | MDS5801-postHMA_CD34 | ✅ Complete |
| 3 | MDS5801-preVEN_342__1 | ✅ Complete |
| 4 | MDS5801-preVEN_342__2 | ✅ Complete |
| 5 | MDS5801-preVEN_342__3 | 🔄 In progress — VEP annotation (Step 6) |
| 6 | MDS5801R-341 | 🔄 In progress — Synthetic quality preprocessing (Step 1) |
| 7 | WM71-0121-5801_CD34 | ⏳ Queued |
| 8 | XTX14_2_MDS_5801_CD34 | ⏳ Queued |

---

## 2. Pipeline Steps

Each sample passes through 10 checkpoint-aware steps. If a step's primary output already exists, it is skipped on re-run.

| Step | Tool | Description |
|------|------|-------------|
| 1 | Custom Python | Synthetic Q30 preprocessing — assigns QUAL=30 to all reads (required by GATK's WellformedReadFilter, which rejects PacBio QUAL=* reads) |
| 2 | GATK HaplotypeCaller | SNV/indel calling restricted to canonical chromosomes (chr1–chrM); output normalized with `bcftools norm` |
| 3 | DeepVariant (PACBIO model) | Deep learning–based SNV/indel calling |
| 4 | Clair3-RNA (hifi model) | Long-read RNA-seq SNV/indel calling; bundled `/opt/models/hifi` model inside container |
| 5 | LongcallR (hifi-isoseq preset) | Long-read isoform-aware variant calling; output sorted before indexing |
| 6 | VEP (Ensembl release 112) | Functional annotation per caller: gnomAD AF, ClinVar, COSMIC, HGVS, canonical transcript, consequence |
| 7 | Custom Python | Multi-caller merge: variants stratified into high (≥3 callers), medium (2 callers), low (1 caller) confidence tiers |
| 8 | mosdepth | BAM-level depth estimation; variants with `BAM_DP ≤ 1500` excluded from final outputs |
| 9 | Custom Python (`add_bam_annotations.py`) | Per-variant BAM re-interrogation: computes `BAM_REF`, `BAM_ALT`, `BAM_VAF` by querying the raw BAM at each variant position via pysam |
| 10 | Custom Python (`vcf_to_table.py`) | Export annotated VCF to TSV with best-consequence GENE assignment and all annotation fields |

### Why BAM-level allele counts (Steps 8–9)?

PacBio long reads frequently have MAPQ=255 (spliced junctions), which GATK excludes when computing FORMAT/DP. This causes GATK's native depth to represent only ~3% of the true read depth. All depth values and allele frequencies in the final outputs are therefore derived directly from the BAM using mosdepth (depth) and pysam (allele counts), not from the VCF FORMAT fields.

---

## 3. Output File Structure

All outputs are under the pipeline working directory:
`/data/salomonis-archive/FASTQs/NCI-R01/rna_seq_varcall/`

### Per-Sample Layout

```
results/{SAMPLE}/
├── {SAMPLE}_qual.bam                         Step 1: synthetic-quality BAM
├── haplotypecaller/
│   └── variants.vcf.gz                       Step 2: GATK normalized VCF
├── deepvariant/
│   └── variants.vcf.gz                       Step 3: DeepVariant VCF
├── clair3_rna/
│   └── variants.vcf.gz                       Step 4: Clair3 VCF
├── longcallr/
│   └── variants.vcf.gz                       Step 5: LongcallR VCF (sorted)
└── merged/
    ├── high_confidence_final.vcf.gz          ≥3 callers, BAM_DP>1500, annotated
    ├── high_confidence_final.tsv             ← PRIMARY OUTPUT
    ├── medium_confidence_final.vcf.gz        2 callers, BAM_DP>1500, annotated
    ├── medium_confidence_final.tsv
    ├── low_confidence_final.vcf.gz           1 caller, BAM_DP>1500, annotated
    └── low_confidence_final.tsv

annotation/{SAMPLE}/gatk/
├── germline.vcf.gz                           GATK germline variants (ClinVar/gnomAD classified)
├── germline.vcf.gz.tbi
├── somatic_known.vcf.gz                      Known somatic: ClinVar pathogenic / COSMIC hits
├── novel_candidates.vcf.gz                   Novel non-germline candidates
├── rare_candidates.vcf.gz                    Rare (gnomAD AF < 0.01) non-germline variants
├── uncertain.vcf.gz                          Uncertain significance
├── non_germline_final.vcf.gz                 All non-germline GATK variants, BAM_DP>1500, annotated
└── non_germline_final.tsv                    ← PRIMARY OUTPUT (GATK-only somatic candidates)
```

### TSV Column Descriptions

| Column | Description |
|--------|-------------|
| `CHROM`, `POS`, `REF`, `ALT` | Variant coordinates (GRCh38, 1-based) |
| `GENE` | Best-consequence gene assignment (canonical transcript preferred; sub-ranked by consequence severity within MODIFIER tier to correctly assign non-coding RNA genes) |
| `CONSEQUENCE` | Most severe VEP consequence for selected transcript |
| `HGVSp` | Protein-level change (if applicable) |
| `HGVSc` | cDNA-level change |
| `COSMIC` | COSMIC variant ID(s) |
| `CLINVAR` | ClinVar significance |
| `gnomAD_AF_max` | Maximum population AF across gnomAD genomes + exomes |
| `SOMATIC_CLASS` | Classification: `SOMATIC_KNOWN`, `NOVEL_CANDIDATE`, `RARE_CANDIDATE`, `UNCERTAIN`, `GERMLINE` |
| `SOMATIC_REASON` | Evidence supporting classification |
| `BAM_REF` | Reference allele count from raw BAM (pysam) |
| `BAM_ALT` | Alternate allele count from raw BAM (pysam) |
| `BAM_VAF` | BAM_ALT / (BAM_REF + BAM_ALT) — recalculated allele frequency |
| `BAM_DP` | Total raw depth at position (mosdepth) |
| `CALLER_COUNT` | Number of callers supporting the variant (merged outputs only) |
| `CALLERS` | Comma-separated list of supporting callers (merged outputs only) |

---

## 4. Results — Completed Samples

### 4.1 MDS5801-1229IR_CD34

| Output | Variants |
|--------|----------|
| Germline (GATK) | 1,101,834 |
| Non-germline GATK (BAM_DP > 1500) | 1,096 |
| High-confidence merged (≥3 callers) | **834** |
| Medium-confidence merged (2 callers) | 437 |
| Low-confidence merged (1 caller) | 1,314 |

**High-confidence somatic classification:**
- SOMATIC_KNOWN: 788 (94.5%) — ClinVar/COSMIC-supported
- NOVEL_CANDIDATE: 20 (2.4%)
- RARE_CANDIDATE: 15 (1.8%)
- UNCERTAIN: 11 (1.3%)

**Caller support (high-confidence):** 570 variants supported by 3 callers; 264 by all 4 callers

---

### 4.2 MDS5801-postHMA_CD34

| Output | Variants |
|--------|----------|
| Germline (GATK) | 1,180,063 |
| Non-germline GATK (BAM_DP > 1500) | 1,146 |
| High-confidence merged (≥3 callers) | **877** |
| Medium-confidence merged (2 callers) | 296 |
| Low-confidence merged (1 caller) | 576 |

**High-confidence somatic classification:**
- SOMATIC_KNOWN: 837 (95.4%)
- RARE_CANDIDATE: 17 (1.9%)
- NOVEL_CANDIDATE: 15 (1.7%)
- UNCERTAIN: 8 (0.9%)

**Caller support (high-confidence):** 592 variants by 3 callers; 285 by all 4 callers

---

### 4.3 MDS5801-preVEN_342__1

| Output | Variants |
|--------|----------|
| Germline (GATK) | 829,004 |
| Non-germline GATK (BAM_DP > 1500) | 563 |
| High-confidence merged (≥3 callers) | **457** |
| Medium-confidence merged (2 callers) | 180 |
| Low-confidence merged (1 caller) | 236 |

**High-confidence somatic classification:**
- SOMATIC_KNOWN: 446 (97.6%)
- NOVEL_CANDIDATE: 6 (1.3%)
- RARE_CANDIDATE: 5 (1.1%)

**Caller support (high-confidence):** 321 variants by 3 callers; 136 by all 4 callers

---

### 4.4 MDS5801-preVEN_342__2

| Output | Variants |
|--------|----------|
| Germline (GATK) | 739,255 |
| Non-germline GATK (BAM_DP > 1500) | 458 |
| High-confidence merged (≥3 callers) | **347** |
| Medium-confidence merged (2 callers) | 167 |
| Low-confidence merged (1 caller) | 217 |

**High-confidence somatic classification:**
- SOMATIC_KNOWN: 338 (97.4%)
- RARE_CANDIDATE: 4 (1.2%)
- NOVEL_CANDIDATE: 4 (1.2%)
- UNCERTAIN: 1 (0.3%)

**Caller support (high-confidence):** 245 variants by 3 callers; 102 by all 4 callers

---

## 5. Samples In Progress

### MDS5801-preVEN_342__3 (Sample 5)

All 4 caller VCFs are complete (Steps 1–5 done). Currently in **Step 6 (VEP annotation)** — annotating the Clair3 VCF (largest caller output, ~40M compressed). GATK and DeepVariant annotation already complete. Expected to finish within ~3–4 hours.

**Note:** This sample's VEP run was interrupted on 2026-03-01 (server process killed after ~45 min). The incomplete output was cleaned up and VEP was restarted with `--force_overwrite` (now default in the pipeline).

### MDS5801R-341 (Sample 6)

Started from Step 1 (synthetic quality preprocessing) on 2026-03-03 13:12. Running sequentially after sample 5 completes in the batch runner. Estimated total runtime: 20–24 hours.

---

## 6. Queued Samples

| Sample | Status |
|--------|--------|
| WM71-0121-5801_CD34 | Queued (sample 7) |
| XTX14_2_MDS_5801_CD34 | Queued (sample 8) |

These will be processed sequentially by `run_all_samples.sh` after samples 5 and 6 complete.

---

## 7. Key Design Notes

### Germline vs. Somatic Separation
GATK variants are classified by `classify_variants.py` using gnomAD population AF, ClinVar significance, and COSMIC annotations:
- **GERMLINE**: gnomAD AF ≥ 0.01 or ClinVar benign/likely-benign
- **SOMATIC_KNOWN**: COSMIC hit or ClinVar pathogenic/likely-pathogenic
- **RARE_CANDIDATE**: gnomAD AF present but < 0.01
- **NOVEL_CANDIDATE**: not in gnomAD (AF = 0 or absent)
- **UNCERTAIN**: conflicting or unknown significance

### Multi-caller Confidence Tiers
Variants are merged across all 4 callers using `bcftools isec`. The confidence tier reflects reproducibility:
- **High** (≥3 callers): highest specificity — primary analysis set
- **Medium** (2 callers): moderate confidence
- **Low** (1 caller): sensitive but less specific — for exploratory use

### BAM_DP Filter (mosdepth > 1500)
Applied to all final outputs. Removes variants at low-coverage positions where allele frequency estimation is unreliable. The 1500-read threshold is appropriate for PacBio KINNEX data, which typically yields 5,000–200,000× coverage at expressed loci.

### GENE Assignment
The best gene is selected from VEP's multi-transcript output using a two-tier ranking system:
1. VEP impact tier: HIGH > MODERATE > LOW > MODIFIER
2. Within MODIFIER: consequence sub-rank (e.g., `non_coding_transcript_exon_variant` > `upstream_gene_variant`) — this correctly assigns mitochondrial rRNA genes (MT-RNR1, MT-RNR2) that would otherwise be mis-assigned to neighboring protein-coding genes.

---

## 8. Software Versions

| Tool | Version / Notes |
|------|----------------|
| GATK HaplotypeCaller | conda env `gatk-env` |
| DeepVariant | v1.6.1 (PACBIO model), `deepvariant_1.6.1.sif` |
| Clair3 | `clair3_latest.sif`, model `/opt/models/hifi` |
| LongcallR | v1.12.0, `longcallr_1.12.0.sif` |
| VEP | Ensembl release 112, `ensembl-vep_release_112.0.sif` |
| mosdepth | `/usr/local/anaconda3-2020/envs/mosdepth/bin/mosdepth` |
| samtools / bcftools / bgzip | conda env `bio-cli` |
| Python | 3.x, conda env `bio-cli` (pysam, pandas) |
| Reference | GRCh38 (chr-prefixed, symlinked from spaceranger) |

---

*Report generated: 2026-03-03*
