# RNA-seq Somatic Variant Calling Pipeline

Detection of somatic mutations in hematologic malignancies using PacBio Iso-Seq (KINNEX) long-read single-cell RNA-seq data. Four independent variant callers are run per sample, variants are annotated with population databases and functional consequences via VEP, and a multi-caller concordance filter isolates high-confidence somatic candidates with BAM-level allele frequency estimates.

---

## Cohort — 8 Samples

MDS/AML patients from the Grimes KINNEX dataset, processed against GRCh38 (chr-prefixed).

| # | Sample | Status |
|---|--------|--------|
| 1 | MDS5801-1229IR_CD34 | ✅ Complete |
| 2 | MDS5801-postHMA_CD34 | ✅ Complete |
| 3 | MDS5801-preVEN_342__1 | ✅ Complete |
| 4 | MDS5801-preVEN_342__2 | ✅ Complete |
| 5 | MDS5801-preVEN_342__3 | ✅ Complete |
| 6 | MDS5801R-341 | ✅ Complete |
| 7 | WM71-0121-5801_CD34 | ✅ Complete |
| 8 | XTX14_2_MDS_5801_CD34 | ✅ Complete |

For per-sample variant counts, somatic classification breakdowns, and run notes see **[RESULTS_REPORT.md](RESULTS_REPORT.md)**.

---

## Pipeline Overview

Each sample passes through a 10-step checkpoint-aware pipeline (`scripts/run_sample_pipeline.sh`). Steps are automatically skipped if their primary output already exists, enabling safe crash recovery and reruns.

```
BAM  (bams/{SAMPLE}.bam)
  │
  ▼ Step 1 — add_synthetic_quality.py
  │   Assigns QUAL=30 to all reads.  Required by GATK WellformedReadFilter;
  │   PacBio Iso-Seq reads carry QUAL=* which GATK rejects.
  │
  ├─► Step 2 — GATK HaplotypeCaller   (canonical chr1–chrM)
  ├─► Step 3 — DeepVariant 1.6.1      (PACBIO model)
  ├─► Step 4 — Clair3-RNA             (hifi model bundled in container)
  └─► Step 5 — LongcallR 1.12.0       (hifi-isoseq preset, SNP-only)
              │  VCFs sorted + indexed after each caller
              ▼
  Step 6 — VEP annotation (Ensembl 112) per caller
              │  gnomAD AF · ClinVar · COSMIC · HGVS · canonical transcript
              │  consequence · SIFT · PolyPhen
              ▼
  Step 7 — classify_variants.py + merge_callers.py
              │  Per-caller germline/somatic separation
              │  Multi-caller merge → high / medium / low confidence tiers
              ▼
  Step 8 — mosdepth
              │  True BAM read depth at each variant position
              │  (bypasses GATK's MAPQ=255 exclusion — see note below)
              ▼
  Step 9 — add_bam_annotations.py
              │  Per-variant BAM re-interrogation via pysam:
              │  BAM_REF · BAM_ALT · BAM_VAF
              ▼
  Step 10 — vcf_to_table.py
              TSV export with best-consequence GENE assignment
```

### Why BAM-level allele counts (Steps 8–9)?

PacBio Iso-Seq reads at splice junctions carry `MAPQ=255` ("mapping quality not computed"), which GATK silently excludes when computing `FORMAT/DP`. For KINNEX data this means GATK's native depth typically represents **2–5% of the true BAM depth**. All depth values and allele frequencies in final outputs are derived directly from the BAM using mosdepth (depth) and pysam (allele counts).

**Example — SRSF2 chr17:76736877 (WM71-0121-5801_CD34):**

| Source | Depth | VAF |
|---|---|---|
| Raw BAM (mosdepth) | 5,548 | 48.5% |
| GATK FORMAT/DP | 162 (2.9%) | — |

---

## Running the Pipeline

### Full cohort (batch mode)

```bash
cd /data/salomonis-archive/FASTQs/NCI-R01/rna_seq_varcall
bash scripts/run_all_samples.sh
```

`run_all_samples.sh` implements **pipeline-stage parallelism**: as soon as sample N's callers (Steps 1–5) finish, sample N+1's callers start in the background — VEP + annotation of sample N overlaps with callers of sample N+1. Controlled by `MAX_CALLER_JOBS=1` (line 34 of the script); increase to 2 for more concurrency on the 753 GB / 56-core server.

### Single sample

```bash
bash scripts/run_sample_pipeline.sh bams/{SAMPLE}.bam
```

### Logs

```
logs/run_all_samples.log         batch runner master log
logs/{SAMPLE}/pipeline.log       per-sample step-by-step log
```

---

## Output Structure

```
results/{SAMPLE}/
├── {SAMPLE}_qual.bam                       Step 1: synthetic-quality BAM
├── haplotypecaller/variants.vcf.gz         Step 2: GATK normalized VCF
├── deepvariant/variants.vcf.gz             Step 3: DeepVariant VCF
├── clair3_rna/variants.vcf.gz              Step 4: Clair3 VCF
├── longcallr/variants.vcf.gz               Step 5: LongcallR VCF (sorted)
└── merged/
    ├── high_confidence_final.vcf.gz        ≥3 callers, BAM_DP>1500, annotated
    ├── high_confidence_final.tsv           ← PRIMARY OUTPUT
    ├── medium_confidence_final.vcf.gz      2 callers, BAM_DP>1500, annotated
    ├── medium_confidence_final.tsv
    ├── low_confidence_final.vcf.gz         1 caller, BAM_DP>1500, annotated
    └── low_confidence_final.tsv

annotation/{SAMPLE}/gatk/
├── germline.vcf.gz                         GATK germline variants
├── somatic_known.vcf.gz                    ClinVar pathogenic / COSMIC hits
├── novel_candidates.vcf.gz                 Not in gnomAD
├── rare_candidates.vcf.gz                  gnomAD AF < 0.01
├── uncertain.vcf.gz                        Conflicting / unknown significance
├── non_germline_final.vcf.gz               All non-germline, BAM_DP>1500
└── non_germline_final.tsv                  ← PRIMARY GATK OUTPUT
```

---

## Variant Callers

| Caller | Environment | Model | SNPs | Indels |
|---|---|---|---|---|
| GATK HaplotypeCaller | `gatk-env` conda | RNA mode (`--dont-use-soft-clipped-bases`) | ✓ | ✓ |
| DeepVariant 1.6.1 | `containers/deepvariant_1.6.1.sif` | PACBIO | ✓ | ✓ |
| Clair3-RNA | `containers/clair3_latest.sif` | `/opt/models/hifi` (bundled) | ✓ | ✓ |
| LongcallR 1.12.0 | `containers/longcallr_1.12.0.sif` | hifi-isoseq preset | ✓ | ✗ |

> **GATK** is restricted to canonical chromosomes (`chr1–chrM`) to avoid crashes on the alt/random/patch contigs present in PacBio BAM headers that are absent from the reference dict.

---

## Somatic Classification (Step 7)

`scripts/classify_variants.py` assigns each GATK variant to one of five tiers:

| Class | Criteria |
|---|---|
| **GERMLINE** | gnomAD AF ≥ 0.01 or ClinVar benign/likely-benign |
| **SOMATIC_KNOWN** | COSMIC hit or ClinVar pathogenic/likely-pathogenic |
| **RARE_CANDIDATE** | gnomAD AF present but < 0.01 |
| **NOVEL_CANDIDATE** | Not in gnomAD (AF = 0 or absent) |
| **UNCERTAIN** | Conflicting or unknown significance |

---

## Multi-Caller Confidence Tiers (Step 7)

`scripts/merge_callers.py` merges variants across all four callers using `bcftools isec`:

| Tier | Criterion | Use |
|---|---|---|
| **High** | ≥ 3 callers | Highest specificity — primary analysis set |
| **Medium** | 2 callers | Moderate confidence |
| **Low** | 1 caller | Sensitive, exploratory |

All final outputs require `BAM_DP > 1500`. The 1,500-read threshold is appropriate for KINNEX data, which typically yields 5,000–200,000× coverage at expressed loci.

---

## TSV Column Descriptions

| Column | Description |
|---|---|
| `CHROM`, `POS`, `REF`, `ALT` | Variant coordinates (GRCh38, 1-based) |
| `GENE` | Best-consequence gene (canonical transcript preferred; consequence sub-ranked within MODIFIER tier to correctly assign non-coding RNA genes) |
| `CONSEQUENCE` | Most severe VEP consequence for selected transcript |
| `HGVSp` | Protein-level change |
| `HGVSc` | cDNA-level change |
| `COSMIC` | COSMIC variant ID(s) |
| `CLINVAR` | ClinVar significance |
| `gnomAD_AF_max` | Maximum population AF across gnomAD genomes + exomes |
| `SOMATIC_CLASS` | `SOMATIC_KNOWN`, `NOVEL_CANDIDATE`, `RARE_CANDIDATE`, `UNCERTAIN`, `GERMLINE` |
| `SOMATIC_REASON` | Evidence supporting classification |
| `BAM_REF` | Reference allele count from raw BAM (pysam) |
| `BAM_ALT` | Alternate allele count from raw BAM (pysam) |
| `BAM_VAF` | BAM_ALT / (BAM_REF + BAM_ALT) |
| `BAM_DP` | Total read depth at position (mosdepth) |
| `CALLER_COUNT` | Number of callers supporting the variant (merged outputs only) |
| `CALLERS` | Comma-separated list of supporting callers (merged outputs only) |

---

## Scripts Reference

| Script | Purpose |
|---|---|
| `run_all_samples.sh` | Batch runner with pipeline-stage parallelism over `bams/*.bam` |
| `run_sample_pipeline.sh` | 10-step checkpoint-aware single-sample pipeline |
| `add_synthetic_quality.py` | Assigns QUAL=30 to all reads (PacBio → GATK compatibility) |
| `classify_variants.py` | Tiered somatic classification using gnomAD + ClinVar + COSMIC |
| `merge_callers.py` | Multi-caller concordance merge → high/medium/low confidence VCFs |
| `add_bam_annotations.py` | Per-variant BAM re-interrogation: BAM_REF, BAM_ALT, BAM_VAF via pysam |
| `vcf_to_table.py` | Exports annotated VCF to TSV with best-consequence GENE selection |
| `annotate_and_filter.sh` | VEP annotation wrapper (called by run_sample_pipeline.sh) |
| `run_callers.sh` | Standalone caller launcher (prototype / testing) |
| `evaluate_callers.py` | Caller performance evaluation against truth set |
| `compare_callers.py` | Cross-caller concordance and overlap analysis |

---

## Key Biological Caveats

1. **No matched normal**: Without paired germline DNA, rare germline variants are indistinguishable from somatic mutations by population frequency alone. COSMIC and ClinVar are the primary evidence for somatic origin.

2. **RNA editing (A→I)**: ADAR-mediated RNA editing produces A→G variants that appear as real mutations. A future filter using REDIportal (~4.5 M known editing sites) would substantially reduce this false-positive class, particularly in NOVEL_CANDIDATE calls.

3. **Intronic variants undetectable**: Iso-Seq reads skip introns (CIGAR `N` operations). Mutations in purely intronic regions — including splice-site variants that do not also affect an exon — cannot be detected from RNA-seq data regardless of caller or depth.

4. **Mitochondrial variants**: gnomAD does not comprehensively index the mitochondrial genome. High-VAF chrM variants (VAF near 1.0) classified as NOVEL_CANDIDATE are most likely haplogroup germline polymorphisms, not somatic mutations.

5. **MAPQ=255 and GATK depth**: All depth and VAF values in final outputs are from mosdepth / pysam. Do not use GATK `FORMAT/DP` or `FORMAT/AD` for PacBio data — they reflect only a small fraction of the true read depth.

---

## Software Versions

| Tool | Version | Location |
|---|---|---|
| GATK HaplotypeCaller | 4.6.2.0 | `gatk-env` conda env |
| DeepVariant | 1.6.1 (PACBIO model) | `containers/deepvariant_1.6.1.sif` |
| Clair3 | latest | `containers/clair3_latest.sif` — model `/opt/models/hifi` |
| LongcallR | 1.12.0 | `containers/longcallr_1.12.0.sif` |
| VEP | Ensembl release 112 | `containers/ensembl-vep_release_112.0.sif` |
| mosdepth | — | `/usr/local/anaconda3-2020/envs/mosdepth/bin/mosdepth` |
| samtools / bcftools / bgzip | 1.23 | `bio-cli` conda env |
| singularity | — | `bio-cli` conda env |
| Python | 3.x (pysam, pandas) | `bio-cli` conda env |
| Reference | GRCh38 chr-prefixed | `reference/genome.fa` → symlink from spaceranger |

---

## Future Directions

### AlphaMissense

[AlphaMissense](https://www.science.org/doi/10.1126/science.adg7492) (DeepMind, 2023) provides pre-computed pathogenicity scores for all ~71 M possible human missense variants using AlphaFold2 protein structure context. Adding an `AlphaMissense_score` and `AlphaMissense_class` (likely_benign / ambiguous / likely_pathogenic) column to the final TSVs would substantially improve somatic candidate prioritization for the ~30% of high-confidence variants that are missense — providing structural evidence independent of COSMIC / ClinVar curation. The pre-computed hg38 scores (2.1 GB TSV) are publicly available and can be integrated as a VEP plugin or a post-annotation join in `vcf_to_table.py`.

### AlphaGenome

[AlphaGenome](https://deepmind.google/discover/blog/alphagenome-a-new-model-for-interpreting-the-genome/) (Google DeepMind, 2025) predicts genome-wide molecular phenotypes — gene expression levels, splicing, chromatin accessibility, and transcription factor occupancy — directly from sequence. Relevant integration points for this pipeline:

- **In silico functional impact**: Score each NOVEL_CANDIDATE and RARE_CANDIDATE for predicted effects on splicing (exon skipping, cryptic splice-site activation) and regulatory activity, providing functional evidence where COSMIC and ClinVar have no annotation.
- **Context-aware variant prioritization**: Variants in predicted regulatory regions (active enhancers, CTCF binding sites) or at high predicted splice-strength positions can be up-ranked; variants in predicted inactive chromatin can be down-ranked to reduce false positives, especially relevant for the high rate of NOVEL_CANDIDATE calls in this RNA-seq dataset.
- **Expression-aware interpretation**: AlphaGenome's expression predictions can flag variants likely to alter the expression level of target genes — a dimension of functional impact not captured by consequence annotation alone.

### RNA editing filter (REDIportal)

A lookup against [REDIportal v2](http://reditools.com/reditoolsDB/) (~4.5 M A→I editing sites catalogued across human tissues) would flag all A→G (and T→C on minus strand) variants overlapping known RNA editing sites. This would substantially reduce the NOVEL_CANDIDATE burden, which is currently enriched for edited sites in highly expressed genes.

### Splice-variant scoring (SpliceAI)

[SpliceAI](https://github.com/Illumina/SpliceAI) delta scores (pre-computed for all positions ±50 nt from annotated splice sites in GRCh38) would flag splice-disrupting SNVs and small indels that the current pipeline classifies as synonymous or intronic. This is particularly relevant for long-read RNA-seq where splicing isoforms are directly observed — a splice-disrupting somatic variant may produce a novel isoform detectable in the same data.

### Clonal evolution across timepoints

Multiple samples from the same patient (MDS5801 pre-VEN ×3 samples, post-HMA, and diagnosis) provide a natural framework for tracking clonal dynamics. A cross-sample merge on shared variant positions with VAF trajectories would identify clones expanding or contracting under venetoclax and HMA therapy — directly supporting the clinical narrative of the cohort.

### Tumor cell fraction–aware VAF correction

VAF in bulk RNA-seq reflects both tumor cell fraction and allele-specific expression. Integration with cell-type deconvolution using the paired scRNA-seq barcodes available in the KINNEX data would allow estimation of tumor cell–specific VAF, separating clonal from subclonal mutations and improving comparability of VAF across timepoints and samples with different tumor purity.

---

*Last updated: 2026-03-03*
