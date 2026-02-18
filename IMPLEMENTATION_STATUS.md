# Variant Caller Comparison - Implementation Status

## Summary
Implementation of a four-caller variant calling comparison pipeline for PacBio Iso-Seq RNA-seq data from patient 5801-diagnosis.

**Goal**: Identify best variant caller (HaplotypeCaller, DeepVariant, Clair3-RNA, LongcallR) for production use on long-read RNA-seq data.

---

## Completed Phases

### Phase 1: Environment Setup ✓
- ✓ samtools 1.23 (bio-cli)
- ✓ bcftools 1.23 (bio-cli)
- ✓ singularity (bio-cli: `/users/pavb5f/.conda/envs/bio-cli/bin/singularity`)
- ✓ Reference genome symlink: `$WORKDIR/reference/genome.fa` → spaceranger GRCh38 (chr-prefixed)
- ✓ Test BAM (11GB, indexed): `/data/salomonis-archive/BAMs/Grimes/scRNA-Seq/KINNEX/5801-diagnosis/scisoseq.mapped.bam`
- ✓ DeepVariant container (2.7 GB): `$WORKDIR/containers/deepvariant_1.6.1.sif`

### Phase 2: Truth Set Preparation ✓
- ✓ `truth_set/truth_set.vcf.gz` — 8 somatic variants, compressed + indexed
- ✓ `truth_set/truth_regions.bed` — target regions ±1kb around variants (chr-prefixed)
- ✓ `truth_set/test_variants_formatted.txt` — formatted for lab scripts

**Truth Set Details** (8 known variants):
| Gene | Position | Protein Change | VAF | Type |
|------|----------|----------------|-----|------|
| RUNX1 | chr21:34799432 | p.W279* | 35% | Nonsense |
| ASXL1 | chr20:32434638 | p.G643fs | 8% | Frameshift |
| SETBP1 | chr18:44951948 | p.G870S | 28% | Missense |
| SRSF2 | chr17:76736877 | p.P95R | 37% | Missense |
| CBL | chr11:119242488 | p.S606fs | UNK | Frameshift |
| MEF2B | chr19:19145882 | p.L341V | UNK | Missense |
| BCOR | chrX:40057210 | p.R1480* | UNK | Nonsense |
| PHF6 | chrX:134415107 | p.R274P | UNK | Missense |

---

## Current State (as of 2026-02-18)

### Tool Status
| Tool | Status | Notes |
|------|--------|-------|
| samtools 1.23 | ✅ working | bio-cli |
| bcftools 1.23 | ✅ working | bio-cli |
| GATK 4.6.2.0 | ❌ broken in bio-cli | `java: undefined symbol: JLI_StringDup` — Java/glibc mismatch |
| DeepVariant SIF | ✅ present | Container at `$WORKDIR/containers/deepvariant_1.6.1.sif` |
| Clair3 | ❌ not installed | Needs `clair3-rna` conda env |
| LongcallR | ❌ not installed | Needs `longcallr` conda env |
| singularity | ✅ working | `/users/pavb5f/.conda/envs/bio-cli/bin/singularity` |

---

## In-Progress Phases

### Phase 3: Install Tools & Run Callers

#### Step 1: Fix GATK — Create `gatk-env`
```bash
conda create -n gatk-env -c conda-forge -c bioconda openjdk=17 gatk4=4.6.2.0 -y
/users/pavb5f/.conda/envs/gatk-env/bin/gatk --version
```

Run HaplotypeCaller:
```bash
bash scripts/run_callers.sh haplotypecaller
```
Output: `results/haplotypecaller/variants.vcf.gz`

#### Step 2: DeepVariant (container already present)
```bash
bash scripts/run_callers.sh deepvariant
```
- Uses `--model_type=PACBIO` (correct for PacBio HiFi reads)
- Output: `results/deepvariant/variants.vcf.gz`

#### Step 3: Install Clair3-RNA
```bash
bash scripts/setup_envs.sh clair3
```
Then update `CLAIR3_MODEL` in `run_callers.sh` to the downloaded model path, then:
```bash
bash scripts/run_callers.sh clair3_rna
```
Output: `results/clair3_rna/merge_output.vcf.gz`

#### Step 4: Install LongcallR
```bash
bash scripts/setup_envs.sh longcallr
bash scripts/run_callers.sh longcallr
```
Output: `results/longcallr/variants.vcf.gz`

#### Or run all at once (after all envs are ready):
```bash
bash scripts/setup_envs.sh all
bash scripts/run_callers.sh all
```

---

## Planned Phases

### Phase 4: Evaluate All 4 Callers
```bash
python scripts/evaluate_callers.py
```
Calculates: sensitivity, precision, F1, VAF correlation, MAE.
Output: `results/comparison/summary.json`, `results/comparison/comparison.csv`

### Phase 5: Generate Report
```bash
python scripts/generate_report.py
```
Output: `results/comparison/COMPARISON_REPORT.md`

---

## Project Structure

```
WORKDIR=/data/salomonis-archive/FASTQs/NCI-R01/rna_seq_varcall/
├── reference/
│   └── genome.fa                       # symlink → spaceranger GRCh38 (chr-prefixed)
├── truth_set/
│   ├── truth_set.vcf.gz                # ✓ 8 somatic variants
│   ├── truth_set.vcf.gz.tbi            # ✓ index
│   ├── truth_regions.bed               # ✓ ±1kb target regions
│   └── test_variants_formatted.txt     # ✓ formatted variants
├── models/
│   └── clair3_rna/                     # PacBio HiFi Clair3 model (to download)
├── containers/
│   └── deepvariant_1.6.1.sif           # ✓ present (2.7 GB)
├── results/
│   ├── haplotypecaller/                # GATK output
│   ├── deepvariant/                    # DeepVariant output
│   ├── clair3_rna/                     # Clair3-RNA output
│   ├── longcallr/                      # LongcallR output
│   └── comparison/                     # Evaluation results
└── logs/                               # Per-caller log files

GitHub repo (scripts tracked here):
├── scripts/
│   ├── run_callers.sh                  # ✓ Updated for 4 new callers
│   ├── setup_envs.sh                   # ✓ Conda env setup
│   ├── evaluate_callers.py             # ✓ Updated caller paths
│   ├── generate_report.py              # ✓ Fixed f-string bug, updated callers
│   ├── tsv_to_vcf.py
│   ├── compress_vcf.py
│   └── prepare_truth_vcf.py
```

---

## Conda Environments

| Env | Tool | Create Command |
|-----|------|----------------|
| `gatk-env` | GATK 4.6.2.0 | `conda create -n gatk-env -c conda-forge -c bioconda openjdk=17 gatk4=4.6.2.0 -y` |
| `bio-cli` | samtools, bcftools, singularity | Already exists |
| `clair3-rna` | Clair3 | `conda create -n clair3-rna -c bioconda -c conda-forge python=3.9 clair3 -y` |
| `longcallr` | LongcallR | `conda create -n longcallr -c bioconda -c conda-forge longcallr -y` |

---

## Key Notes

- **MAPQ=0**: 43% of reads are multi-mappers → reduces SNV sensitivity
- **Read Length**: ~802bp average (PacBio long reads)
- **Cell Count**: ~8-10K cells per sample
- **VAF Known**: 4 of 8 variants (8%, 28%, 35%, 37%)
- **Frameshifts**: 2 variants (ASXL1, CBL) — require indel-capable callers
- **Reference**: Must use chr-prefixed GRCh38 to match BAM headers

---

**Last Updated**: 2026-02-18
**Status**: Phase 2 Complete, Phase 3 In Progress (installing envs)
