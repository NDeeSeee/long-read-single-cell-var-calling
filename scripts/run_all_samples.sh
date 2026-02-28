#!/usr/bin/env bash
# =============================================================================
# run_all_samples.sh
#
# Run the full variant calling pipeline on every BAM in bams/.
# Samples are processed sequentially to avoid overloading the server.
# Each sample is independently checkpointed — re-running is safe.
#
# Usage:
#   bash scripts/run_all_samples.sh [bams_dir]
#
# Default bams_dir: /data/salomonis-archive/FASTQs/NCI-R01/rna_seq_varcall/bams
# =============================================================================

set -euo pipefail

WORKDIR="/data/salomonis-archive/FASTQs/NCI-R01/rna_seq_varcall"
BAMS_DIR="${1:-$WORKDIR/bams}"
PIPELINE="$WORKDIR/scripts/run_sample_pipeline.sh"
MASTER_LOG="$WORKDIR/logs/run_all_samples.log"

mkdir -p "$WORKDIR/logs"

echo "==========================================================" | tee -a "$MASTER_LOG"
echo "run_all_samples.sh started: $(date)" | tee -a "$MASTER_LOG"
echo "BAMs directory: $BAMS_DIR" | tee -a "$MASTER_LOG"
echo "==========================================================" | tee -a "$MASTER_LOG"

BAM_LIST=("$BAMS_DIR"/*.bam)
TOTAL=${#BAM_LIST[@]}
echo "Found $TOTAL BAMs:" | tee -a "$MASTER_LOG"
for bam in "${BAM_LIST[@]}"; do echo "  $(basename $bam)"; done | tee -a "$MASTER_LOG"
echo "" | tee -a "$MASTER_LOG"

PASS=0
FAIL=0
FAILED_SAMPLES=()

for i in "${!BAM_LIST[@]}"; do
    bam="${BAM_LIST[$i]}"
    sample="$(basename "$bam" .bam)"
    n=$((i + 1))

    echo "----------------------------------------------------------" | tee -a "$MASTER_LOG"
    echo "[$n/$TOTAL] Starting: $sample  ($(date))" | tee -a "$MASTER_LOG"
    echo "----------------------------------------------------------" | tee -a "$MASTER_LOG"

    mkdir -p "$WORKDIR/logs/${sample}"
    if bash "$PIPELINE" "$bam" 2>&1 | tee -a "$WORKDIR/logs/${sample}/pipeline.log"; then
        echo "[$n/$TOTAL] SUCCESS: $sample  ($(date))" | tee -a "$MASTER_LOG"
        PASS=$((PASS + 1))
    else
        echo "[$n/$TOTAL] FAILED:  $sample  ($(date))" | tee -a "$MASTER_LOG"
        FAIL=$((FAIL + 1))
        FAILED_SAMPLES+=("$sample")
        echo "  See: $WORKDIR/logs/${sample}/" | tee -a "$MASTER_LOG"
        echo "  Continuing with next sample..." | tee -a "$MASTER_LOG"
    fi
    echo "" | tee -a "$MASTER_LOG"
done

echo "==========================================================" | tee -a "$MASTER_LOG"
echo "run_all_samples.sh finished: $(date)" | tee -a "$MASTER_LOG"
echo "  Completed: $PASS / $TOTAL" | tee -a "$MASTER_LOG"
echo "  Failed:    $FAIL / $TOTAL" | tee -a "$MASTER_LOG"
if [ ${#FAILED_SAMPLES[@]} -gt 0 ]; then
    echo "  Failed samples:" | tee -a "$MASTER_LOG"
    for s in "${FAILED_SAMPLES[@]}"; do echo "    - $s"; done | tee -a "$MASTER_LOG"
fi
echo "==========================================================" | tee -a "$MASTER_LOG"
