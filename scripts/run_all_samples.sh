#!/usr/bin/env bash
# =============================================================================
# run_all_samples.sh
#
# Run the full variant calling pipeline on every BAM in bams/.
# Uses pipeline-stage parallelism: the next sample's callers (Steps 1-5)
# start as soon as the current sample's callers are complete, so VEP+
# annotation (Steps 6-10) of sample N overlaps with callers of sample N+1.
#
# Timeline:
#   Sample N:   [──callers──][──VEP+annotation──]
#   Sample N+1:             [──callers──][──VEP+annotation──]
#
# Resource profile (753 GB RAM, 56 cores):
#   Caller overlap peak: ~40-50 GB RAM — well within server capacity.
#   Max concurrent pipelines: controlled by MAX_CALLER_JOBS (default 1).
#   Increase to 2 to allow 2 samples in the caller phase simultaneously.
#
# Usage:
#   bash scripts/run_all_samples.sh [bams_dir]
#
# Default bams_dir: /data/salomonis-archive/FASTQs/NCI-R01/rna_seq_varcall/bams
# =============================================================================

set -uo pipefail

WORKDIR="/data/salomonis-archive/FASTQs/NCI-R01/rna_seq_varcall"
BAMS_DIR="${1:-$WORKDIR/bams}"
PIPELINE="$WORKDIR/scripts/run_sample_pipeline.sh"
MASTER_LOG="$WORKDIR/logs/run_all_samples.log"

# Maximum number of samples allowed in the caller phase (Steps 1-5) at once.
# 1 = conservative (safe on any server); 2 = allows full pipeline overlap.
MAX_CALLER_JOBS=1

# Poll interval (seconds) when waiting for caller steps to complete
POLL=60

mkdir -p "$WORKDIR/logs"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$MASTER_LOG"; }

# ---------------------------------------------------------------------------
# callers_done: returns 0 if all 4 caller VCFs exist for a sample
# ---------------------------------------------------------------------------
callers_done() {
    local sample="$1"
    local R="$WORKDIR/results/$sample"
    [ -f "$R/haplotypecaller/variants.vcf.gz" ] && \
    [ -f "$R/deepvariant/variants.vcf.gz"     ] && \
    [ -f "$R/clair3_rna/variants.vcf.gz"      ] && \
    [ -f "$R/longcallr/variants.vcf.gz"       ]
}

# ---------------------------------------------------------------------------
# pipeline_done: returns 0 if the final TSV exists for a sample
# ---------------------------------------------------------------------------
pipeline_done() {
    local sample="$1"
    [ -f "$WORKDIR/results/$sample/merged/high_confidence_final.tsv" ]
}

# ---------------------------------------------------------------------------
# already_running: returns 0 if a pipeline process for this BAM is alive
# ---------------------------------------------------------------------------
already_running() {
    local bam="$1"
    pgrep -f "run_sample_pipeline.sh.*$(basename "$bam")" > /dev/null 2>&1
}

log "=========================================================="
log "run_all_samples.sh started (pipeline-parallel mode)"
log "BAMs directory: $BAMS_DIR"
log "MAX_CALLER_JOBS: $MAX_CALLER_JOBS"
log "=========================================================="

BAM_LIST=("$BAMS_DIR"/*.bam)
TOTAL=${#BAM_LIST[@]}
log "Found $TOTAL BAMs:"
for bam in "${BAM_LIST[@]}"; do log "  $(basename "$bam")"; done
log ""

# Track background PIDs and sample names
declare -a BG_PIDS=()
declare -a BG_SAMPLES=()
declare -a CALLER_SAMPLES=()   # samples currently in caller phase
PASS=0
FAIL=0
declare -a FAILED_SAMPLES=()

for i in "${!BAM_LIST[@]}"; do
    bam="${BAM_LIST[$i]}"
    sample="$(basename "$bam" .bam)"
    n=$((i + 1))

    # ------------------------------------------------------------------
    # Skip if already fully complete (checkpoint)
    # ------------------------------------------------------------------
    if pipeline_done "$sample"; then
        log "[$n/$TOTAL] SKIP (complete): $sample"
        PASS=$((PASS + 1))
        continue
    fi

    # ------------------------------------------------------------------
    # If a pipeline is already running for this sample, attach to it:
    # wait until its callers are done, then move on to the next sample.
    # ------------------------------------------------------------------
    if already_running "$bam"; then
        log "[$n/$TOTAL] Already running: $sample — waiting for callers to finish"
        while ! callers_done "$sample" && ! pipeline_done "$sample"; do
            sleep "$POLL"
        done
        log "[$n/$TOTAL] Callers done (or pipeline complete) for $sample — continuing"
        # Register it as a "background" job we need to wait on at the end
        existing_pid=$(pgrep -f "run_sample_pipeline.sh.*$(basename "$bam")" | head -1)
        if [ -n "$existing_pid" ]; then
            BG_PIDS+=("$existing_pid")
            BG_SAMPLES+=("$sample")
        fi
        continue
    fi

    # ------------------------------------------------------------------
    # Wait until there is a free caller slot
    # ------------------------------------------------------------------
    while true; do
        # Remove finished callers from the active list
        active=()
        for s in "${CALLER_SAMPLES[@]+"${CALLER_SAMPLES[@]}"}"; do
            if callers_done "$s"; then
                log "  Callers done for $s — annotation continues in background"
            else
                active+=("$s")
            fi
        done
        CALLER_SAMPLES=("${active[@]+"${active[@]}"}")

        if [ "${#CALLER_SAMPLES[@]}" -lt "$MAX_CALLER_JOBS" ]; then
            break
        fi
        sleep "$POLL"
    done

    # ------------------------------------------------------------------
    # Start pipeline in background
    # ------------------------------------------------------------------
    log "----------------------------------------------------------"
    log "[$n/$TOTAL] Starting: $sample  ($(date))"
    log "----------------------------------------------------------"

    mkdir -p "$WORKDIR/logs/${sample}"
    bash "$PIPELINE" "$bam" \
        >> "$WORKDIR/logs/${sample}/pipeline.log" 2>&1 &
    pid=$!
    BG_PIDS+=("$pid")
    BG_SAMPLES+=("$sample")
    CALLER_SAMPLES+=("$sample")
    log "[$n/$TOTAL] Pipeline running in background (PID $pid)"
done

# ---------------------------------------------------------------------------
# Wait for all background pipelines to finish and collect results
# ---------------------------------------------------------------------------
log ""
log "All samples launched — waiting for remaining pipelines to finish..."
log ""

for j in "${!BG_PIDS[@]}"; do
    pid="${BG_PIDS[$j]}"
    sample="${BG_SAMPLES[$j]}"
    n=$((j + 1))

    if wait "$pid"; then
        log "[$n] SUCCESS: $sample  ($(date))"
        PASS=$((PASS + 1))
    else
        log "[$n] FAILED:  $sample  ($(date)) — see logs/${sample}/pipeline.log"
        FAIL=$((FAIL + 1))
        FAILED_SAMPLES+=("$sample")
    fi
done

log "=========================================================="
log "run_all_samples.sh finished: $(date)"
log "  Completed: $PASS / $TOTAL"
log "  Failed:    $FAIL / $TOTAL"
if [ "${#FAILED_SAMPLES[@]}" -gt 0 ]; then
    log "  Failed samples:"
    for s in "${FAILED_SAMPLES[@]}"; do log "    - $s"; done
fi
log "=========================================================="
