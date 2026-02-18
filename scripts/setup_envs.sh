#!/bin/bash
#
# Create isolated conda environments and pull containers for each variant caller.
# Run this once before executing run_callers.sh.
#
# Usage: bash scripts/setup_envs.sh [all|bio-cli|gatk|clair3|longcallr]
#

set -e

ENV="${1:-all}"

WORKDIR="/data/salomonis-archive/FASTQs/NCI-R01/rna_seq_varcall"
MODEL_DIR="$WORKDIR/models/clair3_rna"

# Detect conda binary and envs root (works on Linux server and macOS)
_CONDA_BIN="${CONDA_EXE:-$(command -v conda 2>/dev/null)}"
CONDA_ENVS_DIR="$("$_CONDA_BIN" info --base 2>/dev/null)/envs"
echo "[INFO] conda envs dir: $CONDA_ENVS_DIR"

# Use mamba if available for faster solves
_INSTALL_CMD="conda"
if command -v mamba &>/dev/null; then
    _INSTALL_CMD="mamba"
fi

# ── bio-cli ───────────────────────────────────────────────────────────────
setup_bio_cli() {
    echo "[INFO] Creating bio-cli environment (samtools, bcftools, bgzip)..."
    "$_INSTALL_CMD" create -n bio-cli -c bioconda -c conda-forge \
        samtools bcftools htslib -y

    echo "[INFO] Verifying..."
    "$CONDA_ENVS_DIR/bio-cli/bin/samtools" --version 2>/dev/null | head -1
    "$CONDA_ENVS_DIR/bio-cli/bin/bcftools" --version 2>/dev/null | head -1
    echo "[SUCCESS] bio-cli ready"
}

# ── gatk-env ──────────────────────────────────────────────────────────────
setup_gatk() {
    echo "[INFO] Creating gatk-env (OpenJDK 17 + GATK 4.6.2.0)..."
    "$_INSTALL_CMD" create -n gatk-env -c conda-forge -c bioconda \
        openjdk=17 gatk4=4.6.2.0 -y

    echo "[INFO] Verifying GATK..."
    # GATK wrapper needs its env bin on PATH to find 'python'
    PATH="$CONDA_ENVS_DIR/gatk-env/bin:$PATH" \
        "$CONDA_ENVS_DIR/gatk-env/bin/gatk" --version
    echo "[SUCCESS] gatk-env ready"
}

# ── clair3 ────────────────────────────────────────────────────────────────
# Clair3 conda package requires Linux x86_64; use Docker (macOS) or Singularity (Linux)
setup_clair3() {
    CLAIR3_IMAGE="hkubal/clair3:latest"
    CLAIR3_SIF="$WORKDIR/containers/clair3_latest.sif"
    mkdir -p "$WORKDIR/containers" 2>/dev/null || true

    if command -v singularity &>/dev/null; then
        echo "[INFO] Pulling Clair3 Singularity SIF..."
        singularity pull --dir "$WORKDIR/containers" "docker://$CLAIR3_IMAGE"
        echo "[SUCCESS] Clair3 SIF ready at $CLAIR3_SIF"
    elif command -v docker &>/dev/null; then
        echo "[INFO] Pulling Clair3 Docker image: $CLAIR3_IMAGE"
        docker pull "$CLAIR3_IMAGE"
        echo "[SUCCESS] Clair3 Docker image ready"
    else
        echo "[WARN] Neither Singularity nor Docker available."
        echo "[INFO] macOS: start Docker Desktop, then re-run."
        echo "[INFO] Linux: module load singularity or conda install singularity"
    fi

    if [ -d "$MODEL_DIR" ] && [ "$(ls -A "$MODEL_DIR" 2>/dev/null)" ]; then
        echo "[INFO] Clair3 models already present at $MODEL_DIR, skipping download."
    else
        echo "[INFO] Downloading Clair3 PacBio HiFi models..."
        mkdir -p "$MODEL_DIR"
        MODEL_URL="http://www.bio8.cs.hku.hk/clair3_models/clair3_models.tar.gz"
        wget -q --show-progress "$MODEL_URL" -O /tmp/clair3_models.tar.gz
        tar -xzf /tmp/clair3_models.tar.gz -C "$MODEL_DIR" --strip-components=1
        rm /tmp/clair3_models.tar.gz
        echo "[SUCCESS] Clair3 models downloaded"
        echo "[INFO] Available models:"
        ls "$MODEL_DIR"
        echo "[INFO] Recommended: hifi_revio or hifi_sequel2"
        echo "[INFO] Set CLAIR3_MODEL in run_callers.sh to the chosen subdirectory."
    fi
}

# ── longcallr ─────────────────────────────────────────────────────────────
# longcallr (C++ BAM→VCF) is available as a BioContainers image
setup_longcallr() {
    LONGCALLR_IMAGE="quay.io/biocontainers/longcallr:1.12.0--py312h67e1f27_0"
    LONGCALLR_SIF="$WORKDIR/containers/longcallr_1.12.0.sif"
    mkdir -p "$WORKDIR/containers" 2>/dev/null || true

    if command -v singularity &>/dev/null; then
        echo "[INFO] Pulling LongcallR Singularity SIF..."
        singularity pull --dir "$WORKDIR/containers" "docker://$LONGCALLR_IMAGE"
        echo "[SUCCESS] LongcallR SIF ready at $LONGCALLR_SIF"
    elif command -v docker &>/dev/null; then
        echo "[INFO] Pulling LongcallR Docker image: $LONGCALLR_IMAGE"
        docker pull "$LONGCALLR_IMAGE"
        echo "[SUCCESS] LongcallR Docker image ready"
    else
        echo "[WARN] Neither Singularity nor Docker available."
        echo "[INFO] macOS: start Docker Desktop, then re-run."
        echo "[INFO] Linux: singularity pull docker://$LONGCALLR_IMAGE"
    fi
}

# ── main ──────────────────────────────────────────────────────────────────
case "$ENV" in
    bio-cli)
        setup_bio_cli
        ;;
    gatk)
        setup_gatk
        ;;
    clair3)
        setup_clair3
        ;;
    longcallr)
        setup_longcallr
        ;;
    all)
        setup_bio_cli
        setup_gatk
        setup_clair3
        setup_longcallr
        ;;
    *)
        echo "[ERROR] Unknown target: $ENV"
        echo "Usage: $0 [all|bio-cli|gatk|clair3|longcallr]"
        exit 1
        ;;
esac

echo ""
echo "========================================================"
echo "[DONE] Setup complete. Verify with:"
echo "  $CONDA_ENVS_DIR/bio-cli/bin/samtools --version"
echo "  PATH=$CONDA_ENVS_DIR/gatk-env/bin:\$PATH $CONDA_ENVS_DIR/gatk-env/bin/gatk --version"
echo "  docker run --rm hkubal/clair3:latest /usr/local/bin/run_clair3.sh --version"
echo "  singularity exec $WORKDIR/containers/longcallr_1.12.0.sif longcallr --version"
echo "========================================================"
