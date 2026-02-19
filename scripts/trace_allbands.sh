#!/bin/bash
# Flat simulations with even and odd subslits for all bands

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="$(dirname "$SCRIPT_DIR")"
OUTPUT_DIR="${SRC_DIR}/../trace_allbands"

FLUX=0.01
JOBS=6

cd "$SRC_DIR"
mkdir -p "$OUTPUT_DIR"

run_sim() {
    local band=$1
    local subslit=$2
    uv run andes-sim simulate --band "$band" --source flat --subslit "$subslit" \
        --flux $FLUX \
        --output-dir "$OUTPUT_DIR" --output-name "${band}_FLAT_${subslit}.fits"
}
export -f run_sim
export FLUX OUTPUT_DIR

parallel -j$JOBS --bar --retries 2 --delay 0.5 run_sim {1} {2} ::: \
    U B V R IZ Y J H ::: even odd

echo "Done. Output in $OUTPUT_DIR"
