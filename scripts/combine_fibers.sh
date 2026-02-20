#!/bin/bash
# Combine per-fiber LFC outputs: all fibers and even/odd splits for each band

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="$(dirname "$SCRIPT_DIR")"
DATA_DIR="${SRC_DIR}/../lfc_allfib_allbands"

cd "$SRC_DIR"

for BAND in U B V R IZ Y J H; do
    echo "=== $BAND ==="

    uv run andes-sim combine --band "$BAND" \
        --input-pattern "${BAND}_LFC_fiber{fib:02d}.fits" \
        --mode all \
        --output-dir "$DATA_DIR" \
        --output "${BAND}_LFC_combined_all.fits"

    uv run andes-sim combine --band "$BAND" \
        --input-pattern "${BAND}_LFC_fiber{fib:02d}.fits" \
        --mode even_odd \
        --output-dir "$DATA_DIR"
done

echo "Done. Output in $DATA_DIR"
