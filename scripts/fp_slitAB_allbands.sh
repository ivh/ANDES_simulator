#!/bin/bash
# FP and flat simulations with slitA and slitB for all bands,
# plus combined slitAB frames

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="$(dirname "$SCRIPT_DIR")"
OUTPUT_DIR="${SRC_DIR}/../fp_slitAB"

FLUX=0.01
JOBS=6

cd "$SRC_DIR"
mkdir -p "$OUTPUT_DIR"

run_sim() {
    local band=$1
    local source=$2
    local subslit=$3
    local tag=$(echo "$source" | tr '[:lower:]' '[:upper:]')
    uv run andes-sim simulate --band "$band" --source "$source" --subslit "$subslit" \
        --flux $FLUX \
        --output-dir "$OUTPUT_DIR" --output-name "${band}_${tag}_${subslit}.fits"
}
export -f run_sim
export FLUX OUTPUT_DIR

parallel -j$JOBS --bar --retries 2 --delay 0.5 run_sim {1} {2} {3} ::: \
    U B V R IZ Y J H ::: fp flat ::: slitA slitB

echo "Adding slitA + slitB pairs..."
for BAND in U B V R IZ Y J H; do
    for TAG in FP FLAT; do
        A="${OUTPUT_DIR}/${BAND}_${TAG}_slitA.fits"
        B="${OUTPUT_DIR}/${BAND}_${TAG}_slitB.fits"
        if [[ ! -f "$A" || ! -f "$B" ]]; then
            echo "  SKIP ${BAND}_${TAG} (missing slitA or slitB)"
            continue
        fi
        uv run python -c "
from astropy.io import fits
import numpy as np
a = fits.open('$A')
b = fits.open('$B')
a[0].data = np.array(a[0].data, dtype=np.float32) + np.array(b[0].data, dtype=np.float32)
a.writeto('${OUTPUT_DIR}/${BAND}_${TAG}_slitAB.fits', overwrite=True)
a.close(); b.close()
"
        echo "  ${BAND}_${TAG}_slitAB.fits"
    done
done

echo "Done. Output in $OUTPUT_DIR"
