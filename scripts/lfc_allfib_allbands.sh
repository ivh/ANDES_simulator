#!/bin/bash
# LFC per-fiber simulations with velocity offsets for all bands

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="$(dirname "$SCRIPT_DIR")"
OUTPUT_DIR="${SRC_DIR}/../lfc_allfib_allbands"

FLUX=1
JOBS=6

cd "$SRC_DIR"
mkdir -p "$OUTPUT_DIR"

# Fiber counts: UBVRIZ=66, YJH=75
get_nfibers() {
    case $1 in
        Y|J|H) echo 75 ;;
        *)     echo 66 ;;
    esac
}

run_sim() {
    local band=$1
    local fiber=$2
    local src_dir=$3
    local out_dir=$4
    local flux=$5
    uv run andes-sim simulate --band "$band" --source lfc --fiber "$fiber" \
        --flux "$flux" \
        --velocity-shift "${src_dir}/data/vel_shifts_${band}.json" \
        --output-dir "$out_dir" \
        --output-name "${band}_LFC_fiber$(printf '%02d' "$fiber").fits"
}
export -f run_sim

# Build job list, skipping existing outputs
JOBFILE=$(mktemp)
SKIPPED=0
for BAND in U B V R IZ Y J H; do
    NFIB=$(get_nfibers "$BAND")
    for FIB in $(seq 1 "$NFIB"); do
        OUTFILE="${OUTPUT_DIR}/${BAND}_LFC_fiber$(printf '%02d' "$FIB").fits"
        if [[ -f "$OUTFILE" ]]; then
            SKIPPED=$((SKIPPED + 1))
        else
            echo "$BAND $FIB $SRC_DIR $OUTPUT_DIR $FLUX"
        fi
    done
done > "$JOBFILE"

TOTAL=$(wc -l < "$JOBFILE")
if [[ $TOTAL -eq 0 ]]; then
    echo "All files exist ($SKIPPED skipped). Nothing to do."
    rm -f "$JOBFILE"
    exit 0
fi
echo "Running $TOTAL LFC fiber simulations ($SKIPPED already exist)..."
parallel -j$JOBS --bar --retries 100 --delay 0.5 --colsep ' ' \
    run_sim {1} {2} {3} {4} {5} < "$JOBFILE"

rm -f "$JOBFILE"
echo "Done. Output in $OUTPUT_DIR"
