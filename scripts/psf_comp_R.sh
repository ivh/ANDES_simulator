#!/bin/bash
# PSF comparison simulations for R-band
# Generates flat and LFC frames with even/odd fiber illumination
# for all available R-band HDF models

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="$(dirname "$SCRIPT_DIR")"
OUTPUT_DIR="${SRC_DIR}/../psf_comp_R"

WL_MIN=700
WL_MAX=715
FLUX=1
JOBS=6

cd "$SRC_DIR"
mkdir -p "$OUTPUT_DIR"

run_sim() {
    local tag=$1
    local hdf=$2
    local source=$3
    local subslit=$4
    local name="${tag}_${source}_${subslit}_wl${WL_MIN}-${WL_MAX}"
    uv run andes-sim simulate --hdf "$hdf" --source "$source" --subslit "$subslit" \
        --wl-min $WL_MIN --wl-max $WL_MAX --flux $FLUX \
        --output-dir "$OUTPUT_DIR" --output-name "${name}.fits"
}
export -f run_sim
export WL_MIN WL_MAX FLUX OUTPUT_DIR

parallel -j$JOBS --bar --retries 2 --delay 0.5 run_sim {1} {2} {3} {4} ::: \
    R0 R1 R2 :::+ \
    "HDF/ANDES_123_R3.hdf" \
    "HDF/ANDES_full_F18A33_win_jmr_MC_T0019_Rband_p0.hdf" \
    "HDF/Andes_full_F18A33_win_jmr_MC_T0108_Rband_P0_cfg1.hdf" \
    ::: flat lfc ::: even odd

echo "Done. Output in $OUTPUT_DIR"
