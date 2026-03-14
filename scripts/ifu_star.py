# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "astropy"]
# ///
"""
IFU star observation simulation for ANDES YJH bands.

Simulates a star centered on the hexagonal IFU with:
  - Star spectrum (HR1544 * sky transmission) on rings 0-4, brightness
    decreasing outward following the IFU spatial sampling.
  - Sky emission on all IFU fibers (uniform).
  - Fabry-Perot calibration on IFU cal fibers (1, 75).

All frames are combined into a single output per band.
Individual simulations run in parallel (6 workers).
"""

import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from astropy.io import fits

SRC_DIR = Path(__file__).resolve().parent.parent
SED_DIR = SRC_DIR / "SED"
OUT_DIR = SRC_DIR.parent / "ifu_star"

BANDS = ["Y", "J", "H"]
N_WORKERS = 6

# Ring flux weights: PSF profile across IFU (ring0 = center, brightest)
RING_FLUX = {
    "ring0": 10.0,
    "ring1": 0.5,
    "ring2": 0.2,
    "ring3": 0.08,
    "ring4": 0.03,
}

# Overall scaling
STAR_FLUX = 1.2
SKY_FLUX = 0.01
FP_FLUX = 80.0


def make_transmitted_star():
    """Multiply HR1544 by sky transmission, write to SED/."""
    star = np.loadtxt(SED_DIR / "HR1544.csv", delimiter=",", comments="#")
    trans = np.loadtxt(SED_DIR / "sky_transmission_YK.csv", delimiter=",")

    # Interpolate transmission onto star wavelength grid
    transmission = np.interp(star[:, 0], trans[:, 0], trans[:, 1],
                             left=0.0, right=0.0)
    transmitted = star.copy()
    transmitted[:, 1] *= transmission

    out = SED_DIR / "HR1544_transmitted.csv"
    with open(out, "w") as f:
        f.write("# scaling: 850\n")
        for wl, fl in transmitted:
            f.write(f"{wl},{fl}\n")
    print(f"Wrote {out}")
    return out


def run_sim(band, args, label):
    """Run a single andes-sim subprocess. Returns (label, returncode, stderr)."""
    cmd = ["uv", "run", "andes-sim", "simulate", "--band", band,
           "--output-dir", str(OUT_DIR)] + args
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(SRC_DIR))
    return label, result.returncode, result.stderr


def build_jobs(spectrum):
    """Build list of (band, args, label) for all simulations."""
    jobs = []
    for band in BANDS:
        for ring, weight in RING_FLUX.items():
            flux = STAR_FLUX * weight
            jobs.append((band, [
                "--source", str(spectrum),
                "--fiber", ring,
                "--flux", f"{flux}",
                "--output-name", f"{band}_star_{ring}.fits",
            ], f"{band} star {ring}"))

        jobs.append((band, [
            "--source", str(SED_DIR / "sky_emission_YK.csv"),
            "--fiber", "ifu",
            "--flux", f"{SKY_FLUX}",
            "--output-name", f"{band}_sky_ifu.fits",
        ], f"{band} sky emission"))

        jobs.append((band, [
            "--source", "fp",
            "--fiber", "cal_ifu",
            "--flux", f"{FP_FLUX}",
            "--output-name", f"{band}_fp_cal.fits",
        ], f"{band} FP cal"))

    return jobs


def combine_frames(band):
    """Add all component frames into a single FITS image."""
    parts = (
        [f"{band}_star_{ring}.fits" for ring in RING_FLUX]
        + [f"{band}_sky_ifu.fits", f"{band}_fp_cal.fits"]
    )
    combined = None
    header = None
    for name in parts:
        path = OUT_DIR / name
        if not path.exists():
            print(f"  WARNING: {path} missing, skipping")
            continue
        with fits.open(path) as hdul:
            data = hdul[0].data.astype(np.float64)
            if combined is None:
                combined = data
                header = hdul[0].header.copy()
            else:
                combined += data

    out = OUT_DIR / f"{band}_ifu_star_combined.fits"
    fits.writeto(out, combined, header=header, overwrite=True)
    print(f"  Combined -> {out}")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Creating transmitted star spectrum")
    spectrum = make_transmitted_star()

    jobs = build_jobs(spectrum)
    n_total = len(jobs)
    n_done = 0
    failed = []

    print(f"Running {n_total} simulations with {N_WORKERS} workers")
    with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
        futures = {pool.submit(run_sim, *job): job[2] for job in jobs}
        for future in as_completed(futures):
            label, rc, stderr = future.result()
            n_done += 1
            if rc != 0:
                err = stderr.splitlines()[-1] if stderr else "unknown error"
                print(f"  [{n_done}/{n_total}] FAILED {label}: {err}")
                failed.append(label)
            else:
                print(f"  [{n_done}/{n_total}] {label}")

    if failed:
        print(f"\n{len(failed)} simulations failed: {failed}", file=sys.stderr)
        sys.exit(1)

    print("Combining frames")
    for band in BANDS:
        combine_frames(band)

    print("Done.")


if __name__ == "__main__":
    main()
