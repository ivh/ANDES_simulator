# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "astropy>=7.0", "skycalc_ipy"]
# ///
"""
R-band star+sky observation simulation for ANDES.

Simulates a stellar observation with:
  - Star spectrum (attenuated by sky transmission) + sky emission on slit A.
  - Sky emission only on slit B (for sky subtraction).
  - Fabry-Perot calibration on cal fibers (33, 34).
  - Per-fiber velocity shifts from data/vel_shifts_R.json.

Sky transmission and emission are fetched from ESO SkyCalc and cached
in SED/sky_transmission_R.csv and SED/sky_emission_R.csv.

All component frames are combined into a single output.
Individual simulations run in parallel.
"""

import os
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from astropy.io import fits

SRC_DIR = Path(__file__).resolve().parent.parent
SED_DIR = SRC_DIR / "SED"
DATA_DIR = SRC_DIR / "data"
OUT_DIR = SRC_DIR.parent / "R_starsky"

BAND = "R"
N_WORKERS = 6
FIB_EFF = "0.75-0.95"

# --- Input spectrum (must cover R-band: ~621-764 nm) ---
STAR_SPECTRUM = SED_DIR / "phoenix.csv"

# --- Flux scaling ---
STAR_FLUX = 10.0
SKY_FLUX = 1.0
FP_FLUX = 10.0


def make_sky_files():
    """Fetch R-band sky transmission and emission from ESO SkyCalc.

    Results are cached; delete the CSV files to regenerate.
    Returns (transmission_path, emission_path).
    """
    trans_path = SED_DIR / "sky_transmission_R.csv"
    emiss_path = SED_DIR / "sky_emission_R.csv"

    if trans_path.exists() and emiss_path.exists():
        print(f"Using cached sky files: {trans_path.name}, {emiss_path.name}")
        return trans_path, emiss_path

    import skycalc_ipy

    print("Fetching R-band sky model from ESO SkyCalc...")
    sc = skycalc_ipy.SkyCalc()
    sc["wmin"] = 600
    sc["wmax"] = 780
    sc["wgrid_mode"] = "fixed_wavelength_step"
    sc["wdelta"] = 0.01
    tbl = sc.get_sky_spectrum()

    with open(trans_path, "w") as f:
        for row in tbl:
            f.write(f"{row['lam']:.4f},{row['trans']:.6e}\n")
    print(f"Wrote {trans_path}")

    with open(emiss_path, "w") as f:
        f.write("# scaling: 5\n")
        for row in tbl:
            f.write(f"{row['lam']:.4f},{row['flux']:.6e}\n")
    print(f"Wrote {emiss_path}")

    return trans_path, emiss_path


def make_transmitted_star(sky_trans_path):
    """Multiply star spectrum by sky transmission, write to SED/."""
    star = np.loadtxt(STAR_SPECTRUM, delimiter=",", comments="#")
    trans = np.loadtxt(sky_trans_path, delimiter=",", comments="#")

    transmission = np.interp(star[:, 0], trans[:, 0], trans[:, 1],
                             left=0.0, right=0.0)
    transmitted = star.copy()
    transmitted[:, 1] *= transmission

    out = SED_DIR / "star_transmitted_R.csv"
    with open(out, "w") as f:
        with open(STAR_SPECTRUM) as orig:
            for line in orig:
                if line.startswith("# scaling:"):
                    f.write(line)
                    break
        for wl, fl in transmitted:
            f.write(f"{wl},{fl}\n")
    print(f"Wrote {out}")
    return out


def run_sim(args, label):
    """Run a single andes-sim subprocess with isolated numba cache."""
    cache_dir = tempfile.mkdtemp(prefix="numba_")
    env = {**os.environ, "NUMBA_CACHE_DIR": cache_dir}
    cmd = ["uv", "run", "andes-sim", "simulate", "--band", BAND,
           "--output-dir", str(OUT_DIR), "--fib-eff", FIB_EFF] + args
    try:
        result = subprocess.run(cmd, capture_output=True, text=True,
                                cwd=str(SRC_DIR), env=env)
        return label, result.returncode, result.stderr
    finally:
        import shutil
        shutil.rmtree(cache_dir, ignore_errors=True)


def build_jobs(transmitted_star, sky_emission):
    """Build simulation jobs for star+sky scene."""
    vshift_file = DATA_DIR / f"vel_shifts_{BAND}.json"
    vshift_arg = (["--velocity-shift", str(vshift_file)]
                  if vshift_file.exists() else [])

    jobs = []

    # Slit A: transmitted star
    jobs.append(([
        "--source", str(transmitted_star),
        "--subslit", "slitA",
        "--flux", f"{STAR_FLUX}",
        "--output-name", f"{BAND}_star_slitA.fits",
    ] + vshift_arg, "star slitA"))

    # Slit A: sky emission
    jobs.append(([
        "--source", str(sky_emission),
        "--subslit", "slitA",
        "--flux", f"{SKY_FLUX}",
        "--output-name", f"{BAND}_sky_slitA.fits",
    ] + vshift_arg, "sky slitA"))

    # Slit B: sky emission only
    jobs.append(([
        "--source", str(sky_emission),
        "--subslit", "slitB",
        "--flux", f"{SKY_FLUX}",
        "--output-name", f"{BAND}_sky_slitB.fits",
    ] + vshift_arg, "sky slitB"))

    # Cal fibers: Fabry-Perot
    jobs.append(([
        "--source", "fp",
        "--subslit", "cal_sl",
        "--flux", f"{FP_FLUX}",
        "--output-name", f"{BAND}_fp_cal.fits",
    ] + vshift_arg, "FP cal"))

    return jobs


FRAME_NAMES = [
    f"{BAND}_star_slitA.fits",
    f"{BAND}_sky_slitA.fits",
    f"{BAND}_sky_slitB.fits",
    f"{BAND}_fp_cal.fits",
]


def combine_frames():
    """Add all component frames into a single FITS image."""
    combined = None
    header = None
    for name in FRAME_NAMES:
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

    out = OUT_DIR / f"{BAND}_starsky_combined.fits"
    fits.writeto(out, combined, header=header, overwrite=True)
    print(f"  Combined -> {out}")


def main():
    if not STAR_SPECTRUM.exists():
        print(f"ERROR: star spectrum not found: {STAR_SPECTRUM}", file=sys.stderr)
        sys.exit(1)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Preparing sky model")
    sky_trans, sky_emiss = make_sky_files()

    print("Creating transmitted star spectrum")
    transmitted = make_transmitted_star(sky_trans)

    jobs = build_jobs(transmitted, sky_emiss)
    n_total = len(jobs)
    n_done = 0
    failed = []

    print(f"Running {n_total} simulations with {N_WORKERS} workers")
    with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
        futures = {pool.submit(run_sim, *job): job[1] for job in jobs}
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
    combine_frames()

    print("Done.")


if __name__ == "__main__":
    main()
