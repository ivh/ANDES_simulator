# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "astropy>=7.0", "skycalc_ipy"]
# ///
"""
LR-J band star+sky observation simulation for MOSAIC.

Simulates a stellar observation on a single MOSAIC NIR fiber bundle in the
middle of the detector:
  - Star spectrum (attenuated by sky transmission) + sky emission on the
    chosen bundle.
  - All other bundles are dark.

Sky transmission and emission are fetched from ESO SkyCalc and cached in
SED/sky_transmission_YK.csv / SED/sky_emission_YK.csv (LR-J 950-1340 nm
falls inside the cached YK range, so the cache is reused if present).

Star and sky frames are added into a single output. Both component
simulations run in parallel.
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
OUT_DIR = SRC_DIR.parent / "mosaic_LR-J_starsky"

BAND = "LR-J"
N_WORKERS = 2
FIB_EFF = "0.75-0.95"

# MOSAIC NIR has 90 bundles x 7 fibers = 630 fibers.
# Bundle 45 sits roughly in the middle of the spatial axis.
BUNDLE = 45

STAR_SPECTRUM = SED_DIR / "phoenix.csv"

# MOSAIC NIR has 630 fibers (vs 66 for ANDES R), so PyEchelle iterates
# many more fibers; keep the per-fiber photon count modest so a full
# bundle:45 raytrace finishes in roughly a minute. Scale up later if
# higher-SNR frames are needed.
STAR_FLUX = 0.5
SKY_FLUX = 0.005


def make_sky_files():
    """Fetch LR-J sky transmission and emission from ESO SkyCalc.

    Reuses the YK cache if present (it spans 950-2400 nm which covers LR-J).
    Returns (transmission_path, emission_path).
    """
    trans_path = SED_DIR / "sky_transmission_YK.csv"
    emiss_path = SED_DIR / "sky_emission_YK.csv"

    if trans_path.exists() and emiss_path.exists():
        print(f"Using cached sky files: {trans_path.name}, {emiss_path.name}")
        return trans_path, emiss_path

    import skycalc_ipy

    print("Fetching YK sky model from ESO SkyCalc...")
    sc = skycalc_ipy.SkyCalc()
    sc["wmin"] = 950
    sc["wmax"] = 2400
    sc["wgrid_mode"] = "fixed_wavelength_step"
    sc["wdelta"] = 0.05
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

    out = SED_DIR / "star_transmitted_LR-J.csv"
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
    """Run a single mosaic-sim subprocess with isolated numba cache.

    Per-job stdout/stderr are captured into OUT_DIR/<label>.log so PyEchelle's
    per-fiber photon prints are visible while the job runs, instead of
    buffered in a pipe until the process exits.
    """
    cache_dir = tempfile.mkdtemp(prefix="numba_")
    env = {**os.environ, "NUMBA_CACHE_DIR": cache_dir}
    cmd = ["uv", "run", "mosaic-sim", "simulate", "--band", BAND,
           "--output-dir", str(OUT_DIR), "--fib-eff", FIB_EFF] + args
    log_path = OUT_DIR / f"{label.replace(' ', '_')}.log"
    try:
        with open(log_path, "w") as logf:
            result = subprocess.run(cmd, stdout=logf, stderr=subprocess.STDOUT,
                                    text=True, cwd=str(SRC_DIR), env=env)
        return label, result.returncode, log_path
    finally:
        import shutil
        shutil.rmtree(cache_dir, ignore_errors=True)


def build_jobs(transmitted_star, sky_emission):
    """Build star and sky simulation jobs for the chosen bundle."""
    bundle_arg = f"bundle:{BUNDLE}"
    jobs = [
        ([
            "--source", str(transmitted_star),
            "--fiber", bundle_arg,
            "--flux", f"{STAR_FLUX}",
            "--output-name", f"{BAND}_star_bundle{BUNDLE}.fits",
        ], f"star bundle{BUNDLE}"),
        ([
            "--source", str(sky_emission),
            "--fiber", bundle_arg,
            "--flux", f"{SKY_FLUX}",
            "--output-name", f"{BAND}_sky_bundle{BUNDLE}.fits",
        ], f"sky bundle{BUNDLE}"),
    ]
    return jobs


FRAME_NAMES = [
    f"{BAND}_star_bundle{BUNDLE}.fits",
    f"{BAND}_sky_bundle{BUNDLE}.fits",
]


def combine_frames():
    """Add the star and sky frames into a single combined image."""
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
            label, rc, log_path = future.result()
            n_done += 1
            if rc != 0:
                print(f"  [{n_done}/{n_total}] FAILED {label} (see {log_path})")
                failed.append(label)
            else:
                print(f"  [{n_done}/{n_total}] {label} (log: {log_path})")

    if failed:
        print(f"\n{len(failed)} simulations failed: {failed}", file=sys.stderr)
        sys.exit(1)

    print("Combining frames")
    combine_frames()

    print("Done.")


if __name__ == "__main__":
    main()
