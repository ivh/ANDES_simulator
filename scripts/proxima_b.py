# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "astropy", "scipy"]
# ///
"""Proxima b reflected-light observing sequence on the ANDES SCAO-IFU.

Reproduces the Bugatti et al. 2025 / ANDES IFS WG acquisition strategy:

    SHORT_A   star scanned over every IFU spaxel with the tip-tilt mirror,
              giving a high-S/N star-only spectrum in each spaxel
    LONG      star centred on the central spaxel behind a neutral density
              filter; the AO residual halo falls into all 61 spaxels; the
              planet sits in one of them with its own Doppler shift
    SHORT_B   repeat of SHORT_A

The two short exposures bracket the long one and serve as the halo+telluric
template.  See ../../../exopl_refl/PLAN.md.

=============================================================================
THE HALO IS NOT PHYSICS YET -- READ THIS
=============================================================================
The E2E simulator has NO telescope-side PSF model.  Grepping andes_simulator
for strehl/halo/coupling/seeing/coronagraph/adaptive returns nothing; the only
"Airy" in the package is the Fabry-Perot transmission function, and every
"PSF" is the ZEMAX *spectrograph* PSF at the detector, which is a different
thing entirely.  scripts/ifu_star.py fakes it with five hand-written ring
weights (10.0/0.5/0.2/0.08/0.03) and no physical basis.

So the focal-plane PSF -- the thing that decides how much starlight couples
into each spaxel, i.e. the entire contrast problem -- has to be supplied here.
The focal-plane coupling therefore comes from `andes_coupling.py`, which is
anchored to the three numbers the PDR package does give: Strehl 0.6/0.3 in H/Y,
raw contrast ~1.5e-3 at 25-45 mas without a coronagraph, and the R-AND-102.0a
requirement of 3.0e-3 at 20 mas with one.  It is still an analytic model -- the
AO subsystem documents holding the actual contrast curves (E-AND-AO-TNO-04-00-001,
E-AND-AO-ANR-14-00-003) are referenced but not delivered in the PDR package -- but
it is no longer arbitrary.  See exopl_refl/phase1/RESULTS.md.
=============================================================================
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).resolve().parent))
from andes_coupling import Coupling, T_SCAO, STREHL

SRC = Path(__file__).resolve().parent.parent
SED = SRC / "SED"
OUT_DEFAULT = SRC.parent / "proxima_b"

# --- IFU: 61 hexagonal spaxels in 5 rings, mapped to pseudo-slit fibres ------
# from andes_simulator/core/andes.py : FIBER_CONFIGS['YJH_IFU']
RING_FIBERS = {
    0: [3],
    1: [6, 8, 10, 12, 14, 16],
    2: list(range(18, 30)),
    3: list(range(31, 49)),
    4: list(range(50, 74)),
}
CAL_FIBERS = [1, 75]

LAMBDA_D_MAS = {"Y": 5.47, "J": 6.64, "H": 8.65}   # ELT, D=39 m, band centre
SEP_PROXB_MAS = 37.3                                # max elongation


def hex_ring(n: int):
    """Axial coordinates of ring n of a hexagonal lattice (6n cells)."""
    if n == 0:
        return [(0, 0)]
    dirs = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
    q, r, out = n, 0, []
    for d in range(6):
        for _ in range(n):
            out.append((q, r))
            q += dirs[(d + 2) % 6][0]
            r += dirs[(d + 2) % 6][1]
    return out


def build_ifu(pitch_mas: float):
    """-> list of (fiber, x_mas, y_mas, ring), one per spaxel."""
    s = pitch_mas / np.sqrt(3)          # circumradius
    spaxels = []
    for ring, fibers in RING_FIBERS.items():
        cells = hex_ring(ring)
        if len(cells) != len(fibers):
            raise ValueError(f"ring {ring}: {len(cells)} lattice cells vs {len(fibers)} fibres")
        for (q, r), fib in zip(cells, fibers):
            x = s * np.sqrt(3) * (q + r / 2.0)
            y = s * 1.5 * r
            spaxels.append((fib, x, y, ring))
    return spaxels


def run_sim(band, args, label, out_dir, fib_eff):
    """One andes-sim subprocess with an isolated numba cache."""
    cache = tempfile.mkdtemp(prefix="numba_")
    env = {**os.environ, "NUMBA_CACHE_DIR": cache}
    cmd = ["uv", "run", "andes-sim", "simulate", "--band", band,
           "--output-dir", str(out_dir), "--fib-eff", fib_eff] + args
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(SRC), env=env)
        return label, r.returncode, r.stderr[-500:] if r.returncode else ""
    finally:
        shutil.rmtree(cache, ignore_errors=True)


def combine(out_dir: Path, pattern: str, output: Path):
    """Sum every frame matching `pattern` into one detector image."""
    files = sorted(out_dir.glob(pattern))
    if not files:
        raise FileNotFoundError(pattern)
    acc, hdr = None, None
    for fn in files:
        d = fits.getdata(fn).astype(np.float64)
        acc = d if acc is None else acc + d
        if hdr is None:
            hdr = fits.getheader(fn)
    hdr["NCOMB"] = (len(files), "frames summed")
    fits.writeto(output, acc.astype(np.float32), hdr, overwrite=True)
    return output, len(files), float(acc.sum())


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--band", default="J", choices=["Y", "J", "H"])
    p.add_argument("--ifu-scale", type=float, default=16.0, help="spaxel pitch [mas]")
    p.add_argument("--exptime-long", type=float, default=3600.0)
    p.add_argument("--exptime-short", type=float, default=120.0,
                   help="TOTAL short-exposure time; the scan dwells t/61 per spaxel")
    p.add_argument("--coronagraph", action="store_true",
                   help="use the coronagraph (R-AND-102.0a contrast, 0.54 transmission)")
    p.add_argument("--nd", type=float, default=1e-3,
                   help="neutral-density transmission on the central spaxel")
    p.add_argument("--planet-pa", type=float, default=35.0, help="planet position angle [deg]")
    p.add_argument("--albedo", type=float, default=0.3)
    p.add_argument("--planet-rv", type=float, default=30000.0, help="planet RV [m/s]")
    p.add_argument("--flux-scale", type=float, default=1.0,
                   help="global multiplier on all fluxes (see CALIBRATION note)")
    p.add_argument("--jobs", type=int, default=6)
    p.add_argument("--out", type=Path, default=OUT_DEFAULT)
    p.add_argument("--only", choices=["shortA", "long", "shortB"], action="append")
    p.add_argument("--dry-run", action="store_true")
    a = p.parse_args()

    spectrum = SED / "proxima.csv"
    if not spectrum.exists():
        sys.exit(f"{spectrum} missing -- run scripts/make_proxima_sed.py first")

    lam_d = LAMBDA_D_MAS[a.band]
    spaxels = build_ifu(a.ifu_scale)
    cm = Coupling(a.band, coronagraph=a.coronagraph)

    # ---- geometry -----------------------------------------------------------
    pa = np.radians(a.planet_pa)
    px, py = SEP_PROXB_MAS * np.cos(pa), SEP_PROXB_MAS * np.sin(pa)
    # planet's spaxel = nearest centre
    pfib, pd = None, 1e9
    for fib, x, y, ring in spaxels:
        dd = np.hypot(x - px, y - py)
        if dd < pd:
            pd, pfib, pring = dd, fib, ring

    rho_star = cm.rho(spaxels, 0.0, 0.0, a.ifu_scale)   # star on centre -> halo everywhere
    rho_planet = cm.rho(spaxels, px, py, a.ifu_scale)   # planet -> mostly its own spaxel

    contrast = a.albedo * (1.07 * 6.371e6 / (0.04848 * 1.496e11)) ** 2 / np.pi

    # ---- report -------------------------------------------------------------
    print(f"band {a.band}   lam/D = {lam_d:.2f} mas   pitch = {a.ifu_scale:.0f} mas")
    print(f"planet at PA {a.planet_pa:.0f} deg, {SEP_PROXB_MAS:.1f} mas "
          f"= {SEP_PROXB_MAS/lam_d:.2f} lam/D -> fibre {pfib} (ring {pring}), "
          f"{pd:.1f} mas off centre")
    print(f"planet/star contrast = {contrast:.3e}")
    print(f"Strehl {cm.S:.2f}, halo anchored to C={cm.contrast:.1e} @ "
          f"{cm.r_anchor:.0f} mas -> control radius {cm.r_c:.0f} mas; "
          f"T_SCAO {T_SCAO['coro' if a.coronagraph else 'nocoro']:.2f}")
    print(f"coupling: central spaxel {rho_star[RING_FIBERS[0][0]]:.4f}, "
          f"planet's spaxel (halo) {rho_star[pfib]:.3e}, "
          f"planet into it {rho_planet[pfib]:.4f}")
    print(f"planet/halo in fibre {pfib} = "
          f"{contrast*rho_planet[pfib]/rho_star[pfib]:.3e}")
    print(f"slit separation planet<->centre = {abs(pfib - RING_FIBERS[0][0])} fibres "
          f"(>=3 needed, see exopl_refl/phase5/RESULTS.md)")

    # ---- build the job list -------------------------------------------------
    want = set(a.only) if a.only else {"shortA", "long", "shortB"}
    jobs, dwell = [], a.exptime_short / len(spaxels)

    for tag in ("shortA", "shortB"):
        if tag not in want:
            continue
        for fib, x, y, ring in spaxels:
            # tip-tilt scan: star centred on THIS spaxel for `dwell` seconds
            rho = cm.rho([(fib, x, y, ring)], x, y, a.ifu_scale)[fib]
            jobs.append((a.band, [
                "--source", str(spectrum), "--fiber", str(fib),
                "--exposure", f"{dwell:.4f}",
                "--flux", f"{rho*a.flux_scale:.6e}",
                "--output-name", f"{a.band}_{tag}_f{fib:02d}.fits",
            ], f"{tag} fibre {fib}"))

    if "long" in want:
        for fib, x, y, ring in spaxels:
            nd = a.nd if ring == 0 else 1.0
            jobs.append((a.band, [
                "--source", str(spectrum), "--fiber", str(fib),
                "--exposure", f"{a.exptime_long:.1f}",
                "--flux", f"{rho_star[fib]*nd*a.flux_scale:.6e}",
                "--output-name", f"{a.band}_long_halo_f{fib:02d}.fits",
            ], f"long halo fibre {fib}"))
        # the planet: its own Doppler shift, own coupling
        jobs.append((a.band, [
            "--source", str(spectrum), "--fiber", str(pfib),
            "--exposure", f"{a.exptime_long:.1f}",
            "--flux", f"{contrast*rho_planet[pfib]*a.flux_scale:.6e}",
            "--velocity-shift", f"{a.planet_rv:.1f}",
            "--output-name", f"{a.band}_long_planet_f{pfib:02d}.fits",
        ], f"long PLANET fibre {pfib}"))
        # sky emission on every IFU fibre
        sky = SED / "sky_emission_YK.csv"
        if sky.exists():
            jobs.append((a.band, [
                "--source", str(sky), "--fiber", "ifu",
                "--exposure", f"{a.exptime_long:.1f}", "--flux", "0.01",
                "--output-name", f"{a.band}_long_sky.fits",
            ], "long sky"))
        # simultaneous FP on the IFU calibration fibres
        jobs.append((a.band, [
            "--source", "fp", "--fiber", "cal_ifu",
            "--exposure", f"{a.exptime_long:.1f}", "--flux", "80",
            "--output-name", f"{a.band}_long_fp.fits",
        ], "long FP cal"))

    print(f"\n{len(jobs)} simulations -> {a.out}")
    if a.dry_run:
        for _, args_, lab in jobs[:8]:
            print("  ", lab, " ".join(args_[:8]))
        print(f"   ... ({len(jobs)} total)")
        return

    a.out.mkdir(parents=True, exist_ok=True)
    fib_eff = "1.0"
    done = fail = 0
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        futs = {ex.submit(run_sim, b, ar, lb, a.out, fib_eff): lb for b, ar, lb in jobs}
        for fut in as_completed(futs):
            lab, rc, err = fut.result()
            done += 1
            if rc:
                fail += 1
                print(f"  FAILED {lab}: {err[:200]}")
            if done % 10 == 0:
                print(f"  {done}/{len(jobs)}")
    print(f"{done - fail}/{len(jobs)} succeeded")

    for tag, pat in (("shortA", f"{a.band}_shortA_f*.fits"),
                     ("long", f"{a.band}_long_*.fits"),
                     ("shortB", f"{a.band}_shortB_f*.fits")):
        if tag not in want:
            continue
        try:
            out, n, tot = combine(a.out, pat, a.out / f"{a.band}_proxb_{tag}.fits")
            print(f"  {out.name}: {n} frames, {tot:.4e} e-")
        except FileNotFoundError:
            print(f"  {tag}: nothing to combine")


if __name__ == "__main__":
    main()
