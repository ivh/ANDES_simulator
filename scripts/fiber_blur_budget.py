"""How much extra image blur the R requirement allows, and what it does to FWHM/pitch.

The HDF models contain the ZEMAX design PSFs but no as-built degradation
(tolerances, alignment, thermal, charge diffusion/IPC). The TRS (ESO-391757)
has no cross-dispersion or crosstalk requirement; the only relevant one is
R-AND-88 (average R=100k, variation within -10%/+30%). This script adds an
isotropic Gaussian blur sigma (0..2 px, 0.05 px steps) to the simulated fiber
image and records, per band/order/column, both the resolving power and the
cross-slit FWHM/pitch, so fiber_blur_eval.py can find the largest sigma that
still meets R-AND-88.

Per order at the 20/50/80% columns, middle fiber: disk (unit circular slit
through the HDF affine transform) + ZEMAX PSF + 1-px pixel + sigma*N(0,1) in
x and y. Points are decomposed in oblique coordinates: along the slit axis
(lines of constant wavelength -> cross-slit FWHM, compared with the pitch from
fiber_pitch.py) and along the trace (tilt-corrected LSF -> R = lambda/FWHM with
the local dispersion from a cubic spline of the trace). Also stored: sampling
(LSF FWHM projected on x, px) and the flux fraction beyond +-pitch/2.

Result (2026-09): model R ~116-129k; sigma allowed ~1.45-1.85 px (UBVRIZ),
~0.7-0.85 px (YJH), giving FWHM/pitch ~1.0-1.15 median (<=1.21 max) and
21-33% of the flux beyond +-pitch/2. R-AND-88 does not constrain crosstalk.

    uv run python scripts/fiber_blur_budget.py -o blur_budget.npy   # ~10 min
    uv run python scripts/fiber_blur_eval.py blur_budget.npy
"""
import argparse
import sys

import h5py
import numpy as np
from scipy.interpolate import CubicSpline

from fiber_pitch import FRACS, MODELS, analyse, hdf_path
from fiber_profile_fwhm import fwhm

SIG = np.round(np.arange(0, 2.01, 0.05), 2)


def splines(f, fib, o):
    d = f[f'CCD_1/fiber_{fib}/order{o}'][()]
    w = d['wavelength'].astype(float)
    s = np.argsort(w)
    return {k: CubicSpline(w[s], d[k].astype(float)[s]) for k in
            ['translation_x', 'translation_y', 'rotation', 'scale_x', 'scale_y', 'shear']}


def band_rows(b, n, rng):
    Nx, Ny, fibs, orders, out = analyse(b)
    fib = fibs[len(fibs) // 2]
    rows = []
    with h5py.File(hdf_path(b), 'r') as f:
        pix = float(f['CCD_1'].attrs['pixelsize'])
        for fr in FRACS:
            for (o, wl, d, _, _) in out[fr]['per']:
                sp = splines(f, fib, o)
                p = {k: float(v(wl)) for k, v in sp.items()}
                t = np.array([sp['translation_x'](wl, 1), sp['translation_y'](wl, 1)])  # px/um
                D = np.hypot(*t)
                t /= D
                s0, s1 = splines(f, fibs[0], o), splines(f, fibs[-1], o)
                ax = np.array([s1['translation_x'](wl) - s0['translation_x'](wl),
                               s1['translation_y'](wl) - s0['translation_y'](wl)])
                ax /= np.hypot(*ax)
                Minv = np.linalg.inv(np.c_[ax, t])

                r = np.sqrt(rng.random(n)) / 2
                ph = rng.random(n) * 2 * np.pi
                x, y = r * np.cos(ph), r * np.sin(ph)
                X = p['scale_x'] * np.cos(p['rotation']) * x - p['scale_y'] * np.sin(p['rotation'] + p['shear']) * y
                Y = p['scale_x'] * np.sin(p['rotation']) * x + p['scale_y'] * np.cos(p['rotation'] + p['shear']) * y
                g = f[f'CCD_1/fiber_{fib}/psf_order_{o}']
                ks = list(g)
                ds = g[ks[np.argmin([abs(g[k].attrs['wavelength'] - wl) for k in ks])]]
                psf = np.clip(ds[()], 0, None)
                psf /= psf.sum()
                dsp = float(ds.attrs['dataSpacing']) / pix
                i0, i1 = np.unravel_index(rng.choice(psf.size, n, p=psf.ravel()), psf.shape)
                # pyechelle maps PSF axis 0 -> detector x
                X = X + (i0 + rng.random(n)) * dsp + rng.random(n)
                Y = Y + (i1 + rng.random(n)) * dsp + rng.random(n)
                gx, gy = rng.standard_normal(n), rng.standard_normal(n)
                pitch = np.median(d)
                for s in SIG:
                    A, B = Minv @ np.vstack([X + s * gx, Y + s * gy])
                    fb = fwhm(B)
                    rows.append(dict(o=o, fr=fr, wl=wl, s=s, R=wl / (fb / D), samp=fb * abs(t[0]),
                                     ratio=fwhm(A) / pitch,
                                     leak=np.mean(abs(A - np.median(A)) > pitch / 2)))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('-o', '--output', default='blur_budget.npy')
    ap.add_argument('-n', type=int, default=200_000, help='photons per order/column point')
    ap.add_argument('--seed', type=int, default=3)
    ap.add_argument('--bands', nargs='+', default=list(MODELS))
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    res = {}
    for b in args.bands:
        res[b] = band_rows(b, args.n, rng)
        print(b, 'done', file=sys.stderr)
    np.save(args.output, res, allow_pickle=True)


if __name__ == '__main__':
    main()
