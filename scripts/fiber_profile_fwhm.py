"""Cross-slit FWHM of a single fiber image vs. the fiber pitch, from the HDF models.

Question answered: is the fiber pitch along the pseudo-slit (scripts/fiber_pitch.py)
the same as the FWHM of each fiber's spatial profile? It is not.

Monte-Carlo of PyEchelle's own photon chain for the middle fiber: uniform
points in the unit circular slit (pyechelle.slit.circular) -> HDF affine
transform (rotation/scale/shear) -> offset drawn from the nearest ZEMAX PSF
(dataSpacing um per PSF pixel; PSF axis 0 maps to detector x, as in
pyechelle.raytracing) -> optional 1-px boxcar. The profile is projected on the
slit axis, at the 20/50/80% columns in every order. Columns printed: pitch,
scale_y, FWHM of disk only / PSF only / disk*PSF / disk*PSF*pixel, and the
FWHM/pitch ratio range.

Result (2026-09): FWHM ~0.82 x pitch in all bands. The profile is the
projected disk (semi-ellipse, FWHM = sqrt(3)/2 D) slightly narrowed by the
0.4-0.9 px ZEMAX PSF -- blurring a strongly concave edge lowers the half-max
points faster than the peak. Flat-topped, not Gaussian. --leak prints the flux
fraction beyond +-pitch/2 (2.4-5% median, up to 8% worst order), which is
pure PSF wing because the disks touch but do not overlap.

    uv run python scripts/fiber_profile_fwhm.py [--leak]
"""
import argparse

import h5py
import numpy as np
from scipy.ndimage import gaussian_filter1d

from fiber_pitch import FRACS, MODELS, analyse, hdf_path

N = 400_000


def fwhm(v, bw=0.01):
    h, e = np.histogram(v, bins=np.arange(v.min() - 0.1, v.max() + 0.1, bw))
    c = 0.5 * (e[1:] + e[:-1])
    h = gaussian_filter1d(h.astype(float), 2)
    m = h.max() / 2
    i = np.where(h >= m)[0]
    a, b = i[0], i[-1]
    return (np.interp(m, [h[b + 1], h[b]], [c[b + 1], c[b]])
            - np.interp(m, [h[a - 1], h[a]], [c[a - 1], c[a]]))


def params(f, fib, o, wl, keys):
    d = f[f'CCD_1/fiber_{fib}/order{o}'][()]
    w = d['wavelength'].astype(float)
    s = np.argsort(w)
    return {k: np.interp(wl, w[s], d[k].astype(float)[s]) for k in keys}


def slit_axis(f, fibs, o, wl):
    P = [list(params(f, fi, o, wl, ['translation_x', 'translation_y']).values())
         for fi in (fibs[0], fibs[-1])]
    ax = np.subtract(P[1], P[0])
    return ax / np.hypot(*ax)


def profiles(f, fib, o, wl, pix, ax, rng):
    """Projected samples: disk only, PSF only, disk*PSF, disk*PSF*pixel."""
    p = params(f, fib, o, wl, ['rotation', 'scale_x', 'scale_y', 'shear'])
    r = np.sqrt(rng.random(N)) / 2
    ph = rng.random(N) * 2 * np.pi
    x, y = r * np.cos(ph), r * np.sin(ph)
    X = p['scale_x'] * np.cos(p['rotation']) * x - p['scale_y'] * np.sin(p['rotation'] + p['shear']) * y
    Y = p['scale_x'] * np.sin(p['rotation']) * x + p['scale_y'] * np.cos(p['rotation'] + p['shear']) * y
    geo = X * ax[0] + Y * ax[1]
    g = f[f'CCD_1/fiber_{fib}/psf_order_{o}']
    ks = list(g)
    pw = np.array([g[k].attrs['wavelength'] for k in ks])
    ds = g[ks[np.argmin(abs(pw - wl))]]
    psf = np.clip(ds[()], 0, None)
    psf /= psf.sum()
    sp = float(ds.attrs['dataSpacing']) / pix
    dx, dy = np.unravel_index(rng.choice(psf.size, N, p=psf.ravel()), psf.shape)
    dx = (dx + rng.random(N)) * sp
    dy = (dy + rng.random(N)) * sp
    ps = dx * ax[0] + dy * ax[1]
    return geo, ps, geo + ps, geo + ps + rng.random(N)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--leak', action='store_true',
                    help='also print flux fraction beyond +-pitch/2 (50%% column)')
    ap.add_argument('--seed', type=int, default=1)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    print(f'{"band":4s} {"pitch":>5s} {"sy":>5s} | {"disk":>5s} {"PSF":>5s} {"d*PSF":>5s} {"+pix":>5s}'
          f' | FWHM/pitch min  med  max' + ('  | leak med  max' if args.leak else ''))
    for b in MODELS:
        Nx, Ny, fibs, orders, out = analyse(b)
        fib = fibs[len(fibs) // 2]
        rows, leak = [], []
        with h5py.File(hdf_path(b), 'r') as f:
            pix = float(f['CCD_1'].attrs['pixelsize'])
            for fr in FRACS:
                for (o, wl, d, _, _) in out[fr]['per']:
                    ax = slit_axis(f, fibs, o, wl)
                    sy = params(f, fib, o, wl, ['scale_y'])['scale_y']
                    prof = profiles(f, fib, o, wl, pix, ax, rng)
                    rows.append([np.median(d), sy] + [fwhm(v) for v in prof])
                    if args.leak and fr == 0.5:
                        t = prof[2] - np.median(prof[2])
                        leak.append(np.mean(abs(t) > np.median(d) / 2))
        R = np.array(rows)
        rat = R[:, 5] / R[:, 0]
        M = np.median(R, axis=0)
        line = (f'{b:4s} {M[0]:5.2f} {M[1]:5.2f} | {M[2]:5.2f} {M[3]:5.2f} {M[4]:5.2f} {M[5]:5.2f}'
                f' |            {rat.min():.2f} {np.median(rat):.2f} {rat.max():.2f}')
        if args.leak:
            line += f'  | {np.median(leak) * 100:4.1f}% {max(leak) * 100:4.1f}%'
        print(line)


if __name__ == '__main__':
    main()
