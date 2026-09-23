"""Evaluate fiber_blur_budget.py output against R-AND-88 (TRS ESO-391757).

Prints (1) the models as they are (sigma=0): R min/mean/max, minimum sampling
in px per resolution element (compare R-AND-89: >=4 px below 950 nm, R-AND-90:
>=2.5 px above), FWHM/pitch; and (2) for two readings of R-AND-88 -- every
point R>=90k ("min90"), mean R>=100k ("avg100") -- the largest added isotropic
Gaussian sigma that still passes, with the resulting cross-slit FWHM/pitch and
flux fraction beyond +-pitch/2. The mean is over the sampled order/column
points, not flux-weighted. A '>' marks sigma hitting the top of the grid.

    uv run python scripts/fiber_blur_eval.py blur_budget.npy
"""
import argparse

import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('input', help='.npy written by fiber_blur_budget.py')
    res = np.load(ap.parse_args().input, allow_pickle=True).item()

    print('sigma=0 (model as is)')
    print(f'{"band":4s} {"Rmin":>6s} {"Rmean":>6s} {"Rmax":>6s} {"samp_min":>8s} {"FWHM/p med":>10s}')
    for b, rows in res.items():
        r0 = [r for r in rows if r['s'] == 0]
        R = np.array([r['R'] for r in r0])
        print(f'{b:4s} {R.min() / 1e3:6.1f} {R.mean() / 1e3:6.1f} {R.max() / 1e3:6.1f} '
              f'{min(r["samp"] for r in r0):8.2f} {np.median([r["ratio"] for r in r0]):10.2f}')

    print('\nlargest added isotropic Gaussian sigma [px] meeting R-AND-88, and resulting cross-slit numbers')
    print(f'{"band":4s} {"crit":6s} {"sigma":>5s} {"FWHMg":>5s} {"Rmin":>6s} {"Rmean":>6s} | '
          f'{"FWHM/p med":>10s} {"max":>5s} | {"leak med":>8s} {"max":>5s}')
    crits = [('min90', lambda R: R.min() >= 90e3), ('avg100', lambda R: R.mean() >= 100e3)]
    for b, rows in res.items():
        sig = sorted({r['s'] for r in rows})
        by = {s: [r for r in rows if r['s'] == s] for s in sig}
        for crit, ok in crits:
            good = [s for s in sig if ok(np.array([r['R'] for r in by[s]]))]
            if not good:
                print(f'{b:4s} {crit:6s}  none (already violated at sigma=0)')
                continue
            s = max(good)
            rr = by[s]
            R = np.array([r['R'] for r in rr])
            ra = np.array([r['ratio'] for r in rr])
            lk = np.array([r['leak'] for r in rr])
            edge = '>' if s == sig[-1] else ' '
            print(f'{b:4s} {crit:6s} {edge}{s:4.2f} {2.3548 * s:5.2f} {R.min() / 1e3:6.1f} {R.mean() / 1e3:6.1f} | '
                  f'{np.median(ra):10.2f} {ra.max():5.2f} | {np.median(lk) * 100:7.1f}% {lk.max() * 100:4.1f}%')


if __name__ == '__main__':
    main()
