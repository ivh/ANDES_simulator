"""Fiber pitch along the pseudo-slit, measured from the PyEchelle HDF models.

For each ANDES band (U B V R IZ Y J H, default HDF model per band) this finds,
in every order, the wavelength at which the fiber traces cross detector columns
at 20/50/80% of the detector width, evaluates all fiber-image centres
(translation_x/y of the HDF affine transforms) there, projects them onto the
principal slit axis and differences adjacent fibers. Prints min/median/max
pitch in detector pixels per band and column, pooled over orders and fiber
pairs.

Result (2026-09): pitch ~4.5-5.8 px (UBVRIZ), ~2.25-2.46 px (YJH); <1% change
with column, +-4-8% order-to-order, no designed gaps along the slit. Pitch
equals the fiber image diameter (scale_y) to 0.1% in every order, i.e. the
models have contiguous fibers with no cladding/buffer between them.

    uv run python scripts/fiber_pitch.py
"""
from pathlib import Path

import h5py
import numpy as np

MODELS = {
    'U': 'ANDES_U_v88', 'B': 'ANDES_B_v88', 'V': 'ANDES_V_v88',
    'R': 'ANDES_123_R3', 'IZ': 'ANDES_123_IZ3',
    'Y': 'ANDES_75fibre_Y', 'J': 'ANDES_75fibre_J', 'H': 'ANDES_75fibre_H',
}
BASE = Path(__file__).resolve().parent.parent / 'HDF'
FRACS = (0.2, 0.5, 0.8)


def hdf_path(band):
    return BASE / f'{MODELS[band]}.hdf'


def band_data(band):
    with h5py.File(hdf_path(band), 'r') as f:
        ccd = f['CCD_1']
        Nx, Ny = int(ccd.attrs['Nx']), int(ccd.attrs['Ny'])
        fibs = sorted(int(k.split('_')[1]) for k in ccd if k.startswith('fiber_'))
        orders = sorted(int(k[5:]) for k in ccd[f'fiber_{fibs[0]}'] if k.startswith('order'))
        data = {}
        for fi in fibs:
            for o in orders:
                d = ccd[f'fiber_{fi}/order{o}'][()]
                data[(fi, o)] = (np.asarray(d['wavelength'], float),
                                 np.asarray(d['translation_x'], float),
                                 np.asarray(d['translation_y'], float))
    return Nx, Ny, fibs, orders, data


def interp(w, v, wt):
    if w[0] > w[-1]:
        w, v = w[::-1], v[::-1]
    return np.interp(wt, w, v)


def analyse(band, fracs=FRACS):
    """Per column fraction: pooled pitch stats plus per-order (order, wl, gaps, fibers, off-axis ratio)."""
    Nx, Ny, fibs, orders, data = band_data(band)
    out = {}
    for fr in fracs:
        col = fr * Nx
        gaps, per = [], []
        for o in orders:
            ws = np.sort(data[(fibs[0], o)][0])
            tx_med = np.median([interp(data[(fi, o)][0], data[(fi, o)][1], ws) for fi in fibs], axis=0)
            if not (tx_med.min() <= col <= tx_med.max()):
                continue
            s = np.argsort(tx_med)
            wl = np.interp(col, tx_med[s], ws[s])
            x = np.array([interp(data[(fi, o)][0], data[(fi, o)][1], wl) for fi in fibs])
            y = np.array([interp(data[(fi, o)][0], data[(fi, o)][2], wl) for fi in fibs])
            pts = np.c_[x - x.mean(), y - y.mean()]
            # slit is slightly bowed, so project on its principal axis rather than on y
            _, sv, vt = np.linalg.svd(pts, full_matrices=False)
            p = pts @ vt[0]
            idx = np.argsort(p)
            d = np.diff(p[idx])
            gaps.append(d)
            per.append((o, wl, d, np.array(fibs)[idx], sv[1] / sv[0]))
        if not gaps:
            out[fr] = None
            continue
        allg = np.concatenate(gaps)
        out[fr] = dict(n_orders=len(gaps), med=np.median(allg), mx=allg.max(), mn=allg.min(),
                       per=per, col=col)
    return Nx, Ny, fibs, orders, out


def main():
    print(f'{"band":4s} {"nfib":>4s} {"nord":>4s}  ' +
          '  '.join(f'col {fr:.0%} min/med/max' for fr in FRACS))
    for b in MODELS:
        Nx, Ny, fibs, orders, out = analyse(b)
        cells = []
        for fr in FRACS:
            r = out[fr]
            cells.append('   (no orders)   ' if r is None else
                         f'{r["mn"]:5.3f}/{r["med"]:5.3f}/{r["mx"]:5.3f}')
        print(f'{b:4s} {len(fibs):4d} {len(orders):4d}  ' + '  '.join(f'{c:>22s}' for c in cells))


if __name__ == '__main__':
    main()
