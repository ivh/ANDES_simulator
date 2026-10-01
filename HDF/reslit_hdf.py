"""Build a PyEchelle HDF for a new slit layout from an existing one, without Zemax.

The source HDF's fibers are treated as samples of the spectrograph's field at known
slit positions. Per order and wavelength sample, the affine transformation is split
into the image of the fiber centre and the linear part per micron of fiber size;
both are interpolated along the slit (cubic spline) and re-assembled for each new
fiber position and size. PSFs are taken from the nearest source fiber and stored as
HDF5 hard links when identical, so the output is small.

Source slit positions are read from the `slit_position_um` attributes if the source
was made by this script, otherwise they must be given via the same layout options
make_hdf.py was run with (--src-*).

Examples:
    # MOSAIC NIR LR-J, real 150 um cores at the same 203 um pitch
    uv run python HDF/reslit_hdf.py HDF/MOSAIC_NIR_LR_J.hdf out.hdf \\
        --src-nbundles 90 --src-fibers-per-bundle 7 --src-pitch 203 --src-gap 203 \\
        --nbundles 90 --fibers-per-bundle 7 --pitch 203 --gap 203 --fiber-size 150

    # arbitrary slit: one line per fiber, "position_um [size_um]"
    uv run python HDF/reslit_hdf.py HDF/MOSAIC_NIR_LR_J.hdf out.hdf \\
        --src-nbundles 90 --src-fibers-per-bundle 7 --src-pitch 203 --src-gap 203 \\
        --positions my_slit.txt --fiber-size 150
"""

import argparse
import hashlib
import sys

import h5py
import numpy as np
from scipy.interpolate import CubicSpline

TF_FIELDS = ["rotation", "scale_x", "scale_y", "shear", "translation_x", "translation_y", "wavelength"]


def bundle_positions(nbundles, fibers_per_bundle, pitch, gap):
    """Same convention as make_hdf.py: centred, gap added after each bundle."""
    pos, y = [], 0.0
    for _ in range(nbundles):
        for _ in range(fibers_per_bundle):
            pos.append(y)
            y += pitch
        y += gap
    pos = np.array(pos)
    return pos - pos.mean()


def to_matrix(tf):
    """(n, ..., 7) structured fields -> linear part m (..., 2, 2) and translation t (..., 2)."""
    rot, sx, sy, sh = tf["rotation"], tf["scale_x"], tf["scale_y"], tf["shear"]
    m = np.stack([
        np.stack([sx * np.cos(rot), -sy * np.sin(rot + sh)], -1),
        np.stack([sx * np.sin(rot), sy * np.cos(rot + sh)], -1),
    ], -2)
    t = np.stack([tf["translation_x"], tf["translation_y"]], -1)
    return m.astype(np.float64), t.astype(np.float64)


def from_matrix(m, t, rot_ref, beta_ref):
    """Inverse of to_matrix, matching pyechelle.optics.decompose_affine_matrix.

    Angles are unwrapped towards the reference values so the output keeps the source's
    branch (make_hdf.py shifts rotation by 2pi to avoid wrapping in wavelength splines).
    """
    sx = np.hypot(m[..., 0, 0], m[..., 1, 0])
    sy = np.hypot(m[..., 0, 1], m[..., 1, 1])
    rot = np.arctan2(m[..., 1, 0], m[..., 0, 0])
    beta = np.arctan2(-m[..., 0, 1], m[..., 1, 1])
    rot += 2 * np.pi * np.round((rot_ref - rot) / (2 * np.pi))
    beta += 2 * np.pi * np.round((beta_ref - beta) / (2 * np.pi))
    return rot, sx, sy, beta - rot, t[..., 0], t[..., 1]


def read_source(h5, src_positions, src_size):
    ccd = h5["CCD_1"]
    fibers = sorted(int(k[6:]) for k in ccd if k.startswith("fiber_"))
    attrs0 = ccd[f"fiber_{fibers[0]}"].attrs
    if "slit_position_um" in attrs0:
        pos = np.array([ccd[f"fiber_{i}"].attrs["slit_position_um"] for i in fibers])
        size = np.array([ccd[f"fiber_{i}"].attrs["fiber_size_um"] for i in fibers])
    else:
        if src_positions is None:
            sys.exit("Source HDF has no slit_position_um attributes; give its layout via --src-*")
        if len(src_positions) != len(fibers):
            sys.exit(f"Source layout has {len(src_positions)} fibers, HDF has {len(fibers)}")
        pos, size = src_positions, np.full(len(fibers), src_size)
    if np.any(np.diff(pos) <= 0):
        sys.exit("Source slit positions must increase with fiber number")
    return fibers, pos, size


def interpolate_order(ccd, fibers, src_pos, src_size, order, new_pos, new_size):
    tf = np.stack([ccd[f"fiber_{i}/order{order}"][()] for i in fibers])  # (nfib, nwl) structured
    wl = tf["wavelength"]
    if not np.allclose(wl, wl[0]):
        sys.exit(f"order {order}: wavelength samples differ between source fibers; not supported yet")

    m, t = to_matrix(tf)
    centre = t + 0.5 * (m[..., 0] + m[..., 1])  # image of the box centre (0.5, 0.5)
    m_per_um = m / src_size[:, None, None, None]

    def spline(y):
        return CubicSpline(src_pos, y, axis=0)(new_pos)

    m_new = spline(m_per_um) * new_size[:, None, None, None]
    c_new = spline(centre)
    t_new = c_new - 0.5 * (m_new[..., 0] + m_new[..., 1])

    rot_ref = spline(tf["rotation"].astype(np.float64))
    beta_ref = spline((tf["rotation"] + tf["shear"]).astype(np.float64))
    params = from_matrix(m_new, t_new, rot_ref, beta_ref)

    out = np.empty((len(new_pos), wl.shape[1]), dtype=tf.dtype)
    for name, val in zip(TF_FIELDS[:6], params):
        out[name] = val
    out["wavelength"] = wl[0]
    return out


def psf_digest(group):
    h = hashlib.sha1()
    for k in sorted(group):
        d = group[k]
        h.update(k.encode())
        h.update(d[()].tobytes())
        h.update(repr(sorted(d.attrs.items())).encode())
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("source")
    p.add_argument("output")

    g = p.add_argument_group("source layout (as given to make_hdf.py; ignored if the HDF carries positions)")
    g.add_argument("--src-nbundles", type=int)
    g.add_argument("--src-fibers-per-bundle", type=int)
    g.add_argument("--src-pitch", type=float, help="make_hdf.py --fiber-size (used as pitch and field size)")
    g.add_argument("--src-gap", type=float, default=0.0, help="make_hdf.py --bundle-gap")
    g.add_argument("--src-nfibers", type=int, help="for make_hdf.py --layout linear")

    g = p.add_argument_group("new layout")
    g.add_argument("--nbundles", type=int)
    g.add_argument("--fibers-per-bundle", type=int)
    g.add_argument("--nfibers", type=int, help="linear layout: number of fibers at --pitch")
    g.add_argument("--pitch", type=float)
    g.add_argument("--gap", type=float, default=0.0)
    g.add_argument("--positions", help="text file, one fiber per line: position_um [size_um]")
    g.add_argument("--fiber-size", type=float, help="field size in microns (default: --pitch)")
    g.add_argument("--field-shape", help="override field_shape attribute (default: from source)")
    g.add_argument("--extrapolate", type=float, default=0.0,
                   help="allowed distance beyond the outermost source fibers, in microns (default: 0)")
    args = p.parse_args()

    src_positions, src_size = None, args.src_pitch
    if args.src_nbundles:
        src_positions = bundle_positions(args.src_nbundles, args.src_fibers_per_bundle, args.src_pitch, args.src_gap)
    elif args.src_nfibers:
        src_positions = bundle_positions(1, args.src_nfibers, args.src_pitch, 0.0)

    if args.positions:
        data = np.loadtxt(args.positions, ndmin=2)
        new_pos = data[:, 0]
        new_size = data[:, 1] if data.shape[1] > 1 else np.full(len(new_pos), args.fiber_size or np.nan)
    elif args.nbundles:
        new_pos = bundle_positions(args.nbundles, args.fibers_per_bundle, args.pitch, args.gap)
        new_size = np.full(len(new_pos), args.fiber_size or args.pitch)
    elif args.nfibers:
        new_pos = bundle_positions(1, args.nfibers, args.pitch, 0.0)
        new_size = np.full(len(new_pos), args.fiber_size or args.pitch)
    else:
        sys.exit("Give the new layout via --nbundles, --nfibers or --positions")
    if np.any(~np.isfinite(new_size)):
        sys.exit("Fiber size missing: give --fiber-size or a size column in --positions")

    with h5py.File(args.source, "r") as src, h5py.File(args.output, "w") as dst:
        ccd = src["CCD_1"]
        fibers, src_pos, src_sizes = read_source(src, src_positions, src_size)

        lo, hi = src_pos[0] - args.extrapolate, src_pos[-1] + args.extrapolate
        outside = (new_pos < lo) | (new_pos > hi)
        if outside.any():
            sys.exit(f"{outside.sum()} fibers outside the sampled slit [{src_pos[0]:.1f}, {src_pos[-1]:.1f}] um "
                     f"(use --extrapolate to allow)")

        orders = sorted(int(k[5:]) for k in ccd[f"fiber_{fibers[0]}"] if k.startswith("order"))
        nearest = np.abs(new_pos[:, None] - src_pos[None, :]).argmin(axis=1)
        print(f"{len(fibers)} source fibers ({src_pos[0]:.1f}..{src_pos[-1]:.1f} um) -> "
              f"{len(new_pos)} fibers ({new_pos.min():.1f}..{new_pos.max():.1f} um), orders {orders}")

        dccd = dst.create_group("CCD_1")
        for k, v in ccd.attrs.items():
            dccd.attrs[k] = v
        src.copy(ccd["Spectrograph"], dccd, "Spectrograph")
        dst.attrs["reslit_source"] = str(args.source)

        tfs = {o: interpolate_order(ccd, fibers, src_pos, src_sizes, o, new_pos, new_size) for o in orders}

        psf_store = {}  # digest -> path of first copy in dst
        digests = {}
        for j in range(len(new_pos)):
            sfib = fibers[nearest[j]]
            g = dccd.create_group(f"fiber_{j + 1}")
            for k, v in ccd[f"fiber_{sfib}"].attrs.items():
                if k not in ("slit_position_um", "fiber_size_um"):
                    g.attrs[k] = v
            if args.field_shape:
                g.attrs["field_shape"] = args.field_shape
            g.attrs["slit_position_um"] = new_pos[j]
            g.attrs["fiber_size_um"] = new_size[j]
            for o in orders:
                g.create_dataset(f"order{o}", data=tfs[o][j])
                key = f"psf_order_{o}"
                sgroup = ccd[f"fiber_{sfib}/{key}"]
                if (sfib, o) not in digests:
                    digests[sfib, o] = psf_digest(sgroup)
                digest = digests[sfib, o]
                if digest in psf_store:
                    g[key] = dst[psf_store[digest]]  # hard link
                else:
                    src.copy(sgroup, g, key)
                    psf_store[digest] = g[key].name
        print(f"wrote {args.output}: {len(psf_store)} distinct PSF sets")


if __name__ == "__main__":
    main()
