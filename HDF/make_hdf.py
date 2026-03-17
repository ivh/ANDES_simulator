# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "pyechelle>=0.4",
#     "zospy",
#     "h5py",
#     "numpy",
# ]
# ///
"""Convert a Zemax optical model to a PyEchelle HDF spectrograph model.

Connects to OpticStudio via ZOS-API, configures fibers and diffraction orders,
then samples transformations and PSFs across each order and writes them to HDF.

Usage:
    # Simple linear slit (e.g. ANDES)
    uv run make_hdf.py model.zmx output.hdf --orders 90 108 --layout linear --nfibers 75

    # Bundled slit (e.g. MOSAIC: 90 mini-IFU bundles of 7 fibers)
    uv run make_hdf.py model.zmx output.hdf --orders 50 80 --layout bundles \
        --nbundles 90 --fibers-per-bundle 7 --bundle-gap 200

    # Quick API connectivity test
    uv run make_hdf.py model.zmx output.hdf --orders 90 108 --test-api
"""

import argparse
import sys
import time

import numpy as np
from pyechelle.CCD import CCD
from pyechelle.hdfbuilder import HDFBuilder
from pyechelle.spectrograph import InteractiveZEMAX


def connect(name, zmx_path):
    print(f"Connecting to OpticStudio, loading {zmx_path} ...")
    zmx = InteractiveZEMAX(name=name, zemax_filepath=zmx_path)
    print("Connected.")
    return zmx


def fiber_positions_linear(nfibers, fiber_pitch):
    """Evenly spaced fibers in a vertical slit (e.g. ANDES)."""
    y = (np.arange(nfibers) - nfibers / 2 + 0.5) * fiber_pitch
    return y


def fiber_positions_bundles(nbundles, fibers_per_bundle, fiber_pitch, bundle_gap):
    """Bundles of fibers with gaps between them (e.g. MOSAIC).

    Each bundle is a contiguous group of fibers at fiber_pitch spacing.
    Between bundles there is an additional bundle_gap (in microns).

    TODO: fill in actual values from MOSAIC fiber link design:
    - fiber_pitch: center-to-center within bundle (fiber diameter + cladding)
    - bundle_gap: extra gap between last fiber of bundle N and first of bundle N+1
    """
    positions = []
    y = 0.0
    for b in range(nbundles):
        for f in range(fibers_per_bundle):
            positions.append(y)
            y += fiber_pitch
        y += bundle_gap
    positions = np.array(positions)
    positions -= positions.mean()  # center on slit
    return positions


def ensure_detector_aperture(zmx, nx, ny, pixelsize):
    """Add a rectangular aperture to the image surface if it doesn't have one.

    PyEchelle's wavelength walk algorithm needs the detector surface to vignette
    rays that fall off the CCD edge. Without it, the walk overshoots.
    """
    import zospy
    lde = zmx._oss.LDE
    det = lde.GetSurfaceAt(lde.NumberOfSurfaces)
    ap = det.ApertureData
    if str(ap.CurrentType) == "None":
        xhw = nx * pixelsize / 2000.0   # half-width in mm
        yhw = ny * pixelsize / 2000.0
        rect_type = zospy.constants.Editors.LDE.SurfaceApertureTypes.RectangularAperture
        settings = ap.CreateApertureTypeSettings(rect_type)
        settings._S_RectangularAperture.XHalfWidth = xhw
        settings._S_RectangularAperture.YHalfWidth = yhw
        ap.ChangeApertureTypeSettings(settings)
        print(f"Added rectangular aperture to detector: {xhw:.1f} x {yhw:.1f} mm half-width")


def configure(zmx, args):
    if args.config is not None:
        mce = zmx._oss.MCE
        n_configs = mce.NumberOfConfigurations
        if args.config < 1 or args.config > n_configs:
            raise ValueError(f"Config {args.config} out of range (1-{n_configs})")
        # Navigate to the target config
        while mce.CurrentConfiguration != args.config:
            mce.NextConfiguration()
        print(f"Switched to configuration {args.config} "
              f"(of {n_configs} total)")

    # Find the active grating surface: pick the diffraction grating with
    # nonzero groove density (MCE sets inactive gratings to 0).
    # This handles multiple gratings with the same or different comment names.
    grating_surface = args.grating_surface
    if isinstance(grating_surface, str):
        # Scan all surfaces for diffraction gratings with nonzero groove density
        lde = zmx._lde
        active = None
        for i in range(lde.NumberOfSurfaces):
            s = lde.GetSurfaceAt(i)
            if "Diffraction" in s.TypeName:
                gd = s.GetCellAt(12).DoubleValue
                if gd > 0:
                    active = i
                    print(f"Found active grating at surface {i} "
                          f"({gd:.4f} l/um, \"{s.Comment}\")")
                    break
        if active is not None:
            grating_surface = active
    zmx.set_grating(surface=grating_surface, blaze=args.blaze)
    if args.groove_density is not None:
        zmx._groves_per_micron = args.groove_density / 1000.0
        zmx._zmx_grating.GetCellAt(12).DoubleValue = args.groove_density / 1000.0
        print(f"Overrode groove density to {args.groove_density} l/mm "
              f"({args.groove_density / 1000.0} l/um)")
    ensure_detector_aperture(zmx, args.nx, args.ny, args.pixelsize)
    if args.skip_ccd_check:
        # Bypass _check_ccd for models without a rectangular detector aperture
        zmx._ccd.update({1: CCD(args.nx, args.ny, pixelsize=args.pixelsize)})
    else:
        zmx.add_ccd(1, CCD(args.nx, args.ny, pixelsize=args.pixelsize))

    orders = list(range(args.orders[0], args.orders[1]))

    if args.layout == "linear":
        yfield = fiber_positions_linear(args.nfibers, args.fiber_size)
        nfibers = args.nfibers
    elif args.layout == "bundles":
        yfield = fiber_positions_bundles(
            args.nbundles, args.fibers_per_bundle,
            args.fiber_size, args.bundle_gap,
        )
        nfibers = len(yfield)

    # Auto-detect slit axis from native field coordinates: if the model's
    # fields span more in X than Y, the slit runs along X.
    fields = zmx._oss.SystemData.Fields
    nf = fields.NumberOfFields
    if nf > 1:
        xs = [fields.GetField(i).X for i in range(1, nf + 1)]
        ys = [fields.GetField(i).Y for i in range(1, nf + 1)]
        x_span = max(xs) - min(xs)
        y_span = max(ys) - min(ys)
        slit_along_x = x_span > y_span
    else:
        slit_along_x = False
    if slit_along_x:
        print(f"Slit axis: X (native field X span {x_span:.1f} > Y span {y_span:.1f})")
    else:
        print(f"Slit axis: Y (native field Y span {y_span:.1f} >= X span {x_span:.1f})")

    # positions are in microns, add_field expects mm
    for i, y in enumerate(yfield):
        pos = y / 1000
        if slit_along_x:
            zmx.add_field(pos, 0.0, args.fiber_size, args.fiber_size,
                           shape=args.fiber_shape, name="fiber")
        else:
            zmx.add_field(0.0, pos, args.fiber_size, args.fiber_size,
                           shape=args.fiber_shape, name="fiber")
        zmx.set_orders(1, i + 1, orders)

    zmx.psf_settings(
        image_delta=args.psf_delta,
        image_sampling=args.psf_image_sampling,
        pupil_sampling=args.psf_pupil_sampling,
    )

    n_calls = nfibers * len(orders) * (args.n_samples + args.n_psf_samples)
    print(f"Configured: {nfibers} fibers, orders {args.orders[0]}-{args.orders[1]-1} "
          f"({len(orders)} orders), ~{n_calls} Zemax API calls")


def patch_wavelength_range(zmx, wl_min, wl_max):
    """Override get_wavelength_range to return a fixed range.

    Needed when SingleRayNormUnpol doesn't work (old OpticStudio versions),
    which breaks pyechelle's walk-to-detector-edge algorithm.
    """
    def fixed_get_wavelength_range(order=None, fiber=None, ccd_index=None):
        if order is not None:
            zmx._zos_set_current_order(order)
        if fiber is not None:
            zmx._zos_set_current_field(fiber)
        return wl_min, wl_max

    zmx.get_wavelength_range = fixed_get_wavelength_range
    print(f"Patched wavelength range to {wl_min:.4f} - {wl_max:.4f} um")


def test_api(zmx, args):
    """Quick sanity check: get one transformation and one PSF."""
    orders = list(range(args.orders[0], args.orders[1]))
    mid_order = orders[len(orders) // 2]
    wl_range = zmx.get_wavelength_range(mid_order, 1, 1)
    mid_wl = (wl_range[0] + wl_range[1]) / 2
    print(f"Order {mid_order}: wavelength range {wl_range[0]:.4f} - {wl_range[1]:.4f} um")

    t = zmx.get_transformation(mid_wl, mid_order, 1, 1)
    print(f"Transformation at {mid_wl:.4f} um: {t}")

    psf = zmx.get_psf(mid_wl, mid_order, 1, 1)
    print(f"PSF at {mid_wl:.4f} um: shape {np.array(psf.data).shape}, "
          f"sampling {psf.sampling} um")


def fix_dispersion_axis(output_path):
    """Ensure dispersion is along tx (horizontal).

    Some Zemax models (e.g. NIR) disperse along ty instead of tx. Detect this
    by checking which translation component varies more with wavelength in the
    first fiber, and swap the x/y axes across the entire HDF if needed.
    """
    import h5py

    with h5py.File(output_path, "r+") as f:
        ccd = f["CCD_1"]
        fiber_key = next(k for k in ccd.keys() if k.startswith("fiber_"))
        order_key = next(
            k for k in ccd[fiber_key].keys()
            if k.startswith("order") and not k.startswith("psf_")
        )
        data = ccd[fiber_key][order_key][:]
        tx = data["translation_x"]
        ty = data["translation_y"]
        tx_range = float(tx.max() - tx.min())
        ty_range = float(ty.max() - ty.min())

        if ty_range > tx_range:
            print(f"Dispersion along ty (range {ty_range:.1f} vs tx {tx_range:.1f}) "
                  f"— swapping x/y axes")
            for fk in ccd.keys():
                if not fk.startswith("fiber_"):
                    continue
                for ok in ccd[fk].keys():
                    if not ok.startswith("order") or ok.startswith("psf_"):
                        continue
                    d = ccd[fk][ok][:]
                    d["translation_x"], d["translation_y"] = (
                        d["translation_y"].copy(),
                        d["translation_x"].copy(),
                    )
                    d["scale_x"], d["scale_y"] = (
                        d["scale_y"].copy(),
                        d["scale_x"].copy(),
                    )
                    ccd[fk][ok][:] = d
        else:
            print("Dispersion already along tx — no swap needed.")


def unwrap_rotation(output_path):
    """Fix rotation angle wrapping in the HDF file.

    PyEchelle stores rotation values that can wrap between -pi and +pi across
    wavelength samples.  Linear interpolation then passes through 0, producing
    wrong trace positions.  Shifting negative values by +2pi keeps all angles
    near pi and prevents this.
    """
    import h5py

    fixed = 0
    with h5py.File(output_path, "r+") as f:
        ccd = f["CCD_1"]
        for key in ccd.keys():
            if not key.startswith("fiber_"):
                continue
            for order_key in ccd[key].keys():
                if not order_key.startswith("order") or order_key.startswith("psf_"):
                    continue
                data = ccd[key][order_key][:]
                rot = data["rotation"]
                mask = rot < 0
                if mask.any():
                    rot[mask] += 2 * np.pi
                    data["rotation"] = rot
                    ccd[key][order_key][:] = data
                    fixed += int(mask.sum())
    if fixed:
        print(f"Unwrapped {fixed} rotation values in {output_path}")
    else:
        print("No rotation unwrapping needed.")


def build_hdf(zmx, args):
    print(f"Writing {args.output} ...")
    t0 = time.time()
    hdf = HDFBuilder(zmx, args.output)
    hdf.save_to_hdf(
        n_transformation_per_order=args.n_samples,
        n_psfs_per_order=args.n_psf_samples,
    )
    hdf.close()
    elapsed = time.time() - t0
    print(f"Done in {elapsed / 3600:.1f} hours. Wrote {args.output}")
    fix_dispersion_axis(args.output)
    unwrap_rotation(args.output)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("zmx_file", help="Path to .zmx/.zos Zemax file")
    p.add_argument("output", help="Output .hdf path")

    g = p.add_argument_group("instrument")
    g.add_argument("--name", default="spectrograph", help="Instrument name (default: spectrograph)")
    g.add_argument("--orders", type=int, nargs=2, required=True, metavar=("FIRST", "LAST"),
                   help="Order range [first, last) -- last is exclusive, as in Python range()")
    g.add_argument("--fiber-size", type=float, default=474, help="Fiber size in microns (default: 474)")
    g.add_argument("--fiber-shape", default="circular", choices=["circular", "rectangular"],
                   help="Fiber shape (default: circular)")
    g.add_argument("--blaze", type=float, default=76, help="Grating blaze angle in degrees (default: 76)")
    g.add_argument("--grating-surface", default="ECHELLE", help="Grating surface name in Zemax (default: ECHELLE)")
    g.add_argument("--groove-density", type=float, default=None,
                   help="Override grating groove density in lines/mm (reads from Zemax if not set)")

    g = p.add_argument_group("fiber layout")
    g.add_argument("--layout", default="linear", choices=["linear", "bundles"],
                   help="Fiber layout: 'linear' for evenly spaced, 'bundles' for grouped (default: linear)")
    # linear layout
    g.add_argument("--nfibers", type=int, default=75, help="Number of fibers for linear layout (default: 75)")
    # bundle layout
    g.add_argument("--nbundles", type=int, default=90, help="Number of bundles (default: 90)")
    g.add_argument("--fibers-per-bundle", type=int, default=7, help="Fibers per bundle (default: 7)")
    g.add_argument("--bundle-gap", type=float, default=200,
                   help="Extra gap between bundles in microns (default: 200, TODO: get real value)")

    g = p.add_argument_group("detector")
    g.add_argument("--nx", type=int, default=4096, help="Detector X pixels (default: 4096)")
    g.add_argument("--ny", type=int, default=4096, help="Detector Y pixels (default: 4096)")
    g.add_argument("--pixelsize", type=float, default=15, help="Pixel size in microns (default: 15)")

    g = p.add_argument_group("PSF settings")
    g.add_argument("--psf-delta", type=float, default=3, help="PSF image delta in microns (default: 3)")
    g.add_argument("--psf-image-sampling", default="128x128", help="PSF image sampling (default: 128x128)")
    g.add_argument("--psf-pupil-sampling", default="64x64", help="PSF pupil sampling (default: 64x64)")

    g = p.add_argument_group("sampling")
    g.add_argument("--n-samples", type=int, default=15,
                   help="Transformation samples per order (default: 15)")
    g.add_argument("--n-psf-samples", type=int, default=15,
                   help="PSF samples per order (default: 15)")

    p.add_argument("--config", type=int, default=None,
                   help="Zemax multi-configuration number (e.g. 2 for LR-R)")
    p.add_argument("--skip-ccd-check", action="store_true",
                   help="Skip CCD size vs detector aperture check (for models without rect aperture)")
    p.add_argument("--wl-range", type=float, nargs=2, metavar=("MIN", "MAX"),
                   help="Override wavelength range in microns (bypasses walk-to-edge algorithm)")
    p.add_argument("--test-api", action="store_true",
                   help="Just connect and test one transformation + PSF, then exit")

    args = p.parse_args()

    zmx = connect(args.name, args.zmx_file)
    configure(zmx, args)

    if args.wl_range:
        patch_wavelength_range(zmx, args.wl_range[0], args.wl_range[1])

    if args.test_api:
        test_api(zmx, args)
        print("API test passed.")
        sys.exit(0)

    build_hdf(zmx, args)


if __name__ == "__main__":
    main()
