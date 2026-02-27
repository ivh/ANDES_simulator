# /// script
# requires-python = ">=3.10"
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


def configure(zmx, args):
    zmx.set_grating(surface=args.grating_surface, blaze=args.blaze)
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

    # positions are in microns, add_field expects mm
    for i, y in enumerate(yfield):
        zmx.add_field(0.0, y / 1000, args.fiber_size, args.fiber_size,
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

    p.add_argument("--test-api", action="store_true",
                   help="Just connect and test one transformation + PSF, then exit")

    args = p.parse_args()

    zmx = connect(args.name, args.zmx_file)
    configure(zmx, args)

    if args.test_api:
        test_api(zmx, args)
        print("API test passed.")
        sys.exit(0)

    build_hdf(zmx, args)


if __name__ == "__main__":
    main()
