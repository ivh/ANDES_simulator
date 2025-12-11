# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "astropy",
#     "matplotlib",
#     "scipy",
# ]
# ///

import numpy as np
import random
import sys
from glob import glob
from astropy.io import fits
from scipy.signal import convolve2d
import concurrent.futures
import functools
from pathlib import Path

script_dir = Path(__file__).parent
project_root = script_dir.parent
BASEDIR = str(project_root.parent) + '/'
#fnschema = '{arm}_T0108_fib{fib:02d}.fits'
fnschema = '{arm}_FP_fiber{fib:02d}_shift*.fits'

def custom_gaussian_kernel(dimx=4, dimy=4, sigma=1.0, blank='random'):
    """
    Creates a 2D Gaussian kernel with specified dimensions.
    Then randomly zeros out one outer edge (top/bottom row or left/right column)
    and renormalizes the kernel.
    
    Parameters:
    -----------
    dimx : int
        Width of the kernel
    dimy : int
        Height of the kernel
    sigma : float, optional
        Standard deviation of the Gaussian distribution (default: 1.0)
    blank : str, optional
        Edge to zero out ('top', 'bottom', 'left', 'right') (default: 'random')
    
    Returns:
    --------
    numpy.ndarray
        The modified Gaussian kernel with dimensions (dimy, dimx)
    """
    # Create meshgrid for x and y coordinates, centered even for even-sized kernels
    y, x = np.ogrid[:dimy, :dimx]
    y = y - (dimy - 1) / 2
    x = x - (dimx - 1) / 2

    # Calculate the 2D Gaussian
    kernel = np.exp(-(x**2 + y**2) / (2 * sigma**2))
    
    # Randomly choose an edge to zero out
    if blank == 'random':
        edge = random.choice(['none', 'top', 'bottom', 'left', 'right'])
    else:
        edge = blank
    
    # Zero out the chosen edge
    if edge == 'none':
        pass
    elif edge == 'top':
        kernel[0, :] = 0
    elif edge == 'bottom':
        kernel[-1, :] = 0
    elif edge == 'left':
        kernel[:, 0] = 0
    elif edge == 'right':
        kernel[:, -1] = 0
    
    # Normalize after zeroing
    if np.sum(kernel) > 0:  # Avoid division by zero
        kernel = kernel / np.sum(kernel)
    
    return kernel


# Define the process_fiber functions at module level
def process_fiber_with_kernel(fib, skipfib, arm, detsize, kerpara):
    if fib in skipfib:
        return np.zeros(detsize)
    # Create the pattern for glob
    pattern = f'{BASEDIR}/{arm}/'+fnschema.format(arm=arm,fib=fib)
    matching_files = glob(pattern)
    
    # Check if any files were found
    if not matching_files:
        print(f"Warning: No files found matching pattern: {pattern}")
        return np.zeros(detsize)
        
    # Process the first matching file
    kernel = custom_gaussian_kernel(*kerpara)
    d = fits.open(matching_files[0])[0].data
    d = convolve2d(d, kernel, mode='same')
    return d

def process_fiber_without_kernel(fib, skipfib, arm, detsize):
    if fib in skipfib:
        return np.zeros(detsize)
    # Create the pattern for glob
    pattern = f'{BASEDIR}/{arm}/'+fnschema.format(arm=arm,fib=fib)
    matching_files = glob(pattern)
    
    # Check if any files were found
    if not matching_files:
        print(f"Warning: No files found matching pattern: {pattern}")
        return np.zeros(detsize)
        
    # Process the first matching file
    d = fits.open(matching_files[0])[0].data
    return d


########################################################
if __name__ == "__main__":

    if len(sys.argv) < 2:
        print('Usage: python Dkernel.py <arm> [kerpara]')
        print('  <arm>      = Y, J, H, U, B, V, R, or IZ')
        print('  [kerpara]  = (optional) kernel parameters as comma-separated values, e.g., 5,5,3.2,left')
        sys.exit(1)
    arm = sys.argv[1]
    kerpara_given = len(sys.argv) > 2
    kerpara_str = sys.argv[2] if kerpara_given else None

    if arm in ("Y", "J", "H"):
        n_fibers = 75
        detsize = (4096, 4096)
        skipfib = [2, 4, 5, 7, 9, 11, 13, 15, 17, 30, 49, 74]
    else:
        print(f'Error: Unknown arm {arm}')
        sys.exit(1)

    if kerpara_given:
        try:
            kerpara = kerpara_str.split(',')
            if len(kerpara) != 4:
                raise ValueError
            # convert from FWHM to sigma
            kerpara[2] = float(kerpara[2]) / 2.35
            kerpara[0] = int(kerpara[0])
            kerpara[1] = int(kerpara[1])
            # kerpara[3] remains string
        except Exception:
            print('Error: <kerpara> must be four comma-separated values, e.g., 5,5,3.2,left')
            sys.exit(1)

    dsum = np.zeros(detsize)
    with concurrent.futures.ProcessPoolExecutor(max_workers=6) as executor:
        if kerpara_given:
            # Use partial functions to handle the extra parameters
            process_func = functools.partial(process_fiber_with_kernel, arm=arm, detsize=detsize, kerpara=kerpara)
            for d in executor.map(process_func, range(1, n_fibers + 1), [skipfib]*n_fibers):
                dsum += d
        else:
            process_func = functools.partial(process_fiber_without_kernel, arm=arm, detsize=detsize)
            for d in executor.map(process_func, range(1, n_fibers + 1), [skipfib]*n_fibers):
                dsum += d

    # Save the result
    HDUL = fits.PrimaryHDU(dsum)
    if kerpara_given:
        fname = str(project_root.parent / arm / f'{arm}_FP_IFUsum_kern{kerpara[0]}x{kerpara[1]}s{float(kerpara[2])*2.35:.1f}{kerpara[3]}.fits')
    else:
        fname = str(project_root.parent / arm / f'{arm}_FP_IFUsum.fits')
    print(f'Saving to {fname}')
    HDUL.writeto(fname, overwrite=True)
    
