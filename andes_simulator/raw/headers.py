"""ESO header stamping and MEF assembly for raw frames.

DRL spec v1.2 Ch. 4.1: one FITS file per spectrograph, one image extension
per detector arm, ESO-standard classification headers on the primary HDU.
The keyword set matches the edps workflow's expectations (the calibration
plan YAML is the shared spec; ~/ANDES/edps/tests/make_test_data.py is the
header-only counterpart).
"""

from datetime import datetime, timezone
from typing import Dict, Optional

import numpy as np
from astropy.io import fits

MJD_EPOCH = datetime(1858, 11, 17, tzinfo=timezone.utc)


def long_key(key: str) -> str:
    """Dotted lowercase keyword to FITS form ('dpr.type' -> HIERARCH ESO ...)."""
    return "HIERARCH ESO " + key.upper().replace(".", " ") if "." in key else key.upper()


def mjd_of(dt: datetime) -> float:
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    return (dt - MJD_EPOCH).total_seconds() / 86400.0


def build_raw_hdul(arm: str,
                   images: Dict[str, np.ndarray],
                   ext_meta: Dict[str, Dict],
                   obs_time: datetime,
                   exptime: float,
                   keywords: Dict[str, object],
                   origin_note: Optional[str] = None) -> fits.HDUList:
    """Assemble the MEF for one exposure of one spectrograph arm.

    images: band -> uint16 frame, in the arm's band order.
    ext_meta: band -> {'gain', 'ron', 'readout', 'simulated'} for the
    extension headers. keywords: dotted lowercase ESO keys for the primary.
    """
    primary = fits.PrimaryHDU()
    hdr = primary.header
    hdr['INSTRUME'] = ('ANDES', 'Instrument name')
    hdr['TELESCOP'] = ('ESO-ELT', 'Telescope name')
    hdr['ORIGIN'] = ('E2E-SIM', origin_note or 'ANDES E2E simulator raw frame')
    if obs_time.tzinfo is None:
        obs_time = obs_time.replace(tzinfo=timezone.utc)
    hdr['DATE-OBS'] = (obs_time.strftime('%Y-%m-%dT%H:%M:%S.%f')[:-3],
                       'Observation start')
    hdr['MJD-OBS'] = (mjd_of(obs_time), 'Observation start (MJD)')
    hdr['EXPTIME'] = (exptime, '[s] Exposure time')
    hdr[long_key('seq.arm')] = arm
    for key, value in keywords.items():
        if value is not None:
            hdr[long_key(key)] = value

    hdulist = fits.HDUList([primary])
    for band, image in images.items():
        ext = fits.ImageHDU(data=image, name=band)
        meta = ext_meta.get(band, {})
        ext.header['BUNIT'] = 'ADU'
        if 'gain' in meta:
            ext.header[long_key('det.chip.gain')] = (
                meta['gain'], '[e-/ADU] placeholder value')
        if 'ron' in meta:
            ext.header[long_key('det.chip.ron')] = (
                meta['ron'], '[e-] placeholder value')
        ext.header[long_key('sim.simulated')] = (
            bool(meta.get('simulated', True)),
            'F: detector-noise-only extension')
        hdulist.append(ext)
    return hdulist
