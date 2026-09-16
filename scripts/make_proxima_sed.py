# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy"]
# ///
"""Build SED/proxima.csv: a physically calibrated Proxima Centauri spectrum.

SED/phoenix.csv is already the right model -- PHOENIX Teff=3000 / logg=5.0 /
[M/H]=0 (Husser+ 2013), the same one Bugatti et al. 2025 used for Proxima -- but
it carries arbitrary flux units.  This script normalises it to **photons per
second per wavelength sample, detected, at 100% coupling into one fibre**, so
that `andes-sim --exposure T --flux RHO` yields physically meaningful counts:

    N_i = F_lam(lam_i) * dlam_i * lam_i/(hc) * A_tel * T_optics(lam_i) * T_tell(lam_i)

Units matter here.  The default flux_unit is "ph/s", under which PyEchelle's
CSVSource returns the rows *verbatim* as discrete lines carrying flux*exptime
photons each -- it does NOT integrate a density.  So the per-sample bin width
dlam_i is folded in below, and the file is photons-per-sample.

Pass --flux-unit ph/s/AA to use the density branch instead (which interpolates
onto PyEchelle's own per-order grid and is therefore insensitive to the sampling
of this file); the file would then have to be written without the dlam factor.
Getting this backwards costs exactly the sampling interval in Angstrom -- 50x
for this 0.02 A grid.

T_optics deliberately EXCLUDES the blaze: PyEchelle applies that itself.
Verified -- ZEMAX.get_efficiency() returns SystemEfficiency[GratingEfficiency],
i.e. grating blaze only, with no coatings and no detector QE.
"""
import numpy as np
from pathlib import Path

SRC = Path(__file__).resolve().parent.parent
SED = SRC / "SED"
ETC = Path.home() / "ANDES/ANDES_ETC_April2025"

h, c = 6.62607015e-34, 2.99792458e8
R_SUN, PC = 6.957e8, 3.0856775814913673e16

# Proxima Centauri (Bugatti et al. 2025, Table 1)
R_STAR = 0.141 * R_SUN
DIST = 1.3012 * PC
# ELT (ANDES ETC April 2025)
D_TEL, COBS = 38.5, 0.28
A_TEL = np.pi / 4 * D_TEL**2 * (1 - COBS**2)          # 1072.9 m^2

WL_MIN, WL_MAX = 940.0, 1810.0                         # nm, YJH + margin
AIRMASS_REF = 1.0

print(f"A_tel = {A_TEL:.1f} m^2")

# --- stellar model -----------------------------------------------------------
d = np.loadtxt(SED / "phoenix.csv", delimiter=",", comments="#")
m = (d[:, 0] >= WL_MIN) & (d[:, 0] <= WL_MAX)
wl = d[m, 0]                                           # nm
# PHOENIX surface flux, erg/s/cm2/cm -> W/m2/um is x1e-7; then dilute by (R/d)^2
flam = d[m, 1] * 1e-7 * (R_STAR / DIST) ** 2           # W/m2/um at Earth
print(f"stellar grid: {len(wl)} samples, {wl[0]:.1f}-{wl[-1]:.1f} nm")

# per-sample bin width (trapezoidal midpoints), in um
edges = np.empty(len(wl) + 1)
edges[1:-1] = 0.5 * (wl[1:] + wl[:-1])
edges[0] = wl[0] - 0.5 * (wl[1] - wl[0])
edges[-1] = wl[-1] + 0.5 * (wl[-1] - wl[-2])
dlam = np.diff(edges) / 1000.0                         # um

# --- instrument throughput, WITHOUT the blaze --------------------------------
e = np.loadtxt(ETC / "Efficiencies/efficiencies_interpol_YJH_Apr2025.dat")
# cols: wl_nm, telescope, fibre-link, front-end, back-end, detector, blaze
t_opt_tab = e[:, 1] * e[:, 2] * e[:, 3] * e[:, 4] * e[:, 5]
t_opt = np.interp(wl, e[:, 0], t_opt_tab)
print(f"optics (no blaze): {t_opt.min():.4f} .. {t_opt.max():.4f}, mean {t_opt.mean():.4f}")

# --- telluric transmission ---------------------------------------------------
tt = np.loadtxt(SED / "sky_transmission_YK.csv", delimiter=",")
t_tell = np.interp(wl, tt[:, 0], tt[:, 1], left=1.0, right=1.0)
t_tell = np.clip(t_tell, 0.0, 1.0)
print(f"telluric  (airmass {AIRMASS_REF:.1f}): median {np.median(t_tell):.3f}")

# --- photons -----------------------------------------------------------------
# photons/s in each wavelength sample's bin (flux_unit "ph/s", the default)
nph = flam * dlam * (wl * 1e-9) / (h * c) * A_TEL * t_opt * t_tell

out = SED / "proxima.csv"
with open(out, "w") as f:
    f.write("# Proxima Centauri, physically calibrated for andes-sim\n")
    f.write("# PHOENIX Teff=3000 logg=5.0 [M/H]=0.0 (Husser+ 2013), same model as Bugatti+2025\n")
    f.write(f"# normalisation: (R*/d)^2 with R*={R_STAR/R_SUN:.3f} Rsun, d={DIST/PC:.4f} pc\n")
    f.write(f"# x A_tel={A_TEL:.1f} m2 (D={D_TEL} m, COBS={COBS})\n")
    f.write("# x optics WITHOUT blaze (ETC Apr2025 cols tel*FL*FE*BE*det) -- PyEchelle applies the blaze\n")
    f.write(f"# x telluric transmission at airmass {AIRMASS_REF:.1f} (SED/sky_transmission_YK.csv)\n")
    f.write("# column 2 = photons/s in that sample's bin, detected, 100% fibre coupling\n")
    f.write("# matches the default flux_unit 'ph/s' (photons per sample, not a density)\n")
    f.write("# for a density instead, divide by the bin width and pass --flux-unit ph/s/AA\n")
    f.write("# multiply by the per-spaxel coupling rho via andes-sim --flux\n")
    f.write("# scaling: 1.0\n")
    for a, b in zip(wl, nph):
        f.write(f"{a:.4f},{b:.6e}\n")

print(f"\nwrote {out}  ({len(wl)} rows, ph/s per sample)")
print(f"  total {nph.sum():.3e} ph/s over {WL_MIN:.0f}-{WL_MAX:.0f} nm")
for lo, hi, nm in ((960, 1110, "Y"), (1160, 1350, "J"), (1470, 1800, "H")):
    k = (wl >= lo) & (wl <= hi)
    print(f"  {nm}: {nph[k].sum():.3e} ph/s")
