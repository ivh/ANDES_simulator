# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "scipy"]
# ///
"""ANDES SCAO+IFU coupling model, anchored to the PDR numbers.

The PDR package contains no contrast curve or PSF profile -- the AO subsystem
documents (E-AND-AO-*) that hold them are referenced but not delivered.  What it
does give is three numbers, and they are enough to pin a two-component model:

  Strehl              0.6 in H, 0.3 in Y, median seeing
                      (E-AND-AN-PSR-00-00-003_1 Science Cases, S2.4)
  raw contrast        ~1.5e-3 at 25-45 mas, low-piston, NO coronagraph
                      (same, S2.4; the only statement on halo shape in the
                      package -- below 25 mas it is set by the Airy rings)
  requirement         3.0e-3 baseline / 1.0e-3 goal at 20 mas +- 3.5 mas,
                      1000-1700 nm, WITH coronagraph (R-AND-102.0a)

CONTRAST HERE MEANS I(r)/I(peak) -- a radial profile of PSF intensity
normalised to the peak.  That is NOT the same as rho, the fraction of total
stellar flux landing in a spaxel, which is what the simulator needs.  The two
differ by the peak intensity times the spaxel solid angle; conflating them is
worth roughly a factor of a few, and is a plausible source of the factor-10
disagreement in Q4 of PLAN.md.

Construction
------------
Normalised PSF, unit total flux:

    I(r) = S * Airy(r)  +  (1-S) * halo(r)

Airy is the obstructed-aperture diffraction pattern; its peak is A/lambda^2.
halo(r) is flat inside an AO control radius r_c and falls as r^(-11/3)
(Kolmogorov) outside.  r_c is not guessed -- it is SOLVED so that the model
reproduces the documented contrast at the documented separation.  Everything
then follows with flux conserved by construction.
"""
import numpy as np
from scipy.optimize import brentq
from scipy.special import j1

MAS = np.pi / 180 / 3600 / 1000          # rad per mas
D_TEL, COBS = 38.5, 0.28
A_TEL = np.pi / 4 * D_TEL ** 2 * (1 - COBS ** 2)

# PDR anchors
STREHL = {"Y": 0.30, "J": 0.45, "H": 0.60}      # J interpolated; Y/H documented
LAM_UM = {"Y": 1.035, "J": 1.255, "H": 1.635}
CONTRAST_NOCORO = 1.5e-3                         # at 25-45 mas, low piston
R_ANCHOR_MAS = 35.0
CONTRAST_CORO = {"baseline": 3.0e-3, "goal": 1.0e-3}   # at 20 mas, requirement
R_CORO_MAS = 20.0
# SCAO optics transmission in the IFU science path (R-AND-86.2)
T_SCAO = {"nocoro": 0.90, "coro": 0.54}


def airy(r_mas, lam_um):
    """Obstructed-aperture Airy pattern, unit total flux, per steradian."""
    r = np.atleast_1d(r_mas) * MAS
    x = np.pi * D_TEL * r / (lam_um * 1e-6)
    e = COBS
    with np.errstate(invalid="ignore", divide="ignore"):
        t = (2 * j1(x) / x - 2 * e * j1(e * x) / x) / (1 - e ** 2)
    t = np.where(x < 1e-9, 1.0, t)
    peak = A_TEL / (lam_um * 1e-6) ** 2          # sr^-1, unit total flux
    return peak * t ** 2


def halo(r_mas, r_c_mas):
    """Flat to r_c then r^-11/3, unit total flux, per steradian."""
    r = np.atleast_1d(r_mas, )
    # normalisation: int_0^inf h(r) 2 pi r dr = 1, with h=h0 inside r_c
    rc = r_c_mas * MAS
    # inner disk: h0 * pi rc^2 ; outer: h0 * 2 pi rc^2/(11/3-2) = h0*2 pi rc^2/(5/3)
    norm = np.pi * rc ** 2 * (1 + 2 / (5 / 3))
    h0 = 1.0 / norm
    out = np.where(r <= r_c_mas, h0, h0 * (np.maximum(r, 1e-9) / r_c_mas) ** (-11 / 3))
    return out


def contrast_at(band, r_anchor, rc, coro_rejection=1.0):
    """Model I(r_anchor) / I_unocculted(0), the standard contrast metric.

    The peak is deliberately the UNOCCULTED one.  Normalising to the occulted
    peak would let the coronagraph cancel out of its own figure of merit.
    """
    S, lam = STREHL[band], LAM_UM[band]
    ia = coro_rejection * S * airy(r_anchor, lam)[0] + (1 - S) * halo(r_anchor, rc)[0]
    pk = S * airy(0.0, lam)[0] + (1 - S) * halo(0.0, rc)[0]
    return ia / pk


def solve_control_radius(band, contrast, r_anchor, coro_rejection=1.0):
    """Find r_c reproducing the documented I(r_anchor)/I(0).

    C(r_anchor) is NOT monotonic in r_c: a tight halo has fallen off by the
    anchor, a very wide one is diluted, so the curve peaks in between and there
    are two roots.  We take the LARGE-r_c branch -- the small one implies a
    control radius of a few lambda/D, which no ELT SCAO has (M4 gives of order
    40 lambda/D).  If even the peak falls short of the target, the target is
    below the diffraction floor and cannot be met without suppressing the Airy
    pattern, i.e. without a coronagraph; we report that rather than hide it.
    """
    grid = np.logspace(np.log10(5.0), np.log10(5000.0), 400)
    vals = np.array([contrast_at(band, r_anchor, rc, coro_rejection) for rc in grid])
    i_pk = int(np.argmax(vals))
    if vals[-1] > contrast:
        raise ValueError(
            f"{band}: even an infinitely wide halo leaves C({r_anchor:.0f} mas) = "
            f"{vals[-1]:.2e} > target {contrast:.2e} -- the residual diffraction "
            f"floor is above the requirement; more core suppression is needed")
    if vals[i_pk] < contrast:
        raise ValueError(
            f"{band}: max achievable C({r_anchor:.0f} mas) = {vals[i_pk]:.2e} "
            f"< target {contrast:.2e} -- the Airy floor is "
            f"{contrast_at(band, r_anchor, 1e6, coro_rejection):.2e}; "
            "needs core suppression")
    f = lambda rc: contrast_at(band, r_anchor, rc, coro_rejection) - contrast
    return brentq(f, grid[i_pk], grid[-1])


class Coupling:
    """Per-spaxel flux fraction rho for a point source anywhere in the field."""

    def __init__(self, band="J", coronagraph=False, contrast=None, n_sub=61):
        self.band, self.lam, self.S = band, LAM_UM[band], STREHL[band]
        self.coro = coronagraph
        if contrast is None:
            contrast = CONTRAST_CORO["baseline"] if coronagraph else CONTRAST_NOCORO
        self.r_anchor = R_CORO_MAS if coronagraph else R_ANCHOR_MAS
        self.contrast = contrast
        self.coro_rejection = 0.01 if coronagraph else 1.0
        self.r_c = solve_control_radius(band, contrast, self.r_anchor,
                                        self.coro_rejection)
        self.n_sub = n_sub

    IWA_MAS = 20.0        # coronagraph inner working angle (R-AND-102.0a bin)

    def intensity(self, r_mas, occulted=True):
        """PSF of a source, per steradian, unit total flux at the entrance.

        `occulted` says whether the coronagraph acts on THIS source.  It acts on
        the on-axis star; an off-axis companion beyond the IWA passes through.
        Suppressing the planet as well -- which a field-independent rejection
        factor would do -- is wrong and was a bug here.
        """
        rej = self.coro_rejection if (self.coro and occulted) else 1.0
        return (rej * self.S * airy(r_mas, self.lam)
                + (1 - self.S) * halo(r_mas, self.r_c))

    def offaxis_throughput(self, d_mas):
        """Coronagraph throughput for a source d from the axis (0 -> 1 at IWA)."""
        if not self.coro:
            return 1.0
        return float(np.clip((d_mas / self.IWA_MAS) ** 2, 0.0, 1.0))

    def rho(self, spaxels, src_x, src_y, pitch_mas):
        """Fraction of a point source at (src_x,src_y) landing in each spaxel."""
        s = pitch_mas / np.sqrt(3)
        u = np.linspace(-s, s, self.n_sub)
        gx, gy = np.meshgrid(u, u)
        m = _in_hex(gx, gy, s)
        dA = (2 * s / (self.n_sub - 1)) ** 2 * MAS ** 2      # steradian per cell
        d_src = np.hypot(src_x, src_y)
        occ = d_src < self.IWA_MAS
        # the occulted star's suppression is already in coro_rejection; the
        # off-axis ramp applies only to sources outside the IWA
        thr = 1.0 if occ else self.offaxis_throughput(d_src)
        out = {}
        for fib, cx, cy, ring in spaxels:
            rr = np.hypot(gx[m] + cx - src_x, gy[m] + cy - src_y)
            out[fib] = float(np.sum(self.intensity(rr, occulted=occ) * dA) * thr)
        return out


def _in_hex(x, y, s):
    """Pointy-top hexagon, circumradius s: flats at x = +-s*sqrt(3)/2.

    Adjacent centres are s*sqrt(3) apart, so pitch = s*sqrt(3) and the
    inradius (half the flat-to-flat width) is pitch/2.
    """
    x, y = np.abs(x), np.abs(y)
    return (x <= s * np.sqrt(3) / 2) & (np.sqrt(3) * y + x <= np.sqrt(3) * s)


if __name__ == "__main__":
    print("PDR-anchored SCAO+IFU coupling model\n")
    print(f"{'band':>5} {'mode':>9} {'S':>5} {'lam/D':>6} {'anchor':>15} "
          f"{'r_ctrl mas':>11} {'lam/D':>7} {'Airy floor':>11}")
    for band in ("Y", "J", "H"):
        ld = LAM_UM[band] * 1e-6 / D_TEL / MAS
        for coro, tag in ((False, "no coro"), (True, "coro")):
            try:
                c = Coupling(band, coronagraph=coro)
                fl = contrast_at(band, c.r_anchor, 1e6, c.coro_rejection)
                print(f"{band:>5} {tag:>9} {c.S:5.2f} {ld:6.2f} "
                      f"{c.contrast:.1e}@{c.r_anchor:.0f}mas {c.r_c:11.1f} "
                      f"{c.r_c/ld:7.1f} {fl:11.2e}")
            except ValueError as e:
                print(f"{band:>5} {tag:>9}  UNREACHABLE: {e}")

    print("\nflux conservation check (integral of I over the field):")
    for band in ("Y", "J", "H"):
        c = Coupling(band)
        r = np.logspace(-3, 5, 200000)
        tot = np.trapezoid(c.intensity(r) * 2 * np.pi * r * MAS ** 2, r)
        print(f"  {band}: {tot:.4f}")


def hex_ring(n):
    if n == 0:
        return [(0, 0)]
    dirs = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
    q, r, out = n, 0, []
    for d in range(6):
        for _ in range(n):
            out.append((q, r))
            q += dirs[(d + 2) % 6][0]
            r += dirs[(d + 2) % 6][1]
    return out


def build_ifu(pitch):
    s = pitch / np.sqrt(3)
    out, fib = [], 0
    for ring in range(5):
        for (q, r) in hex_ring(ring):
            fib += 1
            out.append((fib, s * np.sqrt(3) * (q + r / 2.0), s * 1.5 * r, ring))
    return out


def report(pitch=16.0, sep=37.3, pa=35.0, albedo=0.3):
    cp = albedo * (1.07 * 6.371e6 / (0.04848 * 1.496e11)) ** 2 / np.pi
    sp = build_ifu(pitch)
    px, py = sep * np.cos(np.radians(pa)), sep * np.sin(np.radians(pa))
    pfib = min(sp, key=lambda t: np.hypot(t[1] - px, t[2] - py))[0]
    print(f"\nProxima b: sep {sep} mas, PA {pa} deg, {pitch:.0f} mas spaxels, "
          f"planet/star contrast {cp:.2e}")
    print(f"{'band':>5} {'mode':>8} {'rho_halo':>10} {'rho_planet':>11} "
          f"{'planet/halo':>12} {'T_scao':>7}")
    res = {}
    for band in ("Y", "J", "H"):
        for coro, tag in ((False, "no coro"), (True, "coro")):
            c = Coupling(band, coronagraph=coro)
            rh = c.rho(sp, 0.0, 0.0, pitch)[pfib]
            rp = c.rho(sp, px, py, pitch)[pfib]
            t = T_SCAO["coro" if coro else "nocoro"]
            print(f"{band:>5} {tag:>8} {rh:10.3e} {rp:11.4f} "
                  f"{cp*rp/rh:12.3e} {t:7.2f}")
            res[(band, tag)] = (rh, rp, cp * rp / rh)
    return res


if __name__ == "__main__":
    report()
