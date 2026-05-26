# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "astropy", "scipy", "matplotlib"]
# ///
"""Measure the inter-fiber pitch in H-band using two flat-field exposures
(fibers 10 and 20) and estimate the cross-talk of the fiber profile at
1..4 fiber distances.

Background is estimated as the mean of the bin-sums at the +/-4 pitch
positions (converted to a per-pixel level) and subtracted from each
profile. The Gaussian is then re-fit with the baseline fixed to zero."""

from pathlib import Path
import numpy as np
from astropy.io import fits
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent.parent
F10 = HERE / "H_FF_fib10_wl1610-1620.fits"
F20 = HERE / "H_FF_fib20_wl1610-1620.fits"
COL_LO, COL_HI = 2040, 2140
FIBER_DELTA = 10  # fibers between the two exposures


def gauss(y, amp, mu, sigma, base):
    return base + amp * np.exp(-0.5 * ((y - mu) / sigma) ** 2)


def gauss_nobg(y, amp, mu, sigma):
    return amp * np.exp(-0.5 * ((y - mu) / sigma) ** 2)


def median_profile(path):
    data = fits.getdata(path).astype(float)
    return np.median(data[:, COL_LO:COL_HI], axis=1)


def fit_peak(profile, fixed_base=None):
    y = np.arange(profile.size)
    i_peak = int(np.argmax(profile))
    half = 15
    lo, hi = max(0, i_peak - half), min(profile.size, i_peak + half + 1)
    ys, ps = y[lo:hi], profile[lo:hi]
    if fixed_base is None:
        base0 = float(np.median(profile))
        amp0 = float(ps.max() - base0)
        p0 = (amp0, float(i_peak), 2.0, base0)
        popt, _ = curve_fit(gauss, ys, ps, p0=p0)
        return popt  # amp, mu, sigma, base
    amp0 = float(ps.max())
    p0 = (amp0, float(i_peak), 2.0)
    popt, _ = curve_fit(gauss_nobg, ys, ps, p0=p0)
    amp, mu, sigma = popt
    return amp, mu, sigma, 0.0


def bin_sum(profile, center, width):
    """Sum profile flux in [center-width/2, center+width/2], with
    fractional weighting on the edge pixels."""
    lo = center - width / 2.0
    hi = center + width / 2.0
    i0 = int(np.floor(lo))
    i1 = int(np.floor(hi))
    total = 0.0
    for i in range(i0, i1 + 1):
        if i < 0 or i >= profile.size:
            continue
        left = max(lo, i)
        right = min(hi, i + 1)
        w = max(0.0, right - left)
        total += w * profile[i]
    return total


def main():
    p10_raw = median_profile(F10)
    p20_raw = median_profile(F20)

    # First pass: find peak positions and pitch with floating baseline.
    a10_0, mu10_0, s10_0, b10_0 = fit_peak(p10_raw)
    a20_0, mu20_0, s20_0, b20_0 = fit_peak(p20_raw)
    pitch = abs(mu20_0 - mu10_0) / FIBER_DELTA

    # Background = mean of +/-4-pitch bin sums, converted to per-pixel level.
    def bg_level(prof, mu):
        bp = bin_sum(prof, mu + 4 * pitch, pitch)
        bm = bin_sum(prof, mu - 4 * pitch, pitch)
        return 0.5 * (bp + bm) / pitch

    bg10 = bg_level(p10_raw, mu10_0)
    bg20 = bg_level(p20_raw, mu20_0)
    p10 = p10_raw - bg10
    p20 = p20_raw - bg20

    # Re-fit with baseline fixed to zero.
    a10, mu10, s10, _ = fit_peak(p10, fixed_base=0.0)
    a20, mu20, s20, _ = fit_peak(p20, fixed_base=0.0)
    pitch = abs(mu20 - mu10) / FIBER_DELTA
    sigma_mean = 0.5 * (s10 + s20)
    fwhm = 2.0 * np.sqrt(2.0 * np.log(2.0)) * sigma_mean

    print(f"bg level (per pix): fiber 10 = {bg10:.4g}, fiber 20 = {bg20:.4g}")
    print(f"fiber 10: mu={mu10:.3f} px, sigma={s10:.3f} px, amp={a10:.3g}")
    print(f"fiber 20: mu={mu20:.3f} px, sigma={s20:.3f} px, amp={a20:.3g}")
    print(f"per-fiber pitch: {pitch:.3f} px")
    print(f"mean sigma: {sigma_mean:.3f} px  (FWHM = {fwhm:.3f} px)")
    print()
    print(f"Pixel-flux cross-talk: bins of width = pitch ({pitch:.3f} px), "
          "bg-subtracted profiles, baseline fixed at 0.")
    for label, profile, mu in [
        ("fiber 10", p10, mu10),
        ("fiber 20", p20, mu20),
    ]:
        peak_bin = bin_sum(profile, mu, pitch)
        print(f"\n  {label}  (peak-bin sum = {peak_bin:.4g})")
        print(f"  {'N':>2}  {'frac (+side)':>14}  {'frac (-side)':>14}  "
              f"{'frac (mean)':>14}")
        for n in range(1, 4):
            plus = bin_sum(profile, mu + n * pitch, pitch) / peak_bin
            minus = bin_sum(profile, mu - n * pitch, pitch) / peak_bin
            mean = 0.5 * (plus + minus)
            print(f"  {n:>2}  {plus:>14.4e}  {minus:>14.4e}  {mean:>14.4e}")

    fig, ax = plt.subplots(figsize=(9, 5.5))
    y_full = np.arange(p10.size)

    for profile, mu, s, a, col, lab in [
        (p10, mu10, s10, a10, "C0", "fiber 10"),
        (p20, mu20, s20, a20, "C3", "fiber 20"),
    ]:
        rel = y_full - mu
        m = np.abs(rel) <= 12
        ax.plot(rel[m], profile[m], color=col, lw=1.0,
                label=f"{lab} median (bg sub)", drawstyle="steps-mid")
        yy = np.linspace(-12, 12, 600)
        ax.plot(yy, gauss_nobg(yy + mu, a, mu, s),
                color=col, lw=0.8, ls="--", label=f"{lab} gauss fit (bg=0)")

    for n in range(-4, 5):
        ax.axvline((n + 0.5) * pitch, color="grey", ls=":", lw=0.6)
    ymax = max(p10.max(), p20.max()) * 1.1
    for n in range(-4, 5):
        ax.text(n * pitch, ymax, f"{n:+d}" if n else "0",
                color="grey", fontsize=8, ha="center", va="bottom")
    ax.set_xlabel("y pixel relative to peak")
    ax.set_ylabel(f"median over cols {COL_LO}-{COL_HI} (bg subtracted)")
    ax.set_title(
        f"H-band fiber profile  |  sigma~{sigma_mean:.2f} px, "
        f"pitch~{pitch:.2f} px  (dotted = bin borders, ticks = bin index)"
    )
    ax.legend(loc="upper right", fontsize=8)

    out = HERE / "H_fiber_profile.png"
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    print(f"\nplot saved -> {out}")


if __name__ == "__main__":
    main()
