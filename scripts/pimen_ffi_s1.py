"""Rediscover Pi Mensae c from 30-minute FFI pixels instead of the ready-made light curve.

Same star, same sector, but this time we build the light curve ourselves:
download a 25x25 pixel cutout, pick the star's pixels, sum them, regress out
the background, then run the five moves. Compare against pimen_s1.py.
Run from the project root:  python scripts/pimen_ffi_s1.py
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import lightkurve as lk

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tesshunt import pipeline as P

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PLOTS = os.path.join(ROOT, "plots")
DATA = os.path.join(ROOT, "data")
os.makedirs(PLOTS, exist_ok=True)
os.makedirs(DATA, exist_ok=True)

STAR = "Pi Men"
STAR_RADIUS_SUNS = 1.10
STAR_DENSITY_SUNS = 0.85

# 0. Pixels. TESScut carves a small square out of every FFI in the sector.
tpf = lk.search_tesscut(STAR, sector=1).download(cutout_size=25, download_dir=DATA)
print(f"cutout: {tpf.flux.shape[0]} frames of {tpf.flux.shape[1]}x{tpf.flux.shape[2]} pixels, "
      f"cadence {np.nanmedian(np.diff(tpf.time.value))*24*60:.0f} min")

# 0b. Pixels -> light curve.
lc, raw, aper = P.lightcurve_from_cutout(tpf, threshold=10, n_pca=5)
print(f"aperture: {aper.sum()} pixels")
tpf.plot(aperture_mask=aper, mask_color="red")
plt.savefig(os.path.join(PLOTS, "ffi_0_aperture.png"), dpi=120)
ax = raw.normalize().plot(label="raw sum")
lc.normalize().plot(ax=ax, label="background regressed")
ax.set_title("FFI: pixels summed, then background removed")
plt.savefig(os.path.join(PLOTS, "ffi_1_raw_vs_corrected.png"), dpi=120)

# 1-2. clean, flatten
lc = P.clean(lc)
flat = P.flatten(lc, window_hours=33)
n_before = len(flat)
flat = P.mask_jitter(flat, window_hours=6, factor=3)
print(f"{len(flat)} good cadences over {flat.time.max().value - flat.time.min().value:.1f} days "
      f"({n_before - len(flat)} dropped as jitter)")
ax = flat.plot(); ax.set_title("FFI: flattened")
plt.savefig(os.path.join(PLOTS, "ffi_2_flat.png"), dpi=120)

# 3. search
r = P.search(flat, min_period=1.0, max_period=13.0)
print(f"best period {r.period:.4f} d, duration {r.duration*24:.2f} h, depth {r.depth*1e6:.0f} ppm, SNR {r.snr:.1f}")
ax = r.periodogram.plot(); ax.set_title("FFI: periodogram")
plt.savefig(os.path.join(PLOTS, "ffi_3_periodogram.png"), dpi=120)

# 4. fold
folded = P.fold(flat, r.period, r.t0)
ax = folded.scatter(s=4, alpha=0.4)
folded.bin(time_bin_size=0.02).plot(ax=ax, color="red", lw=2)
ax.set_xlim(-0.3, 0.3); ax.set_title(f"FFI: folded on {r.period:.3f} d")
plt.savefig(os.path.join(PLOTS, "ffi_4_folded.png"), dpi=120)

# 5. vet + physics
v = P.vet(flat, r.period, r.t0, r.duration)
print(f"transits stacked: {v['n_transits']}")
print(f"odd {v['odd_depth_ppm']:.0f} ppm, even {v['even_depth_ppm']:.0f} ppm, +/- {v['depth_err_ppm']:.0f} ppm each")
print(f"secondary {v['secondary_depth_ppm']:.0f} ppm (want ~0)")
print(f"planet radius ~{P.planet_radius_earths(r.depth, STAR_RADIUS_SUNS):.2f} Earths (published 2.04)")
print(f"duration {r.duration*24:.2f} h vs max allowed {P.max_duration_hours(r.period, STAR_DENSITY_SUNS):.2f} h")
print(f"plots written to {PLOTS}")
