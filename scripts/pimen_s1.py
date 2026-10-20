"""Rediscover Pi Mensae c in TESS Sector 1 two-minute data.

Five moves: clean, flatten, search, fold, vet. Writes three plots to plots/.
Run from the project root:  python scripts/pimen_s1.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import lightkurve as lk

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PLOTS = os.path.join(ROOT, "plots")
DATA = os.path.join(ROOT, "data")
os.makedirs(PLOTS, exist_ok=True)
os.makedirs(DATA, exist_ok=True)

# 1. Clean: download the SPOC light curve, drop cadences the spacecraft flagged as bad.
sr = lk.search_lightcurve("Pi Men", mission="TESS", sector=1, author="SPOC")
lc = sr.download(quality_bitmask="default", download_dir=DATA)
lc = lc.remove_nans().normalize()
span = lc.time.max().value - lc.time.min().value
print(f"{len(lc)} good cadences over {span:.1f} days")

# 2. Flatten: remove anything slower than ~33 hours (1001 cadences x 2 min).
flat = lc.flatten(window_length=1001)
ax = flat.plot()
ax.set_title("Step 2: flattened")
plt.savefig(os.path.join(PLOTS, "1_flat.png"), dpi=120)

# 3. Search: box least squares over 1 to 13 day periods.
periods = np.linspace(1, 13, 20000)
pg = flat.to_periodogram(method="bls", period=periods, duration=[0.05, 0.1, 0.15, 0.2])
P = pg.period_at_max_power.value
t0 = pg.transit_time_at_max_power.value
dur = pg.duration_at_max_power.value
depth = pg.depth_at_max_power
depth = depth.value if hasattr(depth, "value") else depth
print(f"best period {P:.4f} d, epoch {t0:.3f}, duration {dur*24:.1f} h, depth {depth*1e6:.0f} ppm")
ax = pg.plot()
ax.set_title("Step 3: periodogram")
plt.savefig(os.path.join(PLOTS, "2_periodogram.png"), dpi=120)

# 4. Fold on the winning period; bin to show the dip shape.
folded = flat.fold(period=P, epoch_time=t0)
ax = folded.scatter(s=2, alpha=0.3)
folded.bin(time_bin_size=0.005).plot(ax=ax, color="red", lw=2)
ax.set_xlim(-0.3, 0.3)
ax.set_title(f"Step 4: folded on {P:.3f} d")
plt.savefig(os.path.join(PLOTS, "3_folded.png"), dpi=120)

# 5. Vet: odd vs even depths should match; nothing at phase 0.5.
t = flat.time.value
f = flat.flux.value
phase_t = ((t - t0 + P / 2) % P) - P / 2          # time from nearest transit center
cycle = np.round((t - t0) / P).astype(int)
in_transit = np.abs(phase_t) < dur / 2
odd = f[in_transit & (cycle % 2 == 1)]
even = f[in_transit & (cycle % 2 == 0)]
print(f"odd depth {(1-odd.mean())*1e6:.0f} ppm, even depth {(1-even.mean())*1e6:.0f} ppm (should match)")
phase_sec = ((t - t0) % P) - P / 2                # time from the half-way point
near_sec = np.abs(phase_sec) < dur / 2
print(f"secondary-eclipse depth {(1-f[near_sec].mean())*1e6:.0f} ppm (should be ~0)")
print(f"plots written to {PLOTS}")
