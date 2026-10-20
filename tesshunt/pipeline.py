"""The five moves of a transit search, as plain functions.

    clean   -> drop bad cadences, scale so "normal" = 1.0
    flatten -> remove slow drift so only fast dips remain
    search  -> try every period, score how well dips stack (BLS)
    fold    -> wrap time onto the winning period so dips overlap
    vet     -> three cheap falsification tests

Each function takes a lightkurve LightCurve and returns one, except
search (returns a dict of numbers) and vet (returns a dict of numbers).
"""
from dataclasses import dataclass
import numpy as np
import lightkurve as lk

# ---------- 1. clean ----------

def clean(lc, sigma_clip=5.0):
    """Drop NaNs and quality-flagged cadences, scale to 1.0, clip wild outliers.

    Upward outliers (cosmic rays, flares) are clipped. Downward ones are kept,
    because a transit is a downward outlier and we must not throw it away.
    """
    lc = lc.remove_nans().normalize()
    med = np.nanmedian(lc.flux.value)
    mad = np.nanmedian(np.abs(lc.flux.value - med)) * 1.4826  # robust sigma
    keep = lc.flux.value < med + sigma_clip * mad
    return lc[keep]

def mask_jitter(flat, window_hours=6.0, factor=3.0):
    """Drop stretches where the scatter is far above normal (bad pointing, glitches).

    FFI cutouts carry no quality flags, so we make our own. We look at the
    size of the jump between consecutive points. Inside a transit the points
    still sit close to each other (the dip is smooth; only entry and exit are
    steps), so a rolling MEDIAN of the jumps barely notices it. During bad
    pointing every consecutive pair jumps by far more than the noise, so the
    rolling median shoots up. Drop cadences where it exceeds `factor` x normal.

    Do NOT use rolling standard deviation here: for a bright, quiet star a
    transit raises local scatter several-fold and gets deleted (we hit this).
    Run this AFTER flatten.
    """
    f = flat.flux.value
    jump = np.abs(np.diff(f, prepend=f[0]))
    cadence_days = np.nanmedian(np.diff(flat.time.value))
    w = max(int(round(window_hours / 24.0 / cadence_days)), 3)
    half = w // 2
    local = np.full(len(f), np.nan)
    for i in range(len(f)):
        local[i] = np.nanmedian(jump[max(0, i - half): i + half + 1])
    typical = np.nanmedian(local)
    keep = local < factor * typical
    return flat[keep]

# ---------- 2. flatten ----------

def flatten(lc, window_hours=33.0):
    """Remove everything slower than window_hours. Transits (a few hours) pass through.

    The window must be several times longer than the longest transit you want
    to keep, or the smoother will follow the dip down and erase it.
    """
    cadence_days = np.nanmedian(np.diff(lc.time.value))
    window = int(round(window_hours / 24.0 / cadence_days))
    if window % 2 == 0:
        window += 1
    return lc.flatten(window_length=max(window, 5))

# ---------- 3. search ----------

@dataclass
class SearchResult:
    period: float          # days
    t0: float              # time of a transit center, same units as lc.time
    duration: float        # days
    depth: float           # fraction (0.00025 = 250 ppm)
    power: float           # BLS power at the peak
    snr: float             # peak power / robust scatter of the periodogram
    periodogram: object    # the lightkurve periodogram, for plotting

def search(flat, min_period=0.5, max_period=None, n_periods=20000,
           min_duration_hours=0.5, max_duration_hours=8.0, duration_step_hours=0.25):
    """Box least squares over a dense grid of periods and durations.

    max_period defaults to half the baseline, so at least two transits fit.
    snr is the peak height in units of the periodogram's own scatter; values
    above ~7-9 are conventionally worth a human look.
    """
    baseline = flat.time.max().value - flat.time.min().value
    if max_period is None:
        max_period = baseline / 2.0
    periods = np.linspace(min_period, max_period, n_periods)
    durations = np.arange(min_duration_hours, max_duration_hours + 1e-9, duration_step_hours) / 24.0
    pg = flat.to_periodogram(method="bls", period=periods, duration=durations)
    power = pg.power.value
    med = np.nanmedian(power)
    mad = np.nanmedian(np.abs(power - med)) * 1.4826
    depth = pg.depth_at_max_power
    depth = float(depth.value) if hasattr(depth, "value") else float(depth)
    return SearchResult(
        period=float(pg.period_at_max_power.value),
        t0=float(pg.transit_time_at_max_power.value),
        duration=float(pg.duration_at_max_power.value),
        depth=depth,
        power=float(pg.max_power.value) if hasattr(pg.max_power, "value") else float(pg.max_power),
        snr=float((np.nanmax(power) - med) / mad) if mad > 0 else float("nan"),
        periodogram=pg,
    )

# ---------- 4. fold ----------

def fold(flat, period, t0):
    """Wrap time onto the period so every transit lands at phase 0."""
    return flat.fold(period=period, epoch_time=t0)

# ---------- 5. vet ----------

def vet(flat, period, t0, duration):
    """Three cheap tests. All are in ppm so you can compare by eye.

    odd_depth, even_depth : should agree within noise. Mismatch -> eclipsing binary at 2x period.
    secondary_depth       : dip half a period away. Nonzero -> a star, not a planet.
    depth_err             : rough 1-sigma on a single depth, for judging "agree".
    """
    t = flat.time.value
    f = flat.flux.value
    from_center = ((t - t0 + period / 2) % period) - period / 2
    in_tr = np.abs(from_center) < duration / 2
    cyc = np.round((t - t0) / period).astype(int)
    odd = f[in_tr & (cyc % 2 == 1)]
    even = f[in_tr & (cyc % 2 == 0)]
    out = f[~in_tr]
    scatter = np.nanstd(out)
    from_half = ((t - t0) % period) - period / 2
    sec = f[np.abs(from_half) < duration / 2]
    n_in = max(in_tr.sum(), 1)
    return {
        "odd_depth_ppm": (1 - np.nanmean(odd)) * 1e6 if len(odd) else float("nan"),
        "even_depth_ppm": (1 - np.nanmean(even)) * 1e6 if len(even) else float("nan"),
        "secondary_depth_ppm": (1 - np.nanmean(sec)) * 1e6 if len(sec) else float("nan"),
        "depth_err_ppm": scatter / np.sqrt(n_in / 2) * 1e6,
        "n_transits": int(len(np.unique(cyc[in_tr]))),
    }

# ---------- physics helpers ----------

EARTH_RADII_PER_SUN = 109.1

def planet_radius_earths(depth, star_radius_suns):
    """Depth is an area ratio, so radius ratio is its square root."""
    return np.sqrt(depth) * star_radius_suns * EARTH_RADII_PER_SUN

def max_duration_hours(period_days, star_density_suns=1.0):
    """Longest possible transit (dead-center crossing, circular orbit).

    A dip much longer than this is a red flag: the star is bigger than the
    catalog says, or the dip belongs to a different, puffier star in the pixel.
    """
    return 13.0 * (period_days / 365.25) ** (1 / 3) * star_density_suns ** (-1 / 3)

# ---------- FFI pixels -> light curve ----------

def lightcurve_from_cutout(tpf, threshold=10, n_pca=5):
    """Build a light curve from a TESScut pixel cube.

    1. aperture: pixels brighter than `threshold` sigma above the cutout median,
       connected to the center pixel. Those are the star.
    2. raw light curve: sum the aperture pixels each cadence.
    3. background: the non-aperture pixels see the same scattered light and
       drift as the star. Take their top n_pca patterns and regress them out.
    Returns (corrected_lc, raw_lc, aperture_mask).
    """
    aper = tpf.create_threshold_mask(threshold=threshold, reference_pixel="center")
    raw = tpf.to_lightcurve(aperture_mask=aper)
    regressors = tpf.flux[:, ~aper]
    regressors = regressors.value if hasattr(regressors, "value") else regressors
    dm = lk.DesignMatrix(regressors, name="background").pca(n_pca).append_constant()
    corrector = lk.RegressionCorrector(raw.remove_nans())
    corrected = corrector.correct(dm)
    return corrected, raw, aper
