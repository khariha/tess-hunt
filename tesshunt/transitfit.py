"""Physical transit fitting: gather near-transit data from every sector, then
least-squares and MCMC with batman. Used by scripts/fit_transit.py and
scripts/mcmc_fit.py.
"""
import numpy as np
import batman
from tess_stars2px import tess_stars2px_function_entry
from tesshunt import pipeline as P
from tesshunt.sweep_core import get_tpf, _central_aperture

G_SI, RSUN, MSUN, REARTH, RJUP, DAY, AU = 6.674e-11, 6.957e8, 1.989e30, 6.371e6, 7.1492e7, 86400.0, 1.496e11
RHO_SUN = 1408.0

def gather(tic, ra, dec, P0, t00, dur0, window=3.0, cutout=15, log=print, sectors=None):
    """Return dict of arrays (t, f, e, exp, sector) for points within `window` durations of a transit."""
    if sectors is None:
        sectors = sorted(set(int(s) for s in tess_stars2px_function_entry(tic, ra, dec)[3]))
    T, F, E, X, S = [], [], [], [], []
    for s in sectors:
        try:
            tpf = get_tpf(ra, dec, s, cutout, "s3")
            if tpf is None:
                continue
            try:
                lc, raw, aper = P.lightcurve_from_cutout(tpf, threshold=5)
                if aper.sum() == 0 or aper.sum() > 0.5 * aper.size:
                    raise ValueError
            except Exception:
                lc, raw, aper = _central_aperture(tpf)
            flat = P.mask_jitter(P.flatten(P.clean(lc), 33), 6, 3)
            t, f = flat.time.value, flat.flux.value
            exp = float(np.nanmedian(np.diff(t)))
            ph = ((t - t00 + P0 / 2) % P0) - P0 / 2
            keep = np.abs(ph) < window * dur0
            oot = (np.abs(ph) > dur0) & keep
            sig = np.nanstd(f[oot]) if oot.sum() > 10 else np.nanstd(f)
            T.append(t[keep]); F.append(f[keep]); E.append(np.full(keep.sum(), sig)); X.append(np.full(keep.sum(), exp)); S.append(np.full(keep.sum(), s))
            log(f"sector {s}: {keep.sum()} points near transit, {exp*1440:.0f}-min cadence, scatter {sig*1e6:.0f} ppm")
        except Exception as e:
            log(f"sector {s}: {type(e).__name__}: {e}")
    return dict(t=np.concatenate(T), f=np.concatenate(F), e=np.concatenate(E), exp=np.concatenate(X), sector=np.concatenate(S), sectors=sectors)

class Model:
    """batman light curve with per-cadence super-sampling, cached per exposure group."""
    def __init__(self, t, exp):
        self.groups = []
        for ex in np.unique(exp):
            m = exp == ex
            pr = batman.TransitParams()
            pr.t0, pr.per, pr.rp, pr.a, pr.inc, pr.ecc, pr.w = 0.0, 1.0, 0.1, 10.0, 90.0, 0.0, 90.0
            pr.u, pr.limb_dark = [0.3, 0.3], "quadratic"
            ss = max(int(round(ex * 1440 / 2)), 1)
            self.groups.append((m, pr, batman.TransitModel(pr, t[m], supersample_factor=ss, exp_time=ex)))
        self.n = len(t)

    def __call__(self, t0, per, rp, ars, b, u1, u2):
        inc = np.degrees(np.arccos(np.clip(b / ars, 0, 1)))
        out = np.empty(self.n)
        for m, pr, tm in self.groups:
            pr.t0, pr.per, pr.rp, pr.a, pr.inc, pr.u = t0, per, rp, ars, inc, [u1, u2]
            out[m] = tm.light_curve(pr)
        return out

def derived(theta, st_rad, st_mass, st_teff):
    t0, per, rp, ars, b, u1, u2 = theta
    rho = 3 * np.pi * ars ** 3 / (G_SI * (per * DAY) ** 2) / RHO_SUN
    rho_cat = (st_mass * MSUN / (4 / 3 * np.pi * (st_rad * RSUN) ** 3)) / RHO_SUN if st_mass > 0 and st_rad > 0 else np.nan
    rp_re = rp * st_rad * RSUN / REARTH
    dur = per / np.pi * np.arcsin(np.clip(np.sqrt(max((1 + rp) ** 2 - b ** 2, 0)) / ars, 0, 1)) * 24
    return dict(rp_re=rp_re, rp_rjup=rp * st_rad * RSUN / RJUP, depth_ppm=rp ** 2 * 1e6, dur_h=dur,
                a_au=ars * st_rad * RSUN / AU, teq=st_teff * np.sqrt(1 / (2 * ars)), rho_star=rho, rho_cat=rho_cat, inc=np.degrees(np.arccos(np.clip(b / ars, 0, 1))))
