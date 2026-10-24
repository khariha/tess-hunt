"""Fit a physical transit model to all sectors of one candidate.

Free parameters: t0, period, Rp/Rs, a/Rs, impact parameter b.
Limb darkening fixed (quadratic, TESS band, cool dwarf).
Each sector is super-sampled to its own exposure time so 30-min and 200-s
data are modeled consistently. Reports planet radius, orbit geometry, and the
stellar density implied by the transit vs the catalog density (a planet on a
circular orbit must agree; a blended binary usually does not).

    python -W ignore scripts/fit_transit.py --sector 1 --tic 389723320
"""
import argparse, os, sys, time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import batman
from scipy.optimize import least_squares
from tess_stars2px import tess_stars2px_function_entry

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tesshunt import pipeline as P
from tesshunt.sweep_core import get_tpf, _central_aperture

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "plots", "candidates", "followup")
os.makedirs(OUT, exist_ok=True)
G_SI, RSUN, MSUN, REARTH, DAY = 6.674e-11, 6.957e8, 1.989e30, 6.371e6, 86400.0

ap = argparse.ArgumentParser()
ap.add_argument("--sector", type=int, default=1); ap.add_argument("--tic", type=int, required=True)
ap.add_argument("--period", type=float, default=None); ap.add_argument("--t0", type=float, default=None)
ap.add_argument("--u1", type=float, default=0.35); ap.add_argument("--u2", type=float, default=0.25)
args = ap.parse_args()
S0, tic = args.sector, args.tic
def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)

sw = pd.read_csv(os.path.join(DATA, f"sweep_s{S0}.csv")).set_index("ID").loc[tic]
tg = pd.read_csv(os.path.join(DATA, f"targets_s{S0}.csv")).drop_duplicates("ID").set_index("ID").loc[tic]
P0 = args.period or float(sw.period_d); t00 = args.t0 or float(sw.t0); dur0 = float(sw.duration_h) / 24
sectors = sorted(set(int(s) for s in tess_stars2px_function_entry(tic, tg.ra, tg.dec)[3]))

# ---- gather data: keep only points within 3 durations of a predicted transit ----
T, F, E, EXP = [], [], [], []
for s in sectors:
    try:
        tpf = get_tpf(tg.ra, tg.dec, s, 15, "s3")
        if tpf is None: continue
        try:
            lc, raw, aper = P.lightcurve_from_cutout(tpf, threshold=5)
            if aper.sum() == 0 or aper.sum() > 0.5 * aper.size: raise ValueError
        except Exception:
            lc, raw, aper = _central_aperture(tpf)
        flat = P.mask_jitter(P.flatten(P.clean(lc), 33), 6, 3)
        t, f = flat.time.value, flat.flux.value
        exp = float(np.nanmedian(np.diff(t)))
        ph = ((t - t00 + P0 / 2) % P0) - P0 / 2
        keep = np.abs(ph) < 3 * dur0
        oot = (np.abs(ph) > dur0) & keep
        sig = np.nanstd(f[oot]) if oot.sum() > 10 else np.nanstd(f)
        T.append(t[keep]); F.append(f[keep]); E.append(np.full(keep.sum(), sig)); EXP.append(np.full(keep.sum(), exp))
        log(f"sector {s}: {keep.sum()} points near transit, {exp*1440:.0f}-min cadence, scatter {sig*1e6:.0f} ppm")
    except Exception as e:
        log(f"sector {s}: {type(e).__name__}: {e}")
t = np.concatenate(T); f = np.concatenate(F); e = np.concatenate(E); exp = np.concatenate(EXP)
base = t.max() - t.min()

# ---- model ----
def model(theta, t, exp):
    t0, per, rp, ars, b = theta
    inc = np.degrees(np.arccos(np.clip(b / ars, 0, 1)))
    out = np.empty_like(t)
    for ex in np.unique(exp):
        m = exp == ex
        pr = batman.TransitParams()
        pr.t0, pr.per, pr.rp, pr.a, pr.inc, pr.ecc, pr.w = t0, per, rp, ars, inc, 0.0, 90.0
        pr.u, pr.limb_dark = [args.u1, args.u2], "quadratic"
        ss = max(int(round(ex * 1440 / 2)), 1)        # ~2-min sub-steps
        out[m] = batman.TransitModel(pr, t[m], supersample_factor=ss, exp_time=ex).light_curve(pr)
    return out

def resid(theta): return (f - model(theta, t, exp)) / e
# initial guess from the sweep numbers
rp0 = np.sqrt(float(sw.depth_ppm) / 1e6)
# refine t0 to the middle of the baseline so period and epoch are not correlated
n_mid = np.round((np.median(t) - t00) / P0); t0_mid = t00 + n_mid * P0
ars0 = (P0 * DAY) / (np.pi * dur0 * DAY) if dur0 > 0 else 15.0
x0 = [t0_mid, P0, rp0, max(ars0, 3.0), 0.3]
lb = [t0_mid - 0.1, P0 * 0.995, 0.005, 1.5, 0.0]; ub = [t0_mid + 0.1, P0 * 1.005, 0.6, 200.0, 1.2]
fit = least_squares(resid, x0, bounds=(lb, ub), x_scale=[0.001, 1e-5, 0.01, 1.0, 0.1])
theta = fit.x
chi2 = float(np.sum(fit.fun ** 2)); dof = len(f) - len(theta)
# rough uncertainties from the Jacobian, scaled by reduced chi2
try:
    cov = np.linalg.inv(fit.jac.T @ fit.jac) * max(chi2 / dof, 1.0)
    err = np.sqrt(np.diag(cov))
except Exception:
    err = np.full(len(theta), np.nan)
t0f, Pf, rpf, arsf, bf = theta
rho_transit = 3 * np.pi * arsf ** 3 / (G_SI * (Pf * DAY) ** 2)                 # kg/m^3
rho_sun = 1408.0
rho_cat = float(tg.mass) * MSUN / (4 / 3 * np.pi * (float(tg.rad) * RSUN) ** 3) if tg.mass > 0 and tg.rad > 0 else np.nan
rp_re = rpf * float(tg.rad) * RSUN / REARTH
dur_full = Pf / np.pi * np.arcsin(np.sqrt((1 + rpf) ** 2 - bf ** 2) / arsf) * 24 if arsf > 1 + rpf else np.nan
a_au = arsf * float(tg.rad) * RSUN / 1.496e11
teq = float(tg.Teff) * np.sqrt(1 / (2 * arsf)) if arsf > 0 else np.nan

res = {
    "t0_btjd": (t0f, err[0]), "period_d": (Pf, err[1]), "rp_rs": (rpf, err[2]), "a_rs": (arsf, err[3]), "impact_b": (bf, err[4]),
    "planet_radius_earth": rp_re, "planet_radius_jup": rp_re / 11.2, "depth_ppm": rpf ** 2 * 1e6,
    "duration_h": dur_full, "a_au": a_au, "T_eq_K_(albedo0)": teq,
    "stellar_density_from_transit_(solar)": rho_transit / rho_sun, "stellar_density_catalog_(solar)": rho_cat / rho_sun if np.isfinite(rho_cat) else np.nan,
    "chi2_per_dof": chi2 / dof, "n_points": len(f), "baseline_d": base, "sectors": sectors,
}
log("fit results:")
for k, v in res.items():
    print(f"  {k:40s} {v}")

# ---- plot ----
ph = ((t - t0f + Pf / 2) % Pf) - Pf / 2
order = np.argsort(ph)
tt = np.linspace(-3 * dur0, 3 * dur0, 2000)
fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
for ex in np.unique(exp):
    m = exp == ex
    ax[0].scatter(ph[m] * 24, f[m], s=3, alpha=0.4, label=f"{ex*1440:.0f}-min data")
ax[0].plot(tt * 24, model([0, Pf, rpf, arsf, bf], tt, np.full(len(tt), 120 / 86400)), "r-", lw=2, label="fit (200-s sampling)")
lo, hi = np.nanpercentile(f, [0.5, 99.5]); ax[0].set_ylim(lo - 0.3 * (hi - lo), hi + 0.3 * (hi - lo))
ax[0].legend(fontsize=8); ax[0].set_ylabel("relative flux")
ax[0].set_title(f"TIC {tic}: P={Pf:.6f} d  Rp/Rs={rpf:.4f}  b={bf:.2f}  a/Rs={arsf:.1f}  ->  Rp={rp_re:.2f} Re  rho*={rho_transit/rho_sun:.1f} (catalog {rho_cat/rho_sun:.1f}) solar")
bins = np.linspace(-3 * dur0, 3 * dur0, 60)
idx = np.digitize(ph, bins); bm = [np.nanmean(f[idx == i]) if (idx == i).sum() > 2 else np.nan for i in range(1, len(bins))]
ax[1].plot((bins[:-1] + bins[1:]) / 2 * 24, bm, "ko", ms=4, label="binned data")
ax[1].plot(tt * 24, model([0, Pf, rpf, arsf, bf], tt, np.full(len(tt), 120 / 86400)), "r-", lw=2)
ax[1].set_xlabel("hours from mid-transit"); ax[1].set_ylabel("relative flux"); ax[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig(os.path.join(OUT, f"s{S0}_tic{tic}_fit.png"), dpi=110)
pd.Series({k: (v[0] if isinstance(v, tuple) else v) for k, v in res.items()}).to_csv(os.path.join(OUT, f"s{S0}_tic{tic}_fit.csv"))
log(f"wrote {OUT}/s{S0}_tic{tic}_fit.png")
