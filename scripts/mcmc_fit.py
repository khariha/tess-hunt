"""MCMC transit fit with free quadratic limb darkening, for the submission.

    nohup python -W ignore scripts/mcmc_fit.py --sector 1 --tic 389723320 --period 3.927167 > data/mcmc_389723320.log 2>&1 &

Writes plots/candidates/followup/s<N>_tic<ID>_mcmc.{csv,png,npz}
"""
import argparse, os, sys, time
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import emcee
from scipy.optimize import least_squares
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tesshunt.transitfit import gather, Model, derived

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data"); OUT = os.path.join(ROOT, "plots", "candidates", "followup"); os.makedirs(OUT, exist_ok=True)
ap = argparse.ArgumentParser()
ap.add_argument("--sector", type=int, default=1); ap.add_argument("--tic", type=int, required=True)
ap.add_argument("--period", type=float, default=None); ap.add_argument("--walkers", type=int, default=32)
ap.add_argument("--steps", type=int, default=3000); ap.add_argument("--burn", type=int, default=1000)
ap.add_argument("--sectors", type=int, nargs="*", default=None, help="restrict to these sectors (e.g. undiluted apertures)")
ap.add_argument("--max-b", type=float, default=None, help="cap impact parameter (e.g. 1.0 to forbid grazing-only solutions)")
ap.add_argument("--u1", type=float, default=0.30); ap.add_argument("--u2", type=float, default=0.30); ap.add_argument("--u-sigma", type=float, default=0.10)
args = ap.parse_args(); S0, tic = args.sector, args.tic
def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)

sw = pd.read_csv(os.path.join(DATA, f"sweep_s{S0}.csv")).set_index("ID").loc[tic]
tg = pd.read_csv(os.path.join(DATA, f"targets_s{S0}.csv")).drop_duplicates("ID").set_index("ID").loc[tic]
P0 = args.period or float(sw.period_d); t00 = float(sw.t0); dur0 = float(sw.duration_h) / 24
d = gather(tic, tg.ra, tg.dec, P0, t00, dur0, log=log, sectors=args.sectors)
t, f, e, exp = d["t"], d["f"], d["e"], d["exp"]
model = Model(t, exp)
n_mid = np.round((np.median(t) - t00) / P0); t0_mid = t00 + n_mid * P0
rp0 = np.sqrt(float(sw.depth_ppm) / 1e6)

def lnprior(th):
    t0, per, rp, ars, b, u1, u2 = th
    bmax = args.max_b if args.max_b is not None else 1 + rp
    if not (abs(t0 - t0_mid) < 0.1 and P0 * 0.995 < per < P0 * 1.005 and 0.005 < rp < 0.6 and 1.5 < ars < 300 and 0 <= b < bmax): return -np.inf
    if not (0 <= u1 <= 1 and 0 <= u1 + u2 <= 1 and u2 > -1): return -np.inf
    return -0.5 * (((u1 - args.u1) / args.u_sigma) ** 2 + ((u2 - args.u2) / args.u_sigma) ** 2)
def lnprob(th):
    lp = lnprior(th)
    if not np.isfinite(lp): return -np.inf
    return lp - 0.5 * np.sum(((f - model(*th)) / e) ** 2)

x0 = np.array([t0_mid, P0, rp0, 15.0, 0.5, args.u1, args.u2])
ls = least_squares(lambda th: np.concatenate([(f - model(*th)) / e, [(th[5] - args.u1) / args.u_sigma, (th[6] - args.u2) / args.u_sigma]]), x0,
                   bounds=([t0_mid - 0.1, P0 * 0.995, 0.005, 1.5, 0, 0, -0.5], [t0_mid + 0.1, P0 * 1.005, 0.6, 300, 1.2, 1, 1]))
x0 = ls.x; log(f"least squares start: {x0}")
ndim = len(x0); p0 = x0 + np.array([1e-4, 1e-6, 1e-3, 0.3, 0.02, 0.02, 0.02]) * np.random.randn(args.walkers, ndim)
p0[:, 4] = np.clip(p0[:, 4], 0, 0.95)
sampler = emcee.EnsembleSampler(args.walkers, ndim, lnprob)
t_start = time.time()
sampler.run_mcmc(p0, args.steps, progress=False)
log(f"mcmc done in {(time.time()-t_start)/60:.1f} min; mean acceptance {np.mean(sampler.acceptance_fraction):.2f}")
chain = sampler.get_chain(discard=args.burn, flat=True)
names = ["t0", "period", "rp_rs", "a_rs", "b", "u1", "u2"]
rows = []
for i, n in enumerate(names):
    lo, md, hi = np.percentile(chain[:, i], [16, 50, 84]); rows.append(dict(param=n, median=md, minus=md - lo, plus=hi - md))
der = np.array([[v for v in derived(th, float(tg.rad), float(tg.mass), float(tg.Teff)).values()] for th in chain[::20]])
dnames = list(derived(chain[0], float(tg.rad), float(tg.mass), float(tg.Teff)).keys())
for i, n in enumerate(dnames):
    lo, md, hi = np.percentile(der[:, i], [16, 50, 84]); rows.append(dict(param=n, median=md, minus=md - lo, plus=hi - md))
res = pd.DataFrame(rows); res.to_csv(os.path.join(OUT, f"s{S0}_tic{tic}_mcmc.csv"), index=False)
np.savez(os.path.join(OUT, f"s{S0}_tic{tic}_mcmc.npz"), chain=chain, names=names)
print(res.to_string(index=False))
med = res.set_index("param").loc[names, "median"].values
chi2 = np.sum(((f - model(*med)) / e) ** 2); log(f"chi2/dof at median = {chi2/(len(f)-ndim):.3f}, n={len(f)}")
# plot
ph = ((t - med[0] + med[1] / 2) % med[1]) - med[1] / 2
tt = np.linspace(-3 * dur0, 3 * dur0, 1500); mm = Model(tt + med[0], np.full(len(tt), 120 / 86400))
fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
ax[0].scatter(ph * 24, f, s=3, alpha=0.35, c=exp * 1440, cmap="viridis"); ax[0].plot(tt * 24, mm(*med), "r-", lw=2)
lo_, hi_ = np.nanpercentile(f, [0.5, 99.5]); ax[0].set_ylim(lo_ - 0.3 * (hi_ - lo_), hi_ + 0.3 * (hi_ - lo_)); ax[0].set_ylabel("relative flux")
r = res.set_index("param")
ax[0].set_title(f"TIC {tic} MCMC: P={r.loc['period','median']:.6f} d  Rp={r.loc['rp_re','median']:.2f}+{r.loc['rp_re','plus']:.2f}-{r.loc['rp_re','minus']:.2f} Re  b={r.loc['b','median']:.2f}  rho*={r.loc['rho_star','median']:.1f} (cat {r.loc['rho_cat','median']:.1f})")
bins = np.linspace(-3 * dur0, 3 * dur0, 60); idx = np.digitize(ph, bins)
bm = [np.nanmean(f[idx == i]) if (idx == i).sum() > 2 else np.nan for i in range(1, len(bins))]
ax[1].plot((bins[:-1] + bins[1:]) / 2 * 24, bm, "ko", ms=4); ax[1].plot(tt * 24, mm(*med), "r-", lw=2); ax[1].set_xlabel("hours from mid-transit"); ax[1].set_ylabel("relative flux (binned)")
plt.tight_layout(); plt.savefig(os.path.join(OUT, f"s{S0}_tic{tic}_mcmc.png"), dpi=110)
try:
    import corner
    fig = corner.corner(chain, labels=names, quantiles=[0.16, 0.5, 0.84], show_titles=True, title_fmt=".4f")
    fig.savefig(os.path.join(OUT, f"s{S0}_tic{tic}_corner.png"), dpi=80)
except Exception as ex:
    log(f"corner plot skipped: {ex}")
log("done")
