"""Transit shape and stellar density: posterior T12/T23 from the existing MCMC chain, the eccentricity or
radius error that reconciles the transit density with the TIC, a grazing (V-shaped) alternative fit, and a
re-fit with the TIC density as a prior. Uses cached cubes and the undiluted sectors.

    python -W ignore scripts/deep_vet3.py --tic 259719966 --period 11.077207 --t0 3877.823115 --dur-h 2.372 --sectors 1 27 67 94 95 103
"""
import argparse, os, sys, time, json
import numpy as np, pandas as pd, emcee
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import least_squares
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tesshunt import pipeline as P
from tesshunt.cubes import load_tpf, ROOT
from tesshunt.transitfit import Model, derived, G_SI, DAY, RHO_SUN
from tesshunt.sweep_core import _central_aperture
ap = argparse.ArgumentParser(); ap.add_argument("--tic", type=int, required=True); ap.add_argument("--sector", type=int, default=1)
ap.add_argument("--period", type=float, required=True); ap.add_argument("--t0", type=float, required=True); ap.add_argument("--dur-h", type=float, required=True)
ap.add_argument("--sectors", type=int, nargs="*", required=True); ap.add_argument("--steps", type=int, default=2500); ap.add_argument("--burn", type=int, default=800)
a = ap.parse_args(); tic = a.tic; PER, T00, DUR = a.period, a.t0, a.dur_h / 24
OUT = os.path.join(ROOT, "plots", "candidates", "submission", f"s{a.sector}_tic{tic}"); DV = os.path.join(ROOT, "data", f"deepvet_{tic}"); FU = os.path.join(ROOT, "plots", "candidates", "followup")
def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)
tg = pd.read_csv(os.path.join(ROOT, "data", f"targets_s{a.sector}.csv")).drop_duplicates("ID").set_index("ID").loc[tic]
R_, M_, TEFF = float(tg.rad), float(tg.mass), float(tg.Teff); rho_cat = M_ / R_ ** 3
res = {}
# ---------- 1. posterior shape from the existing chain ----------
ch = np.load(os.path.join(FU, f"s{a.sector}_tic{tic}_mcmc.npz")); chain = ch["chain"]; names = list(ch["names"])
t0c, perc, rp, ars, b = [chain[:, names.index(n)] for n in ["t0", "period", "rp_rs", "a_rs", "b"]]
T14 = perc / np.pi * np.arcsin(np.sqrt(np.clip((1 + rp) ** 2 - b ** 2, 0, None)) / ars) * 24
T23 = perc / np.pi * np.arcsin(np.sqrt(np.clip((1 - rp) ** 2 - b ** 2, 0, None)) / ars) * 24
rho = 3 * np.pi * ars ** 3 / (G_SI * (perc * DAY) ** 2) / RHO_SUN
def q(x): lo, md, hi = np.percentile(x, [16, 50, 84]); return [round(float(md), 4), round(float(md - lo), 4), round(float(hi - md), 4)]
g = (rho / rho_cat) ** (1 / 3)                       # photo-eccentric factor (1 + e sin w)/sqrt(1-e^2)
e_min = (g ** 2 - 1) / (g ** 2 + 1)                   # minimum e (w = 90 deg) that reconciles the densities
R_needed = R_ * (rho_cat / rho) ** (1 / 3)            # or: the radius (at TIC mass) that would make the star this dense
res["shape"] = dict(T14_h=q(T14), T23_h=q(T23), T12_h=q((T14 - T23) / 2), T23_over_T14=q(T23 / T14), b=q(b), rp_rs=q(rp), a_rs=q(ars), rho_transit=q(rho), rho_tic=round(rho_cat, 3),
                    rho_ratio=q(rho / rho_cat), e_min_if_density_is_tic=q(e_min), R_star_if_circular_and_TIC_mass=q(R_needed), R_tic=R_,
                    frac_posterior_grazing=float(np.mean(b > 1 - rp)), frac_posterior_b_gt_0p8=float(np.mean(b > 0.8)), rho_sigma_tension=round(float((np.median(rho) - rho_cat) / np.std(rho)), 2))
log("shape: " + json.dumps(res["shape"]))
# ---------- 2. gather near-transit data from the cached cubes (same recipe as transitfit.gather) ----------
T, F, E, X, S = [], [], [], [], []
for s in a.sectors:
    tpf = load_tpf(tic, s, float(tg.ra), float(tg.dec))
    try:
        lc, raw, aper = P.lightcurve_from_cutout(tpf, threshold=5)
        if aper.sum() == 0 or aper.sum() > 0.5 * aper.size: raise ValueError
    except Exception: lc, raw, aper = _central_aperture(tpf)
    flat = P.mask_jitter(P.flatten(P.clean(lc), 33), 6, 3); t, f = flat.time.value, flat.flux.value; exp = float(np.nanmedian(np.diff(t)))
    ph = ((t - T00 + PER / 2) % PER) - PER / 2; keep = np.abs(ph) < 3 * DUR; oot = (np.abs(ph) > DUR) & keep; sig = np.nanstd(f[oot])
    T.append(t[keep]); F.append(f[keep]); E.append(np.full(keep.sum(), sig)); X.append(np.full(keep.sum(), exp)); S.append(np.full(keep.sum(), s))
t, f, e, exp, sec = map(np.concatenate, (T, F, E, X, S)); model = Model(t, exp); log(f"{len(t)} points from sectors {a.sectors}")
u1, u2 = 0.355, 0.332
def chi2(th): return float(np.sum(((f - model(*th)) / e) ** 2))
# ---------- 3. free fit vs grazing-forced (V-shaped) fit ----------
x0 = [T00, PER, 0.196, 37.8, 0.55]
free = least_squares(lambda th: (f - model(th[0], th[1], th[2], th[3], th[4], u1, u2)) / e, x0, bounds=([T00 - .05, PER * .999, .01, 2, 0], [T00 + .05, PER * 1.001, 0.9, 300, 1.5]))
graz = []
for rp_g in [0.25, 0.35, 0.5, 0.7]:
    for b_g in [1 - rp_g + 0.02, 1.0]:
        r = least_squares(lambda th: (f - model(th[0], th[1], rp_g, th[2], b_g, u1, u2)) / e, [T00, PER, 20.0], bounds=([T00 - .05, PER * .999, 2], [T00 + .05, PER * 1.001, 300]))
        graz.append(dict(rp_rs=rp_g, b=b_g, a_rs=round(r.x[2], 2), chi2=round(float(np.sum(r.fun ** 2)), 1)))
graz = pd.DataFrame(graz); chi2_free = float(np.sum(free.fun ** 2))
res["grazing"] = dict(free_fit=dict(rp_rs=round(free.x[2], 4), a_rs=round(free.x[3], 2), b=round(free.x[4], 3), chi2=round(chi2_free, 1), n=int(len(f))),
                      grazing_fits=graz.to_dict("records"), best_grazing_delta_chi2=round(float(graz.chi2.min() - chi2_free), 1))
log("grazing: " + json.dumps(res["grazing"]))
# ---------- 4. MCMC with the TIC density as a prior (rho = 4.08 +/- 25%), to see what b and Rp/Rs the catalog star implies ----------
def lnprob(th, prior):
    t0, per, rp_, ars_, b_, u1_, u2_ = th
    if not (abs(t0 - T00) < .05 and PER * .999 < per < PER * 1.001 and .01 < rp_ < .6 and 2 < ars_ < 300 and 0 <= b_ < 1 + rp_ and 0 <= u1_ <= 1 and 0 <= u1_ + u2_ <= 1 and u2_ > -1): return -np.inf
    lp = -0.5 * (((u1_ - u1) / .1) ** 2 + ((u2_ - u2) / .1) ** 2)
    if prior:
        rho_ = 3 * np.pi * ars_ ** 3 / (G_SI * (per * DAY) ** 2) / RHO_SUN; lp += -0.5 * ((rho_ - rho_cat) / (0.25 * rho_cat)) ** 2
    return lp - 0.5 * np.sum(((f - model(*th)) / e) ** 2)
p0 = np.array([free.x[0], free.x[1], free.x[2], free.x[3], free.x[4], u1, u2]); nw = 28
walkers = p0 + np.array([1e-4, 1e-6, 1e-3, 0.3, 0.02, 0.02, 0.02]) * np.random.randn(nw, 7); walkers[:, 4] = np.clip(walkers[:, 4], 0, 0.9)
sam = emcee.EnsembleSampler(nw, 7, lnprob, args=(True,)); sam.run_mcmc(walkers, a.steps, progress=False); cp = sam.get_chain(discard=a.burn, flat=True)
rp_p, ars_p, b_p = cp[:, 2], cp[:, 3], cp[:, 4]; rho_p = 3 * np.pi * ars_p ** 3 / (G_SI * (cp[:, 1] * DAY) ** 2) / RHO_SUN
res["density_prior_fit"] = dict(prior="rho = TIC +/- 25%", rp_rs=q(rp_p), a_rs=q(ars_p), b=q(b_p), rho=q(rho_p), rp_re=q(rp_p * R_ * 109.1), depth_ppm=q(rp_p ** 2 * 1e6),
                                chi2_at_median=round(chi2([np.median(cp[:, i]) for i in range(7)]), 1), chi2_free=round(chi2_free, 1), acceptance=round(float(np.mean(sam.acceptance_fraction)), 2))
log("density prior: " + json.dumps(res["density_prior_fit"]))
json.dump(res, open(os.path.join(DV, "shape_density.json"), "w"), indent=2)
# ---------- figure ----------
fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
ax[0].hist(T23 / T14, 60, color="0.5"); ax[0].set_xlabel("T23 / T14 (flat-bottom fraction)"); ax[0].set_title(f"posterior: {np.median(T23/T14):.2f}; V-shape = 0; grazing fraction {res['shape']['frac_posterior_grazing']:.3f}", fontsize=9)
ax[1].hist(rho, 60, color="0.5", label="transit (free fit)"); ax[1].hist(rho_p, 60, color="tab:blue", alpha=0.6, label="with TIC prior"); ax[1].axvline(rho_cat, color="r", label=f"TIC {rho_cat:.2f}"); ax[1].set_xlabel("stellar density (solar)"); ax[1].legend(fontsize=8)
ax[1].set_title(f"e >= {np.median(e_min):.2f} or R* = {np.median(R_needed):.3f} Rsun reconciles", fontsize=9)
ph = ((t - free.x[0] + PER / 2) % PER) - PER / 2; o = np.argsort(ph); tt_ = np.linspace(-2 * DUR, 2 * DUR, 1000); mm = Model(tt_ + free.x[0], np.full(len(tt_), 200 / 86400))
bins = np.linspace(-2 * DUR, 2 * DUR, 70); idx = np.digitize(ph, bins); bm = [np.nanmean(f[idx == i]) if (idx == i).sum() > 2 else np.nan for i in range(1, len(bins))]
ax[2].plot((bins[:-1] + bins[1:]) / 2 * 24, bm, "ko", ms=3); ax[2].plot(tt_ * 24, mm(free.x[0], free.x[1], free.x[2], free.x[3], free.x[4], u1, u2), "r-", label=f"free: b={free.x[4]:.2f} chi2={chi2_free:.0f}")
gb = graz.sort_values("chi2").iloc[0]; r = least_squares(lambda th: (f - model(th[0], th[1], gb.rp_rs, th[2], gb.b, u1, u2)) / e, [T00, PER, 20.0], bounds=([T00 - .05, PER * .999, 2], [T00 + .05, PER * 1.001, 300]))
ax[2].plot(tt_ * 24, mm(r.x[0], r.x[1], gb.rp_rs, r.x[2], gb.b, u1, u2), "b--", label=f"best grazing: rp={gb.rp_rs} b={gb.b} dchi2=+{gb.chi2-chi2_free:.0f}")
ax[2].legend(fontsize=8); ax[2].set_xlabel("hours from mid-transit"); ax[2].set_title("binned fold (undiluted sectors) vs flat-bottom and V-shaped models", fontsize=9)
plt.tight_layout(); fig.savefig(os.path.join(OUT, "addendum_shape_density.png"), dpi=110); log("done")
