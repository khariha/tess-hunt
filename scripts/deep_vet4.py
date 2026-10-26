"""Joint primary + secondary eclipse fit (batman, eccentric orbit, emcee) from cached cubes.
Both events are masked before flattening so the smoother does not put shoulders on them.

    python -W ignore scripts/deep_vet4.py --tic 259719966 --period 11.077207 --t0 3877.823115 --dur-h 2.372 --sec-phase 0.54 --sectors 1 27 67 94 95 102 103 104
"""
import argparse, os, sys, time, json
import numpy as np, pandas as pd, emcee, batman
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import least_squares
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tesshunt import pipeline as P
from tesshunt.cubes import load_tpf, ROOT
from tesshunt.transitfit import G_SI, DAY, RHO_SUN
import lightkurve as lk
ap = argparse.ArgumentParser(); ap.add_argument("--tic", type=int, required=True); ap.add_argument("--sector", type=int, default=1)
ap.add_argument("--period", type=float, required=True); ap.add_argument("--t0", type=float, required=True); ap.add_argument("--dur-h", type=float, required=True)
ap.add_argument("--sec-phase", type=float, required=True); ap.add_argument("--sectors", type=int, nargs="*", required=True)
ap.add_argument("--steps", type=int, default=2500); ap.add_argument("--burn", type=int, default=800)
a = ap.parse_args(); tic = a.tic; PER, T00, DUR, SPH = a.period, a.t0, a.dur_h / 24, a.sec_phase
OUT = os.path.join(ROOT, "plots", "candidates", "submission", f"s{a.sector}_tic{tic}"); DV = os.path.join(ROOT, "data", f"deepvet_{tic}")
def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)
tg = pd.read_csv(os.path.join(ROOT, "data", f"targets_s{a.sector}.csv")).drop_duplicates("ID").set_index("ID").loc[tic]
gaia = pd.read_csv(os.path.join(ROOT, "data", f"gaia_neighbors_{tic}.csv")); gaia["light_frac"] = 10 ** (-0.4 * (gaia.phot_g_mean_mag - gaia.phot_g_mean_mag.iloc[0]))
R_, M_, TEFF = float(tg.rad), float(tg.mass), float(tg.Teff); rho_cat = M_ / R_ ** 3
def phase(t): return ((t - T00) / PER) % 1.0
# ---------- gather both windows, target-only aperture, both events masked before flattening ----------
T, F, E, X, S = [], [], [], [], []
for s in a.sectors:
    tpf = load_tpf(tic, s, float(tg.ra), float(tg.dec))
    tx, ty = [float(v[0]) for v in tpf.wcs.all_world2pix([tg.ra], [tg.dec], 0)]; gx, gy = tpf.wcs.all_world2pix(gaia.ra.values, gaia.dec.values, 0)
    nb = gaia.assign(px=gx, py=gy).iloc[1:].sort_values("light_frac", ascending=False).iloc[0]
    aper = tpf.create_threshold_mask(threshold=5, reference_pixel="center"); ny, nx = aper.shape; yy, xx = np.mgrid[0:ny, 0:nx]
    aper = aper & ~(np.hypot(xx - nb.px, yy - nb.py) <= 1.0)
    raw = tpf.to_lightcurve(aperture_mask=aper); reg = tpf.flux[:, ~aper].value; good = np.isfinite(raw.flux.value) & np.isfinite(reg).all(axis=1)
    lc = P.clean(lk.RegressionCorrector(raw[good]).correct(lk.DesignMatrix(reg[good], name="bg").pca(5).append_constant()))
    t = lc.time.value; ph = phase(t); dpri = np.abs(((ph + 0.5) % 1) - 0.5) * PER; dsec = np.abs(((ph - SPH + 0.5) % 1) - 0.5) * PER
    mask = (dpri < 1.0 * DUR) | (dsec < 1.2 * DUR)
    # flatten with the events masked: fit the trend on unmasked points only (lightkurve flatten mask argument)
    flat = lc.flatten(window_length=P.flatten(lc, 33).meta.get("WINDOW", None) or max(int(round(33 / 24 / np.nanmedian(np.diff(t)))) | 1, 5), mask=mask)
    flat = P.mask_jitter(flat, 6, 3); t, f = flat.time.value, flat.flux.value; ph = phase(t)
    dpri = np.abs(((ph + 0.5) % 1) - 0.5) * PER; dsec = np.abs(((ph - SPH + 0.5) % 1) - 0.5) * PER
    keep = (dpri < 3 * DUR) | (dsec < 3 * DUR); oot = keep & (dpri > DUR) & (dsec > 1.2 * DUR); sig = np.nanstd(f[oot]); exp = float(np.nanmedian(np.diff(t)))
    T.append(t[keep]); F.append(f[keep]); E.append(np.full(keep.sum(), sig)); X.append(np.full(keep.sum(), exp)); S.append(np.full(keep.sum(), s))
    log(f"sector {s}: {keep.sum()} points in the two windows, scatter {sig*1e6:.0f} ppm")
t, f, e, exp, sec = map(np.concatenate, (T, F, E, X, S))
# ---------- eccentric batman model with primary transit and secondary eclipse ----------
class EccModel:
    def __init__(self, t, exp):
        self.groups = []
        for ex in np.unique(exp):
            m = exp == ex; pr = batman.TransitParams(); pr.t0, pr.per, pr.rp, pr.a, pr.inc, pr.ecc, pr.w, pr.fp, pr.t_secondary = T00, PER, 0.2, 38., 89., 0., 90., 0.01, T00 + PER / 2
            pr.u, pr.limb_dark = [0.35, 0.33], "quadratic"; ss = max(int(round(ex * 1440 / 2)), 1)
            self.groups.append((m, pr, batman.TransitModel(pr, t[m], supersample_factor=ss, exp_time=ex), batman.TransitModel(pr, t[m], transittype="secondary", supersample_factor=ss, exp_time=ex)))
        self.n = len(t)
    def __call__(self, t0, rp, ars, b, fp, h, k, u1, u2):
        ecc = np.hypot(h, k); w = np.degrees(np.arctan2(k, h)) % 360.0
        out = np.empty(self.n)
        for m, pr, tm, ts in self.groups:
            # b is the transit impact parameter: b = a cos i (1-e^2)/(1+e sin w)
            cosi = b / ars * (1 + ecc * np.sin(np.radians(w))) / (1 - ecc ** 2)
            if not (0 <= cosi < 1): return np.full(self.n, np.inf)
            pr.t0, pr.rp, pr.a, pr.inc, pr.fp, pr.ecc, pr.w, pr.u = t0, rp, ars, np.degrees(np.arccos(cosi)), fp, ecc, w, [u1, u2]
            pr.t_secondary = tm.get_t_secondary(pr)
            out[m] = (tm.light_curve(pr) + ts.light_curve(pr) - 1.0) / (1.0 + fp)
        return out
model = EccModel(t, exp)
names = ["t0", "rp_rs", "a_rs", "b", "fp", "ecosw", "esinw", "u1", "u2"]
def lnprior(th):
    t0, rp, ars, b, fp, h, k, u1, u2 = th
    if not (abs(t0 - T00) < .05 and .01 < rp < .6 and 2 < ars < 300 and 0 <= b < 1 and 0 <= fp < .2 and np.hypot(h, k) < .9 and 0 <= u1 <= 1 and 0 <= u1 + u2 <= 1): return -np.inf
    return -0.5 * (((u1 - .355) / .1) ** 2 + ((u2 - .332) / .1) ** 2)
def lnprob(th):
    lp = lnprior(th)
    if not np.isfinite(lp): return -np.inf
    mm = model(*th)
    if not np.all(np.isfinite(mm)): return -np.inf
    return lp - 0.5 * np.sum(((f - mm) / e) ** 2)
# starting point: e cos w from the secondary phase, e sin w ~ 0
ecosw0 = (SPH - 0.5) * np.pi / 2; w0 = 0.0; ecc0 = abs(ecosw0)
best = None
for h0 in [0.04, 0.063, 0.09]:
    for k0 in [-0.15, 0.0, 0.15]:
        x0 = [T00, 0.196, 37.8, 0.55, 0.009, h0, k0, .355, .332]
        r = least_squares(lambda th: np.concatenate([(f - model(*th)) / e, [(th[7] - .355) / .1, (th[8] - .332) / .1]]), x0, bounds=([T00 - .05, .01, 2, 0, 0, -.5, -.5, 0, -.5], [T00 + .05, .6, 300, .99, .2, .5, .5, 1, 1]))
        c2 = float(np.sum(r.fun ** 2))
        if best is None or c2 < best[0]: best = (c2, r.x)
x0 = best[1].copy(); x0[7], x0[8] = .355, .332; log(f"least squares: chi2={best[0]:.1f} " + json.dumps(dict(zip(names, [round(float(v), 5) for v in x0]))))
nw = 36; p0 = x0 + np.array([1e-4, 1e-3, .3, .02, 5e-4, .005, .01, .02, .02]) * np.random.randn(nw, 9); p0[:, 3] = np.clip(p0[:, 3], 0, .95); p0[:, 4] = np.abs(p0[:, 4]); p0[:, 7] = np.clip(p0[:, 7], 0.02, 0.6); p0[:, 8] = np.clip(p0[:, 8], 0.0, 0.5)
assert all(np.isfinite(lnprob(w_)) for w_ in p0), 'bad start'
sam = emcee.EnsembleSampler(nw, 9, lnprob); t1 = time.time(); sam.run_mcmc(p0, a.steps, progress=False); log(f"mcmc {(time.time()-t1)/60:.1f} min, acceptance {np.mean(sam.acceptance_fraction):.2f}")
ch = sam.get_chain(discard=a.burn, flat=True)
def q(x): lo, md, hi = np.percentile(x, [16, 50, 84]); return [round(float(md), 5), round(float(md - lo), 5), round(float(hi - md), 5)]
t0, rp, ars, b, fp, ecosw, esinw = [ch[:, i] for i in range(7)]
ecc = np.hypot(ecosw, esinw); w = np.degrees(np.arctan2(esinw, ecosw)) % 360
rho_circ = 3 * np.pi * ars ** 3 / (G_SI * (PER * DAY) ** 2) / RHO_SUN     # a/R* -> density is e-independent (Kepler), the earlier 5.9 came from forcing e=0 on the duration
depth_sec = fp / (1 + fp)                                               # eclipse depth = companion flux / total
J = fp / rp ** 2                                                       # TESS-band surface-brightness ratio companion/host
# secondary phase & durations
pr = batman.TransitParams(); pr.t0, pr.per, pr.rp, pr.a, pr.ecc, pr.w, pr.u, pr.limb_dark, pr.fp = T00, PER, .2, 38, 0, 90, [.3, .3], "quadratic", .01
tsec = []
for th in ch[::50]:
    e_, w_ = np.hypot(th[5], th[6]), np.degrees(np.arctan2(th[6], th[5])) % 360
    cosi = th[3] / th[2] * (1 + e_ * np.sin(np.radians(w_))) / (1 - e_ ** 2); pr.t0, pr.rp, pr.a, pr.inc, pr.ecc, pr.w = th[0], th[1], th[2], np.degrees(np.arccos(cosi)), e_, w_
    tsec.append(batman.TransitModel(pr, np.array([th[0]])).get_t_secondary(pr))
tsec = np.array(tsec); phsec = ((tsec - ch[::50, 0]) / PER) % 1
# companion temperature from J: blackbody ratio integrated over the TESS band (600-1000 nm) at the host Teff
lam = np.linspace(600e-9, 1000e-9, 200); h, c, kB = 6.626e-34, 2.998e8, 1.381e-23
def bb(Tk): return np.trapezoid(1 / (np.exp(h * c / (lam * kB * Tk)) - 1) / lam ** 5, lam)
Ts = np.arange(1200, 3500, 5); ratios = np.array([bb(T_) / bb(TEFF) for T_ in Ts])
T2 = np.interp(np.percentile(J, [16, 50, 84]), ratios, Ts)
res = dict(sectors=a.sectors, n_points=int(len(t)), chi2_at_median=round(float(np.sum(((f - model(*np.median(ch, axis=0))) / e) ** 2)), 1),
           params={n: q(ch[:, i]) for i, n in enumerate(names)}, e_cos_w=q(ecosw), e_sin_w=q(esinw), ecc=q(ecc), omega_deg=q(w), secondary_phase=q(phsec),
           secondary_depth_ppm=q(depth_sec * 1e6), primary_depth_ppm=q(rp ** 2 * 1e6), surface_brightness_ratio_TESS=q(J),
           companion_Teff_blackbody_K=[round(float(T2[1])), round(float(T2[1] - T2[0])), round(float(T2[2] - T2[1]))], host_Teff_K=TEFF,
           rho_star_from_a_rs=q(rho_circ), rho_tic=round(rho_cat, 3), companion_radius_rsun=q(rp * R_), companion_radius_rjup=q(rp * R_ * 9.73),
           note="rho_star_from_a_rs uses Kepler's law on the fitted a/R*, which with e free is no longer the circular-duration density")
json.dump(res, open(os.path.join(DV, "joint_eccentric_fit.json"), "w"), indent=2); np.savez(os.path.join(DV, "joint_eccentric_chain.npz"), chain=ch, names=names)
log(json.dumps(res, indent=1))
# ---------- figure ----------
med = np.median(ch, axis=0); ph = phase(t); fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
for k, (center, title, span) in enumerate([(0.0, "primary transit", 3 * DUR), (float(np.median(phsec)), "secondary eclipse", 3 * DUR)]):
    d = (((ph - center + 0.5) % 1) - 0.5) * PER; m = np.abs(d) < span
    ax[k].plot(d[m] * 24, f[m], ".", ms=1.5, alpha=0.25, color="0.6")
    bins = np.linspace(-span, span, 50); idx = np.digitize(d[m], bins); bm = [np.nanmean(f[m][idx == i]) if (idx == i).sum() > 2 else np.nan for i in range(1, len(bins))]
    be = [np.nanstd(f[m][idx == i]) / np.sqrt((idx == i).sum()) if (idx == i).sum() > 2 else np.nan for i in range(1, len(bins))]
    ax[k].errorbar((bins[:-1] + bins[1:]) / 2 * 24, bm, be, fmt="ko", ms=3, capsize=1.5)
    tt_ = T00 + center * PER + np.linspace(-span, span, 800); mm = EccModel(tt_, np.full(len(tt_), 120 / 86400)); ax[k].plot((tt_ - T00 - center * PER) * 24, mm(*med), "r-", lw=1.5)
    ax[k].set_xlabel(f"hours from {title} centre"); ax[k].set_title(title, fontsize=10)
ax[0].set_ylim(0.945, 1.015); ax[1].set_ylim(0.985, 1.006)
ax[2].plot(ecosw[::20], esinw[::20], ".", ms=1, alpha=0.2); ax[2].set_xlabel("e cos w"); ax[2].set_ylabel("e sin w"); ax[2].axhline(0, color="0.6", lw=.5); ax[2].axvline(0, color="0.6", lw=.5)
ax[2].set_title(f"e = {np.median(ecc):.3f}, secondary at phase {np.median(phsec):.4f}, depth {np.median(depth_sec)*1e6:.0f} ppm, J = {np.median(J):.2f} -> T2 ~ {T2[1]:.0f} K", fontsize=9)
plt.tight_layout(); fig.savefig(os.path.join(OUT, "addendum_secondary_eclipse.png"), dpi=110); log("done")
