"""Second-stage deep vetting (after deep_vet.py): PSF-based dilution and the two blend hypotheses, predicted vs
observed centroid-shift vectors, per-transit difference-image centroids, sector-by-sector secondary scan,
short-period rotation search, flares in transit, jitter-mask / quality-flag audit near transits, sector-103 deep dive.

    python -W ignore scripts/deep_vet2.py --tic 259719966 --period 11.077207 --t0 3877.823115 --dur-h 2.372
"""
import argparse, os, sys, time, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import lightkurve as lk
from astropy.timeseries import LombScargle
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tesshunt import pipeline as P
from tesshunt.cubes import load_tpf, ROOT
from tesshunt.sweep_core import _central_aperture
from tess_stars2px import tess_stars2px_function_entry

ap = argparse.ArgumentParser()
ap.add_argument("--tic", type=int, required=True); ap.add_argument("--sector", type=int, default=1)
ap.add_argument("--period", type=float, required=True); ap.add_argument("--t0", type=float, required=True); ap.add_argument("--dur-h", type=float, required=True)
ap.add_argument("--psf-sigma", type=float, default=0.75, help="Gaussian PSF sigma in px (TESS ~ 0.7-0.8 px)")
a = ap.parse_args(); tic = a.tic; PER, T0, DUR = a.period, a.t0, a.dur_h / 24
OUT = os.path.join(ROOT, "plots", "candidates", "submission", f"s{a.sector}_tic{tic}"); DV = os.path.join(ROOT, "data", f"deepvet_{tic}")
def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)
tg = pd.read_csv(os.path.join(ROOT, "data", f"targets_s{a.sector}.csv")).drop_duplicates("ID").set_index("ID").loc[tic]
gaia = pd.read_csv(os.path.join(ROOT, "data", f"gaia_neighbors_{tic}.csv")); gaia["light_frac"] = 10 ** (-0.4 * (gaia.phot_g_mean_mag - gaia.phot_g_mean_mag.iloc[0]))
sectors = sorted(set(int(s) for s in tess_stars2px_function_entry(tic, tg.ra, tg.dec)[3]))
dep = pd.read_csv(os.path.join(DV, "depth_per_sector.csv"))
def phase_of(t): return ((t - T0 + PER / 2) % PER) - PER / 2
def build_lc(tpf, aper):
    raw = tpf.to_lightcurve(aperture_mask=aper); reg = tpf.flux[:, ~aper]; reg = reg.value if hasattr(reg, "value") else reg
    good = np.isfinite(raw.flux.value) & np.isfinite(reg).all(axis=1)
    dm = lk.DesignMatrix(reg[good], name="bg").pca(5).append_constant(); return lk.RegressionCorrector(raw[good]).correct(dm)
def psf_frac(px, py, aper, sig):
    """fraction of a Gaussian star's flux falling inside the aperture pixels (pixel-integrated via erf)"""
    from scipy.special import erf
    ys, xs = np.nonzero(aper); tot = 0.0
    for y, x in zip(ys, xs):
        fx = 0.5 * (erf((x + 0.5 - px) / (np.sqrt(2) * sig)) - erf((x - 0.5 - px) / (np.sqrt(2) * sig)))
        fy = 0.5 * (erf((y + 0.5 - py) / (np.sqrt(2) * sig)) - erf((y - 0.5 - py) / (np.sqrt(2) * sig)))
        tot += fx * fy
    return tot

rows_dil, rows_cen, rows_pt, rows_sec, rows_rot, rows_fl, rows_q = [], [], [], [], [], [], []
fig103 = None
for s in sectors:
    tpf = load_tpf(tic, s, float(tg.ra), float(tg.dec)); flux = tpf.flux.value; t = tpf.time.value; cad = float(np.nanmedian(np.diff(t)) * 1440)
    tx, ty = [float(v[0]) for v in tpf.wcs.all_world2pix([tg.ra], [tg.dec], 0)]
    gx, gy = tpf.wcs.all_world2pix(gaia.ra.values, gaia.dec.values, 0); G = gaia.assign(px=gx, py=gy)
    nb = G.iloc[1:].sort_values("light_frac", ascending=False).iloc[0]
    try:
        aper = tpf.create_threshold_mask(threshold=5, reference_pixel="center")
        if aper.sum() == 0 or aper.sum() > 0.5 * aper.size: raise ValueError
    except Exception: aper = _central_aperture(tpf)[2]
    ny, nx = aper.shape; yy, xx = np.mgrid[0:ny, 0:nx]
    aper_tgt = aper & ~(np.hypot(xx - nb.px, yy - nb.py) <= 1.0)
    if aper_tgt.sum() == 0: aper_tgt = np.hypot(xx - tx, yy - ty) <= 1.0
    aper_nb = np.hypot(xx - nb.px, yy - nb.py) <= 1.0
    # ---- PSF-based light budget in each aperture: L_i = f_i * A_i (A = aperture fraction of star i)
    for name, ap_ in [("pipeline", aper), ("target_only", aper_tgt), ("neighbour", aper_nb)]:
        A = np.array([psf_frac(r.px, r.py, ap_, a.psf_sigma) for _, r in G.iterrows()]); L = G.light_frac.values * A
        L_t, L_nb, L_other = L[0], L[nb.name], L.sum() - L[0] - L[nb.name]
        tot = L.sum()
        rows_dil.append(dict(sector=s, aperture=name, aper_px=int(ap_.sum()), frac_target_in_aper=round(A[0], 3), frac_nb_in_aper=round(A[nb.name], 3),
                             light_target=round(L_t, 4), light_nb=round(L_nb, 4), light_other=round(L_other, 4), dilution_psf=round((tot - L_t) / L_t, 4),
                             # H_target: dip on target with true depth D -> observed D * L_t/tot ; H_nb: dip on nb with eclipse E -> observed E * L_nb/tot
                             share_target=round(L_t / tot, 4), share_nb=round(L_nb / tot, 4)))
    # ---- predicted centroid shift vector (whole frame, after plane subtraction, as centroid_test measures it) for a dip on target vs on the neighbour
    good = np.isfinite(flux).all(axis=(1, 2)); fl, tt = flux[good], t[good]; ph = phase_of(tt)
    in_tr = np.abs(ph) < DUR / 2; oot = (np.abs(ph) > DUR) & (np.abs(ph) < 3 * DUR)
    fb = np.clip(P._subtract_plane_background(fl), 0, None)
    img = fb[oot].mean(axis=0); tot_img = img.sum(); cx0 = (img * xx).sum() / tot_img; cy0 = (img * yy).sum() / tot_img
    d_obs = float(dep[(dep.sector == s) & (dep.aperture == "pipeline")].depth_ppm.iloc[0]) / 1e6
    # flux lost in transit (in e-/s units of the frame): on target: the dip is d_obs of the aperture flux; model the lost light as a Gaussian at the dip star
    ap_flux = img[aper].sum(); lost = d_obs * ap_flux   # approximately the total flux removed from the frame during transit (PSF mostly inside aperture)
    def shift_if_dip_at(px, py):
        return np.array([(cx0 * tot_img - px * lost) / (tot_img - lost) - cx0, (cy0 * tot_img - py * lost) / (tot_img - lost) - cy0])
    # measured in-transit centroid of the same frames
    totf = fb.sum(axis=(1, 2)) + 1e-9; cx = (fb * xx).sum(axis=(1, 2)) / totf; cy = (fb * yy).sum(axis=(1, 2)) / totf
    sh = np.array([cx[in_tr].mean() - cx[oot].mean(), cy[in_tr].mean() - cy[oot].mean()])
    sh_err = np.array([cx[oot].std(), cy[oot].std()]) / np.sqrt(in_tr.sum())
    # direct-image star position measured the same way as the difference-image centroid (so WCS offsets cancel)
    def wcen(im, px, py, rad=1.5):
        w = np.clip(im, 0, None) * (np.hypot(xx - px, yy - py) <= rad); return (w * xx).sum() / w.sum(), (w * yy).sum() / w.sum()
    diff = fb[oot].mean(axis=0) - fb[in_tr].mean(axis=0); noise = np.sqrt(fb[oot].var(axis=0) / oot.sum() + fb[in_tr].var(axis=0) / in_tr.sum()) + 1e-9
    snr = diff / noise; py_, px_ = np.unravel_index(np.nanargmax(snr), snr.shape)
    dcx, dcy = wcen(diff, px_, py_); icx, icy = wcen(img, tx, ty)
    pred_t, pred_nb = shift_if_dip_at(icx, icy), shift_if_dip_at(nb.px, nb.py)
    chi_t = float(np.sum(((sh - pred_t) / sh_err) ** 2)); chi_nb = float(np.sum(((sh - pred_nb) / sh_err) ** 2))
    rows_cen.append(dict(sector=s, cadence_min=round(cad, 2), obs_dx=round(sh[0], 4), obs_dy=round(sh[1], 4), err=round(float(np.hypot(*sh_err)), 4),
                         pred_target_dx=round(pred_t[0], 4), pred_target_dy=round(pred_t[1], 4), pred_nb_dx=round(pred_nb[0], 4), pred_nb_dy=round(pred_nb[1], 4),
                         obs_mag_px=round(float(np.hypot(*sh)), 4), pred_target_mag_px=round(float(np.hypot(*pred_t)), 4), pred_nb_mag_px=round(float(np.hypot(*pred_nb)), 4),
                         chi2_target=round(chi_t, 1), chi2_nb=round(chi_nb, 1),
                         diffimg_minus_directimg_px=round(float(np.hypot(dcx - icx, dcy - icy)), 3), diffimg_minus_directimg_arcsec=round(float(np.hypot(dcx - icx, dcy - icy)) * 21, 1),
                         directimg_minus_wcs_px=round(float(np.hypot(icx - tx, icy - ty)), 3), peak_snr=round(float(snr[py_, px_]), 1)))
    log(f"S{s}: shift obs ({sh[0]:+.4f},{sh[1]:+.4f})+/-{np.hypot(*sh_err):.4f} | pred on-target ({pred_t[0]:+.4f},{pred_t[1]:+.4f}) chi2 {chi_t:.1f} | pred on-nb ({pred_nb[0]:+.4f},{pred_nb[1]:+.4f}) chi2 {chi_nb:.1f} | diff-direct {np.hypot(dcx-icx, dcy-icy)*21:.1f}\"")
    # ---- per-transit difference-image centroids
    epochs = np.unique(np.round((tt[in_tr] - T0) / PER).astype(int))
    for ep in epochs:
        tc = T0 + ep * PER; m_in = np.abs(tt - tc) < DUR / 2; m_out = (np.abs(tt - tc) > DUR) & (np.abs(tt - tc) < 3 * DUR)
        if m_in.sum() < 2 or m_out.sum() < 5: continue
        d_ = fb[m_out].mean(axis=0) - fb[m_in].mean(axis=0); n_ = np.sqrt(fb[m_out].var(axis=0) / m_out.sum() + fb[m_in].var(axis=0) / m_in.sum()) + 1e-9
        s_ = d_ / n_; py2, px2 = np.unravel_index(np.nanargmax(s_), s_.shape); cxe, cye = wcen(d_, px2, py2)
        q_in = int((tpf.quality[good][m_in] != 0).sum()) if hasattr(tpf, "quality") else -1
        rows_pt.append(dict(sector=s, epoch=int(ep), n_in=int(m_in.sum()), n_out=int(m_out.sum()), peak_snr=round(float(s_[py2, px2]), 1), peak_dist_target_px=round(float(np.hypot(px2 - icx, py2 - icy)), 2),
                            peak_dist_nb_px=round(float(np.hypot(px2 - nb.px, py2 - nb.py)), 2), cen_offset_px=round(float(np.hypot(cxe - icx, cye - icy)), 3), cen_offset_arcsec=round(float(np.hypot(cxe - icx, cye - icy)) * 21, 1),
                            cen_dist_nb_px=round(float(np.hypot(cxe - nb.px, cye - nb.py)), 2), quality_flagged_in_transit=q_in))
    # ---- quality flags / jitter mask audit near transits; sector-1 bad pointing stretch
    lc = P.clean(build_lc(tpf, aper)); flat_pre = P.flatten(lc, 33); flat = P.mask_jitter(flat_pre, 6, 3)
    for ep in np.unique(np.round((t - T0) / PER).astype(int)):
        tc = T0 + ep * PER; w_raw = np.abs(t - tc) < DUR / 2; w_pre = np.abs(flat_pre.time.value - tc) < DUR / 2; w_post = np.abs(flat.time.value - tc) < DUR / 2
        if w_raw.sum() == 0: continue
        rows_q.append(dict(sector=s, epoch=int(ep), cadences_in_cube=int(w_raw.sum()), after_clean=int(w_pre.sum()), after_jitter_mask=int(w_post.sum()),
                           quality_nonzero=int((tpf.quality[w_raw] != 0).sum()), in_s1_bad_pointing=bool(s == 1 and 1347 < tc < 1349.6)))
    # ---- secondary scan per sector at the phase of the stitched spike, and per-sector phase-0.5 depth
    f_ = flat.flux.value; t_ = flat.time.value; phi = ((t_ - T0) / PER) % 1.0; sig = np.nanstd(f_[np.abs(phase_of(t_)) > DUR])
    for p0 in [0.50, 0.52, 0.54, 0.56]:
        inb = np.abs(((phi - p0 + 0.5) % 1.0) - 0.5) < (DUR / PER) / 2
        rows_sec.append(dict(sector=s, phase=p0, n=int(inb.sum()), depth_ppm=round((1 - np.nanmean(f_[inb])) * 1e6) if inb.sum() else np.nan, err_ppm=round(sig / np.sqrt(max(inb.sum(), 1)) * 1e6),
                             times=",".join(f"{x:.2f}" for x in np.unique(np.round(t_[inb], 1))[:6])))
    # ---- rotation search at short periods (0.2-5 d) where scattered-light systematics are weaker; use flattened? no: flatten kills >1.4 d. Use raw corrected, transits masked
    tl, fl_ = lc.time.value, lc.flux.value; m = np.abs(phase_of(tl)) > DUR; tl, fl_ = tl[m], fl_[m] / np.nanmedian(fl_[m])
    for lo, hi in [(0.2, 2.0), (2.0, 5.0)]:
        fr, pw = LombScargle(tl, fl_).autopower(minimum_frequency=1 / hi, maximum_frequency=1 / lo, samples_per_peak=10); i = np.argmax(pw)
        A_ = np.column_stack([np.sin(2 * np.pi * fr[i] * tl), np.cos(2 * np.pi * fr[i] * tl), np.ones_like(tl)]); cf = np.linalg.lstsq(A_, fl_, rcond=None)[0]
        rows_rot.append(dict(sector=s, range_d=f"{lo}-{hi}", peak_period_d=round(1 / fr[i], 3), amp_ppm=round(np.hypot(cf[0], cf[1]) * 1e6), power=round(float(pw[i]), 3)))
    # ---- flares: upward runs in the flattened curve, and whether any falls inside a transit
    fv = flat.flux.value; mad = np.nanmedian(np.abs(fv - np.nanmedian(fv))) * 1.4826; hi_ = fv > 1 + 4 * mad; runs = hi_ & (np.roll(hi_, 1) | np.roll(hi_, -1))
    starts = np.nonzero(runs & ~np.roll(runs, 1))[0]
    for i0 in starts:
        j = i0
        while j + 1 < len(runs) and runs[j + 1]: j += 1
        rows_fl.append(dict(sector=s, t_start=round(float(flat.time.value[i0]), 3), n_points=int(j - i0 + 1), peak_ppm=round(float(fv[i0:j + 1].max() - 1) * 1e6), in_transit=bool(np.abs(phase_of(flat.time.value[i0])) < DUR / 2),
                            amp_vs_depth=round(float(fv[i0:j + 1].max() - 1) / 0.038, 3)))
    # ---- sector 103 deep dive figure: per-transit diff images, pixel light curves, background
    if s == 103:
        fig103, ax = plt.subplots(2, 4, figsize=(16, 7.5))
        for k, ep in enumerate(epochs[:2]):
            tc = T0 + ep * PER; m_in = np.abs(tt - tc) < DUR / 2; m_out = (np.abs(tt - tc) > DUR) & (np.abs(tt - tc) < 3 * DUR)
            d_ = fb[m_out].mean(axis=0) - fb[m_in].mean(axis=0); n_ = np.sqrt(fb[m_out].var(axis=0) / m_out.sum() + fb[m_in].var(axis=0) / m_in.sum()) + 1e-9
            im = ax[0, k].imshow(d_ / n_, origin="lower", cmap="magma"); ax[0, k].contour(aper, levels=[0.5], colors="cyan", linewidths=0.7); ax[0, k].plot(icx, icy, "g+", ms=10); ax[0, k].plot(nb.px, nb.py, "cx")
            ax[0, k].set_xlim(icx - 4, icx + 4); ax[0, k].set_ylim(icy - 4, icy + 4); ax[0, k].set_title(f"S103 epoch {ep}: diff SNR, n_in={m_in.sum()} n_out={m_out.sum()}", fontsize=9); plt.colorbar(im, ax=ax[0, k], fraction=0.046)
            raw_ = fl[:, aper].sum(axis=1); w = np.abs(tt - tc) < 4 * DUR
            ax[1, k].plot((tt[w] - tc) * 24, raw_[w] / np.nanmedian(raw_[m_out]), "k.", ms=2); ax[1, k].set_title(f"raw aperture flux around epoch {ep}", fontsize=9); ax[1, k].set_xlabel("hours")
            bgm = np.nanmedian(fl[:, ~aper], axis=1); ax[1, k].plot((tt[w] - tc) * 24, bgm[w] / np.nanmedian(bgm[m_out]) - 0.1, "r.", ms=2); ax[1, k].text(0.02, 0.05, "red: median background pixel, -0.1", transform=ax[1, k].transAxes, fontsize=8)
        d_all = fb[oot].mean(axis=0) - fb[in_tr].mean(axis=0); im = ax[0, 2].imshow(d_all, origin="lower", cmap="RdBu_r", vmin=-np.nanpercentile(np.abs(d_all), 99), vmax=np.nanpercentile(np.abs(d_all), 99)); plt.colorbar(im, ax=ax[0, 2], fraction=0.046)
        ax[0, 2].contour(aper, levels=[0.5], colors="k", linewidths=0.7); ax[0, 2].plot(icx, icy, "g+", ms=10); ax[0, 2].set_title("S103 raw difference image (e-/s), both transits", fontsize=9)
        ax[0, 3].imshow(np.log10(np.clip(img, 1, None)), origin="lower", cmap="gray"); ax[0, 3].contour(aper, levels=[0.5], colors="red", linewidths=0.7); ax[0, 3].set_title("S103 mean image (log)", fontsize=9)
        ax[1, 2].plot(tt, bgm, "k.", ms=1); [ax[1, 2].axvline(T0 + e * PER, color="r", lw=0.8) for e in epochs]; ax[1, 2].set_title("S103 median background pixel vs time (red = transits)", fontsize=9)
        ax[1, 3].plot(tt, cx - np.nanmedian(cx), "b.", ms=1, label="x"); ax[1, 3].plot(tt, cy - np.nanmedian(cy), "g.", ms=1, label="y"); [ax[1, 3].axvline(T0 + e * PER, color="r", lw=0.8) for e in epochs]
        ax[1, 3].set_ylim(-0.5, 0.5); ax[1, 3].legend(fontsize=8); ax[1, 3].set_title("S103 frame centroid vs time (px)", fontsize=9)
        plt.tight_layout(); fig103.savefig(os.path.join(OUT, "addendum_sector103.png"), dpi=100); plt.close(fig103)

D = pd.DataFrame(rows_dil); C = pd.DataFrame(rows_cen); PT = pd.DataFrame(rows_pt); SEC = pd.DataFrame(rows_sec); ROT = pd.DataFrame(rows_rot); FL = pd.DataFrame(rows_fl); Q = pd.DataFrame(rows_q)
for name, df in [("dilution_psf", D), ("centroid_vectors", C), ("per_transit_diffimage", PT), ("secondary_per_sector", SEC), ("rotation_short", ROT), ("flares", FL), ("transit_cadence_audit", Q)]:
    df.to_csv(os.path.join(DV, f"{name}.csv"), index=False)

# ---- two-hypothesis fit on the per-sector depths: dip on target (D_true * share_target) vs dip on neighbour (E * share_nb), both apertures jointly
dd = dep.merge(D[["sector", "aperture", "share_target", "share_nb"]], on=["sector", "aperture"])
use = dd[dd.aperture.isin(["pipeline", "target_only", "neighbour"])]
def wfit(x, y, e): w = 1 / e ** 2; k = np.sum(w * x * y) / np.sum(w * x * x); return k, np.sqrt(1 / np.sum(w * x * x)), np.sum(((y - k * x) / e) ** 2)
y, e = use.depth_ppm.values / 1e6, use.depth_err_ppm.values / 1e6
kT, kTe, chiT = wfit(use.share_target.values, y, e); kN, kNe, chiN = wfit(use.share_nb.values, y, e)
# allow a per-sector 10% depth jitter (systematics / spots) in the chi2 as well
e2 = np.hypot(e, 0.10 * y); kT2, _, chiT2 = wfit(use.share_target.values, y, e2); kN2, _, chiN2 = wfit(use.share_nb.values, y, e2)
hyp = dict(n_measurements=int(len(use)), H_target_true_depth_pct=round(kT * 100, 2), H_target_err_pct=round(kTe * 100, 2), H_target_chi2=round(chiT, 1), H_target_chi2_with10pct_jitter=round(chiT2, 1),
           H_neighbour_eclipse_pct=round(kN * 100, 1), H_neighbour_err_pct=round(kNe * 100, 1), H_neighbour_chi2=round(chiN, 1), H_neighbour_chi2_with10pct_jitter=round(chiN2, 1),
           H_neighbour_physically_possible=bool(kN <= 1.0), psf_sigma_px=a.psf_sigma)
use = use.assign(pred_target_ppm=(kT * use.share_target * 1e6).round(), pred_nb_ppm=(kN * use.share_nb * 1e6).round())
use.to_csv(os.path.join(DV, "depth_hypotheses.csv"), index=False); json.dump(hyp, open(os.path.join(DV, "depth_hypotheses.json"), "w"), indent=2)
log("hypotheses: " + json.dumps(hyp)); log("\n" + use.to_string(index=False))
fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
for name, mk in [("pipeline", "o"), ("target_only", "s"), ("neighbour", "^")]:
    u = use[use.aperture == name]; ax[0].errorbar(u.share_target, u.depth_ppm, u.depth_err_ppm, fmt=mk, ms=5, label=name, capsize=2); ax[1].errorbar(u.share_nb, u.depth_ppm, u.depth_err_ppm, fmt=mk, ms=5, label=name, capsize=2)
xx_ = np.linspace(0, 1, 10); ax[0].plot(xx_, kT * xx_ * 1e6, "r-", lw=1, label=f"dip on target: D={kT*100:.2f}%, chi2={chiT:.0f}/{len(use)}"); ax[1].plot(xx_, kN * xx_ * 1e6, "r-", lw=1, label=f"dip on neighbour: E={kN*100:.0f}%, chi2={chiN:.0f}/{len(use)}")
ax[0].set_xlabel("target's share of aperture light (PSF model)"); ax[1].set_xlabel("44-arcsec neighbour's share of aperture light (PSF model)"); ax[0].set_ylabel("fitted depth (ppm)")
for a_ in ax: a_.legend(fontsize=8); a_.set_xlim(0, 1.02)
ax[1].set_ylim(ax[0].get_ylim()); plt.tight_layout(); fig.savefig(os.path.join(OUT, "addendum_blend_hypotheses.png"), dpi=110); plt.close(fig)
# ---- centroid vector figure
fig, ax = plt.subplots(1, 1, figsize=(7, 6.5))
for _, r in C.iterrows():
    ax.errorbar(r.obs_dx * 1e3, r.obs_dy * 1e3, xerr=r.err * 1e3, yerr=r.err * 1e3, fmt="ko", ms=4, capsize=2); ax.annotate(f"S{int(r.sector)}", (r.obs_dx * 1e3, r.obs_dy * 1e3), fontsize=8, xytext=(3, 3), textcoords="offset points")
    ax.plot(r.pred_target_dx * 1e3, r.pred_target_dy * 1e3, "g^", ms=6); ax.plot(r.pred_nb_dx * 1e3, r.pred_nb_dy * 1e3, "cx", ms=7)
ax.plot([], [], "ko", label="observed in-transit centroid shift"); ax.plot([], [], "g^", label="predicted: dip on target"); ax.plot([], [], "cx", label="predicted: dip on 44-arcsec neighbour")
ax.axhline(0, color="0.7", lw=0.5); ax.axvline(0, color="0.7", lw=0.5); ax.set_xlabel("dx (millipixel)"); ax.set_ylabel("dy (millipixel)"); ax.legend(fontsize=8); ax.set_title("centroid shift per sector vs the two blend hypotheses (frame coordinates differ per sector)", fontsize=9)
plt.tight_layout(); fig.savefig(os.path.join(OUT, "addendum_centroid_vectors.png"), dpi=110); plt.close(fig)
pd.set_option("display.width", 250)
log("\n" + C.to_string(index=False)); log("\n" + PT.to_string(index=False)); log("\n" + SEC.to_string(index=False)); log("\n" + ROT.to_string(index=False)); log("\n" + FL.to_string(index=False)); log("\n" + Q.to_string(index=False)); log("done")
