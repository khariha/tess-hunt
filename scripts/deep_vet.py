"""Deep vetting of one candidate from cached FFI cubes: per-sector pixel tests, blend budget,
per-sector depth vs aperture contents, odd/even, secondary at any phase, phase-curve limits,
individual transit times, rotation and flares.

    python -W ignore scripts/deep_vet.py --tic 259719966 --period 11.077207 --t0 3877.82312 --dur-h 2.372 --rp 0.196 --ars 37.84 --b 0.555

Needs data/cubes/s<N>_tic<ID>_15x15.fits (scripts/fetch_cubes.py) and data/gaia_neighbors_<ID>.csv.
Writes data/deepvet_<ID>/*.csv and plots/candidates/submission/s1_tic<ID>/addendum_*.png
"""
import argparse, os, sys, time, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import lightkurve as lk
from scipy.optimize import least_squares
from astropy.timeseries import LombScargle
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tesshunt import pipeline as P
from tesshunt.cubes import load_tpf, ROOT
from tesshunt.transitfit import Model
from tesshunt.sweep_core import _central_aperture

ap = argparse.ArgumentParser()
ap.add_argument("--tic", type=int, required=True); ap.add_argument("--sector", type=int, default=1)
ap.add_argument("--period", type=float, required=True); ap.add_argument("--t0", type=float, required=True)
ap.add_argument("--dur-h", type=float, required=True); ap.add_argument("--rp", type=float, required=True)
ap.add_argument("--ars", type=float, required=True); ap.add_argument("--b", type=float, required=True)
ap.add_argument("--u1", type=float, default=0.355); ap.add_argument("--u2", type=float, default=0.332)
ap.add_argument("--sectors", type=int, nargs="*", default=None)
a = ap.parse_args(); tic = a.tic
PER, T0, DUR, RP, ARS, B, U1, U2 = a.period, a.t0, a.dur_h / 24, a.rp, a.ars, a.b, a.u1, a.u2
OUT = os.path.join(ROOT, "plots", "candidates", "submission", f"s{a.sector}_tic{tic}"); os.makedirs(OUT, exist_ok=True)
DV = os.path.join(ROOT, "data", f"deepvet_{tic}"); os.makedirs(DV, exist_ok=True)
def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)

tg = pd.read_csv(os.path.join(ROOT, "data", f"targets_s{a.sector}.csv")).drop_duplicates("ID").set_index("ID").loc[tic]
gaia = pd.read_csv(os.path.join(ROOT, "data", f"gaia_neighbors_{tic}.csv"))
gaia["light_frac"] = 10 ** (-0.4 * (gaia.phot_g_mean_mag - gaia.phot_g_mean_mag.iloc[0]))
from tess_stars2px import tess_stars2px_function_entry
sectors = a.sectors or sorted(set(int(s) for s in tess_stars2px_function_entry(tic, tg.ra, tg.dec)[3]))

def phase_of(t): return ((t - T0 + PER / 2) % PER) - PER / 2

def build_lc(tpf, aper):
    raw = tpf.to_lightcurve(aperture_mask=aper)
    reg = tpf.flux[:, ~aper]; reg = reg.value if hasattr(reg, "value") else reg
    good = np.isfinite(raw.flux.value) & np.isfinite(reg).all(axis=1)
    dm = lk.DesignMatrix(reg[good], name="bg").pca(5).append_constant()
    return lk.RegressionCorrector(raw[good]).correct(dm)

def fit_depth(t, f, e, exp, label=""):
    """depth with geometry fixed at the joint fit; error from the jacobian scaled by reduced chi2"""
    m = Model(t, exp)
    r = least_squares(lambda th: (f - m(T0, PER, th[0], ARS, B, U1, U2)) / e, [RP], bounds=([0.005], [0.6]))
    J = r.jac; cov = np.linalg.inv(J.T @ J) * max(np.sum(r.fun ** 2) / max(len(f) - 1, 1), 1)
    rp, err = r.x[0], np.sqrt(cov[0, 0])
    return rp ** 2 * 1e6, 2 * rp * err * 1e6, len(f)

def fit_time(t, f, e, exp, tc_guess):
    """mid-time of one transit with shape fixed; returns (tc, err, n, chi2red)"""
    m = Model(t, exp)
    r = least_squares(lambda th: (f - m(th[0], PER, RP, ARS, B, U1, U2)) / e, [tc_guess], bounds=([tc_guess - 0.1], [tc_guess + 0.1]))
    J = r.jac; chi2r = np.sum(r.fun ** 2) / max(len(f) - 1, 1)
    cov = np.linalg.inv(J.T @ J) * max(chi2r, 1)
    return r.x[0], np.sqrt(cov[0, 0]), len(f), chi2r

def diff_centroid(diff, snr_map, cx, cy, rad=1.5):
    """flux-weighted centroid of the difference image within rad px of its peak (positive pixels only)"""
    ny, nx = diff.shape; yy, xx = np.mgrid[0:ny, 0:nx]
    py, px = np.unravel_index(np.nanargmax(snr_map), snr_map.shape)
    w = np.clip(diff, 0, None) * (np.hypot(xx - px, yy - py) <= rad)
    if w.sum() <= 0: return np.nan, np.nan
    return (w * xx).sum() / w.sum(), (w * yy).sum() / w.sum()

pix_rows, dep_rows, tt_rows, act_rows, aper_rows = [], [], [], [], []
flats_pipe, flats_tgt, raws_for_phase = [], [], []
fig_px, axes_px = plt.subplots(3, len(sectors), figsize=(2.6 * len(sectors), 8.2))
for k, s in enumerate(sectors):
    tpf = load_tpf(tic, s, float(tg.ra), float(tg.dec))
    if tpf is None: log(f"sector {s}: no cube"); continue
    flux = tpf.flux.value; t = tpf.time.value
    cad_min = float(np.nanmedian(np.diff(t)) * 1440)
    tx, ty = [float(v[0]) for v in tpf.wcs.all_world2pix([tg.ra], [tg.dec], 0)]
    gx, gy = tpf.wcs.all_world2pix(gaia.ra.values, gaia.dec.values, 0)
    gaia_s = gaia.assign(px=gx, py=gy, sep_px=np.hypot(gx - tx, gy - ty))
    nb = gaia_s.iloc[1:].sort_values("light_frac", ascending=False).iloc[0]  # brightest neighbour
    # --- apertures: pipeline (threshold 5, as the sweep/fit used) and target-only (pipeline minus pixels within 1 px of the bright neighbour)
    try:
        aper = tpf.create_threshold_mask(threshold=5, reference_pixel="center")
        if aper.sum() == 0 or aper.sum() > 0.5 * aper.size: raise ValueError
    except Exception:
        aper = _central_aperture(tpf)[2]
    ny, nx = aper.shape; yy, xx = np.mgrid[0:ny, 0:nx]
    near_nb = np.hypot(xx - nb.px, yy - nb.py) <= 1.0
    aper_tgt = aper & ~near_nb
    if aper_tgt.sum() == 0: aper_tgt = (np.hypot(xx - tx, yy - ty) <= 1.0)
    aper_nb = np.hypot(xx - nb.px, yy - nb.py) <= 1.0            # neighbour-only aperture
    # which Gaia sources fall in each aperture (pixel containing the source)
    def inside(ap_):
        ix, iy = np.round(gaia_s.px).astype(int), np.round(gaia_s.py).astype(int)
        ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
        m = np.zeros(len(gaia_s), bool); m[ok] = ap_[iy[ok], ix[ok]]; return m
    in_pipe, in_tgt = inside(aper), inside(aper_tgt)
    dil_pipe = gaia_s.light_frac[in_pipe & (gaia_s.index > 0)].sum(); dil_tgt = gaia_s.light_frac[in_tgt & (gaia_s.index > 0)].sum()
    aper_rows.append(dict(sector=s, cadence_min=round(cad_min, 2), camera=tpf.camera, ccd=tpf.ccd, aper_px=int(aper.sum()), aper_tgt_px=int(aper_tgt.sum()),
                          nb_sep_px=round(float(nb.sep_px), 2), nb_in_pipeline_aper=bool(in_pipe[nb.name]), nb_dist_to_aper_px=round(float(np.min(np.hypot(xx[aper] - nb.px, yy[aper] - nb.py))), 2),
                          gaia_in_pipeline_aper=",".join(f"{str(sid)[-6:]}:{lf*100:.1f}%" for sid, lf in zip(gaia_s.source_id[in_pipe], gaia_s.light_frac[in_pipe])),
                          dilution_pipeline=round(float(dil_pipe), 4), dilution_target_only=round(float(dil_tgt), 4), tic_contratio=float(tg.contratio)))
    # --- light curves
    lcs = {}
    for name, ap_ in [("pipeline", aper), ("target_only", aper_tgt), ("neighbour", aper_nb)]:
        lc = build_lc(tpf, ap_); lc = P.clean(lc)
        flat = P.mask_jitter(P.flatten(lc, 33), 6, 3)
        flat.meta["SECTOR"] = s; lcs[name] = (lc, flat)
    flats_pipe.append(lcs["pipeline"][1]); flats_tgt.append(lcs["target_only"][1])
    # --- centroid / difference image (pipeline aperture for the target position)
    c = P.centroid_test(tpf, PER, T0, DUR, aperture=None)   # target = cutout centre pixel (7,7); report offsets from the WCS target position below
    snr_map, diff = c["diff_snr_map"], c.get("diff_image", c["diff_snr_map"])
    py_, px_ = np.unravel_index(np.nanargmax(snr_map), snr_map.shape)
    dcx, dcy = diff_centroid(diff, snr_map, tx, ty)
    # jackknife over in-transit cadences for the diff-image centroid error
    good = np.isfinite(flux).all(axis=(1, 2)); fl, tt = flux[good], t[good]
    ph = phase_of(tt); in_tr = np.abs(ph) < DUR / 2; oot = (np.abs(ph) > DUR) & (np.abs(ph) < 3 * DUR)
    fb = P._subtract_plane_background(fl); fb = np.clip(fb, 0, None)
    jk = []
    idx_in = np.nonzero(in_tr)[0]
    for drop in idx_in[:: max(1, len(idx_in) // 20)]:
        m_in = in_tr.copy(); m_in[drop] = False
        d_ = fb[oot].mean(axis=0) - fb[m_in].mean(axis=0)
        n_ = np.sqrt(fb[oot].var(axis=0) / oot.sum() + fb[m_in].var(axis=0) / m_in.sum()) + 1e-9
        jk.append(diff_centroid(d_, d_ / n_, tx, ty))
    jk = np.array(jk); njk = len(jk)
    jk_err = np.sqrt((njk - 1) / njk * ((jk - jk.mean(axis=0)) ** 2).sum(axis=0)) if njk > 2 else np.array([np.nan, np.nan])
    off_vec = np.array([dcx - tx, dcy - ty]); off_px = float(np.hypot(*off_vec)); off_err = float(np.hypot(*jk_err))
    nb_vec = np.array([nb.px - tx, nb.py - ty])
    cos_to_nb = float(off_vec @ nb_vec / (np.linalg.norm(off_vec) * np.linalg.norm(nb_vec) + 1e-12))
    dist_peak_to_nb = float(np.hypot(px_ - nb.px, py_ - nb.py)); dist_peak_to_tgt = float(np.hypot(px_ - tx, py_ - ty))
    # centroid-shift direction (in-transit minus out-of-transit flux-weighted centroid of the whole frame)
    cx_, cy_ = c["centroid_x"], c["centroid_y"]; it_, io_ = c["in_transit"], (np.abs(c["phase"]) > DUR) & (np.abs(c["phase"]) < 3 * DUR)
    sh_vec = np.array([cx_[it_].mean() - cx_[io_].mean(), cy_[it_].mean() - cy_[io_].mean()])
    cos_shift_nb = float(sh_vec @ nb_vec / (np.linalg.norm(sh_vec) * np.linalg.norm(nb_vec) + 1e-12))
    # snr at the target pixel and at the neighbour pixel
    snr_tgt = float(snr_map[int(round(ty)), int(round(tx))]); snr_nb = float(snr_map[int(round(nb.py)), int(round(nb.px))]) if 0 <= round(nb.py) < ny and 0 <= round(nb.px) < nx else np.nan
    pix_rows.append(dict(sector=s, cadence_min=round(cad_min, 2), n_in=c["n_in"], n_out=c["n_out"], shift_px=round(c["shift_px"], 4), shift_sigma=round(c["shift_sigma"], 1),
                         shift_dir_cos_to_neighbour=round(cos_shift_nb, 2),
                         peak_dist_target_px=round(dist_peak_to_tgt, 2), peak_dist_neighbour_px=round(dist_peak_to_nb, 2), peak_snr=round(c["diff_peak_snr"], 1),
                         snr_at_target_px=round(snr_tgt, 1), snr_at_neighbour_px=round(snr_nb, 1),
                         diffcen_offset_px=round(off_px, 3), diffcen_offset_err_px=round(off_err, 3), diffcen_offset_arcsec=round(off_px * 21, 1),
                         diffcen_dir_cos_to_neighbour=round(cos_to_nb, 2), neighbour_sep_px=round(float(nb.sep_px), 2)))
    log(f"sector {s} ({cad_min:.0f}-min): shift {c['shift_px']:.4f} px={c['shift_sigma']:.1f}sig cos_nb={cos_shift_nb:+.2f} | peak {dist_peak_to_tgt:.2f} px from target, {dist_peak_to_nb:.2f} from nb, {c['diff_peak_snr']:.0f} sig | diff centroid {off_px:.3f}+/-{off_err:.3f} px ({off_px*21:.0f}\") cos_nb={cos_to_nb:+.2f}")
    # figure column: mean image, diff snr, fold
    axm = axes_px[0, k]; axm.imshow(np.log10(np.clip(np.nanmean(fl, axis=0), 1, None)), origin="lower", cmap="gray"); axm.contour(aper, levels=[0.5], colors="red", linewidths=0.8)
    for _, n in gaia_s[gaia_s.sep_px < 4].iterrows(): axm.plot(n.px, n.py, "o", mfc="none", mec="cyan", ms=4 + 12 * min(n.light_frac, 1), mew=0.8)
    axm.set_title(f"S{s} {cad_min:.0f}-min", fontsize=9); axm.set_xticks([]); axm.set_yticks([])
    axd = axes_px[1, k]; axd.imshow(snr_map, origin="lower", cmap="magma"); axd.contour(aper, levels=[0.5], colors="cyan", linewidths=0.6)
    axd.plot(tx, ty, "g+", ms=10, mew=1.5); axd.plot(nb.px, nb.py, "cx", ms=7); axd.plot(dcx, dcy, "wo", mfc="none", ms=8)
    axd.set_xlim(tx - 4, tx + 4); axd.set_ylim(ty - 4, ty + 4); axd.set_xticks([]); axd.set_yticks([])
    axd.set_title(f"diff: {off_px*21:.0f}\" off, {c['diff_peak_snr']:.0f}σ", fontsize=9)
    axf = axes_px[2, k]; fp = lcs["pipeline"][1].fold(period=PER, epoch_time=T0); axf.plot(fp.time.value * 24, fp.flux.value, ".", ms=2, alpha=0.3, color="0.5")
    fn_ = lcs["neighbour"][1].fold(period=PER, epoch_time=T0); bn = fn_.bin(time_bin_size=DUR / 6); axf.plot(bn.time.value * 24, bn.flux.value - 0.08, "c.-", ms=3, lw=1)
    bp = fp.bin(time_bin_size=DUR / 6); axf.plot(bp.time.value * 24, bp.flux.value, "r.-", ms=3, lw=1)
    axf.set_xlim(-3 * DUR * 24, 3 * DUR * 24); axf.set_ylim(0.86, 1.03); axf.set_xticks([-6, 0, 6]); axf.tick_params(labelsize=7)
    if k: axf.set_yticks([])
    axf.set_title("fold: pipeline (red), nb aperture (cyan, -0.08)", fontsize=7)
    # --- per-sector depth (both apertures) and a 'neighbour aperture' depth
    for name in ["pipeline", "target_only", "neighbour"]:
        flat = lcs[name][1]; tt_, ff_ = flat.time.value, flat.flux.value; ph_ = phase_of(tt_)
        keep = np.abs(ph_) < 3 * DUR; oo = keep & (np.abs(ph_) > DUR)
        sig = np.nanstd(ff_[oo]); exp_ = np.full(keep.sum(), np.nanmedian(np.diff(tt_)))
        d, de, n = fit_depth(tt_[keep], ff_[keep], np.full(keep.sum(), sig), exp_)
        dep_rows.append(dict(sector=s, aperture=name, depth_ppm=round(d), depth_err_ppm=round(de), n_points=n, scatter_ppm=round(sig * 1e6),
                             dilution=round(float(dil_pipe if name == "pipeline" else dil_tgt if name == "target_only" else np.nan), 4)))
    # --- individual transit times (pipeline aperture)
    flat = lcs["pipeline"][1]; tt_, ff_ = flat.time.value, flat.flux.value
    epochs = np.unique(np.round((tt_ - T0) / PER).astype(int))
    for ep in epochs:
        tc = T0 + ep * PER; w = np.abs(tt_ - tc) < 2.5 * DUR
        if (np.abs(tt_ - tc) < DUR / 2).sum() < max(3, int(0.5 * DUR / (cad_min / 1440))): continue   # need at least half the transit covered
        oo = w & (np.abs(tt_ - tc) > DUR); sig = np.nanstd(ff_[oo]) if oo.sum() > 5 else np.nanstd(ff_)
        tcf, tce, n, chi2r = fit_time(tt_[w], ff_[w], np.full(w.sum(), sig), np.full(w.sum(), cad_min / 1440), tc)
        d, de, _ = fit_depth(tt_[w], ff_[w], np.full(w.sum(), sig), np.full(w.sum(), cad_min / 1440))
        tt_rows.append(dict(sector=s, epoch=int(ep), tc_pred=round(tc, 5), tc_fit=round(tcf, 5), oc_min=round((tcf - tc) * 1440, 2), err_min=round(tce * 1440, 2),
                            n_in_transit=int((np.abs(tt_ - tc) < DUR / 2).sum()), depth_ppm=round(d), depth_err_ppm=round(de), chi2red=round(chi2r, 2)))
    # --- activity: rotation (regression-corrected, unflattened, transits masked, per-sector median normalised) and flares
    lc = lcs["pipeline"][0]; tl, fl_ = lc.time.value, lc.flux.value; m_oot = np.abs(phase_of(tl)) > DUR
    tl, fl_ = tl[m_oot], fl_[m_oot] / np.nanmedian(fl_[m_oot])
    freq, pw = LombScargle(tl, fl_).autopower(minimum_frequency=1 / 15.0, maximum_frequency=1 / 0.2, samples_per_peak=10)
    i = np.argmax(pw); prot = 1 / freq[i]
    A = np.column_stack([np.sin(2 * np.pi * freq[i] * tl), np.cos(2 * np.pi * freq[i] * tl), np.ones_like(tl)]); cf = np.linalg.lstsq(A, fl_, rcond=None)[0]
    amp = np.hypot(cf[0], cf[1]) * 1e6
    # flares: >= 2 consecutive points > 4 sigma above the flattened curve (200-s/10-min only; 30-min can't resolve them)
    flat = lcs["pipeline"][1]; fv = flat.flux.value; mad = np.nanmedian(np.abs(fv - np.nanmedian(fv))) * 1.4826
    hi = fv > 1 + 4 * mad; runs = hi & (np.roll(hi, 1) | np.roll(hi, -1)); nfl = int(np.sum(runs & ~np.roll(runs, 1)))
    span = tl.max() - tl.min()
    act_rows.append(dict(sector=s, cadence_min=round(cad_min, 2), ls_peak_period_d=round(prot, 3), ls_power=round(float(pw[i]), 3), ls_amp_ppm=round(amp), fap=float(LombScargle(tl, fl_).false_alarm_probability(pw[i])),
                         flare_events=nfl, flare_rate_per_day=round(nfl / span, 3), oot_scatter_ppm=round(np.nanstd(fv) * 1e6), span_d=round(span, 1)))
    # keep the unflattened, transit-masked, per-sector normalised series for the phase-curve fit
    raws_for_phase.append(pd.DataFrame(dict(t=tl, f=fl_, sector=s)))
plt.tight_layout(); fig_px.savefig(os.path.join(OUT, "addendum_pixels_all_sectors.png"), dpi=110); plt.close(fig_px)
pd.DataFrame(pix_rows).to_csv(os.path.join(DV, "pixel_tests.csv"), index=False)
pd.DataFrame(aper_rows).to_csv(os.path.join(DV, "aperture_contents.csv"), index=False)
pd.DataFrame(dep_rows).to_csv(os.path.join(DV, "depth_per_sector.csv"), index=False)
pd.DataFrame(tt_rows).to_csv(os.path.join(DV, "transit_times.csv"), index=False)
pd.DataFrame(act_rows).to_csv(os.path.join(DV, "activity.csv"), index=False)

# ---------------- blend budget: eclipse depth each Gaia source would need (G-band flux ratios) ----------------
dep = pd.DataFrame(dep_rows); d_obs = float(dep[(dep.aperture == "pipeline") & (dep.dilution < 0.01)].depth_ppm.mean()) / 1e6 if (dep.dilution < 0.01).any() else RP ** 2
rows = []
for i, n in gaia.iterrows():
    if i == 0: continue
    frac_of_target = n.light_frac
    # if that star is the eclipsing one and all of its light sits in the aperture with the target: needed eclipse depth = d_obs * (1 + f) / f
    needed = d_obs * (1 + frac_of_target) / frac_of_target
    rows.append(dict(source_id=n.source_id, sep_arcsec=round(n.sep_arcsec, 1), sep_px=round(n.sep_arcsec / 21, 2), G=round(n.phot_g_mean_mag, 2), bp_rp=round(n.bp_rp, 2) if np.isfinite(n.bp_rp) else np.nan,
                     parallax_mas=round(n.parallax, 2) if np.isfinite(n.parallax) else np.nan, flux_frac_of_target=round(frac_of_target, 4),
                     needed_eclipse_depth_pct=round(needed * 100, 1), physically_possible=bool(needed <= 1.0), needs_total_eclipse=bool(0.5 < needed <= 1.0)))
bl = pd.DataFrame(rows); bl.to_csv(os.path.join(DV, "blend_budget.csv"), index=False)
log("blend budget:\n" + bl.to_string(index=False))

# ---------------- stitched tests: odd/even, secondary at phase 0.5 and any phase, pipeline vs target-only aperture ----------------
def stitch(flats):
    st = lk.LightCurveCollection(flats).stitch(); st = st[np.isfinite(st.flux.value)]
    return st.time.value, st.flux.value, np.concatenate([np.full(len(f), np.nanmedian(np.diff(f.time.value))) for f in flats]), np.concatenate([np.full(len(f), f.meta["SECTOR"]) for f in flats])
res_stitched = {}
fig_s, axs = plt.subplots(2, 2, figsize=(13, 8))
for j, (name, flats) in enumerate([("pipeline", flats_pipe), ("target_only", flats_tgt)]):
    t, f, exp, sec = stitch(flats); ph = phase_of(t)
    # per-point error = per-sector out-of-transit scatter
    e = np.empty_like(f)
    for s in np.unique(sec):
        m = sec == s; e[m] = np.nanstd(f[m & (np.abs(ph) > DUR)])
    near = np.abs(ph) < 3 * DUR
    ep = np.round((t - T0) / PER).astype(int)
    d_all = fit_depth(t[near], f[near], e[near], exp[near]); d_odd = fit_depth(t[near & (ep % 2 == 1)], f[near & (ep % 2 == 1)], e[near & (ep % 2 == 1)], exp[near & (ep % 2 == 1)])
    d_even = fit_depth(t[near & (ep % 2 == 0)], f[near & (ep % 2 == 0)], e[near & (ep % 2 == 0)], exp[near & (ep % 2 == 0)])
    # secondary scan: weighted mean in a box of the transit duration at every phase, away from the transit
    phases = np.arange(0.06, 0.95, 0.005); dep_scan = []; err_scan = []
    phi = ((t - T0) / PER) % 1.0
    w = 1 / e ** 2
    for p0 in phases:
        dphi = np.abs(((phi - p0 + 0.5) % 1.0) - 0.5); inb = dphi < (DUR / PER) / 2
        if inb.sum() < 20: dep_scan.append(np.nan); err_scan.append(np.nan); continue
        mu_in = np.sum(w[inb] * f[inb]) / np.sum(w[inb]); dep_scan.append((1 - mu_in) * 1e6); err_scan.append(np.sqrt(1 / np.sum(w[inb])) * 1e6)
    dep_scan, err_scan = np.array(dep_scan), np.array(err_scan)
    # the scan's own scatter is the honest noise floor (red noise); quote both
    sig_scan = np.nanstd(dep_scan); i5 = np.nanargmin(np.abs(phases - 0.5)); imax = np.nanargmax(dep_scan / err_scan)
    res_stitched[name] = dict(n_points=int(len(t)), n_transits=int(len(np.unique(ep[np.abs(ph) < DUR / 2]))),
                              depth_all_ppm=round(d_all[0]), depth_all_err=round(d_all[1]), depth_odd_ppm=round(d_odd[0]), depth_odd_err=round(d_odd[1]),
                              depth_even_ppm=round(d_even[0]), depth_even_err=round(d_even[1]), odd_minus_even_ppm=round(d_odd[0] - d_even[0]), odd_minus_even_sigma=round(abs(d_odd[0] - d_even[0]) / np.hypot(d_odd[1], d_even[1]), 2),
                              secondary_phase05_ppm=round(dep_scan[i5]), secondary_phase05_err_ppm=round(err_scan[i5]), secondary_scan_rms_ppm=round(sig_scan),
                              secondary_scan_max_ppm=round(dep_scan[imax]), secondary_scan_max_phase=round(phases[imax], 3), secondary_scan_max_sigma=round(dep_scan[imax] / err_scan[imax], 2),
                              secondary_3sigma_limit_ppm=round(3 * max(sig_scan, np.nanmedian(err_scan))))
    log(f"{name}: " + json.dumps(res_stitched[name]))
    ax = axs[0, j]; ax.plot(ph * 24, f, ".", ms=1.5, alpha=0.2, color="0.6")
    bins = np.linspace(-3 * DUR, 3 * DUR, 61); idx = np.digitize(ph[near], bins); bm = [np.nanmean(f[near][idx == i]) if (idx == i).sum() > 2 else np.nan for i in range(1, len(bins))]
    ax.plot(((bins[:-1] + bins[1:]) / 2) * 24, bm, "ko", ms=3)
    tt_ = np.linspace(-3 * DUR, 3 * DUR, 800); mm = Model(tt_ + T0, np.full(len(tt_), 120 / 86400)); ax.plot(tt_ * 24, mm(T0, PER, np.sqrt(d_all[0] / 1e6), ARS, B, U1, U2), "r-", lw=1.5)
    ax.set_xlim(-3 * DUR * 24, 3 * DUR * 24); ax.set_ylim(0.93, 1.03); ax.set_xlabel("hours from mid-transit"); ax.set_title(f"{name} aperture: depth {d_all[0]:.0f}+/-{d_all[1]:.0f} ppm; odd-even {d_odd[0]-d_even[0]:+.0f}+/-{np.hypot(d_odd[1], d_even[1]):.0f}", fontsize=10)
    ax = axs[1, j]; ax.errorbar(phases, dep_scan, err_scan, fmt=".", ms=3, color="k", elinewidth=0.5); ax.axhline(0, color="r", lw=0.8); ax.axvline(0.5, color="b", ls="--", lw=0.8)
    ax.set_xlabel("orbital phase"); ax.set_ylabel("box depth (ppm)"); ax.set_title(f"secondary scan: phase 0.5 = {dep_scan[i5]:+.0f}+/-{err_scan[i5]:.0f} ppm; scan rms {sig_scan:.0f} ppm; max {dep_scan[imax]:+.0f} at {phases[imax]:.2f}", fontsize=10)
plt.tight_layout(); fig_s.savefig(os.path.join(OUT, "addendum_oddeven_secondary.png"), dpi=110); plt.close(fig_s)
json.dump(res_stitched, open(os.path.join(DV, "stitched_tests.json"), "w"), indent=2)

# ---------------- phase curve: sinusoids at P (beaming/reflection) and P/2 (ellipsoidal) on the unflattened, transit-masked series ----------------
# per sector: cubic-in-time per orbit half (removes scattered light / drift) fitted jointly with the four sinusoid terms shared across sectors
R = pd.concat(raws_for_phase); t, f, sec = R.t.values, R.f.values, R.sector.values
cols = [np.sin(2 * np.pi * (t - T0) / PER), np.cos(2 * np.pi * (t - T0) / PER), np.sin(4 * np.pi * (t - T0) / PER), np.cos(4 * np.pi * (t - T0) / PER)]
blocks = []
for s in np.unique(sec):
    m = sec == s; ts = t[m]; mid = np.nanmedian(ts)
    for half in [ts <= mid, ts > mid]:
        x = (ts - ts.mean()) / 10.0
        for deg in range(4):
            col = np.zeros(len(t)); col[np.nonzero(m)[0][half]] = x[half] ** deg; blocks.append(col)
A = np.column_stack(cols + blocks)
# iterate with 4-sigma clipping
keep = np.ones(len(t), bool)
for _ in range(3):
    cf, *_ = np.linalg.lstsq(A[keep], f[keep], rcond=None); r = f - A @ cf; sd = np.nanstd(r[keep]); keep = np.abs(r) < 4 * sd
cov = np.linalg.inv(A[keep].T @ A[keep]) * sd ** 2
# scale errors for red noise: binned residual rms at 1-day bins vs white expectation
tb = np.floor(t[keep]); rb = pd.Series(r[keep]).groupby(tb).agg(["mean", "count"]); white = sd / np.sqrt(rb["count"]); beta = float(np.sqrt(np.nanmean((rb["mean"] / white) ** 2)))
err = np.sqrt(np.diag(cov))[:4] * beta
A1 = np.hypot(cf[0], cf[1]) * 1e6; A2 = np.hypot(cf[2], cf[3]) * 1e6
phase_res = dict(sin_P_ppm=round(cf[0] * 1e6, 1), cos_P_ppm=round(cf[1] * 1e6, 1), sin_P2_ppm=round(cf[2] * 1e6, 1), cos_P2_ppm=round(cf[3] * 1e6, 1),
                 err_each_ppm=[round(x * 1e6, 1) for x in err], red_noise_beta=round(beta, 2), amp_P_ppm=round(A1, 1), amp_P2_ppm=round(A2, 1),
                 limit_3sigma_P_ppm=round(3 * np.mean(err[:2]) * 1e6, 1), limit_3sigma_P2_ppm=round(3 * np.mean(err[2:]) * 1e6, 1), n_points=int(keep.sum()),
                 note="beaming = sin(P) term, reflection/thermal = -cos(P), ellipsoidal = -cos(P/2); per-sector cubic per orbit-half absorbs trends, so power at P is partly absorbed too: limits are optimistic by up to ~2x")
json.dump(phase_res, open(os.path.join(DV, "phase_curve.json"), "w"), indent=2); log("phase curve: " + json.dumps(phase_res))

# ---------------- figures: depth vs dilution, O-C, rotation ----------------
dep = pd.DataFrame(dep_rows); tt = pd.DataFrame(tt_rows); act = pd.DataFrame(act_rows); apc = pd.DataFrame(aper_rows)
fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
for name, mk, col in [("pipeline", "o", "k"), ("target_only", "s", "tab:blue"), ("neighbour", "^", "tab:cyan")]:
    d = dep[dep.aperture == name]; ax[0].errorbar(d.sector.astype(str), d.depth_ppm, d.depth_err_ppm, fmt=mk, color=col, label=name, ms=5, capsize=2)
d0 = res_stitched["target_only"]["depth_all_ppm"]
ax[0].axhline(d0, color="tab:blue", lw=0.8, ls="--"); ax[0].plot(apc.sector.astype(str), d0 / (1 + apc.dilution_pipeline), "kx", ms=8, label="undiluted depth / (1+dilution)")
ax[0].set_ylabel("fitted depth (ppm)"); ax[0].set_xlabel("sector"); ax[0].legend(fontsize=8); ax[0].set_title("depth per sector and aperture", fontsize=10)
if len(tt):
    ax[1].errorbar(tt.epoch, tt.oc_min, tt.err_min, fmt="o", ms=4, color="k", capsize=2)
    for s in tt.sector.unique(): m = tt.sector == s; ax[1].annotate(f"S{s}", (tt.epoch[m].mean(), tt.oc_min[m].max() + 1), ha="center", fontsize=7)
    ax[1].axhline(0, color="r", lw=0.8); ax[1].set_xlabel("epoch"); ax[1].set_ylabel("O - C (min)"); chi2 = np.sum((tt.oc_min / tt.err_min) ** 2)
    ax[1].set_title(f"transit times: rms {np.std(tt.oc_min):.1f} min, chi2/N = {chi2/len(tt):.2f} (N={len(tt)})", fontsize=10)
ax[2].bar(act.sector.astype(str), act.ls_amp_ppm, color="0.6"); ax[2].set_ylabel("LS peak amplitude (ppm)"); ax[2].set_xlabel("sector")
for i, r_ in act.iterrows(): ax[2].annotate(f"{r_.ls_peak_period_d:.1f} d", (i, r_.ls_amp_ppm), ha="center", va="bottom", fontsize=7)
ax[2].set_title("out-of-transit periodicity per sector (0.2-15 d)", fontsize=10)
plt.tight_layout(); fig.savefig(os.path.join(OUT, "addendum_depth_timing_activity.png"), dpi=110); plt.close(fig)
log("tables:\n" + apc.to_string(index=False) + "\n\n" + dep.pivot(index="sector", columns="aperture", values=["depth_ppm", "depth_err_ppm"]).to_string() + "\n\n" + tt.to_string(index=False) + "\n\n" + act.to_string(index=False))
log("done")
