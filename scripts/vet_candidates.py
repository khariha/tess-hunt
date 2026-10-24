"""Vet every candidate row in data/sweep_s<N>.csv and write a one-page report each.

For each candidate:
    1. catalog check  : is the TIC already a TOI or a known planet host? (Exoplanet Archive TAP)
    2. neighbors      : Gaia DR3 within 2 pixels, light fraction, and the eclipse depth each
                        neighbor would need to fake the dip
    3. centroid test  : does the center of light move during the dip? where does the
                        difference image peak?
    4. verdict        : on_target / blended / off_target / anomaly, with the reasons

Outputs:
    data/vetting_s<N>.csv                       one row per candidate
    plots/candidates/reports/s<N>_tic<ID>.png   light curve, fold, image, difference image
    plots/candidates/reports/s<N>_tic<ID>.md    the report

Run from the project root (network, no paid APIs):
    python -W ignore scripts/vet_candidates.py --sector 1
    python -W ignore scripts/vet_candidates.py --sector 1 --tic 389723320
"""
import argparse, io, os, socket, sys, time
import numpy as np
import pandas as pd
import requests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tesshunt import pipeline as P
from tesshunt.sweep_core import get_tpf

socket.setdefaulttimeout(120)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
REP = os.path.join(ROOT, "plots", "candidates", "reports")
os.makedirs(REP, exist_ok=True)
TAP = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"

ap = argparse.ArgumentParser()
ap.add_argument("--sector", type=int, default=1)
ap.add_argument("--tic", type=int, default=None)
ap.add_argument("--cutout", type=int, default=15)
args = ap.parse_args()
S = args.sector

def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)

def tap(q):
    r = requests.get(TAP, params={"query": q, "format": "csv"}, timeout=180); r.raise_for_status()
    return pd.read_csv(io.StringIO(r.text))

def catalog_check(tic):
    toi = tap(f"select toi, tfopwg_disp, pl_orbper, pl_trandep from toi where tid = {tic}")
    try:
        ps = tap(f"select pl_name, discoverymethod, disc_year from ps where tic_id = 'TIC {tic}'")
    except Exception:
        ps = pd.DataFrame()
    return toi, ps

def gaia_neighbors(ra, dec, tpf, radius_arcsec=45):
    from astroquery.gaia import Gaia
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    Gaia.ROW_LIMIT = 50
    g = Gaia.cone_search_async(SkyCoord(ra, dec, unit="deg"), radius=radius_arcsec * u.arcsec).get_results().to_pandas()
    px, py = tpf.wcs.all_world2pix(g.ra.values, g.dec.values, 0)
    tx, ty = tpf.wcs.all_world2pix([ra], [dec], 0)
    g["dx_px"], g["dy_px"] = px - float(tx[0]), py - float(ty[0])
    g["sep_px"] = np.hypot(g.dx_px, g.dy_px)
    g["sep_arcsec"] = g.sep_px * 21.0
    g = g.sort_values("sep_px").reset_index(drop=True)
    target_G = g.loc[0, "phot_g_mean_mag"]
    g["light_frac"] = 10 ** (-0.4 * (g.phot_g_mean_mag - target_G))   # relative to the target, G band
    return g[["source_id", "phot_g_mean_mag", "parallax", "dx_px", "dy_px", "sep_px", "sep_arcsec", "light_frac"]], (float(tx[0]), float(ty[0]))

def verdict(row, c, g):
    reasons = []
    depth = row.depth_ppm / 1e6
    near = g[(g.sep_px > 0.2) & (g.sep_px <= 2.0)]
    shift_ok = np.isfinite(c["shift_sigma"]) and c["shift_sigma"] < 3
    peak_ok = np.isfinite(c["diff_peak_dpx"]) and c["diff_peak_dpx"] <= 0.75
    # expected centroid shift if the dip is on the target and neighbors stay constant
    exp_shift = float((depth * (near.light_frac * near.sep_px).sum()) / (1 + near.light_frac.sum())) if len(near) else 0.0
    unresolved = g[(g.sep_px <= 0.2) & (g.index > 0)]
    if len(unresolved) and unresolved.light_frac.max() > 0.1:
        reasons.append(f"unresolved companion within 0.2 px ({unresolved.light_frac.max()*100:.0f}% of light): TESS cannot tell which star dips")
    capable = near[(depth / near.light_frac) < 0.6] if len(near) else near
    if not peak_ok and c["diff_peak_snr"] > 4:
        reasons.append(f"difference image peaks {c['diff_peak_dpx']:.1f} px off target at {c['diff_peak_snr']:.0f} sigma")
        v = "off_target"
    elif not shift_ok and c["shift_px"] > 5 * max(exp_shift, 1e-4):
        reasons.append(f"centroid shifts {c['shift_px']:.3f} px ({c['shift_sigma']:.0f} sigma), {c['shift_px']/max(exp_shift,1e-4):.0f}x what an on-target dip predicts")
        v = "anomaly"
    elif len(unresolved) and unresolved.light_frac.max() > 0.1:
        v = "blended"
    elif not np.isfinite(c["diff_peak_snr"]):
        v = "untestable"
        reasons.append("dip occupies too much of the period to define out-of-transit frames (rotator-like)")
    elif c["diff_peak_snr"] < 4 and shift_ok:
        v = "likely_on_target"
        reasons.append(f"difference image too weak to locate the dip ({c['diff_peak_snr']:.0f} sigma); centroid shift {c['shift_sigma']:.1f} sigma is consistent with on-target")
    else:
        v = "on_target"
        reasons.append(f"centroid shift {c['shift_sigma']:.1f} sigma; difference image peak {c['diff_peak_dpx']:.1f} px from target at {c['diff_peak_snr']:.0f} sigma")
    if len(capable):
        reasons.append(f"{len(capable)} neighbor(s) within 2 px bright enough to fake the dip with a <60% eclipse")
    if row.secondary_ppm > 2 * row.depth_err_ppm:
        reasons.append(f"possible secondary eclipse {row.secondary_ppm:.0f} ppm ({row.secondary_ppm/row.depth_err_ppm:.1f} sigma)")
    if row.period_d < 0.5:
        reasons.append("period under 12 h: check for a contact binary at 2x period")
    return v, reasons, exp_shift

sw = pd.read_csv(os.path.join(DATA, f"sweep_s{S}.csv"))
tg = pd.read_csv(os.path.join(DATA, f"targets_s{S}.csv")).drop_duplicates("ID").set_index("ID")
cands = sw[sw.candidate == True]
if args.tic:
    cands = cands[cands.ID == args.tic]
log(f"{len(cands)} candidates to vet")
out_rows = []
for row in cands.itertuples(index=False):
    tic = int(row.ID); t = tg.loc[tic]
    log(f"TIC {tic}")
    try:
        toi, ps = catalog_check(tic)
        tpf = get_tpf(t.ra, t.dec, S, args.cutout, "s3")
        lc, raw, aper = P.lightcurve_from_cutout(tpf, threshold=5)
        flat = P.mask_jitter(P.flatten(P.clean(lc), 33), 6, 3)
        c = P.centroid_test(tpf, row.period_d, row.t0, row.duration_h / 24, aperture=aper)
        g, (tx, ty) = gaia_neighbors(t.ra, t.dec, tpf)
        v, reasons, exp_shift = verdict(row, c, g)
        # ---- figure ----
        fig, ax = plt.subplots(2, 2, figsize=(13, 9))
        flat.plot(ax=ax[0, 0]); ax[0, 0].set_title(f"TIC {tic}  Tmag {t.Tmag:.1f}  R={t.rad:.2f}  Teff={t.Teff:.0f}  d={t.d:.0f} pc")
        folded = P.fold(flat, row.period_d, row.t0)
        folded.scatter(ax=ax[0, 1], s=4, alpha=0.5)
        folded.bin(time_bin_size=max(row.duration_h / 24 / 8, 0.005)).plot(ax=ax[0, 1], color="red", lw=2)
        ax[0, 1].set_xlim(-2 * row.duration_h / 24, 2 * row.duration_h / 24)
        ax[0, 1].set_title(f"P={row.period_d:.4f} d  depth={row.depth_ppm:.0f} ppm  dur={row.duration_h:.1f} h  SNR={row.snr:.0f}  Rp~{row.planet_re:.1f} Re")
        mean = np.nanmean(tpf.flux.value, axis=0)
        ax[1, 0].imshow(np.log10(np.clip(mean, 1, None)), origin="lower", cmap="gray")
        ax[1, 0].contour(aper, levels=[0.5], colors="red", linewidths=0.8)
        for _, n in g[g.sep_px < 3].iterrows():
            ax[1, 0].plot(tx + n.dx_px, ty + n.dy_px, "o", mfc="none", mec="cyan", ms=6 + 10 * min(n.light_frac, 1))
        ax[1, 0].set_title("average image (log); cyan = Gaia stars, size = brightness")
        im = ax[1, 1].imshow(c["diff_snr_map"], origin="lower", cmap="magma"); plt.colorbar(im, ax=ax[1, 1], fraction=0.046)
        ax[1, 1].contour(aper, levels=[0.5], colors="cyan", linewidths=0.8)
        ax[1, 1].set_title(f"where the dip happened (sigma): peak {c['diff_peak_dpx']:.1f} px off, shift {c['shift_sigma']:.0f} sigma -> {v}")
        plt.tight_layout(); plt.savefig(os.path.join(REP, f"s{S}_tic{tic}.png"), dpi=100); plt.close(fig)
        # ---- report ----
        md = [f"# TIC {tic}  (sector {S})  verdict: **{v}**", "",
              f"Star: Tmag {t.Tmag:.2f}, Teff {t.Teff:.0f} K, R {t.rad:.3f} Rsun, M {t.mass:.3f} Msun, d {t.d:.1f} pc, contratio {t.contratio:.3f}, camera {int(t.camera)} ccd {int(t.ccd)} col {t.col:.0f} row {t.row:.0f}", "",
              f"Signal: P {row.period_d:.5f} d, t0 {row.t0:.4f}, depth {row.depth_ppm:.0f} ppm, duration {row.duration_h:.2f} h (max {row.max_duration_h:.2f} h), SNR {row.snr:.1f}, {int(row.n_transits)} transits, Rp ~{row.planet_re:.2f} Re",
              f"Vetting: odd {row.odd_ppm:.0f} / even {row.even_ppm:.0f} +/- {row.depth_err_ppm:.0f} ppm; secondary {row.secondary_ppm:.0f} ppm", "",
              f"Catalog: {'already TOI ' + str(toi.toi.tolist()) if len(toi) else 'not a TOI'}; {'known planet host: ' + str(ps.pl_name.tolist()) if len(ps) else 'not a known planet host'}", "",
              f"Centroid: shift {c['shift_px']:.4f} px = {c['shift_sigma']:.1f} sigma (on-target prediction {exp_shift:.4f} px); difference-image peak {c['diff_peak_dpx']:.1f} px from target at {c['diff_peak_snr']:.1f} sigma; n_in {c['n_in']} n_out {c['n_out']}", "",
              "Gaia DR3 within 2 px:", g[g.sep_px <= 2].round(3).to_string(index=False), "",
              "Reasons:"] + [f"- {r}" for r in reasons]
        open(os.path.join(REP, f"s{S}_tic{tic}.md"), "w").write("\n".join(md))
        out_rows.append(dict(ID=tic, verdict=v, is_toi=len(toi) > 0, known_host=len(ps) > 0,
                             shift_px=c["shift_px"], shift_sigma=c["shift_sigma"], exp_shift_px=exp_shift,
                             diff_peak_dpx=c["diff_peak_dpx"], diff_peak_snr=c["diff_peak_snr"],
                             n_neighbors_2px=int((g.sep_px <= 2).sum() - 1), brightest_neighbor_frac=float(g[g.index > 0].light_frac.max()) if len(g) > 1 else 0.0,
                             reasons=" | ".join(reasons)))
        log(f"  -> {v}: " + "; ".join(reasons))
    except Exception as e:
        log(f"  failed: {type(e).__name__}: {e}")
        out_rows.append(dict(ID=tic, verdict=f"error: {type(e).__name__}", reasons=str(e)[:120]))
vpath = os.path.join(DATA, f"vetting_s{S}.csv")
new_df = pd.DataFrame(out_rows)
if os.path.exists(vpath) and len(new_df):
    old_df = pd.read_csv(vpath)
    new_df = pd.concat([old_df[~old_df["ID"].isin(new_df["ID"])], new_df], ignore_index=True)
new_df.to_csv(vpath, index=False)
log(f"wrote {os.path.join(DATA, f'vetting_s{S}.csv')} and reports in {REP}")
