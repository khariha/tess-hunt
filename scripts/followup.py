"""Multi-sector follow-up for one candidate: does the signal repeat in every sector?

For a TIC with a Sector-N detection (period P, epoch t0, depth):
    1. list every sector that observed the star (tess-point)
    2. for each sector, build the FFI light curve from the S3 cube and search a
       narrow period window around P -> does this sector see the same dip?
    3. stitch all sectors, refine P on the full baseline, fold, and vet again
    4. write plots/candidates/followup/s<N>_tic<ID>.png and a .md summary

    python -W ignore scripts/followup.py --sector 1 --tic 389723320
"""
import argparse, os, sys, time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import lightkurve as lk
from tess_stars2px import tess_stars2px_function_entry

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tesshunt import pipeline as P
from tesshunt.sweep_core import get_tpf, _central_aperture

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "plots", "candidates", "followup")
os.makedirs(OUT, exist_ok=True)

ap = argparse.ArgumentParser()
ap.add_argument("--sector", type=int, default=1, help="sector of the original detection")
ap.add_argument("--tic", type=int, required=True)
ap.add_argument("--cutout", type=int, default=15)
args = ap.parse_args()
S0, tic = args.sector, args.tic

def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)

sw = pd.read_csv(os.path.join(DATA, f"sweep_s{S0}.csv")).set_index("ID").loc[tic]
tg = pd.read_csv(os.path.join(DATA, f"targets_s{S0}.csv")).drop_duplicates("ID").set_index("ID").loc[tic]
P0, t00, dur0, depth0 = float(sw.period_d), float(sw.t0), float(sw.duration_h) / 24, float(sw.depth_ppm)
sectors = sorted(set(int(s) for s in tess_stars2px_function_entry(tic, tg.ra, tg.dec)[3]))
log(f"TIC {tic}: sector-{S0} detection P={P0:.5f} d depth={depth0:.0f} ppm; observed in sectors {sectors}")

per_sector = []
flats = []
for s in sectors:
    try:
        tpf = get_tpf(tg.ra, tg.dec, s, args.cutout, "s3")
        if tpf is None:
            per_sector.append(dict(sector=s, status="no cube")); log(f"  sector {s}: no cube"); continue
        cad_min = float(np.nanmedian(np.diff(tpf.time.value)) * 1440)
        try:
            lc, raw, aper = P.lightcurve_from_cutout(tpf, threshold=5)
            if aper.sum() == 0 or aper.sum() > 0.5 * aper.size: raise ValueError
        except Exception:
            lc, raw, aper = _central_aperture(tpf)
        flat = P.mask_jitter(P.flatten(P.clean(lc), 33), 6, 3)
        if len(flat) < 300:
            per_sector.append(dict(sector=s, status="too few cadences", cadence_min=cad_min)); continue
        r = P.search(flat, min_period=P0 * 0.97, max_period=P0 * 1.03, n_periods=3000,
                     min_duration_hours=0.5, max_duration_hours=max(4.0, dur0 * 24 * 2), duration_step_hours=0.25)
        v = P.vet(flat, r.period, r.t0, r.duration)
        noise = float(np.nanstd(flat.flux.value)) * 1e6
        per_sector.append(dict(sector=s, status="ok", cadence_min=cad_min, n=len(flat), noise_ppm=noise,
                               period=r.period, t0=r.t0, duration_h=r.duration * 24, depth_ppm=r.depth * 1e6,
                               snr=r.snr, n_transits=v["n_transits"], odd=v["odd_depth_ppm"], even=v["even_depth_ppm"],
                               err=v["depth_err_ppm"], secondary=v["secondary_depth_ppm"]))
        flat.meta["SECTOR"] = s
        flats.append(flat)
        log(f"  sector {s} ({cad_min:.0f}-min): P={r.period:.4f} depth={r.depth*1e6:.0f} ppm snr={r.snr:.1f} transits={v['n_transits']}")
    except Exception as e:
        per_sector.append(dict(sector=s, status=f"error: {type(e).__name__}: {str(e)[:60]}")); log(f"  sector {s}: {type(e).__name__}: {e}")

ps = pd.DataFrame(per_sector)
ok = ps[ps.status == "ok"] if "status" in ps else ps.iloc[0:0]
detected = ok[(ok.snr >= 7) & (np.abs(ok.period - P0) / P0 < 0.02) & (ok.depth_ppm > 0.4 * depth0)]
log(f"signal recovered in {len(detected)} of {len(ok)} usable sectors")

# ---- joint refine ----
joint = None
if len(flats) >= 2:
    stitched = lk.LightCurveCollection(flats).stitch()
    stitched = stitched[np.isfinite(stitched.flux.value)]
    base = stitched.time.max().value - stitched.time.min().value
    periods = np.linspace(P0 * 0.995, P0 * 1.005, 40000)
    pg = stitched.to_periodogram(method="bls", period=periods, duration=np.arange(0.5, max(4.0, dur0 * 48) + 1e-9, 0.25) / 24)
    Pj, t0j, durj = float(pg.period_at_max_power.value), float(pg.transit_time_at_max_power.value), float(pg.duration_at_max_power.value)
    dj = pg.depth_at_max_power; dj = float(dj.value) if hasattr(dj, "value") else float(dj)
    vj = P.vet(stitched, Pj, t0j, durj)
    joint = dict(period=Pj, t0=t0j, duration_h=durj * 24, depth_ppm=dj * 1e6, baseline_d=base, n_transits=vj["n_transits"],
                 odd=vj["odd_depth_ppm"], even=vj["even_depth_ppm"], err=vj["depth_err_ppm"], secondary=vj["secondary_depth_ppm"],
                 period_err_est=P0 * durj / base)
    joint["secondary_any_phase_ppm"] = vj.get("secondary_any_phase_ppm", float("nan")); joint["secondary_any_phase_at"] = vj.get("secondary_any_phase_at", float("nan"))
    log(f"joint: P={Pj:.6f} d (+/-~{joint['period_err_est']:.6f}) depth={dj*1e6:.0f} ppm over {base:.0f} days, {vj['n_transits']} transits, odd/even {vj['odd_depth_ppm']:.0f}/{vj['even_depth_ppm']:.0f} +/- {vj['depth_err_ppm']:.0f}, secondary at 0.5 {vj['secondary_depth_ppm']:.0f}, deepest at any phase {joint['secondary_any_phase_ppm']:.0f} ppm (phase {joint['secondary_any_phase_at']:.2f})")
    if np.isfinite(joint["secondary_any_phase_ppm"]) and joint["secondary_any_phase_ppm"] > 4 * vj["depth_err_ppm"]:
        log("WARNING: secondary eclipse away from phase 0.5 -> likely eccentric eclipsing binary")

# ---- figure ----
n = len(flats)
fig, axes = plt.subplots(max(n, 1) + (1 if joint else 0), 1, figsize=(10, 3.2 * (max(n, 1) + (1 if joint else 0))))
axes = np.atleast_1d(axes)
for ax, flat, row in zip(axes, flats, [r for r in per_sector if r["status"] == "ok"]):
    Pp, tp = (joint["period"], joint["t0"]) if joint else (row["period"], row["t0"])
    f = flat.fold(period=Pp, epoch_time=tp)
    f.scatter(ax=ax, s=3, alpha=0.4)
    f.bin(time_bin_size=max(row["duration_h"] / 24 / 6, 0.003)).plot(ax=ax, color="red", lw=1.5)
    ax.set_xlim(-2.5 * dur0, 2.5 * dur0)
    lo, hi = np.nanpercentile(f.flux.value, [0.5, 99.5]); ax.set_ylim(lo - 0.3 * (hi - lo), hi + 0.3 * (hi - lo))
    ax.set_title(f"sector {row['sector']} ({row['cadence_min']:.0f}-min cadence): own search P={row['period']:.4f} depth={row['depth_ppm']:.0f} ppm snr={row['snr']:.1f}", fontsize=10)
if joint:
    ax = axes[-1]
    f = stitched.fold(period=joint["period"], epoch_time=joint["t0"])
    f.scatter(ax=ax, s=3, alpha=0.3)
    f.bin(time_bin_size=max(joint["duration_h"] / 24 / 8, 0.003)).plot(ax=ax, color="red", lw=2)
    ax.set_xlim(-2.5 * dur0, 2.5 * dur0)
    lo, hi = np.nanpercentile(f.flux.value, [0.5, 99.5]); ax.set_ylim(lo - 0.3 * (hi - lo), hi + 0.3 * (hi - lo))
    ax.set_title(f"ALL SECTORS stitched: P={joint['period']:.6f} d depth={joint['depth_ppm']:.0f} ppm, {joint['n_transits']} transits over {joint['baseline_d']:.0f} days", fontsize=10)
plt.tight_layout(); plt.savefig(os.path.join(OUT, f"s{S0}_tic{tic}.png"), dpi=100)

md = [f"# TIC {tic} follow-up", "", f"Sector-{S0} detection: P {P0:.5f} d, depth {depth0:.0f} ppm", "",
      f"Observed in sectors {sectors}; usable {len(ok)}; signal recovered in {len(detected)}: {detected.sector.tolist() if len(detected) else []}", "",
      ps.round(3).to_string(index=False), ""]
if joint:
    md += ["Joint refine:", pd.Series(joint).round(6).to_string()]
open(os.path.join(OUT, f"s{S0}_tic{tic}.md"), "w").write("\n".join(md))
log(f"wrote {OUT}/s{S0}_tic{tic}.png and .md")
