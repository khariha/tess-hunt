"""Search for a secondary eclipse at ANY phase, not just 0.5 (eccentric orbits move it).

    python -W ignore scripts/secondary_scan.py --sector 1 --tic 389723320 --period 3.927167 --t0 3899.38949 --dur-h 1.50
Stitches every sector (transits masked), slides a box of the transit duration across phase 0.05..0.95,
and reports the deepest dip with its significance and the phase-0.5 value.
"""
import argparse, os, sys, time
import numpy as np, pandas as pd
import lightkurve as lk
from tess_stars2px import tess_stars2px_function_entry
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tesshunt import pipeline as P
from tesshunt.sweep_core import get_tpf, _central_aperture
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); DATA = os.path.join(ROOT, "data")
ap = argparse.ArgumentParser(); ap.add_argument("--sector", type=int, default=1); ap.add_argument("--tic", type=int, required=True)
ap.add_argument("--period", type=float, required=True); ap.add_argument("--t0", type=float, required=True); ap.add_argument("--dur-h", type=float, required=True)
a = ap.parse_args()
def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)
tg = pd.read_csv(os.path.join(DATA, f"targets_s{a.sector}.csv")).drop_duplicates("ID").set_index("ID").loc[a.tic]
sectors = sorted(set(int(s) for s in tess_stars2px_function_entry(a.tic, tg.ra, tg.dec)[3]))
P0, t0, dur = a.period, a.t0, a.dur_h / 24
flats = []
for s in sectors:
    try:
        tpf = get_tpf(tg.ra, tg.dec, s, 15, "s3")
        if tpf is None: continue
        try:
            lc, raw, aper = P.lightcurve_from_cutout(tpf, threshold=5)
            if aper.sum() == 0 or aper.sum() > 0.5 * aper.size: raise ValueError
        except Exception:
            lc, raw, aper = _central_aperture(tpf)
        lc = P.clean(lc)
        ph = ((lc.time.value - t0 + P0 / 2) % P0) - P0 / 2
        flat = P.mask_jitter(P.flatten(lc[np.abs(ph) > 1.0 * dur], 33), 6, 3)
        flat.meta["SECTOR"] = s; flats.append(flat); log(f"sector {s}: {len(flat)} points")
    except Exception as e:
        log(f"sector {s}: {type(e).__name__}: {e}")
# robust version: medians, MAD errors, and a real eclipse must appear in EVERY sector at the same phase
secs_used = []
for fl in flats:
    fl.meta["SEC"] = int(fl.meta.get("SECTOR", 0)) if fl.meta.get("SECTOR") else len(secs_used)
    secs_used.append(fl)
def depth_at(fl, c):
    t, f = fl.time.value, fl.flux.value
    phase = ((t - t0) % P0) / P0
    w = np.abs(((phase - c + 0.5) % 1) - 0.5) < (dur / P0) / 2
    if w.sum() < 4: return np.nan, np.nan, 0
    mad = np.nanmedian(np.abs(f - np.nanmedian(f))) * 1.4826
    return (np.nanmedian(f[~w]) - np.nanmedian(f[w])) * 1e6, mad / np.sqrt(w.sum()) * 1e6 * 1.25, int(w.sum())
rows = []
for c in np.arange(0.10, 0.901, 0.005):
    ds = [depth_at(fl, c) for fl in flats]
    d = np.array([x[0] for x in ds]); e = np.array([x[1] for x in ds]); ok = np.isfinite(d) & np.isfinite(e) & (e > 0)
    if ok.sum() == 0: continue
    wmean = np.sum(d[ok] / e[ok] ** 2) / np.sum(1 / e[ok] ** 2); werr = np.sqrt(1 / np.sum(1 / e[ok] ** 2))
    n_pos = int(np.sum(d[ok] > 2 * e[ok]))     # sectors where the dip is individually significant
    rows.append((c, wmean, werr, n_pos, int(ok.sum()), ";".join(f"{x:+.0f}" for x in d[ok])))
df = pd.DataFrame(rows, columns=["phase", "depth_ppm", "err_ppm", "n_sectors_2sig", "n_sectors", "per_sector_ppm"]); df["sigma"] = df.depth_ppm / df.err_ppm
best = df.loc[df.sigma.idxmax()]; p5 = df.iloc[(df.phase - 0.5).abs().argmin()]
log(f"TIC {a.tic}: {len(flats)} sectors")
log(f"  phase 0.50: {p5.depth_ppm:+.0f} +/- {p5.err_ppm:.0f} ppm ({p5.sigma:.1f} sigma), per sector {p5.per_sector_ppm}")
log(f"  deepest anywhere: {best.depth_ppm:+.0f} +/- {best.err_ppm:.0f} ppm ({best.sigma:.1f} sigma) at phase {best.phase:.3f}; individually >2 sigma in {best.n_sectors_2sig}/{best.n_sectors} sectors; per sector {best.per_sector_ppm}")
top = df.sort_values("sigma", ascending=False).head(4)
for _, r in top.iterrows(): log(f"    phase {r.phase:.3f}: {r.depth_ppm:+.0f} +/- {r.err_ppm:.0f} ppm ({r.sigma:.1f} sigma), {r.n_sectors_2sig}/{r.n_sectors} sectors agree")
df.to_csv(os.path.join(DATA, f"secondary_scan_{a.tic}.csv"), index=False)
