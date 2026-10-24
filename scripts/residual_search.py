"""Search for a SECOND planet: mask the known transits, stitch every sector, run BLS down to a low depth floor.

    nohup python -W ignore scripts/residual_search.py --sector 1 --tic 389723320 --period 3.927167 > data/residual_389723320.log 2>&1 &
Writes plots/candidates/followup/s<N>_tic<ID>_residual.png and prints the best residual signal with vetting numbers.
"""
import argparse, os, sys, time
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import lightkurve as lk
from tess_stars2px import tess_stars2px_function_entry
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tesshunt import pipeline as P
from tesshunt.sweep_core import get_tpf, _central_aperture
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); DATA = os.path.join(ROOT, "data"); OUT = os.path.join(ROOT, "plots", "candidates", "followup")
ap = argparse.ArgumentParser(); ap.add_argument("--sector", type=int, default=1); ap.add_argument("--tic", type=int, required=True)
ap.add_argument("--period", type=float, required=True); ap.add_argument("--t0", type=float, default=None); ap.add_argument("--duration-h", type=float, default=None)
ap.add_argument("--min-period", type=float, default=0.5); ap.add_argument("--max-period", type=float, default=30.0)
args = ap.parse_args(); S0, tic = args.sector, args.tic
def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)
sw = pd.read_csv(os.path.join(DATA, f"sweep_s{S0}.csv")).set_index("ID").loc[tic]
tg = pd.read_csv(os.path.join(DATA, f"targets_s{S0}.csv")).drop_duplicates("ID").set_index("ID").loc[tic]
P1 = args.period; t01 = args.t0 or float(sw.t0); dur1 = (args.duration_h or float(sw.duration_h)) / 24
sectors = sorted(set(int(s) for s in tess_stars2px_function_entry(tic, tg.ra, tg.dec)[3]))
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
        # mask the known planet BEFORE flattening so the smoother does not bend around its dips
        ph = ((lc.time.value - t01 + P1 / 2) % P1) - P1 / 2
        keep = np.abs(ph) > 1.0 * dur1
        flat = P.mask_jitter(P.flatten(lc[keep], 33), 6, 3)
        # bin 200-s sectors to 10 min so all sectors weigh in comparably and BLS stays fast
        cad = np.nanmedian(np.diff(flat.time.value)) * 1440
        if cad < 8: flat = flat.bin(time_bin_size=10 / 1440).remove_nans()
        flats.append(flat); log(f"sector {s}: {len(flat)} points after masking planet 1 ({cad:.0f}-min cadence)")
    except Exception as e:
        log(f"sector {s}: {type(e).__name__}: {e}")
st = lk.LightCurveCollection(flats).stitch(); st = st[np.isfinite(st.flux.value)]
# drop isolated single-point dips: a real transit spans several 10-min points, a glitch spans one
fv = st.flux.value; sig = np.nanstd(fv)
low = fv < 1 - 3 * sig
isolated = low & ~np.roll(low, 1) & ~np.roll(low, -1)
st = st[~isolated]; log(f"dropped {isolated.sum()} isolated single-point dips")
base = st.time.max().value - st.time.min().value
log(f"stitched: {len(st)} points over {base:.0f} days, scatter {np.nanstd(st.flux.value)*1e6:.0f} ppm per 10-min point")
# period grid fine enough for the long baseline
q = 0.5 / 24  # allow half an hour of phase drift across the whole baseline
periods = []
p = args.min_period
while p < args.max_period:
    periods.append(p); p += p * p * q / base
periods = np.array(periods); log(f"{len(periods)} trial periods (astropy BLS, direct)")
from astropy.timeseries import BoxLeastSquares
tt_, ff_ = st.time.value, st.flux.value
dy_ = np.full(len(ff_), np.nanstd(ff_))
bls = BoxLeastSquares(tt_, ff_, dy=dy_)
durs = np.array([0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0]) / 24
res = bls.power(periods, durs, objective="snr")
power = np.asarray(res.power); med = np.nanmedian(power); mad = np.nanmedian(np.abs(power - med)) * 1.4826
class _PG:  # minimal stand-in for the lightkurve periodogram used below
    def __init__(self, per, pw): self.period = per * __import__("astropy.units", fromlist=["day"]).day; self.power = pw
    def plot(self, ax=None):
        ax.plot(self.period.value, self.power, "k-", lw=0.5); ax.set_xlabel("period [d]"); ax.set_ylabel("BLS power"); return ax
pg = _PG(periods, power)
order = np.argsort(power)[::-1]
log("top residual peaks (period, depth ppm, SNR):")
seen = []
rows = []
for i in order:
    per = float(pg.period[i].value)
    if any(abs(per / s - round(per / s)) < 0.01 or abs(s / per - round(s / per)) < 0.01 for s in seen): continue
    if any(abs(per / (P1 * k) - 1) < 0.01 for k in (0.5, 1, 2, 3)): continue
    seen.append(per)
    rows.append(dict(period=per, snr=(power[i] - med) / mad))
    if len(rows) >= 12: break
best = None
for k, rrow in enumerate(rows):
    per = rrow["period"]
    r2 = bls.power(np.linspace(per * 0.999, per * 1.001, 400), durs, objective="snr")
    j = int(np.nanargmax(r2.power)); t0 = float(r2.transit_time[j]); dur = float(r2.duration[j]); dep = float(r2.depth[j])
    v = P.vet(st, per, t0, dur)
    if v["n_transits"] < 3 or dur * 24 < 1.0:
        continue
    log(f"  {per:9.4f} d  depth {dep*1e6:6.0f} ppm  SNR {rrow['snr']:5.1f}  transits {v['n_transits']:3d}  odd/even {v['odd_depth_ppm']:.0f}/{v['even_depth_ppm']:.0f} +/- {v['depth_err_ppm']:.0f}  sec {v['secondary_depth_ppm']:.0f}  Rp~{P.planet_radius_earths(max(dep,0), float(tg.rad)):.2f} Re  maxdur {P.max_duration_hours(per, float(tg.mass)/float(tg.rad)**3):.1f} h vs {dur*24:.1f} h")
    if best is None: best = (per, t0, dur, dep, rrow["snr"], v)
fig, ax = plt.subplots(2, 1, figsize=(11, 8))
pg.plot(ax=ax[0]); ax[0].set_title(f"TIC {tic}: residual periodogram after masking P={P1:.4f} d")
if best is None:
    log("RESULT: no residual signal with >= 3 transits and >= 1 h duration among the top peaks")
    ax[1].text(0.5, 0.5, "no residual signal passes (>=3 transits, >=1 h)", ha="center", va="center", transform=ax[1].transAxes)
    plt.tight_layout(); plt.savefig(os.path.join(OUT, f"s{S0}_tic{tic}_residual.png"), dpi=100); log("done"); sys.exit(0)
per, t0, dur, dep, snr, v = best
f = st.fold(period=per, epoch_time=t0); f.scatter(ax=ax[1], s=2, alpha=0.3); f.bin(time_bin_size=max(dur / 6, 0.003)).plot(ax=ax[1], color="red", lw=2)
ax[1].set_xlim(-3 * dur, 3 * dur); lo, hi = np.nanpercentile(f.flux.value, [1, 99]); ax[1].set_ylim(lo, hi)
ax[1].set_title(f"best residual: P={per:.4f} d depth={dep*1e6:.0f} ppm SNR={snr:.1f} transits={v['n_transits']}")
plt.tight_layout(); plt.savefig(os.path.join(OUT, f"s{S0}_tic{tic}_residual.png"), dpi=100); log("done")
