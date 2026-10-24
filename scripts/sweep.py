"""Run the FFI transit search down a target list, one star at a time.

Order: tiny stars first (deepest possible dips), then small, then warm;
brightest first within each tier. Every star gets one row in
data/sweep_s<N>.csv whether it succeeded or failed, so the sweep is
resumable: rerun and it skips stars already done.

Stars that clear the candidate bar get a folded plot in plots/candidates/.

Run from the project root, detached so it survives the session:
    nohup python -W ignore scripts/sweep.py --sector 1 > data/sweep_s1.log 2>&1 &
Smoke test:  python -W ignore scripts/sweep.py --sector 1 --limit 3
"""
import argparse, os, sys, time, traceback
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import lightkurve as lk

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tesshunt import pipeline as P

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
CAND = os.path.join(ROOT, "plots", "candidates")
os.makedirs(CAND, exist_ok=True)

ap = argparse.ArgumentParser()
ap.add_argument("--sector", type=int, default=1)
ap.add_argument("--limit", type=int, default=None, help="stop after this many new stars")
ap.add_argument("--max-contratio", type=float, default=1.0, help="skip stars whose neighbors outshine them")
ap.add_argument("--cutout", type=int, default=15)
ap.add_argument("--min-snr", type=float, default=9.0)
ap.add_argument("--min-depth-ppm", type=float, default=3000.0)
args = ap.parse_args()
S = args.sector

def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

COLS = ["ID", "Tmag", "Teff", "rad", "contratio", "camera", "ccd", "status", "n_cadences", "aperture_px",
        "noise_ppm", "period_d", "t0", "duration_h", "max_duration_h", "depth_ppm", "snr",
        "odd_ppm", "even_ppm", "depth_err_ppm", "secondary_ppm", "n_transits", "planet_re", "candidate", "seconds"]

# ---- queue ----
targets = pd.read_csv(os.path.join(DATA, f"targets_s{S}.csv")).drop_duplicates("ID")
targets = targets[targets["contratio"].fillna(0) <= args.max_contratio]
targets["tier"] = pd.cut(targets["rad"], [0, 0.30, 0.45, 1.0], labels=[0, 1, 2]).astype(int)
targets = targets.sort_values(["tier", "Tmag"]).reset_index(drop=True)

out_path = os.path.join(DATA, f"sweep_s{S}.csv")
done = set()
if os.path.exists(out_path):
    done = set(pd.read_csv(out_path, usecols=["ID"])["ID"].astype(str))
    log(f"resuming: {len(done)} stars already swept")
queue = targets[~targets["ID"].astype(str).isin(done)]
if args.limit:
    queue = queue.head(args.limit)
log(f"sector {S}: {len(queue)} stars to sweep ({len(targets)} eligible, {len(done)} done)")

def write_row(row):
    pd.DataFrame([row], columns=COLS).to_csv(out_path, mode="a", header=not os.path.exists(out_path), index=False)

# ---- one star ----
def sweep_star(t):
    row = {c: np.nan for c in COLS}
    row.update(ID=int(t.ID), Tmag=t.Tmag, Teff=t.Teff, rad=t.rad, contratio=t.contratio,
               camera=int(t.camera), ccd=int(t.ccd), candidate=False)
    t_start = time.time()
    try:
        sr = lk.search_tesscut(f"{t.ra} {t.dec}", sector=S)
        if len(sr) == 0:
            row["status"] = "no_cutout"; return row
        tpf = sr.download(cutout_size=args.cutout, download_dir=DATA)
        try:
            lc, raw, aper = P.lightcurve_from_cutout(tpf, threshold=5, n_pca=5)
        except Exception:
            # faint star: fall back to the central 3x3 pixels
            aper = np.zeros(tpf.flux.shape[1:], dtype=bool)
            c = tpf.flux.shape[1] // 2
            aper[c-1:c+2, c-1:c+2] = True
            raw = tpf.to_lightcurve(aperture_mask=aper)
            reg = tpf.flux[:, ~aper]; reg = reg.value if hasattr(reg, "value") else reg
            dm = lk.DesignMatrix(reg, name="bg").pca(5).append_constant()
            lc = lk.RegressionCorrector(raw.remove_nans()).correct(dm)
        row["aperture_px"] = int(aper.sum())
        lc = P.clean(lc)
        flat = P.flatten(lc, window_hours=33)
        flat = P.mask_jitter(flat, window_hours=6, factor=3)
        row["n_cadences"] = len(flat)
        if len(flat) < 300:
            row["status"] = "too_few_cadences"; return row
        row["noise_ppm"] = float(np.nanstd(flat.flux.value)) * 1e6
        r = P.search(flat, min_period=0.3, max_period=13.0, n_periods=15000,
                     min_duration_hours=0.5, max_duration_hours=6.0, duration_step_hours=0.25)
        v = P.vet(flat, r.period, r.t0, r.duration)
        row.update(status="ok", period_d=r.period, t0=r.t0, duration_h=r.duration * 24,
                   max_duration_h=P.max_duration_hours(r.period, (t.mass / t.rad**3) if t.rad > 0 else 1.0),
                   depth_ppm=r.depth * 1e6, snr=r.snr,
                   odd_ppm=v["odd_depth_ppm"], even_ppm=v["even_depth_ppm"], depth_err_ppm=v["depth_err_ppm"],
                   secondary_ppm=v["secondary_depth_ppm"], n_transits=v["n_transits"],
                   planet_re=P.planet_radius_earths(max(r.depth, 0), t.rad))
        is_cand = (r.snr >= args.min_snr and r.depth * 1e6 >= args.min_depth_ppm and v["n_transits"] >= 2)
        row["candidate"] = bool(is_cand)
        if is_cand:
            folded = P.fold(flat, r.period, r.t0)
            fig, ax = plt.subplots(2, 1, figsize=(9, 7))
            flat.plot(ax=ax[0]); ax[0].set_title(f"TIC {int(t.ID)}  Tmag {t.Tmag:.1f}  R={t.rad:.2f} Rsun")
            folded.scatter(ax=ax[1], s=4, alpha=0.5)
            folded.bin(time_bin_size=max(r.duration / 8, 0.005)).plot(ax=ax[1], color="red", lw=2)
            ax[1].set_xlim(-2 * r.duration, 2 * r.duration)
            ax[1].set_title(f"P={r.period:.3f} d  depth={r.depth*1e6:.0f} ppm  dur={r.duration*24:.1f} h  SNR={r.snr:.1f}  "
                            f"odd/even={v['odd_depth_ppm']:.0f}/{v['even_depth_ppm']:.0f}  sec={v['secondary_depth_ppm']:.0f}")
            plt.tight_layout(); plt.savefig(os.path.join(CAND, f"s{S}_tic{int(t.ID)}.png"), dpi=110); plt.close(fig)
        return row
    except Exception as e:
        row["status"] = f"error: {type(e).__name__}: {str(e)[:80]}"
        return row
    finally:
        row["seconds"] = round(time.time() - t_start, 1)

# ---- loop ----
n_cand = 0
for i, t in enumerate(queue.itertuples(index=False)):
    row = sweep_star(t)
    write_row(row)
    if row.get("candidate"):
        n_cand += 1
    tag = "CANDIDATE" if row.get("candidate") else row["status"]
    log(f"{i+1}/{len(queue)} TIC {int(t.ID)} Tmag {t.Tmag:.1f} R {t.rad:.2f} | {tag} | "
        f"P={row['period_d']:.3f} depth={row['depth_ppm']:.0f} snr={row['snr']:.1f} | {row['seconds']}s"
        if row["status"] == "ok" else f"{i+1}/{len(queue)} TIC {int(t.ID)} | {tag} | {row['seconds']}s")
log(f"done. {n_cand} candidates this run. results in {out_path}")
