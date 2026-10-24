"""Run the FFI transit search down a target list on this laptop, one star at a time.

Order: tiny stars first, then small, then warm; brightest first within a tier.
Resumable: one row per star in data/sweep_s<N>.csv; rerun skips stars done.
Candidates get a plot in plots/candidates/.

    nohup python -W ignore scripts/sweep.py --sector 1 > data/sweep_s1.log 2>&1 &
    python -W ignore scripts/sweep.py --sector 1 --limit 3          # smoke test

For the whole sector use scripts/modal_sweep.py (parallel, in the cloud).
"""
import argparse, os, sys, time
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tesshunt.sweep_core import sweep_star, COLS

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
CAND = os.path.join(ROOT, "plots", "candidates")
os.makedirs(CAND, exist_ok=True)

ap = argparse.ArgumentParser()
ap.add_argument("--sector", type=int, default=1)
ap.add_argument("--limit", type=int, default=None)
ap.add_argument("--max-contratio", type=float, default=1.0)
ap.add_argument("--cutout", type=int, default=15)
ap.add_argument("--source", choices=["s3", "tesscut"], default="s3")
ap.add_argument("--min-snr", type=float, default=9.0)
ap.add_argument("--min-depth-ppm", type=float, default=3000.0)
args = ap.parse_args()
S = args.sector

def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

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
log(f"sector {S}: {len(queue)} stars to sweep ({len(targets)} eligible, {len(done)} done), source={args.source}")

n_cand = 0
for i, t in enumerate(queue.itertuples(index=False)):
    row, png = sweep_star(t, S, cutout=args.cutout, source=args.source,
                          min_snr=args.min_snr, min_depth_ppm=args.min_depth_ppm)
    pd.DataFrame([row], columns=COLS).to_csv(out_path, mode="a", header=not os.path.exists(out_path), index=False)
    if png:
        n_cand += 1
        open(os.path.join(CAND, f"s{S}_tic{int(t.ID)}.png"), "wb").write(png)
    tag = "CANDIDATE" if row["candidate"] else row["status"]
    if row["status"] == "ok":
        log(f"{i+1}/{len(queue)} TIC {int(t.ID)} Tmag {t.Tmag:.1f} R {t.rad:.2f} | {tag} | "
            f"P={row['period_d']:.3f} depth={row['depth_ppm']:.0f} snr={row['snr']:.1f} | {row['seconds']}s")
    else:
        log(f"{i+1}/{len(queue)} TIC {int(t.ID)} | {tag} | {row['seconds']}s")
log(f"done. {n_cand} candidates this run. results in {out_path}")
