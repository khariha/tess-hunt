"""Fetch and cache the FFI pixel cube of one star for every sector that observed it.

    nohup python -W ignore scripts/fetch_cubes.py --tic 259719966 > data/fetch_259719966.log 2>&1 &

Writes data/cubes/s<N>_tic<ID>_15x15.fits (full astrocut TPF with WCS) once per sector; skips existing files.
Other scripts load them with tesshunt.cubes.load_tpf(tic, sector), which falls back to the S3 fetch when missing.
"""
import argparse, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pandas as pd
from tess_stars2px import tess_stars2px_function_entry
from tesshunt.cubes import cube_path, load_tpf, ROOT
ap = argparse.ArgumentParser(); ap.add_argument("--tic", type=int, required=True); ap.add_argument("--sector", type=int, default=1)
ap.add_argument("--cutout", type=int, default=15); ap.add_argument("--sectors", type=int, nargs="*", default=None)
a = ap.parse_args()
def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)
tg = pd.read_csv(os.path.join(ROOT, "data", f"targets_s{a.sector}.csv")).drop_duplicates("ID").set_index("ID").loc[a.tic]
sectors = a.sectors or sorted(set(int(s) for s in tess_stars2px_function_entry(a.tic, tg.ra, tg.dec)[3]))
log(f"TIC {a.tic}: sectors {sectors}")
for s in sectors:
    p = cube_path(a.tic, s, a.cutout)
    if os.path.exists(p): log(f"sector {s}: cached {p}"); continue
    t0 = time.time()
    tpf = load_tpf(a.tic, s, float(tg.ra), float(tg.dec), a.cutout)
    log(f"sector {s}: {'no cube' if tpf is None else f'{tpf.flux.shape} saved'} in {time.time()-t0:.0f} s")
