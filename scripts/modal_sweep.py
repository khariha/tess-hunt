"""Sweep a whole sector in parallel on Modal CPU containers, slicing MAST's
public S3 sector cubes directly (no TESScut rate limit).

    modal run scripts/modal_sweep.py --sector 1 --limit 200        # trial
    modal run scripts/modal_sweep.py --sector 1                    # whole list
    modal run scripts/modal_sweep.py --sector 1 --tier tiny        # one tier

Results stream back to data/sweep_s<N>.csv (same file and columns as the
laptop sweep, resumable) and candidate plots to plots/candidates/.
Costs money: Modal CPU time only, no GPU. Kishore runs this.
"""
import os, sys, time
import modal

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("numpy<2.5", "pandas", "scipy", "matplotlib", "astropy", "lightkurve>=2.5",
                 "astrocut>=1.0", "s3fs", "fsspec")
    .env({"PYTHONWARNINGS": "ignore", "MPLBACKEND": "Agg"})
    .add_local_python_source("tesshunt")
)
app = modal.App("tess-hunt-sweep", image=image)

@app.function(cpu=1.0, memory=2048, timeout=900, retries=modal.Retries(max_retries=1),
              max_containers=int(os.environ.get("TESS_MAX_CONTAINERS", "50")))
def sweep_one(t: dict, sector: int, cutout: int, min_snr: float, min_depth_ppm: float):
    from types import SimpleNamespace
    from tesshunt.sweep_core import sweep_star
    row, png = sweep_star(SimpleNamespace(**t), sector, cutout=cutout, source="s3",
                          min_snr=min_snr, min_depth_ppm=min_depth_ppm)
    return row, png

@app.local_entrypoint()
def main(sector: int = 1, limit: int = 0, tier: str = "all", max_contratio: float = 1.0,
         cutout: int = 15, min_snr: float = 9.0, min_depth_ppm: float = 3000.0):
    import pandas as pd
    from tesshunt.sweep_core import COLS
    DATA = os.path.join(ROOT, "data")
    CAND = os.path.join(ROOT, "plots", "candidates")
    os.makedirs(CAND, exist_ok=True)

    targets = pd.read_csv(os.path.join(DATA, f"targets_s{sector}.csv")).drop_duplicates("ID")
    targets = targets[targets["contratio"].fillna(0) <= max_contratio]
    targets["tier"] = pd.cut(targets["rad"], [0, 0.30, 0.45, 1.0], labels=["tiny", "small", "warm"]).astype(str)
    if tier != "all":
        targets = targets[targets["tier"] == tier]
    targets = targets.sort_values(["tier", "Tmag"], key=lambda s: s.map({"tiny": 0, "small": 1, "warm": 2}) if s.name == "tier" else s)

    out_path = os.path.join(DATA, f"sweep_s{sector}.csv")
    done = set()
    if os.path.exists(out_path):
        done = set(pd.read_csv(out_path, usecols=["ID"])["ID"].astype(str))
    queue = targets[~targets["ID"].astype(str).isin(done)]
    if limit:
        queue = queue.head(limit)
    rows = [{k: (None if pd.isna(v) else v) for k, v in r.items()}
            for r in queue[["ID", "ra", "dec", "Tmag", "Teff", "rad", "mass", "contratio", "camera", "ccd"]].to_dict("records")]
    print(f"sector {sector}: {len(rows)} stars to sweep ({len(done)} already done)", flush=True)
    if not rows:
        return

    t0 = time.time(); n = 0; n_cand = 0
    for row, png in sweep_one.map(rows, kwargs=dict(sector=sector, cutout=cutout, min_snr=min_snr,
                                                    min_depth_ppm=min_depth_ppm), order_outputs=False):
        pd.DataFrame([row], columns=COLS).to_csv(out_path, mode="a", header=not os.path.exists(out_path), index=False)
        n += 1
        if png:
            n_cand += 1
            open(os.path.join(CAND, f"s{sector}_tic{row['ID']}.png"), "wb").write(png)
            print(f"  CANDIDATE TIC {row['ID']}  P={row['period_d']:.3f} d  depth={row['depth_ppm']:.0f} ppm  snr={row['snr']:.1f}", flush=True)
        if n % 50 == 0:
            rate = n / (time.time() - t0)
            print(f"  {n}/{len(rows)} done, {n_cand} candidates, {rate*60:.0f} stars/min, "
                  f"~{(len(rows)-n)/max(rate,1e-9)/60:.0f} min left", flush=True)
    print(f"done: {n} stars, {n_cand} candidates, {(time.time()-t0)/60:.1f} min. results in {out_path}")
