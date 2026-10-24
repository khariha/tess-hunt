"""Sweep a whole sector on Modal, fully in the cloud. Laptop can be closed.

Everything lives on a Modal Volume ("tess-hunt-data"):
    /vol/targets_s<N>.csv        uploaded from data/ on first launch
    /vol/sweep_s<N>.csv          one row per star, resumable
    /vol/candidates/*.png        folded plots for stars that clear the bar

Commands (from the project root):
    modal run --detach scripts/modal_sweep.py --sector 1              # launch, then close the lid
    modal run --detach scripts/modal_sweep.py --sectors "2 3 4 5 6"   # select targets + sweep, sector by sector, all in the cloud
    modal run --detach scripts/modal_sweep.py --sector 1 --tier tiny  # one tier
    modal run scripts/modal_sweep.py --sector 1 --mode status         # how far along
    modal run scripts/modal_sweep.py --sector 1 --mode pull           # copy results to data/ and plots/
    modal run scripts/modal_sweep.py --sector 1 --mode push           # upload local data/sweep_s<N>.csv progress

Costs Modal CPU time only (no GPU): roughly 8 s per star.
"""
import os, sys, time
import modal

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

VOL_NAME = "tess-hunt-data"
vol = modal.Volume.from_name(VOL_NAME, create_if_missing=True)
V = "/vol"

image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("numpy<2.5", "pandas", "scipy", "matplotlib", "astropy", "lightkurve>=2.5",
                 "astrocut>=1.0", "s3fs", "fsspec", "astroquery", "tess-point", "requests")
    .env({"PYTHONWARNINGS": "ignore", "MPLBACKEND": "Agg"})
    .add_local_python_source("tesshunt")
    .add_local_dir(os.path.join(ROOT, "scripts"), remote_path="/root/scripts")
)
app = modal.App("tess-hunt-sweep", image=image)

# ---------- worker: one star ----------
@app.function(cpu=1.0, memory=2048, timeout=900, retries=modal.Retries(max_retries=1),
              max_containers=int(os.environ.get("TESS_MAX_CONTAINERS", "50")))
def sweep_one(t: dict, sector: int, cutout: int, min_snr: float, min_depth_ppm: float):
    from types import SimpleNamespace
    from tesshunt.sweep_core import sweep_star
    return sweep_star(SimpleNamespace(**t), sector, cutout=cutout, source="s3",
                      min_snr=min_snr, min_depth_ppm=min_depth_ppm)

# ---------- driver: hands out stars, collects rows, lives on the Volume ----------
@app.function(cpu=1.0, memory=1024, timeout=24 * 3600, volumes={V: vol})
def drive(sector: int, limit: int, tier: str, max_contratio: float, cutout: int,
          min_snr: float, min_depth_ppm: float):
    import pandas as pd
    from tesshunt.sweep_core import COLS
    vol.reload()
    tpath = f"{V}/targets_s{sector}.csv"
    out_path = f"{V}/sweep_s{sector}.csv"
    cand_dir = f"{V}/candidates"
    os.makedirs(cand_dir, exist_ok=True)

    targets = pd.read_csv(tpath).drop_duplicates("ID")
    targets = targets[targets["contratio"].fillna(0) <= max_contratio]
    targets["tier"] = pd.cut(targets["rad"], [0, 0.30, 0.45, 1.0], labels=["tiny", "small", "warm"]).astype(str)
    if tier != "all":
        targets = targets[targets["tier"] == tier]
    order = {"tiny": 0, "small": 1, "warm": 2}
    targets = targets.assign(_o=targets["tier"].map(order)).sort_values(["_o", "Tmag"])
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
        return {"done": 0, "candidates": 0}

    t0 = time.time(); n = 0; n_cand = 0
    buf = []
    for row, png in sweep_one.map(rows, kwargs=dict(sector=sector, cutout=cutout, min_snr=min_snr,
                                                    min_depth_ppm=min_depth_ppm), order_outputs=False):
        buf.append(row); n += 1
        if png:
            n_cand += 1
            open(f"{cand_dir}/s{sector}_tic{row['ID']}.png", "wb").write(png)
            print(f"  CANDIDATE TIC {row['ID']}  P={row['period_d']:.3f} d  depth={row['depth_ppm']:.0f} ppm  snr={row['snr']:.1f}", flush=True)
        if len(buf) >= 50 or n == len(rows):
            pd.DataFrame(buf, columns=COLS).to_csv(out_path, mode="a", header=not os.path.exists(out_path), index=False)
            buf = []
            vol.commit()
            rate = n / (time.time() - t0)
            print(f"  {n}/{len(rows)} done, {n_cand} candidates, {rate*60:.0f} stars/min, "
                  f"~{(len(rows)-n)/max(rate,1e-9)/60:.0f} min left", flush=True)
    vol.commit()
    print(f"done: {n} stars, {n_cand} candidates, {(time.time()-t0)/60:.1f} min")
    return {"done": n, "candidates": n_cand}

# ---------- cloud target selection (runs scripts/select_targets.py inside a container) ----------
@app.function(cpu=1.0, memory=2048, timeout=4 * 3600, volumes={V: vol})
def select(sector: int):
    import subprocess, shutil
    vol.reload()
    if os.path.exists(f"{V}/targets_s{sector}.csv"):
        print(f"targets_s{sector}.csv already on volume", flush=True); return
    os.makedirs("/root/data", exist_ok=True)
    cache = f"{V}/ctoi_tics.csv"
    if os.path.exists(cache): shutil.copy(cache, "/root/data/ctoi_tics.csv")
    r = subprocess.run([sys.executable, "-W", "ignore", "/root/scripts/select_targets.py", "--sector", str(sector)],
                       capture_output=True, text=True, cwd="/root")
    print(r.stdout[-4000:], flush=True)
    if r.returncode != 0:
        print(r.stderr[-4000:], flush=True); raise RuntimeError(f"select_targets failed for sector {sector}")
    shutil.copy(f"/root/data/targets_s{sector}.csv", f"{V}/targets_s{sector}.csv")
    vol.commit(); print(f"targets_s{sector}.csv written to volume", flush=True)

# ---------- multi-sector pipeline: select then sweep, one sector at a time ----------
@app.function(cpu=1.0, memory=1024, timeout=24 * 3600, volumes={V: vol})
def run_sectors(sectors: list, limit: int, tier: str, max_contratio: float, cutout: int, min_snr: float, min_depth_ppm: float):
    for s in sectors:
        print(f"=== sector {s}: selecting targets ===", flush=True)
        select.remote(s)
        print(f"=== sector {s}: sweeping ===", flush=True)
        out = drive.remote(s, limit, tier, max_contratio, cutout, min_snr, min_depth_ppm)
        print(f"=== sector {s} done: {out} ===", flush=True)
    return {"sectors": sectors}

# ---------- laptop side ----------
@app.local_entrypoint()
def main(sector: int = 1, mode: str = "launch", limit: int = 0, tier: str = "all",
         max_contratio: float = 1.0, cutout: int = 15, min_snr: float = 9.0, min_depth_ppm: float = 3000.0,
         sectors: str = ""):
    import pandas as pd
    if sectors:
        secs = [int(x) for x in sectors.replace(",", " ").split()]
        call = run_sectors.spawn(secs, limit, tier, max_contratio, cutout, min_snr, min_depth_ppm)
        print(f"launched multi-sector pipeline for sectors {secs}, call id {call.object_id}")
        print("each sector: ~25 min target selection + ~3 h sweep. status/pull per sector with --sector N --mode status|pull")
        return
    DATA = os.path.join(ROOT, "data")
    CAND = os.path.join(ROOT, "plots", "candidates")
    os.makedirs(CAND, exist_ok=True)
    local_targets = os.path.join(DATA, f"targets_s{sector}.csv")
    local_sweep = os.path.join(DATA, f"sweep_s{sector}.csv")
    names = {e.path for e in vol.listdir("/")}

    if mode == "push":
        with vol.batch_upload(force=True) as b:
            b.put_file(local_targets, f"targets_s{sector}.csv")
            if os.path.exists(local_sweep):
                b.put_file(local_sweep, f"sweep_s{sector}.csv")
        print("pushed targets and local progress to the volume")
        return

    if mode == "status":
        if f"sweep_s{sector}.csv" not in names:
            print("no sweep file on the volume yet"); return
        df = pd.read_csv(__import__("io").BytesIO(b"".join(vol.read_file(f"sweep_s{sector}.csv"))))
        print(f"{len(df)} stars swept, {int(df.candidate.sum())} candidates, "
              f"{(df.status=='ok').mean()*100:.1f}% ok")
        return

    if mode == "pull":
        if f"sweep_s{sector}.csv" in names:
            open(local_sweep, "wb").write(b"".join(vol.read_file(f"sweep_s{sector}.csv")))
            df = pd.read_csv(local_sweep)
            print(f"pulled {len(df)} rows -> {local_sweep}")
        n = 0
        if "candidates" in names:
            for e in vol.listdir("/candidates"):
                fn = os.path.basename(e.path)
                if fn.startswith(f"s{sector}_"):
                    open(os.path.join(CAND, fn), "wb").write(b"".join(vol.read_file(e.path))); n += 1
        print(f"pulled {n} candidate plots -> {CAND}")
        return

    # mode == "launch": make sure targets (and any local progress) are on the volume, then detach
    if f"targets_s{sector}.csv" not in names:
        with vol.batch_upload() as b:
            b.put_file(local_targets, f"targets_s{sector}.csv")
            if os.path.exists(local_sweep) and f"sweep_s{sector}.csv" not in names:
                b.put_file(local_sweep, f"sweep_s{sector}.csv")
        print("uploaded targets (and local progress) to the volume")
    call = drive.spawn(sector, limit, tier, max_contratio, cutout, min_snr, min_depth_ppm)
    print(f"launched cloud driver, call id {call.object_id}")
    print("close the lid whenever. check with:  modal run scripts/modal_sweep.py --mode status")
    print("fetch results with:                  modal run scripts/modal_sweep.py --mode pull")
