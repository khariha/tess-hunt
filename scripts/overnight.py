"""Unattended vetting chain for a list of sectors. Resumable; one log; one summary at the end.

    nohup python -W ignore scripts/overnight.py --sectors 2 3 4 5 6 > data/overnight.log 2>&1 &

Per sector: pull sweep results + target list from the Modal volume, vet every flagged star
(vet_candidates.py, now with the any-phase secondary test). For every star that vets on_target
or likely_on_target: multi-sector follow-up; if the signal repeats in every usable sector and
no off-phase secondary appears, odd/even fit, MCMC fit, and a submission package.
Writes data/overnight_summary.md.
"""
import argparse, os, re, subprocess, sys, time
import pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); DATA = os.path.join(ROOT, "data"); SCR = os.path.join(ROOT, "scripts")
ap = argparse.ArgumentParser(); ap.add_argument("--sectors", type=int, nargs="+", required=True); ap.add_argument("--skip-pull", action="store_true")
args = ap.parse_args()
PY = [sys.executable, "-W", "ignore"]
def log(m): print(f"[{time.strftime('%m-%d %H:%M:%S')}] {m}", flush=True)
def run(cmd, timeout=3600):
    r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, timeout=timeout)
    return r.returncode, (r.stdout or "") + (r.stderr or "")

summary = ["# Overnight vetting summary", f"Started {time.strftime('%Y-%m-%d %H:%M')}", ""]
grand = []
for S in args.sectors:
    log(f"=== sector {S} ===")
    if not args.skip_pull:
        rc, out = run(["modal", "run", f"{SCR}/modal_sweep.py", "--sector", str(S), "--mode", "pull"], timeout=1800)
        log("pull: " + " | ".join(l for l in out.splitlines() if "pulled" in l.lower()))
        if not os.path.exists(os.path.join(DATA, f"targets_s{S}.csv")):
            rc, out = run(["modal", "volume", "get", "tess-hunt-data", f"targets_s{S}.csv", os.path.join(DATA, f"targets_s{S}.csv"), "--force"], timeout=600)
            log(f"targets_s{S}.csv fetched (rc {rc})")
    sw = pd.read_csv(os.path.join(DATA, f"sweep_s{S}.csv"))
    n_sw, n_ok, n_c = len(sw), int((sw.status == "ok").sum()), int(sw.candidate.sum())
    log(f"sweep: {n_sw} stars, {n_ok/n_sw*100:.1f}% ok, {n_c} flagged")
    rc, out = run(PY + [f"{SCR}/vet_candidates.py", "--sector", str(S)], timeout=7200)
    log("vet: " + str(len(re.findall(r"-> ", out))) + " verdicts" + ("" if rc == 0 else f" (rc {rc})"))
    vpath = os.path.join(DATA, f"vetting_s{S}.csv")
    v = pd.read_csv(vpath) if os.path.exists(vpath) else pd.DataFrame(columns=["ID", "verdict"])
    counts = v.verdict.value_counts().to_dict(); log(f"verdicts: {counts}")
    keep = v[v.verdict.isin(["on_target", "likely_on_target"])]
    rows = []
    for tic in keep.ID.astype(int):
        swr = sw[sw.ID == tic].iloc[0]
        log(f"  TIC {tic}: follow-up ({swr.period_d:.4f} d, {swr.depth_ppm:.0f} ppm)")
        rc, out = run(PY + [f"{SCR}/followup.py", "--sector", str(S), "--tic", str(tic)], timeout=3600)
        m = re.search(r"signal recovered in (\d+) of (\d+) usable sectors", out)
        mj = re.search(r"joint: P=([\d.]+) d .*?depth=(\d+) ppm.*?(\d+) transits", out)
        anyph = re.search(r"deepest at any phase (-?\d+) ppm \(phase ([\d.]+)\)", out)
        eb = "WARNING: secondary eclipse away from phase 0.5" in out
        rec, usable = (int(m.group(1)), int(m.group(2))) if m else (0, 0)
        row = dict(sector=S, ID=tic, verdict=keep[keep.ID == tic].verdict.iloc[0], period=swr.period_d, depth_ppm=swr.depth_ppm,
                   snr=swr.snr, rad=swr.rad, Tmag=swr.Tmag, recovered=f"{rec}/{usable}",
                   joint_period=float(mj.group(1)) if mj else float("nan"), joint_depth_ppm=int(mj.group(2)) if mj else -1,
                   joint_transits=int(mj.group(3)) if mj else 0, any_phase_secondary_ppm=int(anyph.group(1)) if anyph else -1,
                   eclipsing_binary=eb, stage="followup")
        if usable >= 2 and rec == usable and not eb:
            Pj = row["joint_period"]
            log(f"  TIC {tic}: repeats in {rec}/{usable}; fitting")
            run(PY + [f"{SCR}/oddeven_fit.py", "--sector", str(S), "--tic", str(tic), "--period", f"{Pj:.6f}"], timeout=3600)
            rc, out = run(PY + [f"{SCR}/mcmc_fit.py", "--sector", str(S), "--tic", str(tic), "--period", f"{Pj:.6f}", "--walkers", "40", "--steps", "6000", "--burn", "2000"], timeout=7200)
            mr = re.search(r"^\s*rp_re\s+([\d.]+)", out, re.M); mb = re.search(r"^\s*b\s+([\d.]+)", out, re.M)
            row["rp_re"] = float(mr.group(1)) if mr else float("nan"); row["b"] = float(mb.group(1)) if mb else float("nan")
            run(PY + [f"{SCR}/make_ctoi_package.py", "--sector", str(S), "--tic", str(tic)], timeout=600)
            row["stage"] = "packaged"
        elif usable >= 2 and rec == usable and eb:
            row["stage"] = "eclipsing_binary"
        else:
            row["stage"] = "failed_recovery"
        log(f"  TIC {tic}: {row['stage']} (recovered {row['recovered']}, any-phase secondary {row['any_phase_secondary_ppm']} ppm)")
        rows.append(row)
    df = pd.DataFrame(rows); grand.append(df)
    summary += [f"## Sector {S}", f"- swept {n_sw} stars ({n_ok/n_sw*100:.1f}% ok), {n_c} flagged, verdicts: {counts}",
                f"- carried to follow-up: {len(keep)}; packaged: {int((df.stage == 'packaged').sum()) if len(df) else 0}; eclipsing binaries: {int((df.stage == 'eclipsing_binary').sum()) if len(df) else 0}; failed recovery: {int((df.stage == 'failed_recovery').sum()) if len(df) else 0}", ""]
    if len(df):
        summary += [df[["ID", "verdict", "period", "depth_ppm", "snr", "rad", "Tmag", "recovered", "joint_transits", "any_phase_secondary_ppm", "stage"] + ([c for c in ("rp_re", "b") if c in df.columns])].round(3).to_markdown(index=False), ""]
    open(os.path.join(DATA, "overnight_summary.md"), "w").write("\n".join(summary))
if grand:
    pd.concat(grand).to_csv(os.path.join(DATA, "overnight_results.csv"), index=False)
summary += [f"Finished {time.strftime('%Y-%m-%d %H:%M')}"]
open(os.path.join(DATA, "overnight_summary.md"), "w").write("\n".join(summary)); log("done")
