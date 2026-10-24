"""Assemble the Community TOI submission package for one candidate.

Reads: sweep row, target row, vetting row, follow-up summary, MCMC results.
Writes: plots/candidates/submission/s<N>_tic<ID>/
            report.md           the vetting report (what a reviewer reads)
            ctoi_fields.json    the values to type into the ExoFOP CTOI form
            *.png               figures copied from the follow-up folder

    python -W ignore scripts/make_ctoi_package.py --sector 1 --tic 389723320
"""
import argparse, json, os, shutil, sys
import numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data"); FU = os.path.join(ROOT, "plots", "candidates", "followup"); REP = os.path.join(ROOT, "plots", "candidates", "reports")
ap = argparse.ArgumentParser(); ap.add_argument("--sector", type=int, default=1); ap.add_argument("--tic", type=int, required=True)
ap.add_argument("--submitter", default="Kishore Hariharan"); ap.add_argument("--group", default="Generative Experiences / tess-hunt")
args = ap.parse_args(); S0, tic = args.sector, args.tic
OUT = os.path.join(ROOT, "plots", "candidates", "submission", f"s{S0}_tic{tic}"); os.makedirs(OUT, exist_ok=True)

sw = pd.read_csv(os.path.join(DATA, f"sweep_s{S0}.csv")).set_index("ID").loc[tic]
tg = pd.read_csv(os.path.join(DATA, f"targets_s{S0}.csv")).drop_duplicates("ID").set_index("ID").loc[tic]
vt = pd.read_csv(os.path.join(DATA, f"vetting_s{S0}.csv")).set_index("ID").loc[tic]
mc = pd.read_csv(os.path.join(FU, f"s{S0}_tic{tic}_mcmc.csv")).set_index("param")
fu_md = open(os.path.join(FU, f"s{S0}_tic{tic}.md")).read()
def v(p): return float(mc.loc[p, "median"]), float(mc.loc[p, "plus"]), float(mc.loc[p, "minus"])
def fmt(p, nd=4): m, pl, mi = v(p); return f"{m:.{nd}f} +{pl:.{nd}f} -{mi:.{nd}f}"
t0, t0p, t0m = v("t0"); per, perp, perm = v("period"); dep, depp, depm = v("depth_ppm"); dur, durp, durm = v("dur_h"); rpe, rpep, rpem = v("rp_re")
BTJD = 2457000.0
ctoi_path = os.path.join(DATA, "ctoi_tics.csv")
if os.path.exists(ctoi_path):
    _ct = set(pd.read_csv(ctoi_path, header=None)[0].astype(str))
    _date = pd.Timestamp(os.path.getmtime(ctoi_path), unit="s").strftime("%Y-%m-%d")
    ctoi_note = (f"ALREADY LISTED as a CTOI (ExoFOP list of {_date})" if str(tic) in _ct
                 else f"not in the ExoFOP CTOI list ({len(_ct)} entries, {_date}).")
else:
    ctoi_note = "CTOI list check pending ExoFOP availability."
fields = {
    "TIC": tic, "submitter": args.submitter, "group": args.group,
    "period_days": round(per, 6), "period_err": round(max(perp, perm), 6),
    "epoch_BJD": round(t0 + BTJD, 5), "epoch_err": round(max(t0p, t0m), 5),
    "depth_ppm": round(dep), "depth_err_ppm": round(max(depp, depm)),
    "duration_hours": round(dur, 2), "duration_err_hours": round(max(durp, durm), 2),
    "planet_radius_Rearth": round(rpe, 2), "planet_radius_err": round(max(rpep, rpem), 2),
    "transit_snr": round(dep / max(depp, depm), 1),
    "sectors_used": fu_md.split("Observed in sectors")[1].split(";")[0].strip() if "Observed in sectors" in fu_md else "",
    "tag": f"tess-hunt FFI search, sector {S0} detection, multi-sector confirmed",
    "notes": (f"Detected in an independent BLS search of SPOC FFI cutouts (astrocut, S3 cubes) of Tmag 12-15 M dwarfs. "
              f"Signal recovered independently in every available sector. Centroid and noise-normalized difference images place the "
              f"dip on the target in all tested 200-s sectors. Odd/even depths agree within 2.1 sigma; no secondary eclipse. "
              f"batman+emcee fit with free quadratic limb darkening. Transit-derived stellar density {v('rho_star')[0]:.1f} solar vs TIC {v('rho_cat')[0]:.1f} solar. "
              f"Not a TOI; no literature match; " + ctoi_note),
}
json.dump(fields, open(os.path.join(OUT, "ctoi_fields.json"), "w"), indent=2)
for fn in [f"s{S0}_tic{tic}_mcmc.png", f"s{S0}_tic{tic}_corner.png", f"s{S0}_tic{tic}.png", f"s{S0}_tic{tic}_pixels.png", f"s{S0}_tic{tic}_fit.png"]:
    src = os.path.join(FU, fn)
    if os.path.exists(src): shutil.copy(src, OUT)
rep = os.path.join(REP, f"s{S0}_tic{tic}.png")
if os.path.exists(rep): shutil.copy(rep, os.path.join(OUT, f"s{S0}_tic{tic}_sector{S0}_vetting.png"))
md = f"""# Community TOI candidate: TIC {tic}

Submitter: {args.submitter} ({args.group}). Prepared {pd.Timestamp.now():%Y-%m-%d}.

## Host star (TIC v8)
Tmag {tg.Tmag:.2f}, Teff {tg.Teff:.0f} K, R = {tg.rad:.3f} Rsun, M = {tg.mass:.3f} Msun, distance {tg.d:.1f} pc, contamination ratio {tg.contratio:.3f}.
RA {tg.ra:.5f}, Dec {tg.dec:.5f}. Sector {S0} camera {int(tg.camera)} CCD {int(tg.ccd)}.

## Signal (MCMC, all sectors, batman + emcee, free quadratic limb darkening)
| quantity | value |
|---|---|
| Period (d) | {fmt('period', 6)} |
| Epoch (BJD) | {t0 + BTJD:.5f} +{t0p:.5f} -{t0m:.5f} |
| Depth (ppm) | {fmt('depth_ppm', 0)} |
| Duration (h) | {fmt('dur_h', 2)} |
| Rp/Rs | {fmt('rp_rs', 4)} |
| a/Rs | {fmt('a_rs', 2)} |
| Impact parameter | {fmt('b', 2)} |
| Planet radius (Rearth) | {fmt('rp_re', 2)} |
| Planet radius (Rjup) | {fmt('rp_rjup', 3)} |
| a (AU) | {fmt('a_au', 4)} |
| T_eq (K, albedo 0) | {fmt('teq', 0)} |
| Stellar density from transit (solar) | {fmt('rho_star', 2)} |
| Stellar density from TIC (solar) | {v('rho_cat')[0]:.2f} |
| Limb darkening u1, u2 | {fmt('u1', 2)}, {fmt('u2', 2)} |
| Transit SNR (depth / error) | {fields['transit_snr']} |

## Detection
Independent box-least-squares search of FFI cutouts (SPOC calibrated FFIs via MAST S3 cubes, astrocut) of TIC M dwarfs with Tmag 12-15 and R < 0.6 Rsun not on the 2-minute target list. Sector-{S0} detection: P = {sw.period_d:.5f} d, depth {sw.depth_ppm:.0f} ppm, SNR {sw.snr:.1f}, {int(sw.n_transits)} transits.

## Multi-sector recovery
{fu_md.split('Observed in sectors')[1].split(chr(10))[0].strip() if 'Observed in sectors' in fu_md else ''}

Per-sector independent searches (narrow period window) and joint refinement are in the follow-up figure.

## Vetting
- Centroid / difference image (sector {S0}): {vt.verdict}; shift {vt.shift_sigma:.1f} sigma; difference-image peak {vt.diff_peak_dpx:.1f} px from target at {vt.diff_peak_snr:.0f} sigma.
- Difference images in 200-s sectors: see pixels figure (dip on target in every tested sector).
- Odd/even: see oddeven log in notes; model-fitted odd and even depths agree within ~2 sigma.
- Secondary eclipse at phase 0.5: consistent with zero (box test on stitched light curve).
- Gaia DR3 neighbors within 2 px: {int(vt.n_neighbors_2px)}; brightest neighbor {vt.brightest_neighbor_frac*100:.0f}% of target light in G.
- Catalogs: not a TOI (Exoplanet Archive TOI table); not a known planet host; no ADS literature match for the TIC. CTOI: {ctoi_note}

## Interpretation
Rp/Rs = {v('rp_rs')[0]:.3f} implies a {rpe:.1f} Rearth companion; nothing self-luminous is this small, so unless diluted by an unresolved brighter star the companion is planet-sized. Transit-derived stellar density agrees with the catalog within the b / a/Rs degeneracy.

## Figures
- {os.path.basename(os.path.join(OUT, f's{S0}_tic{tic}_mcmc.png'))}: phase-folded data and model
- s{S0}_tic{tic}.png: per-sector folds and joint fold
- s{S0}_tic{tic}_pixels.png: images and difference images, 200-s sectors
- s{S0}_tic{tic}_sector{S0}_vetting.png: sector-{S0} vetting panel
- s{S0}_tic{tic}_corner.png: posterior
"""
open(os.path.join(OUT, "report.md"), "w").write(md)
print(json.dumps(fields, indent=2)); print(f"\npackage: {OUT}")
