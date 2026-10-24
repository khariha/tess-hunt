"""Fit the transit depth separately to odd- and even-numbered transits, with
all geometry fixed at the joint fit. A planet gives the same depth; an
eclipsing binary at twice the period gives two different depths.

    python -W ignore scripts/oddeven_fit.py --sector 1 --tic 389723320 --period 3.927167
"""
import argparse, os, sys, time
import numpy as np, pandas as pd
from scipy.optimize import least_squares
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tesshunt.transitfit import gather, Model
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); DATA = os.path.join(ROOT, "data")
ap = argparse.ArgumentParser(); ap.add_argument("--sector", type=int, default=1); ap.add_argument("--tic", type=int, required=True); ap.add_argument("--period", type=float, default=None)
args = ap.parse_args(); S0, tic = args.sector, args.tic
def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)
sw = pd.read_csv(os.path.join(DATA, f"sweep_s{S0}.csv")).set_index("ID").loc[tic]
tg = pd.read_csv(os.path.join(DATA, f"targets_s{S0}.csv")).drop_duplicates("ID").set_index("ID").loc[tic]
P0 = args.period or float(sw.period_d); t00 = float(sw.t0); dur0 = float(sw.duration_h) / 24
d = gather(tic, tg.ra, tg.dec, P0, t00, dur0, log=log); t, f, e, exp, sec = d["t"], d["f"], d["e"], d["exp"], d["sector"]
model = Model(t, exp)
n_mid = np.round((np.median(t) - t00) / P0); t0_mid = t00 + n_mid * P0
x0 = [t0_mid, P0, np.sqrt(float(sw.depth_ppm) / 1e6), 15.0, 0.5]
ls = least_squares(lambda th: (f - model(th[0], th[1], th[2], th[3], th[4], 0.3, 0.3)) / e, x0,
                   bounds=([t0_mid - 0.1, P0 * 0.995, 0.005, 1.5, 0], [t0_mid + 0.1, P0 * 1.005, 0.6, 300, 1.2]))
t0, per, rp, ars, b = ls.x
epoch = np.round((t - t0) / per).astype(int)
def fit_depth(mask, label):
    m2 = Model(t[mask], exp[mask])
    r = least_squares(lambda th: (f[mask] - m2(t0, per, th[0], ars, b, 0.3, 0.3)) / e[mask], [rp], bounds=([0.005], [0.6]))
    J = r.jac; cov = np.linalg.inv(J.T @ J) * max(np.sum(r.fun ** 2) / (mask.sum() - 1), 1)
    rp_i, err = r.x[0], np.sqrt(cov[0, 0]); depth, derr = rp_i ** 2 * 1e6, 2 * rp_i * err * 1e6
    log(f"{label:12s} n={mask.sum():5d} transits={len(np.unique(epoch[mask]))}  depth={depth:.0f} +/- {derr:.0f} ppm")
    return depth, derr
log(f"joint: P={per:.6f} rp={rp:.4f} depth={rp**2*1e6:.0f} ppm")
dall = fit_depth(np.ones(len(t), bool), "all")
do = fit_depth(epoch % 2 == 1, "odd"); de = fit_depth(epoch % 2 == 0, "even")
diff = do[0] - de[0]; sig = np.hypot(do[1], de[1])
log(f"odd - even = {diff:+.0f} +/- {sig:.0f} ppm  ({abs(diff)/sig:.1f} sigma)")
for s in np.unique(sec):
    fit_depth(sec == s, f"sector {s}")
# secondary eclipse at phase 0.5: fit a box of the transit duration at t0 + P/2
ph2 = ((t - (t0 + per / 2) + per / 2) % per) - per / 2
in2 = np.abs(ph2) < (dur0 / 2); near2 = (np.abs(ph2) > dur0) & (np.abs(ph2) < 3 * dur0)
if in2.sum() > 5 and near2.sum() > 20:
    sec_depth = (np.average(f[near2], weights=1 / e[near2] ** 2) - np.average(f[in2], weights=1 / e[in2] ** 2)) * 1e6
    sec_err = np.sqrt(1 / np.sum(1 / e[in2] ** 2) + 1 / np.sum(1 / e[near2] ** 2)) * 1e6
    log(f"secondary at phase 0.5: {sec_depth:+.0f} +/- {sec_err:.0f} ppm ({in2.sum()} points)")
else:
    log("secondary: not enough phase-0.5 coverage in the near-transit window (expected; gather() only keeps transit windows)")
