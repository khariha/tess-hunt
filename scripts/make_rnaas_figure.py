"""One figure for the Research Note: three phase-folded transits with the MCMC model.
    python -W ignore scripts/make_rnaas_figure.py
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tesshunt.transitfit import gather, Model
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); DATA = os.path.join(ROOT, "data"); FU = os.path.join(ROOT, "plots", "candidates", "followup")
os.makedirs(os.path.join(ROOT, "docs"), exist_ok=True)
sw = pd.read_csv(os.path.join(DATA, "sweep_s1.csv")).set_index("ID"); tg = pd.read_csv(os.path.join(DATA, "targets_s1.csv")).drop_duplicates("ID").set_index("ID")
cases = [(389723320, None, "TIC 389723320  P = 3.9272 d  4.6 R$_\\oplus$"), (114949540, None, "TIC 114949540  P = 2.7392 d  4.35 R$_\\oplus$")]
fig, axes = plt.subplots(1, 2, figsize=(9, 3.6))
for ax, (tic, secs, title) in zip(axes, cases):
    mc = pd.read_csv(os.path.join(FU, f"s1_tic{tic}_mcmc.csv")).set_index("param")["median"]
    t0, per, rp, ars, b, u1, u2 = [mc[k] for k in ["t0", "period", "rp_rs", "a_rs", "b", "u1", "u2"]]
    dur0 = float(sw.loc[tic].duration_h) / 24
    d = gather(tic, tg.loc[tic].ra, tg.loc[tic].dec, per, t0, dur0, window=2.5, log=lambda m: None, sectors=secs)
    ph = ((d["t"] - t0 + per / 2) % per) - per / 2
    ax.scatter(ph * 24, d["f"], s=2, c="0.7", alpha=0.5, rasterized=True)
    bins = np.linspace(-2.5 * dur0, 2.5 * dur0, 45); idx = np.digitize(ph, bins)
    bm = [np.nanmean(d["f"][idx == i]) if (idx == i).sum() > 2 else np.nan for i in range(1, len(bins))]
    ax.plot((bins[:-1] + bins[1:]) / 2 * 24, bm, "ko", ms=3)
    tt = np.linspace(-2.5 * dur0, 2.5 * dur0, 1200); m = Model(tt + t0, np.full(len(tt), 120 / 86400))
    ax.plot(tt * 24, m(t0, per, rp, ars, b, u1, u2), "r-", lw=1.8)
    lo, hi = np.nanpercentile(d["f"], [1, 99]); ax.set_ylim(lo - 0.2 * (hi - lo), hi + 0.2 * (hi - lo))
    ax.set_title(title, fontsize=10); ax.set_xlabel("hours from mid-transit")
axes[0].set_ylabel("relative flux")
plt.tight_layout(); plt.savefig(os.path.join(ROOT, "docs", "rnaas_figure.png"), dpi=200); plt.savefig(os.path.join(ROOT, "docs", "rnaas_figure.pdf")); print("saved docs/rnaas_figure.png and .pdf")
