"""One star, start to finish: pixels -> light curve -> search -> vet -> result row.

Shared by scripts/sweep.py (laptop, sequential) and scripts/modal_sweep.py
(cloud, parallel). Keep all science here so both paths behave identically.
"""
import io, time
import numpy as np
import lightkurve as lk
from tesshunt import pipeline as P

COLS = ["ID", "Tmag", "Teff", "rad", "contratio", "camera", "ccd", "status", "n_cadences", "aperture_px",
        "noise_ppm", "period_d", "t0", "duration_h", "max_duration_h", "depth_ppm", "snr",
        "odd_ppm", "even_ppm", "depth_err_ppm", "secondary_ppm", "n_transits", "planet_re", "candidate", "seconds"]

def get_tpf(ra, dec, sector, size=15, source="s3"):
    """Fetch the pixel cube for one star.

    source="s3": slice MAST's public sector cubes on AWS directly (no rate limit,
                 fast when run near us-east-1). Needs astrocut >= 1.0.
    source="tesscut": NASA's hosted cutout service (5 req/s, 10 concurrent).
    """
    if source == "s3":
        from astrocut import TessFootprintCutout
        c = TessFootprintCutout(f"{ra} {dec}", cutout_size=size, sequence=sector)
        if not c.tpf_cutouts:
            return None
        return lk.TessTargetPixelFile(c.tpf_cutouts[0])
    sr = lk.search_tesscut(f"{ra} {dec}", sector=sector)
    if len(sr) == 0:
        return None
    return sr.download(cutout_size=size)

def _central_aperture(tpf, half=1):
    aper = np.zeros(tpf.flux.shape[1:], dtype=bool)
    c = tpf.flux.shape[1] // 2
    aper[c-half:c+half+1, c-half:c+half+1] = True
    raw = tpf.to_lightcurve(aperture_mask=aper)
    reg = tpf.flux[:, ~aper]; reg = reg.value if hasattr(reg, "value") else reg
    dm = lk.DesignMatrix(reg, name="bg").pca(5).append_constant()
    return lk.RegressionCorrector(raw.remove_nans()).correct(dm), raw, aper

def sweep_star(t, sector, cutout=15, source="s3", min_snr=9.0, min_depth_ppm=3000.0, make_plot=True):
    """t: anything with .ID .ra .dec .Tmag .Teff .rad .mass .contratio .camera .ccd
    Returns (row dict with COLS keys, png bytes or None)."""
    row = {c: np.nan for c in COLS}
    row.update(ID=int(t.ID), Tmag=t.Tmag, Teff=t.Teff, rad=t.rad, contratio=t.contratio,
               camera=int(t.camera), ccd=int(t.ccd), candidate=False)
    png = None
    t_start = time.time()
    try:
        tpf = get_tpf(t.ra, t.dec, sector, cutout, source)
        if tpf is None:
            row["status"] = "no_cutout"; return row, None
        try:
            lc, raw, aper = P.lightcurve_from_cutout(tpf, threshold=5, n_pca=5)
            if aper.sum() == 0 or aper.sum() > 0.5 * aper.size:
                raise ValueError("bad aperture")
        except Exception:
            lc, raw, aper = _central_aperture(tpf)
        row["aperture_px"] = int(aper.sum())
        lc = P.clean(lc)
        flat = P.flatten(lc, window_hours=33)
        flat = P.mask_jitter(flat, window_hours=6, factor=3)
        row["n_cadences"] = len(flat)
        if len(flat) < 300:
            row["status"] = "too_few_cadences"; return row, None
        row["noise_ppm"] = float(np.nanstd(flat.flux.value)) * 1e6
        r = P.search(flat, min_period=0.3, max_period=13.0, n_periods=15000,
                     min_duration_hours=0.5, max_duration_hours=6.0, duration_step_hours=0.25)
        v = P.vet(flat, r.period, r.t0, r.duration)
        density = (t.mass / t.rad ** 3) if (t.rad and t.rad > 0 and t.mass and t.mass > 0) else 1.0
        row.update(status="ok", period_d=r.period, t0=r.t0, duration_h=r.duration * 24,
                   max_duration_h=P.max_duration_hours(r.period, density),
                   depth_ppm=r.depth * 1e6, snr=r.snr,
                   odd_ppm=v["odd_depth_ppm"], even_ppm=v["even_depth_ppm"], depth_err_ppm=v["depth_err_ppm"],
                   secondary_ppm=v["secondary_depth_ppm"], n_transits=v["n_transits"],
                   planet_re=P.planet_radius_earths(max(r.depth, 0), t.rad))
        is_cand = bool(r.snr >= min_snr and r.depth * 1e6 >= min_depth_ppm and v["n_transits"] >= 2)
        row["candidate"] = is_cand
        if is_cand and make_plot:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            folded = P.fold(flat, r.period, r.t0)
            fig, ax = plt.subplots(2, 1, figsize=(9, 7))
            flat.plot(ax=ax[0]); ax[0].set_title(f"TIC {int(t.ID)}  Tmag {t.Tmag:.1f}  R={t.rad:.2f} Rsun  sector {sector}")
            folded.scatter(ax=ax[1], s=4, alpha=0.5)
            folded.bin(time_bin_size=max(r.duration / 8, 0.005)).plot(ax=ax[1], color="red", lw=2)
            ax[1].set_xlim(-2 * r.duration, 2 * r.duration)
            ax[1].set_title(f"P={r.period:.3f} d  depth={r.depth*1e6:.0f} ppm  dur={r.duration*24:.1f} h  SNR={r.snr:.1f}  "
                            f"odd/even={v['odd_depth_ppm']:.0f}/{v['even_depth_ppm']:.0f}  sec={v['secondary_depth_ppm']:.0f}")
            plt.tight_layout()
            buf = io.BytesIO(); plt.savefig(buf, format="png", dpi=110); plt.close(fig)
            png = buf.getvalue()
        return row, png
    except Exception as e:
        row["status"] = f"error: {type(e).__name__}: {str(e)[:80]}"
        return row, None
    finally:
        row["seconds"] = round(time.time() - t_start, 1)
