# TIC 114949540: all-phase secondary scan and per-sector pixel tests (2026-10-25)

Run after TIC 259719966 turned out to be an eclipsing binary with its secondary at phase 0.54. `scripts/deep_vet.py --tic 114949540`
(ephemeris from the MCMC fit: P = 2.739176 d, T0 = 4212.950524, T14 = 1.55 h, Rp/Rs = 0.0930, a/Rs = 13.12, b = 0.53); tables in
`data/deepvet_114949540/`, log `data/deepvet_114949540.log`; figures `addendum_*.png` here. Three sectors (1, 27, 105), 26 transits.

**Verdict: still a planet candidate.** No secondary at any phase, odd = even, dip on target in all three sectors. One caveat TESS cannot
remove: a Gaia source 13.9 arcsec (0.66 px) away with 12% of the target's flux sits inside every aperture and would need an 8.3% eclipse
to fake the dip; only ground-based photometry that resolves it can close that.

| test | result |
|---|---|
| secondary at phase 0.50 | +175 +/- 190 ppm |
| secondary at any phase | scan rms 248 ppm; largest feature +422 ppm at phase 0.15, per sector 65, -176, -286 (+/- 350-420) ppm: nothing repeats. 3-sigma limit ~750 ppm |
| odd/even (26 transits) | -276 +/- 448 ppm (0.6 sigma) |
| difference image / centroid, per sector | peak 0.26-0.64 px from target at 4, 14, 8 sigma; 2.8-3.5 px from the brightest neighbour; centroid shift 0.5-2.2 sigma; difference-image centroid 5-7 arcsec from target |
| blend budget (Gaia DR3 within 2.5 px) | G 18.34 at 13.9 arcsec / 0.66 px (12% of target): needs 8.3%, unresolvable by TESS. G 16.99 at 45 arcsec / 2.16 px (41%): needs 3.0%; outside every aperture (2.95-2.99 px away) and its own aperture shows 25 ppm in S1/S27, 1,664 ppm in S105 (target wing). Four G ~20 sources at 1.8-2.1 px need 34-83% |
| per-sector depth | 7.9, 8.7, 9.7 kppm (+/- 0.3-0.5) in 2, 1, 2-px apertures |
| transit times (25 fitted) | O-C rms 3.5 min, chi2/N = 0.90, no TTV |
| activity | no coherent rotation; flares rare, none in transit |

Figures: `addendum_pixels_all_sectors.png`, `addendum_oddeven_secondary.png`, `addendum_secondary_scan_peak.png`, `addendum_depth_timing_activity.png`.
