# TIC 389723320: all-phase secondary scan and per-sector pixel tests (2026-10-25)

Run after TIC 259719966 turned out to be an eclipsing binary with its secondary at phase 0.54. `scripts/deep_vet.py --tic 389723320`
(ephemeris from the MCMC fit: P = 3.927167 d, T0 = 3899.38949, T14 = 1.50 h, Rp/Rs = 0.1137, a/Rs = 16.64, b = 0.74); tables in
`data/deepvet_389723320/`, log `data/deepvet_389723320.log`; figures `addendum_*.png` here. Six sectors (1, 28, 68, 95, 105, 106), 37 transits.

**Verdict: still a planet candidate.** No secondary eclipse at any phase; dip on target in every sector; two soft spots noted below.

| test | result |
|---|---|
| secondary at phase 0.50 | +62 +/- 191 ppm |
| secondary at any phase (box of T14, phases 0.06-0.95) | scan rms 289 ppm; largest feature +1,364 ppm at phase 0.32, which per sector is 279, 377, -295, 160, 277, 1,216 (+/- 370-580) ppm: only the noisiest sector (106, 9 ppt scatter) contributes, so it is not a repeating eclipse (`secondary_peak_per_sector.csv`). 3-sigma limit on a secondary at any phase: ~900 ppm |
| odd/even (37 transits) | odd 13,858 +/- 444, even 12,376 +/- 362, difference +1,482 +/- 573 ppm (2.6 sigma). Soft spot: worth re-checking with more sectors; at 2x the period this would be a 1.4 vs 1.2% pair of eclipses, which the flat scan at phase 0.5 of the doubled period (= phase 0.25/0.75 here, both ~0) does not support |
| difference image / centroid, per sector | peak 0.12-0.89 px from target, 6-16 sigma, 2.2-3.0 px from the brightest neighbour; centroid shift 0.5-3.0 sigma; difference-image centroid 3-9 arcsec from target in five sectors, 17 arcsec in the weak sector 106 (6 sigma) |
| blend budget (Gaia DR3 within 2.5 px) | G 17.39 at 31.7 arcsec (1.5 px, 23% of target flux) would need a 7.3% eclipse; its own 1-px aperture shows 25-91 ppm (857 ppm in S105 where apertures touch). G 20.10 at 1.8 px would need 74%. No source inside 1 px |
| per-sector depth (pipeline aperture, 1-3 px) | 14.3, 12.5, 11.9, 11.6, 16.9, 13.4 kppm (+/- 0.4-0.9): a 40% spread, largest in the 1-px aperture of S105. Soft spot: aperture-dependent dilution/PSF losses; the submitted radius should quote the multi-sector MCMC with this spread as a systematic |
| transit times (34 fitted) | O-C rms 2.9 min, chi2/N = 0.90, no TTV |
| activity | no coherent rotation period; 0-2 flares per sector, none in transit |

Figures: `addendum_pixels_all_sectors.png`, `addendum_oddeven_secondary.png`, `addendum_secondary_scan_peak.png`, `addendum_depth_timing_activity.png`.
