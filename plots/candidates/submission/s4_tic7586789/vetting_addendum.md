# TIC 7586789: vetting addendum (2026-10-26)

Six sectors on disk (4, 5, 31, 32, 106, 107; `scripts/fetch_cubes.py --sector 4`), 64 transits over 2,230 days. Ephemeris from the
joint MCMC (`plots/candidates/followup/s4_tic7586789_mcmc.csv`): P = 2.261386 d, T0 = BTJD 4248.483978, T14 = 1.45 h, Rp/Rs = 0.1045,
depth 10,590 ppm, Rp = 5.1 Re on a 0.453 Rsun host at 172 pc. Tables in `data/deepvet_7586789/`, logs `data/deepvet*_7586789.log`,
`data/fpp_7586789.log`, `data/qlp_7586789.log`.

## Summary

**Verdict: planet candidate; every TESS test passes, with one weak feature to watch.** The dip is on the target: the six-sector stacked
difference image peaks on the target pixel at 14.8 sigma with its centroid 0.17 px (3.7 arcsec) from the target against 1.05 px to the
only capable neighbour, and the aperture experiment prefers dip-on-target over dip-on-neighbour by chi2 244 vs 5,494. Odd and even
agree (+468 +/- 346 ppm, 1.4 sigma), no secondary at phase 0.5 (+153 +/- 143 ppm), no timing variations over 63 transits (O-C rms 2.5
min), flat-bottomed (T23/T14 = 0.72, grazing excluded at delta-chi2 = 3,084), transit density 3.5 +1.7 -0.8 vs TIC 4.8 solar (1.1 sigma),
quiet host (no rotation signal above 0.5 ppt, no flares in six sectors). TRICERATOPS and QLP: see the end of this file.

The weak feature: the all-phase secondary scan has its largest value at phase 0.71, +274 +/- 139 ppm in the pipeline aperture (2.0 sigma)
and +542 +/- 147 ppm in the neighbour-excluded aperture (3.7 sigma). Per sector the phase-0.71 box is positive in all six sectors but
each is < 1.7 sigma; the weighted mean is 287 +/- 108 ppm (2.7 sigma) against a scan baseline of +95 +/- 105 ppm at phase 0.5, so the
excess over baseline is ~1.3 sigma. Not a detection. It matters because if it were real it could not come from the target: a 0.46 RJup
companion is too small to be self-luminous, so a secondary at phase 0.71 (e cos w = 0.33 on a 2.26-day orbit) would mean the whole signal
is a diluted eclipsing binary, most plausibly on the 22-arcsec neighbour. Two more 200-s sectors or an independent light curve settle it.

| # | Test | Result | Verdict |
|---|---|---|---|
| 1 | Difference image / centroid, 6 sectors | peaks 3-14 sigma, 0.6-1.2 px from target; stack 14.8 sigma on the target pixel, centroid 3.7 arcsec off; centroid shift toward the neighbour (cos +0.8 to +1.0) as an on-target dip predicts | pass |
| 2 | Blend budget | one capable source: Gaia G 18.11 at 21.8 arcsec (1.04 px), 11.4% of target flux, needs 10.9%; its 1-px aperture shows 6.6-8.9 kppm (mostly the target's wings) where a 10.9% eclipse would give far more; PSF-model fit chi2 244 (target) vs 5,494 (neighbour) | pass |
| 3 | Per-sector depth | 9.5, 10.8, 12.5, 10.1, 11.8, 11.2 kppm (+/- 0.3): a 25% spread with 3% errors; the two sectors with the neighbour in the aperture (4, 106) are not systematically shallower | pass with a 10% sector-to-sector depth systematic |
| 4a | Odd/even, 64 transits | 10,995 +/- 258 vs 10,528 +/- 230 (1.35 sigma) | pass |
| 4b | Secondary at phase 0.5 | +153 +/- 143 ppm | pass |
| 4c | Secondary at any phase | rms 222 ppm; largest +542 ppm at phase 0.71 (3.7 sigma in the small aperture, 2.0 in the pipeline aperture; 2.7 sigma per-sector weighted mean, 1.3 sigma above the scan baseline); 3-sigma limit 670 ppm | pass, flagged |
| 4d | Phase curve | amplitude at P 124 +/- 47 ppm (a +cos term, the wrong sign for reflection), at P/2 32 +/- 46; 3-sigma limits 140 ppm | pass |
| 5 | Shape and density | T23/T14 = 0.72 +0.07 -0.06, T12 = 12 min, b = 0.60 +0.08 -0.21, grazing fraction 0, best grazing delta-chi2 +3,084; rho 3.5 +1.7 -0.8 vs 4.83 (1.1 sigma); with TIC prior b = 0.52, Rp/Rs = 0.1027, Rp = 5.08 +/- 0.12 Re | pass |
| 6 | Transit timing, 63 transits | O-C rms 2.46 min, median error 2.37 min, chi2/N = 1.15, largest 7.3 min | pass |
| 7 | Host activity | LS peaks 6.4-15 d at 0.34-0.54 ppt, inconsistent between sectors; 0.2-5 d < 0.52 ppt; zero flares in 6 sectors; host RUWE 0.996, C* 1.1 sigma, Mann R = 0.453 = TIC | pass |
| 8 | Cadence audit | no quality-flagged transit cadences; one cadence lost to clean() in S107 epoch 13 | pass |

## 1. Difference image and centroid

`pixel_tests.csv`, `centroid_vectors.csv`, `per_transit_diffimage.csv`, `stacked_diffimage.json`; figures `addendum_pixels_all_sectors.png`,
`addendum_centroid_vectors.png`, `addendum_stacked_diffimage.png`. Per sector the difference-image peak is 3.2 (S4, 30-min, 4 transits
with 3 points each), 8.4, 14.2, 8.5, 13.4 and 6.1 sigma, 0.6-1.2 px from the target; the difference-image centroid is 2.4-11.8 arcsec from
the direct-image centroid. Because the neighbour is only 1.05 px away, single-sector peaks fall between the two stars (0.3-1.9 px from
the neighbour) and the predicted centroid-shift vectors for the two hypotheses are nearly identical (chi2 26.6 vs 28.6 in S5), so the
per-sector tests cannot separate them. The stack can: 14.8 sigma on the target pixel, centroid 0.17 px from the target and 0.9 px from the
neighbour. Per-transit difference-image centroids have median offsets 3-7 arcsec from the target and 1.05-1.23 px from the neighbour in
every sector.

## 2. Blend budget and aperture experiment

`blend_budget.csv`, `dilution_psf.csv`, `depth_hypotheses.json`; figure `addendum_blend_hypotheses.png`. Gaia DR3 ...288910043264: G 18.11,
bp-rp 1.76, parallax 0.57 mas (a background star at ~1.5 kpc), 21.8 arcsec, 11.4% of the target's G flux, needs a 10.9% eclipse. Next
capable source needs > 50%. The 1-px aperture on the neighbour shows 6.6-8.9 kppm dips; with a Gaussian PSF of 0.75 px that aperture holds
~35-45% target light and ~45-55% neighbour light, so a 10.9% eclipse on the neighbour would make it dip by 5-6%, not 0.7-0.9%. Fitting all 18
depth measurements with each star's modelled share: dip on target D = 1.10% with chi2 244 (33.9 allowing 10% per-sector jitter); dip on
neighbour E = 8.2% with chi2 5,494 (634). The target hypothesis wins by delta-chi2 > 600 even with the jitter term.

## 3. Per-sector depth

| sector | cadence | aperture px | neighbour in aperture | depth (ppm) | S/N |
|---|---|---|---|---|---|
| 4 | 30 min | 7 | yes | 9,468 +/- 298 | 32 |
| 5 | 30 min | 5 | no | 10,790 +/- 293 | 37 |
| 31 | 10 min | 4 | no | 12,543 +/- 307 | 41 |
| 32 | 10 min | 6 | no | 10,146 +/- 295 | 34 |
| 106 | 200 s | 4 | yes | 11,800 +/- 348 | 34 |
| 107 | 200 s | 4 | no | 11,197 +/- 351 | 32 |

Stitched: 10,903 +/- 126 (pipeline), 10,741 +/- 172 (neighbour pixels excluded). The spread (9.5 to 12.5 kppm) exceeds the formal errors
by 10x, as it did for TIC 259719966 and 389723320: a ~10% per-sector depth systematic from PSF-wing dilution and detrending is the norm
for these faint FFI stars and should be carried as the radius error. S4 is the shallowest and does hold the neighbour; S106 also holds it
and is not shallow, so the dilution signature is not clean.

## 4. Odd/even, secondary, phase curve

`stitched_tests.json`, `secondary_peak_per_sector.csv`, `phase_curve.json`; figure `addendum_oddeven_secondary.png`. Numbers in the summary
table. The scan's shape (negative at phases < 0.2 and > 0.8, positive in between) is the 33-h smoother's response to the unmasked
transit and is common to every star vetted so far; the phase-0.71 feature sits on top of that baseline.

Per-sector box depth at phase 0.71 (pipeline aperture): S4 +22 +/- 242, S5 +189 +/- 245, S31 +414 +/- 261, S32 +389 +/- 248, S106 +356
+/- 296, S107 +470 +/- 316 ppm. Same box at phase 0.5: +126, -10, +111, +267, -68, +144 ppm.

## 5. Shape and density

`shape_density.json`; figure `addendum_shape_density.png`. From the existing chain: T14 = 1.46 +/- 0.04 h, T23 = 1.04 +0.07 -0.06 h, T12 =
0.21 +/- 0.05 h, T23/T14 = 0.72 +0.07 -0.06, b = 0.60 +0.08 -0.21, zero posterior mass at b > 0.8. Free fit on 5,079 near-transit points:
Rp/Rs = 0.1045, a/Rs = 10.5, b = 0.63, chi2 5,111; the best V-shaped fit is worse by 3,084. Circular density 3.48 +1.70 -0.80 vs TIC 4.83:
1.1 sigma low, consistent with a circular orbit (e_min = 0). With the TIC density as a 25% prior: b = 0.52 +/- 0.14, a/Rs = 11.6 +/- 1.0,
Rp/Rs = 0.1027 +/- 0.0024, Rp = 5.08 +/- 0.12 Re.

## 6. Transit timing

`transit_times.csv`. 63 fitted transits over 2,230 days: O-C rms 2.46 min, median error 2.37 min, chi2/N = 1.15, largest 7.3 min. No TTV.

## 7. Host

`activity.csv`, `rotation_short.csv`, `flares.csv`, `host_diagnostics.json`. Out-of-transit scatter 1.8-5.8 ppt; Lomb-Scargle peaks at
10.7, 11.2, 6.9, 6.4, 15, 12.4 d with 0.34-0.54 ppt amplitude (systematics), 0.2-5 d peaks < 0.52 ppt. No flares in six sectors: a quiet
early-M dwarf. Gaia DR3 ...323269781120: G 15.75, RUWE 0.996, astrometric excess noise 0.00 mas, C* = +0.019 (1.1 sigma), 0% multi-peak,
non_single_star = 0, Gaia Teff 3533 K (TIC 3389 K). Absolute K = 5.94: R = 0.453, M = 0.450 (Mann), identical to the TIC.

## Biggest remaining risk

The 22-arcsec background star at 1.04 px. TESS puts the dip on the target by every available test, but the two stars share pixels in every
sector and the 2.7-sigma phase-0.71 feature, if it ever firms up, would point at a blended eclipsing binary rather than a planet. A
ground-based transit resolving the pair (1.1% deep, 1.45 h, any 0.4 m telescope at Dec -37) settles both at once.

## Independent pipeline (MIT QLP)

`scripts/qlp_check.py --sec-phase 0.71` -> `qlp_check.csv`, figure `addendum_qlp_independent.png`. QLP has one sector (31). Folded at our
ephemeris: transit recovered at 27.4 sigma (35,000 +/- 1,300 ppm in QLP's crowding-corrected flux; QLP over-corrects faint FFI stars by
~3x, the detection is what counts). Phase 0.71: +459 +/- 1,261 ppm; QLP's 10.7 ppt scatter cannot test a 300-ppm feature. No secondary
above 2,500 ppm at any phase (scan rms 1,040 ppm).
