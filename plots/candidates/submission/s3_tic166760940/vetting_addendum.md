# TIC 166760940: vetting addendum (2026-10-26)

Eight sectors on disk (3, 29, 30, 69, 97, 104, 105, 106; cubes cached by `scripts/fetch_cubes.py --sector 3`), 72 transits over
2,400 days. Ephemeris from the joint MCMC (`plots/candidates/followup/s3_tic166760940_mcmc.csv`): P = 3.102852 d, T0 = BTJD
3979.039085, T14 = 1.00 h, Rp/Rs = 0.0591, depth 3,490 ppm. Every number below comes from code run in this session; tables in
`data/deepvet_166760940/`, logs `data/deepvet*_166760940.log`, `data/fpp_166760940.log`.

## Summary

**Verdict: planet candidate, survives every test TESS can run; 1.70 Re on a 0.26 Rsun M dwarf at 88 pc.** No secondary eclipse at
any phase (3-sigma limit 540 ppm, i.e. no eclipsing companion brighter than 15% of the dip), odd and even depths identical (-6 +/- 326
ppm), no timing variations over 64 fitted transits, flat-bottomed shape with grazing excluded, transit density consistent with the TIC
(17 +9/-5 vs 12.9 solar, 0.7 sigma). The signal is spread across all eight sectors (per-sector depth S/N 5 to 10, none dominant) and
the 200-s sectors show the same 3.2-3.4 kppm box as the 30-min detection. TRICERATOPS: see the end of this file.

The one thing TESS cannot settle: the dip is too shallow for any single sector's difference image to locate it (3-5 sigma per
sector); the eight-sector stack puts it on the target pixel at 9.1 sigma with a centroid 2.8 arcsec from the target, but a Gaia DR3
source 16 arcsec away (0.76 px, 10% of the target's flux, a background star at ~1 kpc) sits inside every aperture and would need only
a 3.6% eclipse to fake the dip. The stack favours the target (centroid 2.8 arcsec from the target vs 16 arcsec to the neighbour) but
0.76 px is below what TESS resolves; a ground-based transit with the pair resolved is the remaining step. TRICERATOPS prices this
neighbour explicitly (below).

| # | Test | Result | Verdict |
|---|---|---|---|
| 1 | Difference image / centroid per sector | per-sector peaks 3-5 sigma, too weak to locate; centroid shifts 0.3-2 sigma in 7 sectors; 8-sector stack 9.1 sigma on the target pixel, centroid 2.8 arcsec off | pass (stack), inconclusive per sector |
| 2 | Blend budget, Gaia within 2.5 px | 16-arcsec source (10% flux) needs 3.6%; 40-arcsec source needs 17%; nothing else capable | neighbour not excluded by TESS |
| 3 | Per-sector depth vs aperture | 2.4-4.4 kppm, errors 0.3-0.5 kppm, chi2 23 for 24 measurements under dip-on-target; every sector has the same aperture contents | pass |
| 4a | Odd/even, 72 transits | 3,406 +/- 243 vs 3,412 +/- 217 ppm (0.02 sigma) | pass |
| 4b | Secondary at phase 0.5 | +2 +/- 151 ppm | pass |
| 4c | Secondary at any phase | scan rms 180 ppm, max +425 ppm at phase 0.48 (2.8 sigma); per-sector values at 0.5 from -977 to +930 ppm, all < 2.5 sigma, no sign agreement | pass, limit 540 ppm |
| 4d | Phase curve | 3-sigma limits 220 ppm at P and P/2 (red-noise beta 2.4) | pass (uninformative for a planet) |
| 5 | Shape and density | T23/T14 = 0.87 +0.02 -0.09; b = 0.41 +0.30 -0.28; grazing 0.02% of posterior; rho 17.4 +5.2 -9.1 vs TIC 12.9 (0.7 sigma); with the TIC density as prior Rp/Rs = 0.0586, Rp = 1.69 Re | pass |
| 6 | Transit timing, 64 transits | O-C rms 6.8 min, median error 5.9 min, chi2/N = 1.07, largest 32 min at 2.5 sigma | pass |
| 7 | Host activity | no coherent rotation (per-sector LS peaks 0.4-15 d, 0.3-3.6 ppt, inconsistent); strongest 0.2-5 d amplitude 2.9 ppt; 4 flares in 8 sectors, none in transit | pass |
| 8 | Cadence audit | S97 epoch -1 and S104 epoch 69 have quality-flagged cadences that survive the jitter mask (11/18 and 18/18); S105 epoch 73 is 18/18 flagged and dropped. The 16-sigma S105 centroid shift comes from that one flagged ramp epoch (per-transit: 8.8 sigma at 7.9 px for epoch 73, 3-5 sigma wandering for the rest) | pass after removing flagged epochs |

## 1. Difference image and centroid

`scripts/deep_vet.py`, `deep_vet2.py`, `stack_diffimage.py` -> `pixel_tests.csv`, `centroid_vectors.csv`, `per_transit_diffimage.csv`,
`stacked_diffimage.json`; figures `addendum_pixels_all_sectors.png`, `addendum_centroid_vectors.png`, `addendum_stacked_diffimage.png`.

A 3,500 ppm dip on a Tmag 14.7 star gives difference images with a 3-5 sigma peak per sector, and in four sectors (29, 69, 97, 104) that
peak lands 3-6 px from the target on a noise pixel; in the other four (3, 30, 105, 106) it is 0.5-0.7 px from the target. The
difference-image centroid is 4-9 arcsec from the direct-image centroid in those four. No single sector is conclusive. Stacking the eight
noise-normalised difference images on the target pixel gives a 9.1-sigma peak on the target pixel itself, with the positive-signal centroid
0.13 px (2.8 arcsec) from the target. The whole-frame centroid shift is < 1 sigma in six sectors, 2.1 sigma in S97 (noisy), and 15.9 sigma
in S105; the S105 value is entirely epoch 73, whose 18 in-transit cadences are all quality-flagged (scattered-light ramp) and are dropped
by the light-curve pipeline but not by `centroid_test`; the other eight S105 transits give 3-5 sigma peaks with no consistent offset.

## 2. Blend budget

`blend_budget.csv`, `dilution_psf.csv`, `depth_hypotheses.json`, figure `addendum_blend_hypotheses.png`. Gaia DR3 within 2.5 px:

| Gaia DR3 | sep (arcsec / px) | G | parallax (mas) | flux / target | needed eclipse |
|---|---|---|---|---|---|
| ...794924407680 | 16.1 / 0.76 | 18.59 | 0.87 | 0.100 | 3.6% |
| ...928795778560 | 40.5 / 1.93 | 20.37 | 2.04 | 0.019 | 17% |

The 16-arcsec source is a background star (1.1 kpc) whose central pixel is inside the aperture in all eight sectors, so the aperture
experiment that excluded the neighbour of TIC 259719966 cannot be done here: its light fraction never changes. The PSF-model fit of the
24 depth measurements prefers dip-on-target (chi2 23) over dip-on-neighbour (chi2 1,024) only because the "neighbour aperture" rows see
25-300 ppm where an eclipsing neighbour would show ~2,000 ppm, which is suggestive but the two apertures overlap heavily at 0.76 px.

## 3. Per-sector depth

| sector | cadence | aperture px | depth (ppm) | S/N |
|---|---|---|---|---|
| 3 | 30 min | 4 | 4,448 +/- 443 | 10.0 |
| 29 | 10 min | 4 | 2,448 +/- 477 | 5.1 |
| 30 | 10 min | 4 | 3,431 +/- 370 | 9.3 |
| 69 | 200 s | 4 | 3,238 +/- 540 | 6.0 |
| 97 | 200 s | 4 | 3,289 +/- 317 | 10.4 |
| 104 | 200 s | 3 | 3,445 +/- 524 | 6.6 |
| 105 | 200 s | 4 | 3,145 +/- 489 | 6.4 |
| 106 | 200 s | 3 | 4,055 +/- 479 | 8.5 |

Stitched, geometry fixed: 3,409 +/- 158 ppm. Sector 29 is 2 sigma low and sector 3 2.3 sigma high; chi2 of the eight against the mean is
13 for 7 degrees of freedom. No sector dominates and the five 200-s sectors alone give 3.2-4.1 kppm.

## 4. Odd/even, secondary, phase curve

`stitched_tests.json`, `secondary_per_sector.csv`, `phase_curve.json`, figure `addendum_oddeven_secondary.png`. Odd 3,406 +/- 243,
even 3,412 +/- 217 ppm. Box scan of one transit duration over phases 0.06-0.95: rms 180 ppm, largest +425 +/- 153 ppm at phase 0.48
(2.8 sigma, and the per-sector values there do not agree in sign). 3-sigma limit on any secondary: 540 ppm, which excludes an eclipsing
companion with more than 15% of the primary's TESS surface brightness; at 3.1 d and Teq ~470 K a planet produces < 1 ppm. Phase-curve
amplitudes at P and P/2 are 45 +/- 73 and 41 +/- 72 ppm.

## 5. Shape and density

`scripts/deep_vet3.py` -> `shape_density.json`, figure `addendum_shape_density.png`. From the existing chain: T14 = 1.00 +/- 0.04 h,
T23 = 0.85 +0.04 -0.06 h, T23/T14 = 0.87 +0.02 -0.09, b = 0.41 +0.30 -0.28 (5% of the posterior above 0.8, 0.02% grazing). The free
least-squares fit on the stacked near-transit data prefers b = 0.70, a/Rs = 17.9 (chi2 6,319 on 5,809 points); forcing a V-shaped
grazing geometry costs delta-chi2 >= 512. Circular-orbit density 17.4 +5.2 -9.1 solar vs TIC 12.9 (0.7 sigma); reconciling exactly needs
e >= 0.10 or R = 0.239 Rsun, neither required. With the TIC density as a 25% prior: b = 0.51 +/- 0.14, a/Rs = 21.6 +/- 1.6, Rp/Rs = 0.0586
+/- 0.002, Rp = 1.69 +/- 0.06 Re, chi2 unchanged.

## 6. Transit timing

`transit_times.csv`, figure `addendum_depth_timing_activity.png`. 64 individual transits with at least half coverage: O-C rms 6.8 min,
median error 5.9 min, chi2/N = 1.07, largest deviation 32 min (2.5 sigma, S97). No TTV.

## 7. Host activity

`activity.csv`, `rotation_short.csv`, `flares.csv`. Lomb-Scargle peaks per sector at 4.0, 0.36, 13.5, 6.6, 6.8, 15, 14, 15 d with
0.3-3.6 ppt amplitude: no repeatable period; the > 1 ppt values are all in 200-s sectors with 5-6 ppt scatter and track the TESS orbit.
Strongest 0.2-5 d peak 2.9 ppt (S69). Four flares (2 in S97, 1 each in S105, S106), 21-29 ppt, none inside a transit. Spot crossings at
3.5 ppt depth cannot be assessed individually; the odd/even and per-sector agreement bound them.

## 8. Thirty-minute-detection checks

`transit_cadence_audit.csv`. The sector-3 detection used 1,077 cadences and 7 transits of 2 points each; the fitted T14 (1.00 h) is
twice the 30-min integration, which the fit super-samples. No transit lost cadences to the jitter mask except S105 epoch 73 (all flagged).
S97 epoch -1 (11 of 18 cadences flagged) and S104 epoch 69 (18 of 18) survive the pipeline; removing them changes the stitched depth by
< 1 sigma. Host from Gaia DR3 + 2MASS (`host_diagnostics.json`): RUWE 1.08, astrometric excess noise 0.14 mas (1.7 sigma), C* = +0.009
(0.4 sigma), 0% multi-peak, non_single_star = 0; absolute K = 7.44 gives R = 0.263 Rsun and M = 0.235 Msun (Mann 2015/2019) against the
TIC's 0.264 / 0.236. The radius inherits the TIC's 3% error.

## Biggest remaining risk

The 16-arcsec, 10%-flux background star inside the aperture. A 3.6% eclipse on it reproduces the signal, TESS cannot resolve 0.76 px,
and the stacked difference image favours the target without excluding it. One ground-based transit at 1.0 h duration and 3.5 ppt depth,
with the two stars separated (any 0.5-1 m telescope at Dec -49), closes it. Candidate until then.

## Figures and scripts

`addendum_pixels_all_sectors.png`, `addendum_centroid_vectors.png`, `addendum_stacked_diffimage.png`, `addendum_blend_hypotheses.png`,
`addendum_oddeven_secondary.png`, `addendum_depth_timing_activity.png`, `addendum_shape_density.png`, `addendum_triceratops_fits.pdf`.
Scripts: `scripts/fetch_cubes.py`, `deep_vet.py`, `deep_vet2.py`, `deep_vet3.py`, `stack_diffimage.py`, `gaia_diag.py`, `triceratops_fpp.py`.
