# TIC 259719966: vetting addendum (2026-10-25)

Prepared from the nine TESS sectors on disk (1, 27, 67, 94, 95, 101-104; SPOC FFI cutouts via astrocut S3, cached in
`data/cubes/s<N>_tic259719966_15x15.fits` by `scripts/fetch_cubes.py`). Every number below comes from code run in this
session; the script and output path are given with each test. Intermediate tables are in `data/deepvet_259719966/`.

## Summary

**TIC 259719966 is an eclipsing binary, not a planet candidate.** A secondary eclipse of 9,200 +/- 180 ppm (0.92%) sits at
orbital phase 0.5405, not 0.5, in every one of the nine sectors (2018-2026), in the raw regression-corrected flux before any
detrending. All previous vetting (sweep, followup.py, oddeven_fit.py, the RNAAS text) tested for a secondary only at phase
0.50, where a box of the transit duration measures 20 +/- 540 ppm; the eclipse is 0.04 in phase (10.7 h) away because the
orbit is eccentric (e cos w = 0.0614 +/- 0.0004, e sin w = 0.26 +/- 0.02, e = 0.27 +/- 0.02). The secondary-to-primary depth
ratio gives a TESS-band surface-brightness ratio of 0.256 +/- 0.006, i.e. a companion of about 2,730 K (blackbody, 600-1000
nm) next to the 3,443 K host. A 0.92 R_Jup object at 2,700 K is a late-M dwarf (M7-M8, roughly 0.08-0.10 M_sun) or a young
brown dwarf; it is not a planet. The eccentricity also removes the "density tension": with e free, a/R* = 31.4 +/- 0.8 gives
rho_star = 3.4 +/- 0.25 solar, consistent with the TIC's 4.1 (the circular fit's 5.9 was the signature of e sin w > 0).

What this means for the paperwork: the RNAAS manuscript (AAS82771) describes this object as a planet-sized candidate and
reports "no secondary eclipse"; that statement is wrong and the manuscript needs a correction or withdrawal of this object
before acceptance. Shubham Kanodia should be told before he schedules time: the object is still a worthwhile low-mass
eclipsing binary (an M3.5 + M7/8 pair with a measured eccentricity is a useful radius calibrator), but it is not a GEMS
planet. Expected RV semi-amplitude for a 0.09 M_sun companion is ~14 km/s, so one or two spectra settle the mass.

Test-by-test:

| # | Test | Result | Verdict |
|---|---|---|---|
| 1 | Difference image / centroid, all 9 sectors | dip on target, 0.6-4.6 arcsec from the direct-image position in the 200-s sectors; S1 within 3.9 arcsec; S103 anomaly traced to one quality-flagged ramp | pass |
| 2 | Blend budget, Gaia DR3 within 2.5 px | 44-arcsec neighbour would need an 18.6% eclipse and its own aperture shows < 0.3%; dip-on-neighbour chi2 = 107,500 vs 1,790 for dip-on-target | pass (neighbour excluded) |
| 3 | Per-sector depth vs aperture contents | depth tracks the neighbour's share of the aperture; undiluted sectors scatter 36.7-43.9 kppm (more than the 0.5 kppm formal errors) | pass, with a 15% sector-to-sector depth spread unexplained |
| 4a | Odd/even, 19 transits | +339 +/- 399 ppm (0.85 sigma) | pass |
| 4b | Secondary at phase 0.5 | -158 +/- 158 ppm | pass, but irrelevant |
| 4c | Secondary at any phase | **9,231 +/- 180 ppm at phase 0.5405, 45 sigma, in all 9 sectors** | **fail: eclipsing binary** |
| 4e | Independent pipeline (MIT QLP, 6 sectors) | secondary at phase 0.5405 seen at 22 sigma; phase 0.50 consistent with zero | **confirms the eclipsing binary** |
| 4d | Ellipsoidal / beaming | 3-sigma limits 310 ppm at P/2, 910 ppm at P (red-noise beta 4.4) | inconclusive; expected signals are 5 ppm and ~200 ppm |
| 5 | Shape and density | flat-bottomed (T23/T14 = 0.55 +/- 0.03, grazing excluded at delta-chi2 = 443); circular density 5.9 vs TIC 4.1 reconciled by e >= 0.125, and the measured e = 0.27 does it | pass as a shape; the density mismatch was a real eccentricity signal |
| 6 | Transit timing, 19 transits | O-C rms 0.8 min, chi2/N = 1.07 | pass, no TTV |
| 7 | Host activity | no coherent rotation period (per-sector LS peaks 6-15 d at 1.6-3.8 ppt are scattered-light systematics); 16 flares in 7 sectors, 0.04-0.19 per day, none inside a transit | pass |
| 8 | Pitfall checks (jitter mask, quality flags, S1 bad pointing, duration grid) | no transit cadence lost; no transit in the S1 bad-pointing stretch; S103 epoch 26 is 100% quality-flagged and correctly dropped from the light curve | pass |

Only outside observations can settle: the companion mass (radial velocities), the host's spectroscopic parameters, and
whether the companion is a bottom-of-the-main-sequence star or a brown dwarf.

## 1. Difference image and centroid in every sector

`scripts/deep_vet.py` -> `data/deepvet_259719966/pixel_tests.csv`; `scripts/deep_vet2.py` -> `centroid_vectors.csv`,
`per_transit_diffimage.csv`; figures `addendum_pixels_all_sectors.png`, `addendum_centroid_vectors.png`,
`addendum_sector103.png`.

Per sector, whole-cutout test as in `pipeline.centroid_test` (plane background per frame, noise-normalised difference image,
in = |phase| < T14/2, out = 1-3 T14 away). "diff - direct" is the distance between the difference-image centroid and the
direct-image centroid of the target measured the same way (so a WCS offset cancels; the raw WCS position is itself 0.05-0.31 px
from the direct-image centroid).

| sector | cadence | centroid shift (px) | shift sigma | shift toward neighbour? (cos) | peak SNR | diff - direct (arcsec) |
|---|---|---|---|---|---|---|
| 1 | 30 min | 0.021 | 6.2 | +0.55 | 7.7 | 3.9 |
| 27 | 10 min | 0.022 | 7.8 | +0.30 | 10.1 | 3.1 |
| 67 | 200 s | 0.020 | 11.9 | +0.38 | 14.4 | 4.6 |
| 94 | 200 s | 0.008 | 4.3 | +0.61 | 17.7 | 1.4 |
| 95 | 200 s | 0.022 | 11.3 | +0.54 | 17.6 | 0.8 |
| 101 | 200 s | 0.025 | 10.2 | +0.29 | 20.2 | 1.3 |
| 102 | 200 s | 0.021 | 13.7 | +0.32 | 16.4 | 0.6 |
| 103 | 200 s | 0.131 | 35.0 | -0.52 | 6.4 | 23.6 |
| 104 | 200 s | 0.021 | 11.6 | +0.57 | 17.8 | 1.0 |

Reading: the centroid shift is 0.02 px (0.4 arcsec) in eight sectors and points toward the neighbours, which is the
direction an on-target dip produces (the star that dims is the target, so the light centroid moves toward everything else).
A dip on the 44-arcsec neighbour would move the centroid away from it. The whole-frame prediction for an on-target dip is
0.003-0.009 px; the observed 0.02 px is 2-3x that in every sector, so the simple model underpredicts the size (plane
subtraction and the non-linear clip at zero both bias it), but the direction test is model-independent: the dip-on-target
vector fits better than the dip-on-neighbour vector in all eight clean sectors (chi2 23 vs 45 in S1; 11 vs 39 in S94).

Sector 1: the difference-image peak is 0.66 px (14 arcsec) from the target and 1.56 px from the neighbour; the
difference-image centroid is 3.9 arcsec from the direct-image centroid. Per transit (3-5 cadences each) the centroids scatter
by up to 1 px, so the 30-min data alone cannot place the dip better than about one pixel; it does not point at the neighbour
(the neighbour is 2.17 px away), but on its own it would not exclude it either. The 200-s sectors do: 0.6-4.6 arcsec in
seven sectors.

Sector 103: the 35-sigma centroid shift and the 25-arcsec difference-image offset come entirely from "epoch 26", a predicted
transit at BTJD 4165.83 that falls in the scattered-light ramp at the start of the sector's second orbit. All 17 cadences
inside that window carry non-zero quality flags; the light-curve pipeline removed them (jitter mask), `centroid_test` did
not. Per transit: epoch 25 peaks 0.73 px from target at 14.3 sigma, epoch 27 0.63 px at 14.2 sigma, epoch 26 1.56 px at
"60 sigma" (the ramp). With epoch 26 removed the sector is on target (`addendum_sector103.png`).

Secondary eclipse location (`secondary_location.csv`, quality-flagged cadences removed): the difference image at phase 0.54
peaks 0.13-0.66 px from the target at 5.5-18 sigma in seven sectors (S1 3.8 sigma at 1.0 px, S95 5.5 sigma at 1.2 px),
centroid shift < 2.4 sigma everywhere. The secondary is on the same star as the primary.

## 2. Blend analysis with Gaia DR3 neighbours

`data/gaia_neighbors_259719966.csv` (Gaia archive, 80 arcsec cone); `deep_vet.py` -> `blend_budget.csv`; `deep_vet2.py` ->
`dilution_psf.csv`, `depth_hypotheses.csv/json`; figure `addendum_blend_hypotheses.png`.

Needed eclipse depth = D_obs (1 + f) / f with f the G-band flux relative to the target, i.e. assuming all of that star's
light is in the aperture (the most favourable case for a blend). Sources within 2.5 px:

| Gaia DR3 | sep (arcsec / px) | G | flux / target | needed eclipse | possible? |
|---|---|---|---|---|---|
| ...628619630976 | 19.9 / 0.95 | 22.05 | 0.0015 | 2,709% | no |
| ...254959886976 | 26.5 / 1.26 | 18.27 | 0.048 | 87% | only as a total eclipse of a 1,060 pc star |
| ...632915830656 | 28.0 / 1.33 | 20.80 | 0.0047 | 854% | no |
| ...254958709376 | 29.9 / 1.42 | 20.84 | 0.0045 | 888% | no |
| ...834778061568 | 32.7 / 1.56 | 19.61 | 0.014 | 289% | no |
| ...254958709632 | 34.2 / 1.63 | 21.09 | 0.0036 | 1,115% | no |
| ...461117434624 | 35.2 / 1.67 | 21.14 | 0.0034 | 1,165% | no |
| ...250662504576 | 37.7 / 1.79 | 21.32 | 0.0029 | 1,377% | no |
| ...254959887872 | 43.7 / 2.08 | 16.38 | 0.274 | 18.6% | yes |
| ...839074264320 | 49.4 / 2.35 | 21.03 | 0.0038 | 1,054% | no |

The only capable blend sources are the 44-arcsec neighbour (G = 16.38, bp-rp 2.40, parallax 4.96 mas: another M dwarf at
203 pc, unrelated to the target at 147 pc) and, marginally, the G = 18.27 star at 26.5 arcsec. Both are tested directly by
the aperture experiment: a 1-px-radius aperture centred on the 44-arcsec neighbour (PSF model: 48-63% of its light is the
neighbour's, 7-30% the target's wings) shows a dip of 700-3,300 ppm in eight sectors and 9,000 ppm in S104, where the
apertures overlap. If the neighbour were eclipsing by 18.6%, that aperture would dip by 9-12%. Fitting all 27 depth
measurements (three apertures x nine sectors) with each star's PSF-modelled share of the aperture light: dip on target gives
D = 4.12% with chi2 = 1,787 (390 with a 10% per-sector depth jitter); dip on neighbour gives chi2 = 107,500 (1,794 with the
jitter). The phase-folded pipeline aperture with the neighbour's pixels excluded gives 40,086 +/- 197 ppm versus 37,102 +/-
189 ppm with them included: excluding the neighbour deepens the dip, which is what an on-target dip does.

## 3. Per-sector depth and aperture contents

`deep_vet.py` -> `aperture_contents.csv`, `depth_per_sector.csv`; figure `addendum_depth_timing_activity.png` (left).
Depths are fitted with the geometry fixed at the joint fit and only Rp/Rs free; dilution is the summed G-band flux of
non-target Gaia sources whose central pixel lies in the aperture (the PSF-wing version is in `dilution_psf.csv`).

| sector | cadence | aperture px | neighbour in aperture? | dilution (point) | depth, pipeline (ppm) | depth, neighbour excluded (ppm) |
|---|---|---|---|---|---|---|
| 1 | 30 min | 4 | no (1.46 px outside) | 0.000 | 37,634 +/- 466 | same aperture |
| 27 | 10 min | 3 | no | 0.056 | 36,728 +/- 417 | same |
| 67 | 200 s | 5 | no (0.50 px outside) | 0.011 | 33,771 +/- 496 | 40,066 +/- 527 |
| 94 | 200 s | 3 | no | 0.053 | 42,823 +/- 498 | same |
| 95 | 200 s | 2 | no | 0.005 | 40,749 +/- 765 | same |
| 101 | 200 s | 6 | yes | 0.333 | 24,137 +/- 333 | 32,031 +/- 398 |
| 102 | 200 s | 5 | yes | 0.281 | 33,808 +/- 356 | 40,392 +/- 410 |
| 103 | 200 s | 4 | no | 0.002 | 41,306 +/- 516 | same |
| 104 | 200 s | 5 | yes | 0.286 | 33,413 +/- 442 | 43,873 +/- 559 |

Depth tracks dilution, not time: the three shallow pipeline depths are exactly the three sectors with the neighbour inside
the aperture, and S67's 1-px aperture change moves it from 33.8 to 40.1 kppm. Two honest residuals: sector 101 (camera 4,
19 ppt scatter, the noisiest sector) stays 20% shallow even with the neighbour's pixels removed, and the neighbour-excluded
depths of the other eight sectors span 36.7-43.9 kppm against ~0.5 kppm formal errors. A 15% sector-to-sector depth spread is
more than detrending noise usually gives; PSF-wing dilution that changes with the per-sector PSF and pointing is the likely
cause (the TIC contamination ratio of 0.007 badly understates the crowding here), spot coverage is another. It does not
change the verdict either way.

## 4. Odd/even, secondary eclipse at any phase, phase curve

`deep_vet.py` -> `stitched_tests.json`, figure `addendum_oddeven_secondary.png`; `deep_vet2.py` -> `secondary_per_sector.csv`;
`deep_vet4.py` -> `joint_eccentric_fit.json`, `joint_eccentric_chain.npz`, figure `addendum_secondary_eclipse.png`;
`deep_vet.py` -> `phase_curve.json`.

Odd/even on the stitched nine-sector light curve (19 transits), neighbour-excluded aperture: odd 40,229 +/- 258, even 39,890
+/- 304, difference +339 +/- 399 ppm (0.85 sigma). Pipeline aperture: +677 +/- 382 (1.8 sigma). Pass.

Secondary at phase 0.50 (box of one transit duration): -158 +/- 158 ppm. Pass, and meaningless, because:

Secondary at any phase: sliding the same box from phase 0.06 to 0.95 in steps of 0.005 gives a flat scan with rms 660-740
ppm everywhere except a 45-sigma dip at phase 0.540 (box depth 6,700-7,500 ppm, under-estimated because the box is narrower
than the eclipse), flanked by the shoulders the 33-h smoother puts around any unmasked dip. Per sector, the box at phase 0.54
measures 4,600-8,300 ppm in all nine sectors (`secondary_per_sector.csv`), and the dip is visible in the raw
regression-corrected flux of every sector before flattening.

Joint eccentric fit of primary and secondary (batman, both events masked before flattening, sectors 1, 27, 67, 94, 95,
102, 103, 104 with the neighbour's pixels excluded, 6,723 points, emcee 36 walkers x 2,500 steps, free quadratic limb darkening with a 0.1 prior, chi2 = 6,947):

| quantity | value |
|---|---|
| secondary eclipse depth | 9,242 +/- 168 ppm |
| secondary phase | 0.54045 +/- 0.00006 |
| e cos w | 0.0614 +/- 0.0004 |
| e sin w | 0.260 +/- 0.019 |
| e | 0.267 +/- 0.018 |
| w | 76.7 +/- 1.0 deg |
| Rp/Rs | 0.1908 +/- 0.0013 |
| a/Rs | 31.3 +/- 0.75 |
| b | 0.355 +0.013 -0.018 |
| primary depth (Rp/Rs)^2 | 36,390 +/- 510 ppm |
| surface-brightness ratio, TESS band | 0.256 +/- 0.005 |
| companion T from blackbody ratio, 600-1000 nm | 2,730 +/- 10 K (statistical only; +/- 100-200 K systematic from the blackbody assumption) |
| companion radius | 0.0944 R_sun = 0.918 R_Jup |
| rho_star from a/Rs (Kepler's law, e free) | 3.37 +/- 0.25 solar (TIC: 4.07) |

Interpretation: a secondary eclipse of 0.92% from a 0.19 R_* companion requires the companion to radiate a quarter of the
host's surface brightness in the TESS band. Reflected light from a planet at a/R* = 31 is at most (Rp/a)^2 = 4 x 10^-5 of
the stellar flux, 40 ppm; thermal emission at T_eq ~ 400 K is zero in the TESS band. The companion is a star or a brown
dwarf. At 2,700 K and 0.92 R_Jup the natural match is an M7-M8 dwarf of 0.08-0.10 M_sun; a brown dwarf younger than ~1 Gyr
can also be this warm and this large. The host remains a 0.49 M_sun mid-M dwarf and the system is an eccentric (e = 0.27)
11.08-day low-mass eclipsing binary.

Independent check with MIT's QLP light curves (`data/qlp/`, `data/deepvet_259719966/qlp_check.csv`, figure
`addendum_qlp_independent.png`): QLP publishes FFI photometry for this star in sectors 67, 95, 101-104, extracted with a
pipeline that shares no code, apertures or detrending with tess-hunt. Folding QLP's own flux at our ephemeris gives a
secondary at phase 0.5405 in every one of the six sectors (16,500-21,200 ppm per sector; 19,420 +/- 870 ppm combined, 22
sigma) and nothing at phase 0.50 (-1,000 to +2,500 ppm, all within 1 sigma of zero). QLP's depths are larger than ours
(primary 59-72 kppm vs our 36-44 kppm) because QLP applies its own crowding correction, but the secondary-to-primary ratio is
the same to within the errors (0.31 vs 0.25). The secondary eclipse is not an artifact of this repository's photometry.

Phase-curve limits (unflattened, transit-masked, per-sector cubic per orbit half fitted jointly with sinusoids; red-noise
inflation beta = 4.4): amplitude at P = 1,050 +/- 430 ppm (not significant), at P/2 = 360 +/- 150 ppm; 3-sigma limits 910
ppm at P and 310 ppm at P/2, optimistic by up to 2x because the per-sector polynomials absorb 11-day power. Expected
ellipsoidal amplitude for a 0.09 M_sun companion at a/R* = 31 is ~10 ppm and beaming ~200 ppm: TESS FFIs cannot test either
here. Inconclusive.

## 5. Transit shape and stellar density

`scripts/deep_vet3.py` -> `data/deepvet_259719966/shape_density.json`; figure `addendum_shape_density.png`. From the
existing circular-orbit chain (`plots/candidates/followup/s1_tic259719966_mcmc.npz`): T14 = 2.369 +/- 0.023 h, T23 = 1.30
+/- 0.07 h, T12 = 0.53 +/- 0.04 h, T23/T14 = 0.55 +/- 0.03; the posterior puts zero weight on b > 0.8 and zero on grazing (b >
1 - Rp/Rs). Forcing a grazing geometry (Rp/Rs = 0.35, b = 1) costs delta-chi2 = +443 on 2,307 points; the dip is flat-bottomed
with a 32-min ingress, consistent with a 0.19-0.20 radius ratio and not V-shaped. The circular fit's density, 5.93 +0.43
-0.37 solar, is 4.5 sigma above the TIC's 4.07; reconciling them needs either e >= 0.125 +/- 0.02 (w = 90 deg) or a stellar
radius of 0.436 R_sun at the TIC mass. A re-fit with the TIC density as a 25% prior barely moves (rho = 5.7 +/- 0.35, b =
0.57, Rp/Rs = 0.1965, chi2 unchanged), so the data, not the prior, set the density. The joint eccentric fit in section 4
resolves it: e sin w = 0.26 shortens the transit, and the true a/R* = 31.4 gives 3.4 +/- 0.3 solar.

## 6. Transit timing

`deep_vet.py` -> `transit_times.csv`; figure `addendum_depth_timing_activity.png` (middle). Nineteen individual transit
times (shape fixed, per-transit errors from the Jacobian scaled by reduced chi2): O-C rms 0.8 min, largest |O-C| 2.0 +/- 0.8
min (S104 epoch 29), chi2/N = 1.07 for N = 19 over 2,869 days. No TTV. (S103 "epoch 26" has no measurement: its window is
the quality-flagged ramp.)

## 7. Host activity

`deep_vet.py` -> `activity.csv`; `deep_vet2.py` -> `rotation_short.csv`, `flares.csv`; figure
`addendum_depth_timing_activity.png` (right). Lomb-Scargle on the unflattened, transit-masked, regression-corrected flux
peaks at a different period in every sector (6.3, 7.4, 8.6, 9.3, 9.4, 9.7, 9.7, 15, 15 d) with 1.6-3.8 ppt amplitude:
scattered-light systematics at the TESS orbit timescale, not a rotation signal. Between 0.2 and 5 d the strongest
per-sector peaks are 0.4-3.3 ppt and again inconsistent in period. No rotation period is detected; a spotted fast rotator
with > 3 ppt amplitude is excluded. Flares (>= 2 consecutive points more than 4 sigma high): 16 events in the seven 200-s
and 10-min sectors, 0.04-0.19 per day, peak amplitudes 1.3-5.2% (0.3-1.4 transit depths), none inside a transit or an
eclipse window. Spot crossings: the per-transit depths within a sector agree to 1-2 sigma (S1 37.4/36.8/39.3 kppm; S94
43.2/42.5; S104 34.0/32.8), so spots do not visibly bias individual transits; the sector-to-sector spread in section 3 is
the larger effect.

## 8. Pitfall-log checks for a 30-min first detection

`deep_vet2.py` -> `transit_cadence_audit.csv`. For every predicted transit: cadences in the cube, after clean(), after the
jitter mask, and quality-flag count. No transit lost a cadence to the jitter mask except S103 epoch 26 (17 of 17 cadences
quality-flagged and removed, the ramp above). No sector-1 transit falls in the 1347-1349.6 bad-pointing stretch (S1
transits at 1330.07, 1341.14, 1352.22). The sector-1 BLS duration of 1.75 h came from the coarse 0.25-h grid; the fitted T14
is 2.37 h (circular) and the 30-min integration is modelled by super-sampling in every fit. The sector-1 "anomaly" verdict in
`report.md` came from predicting the centroid shift with neighbours inside 2.0 px while the dominant neighbour sits at 2.08
px; the whole-frame prediction and the direction test above replace it.

## What remains for outside observers

- One or two spectra (K ~ 14 km/s for a 0.09 M_sun companion; ~2 km/s if it were a 13 M_Jup brown dwarf) give the companion
  mass and settle star versus brown dwarf.
- A spectroscopic Teff, metallicity and radius for the host, which the companion radius inherits.
- A ground-based secondary eclipse (0.9%, ~3 h, phase 0.5405, ephemeris in `observer_brief.md` with the secondary at
  mid-transit + 5.987 d) confirms the eclipse on the target with the 44-arcsec neighbour excluded.

## New pitfalls (added to CLAUDE.md this session)

- Secondary tests only at phase 0.5 miss eccentric eclipsing binaries; scan all phases.
- centroid_test must drop quality-flagged / jitter-masked cadences and be run per transit.
- Predict centroid shifts from every Gaia star in the frame and compare direction, not only size.

## Scripts and outputs

`scripts/fetch_cubes.py`, `tesshunt/cubes.py` (cube cache); `scripts/deep_vet.py`, `deep_vet2.py`, `deep_vet3.py`,
`deep_vet4.py` (logs in `data/deepvet*_259719966.log`; tables in `data/deepvet_259719966/`). Figures in this folder:
`addendum_pixels_all_sectors.png`, `addendum_centroid_vectors.png`, `addendum_sector103.png`,
`addendum_blend_hypotheses.png`, `addendum_depth_timing_activity.png`, `addendum_oddeven_secondary.png`,
`addendum_shape_density.png`, `addendum_secondary_eclipse.png`, `addendum_qlp_independent.png`.
