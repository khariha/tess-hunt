# Three transiting candidates around faint M dwarfs from an agent-driven search of TESS full-frame images

**Kishore Hariharan** (Generative Experiences, Orlando, FL, USA)

*Draft for Research Notes of the AAS. Limit: 1,500 words, one figure OR one table. Word count target: 1,300. Figure: three phase-folded transits with the fitted model (one panel each).*

## Abstract

We report three transiting planet candidates around M dwarfs with TESS magnitudes 13.8 to 14.8, identified in an independent box-least-squares search of SPOC full-frame image (FFI) cutouts of 38,620 Sector 1 stars with Tmag 12 to 15 and radii below 0.6 Rsun that were not 2-minute targets and not listed as TESS Objects of Interest. Each signal is recovered independently in every available sector (3 to 9 sectors, 2018 to 2026), is located on the target in noise-normalized difference images, shows consistent odd and even depths and no secondary eclipse, and is fit with a physical transit model. TIC 389723320 hosts a 4.6 Rearth candidate at P = 3.927167 d; TIC 114949540 a 4.35 Rearth candidate at P = 2.739176 d; and TIC 259719966 a 0.94 Rjup companion at P = 11.077205 d whose nature (planet, brown dwarf, or late M dwarf) requires a mass measurement. None appears in the TOI, Community Candidate, or confirmed-planet lists. The search, vetting, and fitting pipeline was built and operated with an AI coding agent under human direction; code and per-candidate reports are public.

## 1. Motivation

The TESS mission's own transit searches of FFI data have concentrated on stars brighter than Tmag ~13.5 (Kunimoto et al. 2022; Huang et al. 2020), leaving the fainter M-dwarf population, where a given planet produces the deepest transits, comparatively under-searched. We targeted that population directly: M dwarfs (Teff < 4000 K, R < 0.6 Rsun) with 12 < Tmag < 15, excluding 2-minute targets, TOIs, Community Candidates, and known planet hosts. For Sector 1 this yields 40,328 unique stars, of which 38,620 with contamination ratio below 1 were searched.

## 2. Data and search

For each star we extracted a 15x15 pixel cutout from the SPOC calibrated FFI sector cubes hosted on AWS (astrocut; Brasseur et al. 2019), defined a threshold aperture, summed the aperture flux, and removed shared systematics by regressing against the leading principal components of the non-aperture pixels (lightkurve; Lightkurve Collaboration 2018). After normalization, a 33-hour Savitzky-Golay flatten, and a point-to-point-jump filter for pointing excursions, we ran a box least squares search (Kovács et al. 2002) over 0.3 to 13 days with durations 0.5 to 6 hours. Candidates required SNR >= 9, depth >= 3000 ppm, at least two transits, duration no more than 1.5x the maximum allowed by the TIC stellar density (Stassun et al. 2019), duty cycle below 15%, odd and even depths consistent within 3 sigma, and no secondary eclipse at phase 0.5 above 3 sigma. Twenty-two stars passed.

## 3. Vetting

Each flagged star was tested against Gaia DR3 (Gaia Collaboration 2023) neighbors within 2 pixels, with the eclipse depth each neighbor would need to reproduce the signal; a flux-weighted centroid shift between in- and out-of-transit cadences; and a difference image (out-of-transit mean minus in-transit mean, per-frame plane-subtracted, divided by its per-pixel noise) whose peak locates the source of the dip (Bryson et al. 2013; Twicken et al. 2018). Of the 22, seven were located on a neighbor, one was a spotted rotator, one an unresolvable 1-arcsecond blend, and three showed anomalous centroid behavior. Every survivor was then searched independently in each sector in which TESS observed it (sectors from tess-point; Burke et al. 2020) and the sectors stitched to refine the period. Two high-SNR single-sector detections failed to repeat and were dropped. Three candidates passed all tests.

## 4. Candidates

Transit parameters were obtained with batman (Kreidberg 2015) and emcee (Foreman-Mackey et al. 2013), fitting epoch, period, Rp/Rs, a/Rs, impact parameter, and quadratic limb darkening with Gaussian priors (0.3 ± 0.1), super-sampling each sector to its exposure time. Stellar parameters are from TIC v8.

**TIC 389723320** (Tmag 14.52, Teff 3271 K, R = 0.37 Rsun, d = 126 pc). Six sectors (1, 28, 68, 95, 105, 106), 36 transits. P = 3.927167 ± 0.000001 d, T0 = 2460899.3895 ± 0.0005 BJD, depth 12,850 ± 690 ppm, duration 1.50 ± 0.06 h, b = 0.74, Rp = 4.62 ± 0.13 Rearth, Teq ≈ 565 K. Odd minus even depth 1,100 ± 520 ppm (2.1 sigma). Difference images in sectors 68, 95, and 105 peak on the target at 14, 14, and 18 sigma. Transit-derived stellar density 4.1 (+1.5/-1.0) rho_sun vs 7.0 from TIC, consistent within the b-a/Rs degeneracy.

**TIC 114949540** (Tmag 14.78, Teff 3320 K, R = 0.42 Rsun, d = 168 pc). Three sectors (1, 27, 105), 25 transits. P = 2.739176 ± 0.000001 d, T0 = 2461212.9505 ± 0.0007 BJD, depth 8,930 ± 710 ppm, duration 1.55 ± 0.07 h, b = 0.53, Rp = 4.35 ± 0.17 Rearth, Teq ≈ 670 K. Odd minus even -240 ± 450 ppm. Difference images in sectors 27 and 105 peak on the target at 15 and 8 sigma.

**TIC 259719966** (Tmag 13.80, Teff 3443 K, R = 0.49 Rsun, d = 147 pc). Nine sectors (1, 27, 67, 94, 95, 101, 102, 103, 104), 20 transits. P = 11.077205 ± 0.000001 d, T0 = 2460877.8231 ± 0.0002 BJD, depth 38,300 ± 1,200 ppm (undiluted sectors), duration 2.37 ± 0.03 h, b = 0.55, Rp = 10.58 ± 0.17 Rearth (0.94 Rjup), Teq ≈ 400 K. Odd minus even -230 ± 370 ppm; secondary 20 ± 540 ppm. Difference images in sectors 94 and 104 peak on the target at 33 and 37 sigma. A Gaia neighbor 2.1 px away (27% of the target's G-band flux) enters the automatic aperture in three sectors and dilutes the depth there; the sense of that effect (shallower when included) independently places the dip on the target. With a radius near Jupiter's, the companion may be a giant planet, a brown dwarf, or a star near the hydrogen-burning limit; its radial-velocity semi-amplitude would exceed 100 m/s for any stellar-mass companion.

A residual search of the stitched light curves of the first two systems, with the known transits masked, found no additional signal above ~1,000 ppm with three or more transits at periods below 30 days.

## 5. Caveats and availability

All three hosts are faint for radial-velocity follow-up; TIC 259719966 is the most tractable. None of the candidates is a TOI, Community Candidate (ExoFOP list of 2026-10-24), or known planet host; none has an ADS literature match; and none appears within 30 arcsec of an entry in the DIAmante FFI candidate catalog (Montalto et al. 2020), the QLP faint-star TOI release (Kunimoto et al. 2022), the TESS eclipsing binary catalog (Prša et al. 2022), or the Gaia DR3 variability classification. The pipeline (target selection, parallel FFI sweep, vetting, multi-sector recovery, MCMC fitting, and report generation) was written and operated with Claude, an AI coding agent, directed by the author, who has no prior background in exoplanet science; we note this because the result bears on how such searches can be scaled. Code, pitfalls log, and per-candidate vetting reports: https://github.com/khariha/tess-hunt.

## Figure caption

Phase-folded TESS FFI photometry (all sectors, gray) with the posterior-median transit model (red) for TIC 389723320 (left), TIC 114949540 (center), and TIC 259719966 (right). Binned points in black.

## References (to be formatted)

Brasseur, C. E., et al. 2019, Astrocut, ascl:1905.007. Bryson, S. T., et al. 2013, PASP, 125, 889. Burke, C. J., et al. 2020, tess-point, ascl:2003.001. Foreman-Mackey, D., et al. 2013, PASP, 125, 306. Gaia Collaboration 2023, A&A, 674, A1. Huang, C. X., et al. 2020, RNAAS, 4, 204. Kovács, G., Zucker, S., & Mazeh, T. 2002, A&A, 391, 369. Kreidberg, L. 2015, PASP, 127, 1161. Kunimoto, M., et al. 2022, ApJS, 259, 33. Lightkurve Collaboration 2018, ascl:1812.013. Montalto, M., et al. 2020, MNRAS, 498, 1726. Prša, A., et al. 2022, ApJS, 258, 16. Stassun, K. G., et al. 2019, AJ, 158, 138. Twicken, J. D., et al. 2018, PASP, 130, 064502.
