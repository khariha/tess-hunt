# Community TOI candidate: TIC 7586789

Submitter: Kishore Hariharan (Generative Experiences / tess-hunt). Prepared 2026-10-25.

## Host star (TIC v8)
Tmag 14.54, Teff 3389 K, R = 0.453 Rsun, M = 0.450 Msun, distance 172.2 pc, contamination ratio 0.009.
RA 63.65490, Dec -37.02509. Sector 4 camera 3 CCD 1.

## Signal (MCMC, all sectors, batman + emcee, free quadratic limb darkening)
| quantity | value |
|---|---|
| Period (d) | 2.261386 +0.000000 -0.000000 |
| Epoch (BJD) | 2461248.48398 +0.00038 -0.00036 |
| Depth (ppm) | 10594 +551 -572 |
| Duration (h) | 1.45 +0.04 -0.05 |
| Rp/Rs | 0.1045 +0.0020 -0.0032 |
| a/Rs | 10.98 +1.56 -0.92 |
| Impact parameter | 0.60 +0.08 -0.21 |
| Planet radius (Rearth) | 5.10 +0.13 -0.14 |
| Planet radius (Rjup) | 0.454 +0.012 -0.012 |
| a (AU) | 0.0242 +0.0043 -0.0024 |
| T_eq (K, albedo 0) | 708 +37 -56 |
| Stellar density from transit (solar) | 3.95 +2.51 -1.05 |
| Stellar density from TIC (solar) | 4.84 |
| Limb darkening u1, u2 | 0.27 +0.08 -0.08, 0.28 +0.09 -0.09 |
| Transit SNR (depth / error) | 18.5 |

## Detection
Independent box-least-squares search of FFI cutouts (SPOC calibrated FFIs via MAST S3 cubes, astrocut) of TIC M dwarfs with Tmag 12-15 and R < 0.6 Rsun not on the 2-minute target list. Sector-4 detection: P = 2.26186 d, depth 8898 ppm, SNR 22.9, 9 transits.

## Multi-sector recovery
[4, 5, 31, 32, 106, 107]; usable 6; signal recovered in 6: [4, 5, 31, 32, 106, 107]

Per-sector independent searches (narrow period window) and joint refinement are in the follow-up figure.

## Vetting
- Centroid / difference image (sector 4): likely_on_target; shift 1.0 sigma; difference-image peak 1.3 px from target at 3 sigma.
- Difference images in 200-s sectors: see pixels figure (dip on target in every tested sector).
- Odd/even: see oddeven log in notes; model-fitted odd and even depths agree within ~2 sigma.
- Secondary eclipse at phase 0.5: consistent with zero (box test on stitched light curve).
- Gaia DR3 neighbors within 2 px: 1; brightest neighbor 11% of target light in G.
- Catalogs: not a TOI (Exoplanet Archive TOI table); not a known planet host; no ADS literature match for the TIC. CTOI: not in the ExoFOP CTOI list (4936 entries, 2026-10-26).

## Interpretation
Rp/Rs = 0.104 implies a 5.1 Rearth companion; nothing self-luminous is this small, so unless diluted by an unresolved brighter star the companion is planet-sized. Transit-derived stellar density agrees with the catalog within the b / a/Rs degeneracy.

## Figures
- s4_tic7586789_mcmc.png: phase-folded data and model
- s4_tic7586789.png: per-sector folds and joint fold
- s4_tic7586789_pixels.png: images and difference images, 200-s sectors
- s4_tic7586789_sector4_vetting.png: sector-4 vetting panel
- s4_tic7586789_corner.png: posterior
