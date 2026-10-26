# Community TOI candidate: TIC 166760940

Submitter: Kishore Hariharan (Generative Experiences / tess-hunt). Prepared 2026-10-25.

## Host star (TIC v8)
Tmag 14.74, Teff 3179 K, R = 0.264 Rsun, M = 0.236 Msun, distance 87.6 pc, contamination ratio 0.013.
RA 37.66595, Dec -49.32110. Sector 3 camera 3 CCD 2.

## Signal (MCMC, all sectors, batman + emcee, free quadratic limb darkening)
| quantity | value |
|---|---|
| Period (d) | 3.102852 +0.000002 -0.000002 |
| Epoch (BJD) | 2460979.03908 +0.00065 -0.00071 |
| Depth (ppm) | 3489 +331 -251 |
| Duration (h) | 1.00 +0.05 -0.04 |
| Rp/Rs | 0.0591 +0.0029 -0.0021 |
| a/Rs | 23.18 +2.13 -5.05 |
| Impact parameter | 0.41 +0.30 -0.28 |
| Planet radius (Rearth) | 1.70 +0.08 -0.06 |
| Planet radius (Rjup) | 0.152 +0.007 -0.006 |
| a (AU) | 0.0283 +0.0026 -0.0055 |
| T_eq (K, albedo 0) | 468 +53 -20 |
| Stellar density from transit (solar) | 17.17 +5.16 -8.17 |
| Stellar density from TIC (solar) | 12.88 |
| Limb darkening u1, u2 | 0.31 +0.10 -0.10, 0.30 +0.10 -0.09 |
| Transit SNR (depth / error) | 10.5 |

## Detection
Independent box-least-squares search of FFI cutouts (SPOC calibrated FFIs via MAST S3 cubes, astrocut) of TIC M dwarfs with Tmag 12-15 and R < 0.6 Rsun not on the 2-minute target list. Sector-3 detection: P = 3.10265 d, depth 4995 ppm, SNR 12.6, 7 transits.

## Multi-sector recovery
[3, 29, 30, 69, 97, 104, 105, 106]; usable 8; signal recovered in 8: [3, 29, 30, 69, 97, 104, 105, 106]

Per-sector independent searches (narrow period window) and joint refinement are in the follow-up figure.

## Vetting
- Centroid / difference image (sector 3): likely_on_target; shift 1.1 sigma; difference-image peak 9.3 px from target at 3 sigma.
- Difference images in 200-s sectors: see pixels figure (dip on target in every tested sector).
- Odd/even: see oddeven log in notes; model-fitted odd and even depths agree within ~2 sigma.
- Secondary eclipse at phase 0.5: consistent with zero (box test on stitched light curve).
- Gaia DR3 neighbors within 2 px: 2; brightest neighbor 10% of target light in G.
- Catalogs: not a TOI (Exoplanet Archive TOI table); not a known planet host; no ADS literature match for the TIC. CTOI: not in the ExoFOP CTOI list (4936 entries, 2026-10-26).

## Interpretation
Rp/Rs = 0.059 implies a 1.7 Rearth companion; nothing self-luminous is this small, so unless diluted by an unresolved brighter star the companion is planet-sized. Transit-derived stellar density agrees with the catalog within the b / a/Rs degeneracy.

## Figures
- s3_tic166760940_mcmc.png: phase-folded data and model
- s3_tic166760940.png: per-sector folds and joint fold
- s3_tic166760940_pixels.png: images and difference images, 200-s sectors
- s3_tic166760940_sector3_vetting.png: sector-3 vetting panel
- s3_tic166760940_corner.png: posterior
