# Community TOI candidate: TIC 46432934

Submitter: Kishore Hariharan (Generative Experiences / tess-hunt). Prepared 2026-10-26.

## Host star (TIC v8)
Tmag 13.31, Teff 3418 K, R = 0.413 Rsun, M = 0.406 Msun, distance 90.7 pc, contamination ratio 0.045.
RA 83.87849, Dec -14.60354. Sector 6 camera 2 CCD 4.

## Signal (MCMC, all sectors, batman + emcee, free quadratic limb darkening)
| quantity | value |
|---|---|
| Period (d) | 1.440445 +0.000000 -0.000000 |
| Epoch (BJD) | 2461045.58677 +0.00004 -0.00004 |
| Depth (ppm) | 24132 +2456 -1170 |
| Duration (h) | 1.20 +0.01 -0.01 |
| Rp/Rs | 0.1548 +0.0046 -0.0031 |
| a/Rs | 6.83 +0.09 -0.11 |
| Impact parameter | 0.89 +0.02 -0.01 |
| Planet radius (Rearth) | 7.01 +0.35 -0.17 |
| Planet radius (Rjup) | 0.624 +0.031 -0.015 |
| a (AU) | 0.0131 +0.0002 -0.0003 |
| T_eq (K, albedo 0) | 926 +10 -7 |
| Stellar density from transit (solar) | 2.04 +0.10 -0.12 |
| Stellar density from TIC (solar) | 5.77 |
| Limb darkening u1, u2 | 0.27 +0.09 -0.10, 0.28 +0.09 -0.10 |
| Transit SNR (depth / error) | 9.8 |

## Detection
Independent box-least-squares search of FFI cutouts (SPOC calibrated FFIs via MAST S3 cubes, astrocut) of TIC M dwarfs with Tmag 12-15 and R < 0.6 Rsun not on the 2-minute target list. Sector-6 detection: P = 1.44054 d, depth 14017 ppm, SNR 22.8, 13 transits.

## Multi-sector recovery
[6, 32, 98, 107, 108, 109, 110]; usable 6; signal recovered in 5: [6, 32, 98, 107, 108]

Per-sector independent searches (narrow period window) and joint refinement are in the follow-up figure.

## Vetting
- Centroid / difference image (sector 6): on_target; shift 12.6 sigma; difference-image peak 0.6 px from target at 7 sigma.
- Difference images in 200-s sectors: see pixels figure (dip on target in every tested sector).
- Odd/even: see oddeven log in notes; model-fitted odd and even depths agree within ~2 sigma.
- Secondary eclipse at phase 0.5: consistent with zero (box test on stitched light curve).
- Gaia DR3 neighbors within 2 px: 3; brightest neighbor 272% of target light in G.
- Catalogs: not a TOI (Exoplanet Archive TOI table); not a known planet host; no ADS literature match for the TIC. CTOI: not in the ExoFOP CTOI list (4936 entries, 2026-10-26).

## Interpretation
Rp/Rs = 0.155 implies a 7.0 Rearth companion; nothing self-luminous is this small, so unless diluted by an unresolved brighter star the companion is planet-sized. Transit-derived stellar density agrees with the catalog within the b / a/Rs degeneracy.

## Figures
- s6_tic46432934_mcmc.png: phase-folded data and model
- s6_tic46432934.png: per-sector folds and joint fold
- s6_tic46432934_pixels.png: images and difference images, 200-s sectors
- s6_tic46432934_sector6_vetting.png: sector-6 vetting panel
- s6_tic46432934_corner.png: posterior
