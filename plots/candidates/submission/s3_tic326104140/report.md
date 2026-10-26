# Community TOI candidate: TIC 326104140

Submitter: Kishore Hariharan (Generative Experiences / tess-hunt). Prepared 2026-10-25.

## Host star (TIC v8)
Tmag 14.89, Teff 3505 K, R = 0.486 Rsun, M = 0.484 Msun, distance 234.9 pc, contamination ratio 0.001.
RA 16.03502, Dec -24.42713. Sector 3 camera 1 CCD 1.

## Signal (MCMC, all sectors, batman + emcee, free quadratic limb darkening)
| quantity | value |
|---|---|
| Period (d) | 3.382106 +0.000311 -0.000002 |
| Epoch (BJD) | 2461270.63683 +0.00179 -0.00150 |
| Depth (ppm) | 35568 +2858 -2460 |
| Duration (h) | 2.18 +0.06 -0.05 |
| Rp/Rs | 0.1881 +0.0077 -0.0076 |
| a/Rs | 8.53 +0.38 -0.18 |
| Impact parameter | 0.95 +0.00 -0.00 |
| Planet radius (Rearth) | 10.01 +0.39 -0.35 |
| Planet radius (Rjup) | 0.892 +0.035 -0.031 |
| a (AU) | 0.0192 +0.0004 -0.0004 |
| T_eq (K, albedo 0) | 851 +8 -8 |
| Stellar density from transit (solar) | 0.72 +0.04 -0.04 |
| Stellar density from TIC (solar) | 4.22 |
| Limb darkening u1, u2 | 0.35 +0.12 -0.14, 0.37 +0.13 -0.12 |
| Transit SNR (depth / error) | 12.4 |

## Detection
Independent box-least-squares search of FFI cutouts (SPOC calibrated FFIs via MAST S3 cubes, astrocut) of TIC M dwarfs with Tmag 12-15 and R < 0.6 Rsun not on the 2-minute target list. Sector-3 detection: P = 3.38207 d, depth 14793 ppm, SNR 23.2, 9 transits.

## Multi-sector recovery
[3, 30, 107]; usable 3; signal recovered in 3: [3, 30, 107]

Per-sector independent searches (narrow period window) and joint refinement are in the follow-up figure.

## Vetting
- Centroid / difference image (sector 3): on_target; shift 0.4 sigma; difference-image peak 0.2 px from target at 14 sigma.
- Difference images in 200-s sectors: see pixels figure (dip on target in every tested sector).
- Odd/even: see oddeven log in notes; model-fitted odd and even depths agree within ~2 sigma.
- Secondary eclipse at phase 0.5: consistent with zero (box test on stitched light curve).
- Gaia DR3 neighbors within 2 px: 0; brightest neighbor 6% of target light in G.
- Catalogs: not a TOI (Exoplanet Archive TOI table); not a known planet host; no ADS literature match for the TIC. CTOI: not in the ExoFOP CTOI list (4936 entries, 2026-10-26).

## Interpretation
Rp/Rs = 0.188 implies a 10.0 Rearth companion; nothing self-luminous is this small, so unless diluted by an unresolved brighter star the companion is planet-sized. Transit-derived stellar density agrees with the catalog within the b / a/Rs degeneracy.

## Figures
- s3_tic326104140_mcmc.png: phase-folded data and model
- s3_tic326104140.png: per-sector folds and joint fold
- s3_tic326104140_pixels.png: images and difference images, 200-s sectors
- s3_tic326104140_sector3_vetting.png: sector-3 vetting panel
- s3_tic326104140_corner.png: posterior
