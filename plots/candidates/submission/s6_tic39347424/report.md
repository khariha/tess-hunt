# Community TOI candidate: TIC 39347424

Submitter: Kishore Hariharan (Generative Experiences / tess-hunt). Prepared 2026-10-26.

## Host star (TIC v8)
Tmag 14.67, Teff 3858 K, R = 0.538 Rsun, M = 0.534 Msun, distance 270.1 pc, contamination ratio 0.030.
RA 88.91108, Dec -2.17499. Sector 6 camera 1 CCD 1.

## Signal (MCMC, all sectors, batman + emcee, free quadratic limb darkening)
| quantity | value |
|---|---|
| Period (d) | 0.926975 +0.000001 -0.000001 |
| Epoch (BJD) | 2460671.15542 +0.00251 -0.00238 |
| Depth (ppm) | 7738 +2477 -1793 |
| Duration (h) | 1.54 +0.06 -0.06 |
| Rp/Rs | 0.0890 +0.0137 -0.0133 |
| a/Rs | 2.32 +0.14 -0.10 |
| Impact parameter | 0.97 +0.02 -0.03 |
| Planet radius (Rearth) | 5.17 +0.77 -0.64 |
| Planet radius (Rjup) | 0.461 +0.069 -0.057 |
| a (AU) | 0.0059 +0.0003 -0.0003 |
| T_eq (K, albedo 0) | 1783 +40 -48 |
| Stellar density from transit (solar) | 0.20 +0.04 -0.02 |
| Stellar density from TIC (solar) | 3.43 |
| Limb darkening u1, u2 | 0.32 +0.10 -0.10, 0.32 +0.10 -0.10 |
| Transit SNR (depth / error) | 3.1 |

## Detection
Independent box-least-squares search of FFI cutouts (SPOC calibrated FFIs via MAST S3 cubes, astrocut) of TIC M dwarfs with Tmag 12-15 and R < 0.6 Rsun not on the 2-minute target list. Sector-6 detection: P = 0.92742 d, depth 5261 ppm, SNR 51.4, 23 transits.

## Multi-sector recovery
[6, 33, 87, 110, 111]; usable 3; signal recovered in 3: [6, 33, 87]

Per-sector independent searches (narrow period window) and joint refinement are in the follow-up figure.

## Vetting
- Centroid / difference image (sector 6): on_target; shift 0.2 sigma; difference-image peak 0.7 px from target at 6 sigma.
- Difference images in 200-s sectors: see pixels figure (dip on target in every tested sector).
- Odd/even: see oddeven log in notes; model-fitted odd and even depths agree within ~2 sigma.
- Secondary eclipse at phase 0.5: consistent with zero (box test on stitched light curve).
- Gaia DR3 neighbors within 2 px: 8; brightest neighbor 31% of target light in G.
- Catalogs: not a TOI (Exoplanet Archive TOI table); not a known planet host; no ADS literature match for the TIC. CTOI: not in the ExoFOP CTOI list (4936 entries, 2026-10-26).

## Interpretation
Rp/Rs = 0.089 implies a 5.2 Rearth companion; nothing self-luminous is this small, so unless diluted by an unresolved brighter star the companion is planet-sized. Transit-derived stellar density agrees with the catalog within the b / a/Rs degeneracy.

## Figures
- s6_tic39347424_mcmc.png: phase-folded data and model
- s6_tic39347424.png: per-sector folds and joint fold
- s6_tic39347424_pixels.png: images and difference images, 200-s sectors
- s6_tic39347424_sector6_vetting.png: sector-6 vetting panel
- s6_tic39347424_corner.png: posterior
