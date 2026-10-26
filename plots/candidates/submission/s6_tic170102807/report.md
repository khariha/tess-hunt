# Community TOI candidate: TIC 170102807

Submitter: Kishore Hariharan (Generative Experiences / tess-hunt). Prepared 2026-10-25.

## Host star (TIC v8)
Tmag 13.47, Teff 3535 K, R = 0.511 Rsun, M = 0.509 Msun, distance 135.0 pc, contamination ratio 0.010.
RA 101.50882, Dec -47.80176. Sector 6 camera 3 CCD 4.

## Signal (MCMC, all sectors, batman + emcee, free quadratic limb darkening)
| quantity | value |
|---|---|
| Period (d) | 0.952986 +0.000002 -0.000002 |
| Epoch (BJD) | 2460665.29534 +0.00344 -0.00528 |
| Depth (ppm) | 2365 +2557 -531 |
| Duration (h) | 0.94 +0.12 -0.08 |
| Rp/Rs | 0.0589 +0.0184 -0.0137 |
| a/Rs | 3.42 +4.91 -0.46 |
| Impact parameter | 0.94 +0.04 -0.58 |
| Planet radius (Rearth) | 2.71 +1.20 -0.32 |
| Planet radius (Rjup) | 0.242 +0.107 -0.029 |
| a (AU) | 0.0151 +0.0047 -0.0077 |
| T_eq (K, albedo 0) | 993 +424 -126 |
| Stellar density from transit (solar) | 3.76 +4.72 -3.32 |
| Stellar density from TIC (solar) | 3.82 |
| Limb darkening u1, u2 | 0.32 +0.10 -0.10, 0.32 +0.10 -0.10 |
| Transit SNR (depth / error) | 0.9 |

## Detection
Independent box-least-squares search of FFI cutouts (SPOC calibrated FFIs via MAST S3 cubes, astrocut) of TIC M dwarfs with Tmag 12-15 and R < 0.6 Rsun not on the 2-minute target list. Sector-6 detection: P = 0.95282 d, depth 7528 ppm, SNR 37.3, 21 transits.

## Multi-sector recovery
[6, 7, 8, 33, 34, 61, 87, 88]; usable 8; signal recovered in 8: [6, 7, 8, 33, 34, 61, 87, 88]

Per-sector independent searches (narrow period window) and joint refinement are in the follow-up figure.

## Vetting
- Centroid / difference image (sector 6): on_target; shift 0.9 sigma; difference-image peak 0.5 px from target at 7 sigma.
- Difference images in 200-s sectors: see pixels figure (dip on target in every tested sector).
- Odd/even: see oddeven log in notes; model-fitted odd and even depths agree within ~2 sigma.
- Secondary eclipse at phase 0.5: consistent with zero (box test on stitched light curve).
- Gaia DR3 neighbors within 2 px: 4; brightest neighbor 4% of target light in G.
- Catalogs: not a TOI (Exoplanet Archive TOI table); not a known planet host; no ADS literature match for the TIC. CTOI: not in the ExoFOP CTOI list (4936 entries, 2026-10-26).

## Interpretation
Rp/Rs = 0.059 implies a 2.7 Rearth companion; nothing self-luminous is this small, so unless diluted by an unresolved brighter star the companion is planet-sized. Transit-derived stellar density agrees with the catalog within the b / a/Rs degeneracy.

## Figures
- s6_tic170102807_mcmc.png: phase-folded data and model
- s6_tic170102807.png: per-sector folds and joint fold
- s6_tic170102807_pixels.png: images and difference images, 200-s sectors
- s6_tic170102807_sector6_vetting.png: sector-6 vetting panel
- s6_tic170102807_corner.png: posterior
