# Community TOI candidate: TIC 114949540

Submitter: Kishore Hariharan (Generative Experiences / tess-hunt). Prepared 2026-10-24.

## Host star (TIC v8)
Tmag 14.78, Teff 3320 K, R = 0.422 Rsun, M = 0.416 Msun, distance 168.4 pc, contamination ratio 0.013.
RA 315.13674, Dec -37.41139. Sector 1 camera 1 CCD 1.

## Signal (MCMC, all sectors, batman + emcee, free quadratic limb darkening)
| quantity | value |
|---|---|
| Period (d) | 2.739176 +0.000001 -0.000001 |
| Epoch (BJD) | 2461212.95052 +0.00074 -0.00074 |
| Depth (ppm) | 8927 +713 -628 |
| Duration (h) | 1.55 +0.07 -0.07 |
| Rp/Rs | 0.0930 +0.0045 -0.0035 |
| a/Rs | 13.12 +1.97 -2.31 |
| Impact parameter | 0.53 +0.19 -0.34 |
| Planet radius (Rearth) | 4.35 +0.17 -0.16 |
| Planet radius (Rjup) | 0.388 +0.015 -0.014 |
| a (AU) | 0.0242 +0.0041 -0.0034 |
| T_eq (K, albedo 0) | 669 +52 -50 |
| Stellar density from transit (solar) | 3.35 +1.98 -1.21 |
| Stellar density from TIC (solar) | 5.55 |
| Limb darkening u1, u2 | 0.35 +0.09 -0.10, 0.32 +0.10 -0.10 |
| Transit SNR (depth / error) | 12.5 |

## Detection
Independent box-least-squares search of FFI cutouts (SPOC calibrated FFIs via MAST S3 cubes, astrocut) of TIC M dwarfs with Tmag 12-15 and R < 0.6 Rsun not on the 2-minute target list. Sector-1 detection: P = 2.73856 d, depth 7793 ppm, SNR 22.3, 9 transits.

## Multi-sector recovery
[1, 27, 105]; usable 3; signal recovered in 3: [1, 27, 105]

Per-sector independent searches (narrow period window) and joint refinement are in the follow-up figure.

## Vetting
- Centroid / difference image (sector 1): on_target; shift 1.6 sigma; difference-image peak 0.5 px from target at 6 sigma.
- Difference images in 200-s sectors: see pixels figure (dip on target in every tested sector).
- Odd/even: see oddeven log in notes; model-fitted odd and even depths agree within ~2 sigma.
- Secondary eclipse at phase 0.5: consistent with zero (box test on stitched light curve).
- Gaia DR3 neighbors within 2 px: 4; brightest neighbor 12% of target light in G.
- Catalogs: not a TOI (Exoplanet Archive TOI table); not a known planet host; no ADS literature match for the TIC. CTOI: not in the ExoFOP CTOI list (4936 entries, 2026-10-24).

## Interpretation
Rp/Rs = 0.093 implies a 4.3 Rearth companion; nothing self-luminous is this small, so unless diluted by an unresolved brighter star the companion is planet-sized. Transit-derived stellar density agrees with the catalog within the b / a/Rs degeneracy.

## Figures
- s1_tic114949540_mcmc.png: phase-folded data and model
- s1_tic114949540.png: per-sector folds and joint fold
- s1_tic114949540_pixels.png: images and difference images, 200-s sectors
- s1_tic114949540_sector1_vetting.png: sector-1 vetting panel
- s1_tic114949540_corner.png: posterior
