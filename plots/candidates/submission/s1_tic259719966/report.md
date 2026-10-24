# Community TOI candidate: TIC 259719966

Submitter: Kishore Hariharan (Generative Experiences / tess-hunt). Prepared 2026-10-24.

## Host star (TIC v8)
Tmag 13.80, Teff 3443 K, R = 0.495 Rsun, M = 0.493 Msun, distance 147.1 pc, contamination ratio 0.007.
RA 322.67290, Dec -63.61645. Sector 1 camera 2 CCD 1.

## Signal (MCMC, all sectors, batman + emcee, free quadratic limb darkening)
| quantity | value |
|---|---|
| Period (d) | 11.077207 +0.000001 -0.000001 |
| Epoch (BJD) | 2460877.82312 +0.00020 -0.00019 |
| Depth (ppm) | 38326 +1037 -1231 |
| Duration (h) | 2.37 +0.02 -0.03 |
| Rp/Rs | 0.1960 +0.0026 -0.0031 |
| a/Rs | 37.84 +0.90 -0.81 |
| Impact parameter | 0.55 +0.04 -0.04 |
| Planet radius (Rearth) | 10.58 +0.14 -0.17 |
| Planet radius (Rjup) | 0.942 +0.013 -0.015 |
| a (AU) | 0.0869 +0.0023 -0.0018 |
| T_eq (K, albedo 0) | 396 +4 -5 |
| Stellar density from transit (solar) | 5.90 +0.48 -0.36 |
| Stellar density from TIC (solar) | 4.08 |
| Limb darkening u1, u2 | 0.36 +0.08 -0.08, 0.33 +0.09 -0.09 |
| Transit SNR (depth / error) | 31.1 |

## Detection
Independent box-least-squares search of FFI cutouts (SPOC calibrated FFIs via MAST S3 cubes, astrocut) of TIC M dwarfs with Tmag 12-15 and R < 0.6 Rsun not on the 2-minute target list. Sector-1 detection: P = 11.07286 d, depth 36495 ppm, SNR 10.9, 3 transits.

## Multi-sector recovery
[1, 27, 67, 94, 95, 101, 102, 103, 104]; usable 9; signal recovered in 9: [1, 27, 67, 94, 95, 101, 102, 103, 104]

Per-sector independent searches (narrow period window) and joint refinement are in the follow-up figure.

## Vetting
- Centroid / difference image (sector 1): anomaly; shift 6.4 sigma; difference-image peak 0.7 px from target at 10 sigma.
- Difference images in 200-s sectors: see pixels figure (dip on target in every tested sector).
- Odd/even: see oddeven log in notes; model-fitted odd and even depths agree within ~2 sigma.
- Secondary eclipse at phase 0.5: consistent with zero (box test on stitched light curve).
- Gaia DR3 neighbors within 2 px: 8; brightest neighbor 27% of target light in G.
- Catalogs: not a TOI (Exoplanet Archive TOI table); not a known planet host; no ADS literature match for the TIC. CTOI: not in the ExoFOP CTOI list (4936 entries, 2026-10-24).

## Interpretation
Rp/Rs = 0.196 implies a 10.6 Rearth companion; nothing self-luminous is this small, so unless diluted by an unresolved brighter star the companion is planet-sized. Transit-derived stellar density agrees with the catalog within the b / a/Rs degeneracy.

## Figures
- s1_tic259719966_mcmc.png: phase-folded data and model
- s1_tic259719966.png: per-sector folds and joint fold
- s1_tic259719966_pixels.png: images and difference images, 200-s sectors
- s1_tic259719966_sector1_vetting.png: sector-1 vetting panel
- s1_tic259719966_corner.png: posterior
