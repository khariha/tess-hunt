# Community TOI candidate: TIC 389723320

Submitter: Kishore Hariharan (Generative Experiences / tess-hunt). Prepared 2026-10-24.

## Host star (TIC v8)
Tmag 14.52, Teff 3271 K, R = 0.373 Rsun, M = 0.361 Msun, distance 126.2 pc, contamination ratio 0.008.
RA 329.29429, Dec -44.60839. Sector 1 camera 1 CCD 1.

## Signal (MCMC, all sectors, batman + emcee, free quadratic limb darkening)
| quantity | value |
|---|---|
| Period (d) | 3.927167 +0.000001 -0.000001 |
| Epoch (BJD) | 2460899.38949 +0.00050 -0.00051 |
| Depth (ppm) | 12849 +596 -690 |
| Duration (h) | 1.50 +0.06 -0.05 |
| Rp/Rs | 0.1137 +0.0027 -0.0035 |
| a/Rs | 16.64 +1.88 -1.32 |
| Impact parameter | 0.74 +0.05 -0.09 |
| Planet radius (Rearth) | 4.62 +0.11 -0.13 |
| Planet radius (Rjup) | 0.412 +0.009 -0.011 |
| a (AU) | 0.0291 +0.0033 -0.0025 |
| T_eq (K, albedo 0) | 565 +26 -29 |
| Stellar density from transit (solar) | 4.11 +1.55 -0.98 |
| Stellar density from TIC (solar) | 6.95 |
| Limb darkening u1, u2 | 0.31 +0.10 -0.10, 0.31 +0.10 -0.10 |
| Transit SNR (depth / error) | 18.6 |

## Detection
Independent box-least-squares search of FFI cutouts (SPOC calibrated FFIs via MAST S3 cubes, astrocut) of TIC M dwarfs with Tmag 12-15 and R < 0.6 Rsun not on the 2-minute target list. Sector-1 detection: P = 3.92736 d, depth 11269 ppm, SNR 15.9, 6 transits.

## Multi-sector recovery
[1, 28, 68, 95, 105, 106]; usable 6; signal recovered in 6: [1, 28, 68, 95, 105, 106]

Per-sector independent searches (narrow period window) and joint refinement are in the follow-up figure.

## Vetting
- Centroid / difference image (sector 1): on_target; shift 0.2 sigma; difference-image peak 0.5 px from target at 5 sigma.
- Difference images in 200-s sectors: see pixels figure (dip on target in every tested sector).
- Odd/even: see oddeven log in notes; model-fitted odd and even depths agree within ~2 sigma.
- Secondary eclipse at phase 0.5: consistent with zero (box test on stitched light curve).
- Gaia DR3 neighbors within 2 px: 2; brightest neighbor 23% of target light in G.
- Catalogs: not a TOI (Exoplanet Archive TOI table); not a known planet host; no ADS literature match for the TIC. CTOI list: check pending ExoFOP availability.

## Interpretation
Rp/Rs = 0.114 implies a 4.6 Rearth companion; nothing self-luminous is this small, so unless diluted by an unresolved brighter star the companion is planet-sized. Transit-derived stellar density agrees with the catalog within the b / a/Rs degeneracy.

## Figures
- s1_tic389723320_mcmc.png: phase-folded data and model
- s1_tic389723320.png: per-sector folds and joint fold
- s1_tic389723320_pixels.png: images and difference images, 200-s sectors
- s1_tic389723320_sector1_vetting.png: sector-1 vetting panel
- s1_tic389723320_corner.png: posterior
