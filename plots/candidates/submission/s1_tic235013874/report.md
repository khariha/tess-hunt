# Community TOI candidate: TIC 235013874

Submitter: Kishore Hariharan (Generative Experiences / tess-hunt). Prepared 2026-10-24.

## Host star (TIC v8)
Tmag 14.55, Teff 3259 K, R = 0.447 Rsun, M = 0.443 Msun, distance 166.5 pc, contamination ratio 0.001.
RA 353.39139, Dec -63.74506. Sector 1 camera 2 CCD 2.

## Signal (MCMC, all sectors, batman + emcee, free quadratic limb darkening)
| quantity | value |
|---|---|
| Period (d) | 5.750919 +0.000002 -0.000002 |
| Epoch (BJD) | 2460904.21833 +0.00031 -0.00032 |
| Depth (ppm) | 21265 +7282 -2585 |
| Duration (h) | 1.46 +0.03 -0.03 |
| Rp/Rs | 0.1450 +0.0267 -0.0089 |
| a/Rs | 21.92 +1.20 -1.21 |
| Impact parameter | 0.89 +0.06 -0.03 |
| Planet radius (Rearth) | 7.12 +1.13 -0.45 |
| Planet radius (Rjup) | 0.634 +0.101 -0.040 |
| a (AU) | 0.0454 +0.0021 -0.0022 |
| T_eq (K, albedo 0) | 493 +12 -11 |
| Stellar density from transit (solar) | 4.23 +0.61 -0.59 |
| Stellar density from TIC (solar) | 4.97 |
| Limb darkening u1, u2 | 0.28 +0.10 -0.10, 0.28 +0.10 -0.10 |
| Transit SNR (depth / error) | 2.9 |

## Detection
Independent box-least-squares search of FFI cutouts (SPOC calibrated FFIs via MAST S3 cubes, astrocut) of TIC M dwarfs with Tmag 12-15 and R < 0.6 Rsun not on the 2-minute target list. Sector-1 detection: P = 5.75036 d, depth 10110 ppm, SNR 11.5, 4 transits.

## Multi-sector recovery
[1, 28, 68, 95, 102, 103, 104]; usable 7; signal recovered in 7: [1, 28, 68, 95, 102, 103, 104]

Per-sector independent searches (narrow period window) and joint refinement are in the follow-up figure.

## Vetting
- Centroid / difference image (sector 1): on_target; shift 0.8 sigma; difference-image peak 0.2 px from target at 7 sigma.
- Difference images in 200-s sectors: see pixels figure (dip on target in every tested sector).
- Odd/even: see oddeven log in notes; model-fitted odd and even depths agree within ~2 sigma.
- Secondary eclipse at phase 0.5: consistent with zero (box test on stitched light curve).
- Gaia DR3 neighbors within 2 px: 1; brightest neighbor 1% of target light in G.
- Catalogs: not a TOI (Exoplanet Archive TOI table); not a known planet host; no ADS literature match for the TIC. CTOI list: check pending ExoFOP availability.

## Interpretation
Rp/Rs = 0.145 implies a 7.1 Rearth companion; nothing self-luminous is this small, so unless diluted by an unresolved brighter star the companion is planet-sized. Transit-derived stellar density agrees with the catalog within the b / a/Rs degeneracy.

## Figures
- s1_tic235013874_mcmc.png: phase-folded data and model
- s1_tic235013874.png: per-sector folds and joint fold
- s1_tic235013874_pixels.png: images and difference images, 200-s sectors
- s1_tic235013874_sector1_vetting.png: sector-1 vetting panel
- s1_tic235013874_corner.png: posterior
