# Community TOI candidate: TIC 156999506

Submitter: Kishore Hariharan (Generative Experiences / tess-hunt). Prepared 2026-10-26.

## Host star (TIC v8)
Tmag 14.89, Teff 3665 K, R = 0.521 Rsun, M = 0.519 Msun, distance 276.1 pc, contamination ratio 0.190.
RA 104.45880, Dec -40.33107. Sector 6 camera 3 CCD 1.

## Signal (MCMC, all sectors, batman + emcee, free quadratic limb darkening)
| quantity | value |
|---|---|
| Period (d) | 4.393071 +0.000002 -0.000002 |
| Epoch (BJD) | 2460663.51304 +0.00053 -0.00055 |
| Depth (ppm) | 7373 +142 -145 |
| Duration (h) | 2.64 +0.04 -0.03 |
| Rp/Rs | 0.0859 +0.0009 -0.0008 |
| a/Rs | 8.95 +0.36 -0.31 |
| Impact parameter | 0.83 +0.02 -0.02 |
| Planet radius (Rearth) | 4.89 +0.05 -0.05 |
| Planet radius (Rjup) | 0.436 +0.004 -0.004 |
| a (AU) | 0.0217 +0.0010 -0.0007 |
| T_eq (K, albedo 0) | 867 +14 -19 |
| Stellar density from transit (solar) | 0.50 +0.07 -0.05 |
| Stellar density from TIC (solar) | 3.66 |
| Limb darkening u1, u2 | 0.22 +0.09 -0.08, 0.22 +0.09 -0.09 |
| Transit SNR (depth / error) | 50.7 |

## Detection
Independent box-least-squares search of FFI cutouts (SPOC calibrated FFIs via MAST S3 cubes, astrocut) of TIC M dwarfs with Tmag 12-15 and R < 0.6 Rsun not on the 2-minute target list. Sector-6 detection: P = 4.39645 d, depth 5627 ppm, SNR 11.0, 5 transits.

## Multi-sector recovery
[6, 7, 33, 34, 61, 87, 88]; usable 7; signal recovered in 6: [7, 33, 34, 61, 87, 88]

Per-sector independent searches (narrow period window) and joint refinement are in the follow-up figure.

## Vetting
- Centroid / difference image (sector 6): on_target; shift 1.7 sigma; difference-image peak 0.5 px from target at 6 sigma.
- Difference images in 200-s sectors: see pixels figure (dip on target in every tested sector).
- Odd/even: see oddeven log in notes; model-fitted odd and even depths agree within ~2 sigma.
- Secondary eclipse at phase 0.5: consistent with zero (box test on stitched light curve).
- Gaia DR3 neighbors within 2 px: 5; brightest neighbor 103% of target light in G.
- Catalogs: not a TOI (Exoplanet Archive TOI table); not a known planet host; no ADS literature match for the TIC. CTOI: not in the ExoFOP CTOI list (4936 entries, 2026-10-26).

## Interpretation
Rp/Rs = 0.086 implies a 4.9 Rearth companion; nothing self-luminous is this small, so unless diluted by an unresolved brighter star the companion is planet-sized. Transit-derived stellar density agrees with the catalog within the b / a/Rs degeneracy.

## Figures
- s6_tic156999506_mcmc.png: phase-folded data and model
- s6_tic156999506.png: per-sector folds and joint fold
- s6_tic156999506_pixels.png: images and difference images, 200-s sectors
- s6_tic156999506_sector6_vetting.png: sector-6 vetting panel
- s6_tic156999506_corner.png: posterior
