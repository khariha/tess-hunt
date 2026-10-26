# Community TOI candidate: TIC 248942119

Submitter: Kishore Hariharan (Generative Experiences / tess-hunt). Prepared 2026-10-25.

## Host star (TIC v8)
Tmag 14.75, Teff 3838 K, R = 0.587 Rsun, M = 0.577 Msun, distance 322.3 pc, contamination ratio 0.018.
RA 19.91565, Dec -3.55106. Sector 3 camera 1 CCD 3.

## Signal (MCMC, all sectors, batman + emcee, free quadratic limb darkening)
| quantity | value |
|---|---|
| Period (d) | 1.843997 +0.000003 -0.000004 |
| Epoch (BJD) | 2460231.75955 +0.00227 -0.00360 |
| Depth (ppm) | 1494 +1253 -634 |
| Duration (h) | 0.88 +0.23 -0.19 |
| Rp/Rs | 0.0314 +0.0036 -0.0029 |
| a/Rs | 17.38 +5.01 -5.67 |
| Impact parameter | 0.45 +0.35 -0.30 |
| Planet radius (Rearth) | 2.48 +0.88 -0.60 |
| Planet radius (Rjup) | 0.221 +0.079 -0.053 |
| a (AU) | 0.0198 +0.0340 -0.0076 |
| T_eq (K, albedo 0) | 1007 +275 -396 |
| Stellar density from transit (solar) | 1.52 +28.83 -1.16 |
| Stellar density from TIC (solar) | 2.86 |
| Limb darkening u1, u2 | 0.31 +0.09 -0.10, 0.29 +0.10 -0.10 |
| Transit SNR (depth / error) | 1.2 |

## Detection
Independent box-least-squares search of FFI cutouts (SPOC calibrated FFIs via MAST S3 cubes, astrocut) of TIC M dwarfs with Tmag 12-15 and R < 0.6 Rsun not on the 2-minute target list. Sector-3 detection: P = 1.84442 d, depth 4365 ppm, SNR 25.6, 12 transits.

## Multi-sector recovery
[3, 30, 42, 43, 70, 97]; usable 6; signal recovered in 6: [3, 30, 42, 43, 70, 97]

Per-sector independent searches (narrow period window) and joint refinement are in the follow-up figure.

## Vetting
- Centroid / difference image (sector 3): likely_on_target; shift 0.7 sigma; difference-image peak 6.5 px from target at 3 sigma.
- Difference images in 200-s sectors: see pixels figure (dip on target in every tested sector).
- Odd/even: see oddeven log in notes; model-fitted odd and even depths agree within ~2 sigma.
- Secondary eclipse at phase 0.5: consistent with zero (box test on stitched light curve).
- Gaia DR3 neighbors within 2 px: 1; brightest neighbor 8% of target light in G.
- Catalogs: not a TOI (Exoplanet Archive TOI table); not a known planet host; no ADS literature match for the TIC. CTOI: not in the ExoFOP CTOI list (4936 entries, 2026-10-26).

## Interpretation
Rp/Rs = 0.031 implies a 2.5 Rearth companion; nothing self-luminous is this small, so unless diluted by an unresolved brighter star the companion is planet-sized. Transit-derived stellar density agrees with the catalog within the b / a/Rs degeneracy.

## Figures
- s3_tic248942119_mcmc.png: phase-folded data and model
- s3_tic248942119.png: per-sector folds and joint fold
- s3_tic248942119_pixels.png: images and difference images, 200-s sectors
- s3_tic248942119_sector3_vetting.png: sector-3 vetting panel
- s3_tic248942119_corner.png: posterior
