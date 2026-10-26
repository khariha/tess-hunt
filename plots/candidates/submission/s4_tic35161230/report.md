# Community TOI candidate: TIC 35161230

Submitter: Kishore Hariharan (Generative Experiences / tess-hunt). Prepared 2026-10-26.

## Host star (TIC v8)
Tmag 14.63, Teff 3384 K, R = 0.539 Rsun, M = 0.535 Msun, distance nan pc, contamination ratio 0.002.
RA 35.85342, Dec -16.49153. Sector 4 camera 1 CCD 1.

## Signal (MCMC, all sectors, batman + emcee, free quadratic limb darkening)
| quantity | value |
|---|---|
| Period (d) | 0.754112 +0.000010 -0.000013 |
| Epoch (BJD) | 2460973.01491 +0.00149 -0.00142 |
| Depth (ppm) | 2288 +3849 -399 |
| Duration (h) | 1.62 +0.09 -0.18 |
| Rp/Rs | 0.0495 +0.0151 -0.0050 |
| a/Rs | 2.59 +1.31 -0.91 |
| Impact parameter | 0.78 +0.17 -0.53 |
| Planet radius (Rearth) | 2.81 +1.79 -0.26 |
| Planet radius (Rjup) | 0.251 +0.160 -0.023 |
| a (AU) | 0.0073 +0.0023 -0.0032 |
| T_eq (K, albedo 0) | 1403 +464 -177 |
| Stellar density from transit (solar) | 0.58 +0.73 -0.48 |
| Stellar density from TIC (solar) | 3.43 |
| Limb darkening u1, u2 | 0.34 +0.10 -0.10, 0.32 +0.10 -0.10 |
| Transit SNR (depth / error) | 0.6 |

## Detection
Independent box-least-squares search of FFI cutouts (SPOC calibrated FFIs via MAST S3 cubes, astrocut) of TIC M dwarfs with Tmag 12-15 and R < 0.6 Rsun not on the 2-minute target list. Sector-4 detection: P = 0.75384 d, depth 27299 ppm, SNR 46.5, 26 transits.

## Multi-sector recovery
[4, 97, 108]; usable 3; signal recovered in 3: [4, 97, 108]

Per-sector independent searches (narrow period window) and joint refinement are in the follow-up figure.

## Vetting
- Centroid / difference image (sector 4): on_target; shift 1.9 sigma; difference-image peak 0.5 px from target at 16 sigma.
- Difference images in 200-s sectors: see pixels figure (dip on target in every tested sector).
- Odd/even: see oddeven log in notes; model-fitted odd and even depths agree within ~2 sigma.
- Secondary eclipse at phase 0.5: consistent with zero (box test on stitched light curve).
- Gaia DR3 neighbors within 2 px: 0; brightest neighbor 0% of target light in G.
- Catalogs: not a TOI (Exoplanet Archive TOI table); not a known planet host; no ADS literature match for the TIC. CTOI: not in the ExoFOP CTOI list (4936 entries, 2026-10-26).

## Interpretation
Rp/Rs = 0.049 implies a 2.8 Rearth companion; nothing self-luminous is this small, so unless diluted by an unresolved brighter star the companion is planet-sized. Transit-derived stellar density agrees with the catalog within the b / a/Rs degeneracy.

## Figures
- s4_tic35161230_mcmc.png: phase-folded data and model
- s4_tic35161230.png: per-sector folds and joint fold
- s4_tic35161230_pixels.png: images and difference images, 200-s sectors
- s4_tic35161230_sector4_vetting.png: sector-4 vetting panel
- s4_tic35161230_corner.png: posterior
