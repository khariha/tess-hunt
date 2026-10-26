# WITHDRAWN 2026-10-25: eclipsing binary, not a planet candidate

A deeper vetting pass found a 0.92% secondary eclipse at orbital phase 0.54 (45 sigma, all nine sectors, confirmed in QLP); the orbit is eccentric (e = 0.27) and the companion is a late-M dwarf or brown dwarf. See vetting_addendum.md. The brief below is kept for the record and must not be circulated.

# TIC 259719966: observer brief

Transiting Jupiter-radius companion candidate around a mid-M dwarf. Prepared
2026-10-25 by Kishore Hariharan (Generative Experiences / tess-hunt) from the
multi-sector fit in report.md. Research Note under review at RNAAS
(AAS82771); not a TOI or CTOI as of 2026-10-25. Code and full vetting:
https://github.com/khariha/tess-hunt

## Target

| | |
|---|---|
| TIC | 259719966 |
| Gaia DR3 | 6403128632917009792 |
| RA, Dec (J2000) | 322.67290, -63.61645 (21h30m41.5s, -63d36m59s) |
| Constellation | Indus |
| T / G | 13.80 / 14.98 |
| Host (TIC v8) | Teff 3443 K, R = 0.495 Rsun, M = 0.493 Msun, d = 147.1 pc |
| Gaia DR3 | RUWE 1.08, non_single_star = 0, single source within 3 arcsec |

## Ephemeris (BJD_TDB)

| | |
|---|---|
| P | 11.077207 +/- 0.000001 d |
| T0 | 2460877.82312 +/- 0.00020 |
| Duration (T14) | 2.37 +/- 0.03 h |
| Depth (undiluted sectors) | 38,300 +1,000 -1,200 ppm |
| Rp/Rs | 0.196 +/- 0.003 |
| a/Rs | 37.8 +/- 0.9 |
| b | 0.55 +/- 0.04 |
| Rp | 10.58 Rearth = 0.942 Rjup (+/- 0.015) |
| a | 0.087 AU |
| T_eq (A = 0) | 396 K |
| Transit-derived rho_star | 5.9 +/- 0.4 solar (TIC: 4.1) |

Prediction uncertainty stays under 1 minute through 2027 (P error times
cycle count).

Upcoming mid-transits, UTC (ingress and egress at -71 and +71 min):

| Epoch | Mid-transit |
|---|---|
| 42 | 2026-10-29 13:33 |
| 43 | 2026-11-09 15:24 |
| 44 | 2026-11-20 17:15 |
| 45 | 2026-12-01 19:06 |
| 46 | 2026-12-12 20:58 |
| 47 | 2026-12-23 22:49 |
| 48 | 2027-01-04 00:40 |
| 49 | 2027-01-15 02:31 |
| 50 | 2027-01-26 04:22 |
| 51 | 2027-02-06 06:14 |
| 52 | 2027-02-17 08:05 |
| 53 | 2027-02-28 09:56 |
| 54 | 2027-03-11 11:47 |
| 55 | 2027-03-22 13:38 |

Dec -64: circumpolar from Chile, South Africa and Australia; never rises
north of +26 latitude.

## Evidence summary

- Detected in Sector 1 (30-min FFIs) in an independent BLS search of 38,620
  Tmag 12-15 M dwarfs; recovered independently in all 9 sectors TESS has
  observed the star (1, 27, 67, 94, 95, 101-104), 20 transits, 2018-2026.
- Noise-normalized difference images in the 200-s sectors 94 and 104 peak on
  the target at 33 and 37 sigma. In the 30-min Sector 1 data the centroid
  shows a 6.4 sigma shift and the difference-image peak lies 0.7 px from the
  target at 10 sigma, consistent with the 21-arcsec pixels and the
  neighbour below rather than a displaced source; the 200-s sectors resolve
  it.
- Odd and even depths agree (-230 +/- 370 ppm). Secondary at phase 0.5:
  20 +/- 540 ppm, which disfavours a companion hot enough to be seen in the
  TESS band.
- Duration consistent with the TIC stellar density within the b / a/Rs
  degeneracy.

## Neighbour

Eight Gaia DR3 sources within 2 TESS pixels. The brightest, 2.1 px
(about 44 arcsec) from the target, carries 27% of the target's G flux and
enters the automatic TESS aperture in three sectors, where the depth
measures shallower. The sense of that effect (shallower when the neighbour
is included) independently places the dip on the target. Ground-based
apertures should exclude it; at 44 arcsec that is routine.

## What would settle the nature of the companion

Radius alone spans giant planet, brown dwarf and late-M star. Expected RV
semi-amplitude for a circular, edge-on orbit around 0.49 Msun at 11.08 d:

| Companion mass | K |
|---|---|
| 1 Mjup | ~150 m/s |
| 13 Mjup | ~1.9 km/s |
| 80 Mjup | ~12 km/s |

A single reconnaissance spectrum excludes the stellar case; two epochs near
quadrature (phases 0.25 and 0.75) with a red or near-infrared spectrograph on
a 4 m to 8 m telescope (ESPRESSO, NIRPS, HARPS, MAROON-X) separate planet
from brown dwarf. Host is faint for RV but not prohibitive: T = 13.8, and the
expected signal is large.

A ground-based transit on the target at the predicted time, with the
44-arcsec neighbour excluded, would close the remaining blend question. The
3.8% depth is reachable with a 0.4 m telescope.

## Caveats

Stellar parameters are TIC v8; the companion radius inherits their error.
No spectroscopic characterization of the host exists. Candidate status
until a mass is measured.

## Figures (this directory)

s1_tic259719966_mcmc.png (phase fold and model), s1_tic259719966.png
(per-sector folds), s1_tic259719966_pixels.png (difference images),
s1_tic259719966_corner.png (posterior), s1_tic259719966_sector1_vetting.png.
