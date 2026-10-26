# tess-hunt

Transit photometry search over TESS data. Goal: an agent that finds unreported planet candidates in full-frame images (FFIs), vets them, and submits worthwhile ones as Community TOIs on ExoFOP.

Kishore is new to astronomy. Explain concepts plainly, ramp up gradually, draw a diagram when it helps.

## Status (24 Oct 2026)
- Sector 1 swept (38,620 red dwarfs, Tmag 12-15, R<0.6): 22 flagged, vetted in data/vetting_s1.csv.
- Submission batch, packages in plots/candidates/submission/ (tracked in git):
  1. TIC 389723320  P=3.927167 d  4.6 Re sub-Neptune, 6 sectors, cleanest
  2. TIC 114949540  P=2.739176 d  4.35 Re hot Neptune, 3 sectors
  3. ~~TIC 259719966~~ WITHDRAWN 2026-10-25: eccentric M+M eclipsing binary (secondary 0.9% at phase 0.54); see plots/candidates/submission/s1_tic259719966/vetting_addendum.md
  Held back: TIC 235013874 (grazing, radius lower limit ~6 Re).
- ExoFOP policy (Aug 2026): community candidates must be published first (refereed journal, or an RNAAS citing published methods). CTOI duplicate check done 2026-10-24: none of the four listed. Kishore has an ExoFOP account.
- Research Note submitted to RNAAS 2026-10-25, manuscript AAS82771 (docs/rnaas/). Revision (two candidates, TIC 259719966 described as an eccentric EB) sent to AAS staff via manuscript correspondence 2026-10-25 evening with rnaas.tex, figure1.pdf, rnaas.pdf attached; awaiting their reply. Next: on acceptance, post to arXiv astro-ph.EP, then file the ExoFOP Published Candidate Upload Request with the note's URL and enter the two candidates from plots/candidates/submission/s1_tic{389723320,114949540}/ctoi_fields.json.
- Public mentions of three candidates (Kanodia email, Sentinel pitch, LinkedIn) are corrected only if someone follows up; genexperiences.com/tess-hunt already shows two.
- Public repo: https://github.com/khariha/tess-hunt
- Running: `modal run --detach scripts/modal_sweep.py --sectors "2 3 4 5 6"` (select + sweep per sector, all in cloud).
  Check: `--sector N --mode status`; fetch: `--sector N --mode pull`; then `scripts/vet_candidates.py --sector N`,
  `scripts/followup.py`, `scripts/mcmc_fit.py`, `scripts/make_ctoi_package.py` for survivors.
- Also check our candidates against published faint-star FFI candidate tables (DIAmante, TESS faint-star search) before submitting.

## Pipeline (scripts, in order)
select_targets.py -> modal_sweep.py (or sweep.py) -> vet_candidates.py -> followup.py -> oddeven_fit.py -> mcmc_fit.py -> make_ctoi_package.py; residual_search.py looks for a second planet.

## Layout
- tesshunt/pipeline.py    the five moves + FFI pixel -> light curve
- tesshunt/sweep_core.py  one star end to end (shared by laptop and Modal sweeps)
- scripts/   runnable analysis scripts; modal_sweep.py is the parallel cloud sweep (Kishore runs it)
- data/      downloaded light curves and cutouts (gitignored)
- plots/     output images (gitignored)

## Reference
- TESS Instrument Handbook v0.1 (Dec 2018): ~/Downloads/TESS_Instrument_Handbook_v0.1-2.pdf
- Key handbook facts: 2-min TPFs for ~20k stars, 30-min FFIs for everything (10 min / 200 s in later sectors); 27-day sectors; 21 arcsec pixels; quality bit 8 = manual exclude; momentum dumps every 2.5-3 days; scattered light near end of each orbit; Tmag 10 = 15,000 e-/s; stars brighter than Tmag ~6 saturate.

## Pitfalls
- (symptom → cause → rule, one line each; add in the same session something goes wrong)
- Plots landed in ~ instead of plots/ → inline heredoc was run from the home directory → always run scripts/*.py from the project root; scripts write to plots/ by absolute path.
- Synthetic 300 ppm transit went undetected at 30-min cadence with 200 ppm noise → only 4 transits x 6 points gives SNR ~7, right at the threshold → for FFI targets expect to need either a bright star (low noise), a deeper dip, or multiple sectors; do not trust SNR < 8.
- TESScut FFI cutouts have quality == 0 everywhere, even in Sector 1's known bad-pointing stretch (days 1347-1349.6) → TESScut does not carry SPOC's cadence flags → always run mask_jitter (or a hand cut) on FFI light curves before searching.
- Jitter mask based on rolling std deleted every transit cadence of Pi Men → a 234 ppm dip on a 60 ppm-noise star raises local scatter 4x, same as jitter → detect jitter with a rolling median of point-to-point jumps, never with local scatter.
- select_targets.py crashed after a 25-min catalog pull → Exoplanet Archive `ps.tic_id` is text like "TIC 36734222", not an int → strip the prefix; and cache every slow network pull to data/ before any later step can fail.
- ExoFOP returned HTTP 522 (down) → CTOI filter silently skipped → script now says so and uses data/ctoi_tics.csv cache when present; rerun the CTOI check before trusting a candidate list.
- 25-min catalog pull lost when the Claude session restarted → background job was tied to the session and the cache only wrote at the end → long pulls save a .partial.csv after every chunk and resume from it; launch them with nohup so they outlive the session.
- Sector 1 target list had 44,401 rows but only 40,328 unique TICs → cone overlap survived the merge with tess-point output → drop_duplicates("ID") after every merge; sweep.py dedupes on read too.
- Sweep takes 5-50 s per star, download-bound → 2,900 tiny stars is ~1 day, full sector ~3 weeks sequential → parallelize downloads before sweeping the small/warm tiers.
- TESScut (NASA hosted cutouts) limits each user to 5 requests/s and 10 concurrent, returns 503 beyond that → use astrocut TessFootprintCutout against MAST's public S3 sector cubes instead (no rate limit; sweep_core source="s3"); run it near us-east-1 (Modal) for speed.
- Two sweeps appending to the same data/sweep_s<N>.csv can interleave rows → never run the laptop sweep and the Modal sweep at the same time on the same sector.
- First sweep "candidate" (TIC 52153067, SNR 29) was a spotted rotator: P=0.49 d, 6 h duration, secondary -1574 ppm → SNR+depth alone do not distinguish a sine wave from a box → candidate bar also requires duration <= 1.5x physical max, duty cycle < 15%, secondary ~0, odd = even.
- Killing the local `modal run` client left the ephemeral app and its 50 workers running in the cloud → two apps then hit the 100-container account limit → stop Modal apps with `modal app stop <app-id>` (see `modal app list`), never with pkill. Workspace limit is 100 containers: keep TESS_MAX_CONTAINERS at 50 and run one sweep at a time.
- Vizier (CDS) cone searches time out from this network; SIMBAD hangs too → use astroquery.gaia (ESA archive) for neighbors, set socket.setdefaulttimeout, and run catalog lookups in the background.
- A raw-counts difference image peaks on whichever bright star is in the frame → normalize the difference image by per-pixel noise (sigma map) before taking the peak.
- Difference image for TIC 114949540 sector 105 peaked 8 px off target as a diffuse glow → scattered-light gradient across the cutout, and per-frame median background subtraction cannot remove a tilt → centroid_test now subtracts a per-frame plane fit to the faint pixels; a diffuse (not point-like) diff-image peak means background, not a neighbor.
- ExoFOP returned 522 from home and 525 from a US Modal container → the site itself was down, not the network → probe from Modal before blaming local connectivity.
- TIC 259719966 transit depth varied 2.5% to 4.4% between sectors → the threshold aperture sometimes included a neighbor 2 px away (27% of the light), diluting the dip → per-sector depth tracks aperture contents; fit radius on undiluted sectors only (mcmc_fit --sectors), and note that dilution going the "shallower" way confirms the dip is on the target.
- Two single-sector candidates with SNR 31 and 45 (TIC 306741052, 212347535) failed multi-sector recovery or showed a secondary → single-sector vetting is not enough; always run followup.py before ranking.
- Sector 6 sweep: 76% of stars failed with "boolean index did not match indexed array" → the sector 6 S3 cubes contain an all-NaN frame; remove_nans() dropped it from the light curve but the regressor matrix kept it → mask non-finite cadences in BOTH before RegressionCorrector; check the ok-fraction in `--mode status` after every new sector before trusting its candidate count.
- TIC 259719966 passed every test and went into the RNAAS as a 0.94 Rjup candidate; deeper vetting found a 0.9% secondary eclipse at phase 0.54 (45 sigma, confirmed in MIT QLP data at 22 sigma) → eccentric orbit (e = 0.27) moved the eclipse away from phase 0.5 and every secondary test only looked at 0.5 → vet() now scans all phases 0.15-0.85 and the candidate bar rejects any-phase secondaries > 4 sigma; a transit-derived density that disagrees with the catalog is a hint of eccentricity, not just degeneracy.
- A manuscript sentence like "no secondary eclipse" must name the phase range tested; the RNAAS (AAS82771) said it flatly and needs correction.
- overnight.py called 5/6 and 6/7 sector recoveries "failed_recovery" → the rule was recovered == usable, but one sector with heavy scattered light can miss a real signal → treat recovered >= max(2, usable-1) as "repeats", then look at the missing sector by hand.
- Two "packaged" stars had empty package folders and NaN radii → the MCMC subprocess failed silently inside the chain and the regex just found nothing → the chain must check the mcmc csv exists before calling a star packaged; refits with --max-b 1.0 were needed for three grazing/diverged fits (b~0.95, rho_star nonsense).
- BLS duration came out exactly 2.4 h → the duration grid was coarse (0.05 d steps) → use a fine duration grid (e.g. 0.5 h to 6 h in 0.25 h steps) before trusting duration.
- TIC 259719966 passed "secondary at phase 0.5: 20 +/- 540 ppm" through RNAAS submission, yet a 0.9% secondary eclipse sits at phase 0.54 in all 9 sectors → every secondary test (vet, oddeven_fit, followup) only looked at phase 0.50, and an eccentric orbit moves the eclipse → scan a box of the transit duration across ALL phases (deep_vet.py secondary scan) before calling anything a planet; a flat-bottomed 4% dip on an M dwarf is an EB until the whole phase curve is clean.
- Sector 103 centroid shift 35 sigma, difference image 25" off → one "transit" fell in the mid-sector scattered-light ramp where every cadence is quality-flagged; the light curve pipeline drops them, centroid_test did not → run centroid_test on tpf[tpf.quality == 0] and per transit, and treat a many-sigma result that only one epoch produces as a data artifact until checked.
- Sector-1 centroid verdict "anomaly, 6x the on-target prediction" → the prediction summed Gaia neighbours within 2.0 px and the dominant neighbour sits at 2.08 px → predict centroid shifts from every Gaia star in the frame (whole-frame centroid), and compare direction as well as size: an on-target dip moves the centroid toward the neighbours, an off-target one away from the dipping star.
- TRICERATOPS gave TIC 259719966 (the eclipsing binary) FPP = 0.0000 and P(planet) = 1.00 → statistical validation only models the primary transit and the field, never the secondary eclipse → never quote a TRICERATOPS FPP without the all-phase secondary scan beside it; validation is the last test, not the first. Also: its default field image is a TESScut call that fails for late sectors and near CCD edges, so pass --field-sectors with early sectors.
