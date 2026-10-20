# tess-hunt

Transit photometry search over TESS data. Goal: an agent that finds unreported planet candidates in full-frame images (FFIs), vets them, and submits worthwhile ones as Community TOIs on ExoFOP.

Kishore is new to astronomy. Explain concepts plainly, ramp up gradually, draw a diagram when it helps.

## Plan
1. Learn on a known planet: Pi Mensae c, Sector 1, 2-minute data (scripts/pimen_s1.py).
2. Build the five-step pipeline as reusable functions: clean, flatten, BLS search, fold, vet.
3. Move to FFI cutouts (TESScut) of stars not on the 2-minute target list.
4. Cross-match candidates against TOI / CTOI / confirmed-planet lists.

## Layout
- scripts/   runnable analysis scripts
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
- BLS duration came out exactly 2.4 h → the duration grid was coarse (0.05 d steps) → use a fine duration grid (e.g. 0.5 h to 6 h in 0.25 h steps) before trusting duration.
