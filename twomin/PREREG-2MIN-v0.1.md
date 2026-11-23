# Pre-registration 2MIN v0.1: does Parker Solar Probe's "2-minute dominance" hold in encounters the claim never used?

MV2 of Munim, Inc. (mv2@munimtech.com). Written 2026-11-23, before job `psp2min-v1` has run and
before any f_mid value from any Parker Solar Probe window has been computed by this project.
Made with AI agents working in MV2's name; no human expert has checked it yet. Licence: text
CC BY 4.0, code Apache-2.0.

## The claim

Huang, Velli, Shi et al., *Dominance of 2-Minute Oscillations near the Alfven Surface*, ApJ Letters
977, L12 (2024), doi:10.3847/2041-8213/ad9271, arXiv:2405.15967. From 4-hour windows of the first
19 encounters (E1-E19) they report that 1/f_mid, the period at the middle of the band holding the
central 50 % (25 %-75 % cumulative) of the magnetic fluctuation energy, falls as the advection time
t_adv = (R - Rsun)/V_r falls, reaching about 2-3 minutes in the youngest wind (smallest t_adv
about 5.7 h), and they argue that the true floor is "likely closer to 2 minutes", the launch period
of Alfven waves seen in SDO/AIA swaying motions. They read this as evidence that the waves start in
the solar atmosphere. No significance test or out-of-sample check was made. The same group's 2026
paper (arXiv:2511.10906) reports 1/f_mid = 105 s for one interval at 9.86 Rsun.

## Why it matters

If the energy-containing period of near-Sun turbulence is pinned at a coronal launch period, the
solar-wind fluctuation spectrum carries a direct imprint of the low corona. If instead it keeps
shrinking as the wind gets younger, the "2 minutes" is where turbulence happened to be at the
distances E1-E19 reached, not a floor. Encounters E20-E27 (June 2024 to March 2026; E22-E27 at
9.86 Rsun) reach the deepest and youngest wind ever sampled and were never used for the claim.
They are a clean out-of-sample test.

## What we run (job `psp2min-v1`, code `twomin/psp2min.py` 0.1.0)

Reused unchanged, fetched at pinned commits with a fatal sha256 check: `psp5min/p5min.py` 0.1.0
(commit a025a506..., sha256 310727b8...; positions, perihelion numbering checked against E22 and E24,
FIELDS `PSP_FLD_L2_MAG_RTN_4_SA_PER_CYC` binned to 0.5 s) and `alfven/psp5alf.py` 0.1.0 (commit
f13f8d6d..., sha256 6fc3e1fb...; SPAN-i `PSP_SWP_SPI_SF00_L3_MOM` via HAPI).

- **Windows.** For every encounter E1-E27 whose +-7 day span ends before 2026-04-30 (the end of
  the magnetometer data), tiled (non-overlapping) 4-hour windows over the +-7 days. Tiled windows,
  not the paper's 15-minute sliding steps, so that the windows used in statistics are independent.
- **Magnetometer.** A window is valid if at most 5 % of its 0.5 s bins are empty and no gap exceeds
  60 s (gaps filled linearly).
- **f_mid (our reading of the paper's metric).** Per window: remove a least-squares line from each
  RTN component, Hann taper, FFT, trace PSD. Cumulative energy of P(f) between f_lo and 1 Hz
  (the 0.5 s Nyquist), where f_lo = 0.380 mHz is the frequency at which half of a Morlet (w0 = 6)
  wavelet's points fall outside the cone of influence for a 4-hour window (the paper's lower limit).
  f25 and f75 are the first FFT frequencies where the cumulative energy reaches 25 % and 75 %;
  **f_mid = sqrt(f25 f75)** (the middle on the paper's log-frequency axis; primary). The arithmetic
  middle (f25 + f75)/2 is reported, not graded. Also reported: 1/f_peak, the maximum of P(f) f
  smoothed over a factor-2 band, between f_lo and 50 mHz.
  Check K0 (printed and unit-tested): for pure f^-5/3 noise our 1/f_mid must lie in 600-900 s (the
  paper's Appendix A gives about 13 min, our continuum value is 758 s).
- **Advection time.** V_r = median SPAN-i radial speed (Sun frame, 50-2,000 km/s) in the window;
  usable if at least 60 % of the window's 1-minute slots hold a sample. R from the hourly position
  at the window centre. t_adv = (R - 1 Rsun) x 695,700 km / V_r.
- **Cleaning (our reading of their 0.005 < dt_adv/<t_adv> < 0.2).** dt_adv/<t_adv> = standard
  deviation over mean of t_adv computed from 10-minute median speeds (at least 6 blocks). The cut
  is applied to the graded samples **unless it removes more than 15 % of usable E1-E19 windows**
  (the paper removed 3.7 %), in which case all grades use the uncut samples. Both are reported.
- A window is **usable** if the magnetometer is valid, V_r is usable and f_mid exists.
- Sample A = usable windows of E1-E19 (the paper's encounters); sample B = E20-E27 (new).

## Gates and tests (fixed now)

- **P0 pipeline**: perihelion numbering correct (E22, E24 within 2 days of the expected dates); at
  least 15 of E1-E19 and at least 6 of E20-E27 have >= 50 % usable windows. P0 FAIL: H1 and H2 not
  gradable (numbers reported as indicative).
- **R1 reproduction on E1-E19**: Spearman correlation between t_adv and 1/f_mid > 0.3 AND the
  geometric mean 1/f_mid of windows with t_adv <= 10 h lies in 90-240 s (at least 20 such windows,
  else NOT GRADABLE). R1 FAIL does not stop grading; H1 and H2 are then reported "with our reading
  of the paper's metric".
- **H1 (main test): level in the youngest new wind.** G = geometric mean of 1/f_mid over sample-B
  windows with t_adv <= 10 h, with an encounter-block bootstrap 95 % interval (10,000 resamples,
  seed 20261123).
  - **CONSISTENT** with 2-minute dominance if 90 s <= G <= 180 s;
  - **SHORTER** if G < 90 s; **LONGER** if G > 180 s;
  - NOT GRADABLE if fewer than 30 such windows or fewer than 3 encounters contribute.
- **H2: floor or continued decline.** Pooling A and B: ratio = G(t_adv <= 7 h) / G(7 h < t_adv <= 10 h).
  **FLOOR** if ratio >= 0.85; **DECLINING** if ratio <= 0.70; **INTERMEDIATE** between; NOT GRADABLE if
  either group has fewer than 20 windows. Reported, not graded: the same ratio in B alone.
- Reported, not graded: binned tables of 1/f_mid against t_adv (edges 0, 6, 7, 8, 10, 12, 15, 20, 25,
  30, 40, 60 h) for A and B (to set beside the paper's Fig. 4c), per-encounter values, the fraction of
  young windows whose P(f) f peak sits at f_lo, and 1/f_peak.

## Predictions (written before the run)

- P0 PASS 80 %. K0 PASS 90 %. R1 PASS 55 % (the paper's wavelet averaging, cleaning and "middle"
  are not fully specified).
- H1: CONSISTENT 55 %, SHORTER 25 %, LONGER 20 %. Reasoning: the one published 9.86 Rsun value (105 s)
  sits inside the band, but the slow, dense wind near recent perihelia may sit lower.
- H2: FLOOR 40 %, INTERMEDIATE 30 %, DECLINING 30 %.

## What each outcome means

- CONSISTENT + FLOOR: first out-of-sample support for a roughly 2-minute energy-containing period in
  the youngest wind, in data the claim never saw.
- SHORTER or DECLINING: the period keeps shrinking below 2 minutes in younger wind, so "2 minutes" is
  not a floor; the spectral argument for a coronal launch period would need revisiting (the paper's
  polarisation and wavelet-burst evidence is not tested here).
- LONGER: the youngest new wind does not reach 2-3 minutes; the trend does not carry over.
- Limits fixed now: FFT with a Hann taper is not the paper's cone-averaged wavelet spectrum;
  t_adv uses partial SPAN-i moments (quality flags not used); spacecraft-frame frequencies are not
  Doppler-corrected (the paper argues the correction is ~1 near the Alfven surface).

Disclosure: this project has already seen the results of `psp5alf-v1` (Alfvenicity of E17-E27
windows near the two 5-minute events, results committed in this repository), which contain no
f_mid or t_adv values. No other output of E20-E27 data has been read. No threshold, sample or
definition in this file will change after the data are seen; any unforeseen case gets a numbered
amendment written before results are read.
