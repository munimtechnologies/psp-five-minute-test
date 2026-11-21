# Pre-registration ALF v0.1: does Alfvenicity single out the published 5-minute wave trains?

MV2 of Munim, Inc. (mv2@munimtech.com). Written 2026-11-21, before any SPAN-i plasma data or
magnetometer data were downloaded or looked at by this project (job `psp5min-v1`, the first
test, has not run yet either). Made with AI agents working in MV2's name; no human expert has
checked it yet. Text CC BY 4.0, code Apache-2.0.

## Why a second test

Huang, Velli, Panasenco et al. (ApJ Letters 999, 2026; arXiv:2511.10906) rest their p-mode
reading on two legs. The first is the spectral "~6 sigma" peak; `PREREG-v0.1.md` calibrates it.
The second is Alfvenicity: for the E22 interval (2024-12-24 09:50-12:05 UTC) they report cross
helicity sigma_c "close to unity", both in time (fluctuations after subtracting a 10-minute
rolling mean) and in frequency (sigma_c(f) from FFTs of V and V_A), a nearly constant |B|, and
power concentrated in B_N. The paper shows no comparison intervals for these signatures.

Inside about 20 solar radii the solar wind is known to be highly Alfvenic almost everywhere. If
ordinary fluctuations at nearby frequencies, and chance peaks of the paper's own statistic, are
just as Alfvenic, then Alfvenicity is consistent with p-modes but does not tell them apart from
ordinary turbulence. If the published peaks are far more Alfvenic than both, that is real
support. Nobody has measured this.

## Data and code

- `alfven/psp5alf.py` 0.1.0, importing `psp5min/p5min.py` 0.1.0 unchanged (the job downloads it
  from commit a025a506d2092e982e82e87c2c27e4f2a477fa98 and prints its sha256).
- B: FIELDS `PSP_FLD_L2_MAG_RTN_4_SA_PER_CYC`, binned to 0.5 s exactly as in psp5min, then
  averaged to a 2 s grid. V and density: SWEAP SPAN-i `PSP_SWP_SPI_SF00_L3_MOM` (`VEL_RTN_SUN`,
  `DENS`, about 1.75 s cadence) via CDAWeb HAPI, linearly interpolated to the 2 s grid where
  the two bracketing samples are at most 10 s apart. No quality-flag filtering (the paper
  describes none); non-finite values are dropped.
- Encounters E17-E27, +-3 days around perihelion; only windows whose centre is within
  **20 solar radii** are used. Windows of 135 and 70 minutes, sliding in 10-minute steps, plus
  tiled windows. Windows overlapping either published interval are excluded from every
  reference set.
- A window is usable if the psp5min statistic is valid (its gap rules and baseline rule) and,
  on the 2 s grid, B and V each have at most 5 % missing points and no gap over 60 s, and at
  least 95 % of densities are finite and positive.
- Alfven units: b = B x 21.812 / sqrt(n) km/s, n = the window's median SPAN-i density.
- sigma_c(t) = 2<dv.db> / (<dv^2> + <db^2>) after subtracting a centred 301-point (10-minute)
  rolling mean from each component of v and b.
- Band sigma_c: from FFTs of mean-removed v and b (no taper), 2 sum Re(V.B*) / sum(|V|^2+|B|^2)
  over the same FFT bins psp5min scores (paper band P3.15: 2 bins at 135 min, 1 at 70 min;
  controls at 1.3, 2.1, 4.8, 6.6, 7.4, 8.2, 9.0, 9.8 mHz). **Near controls = 1.3, 2.1, 4.8 mHz.**
  B_N share = |B_N|^2 / |B|^2 over the same bins. Sign is ignored (|sigma_c|), because the sign
  only follows the field polarity.
- Published values: c_pub(E22) = |band sigma_c| at P3.15 in 2024-12-24 09:50-12:05 (135 min);
  c_pub(E24) = the same in 2025-06-19 08:15-09:25 (70 min). c_min = the smaller of the two.

## Gates and tests (fixed now)

- **R1 (reading check):** |sigma_c(t)| >= 0.8 in both published windows (the paper: "close to
  unity"). R1 FAIL does not stop grading; A1/A2 are then reported "with our reading of their
  sigma_c".
- **A1 (are the published peaks more Alfvenic than ordinary near-Sun fluctuations?).** For each
  published event, q = fraction of reference values strictly below c_pub, where the reference is
  |band sigma_c| in the three near-control bands over all usable tiled windows of the same length.
  Needs at least 30 reference values for each event.
  - **NOT DISTINCTIVE** if q < 0.90 for both events.
  - **DISTINCTIVE** if q >= 0.99 for both events.
  - **MIXED** otherwise; **NOT GRADABLE** if a published window or the reference is missing.
- **A2 (are chance peaks as Alfvenic as the published ones?).** A chance peak is a cluster of
  consecutive usable sliding windows (same encounter, same band, same length) in which a
  near-control band has the paper's statistic S >= 4.0; the window with the largest S
  represents it. Both lengths and all three near-control bands are pooled. f = fraction of
  chance peaks with |band sigma_c| >= c_min. Needs at least 5 chance peaks.
  - **AS ALFVENIC AS PUBLISHED** if f >= 0.5.
  - **LESS ALFVENIC** if f < 0.2.
  - **INTERMEDIATE** otherwise; **NOT GRADABLE** if fewer than 5 peaks or c_min is missing.
- Reported, not graded: q for all eight controls; q against the paper band in other windows;
  q of sigma_c(t); B_N share of the published peaks against near-control bands and chance peaks;
  residual energy; dB/B; the same fraction f for chance peaks with S >= 5.8; cluster counts at
  S >= 5.8 and S >= 4 in the paper band and every control within 20 solar radii (a small
  out-of-sample companion to H1 of PREREG-v0.1 at the published events' own distance).

## Predictions (written before the run)

- R1 PASS 60 % (SPAN-i partial moments near perihelion can miss part of the beam; that lowers
  |sigma_c| for every window, published or not).
- A1: NOT DISTINCTIVE 65 %, MIXED 25 %, DISTINCTIVE 10 %.
- A2: AS ALFVENIC AS PUBLISHED 55 %, INTERMEDIATE 25 %, LESS ALFVENIC 10 %, NOT GRADABLE 10 %.

## What each outcome means

- NOT DISTINCTIVE and AS ALFVENIC AS PUBLISHED: the Alfvenicity of the two wave trains is what
  any near-Sun fluctuation, and any chance peak, shows; it supports "Alfven waves" but adds no
  evidence for a 5-minute solar origin. Combined with a CHANCE-LEVEL H1 the claim would then rest
  on no calibrated evidence. It would not show that p-modes are absent.
- DISTINCTIVE and LESS ALFVENIC: the published wave trains stand out in Alfvenicity too, which
  strengthens the claim, whatever H1 says.
- Limits fixed now: one density (the window median of SPAN-i partial moments) per window, where
  the paper used QTN electron density for E22; band sigma_c from one or two FFT bins is noisy for
  every window alike; sliding windows are not tuned to events.

No threshold, band, window length or sample in this file will be changed after the data are
seen. Any unforeseen case is handled in a numbered amendment written before results are read.
