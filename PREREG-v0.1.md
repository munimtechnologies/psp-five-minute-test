# Pre-registration v0.1: how often does Parker Solar Probe show a "6 sigma" 5-minute peak by chance?

MV2 of Munim, Inc. (mv2@munimtech.com). Written 2026-11-21, before any magnetometer data
were downloaded or looked at by this project. Made with AI agents working in MV2's name; no human
expert has checked it yet. Licence: text CC BY 4.0, code Apache-2.0.

## The claim

Huang, Velli, Panasenco et al., *Evidence of 5 minute Oscillations from Parker Solar Probe*,
ApJ Letters 999 (2026), doi:10.3847/2041-8213/ae43ec, arXiv:2511.10906 (v1, 2025-11-14).
They report "statistically significant (~6 sigma)" peaks at 3.1-3.2 mHz in the magnetic-field
power spectrum, which they read as the Sun's 5-minute p-mode oscillations seen in situ for
the first time:

| event | interval (UTC) | peak | "sigma" |
|---|---|---|---|
| E22 | 2024-12-24 09:50-12:05 (135 min) | 3.2 mHz | 5.8 |
| E24 | 2025-06-19 08:15-09:25 (70 min) | 3.1 mHz | 6.2 |

Their method, as given in arXiv v1: FIELDS fluxgate B (RTN) resampled to 0.5 s; NumPy FFT of each
component over the whole interval; trace; spectrum multiplied by f; the "sigma" is
(peak - mean) / std of P(f)*f over "the frequency range containing the central 50% of the total
fluctuation energy". The search covered all data within +-7 days of perihelion in E1-E24.
The paper reports no false-alarm rate, no control intervals and no correction for the number
of intervals or frequencies searched. We found no published comment or replication
(searched 2026-11-21).

## Why it matters

If the Sun's p-modes reach 10 solar radii as Alfven waves, that ties helioseismology to the solar
wind. The published version already says the waves must largely dissipate near the Alfven
surface. Whether two peaks of this size are rare, or are expected in a search of this size, has
never been measured. A single-spectrum "sigma" is not a probability: single FFT bins scatter
like chi-square with few degrees of freedom, and solar-wind turbulence is intermittent.

## What we run (job `psp5min-v1`, code `psp5min/p5min.py` 0.1.0)

Data: NASA CDAWeb `PSP_FLD_L2_MAG_RTN_4_SA_PER_CYC` (public) and `PSP_HELIO1HR_POSITION` for
perihelion times (each perihelion = minimum of hourly distance within a run of r < 0.25 AU,
numbered from launch). Encounters E1-E27: E1-E24 are the authors' search space; **E25-E27
(Sept 2025-Mar 2026, all at 9.86 Rs) were never used by them**. Data end 2026-04-30.

For each encounter, the +-7-day span is binned to 0.5 s (mean per bin). Windows of 135 and
70 minutes (the two published lengths) slide in 10-minute steps; tiled (non-overlapping)
windows are also kept. A window is valid if at most 5 % of its 0.5 s bins are empty and no gap
exceeds 60 s (gaps filled linearly). Windows that overlap a published interval are excluded
from every count below.

Statistic S (the paper's): as above, with no taper, mean removed, bins k = 1..N/2-1; baseline =
bins from where the cumulative energy first reaches 25 % to where it first reaches 75 %
(at least 10 bins, otherwise the window is excluded). We score the n_p FFT bins nearest a
band centre (n_p = 2 for 135 min, 1 for 70 min, i.e. the bins inside 3.05-3.25 mHz), the paper
band P = 3.15 mHz, and eight control bands at 1.3, 2.1, 4.8, 6.6, 7.4, 8.2, 9.0, 9.8 mHz (the
3-minute range 5.0-6.2 mHz and 2.3-4.4 mHz are kept out of the controls). The flanking
controls are 2.1 and 4.8 mHz.

A **cluster** is a run of consecutive sliding windows (10-min steps) in which a band's S is
at least 5.8. Clusters, not windows, are counted, so one event seen by overlapping windows
counts once.

Statistic L (ours, frequency-local): y = P*f at a bin divided by the median of y over bins
2..w away on both sides (w = max(6, ceil(0.3k))). Per tiled window, L is maximised over the
n_l bins nearest 3.2 mHz (n_l = 5 for 135 min, 3 for 70 min) and over the same-size control
bands; the target's rank among the 9 bands (1 = lowest, 9 = highest) is recorded.

## Gates and tests (fixed now)

- **P0 pipeline**: perihelion numbering puts E22 and E24 within 2 days of 2024-12-24 and
  2025-06-19; at least 20 of E1-E24 have >= 70 % valid sliding 135-min windows; all of
  E25-E27 have >= 50 %. P0 FAIL makes H1 and H2 not gradable (numbers reported as indicative).
- **R0 reproduction**: on the two published intervals, the largest P*f bin in 2.5-4.0 mHz lies
  in 3.00-3.35 mHz and its S is within +-1.5 of the published 5.8 / 6.2. R0 FAIL does not stop
  the grading, but H1 is then reported as "with our reading of their statistic".
- **H1 calibration (main test).** kF = clusters with S >= 5.8 in the two flanking control bands
  (2.1 + 4.8 mHz), summed, in E1-E24 (the authors' own search space).
  - **CHANCE-LEVEL** if kF >= 4 at 135 min or at 70 min: at frequencies where no p-mode is
    expected, the paper's statistic produces at least two "~6 sigma" peaks per band in the
    authors' search, as many as the two published detections.
  - **RARE** if kF <= 1 at both lengths.
  - **UNDECIDED** otherwise.
  Also reported, not graded: the mean over all eight controls, the paper-band count (kT),
  counts at S >= 4, and every cluster with time and distance.
- **H2 frequency specificity.** Rank test on tiled 135-min windows, all encounters, published
  intervals excluded: Z = (mean rank - 5) / sqrt(6.667 / n).
  **SUPPORT** (3 mHz special) if Z >= 3; **NO-EXCESS** if Z < 2; **INDETERMINATE** between.
  Reported, not graded: Z for 70 min, for r <= 30 Rs, for E1-E24 and E25-E27 separately; a tail
  test (target exceedances of the pooled control 99th percentile of L, binomial p with 1/9); and
  where the two published windows' L falls among all control-band values.
  Power caveat, fixed now: H2 can only detect a 3 mHz signal present in a sizeable fraction of
  windows (roughly 10 % or more). NO-EXCESS does **not** rule out rare events like the two
  published ones; H1 is the test of those.
- **H3 out-of-sample (reported, not graded):** clusters with S >= 5.8 in the paper band in
  E25-E27, compared with E22-E24.

## Predictions (written before the run)

- P0 PASS 75 %. R0 PASS 45 % (the paper leaves the taper, detrending and band edges open).
- H1: CHANCE-LEVEL 60 %, UNDECIDED 25 %, RARE 15 %. Reasoning: with 2 FFT bins per band and
  3,584 independent 135-min windows in E1-E24 +-7 d, an exponential (2 degree-of-freedom)
  tail gives about 8 chance peaks per band; a 6 degree-of-freedom tail gives about 1.5;
  intermittent turbulence should push toward the heavier tail.
- H2: NO-EXCESS 70 %, INDETERMINATE 20 %, SUPPORT 10 %.

## What each outcome means

- CHANCE-LEVEL: the published "~6 sigma" does not measure how unusual the peaks are; two peaks of
  that size at 3.1-3.2 mHz are expected by chance in a search of E1-E24. That would not prove
  p-modes are absent (their Alfvenic-polarisation argument is untouched); it would mean the
  spectral evidence needs a calibrated test.
- RARE: peaks of that size are rare at other frequencies too, which strengthens the claim.
- Our sliding windows are not tuned to each event, so they find somewhat lower peaks than
  hand-chosen intervals would; this makes CHANCE-LEVEL harder to reach (it favours the claim).

No threshold, band, window length or sample in this file will be changed after the data are
seen. Any unforeseen case is handled in a numbered amendment written before results are read.
