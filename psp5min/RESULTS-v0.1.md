# Results v0.1: how often does Parker Solar Probe show a "6 sigma" 5-minute peak by chance?

MV2 of Munim, Inc. (mv2@munimtech.com). Graded 2026-11-25 against [`PREREG-v0.1.md`](../PREREG-v0.1.md)
(written 2026-11-21, before any data were read) and amendments
[0.1.1](../PREREG-amendment-0.1.1.md), [0.1.2](../PREREG-amendment-0.1.2.md) and
[0.1.3](../PREREG-amendment-0.1.3.md). Made with AI agents working in MV2's name; no human expert
has checked it yet. Text CC BY 4.0, code Apache-2.0. No threshold, band, window length or sample
was changed after the data were seen.

Job `psp5min-v1` (code `psp5min/p5min.py` 0.1.0, sha256 310727b8...; 15/15 unit tests pass) ran
2026-11-25 14:25-15:29 UTC on public NASA CDAWeb data (FIELDS `PSP_FLD_L2_MAG_RTN_4_SA_PER_CYC`,
`PSP_HELIO1HR_POSITION`), encounters E1-E27, +-7 days around each perihelion. Every output is in
[`results/psp5min-v1/`](../results/psp5min-v1/) (`results.json`, `run.txt`, `clusters_5p8.csv`,
`published_windows.csv`, `perihelia.csv`, logs). The per-window tables (`windows/E01-E27.csv`, about
4,100 rows each) are not committed because of their size; the code rebuilds them from the public data.

The claim: Huang, Velli, Panasenco et al., ApJL 999 L4 (2026), arXiv:2511.10906: "~6 sigma" peaks at
3.2 mHz (E22, 2024-12-24 09:50-12:05, 5.8 sigma) and 3.1 mHz (E24, 2025-06-19 08:15-09:25, 6.2 sigma),
read as the Sun's 5-minute p-modes seen in situ.

## Grades

| Test | Measured | Grade | Our prior |
|---|---|---|---|
| P0 pipeline | numbering correct (E22 2024-12-24T11, E24 2025-06-19T09); 24/24 of E1-E24 with >= 70 % valid 135-min windows; 3/3 of E25-E27 >= 50 % | **PASS** | PASS 75 % |
| R0 reproduction | E22: peak 3.210 mHz, S = 5.31 (published 5.8); E24: peak 3.095 mHz, S = 6.15 (published 6.2) | **PASS** | PASS 45 % |
| **H1 calibration (main test)** | clusters with S >= 5.8 in the flanking control bands (2.1 + 4.8 mHz), E1-E24: **kF = 12 at 135 min** (8 + 4), kF = 1 at 70 min. Rule: CHANCE-LEVEL if kF >= 4 at either length | **CHANCE-LEVEL** | CHANCE-LEVEL 60 % |
| H2 frequency specificity | binding grade comes from `psp5curv-v1` (H2c, amendment 0.1.1), still queued. psp5min's own L rank test: Z = 2.69 (indicative only) | **pending** | - |

H1's guard in amendment 0.1.1 applies only to RARE; CHANCE-LEVEL is unchanged by it.

## H1 in detail (S >= 5.8 clusters per band, E1-E24, published intervals excluded)

| band (mHz) | 3.15 (paper) | 1.3 | 2.1 | 4.8 | 6.6 | 7.4 | 8.2 | 9.0 | 9.8 | controls mean |
|---|---|---|---|---|---|---|---|---|---|---|
| 135-min windows | 2 | 2 | 8 | 4 | 0 | 3 | 2 | 4 | 1 | 3.0 |
| 70-min windows | 0 | 0 | 0 | 1 | 1 | 0 | 2 | 3 | 1 | 1.0 |

Sliding windows: 48,072 (135 min; 35,844 eligible) and 48,240 (70 min; 31,885 eligible) in E1-E24.
At frequencies where nobody expects a p-mode, the paper's own statistic reaches 5.8 about three
times per band at the E22 window length and about once per band at the E24 length, within the
authors' own search space. The two published peaks are what a search of this size produces by chance.

Two examples from `results.json` (information only):
- E22, 2024-12-24 21:50 UTC, 12.4 Rsun: a 2.1 mHz cluster peaking at S = 7.23, the same day as the
  published 3.2 mHz event (09:50-12:05) and higher than its published 5.8.
- E24, 2025-06-19 06:30 UTC, 9.99 Rsun: a 70-min window with S = 5.95 at 4.8 mHz, 1 h 45 min before
  the published 3.1 mHz window on the same morning.

## Reported, not graded

- Paper band in E1-E24 outside the published intervals (135 min): 2 clusters, E15 2023-03-17 14:10
  (13.6 Rsun, S 6.25) and E22 2024-12-25 07:40 (16.5 Rsun, S 5.84).
- **H3 out of sample (E25-E27, never used by the authors), 135 min:** paper band 2 clusters (both E27:
  2026-03-10 08:50, 21.7 Rsun, S 6.71; 2026-03-12 05:40, 12.9 Rsun, S 5.85); flanks 2.1 mHz 2 and
  4.8 mHz 2; mean of the eight controls 0.5. At 70 min: paper band 0, controls mean 0.25.
- L rank test (tiled 135-min windows; superseded as the graded H2 by amendment 0.1.1, indicative
  only until `psp5curv-v1` measures L's built-in bias in simulation): all Z = 2.69 (n = 3,931); E1-E24
  2.60; E25-E27 0.70; r <= 30 Rsun 0.95. 70-min windows: all -1.00, E25-E27 -2.22. Tail test (135 min):
  52 target exceedances of the control 99th percentile vs 315 for eight controls, one-sided binomial
  p = 0.041; 70 min p = 0.32.
- The published windows' L sits at the 99.37th (E22) and 99.94th (E24) percentile of control-band L
  values. Both windows were selected for having a peak, so these percentiles are not a calibration.

## What this means

The paper's "~6 sigma" is a single-spectrum z-score (bounded by sqrt(n_b - 1), amendment 0.1.1), not a
false-alarm probability. Measured on the authors' own encounters, the same statistic gives peaks of
that size at control frequencies at least as often as at 3.1-3.2 mHz, so the spectral peaks alone are
not evidence for p-modes. This does not show p-modes are absent: the indicative rank test leans
slightly toward 3 mHz in E1-E24 (Z 2.6, not seen in E25-E27 or within 30 Rsun), and its binding,
bias-corrected version (H2c), the surrogate false-alarm test (H4) and the reading-sensitivity ladder
are still queued. Our earlier Alfvenicity test ([RESULTS-ALF-v0.1](../alfven/RESULTS-ALF-v0.1.md))
found the published wave trains highly Alfvenic but not more so than typical fluctuations within
20 Rsun.

Limits: sliding windows in 10-minute steps are not tuned to each event, which makes chance peaks
slightly harder to find (it favours the claim); band centres and the 25-75 % baseline follow our
reading of the paper (R0 shows it reproduces both published values within 0.5); data end 2026-04-30.
