# Results of pre-registration ALF v0.1: is the Alfvenicity of the published 5-minute wave trains unusual?

MV2 of Munim, Inc. (mv2@munimtech.com). Graded 2026-11-22. Made with AI agents working in MV2's name;
**no human expert has checked it yet**. Text CC BY 4.0, code Apache-2.0.

Job `psp5alf-v1` ran on the Blade from 16:59 to 17:17 UTC on 2026-11-22 (1,106 s, peak 0.19 GB). It ran
psp5alf 0.1.0 (sha256 6fc3e1fb...c18d) with psp5min's p5min.py 0.1.0 downloaded at commit a025a506
(sha256 310727b8...05a4). py_compile and the unit tests passed. Encounters E17-E27 were read, with
694-697 windows per encounter, and every window passed the Alfven-data rules except 43 in E17.
Raw output: `results/psp5alf-v1/` (`run.txt`, `results.json`, `published.csv`, `peaks.csv`, `tests.txt`,
`sha256.txt`, `versions.txt`, `job.json`). The per-window table `windows.csv` (7,652 rows) is not
committed because of its size. It can be rebuilt by re-running the job.

Grades are copied from the job's `PALF_GRADES` line and checked by hand against PREREG-ALF-v0.1.md.
No threshold was changed.

## Reading check

We first re-measured the paper's own statistic at the paper's own band, which psp5min also uses
(`S_P3.15`). It gives **6.15 for E24** (paper: 6.2 sigma) and **5.31 for E22** (paper: 5.8 sigma).
The peak frequencies are 3.10 and 3.21 mHz (paper: 3.1 and 3.2).

**R1 (|sigma_c(t)| >= 0.8 in both published windows): PASS** (prior 60 %). The values are 0.958 (E24)
and 0.939 (E22), matching the paper's "close to unity".

## A1: are the published peaks more Alfvenic than ordinary near-Sun fluctuations?

Here c_pub is the |band cross helicity| at the paper band in the published window. q is the fraction
of tiled windows (same length, r <= 20 solar radii, published windows excluded) whose near-control-band
value (1.3, 2.1, 4.8 mHz) is lower.

| Event | c_pub | reference values | q | median reference |
|---|---|---|---|---|
| E24 (70 min) | 0.977 | 1,452 | 0.801 | 0.943 |
| E22 (135 min) | 0.967 | 753 | 0.781 | 0.935 |

**A1 = NOT DISTINCTIVE** (q < 0.90 for both). The prior for this outcome was 65 %.

Not graded: q is 0.77 / 0.75 against all eight control bands and 0.77 / 0.72 against the paper band in
other windows. For the time-domain sigma_c, q is 0.76 / 0.59. For the B_N power share, q is 0.90 /
0.95 against the near-control bands.

## A2: are chance peaks of the paper's statistic as Alfvenic as the published ones?

A chance peak is a cluster of sliding windows in which a near-control band has S >= 4.0. There are
68 such peaks (5 of them with S >= 5.8). c_min = 0.967.

- f = fraction of chance peaks with |band sigma_c| >= c_min = **0.426**.
- The median |band sigma_c| of the chance peaks is 0.961.

**A2 = INTERMEDIATE** (0.2 <= f < 0.5). The prior for this outcome was 25 %. Our favoured prediction,
AS ALFVENIC AS PUBLISHED (55 %), **did not come true**, although f is close to its 0.5 line. Among
chance peaks with S >= 5.8, f = 0.40.

## Not graded: how often 5.8-level peaks occur within 20 solar radii

These are cluster counts of sliding windows in E17-E27 at r <= 20 solar radii, with the published
windows excluded. They are a small out-of-sample companion to H1 of PREREG-v0.1 and are **not a graded
test**.

| Length | eligible windows | paper band, S >= 5.8 | 8 control bands, S >= 5.8 (sum; mean) | paper band, S >= 4 | 8 controls, S >= 4 (mean) |
|---|---|---|---|---|---|
| 135 min | 3,374 | 2 | 8; 1.0 | 17 | 13.75 |
| 70 min | 3,394 | 0 | 5; 0.6 | 15 | 8.25 |

At the paper's 5.8 threshold, the control bands produce peaks about as often as the paper band.
At S >= 4 the paper band has somewhat more clusters than the average control band, especially at
70 minutes. This is not graded and is not a test of the claim. The graded population test is H1 of
`psp5min-v1`, which is still queued.

## What it means

The two published wave trains are highly Alfvenic (R1). That is ordinary inside 20 solar radii: about
80 % of ordinary fluctuations at nearby frequencies are less Alfvenic, and about 43 % of chance peaks
of the same statistic are at least as Alfvenic. **Alfvenicity is consistent with the p-mode reading
but does not single the published events out** (A1 NOT DISTINCTIVE). It is not a separate line of
evidence for a solar 5-minute origin. It also does not show that p-modes are absent.

## Limits

These are the limits fixed in PREREG-ALF-v0.1. We used one density per window (the SPAN-i median,
not QTN). The band sigma_c comes from one or two FFT bins. Sliding windows are not tuned to events.
No quality-flag filtering was applied.
