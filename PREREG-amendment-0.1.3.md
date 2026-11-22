# Pre-registration amendment 0.1.3: how much does "~6 sigma" depend on unstated choices?

MV2 of Munim, Inc. (mv2@munimtech.com). Written 2026-11-22, before any output of psp5min-v1,
psp5alf-v1 or psp5curv-v1 existed and before any magnetometer data were read by this project.
Made with AI agents working in MV2's name; no human expert has checked it yet. Licence: text
CC BY 4.0, code Apache-2.0. Earlier documents: PREREG-v0.1, PREREG-ALF-v0.1, amendments 0.1.1
and 0.1.2. Nothing in them changes.

## Why

The paper (Huang, Velli et al., ApJL 999 L4, 2026, arXiv:2511.10906) defines its significance as
(peak of P(f)*f - mean) / std over "the frequency range containing the central 50% of the total
fluctuation energy". Re-reading the method on 2026-11-22, the text does **not** state: a taper,
any detrending of the magnetic field, the upper frequency limit of the "total" energy, whether the
peak bin is inside the baseline, the std convention, or how exactly the two intervals' edges were
chosen. psp5min's grading gate R0 (reproduce 5.8 and 6.2 within +-1.5) tests one reading only, and
PREREG-v0.1 left the case "R0 fails" to a later amendment. Choosing a reading **after** seeing which
one matches would be a forking path, so the full list of readings and the rule that picks one are
fixed here, before any data are read.

## What we run (job `psp5ladder-v1`, code `ladder/psp5ladder.py` 0.1.0)

Only the two published intervals (E22 2024-12-24 09:50-12:05, E24 2025-06-19 08:15-09:25 UTC),
from the same public CDAWeb data (`PSP_FLD_L2_MAG_RTN_4_SA_PER_CYC`), binned to 0.5 s with
psp5min 0.1.0's own code (p5min.py at commit a025a506, imported unchanged), same validity rule.
The statistic S is taken at the largest P*f bin in 2.5-4.0 mHz, under seven readings, in this
fixed order:

| | reading |
|---|---|
| V0 | psp5min's: mean removed, no taper, energy = sum of P over all bins to Nyquist, edges where the cumulative first reaches 25 % and 75 %, peak bin inside the baseline, population std |
| V1 | V0 with sample std (ddof 1) |
| V2 | V0 with the peak bin and its two neighbours left out of the mean and std |
| V3 | V0 with a straight line removed from each component |
| V4 | V0 with a Hann taper |
| V5 | V0 with the total energy counted only up to 0.1 Hz |
| V6 | cumulative of P*f over linear bins (a less defensible reading of "energy"; reported, not in RL1) |

Also: V0 on each interval with one edge moved by -10, -5, +5 or +10 minutes (8 windows per event).

## Tests (fixed now)

- **RL0 gate**: both published intervals valid with a usable V0 baseline. FAIL -> RL1-RL3 NOT GRADABLE.
- **RL1 reading spread** (V0-V5): spread = max S - min S per event, at least 5 finite values per
  event (else NOT GRADABLE). **FRAGILE** if the spread is >= 2.0 in either event; **STABLE** if
  < 1.0 in both; **MODERATE** otherwise.
- **RL2 reproduction**: a reading passes if, for both events, the peak lies in 3.00-3.35 mHz and S
  is within +-1.5 of 5.8 / 6.2 (psp5min's R0 rule). **REPRODUCED** if V0 passes;
  **REPRODUCED-OTHER-READING** if only others pass; **NOT REPRODUCED** if none passes.
  **Selection rule:** the selected reading is the first passing one in the order V0..V6; if none
  passes, the reading with the smallest worst-event |S - published| (ties by list order).
- **RL3 interval edges**: drop = S(V0) - min S over the 8 moved-edge windows, at least 6 finite per
  event. **EDGE-SENSITIVE** if the drop is >= 2.0 in either event; **ROBUST** if < 1.0 in both;
  **PARTLY** otherwise.
- Reported, not graded: for E22 under V0, the arithmetic midpoint, geometric midpoint and the
  50 %-energy frequency of the baseline band, against the paper's 1/f_mid = 105 s (9.52 mHz),
  each flagged if within 25 %; the Samuelson ceiling sqrt(n_b - 1) per reading.

## What follows from each outcome (fixed now)

- If psp5min-v1's R0 **passes**, nothing changes: H1-H4 stand as graded.
- If R0 **fails** and RL2 selects a reading other than V0 that passes, H1 is re-run once under
  that reading only (a new module and job, thresholds unchanged), and both versions are reported.
- If R0 fails and RL2 is NOT REPRODUCED, no further reading is tried: H1 is reported "with our
  reading of their statistic", as PREREG-v0.1 already says.
- FRAGILE means the published "~6 sigma" is not a stable property of the data: it moves by two
  units or more under choices the paper leaves open. EDGE-SENSITIVE means a 5-10 minute move of an
  edge takes away two units or more, which is what a peak picked from many possible intervals
  would show. **Neither proves the peaks are chance**: a short real wave train is also
  edge-sensitive. They measure how much weight the single number "~6 sigma" can carry. H1 and
  H4 remain the calibrated tests.

## Predictions (written before the run)

- RL0 PASS 85 %.
- RL1: FRAGILE 30 %, MODERATE 40 %, STABLE 30 %. (A taper and detrending mostly change the lowest
  bins, which move the 25 % edge; if the baseline runs from about 0.7 to 20 mHz it holds about
  160 bins at 135 min and 80 at 70 min, so the mean and std
  should move little. The peak bin itself is shared by all readings.)
- RL2: REPRODUCED 45 % (same as R0), REPRODUCED-OTHER-READING 25 %, NOT REPRODUCED 30 %.
- RL3: EDGE-SENSITIVE 35 %, PARTLY 40 %, ROBUST 25 %.
- Midpoint: arithmetic within 25 % of 9.52 mHz 55 %; geometric 20 %; 50 %-energy frequency 25 %.
  (For a P*f plateau up to ~20 mHz and a steeper fall above it, the arithmetic midpoint of the
  25-75 % band lands near 10 mHz, the geometric one near 3-4 mHz.)

## Disclosure

This job reads the two published intervals, which psp5min-v1 also scores for R0. If it returns
first, those S values are seen before psp5min's results. psp5min's code, thresholds and grading
are frozen (commit a025a506) and mechanical, so nothing there can be adjusted. No other window
of any encounter is read by this job.
