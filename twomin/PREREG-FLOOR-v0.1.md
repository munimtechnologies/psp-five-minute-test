# PREREG-FLOOR v0.1: does the young solar wind's energy-containing period keep falling below 2 minutes after encounter 27?

MV2 of Munim, Inc. (mv2@munimtech.com). Written 2026-11-27, **before any data used here exist in public**:
on 2026-11-27 the CDAWeb inventory of `PSP_FLD_L2_MAG_RTN_4_SA_PER_CYC` ended at 2026-04-30T23:59:59Z
and that of `PSP_SWP_SPI_SF00_L3_MOM` at 2026-08-07, so no magnetic field from encounter 28 (perihelion
about June 2026) or later had been released. Made with AI agents working in MV2's name; no human expert
has checked it. Text CC BY 4.0, code Apache-2.0.

## Why

[RESULTS-2MIN-v0.1](RESULTS-2MIN-v0.1.md) tested the 2-minute energy-containing period of Huang, Velli,
Shi et al. (ApJL 977 L12, 2024, arXiv:2405.15967) on encounters E20-E27. The level held (H1 CONSISTENT,
154.6 s), but the floor did not show up: the pooled ratio was r = 0.701 (INTERMEDIATE, 0.001 above the
DECLINING line). In E20-E27 alone it was 0.676, and the 17 windows with t_adv under 6 h had 98.7 s.
These numbers were seen before this file was written. They motivate the test and set our priors, but they
are **not** part of the sample. The test uses only encounters that were not public when it was written.

Related: Morton, Huang, Schad & Velli, arXiv:2611.34265 (2026-11-25), average 270 sub-Alfvenic
intervals from E10-E24 and find energy per log frequency peaking at 4-8 mHz. They do not test whether the
period levels off with advection time.

## Sample C

Encounters 28 and 29 (numbered by the frozen perihelion finder), each covering +-7 days of perihelion as
in psp2min 0.1.0. **Trigger:** the job is filed once, when both CDAWeb inventories (MAG RTN 4 Sa/cyc and
SPAN-i L3 moments) cover perihelion + 7 days of E29. Before that, nobody looks at E28 or E29 spectra. If
E29 is not covered by 2027-09-30, the trigger becomes E28 alone (with the same minimum counts below).

## Code

psp2min 0.1.0 and its pinned p5min / psp5alf modules, unchanged except for these three things in a new
version (0.2.0):
1. `P.MAG_STOP` is moved to the newest MAG inventory end.
2. `P.EXPECTED_PERI` gains E28 and E29 (and E30 if the extension below is used). The times come from the
   minimum heliocentric distance in PSP_HELIO1HR_POSITION, and a unit test checks them.
3. A `--sample C` option grades only E28+.
Window length (4 h, tiled), the spectrum (detrended, Hann-tapered trace PSD), the f_lo cone-of-influence
cut, f_mid = sqrt(f25 * f75), the cleaning cut, the t_adv = (R - Rsun)/median V_r definition and every
validity rule stay byte-identical. A REVIEW-FLOOR file comparing the diff with this list must be committed
before the job runs. Any other change voids the test.

## Tests (thresholds fixed now)

G(S) = geometric mean of 1/f_mid over valid, cleaned windows in set S, all from sample C only.

**F1, floor or decline (primary).** r_C = G(t_adv <= 7 h) / G(7 h < t_adv <= 10 h), the same split as H2.
- FLOOR if r_C >= 0.85; DECLINING if r_C <= 0.70; INTERMEDIATE otherwise.
- Needs >= 15 windows in each range, else NOT GRADABLE.
- Reported with an encounter-block bootstrap 95 % interval (2000 resamples, seed 20261127). The interval is
  reported only and does not change the grade.

**F2, level of the youngest wind.** G(t_adv <= 6 h) in C.
- SUB-2-MIN if < 120 s; TWO-MIN if 120-180 s; ABOVE if > 180 s.
- Needs >= 8 windows, else NOT GRADABLE.

**F3, slope (secondary).** The least-squares slope of ln(1/f_mid) against ln(t_adv) over C windows with
t_adv <= 10 h, with an encounter-block bootstrap 95 % interval.
- STILL FALLING if the whole interval is > 0.
- FLAT if the interval contains 0 and lies within [-0.3, 0.3].
- UNRESOLVED otherwise.

**Gate G0 (pipeline).** Both encounters must have >= 50 % usable windows, the perihelion numbering check
must pass, and K0 (pure f^-5/3 noise gives 600-900 s) must pass again. If G0 fails, nothing is graded.

**Extension, once only.** If F1 or F2 is NOT GRADABLE on E28+E29 for lack of windows, the job is re-filed
once, when E30 is covered, on E28-E30 with the same rules. If it is still short of windows, the result is
NOT GRADABLE.

## Reported, not graded

- G in the 12-20 h bin for C against B (E20-E27), to show any overall shift with the solar cycle. Solar
  activity is falling from the cycle-25 maximum, and a shift there would mean the young-wind numbers
  cannot be read alone.
- Per-encounter G(t_adv <= 10 h), and binned G on the RESULTS-2MIN bins.
- The arithmetic-middle f_mid.

## Our priors (written before any data exist)

- F1: DECLINING 50 %, INTERMEDIATE 30 %, FLOOR 20 %.
- F2: SUB-2-MIN 50 %, TWO-MIN 40 %, ABOVE 10 %.
- F3: STILL FALLING 45 %, FLAT 20 %, UNRESOLVED 35 %.

## What each outcome would mean

- **DECLINING with SUB-2-MIN:** the 2-minute period is a value passed on the way in, not a property
  fixed at the source. The energy-containing period keeps shortening at least down to t_adv of a few
  hours.
- **FLOOR with TWO-MIN:** the first out-of-sample evidence that the period levels off near 2 minutes,
  which supports a source-set period such as leakage of convectively driven waves.
- **Anything else:** reported as it falls. No threshold, sample or definition may change after any
  E28+ data are read.

Limits: tiled 4-hour Hann-FFT windows, not the paper's wavelets; spacecraft-frame frequencies; t_adv
from SPAN-i median speed without quality flags; at the youngest ages the sample is small.
