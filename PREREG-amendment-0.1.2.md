# Pre-registration amendment 0.1.2: is H2c robust to correlated windows and to one encounter?

MV2 of Munim, Inc. (mv2@munimtech.com). Written 2026-11-22 (about 01:30 UTC), after a third review
of the queued code and **before any output of job `psp5min-v1`, `psp5alf-v1` or `psp5curv-v1`
existed** (all three were still waiting for compute; no magnetometer or plasma data had been read
by this project). Made with AI agents working in MV2's name; no human expert has checked it yet.
Text CC BY 4.0, code Apache-2.0.

This amendment adds one qualifier (H2c-R) to H2c and fixes, now, how a positive H2c may be
described. Nothing in `PREREG-v0.1.md` or amendment 0.1.1 changes: no threshold, band, window,
sample, grade or prediction. `psp5min` 0.1.0, `psp5alf` 0.1.0 and `psp5curv` 0.1.0 stay frozen.
The new code is `jack/psp5jack.py` 0.1.0 (standard-library Python, 8 unit tests).

## 1. Two assumptions H2c still makes

H2c (amendment 0.1.1) grades Z = (mean rank - 5) / sqrt(6.667 / n) of the 3.2 mHz band's L2 rank
among nine bands, over all tiled 135-minute windows (n is expected to be 3,000-4,000).

- **Independent windows.** The formula treats every window as an independent draw. Adjacent
  windows of one encounter lie in the same solar-wind stream, so their spectra, and therefore the
  target's rank, can be correlated. With n = 3,000 the standard error of the mean rank is only
  0.047 (out of 9), so a mean shift of 0.14 already gives Z = 3. Positive correlation inflates Z.
- **Identically shaped scatter in every band.** L2's fit design is identical in every band, but
  the null is only identical if the ordinates scatter the same way at every frequency. The scatter
  of a periodogram ordinate depends on how the power is shared among the three field components
  (anisotropy) and on intermittency, both of which can change with frequency. SIM-F and SIM-K use
  one fixed component ratio and Gaussian statistics, so they cannot detect this.

Our second review (run 186) found the first problem's cousin in psp5min's L. This one is smaller
and only matters if H2c comes back positive, so we fix the rule now, before data, rather than
change any frozen code.

## 2. New qualifier H2c-R (graded with psp5jack)

After `psp5curv-v1` returns, its `windows.csv` is committed to this repository unchanged, and
job `psp5jack-v1` downloads it at that commit's full 40-character sha (sha256 printed) and runs
`psp5jack.py`. It uses the 135-minute rows with a finite `rank_L2`, and computes:

- **Z_cl**: the same mean rank with the cluster-robust (CR1) variance, one cluster per encounter:
  V = E/(E-1) * sum_e (sum_i (r_ei - rbar))^2 / n^2, Z_cl = (rbar - 5) / sqrt(V).
  phi = V / (6.667 / n) is the variance inflation (1 for independent windows).
- **Leave-one-encounter-out**: the plain Z recomputed with each encounter removed; the minimum
  and maximum and which encounter gives them.
- Reported only: per-encounter mean ranks and Z; the lag-1 correlation of ranks between
  time-adjacent windows of the same encounter; Z_cl for E1-E24 and E25-E27 separately.

**Grade (fixed now):**
- **ROBUST** if Z_cl >= 3.0 **and** the minimum leave-one-encounter-out Z >= 3.0.
- **NOT ROBUST** otherwise.
- **NOT GRADABLE** if fewer than 2,000 windows or fewer than 10 encounters have a rank, or if
  `psp5curv-v1` produced no `windows.csv`.

## 3. How H2c may be described (fixed now)

| H2c (0.1.1) | H2c-R | wording we will use |
|---|---|---|
| SUPPORT | ROBUST | "a 3.2 mHz excess across all windows that survives correlated windows and the removal of any one encounter; not yet controlled for frequency-dependent scatter" |
| SUPPORT | NOT ROBUST / NOT GRADABLE | "SUPPORT under the independent-window formula only; indicative, not a result" |
| INDETERMINATE | any | "indeterminate" (H2c-R reported) |
| NO-EXCESS | any | "no excess" (H2c-R reported only) |

**Commitment.** If H2c is SUPPORT, then before any article or message calls it evidence for
5-minute oscillations across the encounters, a band-trend control will be pre-registered and run:
L2 ranks of all nine bands in every window, and the target's mean rank compared with a
straight-line interpolation in log f of the eight controls' mean ranks, with encounter-clustered
errors. Until then the third column above applies.

This does not touch the main question about the two published events (H1 look-elsewhere, H4
local false-alarm probability). The paper itself claims oscillations only in two sub-Alfvenic
intervals, not across all windows, so H2c, positive or negative, is a secondary test.

## 4. Predictions (written before any result)

- phi > 1 (windows positively correlated within encounters): 70 %.
- |Z_cl| < |plain Z|: 75 %.
- Lag-1 rank correlation between adjacent windows > 0.05: 55 %.
- H2c-R, given H2c = SUPPORT: ROBUST 40 %, NOT ROBUST 60 %.
- H2c-R ROBUST unconditionally: 8 % (H2c SUPPORT was given 15 % in 0.1.1, and Z_cl >= 3
  requires the plain Z to be at least about 3 * sqrt(phi)).
