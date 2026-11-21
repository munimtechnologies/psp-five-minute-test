# Pre-registration amendment 0.1.1: calibrating the frequency test and the paper's "sigma"

MV2 of Munim, Inc. (mv2@munimtech.com). Written 2026-11-21 (about 23:30 UTC), after a second
review of `psp5min/p5min.py` 0.1.0 and **before any output of job `psp5min-v1` or `psp5alf-v1`
existed** (both were still waiting for compute; no magnetometer or plasma data had been read by
this project). Made with AI agents working in MV2's name; no human expert has checked it yet.
Text CC BY 4.0, code Apache-2.0.

This amendment changes how H2 is graded, adds one guard to H1, and adds one new graded test (H4).
Nothing else in `PREREG-v0.1.md` changes: no threshold, band, window length or sample.
The new code is `curv/psp5curv.py` 0.1.0 (job `psp5curv-v1`). It imports `p5min.py` unchanged
at commit a025a506d2092e982e82e87c2c27e4f2a477fa98.

## 1. Why H2 needs a calibrated null

H2 ranks the frequency-local ratio L of the 3.2 mHz band among nine bands and grades
Z = (mean rank - 5) / sqrt(6.667 / n), which assumes the target's rank is uniform on 1..9 when
there is no 3 mHz signal. Two things break that assumption:

- **Different flank counts.** L divides by the median of P*f over bins 2..w away, with
  w = max(6, ceil(0.3k)). In 135-minute windows that is 10 flank bins at 1.3 and 2.1 mHz,
  14-16 at the 3.2 mHz target and 22-48 at 4.8-9.8 mHz. A median of fewer bins scatters more and
  is biased differently, so L does not have the same null distribution in every band. Our rough
  asymptotic estimate is a shift of order +0.8 in Z for 3,000-4,000 windows. The sign is not
  certain, because the small-sample bias of the median works the other way.
- **Spectral knees.** Near the Sun the break between the 1/f range and the inertial range
  lies at a few mHz. A local ratio at a concave knee comes out about 2 % high (for a slope
  change of 0.6 in log-log). If every window had a sharp knee exactly at 3.2 mHz, that alone
  would move Z by about +2.

A Z between 2 and 3 from psp5min could therefore come from the method, not the Sun.

**Change.** H2 is graded on **H2c** from `psp5curv`:
L2 = y_k / exp(a), where a is the intercept of a least-squares quadratic fit of log y over
the offsets d = +-2, ..., +-6 bins around k, maximised over the band's bins exactly as L is.
Every band uses the identical fit design. With a locally smooth (quadratic in log) spectrum and
identically distributed ordinates, L2 then has the same null distribution in every band, and
the target's rank is uniform. The windows are the tiled 135-minute windows (all encounters,
published intervals excluded) under the same validity rule as psp5min, and Z uses the same
formula and **the same thresholds**:
SUPPORT if Z >= 3, NO-EXCESS if Z < 2, INDETERMINATE between.

**SIMGATE (fixed now).** psp5curv simulates, in the frequency domain, the observed number of
windows 200 times with no 3 mHz signal:
- SIM-F: flat P*f, three equal components.
- SIM-K: P ~ f^-1 (1 + (f/f_b)^2)^-0.3, with f_b log-uniform in 1-10 mHz per window and
  component powers 0.2 : 1 : 1.
- SIM-K0 (reported only): a sharp knee exactly at 3.2 mHz in every window, slopes -1 to -1.6.
  This is the worst case.

H2c is gradable only if the mean Z of L2 lies within +-0.5 in both SIM-F and SIM-K, and if
n >= 2,000 windows. Otherwise H2c is NOT GRADABLE. psp5min's own H2 Z (from L) is then
reported as **indicative only**, next to the mean Z that L gives in the same simulations
(that mean is its built-in bias). If `psp5curv-v1` does not return, H2 is NOT GRADABLE.

## 2. The paper's "sigma" has a ceiling (Samuelson's inequality)

The paper's statistic is S = (peak - mean) / std over the n_b bins of the central-50 %-energy
range. When the peak bin is one of those n_b bins (likely at 3.1-3.2 mHz, which lies in the
energy-containing range near the Sun; psp5curv checks it), S can never exceed sqrt(n_b - 1),
whatever the data. That is Samuelson's inequality for the population standard deviation.
So S = 5.8 requires n_b >= 35 bins and S = 6.2 requires n_b >= 40. Conversely, a single bin
of any height can give at most sqrt(n_b - 1). S is therefore a bounded z-score, not a
Gaussian significance, and "~6 sigma" cannot be read as a probability of order 1e-9.
psp5curv reports n_b, the ceiling sqrt(n_b - 1) and S for both published windows.

## 3. New graded test H4: what local false-alarm probability does each published peak have?

For each published window, psp5curv fits a smooth spectrum to each magnetic component: a
weighted (1/max(k, 8)) least-squares cubic in log f of log P, over f <= 30 mHz, with
**2.9-3.5 mHz left out of the fit**, plus Euler's gamma (the bias of a log periodogram
ordinate). It then draws Gaussian-noise surrogate spectra around that smooth fit (each
ordinate = smooth x Exp(1) per component; bins above 30 mHz keep their observed values):
20,000 for the 135-minute E22 window and 40,000 for the 70-minute E24 window. On each surrogate
it computes the paper's statistic exactly as psp5min does, at the highest P*f bin in
2.5-4.0 mHz (S_fpk), as the authors picked their peak. The local false-alarm probability is
p = (hits + 1) / (N + 1), where hits counts surrogates with S_fpk >= the observed S_fpk.

- **EXTREME** if p < 1e-4 in both windows. The peaks are far beyond Gaussian noise on a smooth
  spectrum.
- **ORDINARY** if p > 1e-2 in at least one window.
- **MODERATE** otherwise.
- NOT GRADABLE if either window is invalid or has fewer than 10,000 usable surrogates.

This p is per window and ignores how the windows and their edges were chosen, so it favours
the claim. The look-elsewhere side is H1. It assumes Gaussian statistics; real turbulence is
intermittent, which makes large peaks more common, so it also favours the claim.

**Analytic expectation (written now).** For flat P*f and equal transverse power, the
trace ordinate scatters like chi-square with 6 degrees of freedom (sd/mean = 0.577):
P(S >= 5.8) per bin is about e^-13.04 (1 + 13.04 + 85.0) = 2.2e-4, and P(S >= 6.2) is about
1.2e-4. With 4 degrees of freedom (one transverse component dominant) they are about
4e-4 and 2e-4. The search over 2.5-4.0 mHz covers 12 bins at 135 minutes and 6 at 70 minutes,
so p_fpk is about 1e-3 to 5e-3 for E22 and 7e-4 to 1.5e-3 for E24. Those are Gaussian
equivalents of roughly 2.6-3.2 sigma, not 6. A spectral slope across the baseline raises the
std and lowers p.

## 4. Guard added to H1 (fixed now)

H1 = RARE also requires that at least 50 % of the sliding windows in E1-E24 are eligible
(valid and with a usable baseline) at **both** lengths, read from psp5min's printed
`sliding_windows` and `eligible` counts. Otherwise RARE becomes NOT GRADABLE, because a
low count could come from windows dropped for a failed baseline. CHANCE-LEVEL and UNDECIDED
are unchanged.

## 5. Reported, not graded

L and L2 rank Z for 70-minute windows (C1.3 is dropped from L2 there, because its fit would
reach bin 1), for r <= 30 Rs, and for E1-E24 and E25-E27 separately. SIM-K0. A
smooth broken-power-law knee frequency for every tiled 135-minute window: medians by distance,
the fraction of windows with the knee within 0.1 dex of 3.2 mHz, and the mean L2 rank of the
target in those windows against all the others.

## Predictions (written before any result)

- SIMGATE PASS 80 %.
- H2c: NO-EXCESS 70 %, INDETERMINATE 15 %, SUPPORT 15 %.
- H4: MODERATE 55 %, ORDINARY 30 %, EXTREME 15 %.
- psp5min's L-based Z will differ from the L2-based Z by less than 1.5 (70 %).
