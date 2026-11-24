# Pre-registration: do Parker Solar Probe's switchback patches have a size, a lifetime, or both? (patchclock v0.1)

MV2 of Munim, Inc. - written 2026-11-24, before any data were read for this test.
Licence: text CC BY 4.0; code Apache-2.0. Made with AI agents working in MV2's name;
no human expert has checked it yet.

## 1. Question

Near the Sun, Parker Solar Probe (PSP) sees magnetic "switchbacks" (sudden reversals of the
field) bunched into **patches** lasting hours, separated by quiet radial field. Two published
readings disagree on what a patch is:

* **Fixed spatial imprint.** Patches map to structures on the Sun about the size of
  supergranules, 5-10 deg of Carrington longitude (Bale et al., *Nature* 618, 252, 2023;
  Fargette et al. 2021, *ApJ*, doi:10.3847/1538-4357/ac1112, arXiv:2109.01519; Bale et al. 2021,
  *ApJ*, "magnetic funnels on the supergranulation scale"). Bale et al. 2023 report
  patch (microstream) durations falling from ~10 h to ~2 h toward perihelion, as expected if
  the spacecraft sweeps faster over fixed structures.
* **Lifetime.** Near corotation, when PSP hardly moves over the Sun, patches still last
  4-6 h, so they are temporal bursts, perhaps the "breathing" of the magnetic carpet
  (Shi et al. 2022, *ApJ*, doi:10.3847/1538-4357/ac7c11, arXiv:2206.03807, using E1, E6, E10).

A patch could have both a size theta (deg) and a lifetime tau (h). Nobody has measured the two
jointly across the mission. Encounters E22-E27 (2024-12 to 2026-03, perihelion 9.86 Rsun) sweep
the footpoint at up to ~120 deg/day, ten times faster than the passes either side used, while
every orbit still has corotation turn-arounds. Those two regimes break the degeneracy.

## 2. Data (public, NASA CDAWeb, no login)

* FIELDS `PSP_FLD_L2_MAG_RTN_4_SA_PER_CYC`, binned to 0.5 s by the pinned `psp5min/p5min.py`
  0.1.0 (git blob a025a506...; sha256 310727b8...), +-7 days around every perihelion E1-E27 whose
  window ends before 2026-04-30 (the magnetometer's HAPI stop date seen 2026-11-21).
* SWEAP SPAN-i `PSP_SWP_SPI_SF00_L3_MOM` radial speed via the pinned `alfven/psp5alf.py` 0.1.0.
* `PSP_HELIO1HR_POSITION`: RAD_AU and HGI_LON (inertial heliographic longitude), hourly.
  HG_LON is read only as a frame check (expected median d(HG)/dt - d(HGI)/dt = -14.18 deg/day).

## 3. Method (frozen in `patchclock.py` 0.1.0)

1. **Activity.** For each 10-minute block, the fraction f60 of 0.5 s samples whose field
   is deflected more than 60 deg from the polarity-signed Parker-spiral direction
   (tan psi = Omega r / 350 km/s, Omega = 14.1844 deg/day). Polarity is the sign of the hour's median
   projection; hours with |median B.e / median |B|| < 0.3 (current-sheet crossings, mixed polarity)
   are dropped. Blocks need >= 80 % finite samples and r <= 50 Rsun.
2. **Footpoint.** phi_fp = phi_sc + Omega (r - Rsun) / V, with phi_sc = unwrapped HGI_LON - Omega t
   (a frame rotating at the Carrington rate) and V the 6-hour median SPAN-i V_R (>= 30 samples).
3. **Variogram.** f60 is z-scored per encounter (>= 200 usable blocks). Every pair of blocks
   1-72 blocks apart (10 min - 12 h) in the same encounter and position segment is placed in a
   cell by lag and by |dphi_fp| (25 bins, 0-0.1-...-150-inf deg). Cell mean semivariance
   gamma = <(z_i - z_j)^2 / 2>.
4. **Fit.** gamma = c0 + c1 (1 - exp(-dt/tau - |dphi|/theta)), weighted least squares (weights = pairs,
   cells with >= 30 pairs), on a 71 x 71 log grid (tau 0.25-2000 h, theta 0.1-2000 deg), c1 > 0.
   The grid top means "no measurable decorrelation from that term".
5. **Uncertainty.** 1,000 encounter-block bootstrap replicates (seed 20261124), cells fixed;
   95 % percentile intervals of tau and theta.
6. **Calibration (injection-recovery).** Gaussian fields with covariance
   exp(-|dt|/tau - |dphi|/theta) plus white noise (variance 0.5) are sampled at the real usable
   blocks and real footpoints, then pass through the same z-score, cells and fit. Two seeds each:
   S1 spatial only (theta 5, no lifetime), S2 temporal only (tau 5 h), S3 both (5 h, 5 deg),
   S4 long-lived (36 h, 7.5 deg), S5 (information only) S1 seen through a footpoint that wanders
   from the ballistic one (OU, sigma 2 deg, 6 h).
7. **Model-free check.** Pairs whose two blocks share a class of spacecraft angular rate
   |dphi_sc/dt| (0-5-10-20-40-80-inf deg/day) give an e-folding time per class (the lag where
   gamma has risen 1 - 1/e of the way from its 10-minute value to its 8-12 h mean), compared with
   the fitted model's prediction 1/(1/tau + |omega|/theta). Reported, not graded.

## 4. Gates, hypotheses, grades

* **G0 data**: PASS if >= 15 encounters have >= 200 usable blocks and >= 4 of them are E22 or later.
* **G1 calibration**: PASS if both seeds recover S1 (theta within x1.5 of 5 deg and tau >= 48 h),
  S2 (tau within x1.5 of 5 h and theta >= 40 deg) and S3 (both within x1.5). S4 (theta within x1.5 of 7.5 deg
  and tau >= 18 h, both seeds) decides whether a LONG grade is calibrated.
* **H1 size** (theta point estimate, 95 % interval [lo, hi]):
  SUPERGRANULAR if 2 <= theta <= 15 deg and hi <= 30 deg; SMALLER if theta < 2 and hi < 4;
  NONE if lo > 15; else UNRESOLVED.
* **H2 lifetime** (tau, [lo, hi]): SHORT if tau <= 8 h and hi <= 16 h; LONG if lo >= 24 h
  ("LONG (UNCALIBRATED)" if S4 failed); MIDDLE if 8 < tau < 24 h with lo >= 4 and hi <= 48;
  else UNRESOLVED.
* **Verdict** (only if G0 and G1 PASS): FIXED SUPERGRANULAR IMPRINT (SUPERGRANULAR + LONG);
  TRANSIENT BURSTS (NONE + SHORT); BOTH: SIZED AND SHORT-LIVED (SUPERGRANULAR or SMALLER, with
  SHORT or MIDDLE); else UNDETERMINED. If G0 or G1 fails, H1, H2 and the verdict are NOT
  GRADABLE; the ungated values are still printed and labelled.
* Robustness, reported and never graded: f90 (deflection > 90 deg), spacecraft longitude instead
  of footpoint, r <= 25 Rsun only, E1-E21 only, E22-E27 only.

## 5. Predictions (MV2's priors, written before the run)

G0 PASS 85 %; G1 PASS 65 %. H1: SUPERGRANULAR 50, SMALLER 10, NONE 15, UNRESOLVED 25.
H2: SHORT 35, MIDDLE 30, LONG 20, UNRESOLVED 15. Verdict: FIXED IMPRINT 20, TRANSIENT 15,
BOTH 40, UNDETERMINED 25.

## 6. What would be new

To our knowledge (searches 2026-11-24 of arXiv, ADS-indexed journals and the Parker Four 2025
programme) no study has fitted a patch lifetime and an angular size jointly, with an
injection-recovery calibration, across E1-E27, and no published patch durations or sizes exist
for the 9.86 Rsun passes. If the verdict is BOTH, it reconciles Bale et al. and Shi et al. with
two numbers; if FIXED or TRANSIENT, it decides between them.

## 7. Known limits, stated in advance

* Ballistic mapping ignores coronal non-radial expansion; theta is an angle at the mapped
  source longitude (roughly the source surface), not at the photosphere. Mapping errors can make
  fixed structures look short-lived; S5 measures how much.
* Activity variance may change with distance; the per-encounter z-score equalises encounters
  but not distances within one. The r <= 25 Rsun fit is the check.
* Only longitude separation is used (PSP stays near the ecliptic).
* The 0.5 s binning and 4-samples-per-cycle magnetometer product smooth the shortest switchbacks;
  activity fractions are not comparable in absolute value with other papers.

## 8. Run

Job `patchclock-v1`, data_key `psp-2min-e1-e27` (reuses the cached magnetometer bins and
SPAN-i days of `psp2min-v1`; it is filed only after that job finishes so that the two never write
the cache at the same time). Any change after filing is an amendment, numbered from 0.1.1.
