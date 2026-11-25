# Results of job psp5curv-v1, graded against PREREG-amendment-0.1.1

MV2 of Munim, Inc. (mv2@munimtech.com). Graded 2026-11-25 (about 23:40 UTC) from
`results/psp5curv-v1/results.json` (copied unchanged from the job). Every grade below was
re-derived by hand from that file against the rules in `PREREG-amendment-0.1.1.md`; no threshold,
band, window or sample was moved. Made with AI agents working in MV2's name; no human expert has
checked it yet. Text CC BY 4.0, code Apache-2.0.

## Job

- `psp5curv-v1`, Blade, 2026-11-25 21:02-21:34 UTC, status ok, peak memory 0.18 GB.
- `curv/psp5curv.py` 0.1.0 sha256 `58b6fac1...ce52`; imports `psp5min/p5min.py` at commit
  a025a506 with sha256 `310727b8...05a4` (matches the pin). 15/15 unit tests passed.
- Public FIELDS magnetometer data from NASA CDAWeb for all 27 encounters (14-15 files each).
- 3,931 tiled 135-minute windows and 7,622 tiled 70-minute windows (published intervals excluded).

## Graded

| Test | Rule (fixed 2026-11-21) | Number | Grade | Our prior |
|---|---|---|---|---|
| SIMGATE | mean Z of L2 within +-0.5 in SIM-F and SIM-K, n >= 2,000 | SIM-F -0.024 (sd 1.03), SIM-K +0.091 (sd 0.87), n = 3,931 | **PASS** | 80 % |
| H2c | Z >= 3 SUPPORT, < 2 NO-EXCESS | mean rank 5.0315 of 9, Z = 0.031544 / sqrt(6.667/3931) = **0.766** | **NO-EXCESS** | 70 % |
| H4 | EXTREME if p < 1e-4 in both; ORDINARY if p > 1e-2 in at least one | E22 p = (301+1)/(20000+1) = **0.0151**; E24 p = (114+1)/(40000+1) = **0.00287** | **ORDINARY** | 30 % |

Prediction "psp5min's L-based Z will differ from the L2-based Z by less than 1.5" (70 %):
2.687 - 0.766 = 1.92. **Failed.** The reason is in the simulations below: L carries a built-in
bias of about +2 that L2 removes.

## H2c-R (amendment 0.1.2)

With H2c = NO-EXCESS, amendment 0.1.2 section 3 fixes the wording as "no excess" and makes H2c-R
"reported only". No outcome of H2c-R can change a grade or the wording, so job `psp5jack-v1` is
not filed now, to keep the shared queue short. `results/psp5curv-v1/windows.csv` is committed
unchanged, so anyone can run `jack/psp5jack.py` on it.

## Reported, not graded

**The indicative rank test explained.** psp5min-v1 gave Z = 2.69 from the simple local ratio L
(135-minute windows, all encounters). In 200 simulations with the observed number of windows and
**no** 3 mHz signal, L gives a mean Z of 1.88 (sd 0.94, flat spectra) and 2.05 (sd 0.97, random
knees in 1-10 mHz); a sharp knee at exactly 3.2 mHz in every window gives 8.32. The data's 2.69 is
0.63 above the knee simulation, within one simulation sd. The curvature-corrected L2 gives
-0.02, 0.09 and 0.79 in the same three simulations and 0.77 on the data.

| Sample | L Z | L2 Z |
|---|---|---|
| 135 min, all (n 3,931) | 2.69 | 0.77 |
| 135 min, r <= 30 Rs (n 1,069) | 0.95 | 1.10 |
| 135 min, E1-E24 (n 3,487) | 2.60 | 0.52 |
| 135 min, E25-E27 (n 444) | 0.70 | 0.81 |
| 70 min, all (n 7,622; L2 over 8 bands) | -1.00 | 0.96 |
| 70 min, E25-E27 (n 858) | -2.22 | -1.50 |

**The two published windows.**

| | E22 (2024-12-24 09:50-12:05) | E24 (2025-06-19 08:15-09:25) |
|---|---|---|
| Published score | 5.8 at 3.2 mHz | 6.2 at 3.1 mHz |
| Our S (peak bin) | 5.315 at 3.210 mHz | 6.148 at 3.095 mHz |
| Baseline bins n_b, peak inside | 197, yes | 94, yes |
| Samuelson ceiling sqrt(n_b - 1) | 14.0 | 9.64 |
| Surrogates | 20,000 | 40,000 |
| p, peak searched over 2.5-4.0 mHz (graded) | 0.0151 (z 2.17) | 0.00287 (z 2.76) |
| p, fixed band bin only | 0.00265 (z 2.79) | 0.00045 (z 3.32) |
| Surrogate 99 % / 99.9 % of S (searched) | 5.63 / 7.19 | 5.40 / 6.74 |
| L2 rank of 3.2 mHz among 9 bands | 6 | 8 |

The written analytic expectation (amendment 0.1.1 section 3) was p of about 1e-3 to 5e-3 for E22
and 7e-4 to 1.5e-3 for E24. Both observed p are larger (less significant), as expected when a
spectral slope across the baseline raises the scatter. These p ignore how the windows and edges
were chosen and assume Gaussian statistics; both favour the claim.

**Spectral knees.** A smooth broken power law fitted to every 135-minute window puts the median
knee at 3.02 mHz inside 12 Rs (80 windows), 1.82 mHz at 12-20 Rs (411), 1.15 at 20-30 (578), 1.26
at 30-45 (1,198) and 1.05 beyond 45 Rs (1,664). 9.0 % of windows have the knee within 0.1 dex of
3.2 mHz; 28 % hit the edge of the fit grid. In those 352 knee-near windows the 3.2 mHz band's
mean L2 rank is 4.64, against 5.07 in the other 3,579. So closest to the Sun, where the paper
finds its events, the spectrum typically bends right at the claimed frequency; a simple local
contrast there is biased upward (SIM-K0), the curvature-corrected one is not.

## What this means for the 5-minute claim

All three graded tests of the published "~6 sigma" peaks come out at the chance level: H1
CHANCE-LEVEL (psp5min-v1), H2c NO-EXCESS and H4 ORDINARY (this job). The peaks are reproduced
(R0 PASS, RL2 REPRODUCED), are not sensitive to the exact interval edges (RL3 ROBUST), but their
size depends strongly on unstated method choices (RL1 FRAGILE). None of this tests the Alfvenic
polarization of the wave trains, which the paper also relies on, and none of it shows that
5-minute oscillations are absent from the solar wind.
