# RESULTS-LADDER v0.1: how much does the published "~6 sigma" depend on unstated method choices?

MV2 of Munim, Inc. (mv2@munimtech.com). Graded 2026-11-25 with AI agents working in MV2's name;
no solar physicist has checked it yet. Text CC BY 4.0, code Apache-2.0.

Grades job `psp5ladder-v1` against `PREREG-amendment-0.1.3.md` (written 2026-11-22, before any
data were read). Paper tested: Z. Huang, M. Velli et al., ApJL 999 L4 (2026), arXiv:2511.10906.
Outputs copied unchanged to `results/psp5ladder-v1/`.

## The run

- Dell, 2026-11-25 19:46-20:04 UTC, 1 CPU, 1 GB limit, status ok, every step exit 0.
- `p5min.py` fetched at pinned commit a025a506; its sha256 in the job (`310727b8...`) matches the pin.
  `psp5ladder.py` sha256 `a8fd37a4...`, `PREREG-amendment-0.1.3.md` sha256 `d9b0f186...`.
- Unit tests: 17 run, 16 ok, 1 skipped (`test_samuelson_ceiling`: "peak outside baseline in this
  draw", a guard in the test, not a failure).
- Data: one CDAWeb file per event (`psp_fld_l2_mag_rtn_4_sa_per_cyc_20241224_v02.cdf`,
  `..._20250619_v02.cdf`), coverage of the 0.5 s grid 1.000 in both intervals; time check against
  cdflib 0 and -2.4e-7 s.

## Scores S at the largest P*f bin in 2.5-4.0 mHz

| reading | E22 S (published 5.8) | E22 peak mHz | E24 S (published 6.2) | E24 peak mHz |
|---|---|---|---|---|
| V0 psp5min reading | 5.315 | 3.210 | 6.148 | 3.095 |
| V1 sample std | 5.301 | 3.210 | 6.115 | 3.095 |
| V2 peak +-1 bin out of baseline | 5.922 | 3.210 | 8.761 | 3.095 |
| V3 linear detrend | 5.404 | 3.210 | 6.174 | 3.095 |
| V4 Hann taper | 6.004 | 3.210 | 5.784 | 3.095 |
| V5 energy only up to 0.1 Hz | 4.823 | 3.210 | 4.700 | 3.095 |
| V6 P*f cumulative (not in RL1) | 27.611 | 3.210 | 31.239 | 3.095 |

Every reading puts the peak in the same bin. Baseline bands under V0: E22 3.09-27.28 mHz (197 bins),
E24 2.62-24.76 mHz (94 bins).

## Grades (mechanical, as the amendment fixed them)

| test | number | grade | our prior |
|---|---|---|---|
| RL0 gate | both intervals valid, V0 baselines usable | **PASS** | PASS 85 % |
| RL1 reading spread (V0-V5) | E22 1.18 (4.82-6.00); E24 **4.06** (4.70-8.76) | **FRAGILE** (>= 2.0 in either event) | FRAGILE 30 % |
| RL2 reproduction | V0 passes both (worst abs dS 0.485); V0, V1, V3, V4, V5 pass; V2 fails E24 (+2.56); V6 fails | **REPRODUCED** (selected V0) | REPRODUCED 45 % |
| RL3 interval edges | drop E22 0.69 (moved-edge range 4.62-6.73); E24 0.88 (5.26-6.25) | **ROBUST** (< 1.0 in both) | ROBUST 25 % |

Reported, not graded (E22, V0), against the paper's 1/f_mid = 105 s (9.52 mHz), flagged if within 25 %:
arithmetic midpoint of the baseline band 15.19 mHz (not within; we gave 55 % to "within"), geometric
midpoint 9.18 mHz (within; prior 20 %), 50 %-energy frequency 8.02 mHz (within; prior 25 %).

## What follows (as fixed in the amendment)

- psp5min-v1's R0 passed, so **nothing changes**: H1 (CHANCE-LEVEL) stands as graded in
  `psp5min/RESULTS-v0.1.md`; the binding H2c and H4 still wait for psp5curv-v1.
- **FRAGILE**: the single number "~6 sigma" is not a stable property of the data. For the E24 event
  it ranges from 4.70 (total energy counted only to 0.1 Hz) to 8.76 (peak left out of its own
  baseline) under choices the paper does not state. Information only, post hoc: without V2 the E24
  spread would be 1.47.
- **ROBUST**: moving either edge of either published interval by 5-10 minutes changes S by less than
  one unit downward (it rises to 6.73 for one E22 window). The published values are not produced
  by a finely tuned choice of interval.
- Neither grade says the peaks are chance. That question is H1's, which uses the same statistic at
  control frequencies; the ladder says the statistic itself is reproduced but its size carries
  method-dependent uncertainty of up to about +-2 units.

## Information only

- With the peak inside a baseline of n_b bins and population std, S cannot exceed sqrt(n_b - 1)
  (Samuelson's bound): 14.0 for E22 and 9.64 for E24 under V0, 6.86 for E24 under V5. A "sigma" of
  this kind is bounded by the size of the baseline and is not a Gaussian tail probability.
- The paper's 105 s mid frequency matches the geometric midpoint or the 50 %-energy frequency of our
  V0 baseline band, not its arithmetic midpoint, which suggests the paper's band is read on a
  logarithmic frequency axis. This does not change any grade.

## Prediction score

4 graded predictions: RL0 (PASS, 85 %) and RL2 (REPRODUCED, 45 %) came out as our favoured outcome;
RL1 (FRAGILE, 30 %) and RL3 (ROBUST, 25 %) did not. Midpoint flags: 0 of 3 as favoured (we gave
"within 25 %" 55 % for the arithmetic midpoint, which missed, and 20 % and 25 % for the geometric
midpoint and the 50 %-energy frequency, which both came true).
