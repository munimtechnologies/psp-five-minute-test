# Results 2MIN v0.1: does Parker Solar Probe's "2-minute dominance" hold in encounters the claim never used?

MV2 of Munim, Inc. (mv2@munimtech.com). Graded 2026-11-25 against
[`PREREG-2MIN-v0.1.md`](PREREG-2MIN-v0.1.md) (written 2026-11-23, before any f_mid value was computed)
and the pre-run review [`REVIEW-2MIN-v0.1.md`](REVIEW-2MIN-v0.1.md) (no amendment). Made with AI agents
working in MV2's name; no human expert has checked it yet. Text CC BY 4.0, code Apache-2.0. No
threshold, sample or definition was changed after the data were seen.

Job `psp2min-v1` (`twomin/psp2min.py` 0.1.0 with pinned `p5min.py` and `psp5alf.py`; 11/11 unit
tests pass) ran 2026-11-25 15:41-16:48 UTC on public NASA CDAWeb data (FIELDS MAG RTN 4 Sa/cyc,
SWEAP SPAN-i L3 moments, hourly positions), tiled 4-hour windows over +-7 days of E1-E27. Outputs:
[`results/psp2min-v1/`](../results/psp2min-v1/) (`results.json`, `windows.csv`, `run.txt`, logs).

The claim: Huang, Velli, Shi et al., ApJL 977 L12 (2024), arXiv:2405.15967: in E1-E19 the
energy-containing period 1/f_mid falls with the wind's advection time t_adv to about 2-3 minutes in the
youngest wind, with a floor "likely closer to 2 minutes".

## Grades

| Test | Measured | Grade | Our prior |
|---|---|---|---|
| P0 pipeline | numbering correct; 19/19 of E1-E19 and 8/8 of E20-E27 with >= 50 % usable windows | **PASS** | 80 % |
| K0 | pure f^-5/3 noise gives 1/f_mid = 664.2 s (band 600-900 s) | **PASS** | 90 % |
| Cleaning cut | removes 2.4 % of usable E1-E19 windows (paper: 3.7 %; limit 15 %), so it is used | applied | - |
| R1 reproduction (E1-E19) | Spearman(t_adv, 1/f_mid) = 0.548 (> 0.3); geometric mean 1/f_mid for t_adv <= 10 h = 199.6 s from 69 windows (band 90-240 s) | **PASS** | 55 % |
| **H1 level in the youngest new wind (E20-E27, t_adv <= 10 h)** | G = **154.6 s** (2.6 min), 70 windows from 8 encounters, encounter-block bootstrap 95 % 142.4-167.1 s; band 90-180 s | **CONSISTENT** | 55 % |
| H2 floor or decline (A+B pooled) | G(t_adv <= 7 h) = 138.4 s (n = 46) / G(7-10 h) = 197.4 s (n = 93): ratio **0.701**; FLOOR >= 0.85, DECLINING <= 0.70 | **INTERMEDIATE** | 30 % |

H2 missed the DECLINING threshold by 0.001. We grade it INTERMEDIATE as the rule says.

## Reported, not graded

- Ratio in the new encounters alone (B): 0.676 (below the DECLINING line).
- Per-encounter H1 means (s): E20 156.0 (5 windows), E21 170.8 (1), E22 169.4 (5), E23 131.5 (5),
  E24 139.3 (15), E25 140.8 (8), E26 159.8 (16), E27 176.6 (15). Spearman(t_adv, 1/f_mid) in B 0.696.
- Binned geometric means in B (t_adv bin: n, 1/f_mid): 0-6 h: 17, 98.7 s; 6-7 h: 12, 167.8 s;
  7-8 h: 12, 172.5 s; 8-10 h: 29, 185.8 s; 10-12 h: 43, 236.4 s; 12-15 h: 76, 305.2 s;
  15-20 h: 107, 359.9 s; 20-25 h: 95, 413.4 s; 25-30 h: 68, 483.9 s; 30-40 h: 73, 500.1 s;
  40-60 h: 18, 599.4 s.
- In A (E1-E19), the same bins: 0-6 h: 2, 432.9 s; 6-7 h: 15, 149.4 s; 7-8 h: 20, 199.0 s;
  8-10 h: 32, 218.3 s; 10-12 h: 44, 291.1 s; 12-15 h: 114, 334.6 s; 15-20 h: 311, 393.2 s.
- Arithmetic-middle definition (f25 + f75)/2: H1 G = 89.7 s.
- 4.3 % of young B windows have their P(f) f peak at the lowest frequency; geometric mean
  1/f_peak 333.9 s.

## What this means

In the youngest solar wind ever sampled (E20-E27, down to 9.86 Rsun, none of it used for the claim),
the energy-containing period sits at 2.6 minutes (95 % 2.4-2.8), inside the pre-registered 90-180 s band
for "2-minute dominance". This is the first out-of-sample support for the level the paper reported.
Whether 2 minutes is a floor is not settled: the period still falls by 30 % between 7-10 h and
under 7 h of advection time, right at our DECLINING line, and the youngest new bin (under 6 h,
17 windows) sits near 99 s. A floor would show as a flat ratio near 1; we do not see one.

Limits: a Hann-tapered FFT, not the paper's cone-averaged wavelet spectrum; t_adv uses a window
median of SPAN-i radial speed without quality flags; spacecraft-frame frequencies are not
Doppler-corrected; tiled rather than sliding windows.
