# twomin: out-of-sample test of Parker Solar Probe's "2-minute dominance"

Part of MV2's open, pre-registered checks of the Parker Solar Probe near-Sun oscillation claims
(MV2 of Munim, Inc., mv2@munimtech.com). Made with AI agents working in MV2's name; no human expert
has checked it yet. Code Apache-2.0, text CC BY 4.0.

- Claim: Huang, Velli, Shi et al., ApJL 977 L12 (2024), arXiv:2405.15967 - in E1-E19 the
  energy-containing period 1/f_mid of 4-hour magnetic spectra falls with the wind's advection time
  to about 2-3 minutes near the Alfven surface, read as the Alfven-wave launch period.
- Test: the same quantity in E20-E27 (2024-06 to 2026-03, down to 9.86 Rsun), which the claim
  never used, plus a reproduction on E1-E19. Rules, thresholds and predictions:
  [`PREREG-2MIN-v0.1.md`](PREREG-2MIN-v0.1.md), fixed before the job ran.
- Code: `psp2min.py` 0.1.0 (imports `psp5min/p5min.py` and `alfven/psp5alf.py` 0.1.0 unchanged),
  `test_psp2min.py` (synthetic data only). Job: `jobs/psp2min-v1.json`.

Run: `python psp2min.py --data-dir DATA --out-dir OUT` (public NASA CDAWeb data, no login).
Lines `P2MIN_GRADES`, `P2MIN_H1`, `P2MIN_H2`, `P2MIN_R1` carry every number the grades need.
