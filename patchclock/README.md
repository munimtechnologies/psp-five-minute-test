# patchclock: do Parker Solar Probe's switchback patches have a size, a lifetime, or both?

Part of MV2's open Parker Solar Probe tests (MV2 of Munim, Inc., mv2@munimtech.com).
Code Apache-2.0; text CC BY 4.0. Made with AI agents working in MV2's name; no human expert
has checked it yet.

Switchbacks near the Sun come in patches lasting hours. One reading (Bale et al. 2021, 2023;
Fargette et al. 2021) says patches are fixed structures about the size of supergranules
(5-10 deg on the Sun). Another (Shi et al. 2022) says they are bursts with a lifetime of
4-6 h, because they last that long even when the spacecraft corotates with the Sun.
`patchclock.py` measures both at once across encounters E1-E27: a lifetime tau and an angular
size theta, fitted to the space-time variogram of switchback activity, with an
injection-recovery calibration on the real orbit tracks.

* Rules, thresholds and predictions, fixed before any data were read:
  [PREREG-PATCH-v0.1.md](PREREG-PATCH-v0.1.md)
* Code: [patchclock.py](patchclock.py) 0.1.0; tests (synthetic only):
  [test_patchclock.py](test_patchclock.py)
* Reuses, pinned and unchanged: `psp5min/p5min.py` 0.1.0 and `alfven/psp5alf.py` 0.1.0.
* Data: NASA CDAWeb (FIELDS MAG RTN 4 samples/cycle, SWEAP SPAN-i L3 moments, PSP hourly
  positions), public and without login.

Status (2026-11-24): pre-registered; the job runs after the queued `psp2min-v1` has filled the
shared data cache. Results will go to `results/patchclock-v1/`.

Run locally: `python patchclock.py --data-dir data --out-dir out` (downloads several GB).
