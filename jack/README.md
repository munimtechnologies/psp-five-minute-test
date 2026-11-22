# jack/ - encounter-clustered robustness check for H2c (psp5jack 0.1.0)

Part of [psp-five-minute-test](../README.md), MV2 of Munim, Inc. (mv2@munimtech.com).
Made with AI agents working in MV2's name; no human expert has checked it yet.
Code Apache-2.0, text CC BY 4.0.

`psp5jack.py` reads `windows.csv` from job `psp5curv-v1` and recomputes H2c's rank Z
(3.2 mHz band among nine bands, tiled 135-minute windows) with

- a cluster-robust (CR1) variance, one cluster per Parker Solar Probe encounter, because
  adjacent windows in one solar-wind stream need not be independent;
- the plain Z with each encounter left out;
- the lag-1 rank correlation of time-adjacent windows (reported only).

The grade H2c-R (ROBUST / NOT ROBUST / NOT GRADABLE) and the wording each outcome allows are fixed
in [PREREG-amendment-0.1.2.md](../PREREG-amendment-0.1.2.md), written before any data were read.

    python -m unittest -v test_psp5jack
    python psp5jack.py --windows windows.csv --out jack.json

Standard-library Python 3.9+ only.
