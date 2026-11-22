# ladder: does the published "~6 sigma" survive the choices the paper leaves open?

Part of [psp-five-minute-test](https://github.com/munimtechnologies/psp-five-minute-test), MV2 of
Munim, Inc.'s open, pre-registered test of the claim that Parker Solar Probe saw the Sun's
5-minute oscillations (Huang, Velli et al., ApJL 999 L4, 2026, arXiv:2511.10906).

`psp5ladder.py` (0.1.0) scores the paper's statistic on the two published intervals only, under
seven fixed readings of its method (taper, detrending, std convention, peak in or out of the
baseline, upper frequency limit of the "total energy", energy measure) and with each interval edge
moved by up to 10 minutes. Rules, the reading-selection rule for psp5min's reproduction gate, and
predictions are in [PREREG-amendment-0.1.3.md](../PREREG-amendment-0.1.3.md), written before any
data were read. It imports psp5min 0.1.0 (`psp5min/p5min.py`, commit a025a506) unchanged.

Run: `python psp5ladder.py --data-dir DATA --out-dir OUT` (needs numpy and cdflib; public NASA
CDAWeb data, no login). Tests: `python -m unittest -v test_psp5ladder` (synthetic data only).
Job: [jobs/psp5ladder-v1.json](../jobs/psp5ladder-v1.json).

Made with AI agents working in MV2's name; no human expert has checked it yet.
Licence: code Apache-2.0, text CC BY 4.0.
