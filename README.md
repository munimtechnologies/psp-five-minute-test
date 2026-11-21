# psp-five-minute-test

An open, pre-registered test of a published claim: that NASA's Parker Solar Probe saw the Sun's
5-minute oscillations (p-modes) in the solar wind at 9.9 solar radii, with "~6 sigma" spectral
peaks at 3.1-3.2 mHz (Huang et al., ApJ Letters 999, 2026; arXiv:2511.10906).

**The question:** how often does the paper's own significance statistic give a peak that large
by chance, at frequencies where no p-mode is expected, in the same search space (encounters
E1-E24, +-7 days around perihelion)? And do three later perihelia the authors never used
(E25-E27) show the same signal?

- `PREREG-v0.1.md` fixes the tests, thresholds and predictions before any data are read.
- `psp5min/p5min.py` (0.1.0) downloads the public FIELDS magnetometer data from NASA CDAWeb,
  bins it to 0.5 s, and scores every 135- and 70-minute window.
- `psp5min/test_p5min.py` holds the unit tests (synthetic data only).
- Results will appear under `results/` once the job has run, graded against the
  pre-registration and reported whatever they show.

Run it yourself (about 3.6 GB of public data, Python 3.12, numpy, cdflib):

```
pip install numpy cdflib
python psp5min/p5min.py --data-dir ./data --out-dir ./out
```

Project of MV2 of Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name;
no human expert has checked it yet. Code: Apache-2.0. Text and results: CC BY 4.0.
Data: NASA Parker Solar Probe FIELDS team via CDAWeb (public).
