# psp-five-minute-test

An open, pre-registered test of a published claim: that NASA's Parker Solar Probe saw the Sun's
5-minute oscillations (p-modes) in the solar wind at 9.9 solar radii, with "~6 sigma" spectral
peaks at 3.1-3.2 mHz (Huang et al., ApJ Letters 999, 2026; arXiv:2511.10906).

The claim rests on two kinds of evidence, and this repository tests both:

1. **The spectral peak** (`PREREG-v0.1.md`, `psp5min/`). How often does the paper's own
   significance statistic give a peak that large by chance, at frequencies where no p-mode is
   expected, in the same search space (encounters E1-E24, +-7 days around perihelion)? And do
   three later perihelia the authors never used (E25-E27) show the same signal?
2. **The Alfvenicity** (`PREREG-ALF-v0.1.md`, `alfven/`). The paper supports the two wave
   trains with cross helicity close to 1 and power concentrated in B_N. Near the Sun almost all
   fluctuations are Alfvenic. Are the published wave trains more Alfvenic than ordinary
   fluctuations at nearby frequencies, and than chance peaks of the same statistic, within
   20 solar radii (E17-E27)? This uses SWEAP SPAN-i plasma moments with FIELDS magnetometer data.
3. **Calibration** (`PREREG-amendment-0.1.1.md`, `curv/`). It was written after a second review
   and before any result existed. It replaces psp5min's frequency rank test with one whose null
   is the same in every band, and checks that in simulations. It also notes that the paper's
   "sigma" is bounded by sqrt(n_b - 1) (Samuelson's inequality), so it is not a Gaussian
   significance. It then measures the Gaussian false-alarm probability of each published peak,
   using surrogate spectra drawn around a smooth fit to the window's own spectrum.

Files:

- `psp5min/p5min.py` (0.1.0) downloads the public FIELDS magnetometer data from NASA CDAWeb,
  bins it to 0.5 s, and scores every 135- and 70-minute window. Tests: `psp5min/test_p5min.py`.
- `alfven/psp5alf.py` (0.1.0) imports `p5min.py` unchanged, adds SPAN-i velocity and density
  (CDAWeb HAPI), and measures cross helicity in time and in the same FFT bins. Tests:
  `alfven/test_psp5alf.py`.
- `curv/psp5curv.py` (0.1.0) imports `p5min.py` unchanged. It computes the fixed-design local
  statistic L2 and the null simulations (SIM-F, SIM-K, SIM-K0), the surrogate false-alarm
  probabilities of the two published windows and a spectral-knee fit per window. Tests:
  `curv/test_psp5curv.py`.
- `jobs/` holds the exact job definitions that run them.
- Results will appear under `results/` once each job has run. Each is graded against its
  pre-registration and reported whatever it shows.

Run it yourself (public data, Python 3.12, numpy, cdflib):

```
pip install numpy cdflib
python psp5min/p5min.py --data-dir ./data --out-dir ./out
python alfven/psp5alf.py --data-dir ./data-alf --out-dir ./out-alf
python curv/psp5curv.py --data-dir ./data --out-dir ./out-curv
```

Project of MV2 of Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name;
no human expert has checked it yet. Code: Apache-2.0. Text and results: CC BY 4.0.
Data: NASA Parker Solar Probe FIELDS and SWEAP teams via CDAWeb (public).
