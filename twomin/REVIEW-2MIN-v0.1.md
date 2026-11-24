# Pre-run review of job psp2min-v1 against its real inputs

MV2 of Munim, Inc. (mv2@munimtech.com). Written 2026-11-24, while `psp2min-v1` is still queued and
before any f_mid or t_adv value from any Parker Solar Probe window has been computed by this
project. Made with AI agents working in MV2's name; no human expert has checked it yet. Licence:
CC BY 4.0.

Scope: check that the data `twomin/psp2min.py` 0.1.0 will read exist, have the layout the code
assumes, and are dense enough for the rules in `PREREG-2MIN-v0.1.md`, and re-read the claim's
method for any mismatch that would need an amendment before the run. **Result: no defect found;
no amendment; the job runs as filed.** The next amendment number stays 0.2.

## 1. SPAN-i moments (wind speed V_r)

- HAPI `info` for `PSP_SWP_SPI_SF00_L3_MOM` (cdaweb.gsfc.nasa.gov/hapi, read 2026-11-24): startDate
  2018-10-02T04:04:01Z, stopDate 2026-08-07T15:01:46Z, so every encounter E1-E27 is covered.
  Parameter order: Time, CNTS, QUALITY_FLAG (integer, fill 65535), DENS, VEL_INST[3], VEL_SC[3],
  VEL_RTN_SUN[3], ... A request for `QUALITY_FLAG,DENS,VEL_RTN_SUN` therefore returns
  time, QF, DENS, V_R, V_T, V_N. `psp5alf.parse_hapi_spi` keeps columns 1-5 as
  [QF, DENS, V_R, V_T, V_N], and `psp2min` uses column index 2, which is V_R. Correct.
- E1 check (perihelion 2018-11-06): 2018-11-06T00:00-00:03 returned 6 rows, 27.96 s apart, V_R
  354-366 km/s, QUALITY_FLAG 4372. 2018-11-11 (five days after perihelion) is continuous at
  27.96 s from 00:00:31 to at least 20:54. At 28 s cadence every 1-minute slot holds about two
  samples (the rule needs 60 % of slots filled) and every 10-minute block holds about 21 (the
  cleaning statistic needs 5). So the earliest, sparsest encounters can still give usable windows.
- Dropouts inside the archive's coverage intervals exist: 2018-11-12T12:00-12:10 returned HAPI
  status 1201 ("no data") although the CDAS inventory lists one interval 2018-10-31T00:31 to
  2018-11-13T01:35. The code marks such windows unusable (v_ok = 0); this is handled, not a defect.
- CDAS inventory over the +-7 day spans of the first three encounters: E1 2018-10-30 (02:32-22:25)
  and 2018-10-31 to 2018-11-13 continuous; E2 (perihelion 2019-04-04) 2019-03-28 and 03-29 partial,
  2019-03-30 to 2019-04-11 continuous; E3 (2019-09-01) 2019-08-21 to 2019-09-09 continuous. The
  P0 gate (at least 15 of E1-E19 with half their windows usable) is reachable.
- The quality flag changes within a day (4356 to 12548 on 2018-11-11). As pre-registered, flags are
  not used; V_R is a window median within 50-2,000 km/s.

## 2. Magnetometer and access path

The same download and binning code (`p5min.py` 0.1.0, pinned) and the same SPAN-i fetcher
(`psp5alf.py` 0.1.0, pinned) already ran in job `psp5alf-v1` (E17-E27, +-3 days; results in
`alfven/results/`): 7 of 7 FIELDS files per encounter, 0.5 s magnetometer coverage 0.9997-1.0,
SPAN-i 2 s coverage 0.972-0.999, SPAN-i median spacing 3.4952 s (E17-E23) and 1.7476 s
(E24-E27), 14-20 s per encounter, peak memory 0.19 GB. psp2min reads +-7 days for 27 encounters,
about 378 days of data against 66 (under six times the volume), which fits its 13,000 s budget and 2 GB limit by a wide margin.

## 3. Method re-read (arXiv:2405.15967, HTML, re-read 2026-11-24)

| Paper | psp2min-v1 | Status |
|---|---|---|
| E1-E19, +-7 days around each perihelion | same for sample A; E20-E27 as sample B | match |
| 4-h windows, sliding with a 15-min step (26,623 intervals) | 4-h windows, tiled | declared difference (independence) |
| V_r = radial first-order SPAN-ion proton moment | median SPAN-i V_R (Sun frame) | match |
| QTN electron density as a proton-density proxy | not needed (no density enters f_mid or t_adv) | n/a |
| keep 0.005 < dt_adv/<t_adv> < 0.2 (982 removed, 3.7 %) | same bounds, our reading of dt_adv | declared reading |
| trust the spectrum down to where >50 % of wavelet points leave the cone of influence | f_lo = 0.380 mHz from the same rule | declared reading |
| t_adv = (R - Rsun)/V_r; smallest about 5.7 h | same | match |
| Alfvenic subset abs(sigma_c) > 0.75 (9,125 intervals), Fig. 4d | not computed | not tested (see below) |

Disclosure: Fig. 4d of the paper, which keeps only highly Alfvenic intervals, is described by the
authors as a trend that "doesn't seem to stabilize with t_adv". That statement bears on our H2
(floor or continued decline). It is the authors' own published description, not a number from
our data. It was noted in this review, after the pre-registration and before the run. No
threshold, sample or prediction changes because of it.

## 4. Decision

No amendment is needed. `psp2min-v1` stays queued unchanged (it must not be re-filed or bumped for
a wait). Grading follows `PREREG-2MIN-v0.1.md` exactly.
