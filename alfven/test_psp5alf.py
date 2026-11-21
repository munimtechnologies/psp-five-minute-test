# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. Written by AI agents working for MV2 of Munim, Inc.
"""Unit tests for psp5alf (synthetic data only, no network).
Expected values are computed by hand in the comments."""
import math
import os
import sys
import unittest

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(_HERE), "psp5min"))
sys.path.insert(0, _HERE)

try:
    import numpy as np
    import psp5alf as A
    import p5min as P
    HAVE = True
except Exception:  # noqa: BLE001
    HAVE = False


@unittest.skipUnless(HAVE and os.name == "posix", "needs numpy on POSIX")
class TestPsp5alf(unittest.TestCase):

    def test_va_factor(self):
        # 1 nT, 1 cm^-3: 1e-9 / sqrt(4 pi 1e-7 * 1e6 * 1.67262e-27) = 21812 m/s
        self.assertAlmostEqual(A.va_factor(1.0), 21.812, delta=0.005)
        self.assertAlmostEqual(A.va_factor(100.0), 2.1812, delta=0.0005)

    def test_parse_hapi_spi(self):
        txt = ("#header\n"
               "2024-12-24T11:00:01.627192576Z,4336,4.5544e+03,2.7523e+02,-6.6993e+00,-3.8223e+01\n"
               "2024-12-24T11:00:03.374812928Z,4336,NaN,NaN,NaN,NaN\n"
               "\"status\": {\"code\": 1200}\n"
               "garbage,1,2\n")
        t, a = A.parse_hapi_spi(txt)
        self.assertEqual(len(t), 2)
        self.assertAlmostEqual(float(t[0]), P.to_unix("2024-12-24T11:00:01Z") + 0.627192576,
                               places=5)
        self.assertAlmostEqual(float(a[0, 1]), 4554.4)
        self.assertAlmostEqual(float(a[0, 2]), 275.23)
        self.assertAlmostEqual(float(a[0, 4]), -38.223)
        self.assertTrue(math.isnan(a[1, 2]))

    def test_interp_gap_rule(self):
        ts = np.array([0.0, 1.75, 3.5, 20.0, 21.75])
        grid = np.array([1.0, 3.0, 5.0, 11.0, 21.0, 30.0])
        out = A.interp_to_grid(ts, ts, grid)[:, 0]
        # brackets: 0-1.75 ok, 1.75-3.5 ok, 3.5-20 gap 16.5 > 10, same, 20-21.75 ok, beyond end
        self.assertAlmostEqual(out[0], 1.0)
        self.assertAlmostEqual(out[1], 3.0)
        self.assertTrue(math.isnan(out[2]))
        self.assertTrue(math.isnan(out[3]))
        self.assertAlmostEqual(out[4], 21.0)
        self.assertTrue(math.isnan(out[5]))

    def test_block_mean(self):
        x = np.array([1, 2, 3, 4, 5] + [np.nan] * 7, float)[:, None]
        y = A.block_mean(x, 4)[:, 0]
        self.assertAlmostEqual(y[0], 2.5)
        self.assertAlmostEqual(y[1], 5.0)
        self.assertTrue(math.isnan(y[2]))

    def test_rolling_detrend(self):
        x = np.full((1000, 3), 7.0)
        self.assertLess(float(np.abs(A.rolling_detrend(x)).max()), 1e-9)
        lin = np.arange(1000.0)[:, None] * np.ones((1, 3))
        d = A.rolling_detrend(lin)
        # centred 301-point mean of a line equals the line away from the edges
        self.assertLess(float(np.abs(d[150:850]).max()), 1e-9)
        self.assertGreater(float(np.abs(d[0]).max()), 1.0)

    def test_fill_gaps(self):
        x = np.arange(1000.0)[:, None] * np.ones((1, 3))
        y = x.copy()
        y[100:140] = np.nan          # 40-point run > 30
        self.assertIsNone(A.fill_gaps(y))
        y = x.copy()
        y[100:125] = np.nan          # 25-point run, 2.5 %
        f = A.fill_gaps(y)
        self.assertIsNotNone(f)
        self.assertLess(float(np.abs(f - x).max()), 1e-9)

    def test_sigma_c(self):
        i = np.arange(2000)
        s = np.stack([np.sin(2 * np.pi * 7 * i / 2000), np.zeros(2000), np.zeros(2000)], 1)
        c = np.stack([np.cos(2 * np.pi * 7 * i / 2000), np.zeros(2000), np.zeros(2000)], 1)
        self.assertAlmostEqual(A.sigma_c(s, s), 1.0)
        self.assertAlmostEqual(A.sigma_c(s, -s), -1.0)
        self.assertAlmostEqual(A.sigma_c(s, c), 0.0, places=9)

    def test_band_values(self):
        n = 2100
        i = np.arange(n)
        b = np.zeros((n, 3))
        b[:, 2] = np.sin(2 * np.pi * 13 * i / n)
        v = b.copy()
        v[:, 0] += np.sin(2 * np.pi * 20 * i / n)   # velocity-only power at bin 20
        t = A.band_terms(v, b)
        k13 = A.band_values(t, [13])
        self.assertAlmostEqual(k13["sc"], 1.0, places=9)
        self.assertAlmostEqual(k13["sr"], 0.0, places=9)
        self.assertAlmostEqual(k13["bn_share"], 1.0, places=9)
        k20 = A.band_values(t, [20])
        self.assertAlmostEqual(k20["sc"], 0.0, places=6)
        self.assertAlmostEqual(k20["sr"], 1.0, places=6)
        b2 = np.zeros((n, 3))
        b2[:, 1] = np.sin(2 * np.pi * 13 * i / n)
        b2[:, 2] = np.cos(2 * np.pi * 13 * i / n)
        self.assertAlmostEqual(A.band_values(A.band_terms(b2, b2), [13])["bn_share"], 0.5,
                               places=9)

    def test_alf_window(self):
        rng = np.random.default_rng(3)
        n = 2100
        lay = P.layout(4200)
        b = rng.normal(size=(n, 3)) + np.array([50.0, 0.0, 0.0])
        dens = np.full(n, 100.0)
        v = (b - b.mean(axis=0)) * A.va_factor(100.0) + np.array([300.0, 0.0, 0.0])
        rec = A.alf_window(b, v, dens, lay)
        self.assertAlmostEqual(rec["n_med"], 100.0)
        self.assertAlmostEqual(rec["sc_t"], 1.0, places=9)
        self.assertAlmostEqual(rec["sc_P3.15"], 1.0, places=9)
        self.assertAlmostEqual(rec["sc_C4.8"], 1.0, places=9)
        rec2 = A.alf_window(b, 2 * np.array([300.0, 0, 0]) - v, dens, lay)
        self.assertAlmostEqual(rec2["sc_t"], -1.0, places=9)
        d2 = dens.copy()
        d2[:200] = np.nan            # 9.5 % of densities missing
        self.assertIsNone(A.alf_window(b, v, d2, lay))

    def test_frac_below(self):
        self.assertAlmostEqual(A.frac_below([0.1, 0.2, 0.3, 0.4], 0.35), 0.75)
        self.assertTrue(math.isnan(A.frac_below([0.1, 0.2], float("nan"))))

    def _row(self, enc, D, s, S21=0.0, sc=0.5, tiled=0):
        r = {"enc": enc, "D": D, "s": s, "start": "x", "r_rs": 10.0, "valid": 1, "base_ok": 1,
             "alf_ok": 1, "pub": 0, "sliding": 1, "tiled": tiled, "sc_t": 0.9,
             "S_P3.15": 0.0, "sc_P3.15": 0.5, "bn_P3.15": 0.5}
        for fc in P.CONTROL_CENTRES:
            r[f"S_C{fc}"] = 0.0
            r[f"sc_C{fc}"] = sc
            r[f"bn_C{fc}"] = 0.4
        r["S_C2.1"] = S21
        return r

    def test_summarise_and_grade(self):
        rows = [self._row(22, 70, 0, S21=5.0, sc=-0.85, tiled=1),
                self._row(22, 70, 1200, S21=4.5),
                self._row(22, 70, 2400, S21=3.0),
                self._row(23, 70, 0, S21=4.2, sc=0.1, tiled=1)]
        pubs = [{"enc": 22, "D": 135, "sc_P3.15": 0.9, "sc_t": 0.95},
                {"enc": 24, "D": 70, "sc_P3.15": -0.8, "sc_t": -0.85}]
        res = A.summarise(rows, pubs)
        self.assertAlmostEqual(res["c_min"], 0.8)
        # two peaks: E22 cluster (best S 5.0, |sc| 0.85) and E23 (S 4.2, |sc| 0.1)
        self.assertEqual(res["A2"]["n_peaks"], 2)
        self.assertAlmostEqual(res["A2"]["frac_ge_cmin"], 0.5)
        self.assertEqual(res["A2"]["n_peaks_S_ge_5p8"], 0)
        # E24 (D=70) reference: tiled rows' near controls -> |sc| 0.85 x3 and 0.1 x3
        self.assertEqual(res["A1"]["E24"]["n_ref"], 6)
        self.assertAlmostEqual(res["A1"]["E24"]["q_near"], 0.5)
        g = A.grade(res, pubs)
        self.assertEqual(g["R1"], "PASS")
        self.assertEqual(g["A1"], "NOT GRADABLE")    # n_ref < 30
        self.assertEqual(g["A2"], "NOT GRADABLE")    # fewer than 5 peaks

    def test_grade_rules(self):
        def res(q1, q2, n, frac):
            return {"A1": {"E22": {"q_near": q1, "n_ref": 100},
                           "E24": {"q_near": q2, "n_ref": 100}},
                    "A2": {"n_peaks": n, "frac_ge_cmin": frac}}
        pubs = [{"sc_t": 0.9}, {"sc_t": 0.5}]
        g = A.grade(res(0.5, 0.6, 10, 0.6), pubs)
        self.assertEqual(g["A1"], "NOT DISTINCTIVE")
        self.assertEqual(g["A2"], "AS ALFVENIC AS PUBLISHED")
        self.assertEqual(g["R1"], "FAIL")
        g = A.grade(res(0.995, 0.999, 10, 0.1), pubs)
        self.assertEqual(g["A1"], "DISTINCTIVE")
        self.assertEqual(g["A2"], "LESS ALFVENIC")
        g = A.grade(res(0.5, 0.995, 10, 0.3), pubs)
        self.assertEqual(g["A1"], "MIXED")
        self.assertEqual(g["A2"], "INTERMEDIATE")
        g = A.grade(res(0.5, 0.5, 4, 0.9), pubs)
        self.assertEqual(g["A2"], "NOT GRADABLE")


if __name__ == "__main__":
    unittest.main()
