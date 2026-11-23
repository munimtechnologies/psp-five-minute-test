# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. Written by AI agents working for MV2 of Munim, Inc.
"""Unit tests for psp2min (synthetic data only, no network)."""
import math
import os
import sys
import unittest

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in (os.path.join(os.path.dirname(_HERE), "psp5min"), _HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)

try:
    import numpy as np
    import psp2min as M
    HAVE = True
except Exception:  # noqa: BLE001
    HAVE = False


def synth(spec, n, seed):
    """(n,3) random-phase series whose PSD follows spec(f)."""
    rng = np.random.default_rng(seed)
    f = np.fft.rfftfreq(n, d=M.DT)
    amp = np.zeros_like(f)
    amp[1:] = np.sqrt(spec(f[1:]))
    x = np.zeros((n, 3))
    for c in range(3):
        x[:, c] = np.fft.irfft(amp * np.exp(1j * rng.uniform(0, 2 * np.pi, len(f))), n=n)
    return x


@unittest.skipUnless(HAVE and os.name == "posix", "needs numpy on POSIX")
class TestPsp2min(unittest.TestCase):
    N = 28800

    def test_f_lo(self):
        # 4-hour window: 50 % CoI period = 1.033 * 14400 / (4 sqrt 2) = 2630 s
        self.assertAlmostEqual(1.0 / M.f_lo_for(14400), 2629.6, delta=2.0)

    def test_kolmogorov_saturation(self):
        # paper (Appendix A): ~13 min for f^-5/3 in a 4-hour window; continuum value 758 s
        k = M.kolmogorov_check()
        self.assertGreater(k, 600.0)
        self.assertLess(k, 900.0)

    def test_break_spectrum(self):
        fb = 8e-3

        def spec(f):
            return np.where(f < fb, f ** -0.5, fb ** -0.5 * (f / fb) ** (-5.0 / 3.0))
        x = synth(spec, self.N, 3)
        f, Pw = M.trace_psd(x)
        st = M.fmid_stats(f, Pw, M.f_lo_for(14400))
        # continuum calculation gives 1/f_mid = 133 s
        self.assertGreater(1.0 / st["fmid_geo"], 95.0)
        self.assertLess(1.0 / st["fmid_geo"], 190.0)
        # P*f peaks at the break (smoothed by a factor of 2)
        self.assertGreater(st["fpeak"], fb / 2.5)
        self.assertLess(st["fpeak"], fb * 2.5)
        self.assertEqual(st["fpeak_at_floor"], 0)

    def test_trend_detrend_insensitive(self):
        x = synth(lambda f: f ** (-5.0 / 3.0), self.N, 4)
        f, P1 = M.trace_psd(x)
        y = x + np.linspace(0, 50 * x.std(), self.N)[:, None]
        _, P2 = M.trace_psd(y)
        a = M.fmid_stats(f, P1, M.f_lo_for(14400))["fmid_geo"]
        b = M.fmid_stats(f, P2, M.f_lo_for(14400))["fmid_geo"]
        self.assertAlmostEqual(a, b, delta=1e-9 + 1e-6 * a)

    def test_parseval_hann(self):
        rng = np.random.default_rng(1)
        x = rng.normal(size=(self.N, 3))
        f, Pw = M.trace_psd(x)
        df = f[1] - f[0]
        self.assertAlmostEqual(float(Pw.sum() * df) / 3.0, 1.0, delta=0.03)

    def test_fill_gaps(self):
        x = synth(lambda f: f ** -1.5, self.N, 5).astype(np.float32)
        x[100:400] = np.nan          # 150 s > 60 s
        y, frac, run = M.fill_gaps(x)
        self.assertIsNone(y)
        x2 = synth(lambda f: f ** -1.5, self.N, 5).astype(np.float32)
        x2[100:150] = np.nan
        y2, frac2, run2 = M.fill_gaps(x2)
        self.assertIsNotNone(y2)
        self.assertEqual(run2, 50)
        self.assertFalse(np.isnan(y2).any())

    def test_window_plasma(self):
        pos_t = np.array([0.0, 1e6])
        pos_r = np.array([10.0, 10.0]) / M.P.RS_PER_AU   # 10 Rsun
        ts = np.arange(0.0, 14400.0, 2.0)
        vr = np.full(len(ts), 300.0)
        vr[::7] += 30.0
        out = M.window_plasma(ts, vr, 0.0, 14400.0, pos_t, pos_r)
        self.assertAlmostEqual(out["vr"], 300.0)
        self.assertAlmostEqual(out["tadv_h"], 9 * 695700.0 / 300.0 / 3600.0, places=6)
        self.assertAlmostEqual(out["v_cover"], 1.0)
        self.assertLess(out["dtadv_ratio"], 0.005)    # constant 10-min medians
        vr2 = vr.copy()
        vr2[ts >= 7200] = 600.0                        # abrupt speed change
        out2 = M.window_plasma(ts, vr2, 0.0, 14400.0, pos_t, pos_r)
        self.assertGreater(out2["dtadv_ratio"], 0.2)
        out3 = M.window_plasma(ts[:1000], vr[:1000], 0.0, 14400.0, pos_t, pos_r)
        self.assertLess(out3["v_cover"], M.V_COVER_MIN)
        bad = np.full(len(ts), np.nan)
        out4 = M.window_plasma(ts, bad, 0.0, 14400.0, pos_t, pos_r)
        self.assertTrue(math.isnan(out4["tadv_h"]))

    def test_spearman_and_rank(self):
        self.assertAlmostEqual(M.spearman([1, 2, 3, 4, 5], [2, 4, 6, 8, 10]), 1.0)
        self.assertAlmostEqual(M.spearman([1, 2, 3, 4, 5], [5, 4, 3, 2, 1]), -1.0)
        r = M.rankdata([10, 20, 20, 30])
        self.assertEqual(list(r), [1.0, 2.5, 2.5, 4.0])

    def test_gmean_boot(self):
        self.assertAlmostEqual(M.gmean([1, 100]), 10.0)
        g = {1: [100.0] * 10, 2: [100.0] * 5, 3: [100.0] * 7}
        lo, hi = M.block_boot_gmean(g, n_boot=500)
        self.assertAlmostEqual(lo, 100.0)
        self.assertAlmostEqual(hi, 100.0)
        g2 = {1: [60.0] * 10, 2: [240.0] * 10}
        lo2, hi2 = M.block_boot_gmean(g2, n_boot=2000)
        self.assertLess(lo2, 70.0)
        self.assertGreater(hi2, 200.0)
        a = M.block_boot_gmean(g2, n_boot=300)
        b = M.block_boot_gmean(g2, n_boot=300)
        self.assertEqual(a, b)                          # seeded

    def test_grades(self):
        self.assertEqual(M.grade_h1(120.0, 50, 5), "CONSISTENT")
        self.assertEqual(M.grade_h1(89.9, 50, 5), "SHORTER")
        self.assertEqual(M.grade_h1(180.1, 50, 5), "LONGER")
        self.assertEqual(M.grade_h1(120.0, 29, 5), "NOT GRADABLE")
        self.assertEqual(M.grade_h1(120.0, 50, 2), "NOT GRADABLE")
        self.assertEqual(M.grade_h2(100.0, 30, 110.0, 30)[0], "FLOOR")
        self.assertEqual(M.grade_h2(70.0, 30, 110.0, 30)[0], "DECLINING")
        self.assertEqual(M.grade_h2(85.0, 30, 110.0, 30)[0], "INTERMEDIATE")
        self.assertEqual(M.grade_h2(85.0, 19, 110.0, 30)[0], "NOT GRADABLE")
        self.assertEqual(M.grade_r1(0.5, 150.0, 40), "PASS")
        self.assertEqual(M.grade_r1(0.2, 150.0, 40), "FAIL")
        self.assertEqual(M.grade_r1(0.5, 300.0, 40), "FAIL")
        self.assertEqual(M.grade_r1(0.5, 150.0, 10), "NOT GRADABLE")

    def test_analyse_encounter_and_summary(self):
        # two synthetic "encounters" of 3 windows each; break at 8 mHz -> ~2 min
        fb = 8e-3

        def spec(f):
            return np.where(f < fb, f ** -0.5, fb ** -0.5 * (f / fb) ** (-5.0 / 3.0))
        n = 3 * self.N
        Bv = synth(spec, n, 9).astype(np.float32)
        ts = np.arange(0.0, n * M.DT, 3.0)
        pos_t = np.array([-1e7, 1e7])
        pos_r = np.array([10.0, 10.0]) / M.P.RS_PER_AU
        vr = np.full(len(ts), 300.0) + 5.0 * np.sin(ts / 900.0)
        rows = []
        for enc in (5, 22):
            rows += M.analyse_encounter(enc, 0.0, Bv, 0.0, ts, vr, pos_t, pos_r)
        self.assertEqual(len(rows), 6)
        self.assertTrue(all(r["usable"] == 1 for r in rows))
        self.assertTrue(all(r["clean_ok"] == 1 for r in rows))
        res = M.summarise(rows, [5, 22])
        self.assertEqual(res["n_usable"], {"A_E1_E19": 3, "B_E20_E27": 3})
        self.assertGreater(res["H1"]["gmean_s"], 95.0)
        self.assertLess(res["H1"]["gmean_s"], 190.0)
        self.assertEqual(res["grades"]["H1"], "NOT GRADABLE")   # n < 30
        self.assertTrue(res["clean_cut"]["cut_used_for_grades"])


if __name__ == "__main__":
    unittest.main()
