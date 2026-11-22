# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. Written by AI agents working for MV2 of Munim, Inc.
"""Unit tests for psp5ladder (synthetic data only, no network)."""
import math
import os
import sys
import unittest

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
for _sib in ("psp5min", "curv", "jack", "alfven"):
    _p = os.path.join(os.path.dirname(_HERE), _sib)
    if os.path.isdir(_p):
        sys.path.append(_p)

try:
    import numpy as np
    import p5min
    import psp5ladder as L
    HAVE = True
except Exception:  # noqa: BLE001
    HAVE = False


def _pink(n, seed):
    rng = np.random.default_rng(seed)
    f = np.fft.rfftfreq(n, d=0.5)
    amp = np.zeros_like(f)
    amp[1:] = f[1:] ** -0.5
    hi = f > 0.02
    amp[hi] = (0.02 ** -0.5) * (f[hi] / 0.02) ** (-5.0 / 6.0)
    out = np.zeros((n, 3))
    for c in range(3):
        a = amp * rng.normal(size=len(f))
        ph = rng.uniform(0, 2 * np.pi, size=len(f))
        out[:, c] = np.fft.irfft(a * np.exp(1j * ph), n=n)
    return out


def _with_peak(n=16200, seed=7, power=40.0):
    x = _pink(n, seed)
    t = np.arange(n) * 0.5
    dur = n * 0.5
    f0 = 26 / dur  # exactly on bin 26 (3.21 mHz at 8100 s)
    _, P = p5min.spectrum(x)
    pn = float(np.median(P[19:32]))
    a = math.sqrt(2.0 * power * pn / dur)
    y = x.copy()
    y[:, 2] += a * np.sin(2 * np.pi * f0 * t)
    return y


@unittest.skipUnless(HAVE and os.name == "posix", "needs numpy on POSIX")
class TestLadder(unittest.TestCase):

    def test_spec_matches_p5min(self):
        x = _pink(4000, 1)
        f1, P1 = p5min.spectrum(x)
        f2, P2 = L.spec(x, "none")
        self.assertTrue(np.allclose(f1, f2))
        self.assertTrue(np.allclose(P1, P2, rtol=1e-12, atol=0))

    def test_v0_equals_p5min_statistic(self):
        y = _with_peak()
        rec0 = p5min.analyse_window(y.astype(np.float32), p5min.layout(8100))
        rec = L.score(y.astype(np.float32), "V0")
        self.assertEqual(rec["valid"], 1)
        self.assertAlmostEqual(rec["fpk_mhz"], rec0["fpk_mhz"], places=3)
        self.assertAlmostEqual(rec["S_fpk"], float(rec0["S_fpk"]), places=9)
        self.assertGreater(rec["S_fpk"], 5.8)

    def test_baseline_hand_values(self):
        P = np.ones(100)
        f = np.arange(1, 101) * 1e-3
        bl = L.baseline(P, f, P * f)
        # cumulative = (i+1)/100: first >= 0.25 at i = 24, first >= 0.75 at i = 74
        self.assertEqual((bl["lo"], bl["hi"], bl["nb"]), (24, 74, 51))
        self.assertEqual(bl["med"], 49)
        blx = L.baseline(P, f, P * f, excl=(49, 51))
        self.assertEqual(blx["nb"], 48)
        self.assertEqual(blx["nb_full"], 51)
        # total limited to f <= 0.05 Hz: 50 bins, (i+1)/50 >= 0.25 at i = 12, >= 0.75 at i = 37
        blf = L.baseline(P, f, P * f, fmax=0.0505)
        self.assertEqual((blf["lo"], blf["hi"]), (12, 37))
        # cumulative of P*f = (i+1)(i+2)/(100*101): >= 0.25 at i = 49, >= 0.75 at i = 86
        blp = L.baseline(P, f, P * f, measure="Pf")
        self.assertEqual((blp["lo"], blp["hi"]), (49, 86))

    def test_baseline_mean_std(self):
        P = np.ones(100)
        f = np.arange(1, 101) * 1e-3
        y = np.zeros(100)
        y[24:75] = np.tile([1.0, 3.0], 26)[:51]  # 26 ones, 25 threes
        bl = L.baseline(P, f, y)
        mu = (26 * 1.0 + 25 * 3.0) / 51
        self.assertAlmostEqual(bl["mu"], mu, places=12)
        var0 = (26 * (1 - mu) ** 2 + 25 * (3 - mu) ** 2) / 51
        self.assertAlmostEqual(bl["sd"], math.sqrt(var0), places=12)
        bl1 = L.baseline(P, f, y, ddof=1)
        self.assertAlmostEqual(bl1["sd"], math.sqrt(var0 * 51 / 50), places=12)

    def test_too_few_bins(self):
        P = np.zeros(100)
        P[0] = 1.0
        f = np.arange(1, 101) * 1e-3
        bl = L.baseline(P, f, P * f)
        self.assertFalse(bl["ok"])

    def test_detrend(self):
        n = 1000
        t = np.arange(n, dtype=float)
        x = np.stack([3.0 + 0.2 * t, -1.0 * t, np.full(n, 5.0)], axis=1)
        d = L.detrend_linear(x)
        self.assertLess(float(np.abs(d).max()), 1e-9)

    def test_scale_invariance_all_variants(self):
        y = _with_peak()
        for v in L.VNAMES:
            r1 = L.score(y, v)
            r2 = L.score(7.3 * y, v)
            self.assertAlmostEqual(r1["S_fpk"], r2["S_fpk"], places=6, msg=v)

    def test_all_variants_find_injected_peak(self):
        y = _with_peak(power=80.0)
        for v in L.VNAMES:
            r = L.score(y, v)
            self.assertEqual(r["valid"], 1, msg=v)
            self.assertAlmostEqual(r["fpk_mhz"], 26 / 8.1, places=2, msg=v)

    def test_excluding_peak_raises_s(self):
        y = _with_peak(power=80.0)
        r0 = L.score(y, "V0")
        r2 = L.score(y, "V2")
        if not r0.get("peak_inside_base"):
            self.skipTest("peak outside baseline in this draw")
        self.assertGreater(r2["S_fpk"], r0["S_fpk"])
        self.assertEqual(r2["nb"], r0["nb"] - 3)

    def test_samuelson_ceiling(self):
        y = _with_peak(power=1000.0)
        r = L.score(y, "V0")
        if not r.get("peak_inside_base"):
            self.skipTest("peak outside baseline in this draw")
        self.assertLessEqual(r["S_fpk"], r["samuelson_ceiling"] + 1e-9)
        self.assertAlmostEqual(r["samuelson_ceiling"], math.sqrt(r["nb"] - 1))

    def test_gappy_invalid(self):
        x = _pink(8400, 3).astype(np.float32)
        x[1000:1200, :] = np.nan
        self.assertEqual(L.score(x, "V0")["valid"], 0)
        self.assertEqual(L.score(None, "V0")["valid"], 0)

    def test_edge_windows(self):
        w = L.edge_windows(1000.0, 9100.0)
        self.assertEqual(len(w), 8)
        self.assertEqual(w[0], ("start", -10, 400.0, 9100.0))
        self.assertEqual(w[-1], ("stop", 10, 1000.0, 9700.0))

    def test_segment(self):
        B = np.arange(30.0).reshape(10, 3)
        seg = L.segment(B, 100.0, 101.0, 102.0)  # bins 2..3
        self.assertEqual(seg.shape, (2, 3))
        self.assertEqual(float(seg[0, 0]), 6.0)
        self.assertIsNone(L.segment(B, 100.0, 99.0, 101.0))

    def _res(self, s22, s24, f22=3.2, f24=3.1, edges22=None, edges24=None):
        ev = {}
        for name, svals, fp, pub, edges in (("E22", s22, f22, 5.8, edges22),
                                            ("E24", s24, f24, 6.2, edges24)):
            var = {}
            for v, s in zip(L.VNAMES, svals):
                var[v] = {"valid": 1, "base_ok": 1 if math.isfinite(s) else 0,
                          "S_fpk": s, "fpk_mhz": fp}
            ev[name] = {"published_sigma": pub, "variants": var,
                        "edges": [{"S_fpk": s} for s in (edges or [])]}
        return {"events": ev}

    def test_grade_stable_reproduced_robust(self):
        r = self._res([5.9, 5.8, 6.1, 5.7, 5.5, 5.9, 3.0],
                      [6.0, 5.9, 6.3, 6.1, 6.2, 6.0, 9.0],
                      edges22=[5.5, 5.6, 5.7, 5.8, 5.4, 5.6, 5.9, 5.3],
                      edges24=[5.8, 5.9, 6.1, 6.2, 5.6, 5.9, 6.0, 6.3])
        g = L.grade(r)
        self.assertEqual(g["RL0"], "PASS")
        self.assertEqual(g["RL1"], "STABLE")          # V6 ignored in RL1
        self.assertEqual(g["RL2"], "REPRODUCED")
        self.assertEqual(g["RL2_selected"], "V0")
        self.assertEqual(g["RL3"], "ROBUST")          # drops 0.6 and 0.4

    def test_grade_fragile_other_reading_edge(self):
        r = self._res([3.0, 2.9, 3.5, 3.1, 5.2, 3.0, 4.0],
                      [4.0, 3.9, 4.6, 4.2, 6.0, 4.1, 4.0],
                      edges22=[1.0, 2.5, 2.8, 3.0, 2.9, 2.7, 3.1, 2.0],
                      edges24=[3.9, 3.8, 4.0, 4.1, 4.2, 3.7, 3.9, 4.0])
        g = L.grade(r)
        self.assertEqual(g["RL1"], "FRAGILE")         # E22 spread 2.3
        self.assertEqual(g["RL2"], "REPRODUCED-OTHER-READING")
        self.assertEqual(g["RL2_passing"], ["V4"])
        self.assertEqual(g["RL3"], "EDGE-SENSITIVE")  # E22 drop 2.0

    def test_grade_not_reproduced_selection(self):
        # nobody within 1.5: worst |dS| V0 2.8, V1 2.9, V2 2.2, V3 2.6, V4 2.5, V5 2.9, V6 2.2
        r = self._res([3.0, 2.9, 3.6, 3.2, 3.3, 2.9, 3.6],
                      [3.4, 3.3, 4.0, 3.6, 3.7, 3.3, 4.0])
        g = L.grade(r)
        self.assertEqual(g["RL2"], "NOT REPRODUCED")
        self.assertEqual(g["RL2_passing"], [])
        self.assertEqual(g["RL2_selected"], "V2")     # tie with V6 -> list order
        self.assertEqual(g["RL1"], "STABLE")          # spreads 0.7 and 0.7
        self.assertEqual(g["RL3"], "NOT GRADABLE")    # no edge windows given

    def test_grade_missing_value_not_gradable(self):
        nan = float("nan")
        r = self._res([5.9, nan, nan, 5.7, 5.5, 5.9, 3.0],
                      [6.0, 5.9, 6.3, 6.1, 6.2, 6.0, 9.0])
        g = L.grade(r)
        self.assertEqual(g["RL1"], "NOT GRADABLE")    # only 4 finite for E22
