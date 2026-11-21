# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. Written by AI agents working for MV2 of Munim, Inc.
"""Unit tests for psp5curv (synthetic data only, no network).
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
    import p5min as P
    import psp5curv as C
    HAVE = True
except Exception:  # noqa: BLE001
    HAVE = False


def _pink(n, seed):
    """Same synthetic spectrum as psp5min's tests: P ~ 1/f to 20 mHz, then f^-5/3."""
    rng = np.random.default_rng(seed)
    f = np.fft.rfftfreq(n, d=0.5)
    amp = np.zeros_like(f)
    amp[1:] = f[1:] ** -0.5
    hi = f > 0.02
    amp[hi] = (0.02 ** -0.5) * (f[hi] / 0.02) ** (-5.0 / 6.0)
    out = np.zeros((n, 3))
    for c in range(3):
        ph = rng.uniform(0, 2 * np.pi, size=len(f))
        out[:, c] = np.fft.irfft(amp * np.exp(1j * ph), n=n)
    return out


def _pink_noisy(n, seed):
    """Gaussian noise with the _pink spectrum (periodogram scatters like chi-square)."""
    rng = np.random.default_rng(seed)
    f = np.fft.rfftfreq(n, d=0.5)
    amp = np.zeros_like(f)
    amp[1:] = f[1:] ** -0.5
    hi = f > 0.02
    amp[hi] = (0.02 ** -0.5) * (f[hi] / 0.02) ** (-5.0 / 6.0)
    out = np.zeros((n, 3))
    for c in range(3):
        z = (rng.normal(size=len(f)) + 1j * rng.normal(size=len(f))) / math.sqrt(2.0)
        out[:, c] = np.fft.irfft(amp * z, n=n)
    return out


def _inject(x, mult=40.0):
    n = x.shape[0]
    t = np.arange(n) * 0.5
    _, Pt = P.spectrum(x)
    pn = float(np.median(Pt[19:32]))
    a = math.sqrt(2.0 * mult * pn / 8100.0)
    y = x.copy()
    y[:, 2] += a * np.sin(2 * np.pi * (26 / 8100.0) * t)
    return y


@unittest.skipUnless(HAVE and os.name == "posix", "needs numpy on POSIX")
class TestPsp5curv(unittest.TestCase):

    def test_w_int(self):
        # intercept weights of a + b d + c d^2 on d = +-2..+-6:
        # w_d = 1/10 - 18 (d^2 - 18) / 1308   (mean d^2 = 18, sum (d^2-18)^2 = 1308)
        w = C.W_INT
        self.assertEqual(len(w), 10)
        self.assertAlmostEqual(float(w[0]), -0.147706, places=5)   # d = -6
        self.assertAlmostEqual(float(w[1]), 0.003670, places=5)    # d = -5
        self.assertAlmostEqual(float(w[2]), 0.127523, places=5)    # d = -4
        self.assertAlmostEqual(float(w[3]), 0.223853, places=5)    # d = -3
        self.assertAlmostEqual(float(w[4]), 0.292661, places=5)    # d = -2
        self.assertAlmostEqual(float(w.sum()), 1.0, places=9)
        for i in range(5):
            self.assertAlmostEqual(float(w[i]), float(w[9 - i]), places=9)

    def test_l2_exact_on_quadratic_log_spectrum(self):
        k = np.arange(1, C.KMAX + 1, dtype=float)
        Y = np.exp(0.3 + 0.01 * k - 0.0005 * k * k)[None, :]
        for ks in ([24, 25, 26, 27, 28], [9, 10, 11, 12, 13], [77, 78, 79, 80, 81]):
            self.assertAlmostEqual(float(C.band_L2(Y, ks)[0]), 1.0, places=9)

    def test_band_L_matches_p5min(self):
        rng = np.random.default_rng(3)
        y = rng.exponential(size=8099) + rng.exponential(size=8099)
        Y = y[None, :C.KMAX]
        lay = P.layout(8100)["local"]
        for name, ks in lay.items():
            ref = max(P.local_ratio(y, k) for k in ks)
            self.assertAlmostEqual(float(C.band_L(Y, ks, len(y))[0]), ref, places=12, msg=name)

    def test_spectrum_matches_p5min(self):
        rng = np.random.default_rng(5)
        x = rng.normal(size=(16200, 3))
        f, Pc = C.comp_spectrum(x)
        f2, P2 = P.spectrum(x)
        self.assertTrue(np.allclose(f, f2))
        self.assertTrue(np.allclose(Pc.sum(axis=1), P2, rtol=1e-12, atol=0))

    def test_s_batch_matches_analyse_window(self):
        lay = P.layout(8100)
        for x in (_pink_noisy(16200, 7), _inject(_pink_noisy(16200, 7))):
            x32 = x.astype(np.float32)
            rec = P.analyse_window(x32, lay)
            self.assertEqual(rec["base_ok"], 1)
            f, Pc = C.comp_spectrum(C.prep(x32))
            sel = np.flatnonzero((f >= 2.5e-3) & (f <= 4.0e-3))
            r = C.s_batch(Pc.sum(axis=1)[None, :], f, lay["paper"]["P3.15"], sel)
            self.assertAlmostEqual(float(r["S_band"][0]), float(rec["S_P3.15"]), delta=1e-6)
            self.assertAlmostEqual(float(r["S_fpk"][0]), float(rec["S_fpk"]), delta=1e-6)
            self.assertAlmostEqual(float(f[r["ipk"][0]]) * 1e3, rec["fpk_mhz"], places=3)

    def test_samuelson_ceiling_holds(self):
        # (max - mean) / population sd over n values never exceeds sqrt(n - 1)
        rng = np.random.default_rng(11)
        f = np.arange(1, 4200) / 4200.0
        Pm = (f ** -1.0)[None, :] * rng.exponential(size=(400, 4199)) ** 3
        sel = np.flatnonzero((f >= 2.5e-3) & (f <= 4.0e-3))
        r = C.s_batch(Pm, f, [13], sel, min_bins=1)
        inb = (r["lo"] <= r["ipk"]) & (r["ipk"] <= r["hi"]) & (r["nbase"] >= 2)
        self.assertGreater(int(inb.sum()), 10)
        lim = np.sqrt(r["nbase"][inb] - 1.0)
        self.assertTrue(np.all(r["S_fpk"][inb] <= lim + 1e-9))

    def test_smooth_cubic_power_law(self):
        f = np.arange(1, 8100) / 8100.0
        Pc = np.stack([f ** -1.5, 2 * f ** -1.5, 3 * f ** -1.5], axis=1)
        S, rand, rms = C.smooth_cubic(f, Pc)
        # f <= 30 mHz at 135 min: k <= 243 (k = 243 sits on the edge in floating point)
        self.assertIn(int(rand.sum()), (242, 243))
        ratio = S[rand] / Pc[rand]
        self.assertTrue(np.allclose(ratio, math.exp(C.EULER_GAMMA), rtol=1e-8))
        self.assertTrue(np.all(S[~rand] == 0))
        self.assertLess(max(rms), 1e-8)

    def test_surrogate_test_runs(self):
        rng = np.random.default_rng(13)
        lay = P.layout(8100)
        out = {}
        for tag, x in (("plain", _pink_noisy(16200, 7)),
                       ("line", _inject(_pink_noisy(16200, 7)))):
            f, Pc = C.comp_spectrum(x)
            o = C.surrogate_test(rng, f, Pc, lay["paper"]["P3.15"], 2000, chunk=700)
            self.assertEqual(o["n_ok_fpk"], 2000)
            self.assertTrue(0.0 < o["p_fpk"] <= 1.0)
            self.assertAlmostEqual(o["p_fpk"], (o["hits_fpk"] + 1) / 2001.0)
            self.assertAlmostEqual(o["samuelson_ceiling"], math.sqrt(o["n_base"] - 1))
            out[tag] = o
        self.assertTrue(out["line"]["peak_in_base"])
        self.assertLessEqual(out["line"]["S_fpk_obs"], out["line"]["samuelson_ceiling"] + 1e-9)
        self.assertLess(out["line"]["p_fpk"], 0.05)

    def test_gauss_isf(self):
        self.assertAlmostEqual(C.gauss_isf(0.5), 0.0, places=6)
        self.assertAlmostEqual(C.gauss_isf(0.0013498980), 3.0, places=4)
        self.assertAlmostEqual(C.gauss_isf(2.8665157e-7), 5.0, places=4)
        self.assertTrue(math.isnan(C.gauss_isf(float("nan"))))

    def test_ranks_from(self):
        v = np.array([[3.0, 1.0, 2.0, 5.0], [0.5, 1.0, 2.0, 3.0], [9.0, 1.0, 2.0, 3.0]])
        self.assertEqual(C.ranks_from(v).tolist(), [3, 1, 4])

    def test_l2_usable_and_band_table_70(self):
        lay70 = P.layout(4200)["local"]
        self.assertFalse(C.l2_usable(lay70["C1.3"], C.KMAX))   # bins 4-6: 4 - 6 < 2
        self.assertTrue(C.l2_usable(lay70["C2.1"], C.KMAX))    # bins 8-10: 8 - 6 = 2
        rng = np.random.default_rng(2)
        Y = rng.exponential(size=(5, C.KMAX))
        names, L, use2, L2 = C.band_table(Y, lay70, 4199)
        self.assertEqual(L.shape, (5, 9))
        self.assertEqual(len(use2), 8)
        self.assertNotIn("C1.3", use2)
        self.assertEqual(L2.shape, (5, 8))
        lay135 = P.layout(8100)["local"]
        names, L, use2, L2 = C.band_table(Y, lay135, 8099)
        self.assertEqual(len(use2), 9)

    def test_sim_flat_l2_unbiased(self):
        # flat P*f: every band has the same fit design and iid ordinates,
        # so L2 ranks are exactly uniform and Z has mean 0, sd 1.
        rng = np.random.default_rng(17)
        f = np.arange(1, 8100) / 8100.0
        s = C.run_sims(rng, 400, f, P.layout(8100)["local"], 8099, kinds=("F",), nrep=30)
        self.assertLess(abs(s["F"]["Z_L2_mean"]), 0.8)
        self.assertTrue(0.5 < s["F"]["Z_L2_sd"] < 1.6)

    def test_knee_fit(self):
        f = np.arange(1, 8100) / 8100.0
        fk = f[:C.KMAX]
        rows = []
        for ub in (math.log10(3e-3), math.log10(8e-3)):
            up = np.log10(fk) - ub
            lp = -1.0 * up + (-0.6) * C.KNEE_DELTA * np.log10(1 + 10 ** (up / C.KNEE_DELTA))
            rows.append(10 ** lp)
        kn, b1, b2 = C.knee_fit(f, np.stack(rows))
        # grid step 0.04 dex: worst case 0.023 dex from 3 mHz, 0.003 dex from 8 mHz
        self.assertLess(abs(math.log10(kn[0] / 3.0)), 0.045)
        self.assertLess(abs(math.log10(kn[1] / 8.0)), 0.045)
        for i in range(2):
            self.assertAlmostEqual(float(b1[i]), -1.0, delta=0.15)
            self.assertAlmostEqual(float(b2[i]), -1.6, delta=0.15)

    def test_grade(self):
        base = {"sims": {"F": {"Z_L2_mean": 0.1}, "K": {"Z_L2_mean": -0.2}},
                "D135": {"rank": {"L2": {"all": {"z": 3.5, "n": 3000}}}},
                "published": [{"p_fpk": 5e-4, "n_ok_fpk": 20000},
                              {"p_fpk": 2e-3, "n_ok_fpk": 40000}]}
        g = C.grade(base)
        self.assertEqual((g["SIMGATE"], g["H2c"], g["H4"]), ("PASS", "SUPPORT", "MODERATE"))
        b2 = dict(base, sims={"F": {"Z_L2_mean": 0.1}, "K": {"Z_L2_mean": 0.7}})
        g2 = C.grade(b2)
        self.assertEqual((g2["SIMGATE"], g2["H2c"]), ("FAIL", "NOT GRADABLE"))
        b3 = dict(base, D135={"rank": {"L2": {"all": {"z": 1.2, "n": 3000}}}},
                  published=[{"p_fpk": 4e-5, "n_ok_fpk": 20000},
                             {"p_fpk": 5e-5, "n_ok_fpk": 40000}])
        g3 = C.grade(b3)
        self.assertEqual((g3["H2c"], g3["H4"]), ("NO-EXCESS", "EXTREME"))
        b4 = dict(base, D135={"rank": {"L2": {"all": {"z": 2.5, "n": 3000}}}},
                  published=[{"p_fpk": 4e-5, "n_ok_fpk": 20000},
                             {"p_fpk": 0.02, "n_ok_fpk": 40000}])
        g4 = C.grade(b4)
        self.assertEqual((g4["H2c"], g4["H4"]), ("INDETERMINATE", "ORDINARY"))
        b5 = dict(base, D135={"rank": {"L2": {"all": {"z": 3.5, "n": 1500}}}},
                  published=[{"p_fpk": 4e-5, "n_ok_fpk": 20000}])
        g5 = C.grade(b5)
        self.assertEqual((g5["H2c"], g5["H4"]), ("NOT GRADABLE", "NOT GRADABLE"))

    def test_prep_gap_rule(self):
        x = _pink(8400, 3).astype(np.float32)
        x[1000:1200, :] = np.nan        # 100 s gap > 60 s
        self.assertIsNone(C.prep(x))
        x2 = _pink(8400, 3).astype(np.float32)
        x2[1000:1050, :] = np.nan
        y = C.prep(x2)
        self.assertIsNotNone(y)
        self.assertFalse(np.isnan(y).any())


if __name__ == "__main__":
    unittest.main()
