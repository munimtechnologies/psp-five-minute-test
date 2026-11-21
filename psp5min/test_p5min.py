# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. Written by AI agents working for MV2 of Munim, Inc.
"""Unit tests for psp5min (synthetic data only, no network)."""
import math
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    import numpy as np
    import p5min
    HAVE = True
except Exception:  # noqa: BLE001
    HAVE = False


@unittest.skipUnless(HAVE and os.name == "posix", "needs numpy on POSIX")
class TestP5min(unittest.TestCase):

    def test_time_offset(self):
        # 2017-01-01T00:00:00Z <-> TT2000 536500869184000000 ns
        u = p5min.tt2000_to_unix(536500869184000000)
        self.assertAlmostEqual(float(u), 1483228800.0, places=3)

    def test_time_offset_cdflib(self):
        try:
            import cdflib
            tt = cdflib.cdfepoch.compute_tt2000([2025, 1, 1, 0, 0, 0, 0, 0, 0])
            tt = int(np.asarray(tt).ravel()[0])
        except Exception as e:  # noqa: BLE001
            self.skipTest(f"cdflib unavailable or API differs: {e!r}")
        u = float(p5min.tt2000_to_unix(tt))
        self.assertAlmostEqual(u, p5min.to_unix("2025-01-01T00:00:00Z"), places=3)

    def test_nearest_bins(self):
        self.assertEqual(p5min.nearest_bins(3.15, 8100, 2), [25, 26])
        self.assertEqual(p5min.nearest_bins(3.15, 4200, 1), [13])
        self.assertEqual(p5min.nearest_bins(3.2, 8100, 5), [24, 25, 26, 27, 28])
        self.assertEqual(p5min.nearest_bins(3.2, 4200, 3), [12, 13, 14])

    def test_layout(self):
        for d in (135, 70):
            lay = p5min.layout(d * 60)
            self.assertEqual(len(lay["local"]), 9)
            self.assertEqual(len(lay["paper"]), 9)
            sizes = {len(v) for v in lay["local"].values()}
            self.assertEqual(len(sizes), 1)
        self.assertEqual(p5min.layout(8100)["n_local"], 5)
        self.assertEqual(p5min.layout(4200)["n_local"], 3)
        self.assertEqual(p5min.layout(8100)["n_paper"], 2)
        self.assertEqual(p5min.layout(4200)["n_paper"], 1)

    def test_longest_run(self):
        m = np.array([0, 1, 1, 0, 1, 1, 1, 0], bool)
        self.assertEqual(p5min.longest_run(m), 3)
        self.assertEqual(p5min.longest_run(np.zeros(5, bool)), 0)
        self.assertEqual(p5min.longest_run(np.ones(4, bool)), 4)

    def test_clusters(self):
        cl = p5min.clusters_from_starts([0, 1200, 2400, 6000, 9600, 10800], 1200)
        self.assertEqual([len(c) for c in cl], [3, 1, 2])

    def test_binner(self):
        b = p5min.Binner(100.0, 4)
        t = 100.0 + np.arange(20) * 0.1
        v = np.stack([np.arange(20.0)] * 3, axis=1)
        v[3, 0] = -1e31
        b.add(t, v)
        out = b.result()
        self.assertAlmostEqual(float(out[1, 0]), 7.0, places=4)
        # fill value drops the whole sample (index 3) from bin 0: (0+1+2+4)/4
        self.assertAlmostEqual(float(out[0, 0]), 1.75, places=4)
        self.assertAlmostEqual(float(out[0, 1]), 1.75, places=4)

    def test_spectrum_parseval(self):
        rng = np.random.default_rng(1)
        x = rng.normal(size=(4000, 3))
        f, P = p5min.spectrum(x)
        var = (x - x.mean(axis=0)).var(axis=0).sum()
        df = f[1] - f[0]
        self.assertAlmostEqual(float(P.sum() * df) / var, 1.0, delta=0.01)
        self.assertAlmostEqual(float(f[0]), 1.0 / (4000 * 0.5))

    def test_local_ratio(self):
        y = np.ones(200)
        y[49] = 10.0  # bin k = 50
        self.assertAlmostEqual(p5min.local_ratio(y, 50), 10.0)
        self.assertAlmostEqual(p5min.local_ratio(y, 60), 1.0)
        self.assertTrue(math.isfinite(p5min.local_ratio(y, 1)))

    def _pink(self, n, seed):
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

    def test_injected_peak_found(self):
        n = 16200
        x = self._pink(n, 7)
        lay = p5min.layout(8100)
        base = p5min.analyse_window(x.astype(np.float32), lay)
        self.assertEqual(base["valid"], 1)
        t = np.arange(n) * 0.5
        f0 = 26 / 8100.0  # exactly on bin 26 (3.21 mHz)
        _, P = p5min.spectrum(x)
        pn = float(np.median(P[19:32]))  # noise level near bin 26
        a = math.sqrt(2.0 * 40.0 * pn / 8100.0)  # sine power = 40 x noise
        y = x.copy()
        y[:, 2] += a * np.sin(2 * np.pi * f0 * t)
        rec = p5min.analyse_window(y.astype(np.float32), lay)
        self.assertEqual(rec["base_ok"], 1)
        self.assertAlmostEqual(rec["fpk_mhz"], 3.2099, places=3)
        self.assertGreater(rec["S_P3.15"], 5.8)
        self.assertEqual(rec["rank_T"], 9)
        self.assertGreater(rec["L_T3.2"], 20.0)

    def test_gappy_window_rejected(self):
        x = self._pink(8400, 3).astype(np.float32)
        x[1000:1200, :] = np.nan  # 200 bins = 100 s gap > 60 s
        rec = p5min.analyse_window(x, p5min.layout(4200))
        self.assertEqual(rec["valid"], 0)
        x2 = self._pink(8400, 3).astype(np.float32)
        x2[1000:1050, :] = np.nan
        rec2 = p5min.analyse_window(x2, p5min.layout(4200))
        self.assertEqual(rec2["valid"], 1)

    def test_perihelia(self):
        t = np.arange(0, 200 * 86400, 3600.0)
        r = 0.6 + 0.5 * np.cos(2 * np.pi * t / (88 * 86400))
        r = np.abs(r)
        p = p5min.find_perihelia(t, r)
        self.assertEqual(len(p), 2)
        self.assertAlmostEqual(p[0][0] / 86400, 44, delta=0.1)

    def test_hapi_parse(self):
        txt = ("2024-12-24T11:00:00.000Z,0.0461\n"
               "2024-12-24T12:00:00.000Z,0.0459\n"
               "bad,line\n"
               "2024-12-24T13:00:00.000Z,-1.0E31\n")
        t, v = p5min.parse_hapi_csv(txt)
        self.assertEqual(len(t), 2)
        self.assertAlmostEqual(v[1], 0.0459)

    def test_rank_and_binom(self):
        z = p5min.rank_z([5] * 100)
        self.assertAlmostEqual(z["z"], 0.0)
        z9 = p5min.rank_z([9] * 4)
        self.assertGreater(z9["z"], 3.0)
        self.assertAlmostEqual(p5min.binom_sf(0, 10, 0.1), 1.0)
        self.assertAlmostEqual(p5min.binom_sf(10, 10, 0.5), 0.5 ** 10)

    def test_overlap_pub(self):
        a = p5min.to_unix("2024-12-24T11:00:00Z")
        self.assertEqual(p5min.overlaps_pub(a, a + 600), 1)
        self.assertEqual(p5min.overlaps_pub(a - 86400, a - 86000), 0)


if __name__ == "__main__":
    unittest.main()
