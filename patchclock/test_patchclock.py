# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. Written by AI agents working for MV2 of Munim, Inc.
"""Unit tests for patchclock (synthetic data only, no network)."""
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
    import patchclock as M
    HAVE = True
except Exception:  # noqa: BLE001
    HAVE = False


def e_out(r_rs):
    tp = M.OMEGA_RAD_S * r_rs * M.RS_KM / M.V_REF
    v = np.array([1.0, -tp, 0.0])
    return v / np.linalg.norm(v)


def fake_track(n, seed, rate_max=4.0, period_h=48.0):
    """Block times (s), a footpoint longitude whose rate swings between 0 and
    rate_max deg/h (corotation turn-arounds included), one position segment."""
    t = (np.arange(n) + 0.5) * M.BLOCK_S
    th = t / 3600.0
    rate = 0.5 * rate_max * (1.0 - np.cos(2 * np.pi * (th + 7.0 * seed) / period_h))
    phi = np.cumsum(rate) * (M.BLOCK_S / 3600.0)
    return t, phi


@unittest.skipUnless(HAVE and os.name == "posix", "needs numpy on POSIX")
class TestPatchclock(unittest.TestCase):

    def test_interp_gapped(self):
        t = np.array([0.0, 3600.0, 7200.0, 7200.0 + 5 * 3600.0])
        x = np.array([0.0, 1.0, 2.0, 7.0])
        y = M.interp_gapped(t, x, np.array([1800.0, 5400.0, 9000.0, -10.0, 1e9]))
        self.assertAlmostEqual(y[0], 0.5)
        self.assertAlmostEqual(y[1], 1.5)
        self.assertTrue(math.isnan(y[2]))      # inside a 5-hour gap
        self.assertTrue(math.isnan(y[3]))
        self.assertTrue(math.isnan(y[4]))

    def test_prep_positions_rotating_frame(self):
        t = np.arange(0, 49) * 3600.0
        hgi = (100.0 + 20.0 * t / 3600.0) % 360.0
        hg = (hgi - M.OMEGA_DEG_DAY * t / 86400.0) % 360.0
        pos = M.prep_positions(t, np.full(len(t), 0.1), hgi, hg)
        self.assertTrue(np.allclose(pos["omega"], 20.0 * 24 - M.OMEGA_DEG_DAY, atol=1e-6))
        self.assertAlmostEqual(pos["frame_check_deg_day"], -M.OMEGA_DEG_DAY, places=5)
        self.assertAlmostEqual(pos["lon_u"][-1] - pos["lon_u"][0], 48 * 20.0 - 2 * M.OMEGA_DEG_DAY, places=5)
        self.assertAlmostEqual(pos["r_rs"][0], 0.1 * 215.032, places=6)
        self.assertEqual(int(pos["seg"].max()), 0)

    def test_block_activity(self):
        r = 20.0
        e = e_out(r)
        nblk = 18
        B = np.tile(100.0 * e, (nblk * M.NB, 1))
        # hour 1 (blocks 0-5) outward; block 3: 30 % fully reversed samples
        s = 3 * M.NB
        B[s:s + int(0.3 * M.NB)] = -100.0 * e
        # hour 2 (blocks 6-11) inward; block 8: half the samples deflected 70 deg
        B[6 * M.NB:12 * M.NB] = -100.0 * e
        d70 = math.cos(math.radians(70)) * (-e) + math.sin(math.radians(70)) * np.array([0, 0, 1.0])
        s = 8 * M.NB
        B[s:s + M.NB // 2] = 100.0 * d70
        # hour 3 (blocks 12-17): half outward, half inward -> ambiguous polarity
        B[12 * M.NB:15 * M.NB] = 100.0 * e
        B[15 * M.NB:18 * M.NB] = -100.0 * e
        # block 1: too many missing samples
        B[1 * M.NB:1 * M.NB + M.NB // 2] = np.nan
        act = M.block_activity(B.astype(np.float32), np.full(nblk, r))
        self.assertEqual(act["nblk"], nblk)
        self.assertAlmostEqual(act["f60"][3], 0.3, places=3)
        self.assertAlmostEqual(act["f90"][3], 0.3, places=3)
        self.assertAlmostEqual(act["f60"][0], 0.0)
        self.assertAlmostEqual(act["f60"][8], 0.5, places=3)
        self.assertAlmostEqual(act["f90"][8], 0.0, places=3)
        self.assertEqual(act["pol"][7], -1.0)
        self.assertFalse(act["valid"][1])
        self.assertTrue(math.isnan(act["f60"][1]))
        self.assertFalse(act["valid"][12:].any())
        far = M.block_activity(B[:6 * M.NB].astype(np.float32), np.full(6, 60.0))
        self.assertFalse(far["valid"].any())   # beyond R_MAX_RS

    def test_window_median_v_and_footpoint(self):
        t = np.arange(0, 86400, 60.0)
        v = np.where(t < 43200, 300.0, 500.0)
        v[5] = -1e31
        out = M.window_median_v(t, v, np.array([10800.0, 75600.0, 200000.0]))
        self.assertAlmostEqual(out[0], 300.0)
        self.assertAlmostEqual(out[1], 500.0)
        self.assertTrue(math.isnan(out[2]))
        fp = M.footpoint(np.array([10.0, 10.0]), np.array([1.0, 11.0]), np.array([400.0, 400.0]))
        self.assertAlmostEqual(fp[0], 10.0)
        self.assertAlmostEqual(fp[1], 10.0 + M.OMEGA_DEG_S * 10 * M.RS_KM / 400.0, places=9)

    def test_cell_sums(self):
        z = np.array([0.0, 1.0, 3.0, np.nan, 2.0])
        phi = np.array([0.0, 0.2, 0.6, 1.0, 5.0])
        seg = np.array([0, 0, 0, 0, 1])
        ok = np.isfinite(z)
        S = M.cell_sums(z, phi, seg, ok)
        # lag 1: (0,1) dphi 0.2 -> bin [0.1,0.25); (1,2) dphi 0.4 -> [0.25,0.5); (3,4) excluded
        self.assertEqual(S[0, 0].sum(), 2)
        self.assertEqual(S[0, 0, 1], 1)
        self.assertEqual(S[0, 0, 2], 1)
        self.assertAlmostEqual(S[1, 0, 1], 0.5)
        self.assertAlmostEqual(S[1, 0, 2], 2.0)
        self.assertAlmostEqual(S[2, 0].sum(), 2 * M.BLOCK_S / 3600.0)
        # lag 2: (0,2) dphi 0.6 -> [0.5,0.75)
        self.assertEqual(S[0, 1].sum(), 1)
        self.assertAlmostEqual(S[1, 1, 3], 4.5)
        self.assertAlmostEqual(S[3, 1, 3], 0.6)

    def test_fit_recovers_model(self):
        S = np.zeros((4, M.NLAG, M.NPHI))
        mids = 0.5 * (M.DPHI_EDGES[:-1] + np.minimum(M.DPHI_EDGES[1:], 200.0))
        for k in range(M.NLAG):
            dt = (k + 1) * M.BLOCK_S / 3600.0
            for j in range(M.NPHI):
                g = 0.2 + 0.8 * (1.0 - math.exp(-dt / 5.0 - mids[j] / 4.0))
                S[:, k, j] = (100.0, 100.0 * g, 100.0 * dt, 100.0 * mids[j])
        f, D = M.fit_cells(S)
        self.assertIsNotNone(f)
        self.assertTrue(5.0 / 1.2 <= f["tau_h"] <= 5.0 * 1.2, f)
        self.assertTrue(4.0 / 1.2 <= f["theta_deg"] <= 4.0 * 1.2, f)
        self.assertAlmostEqual(f["c0"], 0.2, delta=0.05)
        self.assertGreater(f["r2"], 0.99)

    def test_efold_time(self):
        lags = (np.arange(M.NLAG) + 1) * M.BLOCK_S / 3600.0
        g = 1.0 - np.exp(-lags / 2.0)
        n = np.full(M.NLAG, 100.0)
        te = M.efold_time(100.0 * g, n)
        self.assertTrue(1.9 <= te <= 2.4, te)
        self.assertTrue(math.isnan(M.efold_time(np.full(M.NLAG, 50.0), n)))

    def test_ou_noise(self):
        rng = np.random.default_rng(3)
        x = M.ou_noise(np.arange(20000) * 0.5, 2.0, 6.0, rng)
        self.assertAlmostEqual(float(x.std()), 2.0, delta=0.4)
        self.assertTrue(np.all(M.ou_noise(np.arange(5.0), 0.0, 6.0, rng) == 0))

    def _pipeline(self, tau, theta, seed):
        rng = np.random.default_rng(seed)
        tot = np.zeros((4, M.NLAG, M.NPHI))
        for i in range(3):
            t, phi = fake_track(700, i)
            a = M.sim_activity(t / 3600.0, phi, tau, theta, 0.5, rng)
            ok = np.ones(len(a), bool)
            z, okz = M.zscore(a, ok)
            tot += M.cell_sums(z, phi, np.zeros(len(a), np.int64), okz)
        f, _ = M.fit_cells(tot)
        return f

    def test_injection_spatial(self):
        f = self._pipeline(math.inf, 5.0, 101)
        self.assertTrue(2.5 <= f["theta_deg"] <= 10.0, f)
        self.assertGreaterEqual(f["tau_h"], 15.0, f)

    def test_injection_temporal(self):
        f = self._pipeline(5.0, math.inf, 102)
        self.assertTrue(2.5 <= f["tau_h"] <= 10.0, f)
        self.assertGreaterEqual(f["theta_deg"], 15.0, f)

    def test_grades(self):
        self.assertEqual(M.grade_h1(5.0, [3.0, 9.0]), "SUPERGRANULAR")
        self.assertEqual(M.grade_h1(1.0, [0.5, 3.0]), "SMALLER")
        self.assertEqual(M.grade_h1(500.0, [40.0, 2000.0]), "NONE")
        self.assertEqual(M.grade_h1(12.0, [4.0, 60.0]), "UNRESOLVED")
        self.assertEqual(M.grade_h1(float("nan"), [1, 2]), "NOT GRADABLE")
        self.assertEqual(M.grade_h2(5.0, [3.0, 9.0], True), "SHORT")
        self.assertEqual(M.grade_h2(100.0, [30.0, 2000.0], True), "LONG")
        self.assertEqual(M.grade_h2(100.0, [30.0, 2000.0], False), "LONG (UNCALIBRATED)")
        self.assertEqual(M.grade_h2(12.0, [6.0, 30.0], True), "MIDDLE")
        self.assertEqual(M.grade_h2(12.0, [2.0, 300.0], True), "UNRESOLVED")
        self.assertEqual(M.verdict("PASS", "PASS", "SUPERGRANULAR", "LONG"), "FIXED SUPERGRANULAR IMPRINT")
        self.assertEqual(M.verdict("PASS", "PASS", "NONE", "SHORT"), "TRANSIENT BURSTS")
        self.assertEqual(M.verdict("PASS", "PASS", "SUPERGRANULAR", "SHORT"), "BOTH: SIZED AND SHORT-LIVED")
        self.assertEqual(M.verdict("PASS", "PASS", "UNRESOLVED", "SHORT"), "UNDETERMINED")
        self.assertEqual(M.verdict("FAIL", "PASS", "SUPERGRANULAR", "LONG"), "NOT GRADABLE")
        info = [{"enc": e, "used_fp": 300} for e in range(1, 28)]
        self.assertEqual(M.grade_g0(info)[0], "PASS")
        info2 = [{"enc": e, "used_fp": 300 if e < 22 else 10} for e in range(1, 28)]
        self.assertEqual(M.grade_g0(info2)[0], "FAIL")
        inj = {k: [{"pass": True}, {"pass": True}] for k in ("S1_spatial", "S2_temporal", "S3_mixed")}
        inj["S4_long"] = [{"pass": False}, {"pass": True}]
        self.assertEqual(M.grade_g1(inj), ("PASS", False))
        inj["S2_temporal"][1]["pass"] = False
        self.assertEqual(M.grade_g1(inj)[0], "FAIL")
        self.assertTrue(M.scenario_pass("S1_spatial", {"tau_h": 2000.0, "theta_deg": 5.3}))
        self.assertFalse(M.scenario_pass("S1_spatial", {"tau_h": 20.0, "theta_deg": 5.3}))
        self.assertFalse(M.scenario_pass("S3_mixed", None))


if __name__ == "__main__":
    unittest.main()
