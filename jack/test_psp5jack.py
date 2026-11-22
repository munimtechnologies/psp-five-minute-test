# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. Written by AI agents working for MV2 of Munim, Inc.
"""Unit tests for psp5jack 0.1.0. Expected values worked out by hand."""
import math
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import psp5jack as J  # noqa: E402

T0 = J.parse_time("2024-12-20T00:00:00Z")
STEP = 135 * 60.0


def rows_from(groups):
    """groups: {enc: [rank, ...]} -> rows with time-adjacent 135-min starts."""
    out = []
    for enc, ranks in groups.items():
        for i, r in enumerate(ranks):
            out.append((enc, T0 + enc * 86400 * 30 + i * STEP, float(r)))
    return out


class TestFormulas(unittest.TestCase):
    def test_zero_mean_but_overdispersed(self):
        # A [9,9], B [1,1], C [5,5]: rbar 5, Z 0; S = 8, -8, 0; V = 1.5*128/36 = 5.3333
        rows = rows_from({1: [9, 9], 2: [1, 1], 3: [5, 5]})
        cr = J.cluster_robust(rows)
        self.assertEqual(cr["n"], 6)
        self.assertEqual(cr["E"], 3)
        self.assertAlmostEqual(cr["rbar"], 5.0)
        self.assertAlmostEqual(cr["z_cl"], 0.0)
        self.assertAlmostEqual(cr["var_cl"], 5.333333333, places=6)
        self.assertAlmostEqual(cr["phi"], 4.8, places=6)
        self.assertAlmostEqual(J.z_iid([r for _, _, r in rows]), 0.0)

    def test_hand_example(self):
        # A [9,8], B [7,6], C [5,7]: n 6, rbar 7
        rows = rows_from({1: [9, 8], 2: [7, 6], 3: [5, 7]})
        self.assertAlmostEqual(J.z_iid([r for _, _, r in rows]), 1.897367, places=5)
        cr = J.cluster_robust(rows)
        # S = 3, -1, -2 -> 14; V = 1.5 * 14 / 36 = 0.583333
        self.assertAlmostEqual(cr["var_cl"], 0.5833333, places=6)
        self.assertAlmostEqual(cr["z_cl"], 2.618615, places=5)
        self.assertAlmostEqual(cr["phi"], 0.525, places=6)
        loo = J.leave_one_out(rows)
        self.assertEqual(loo["argmin"], 1)
        self.assertAlmostEqual(loo["min"], 0.968246, places=5)   # [7,6,5,7]
        self.assertEqual(loo["argmax"], 3)
        self.assertAlmostEqual(loo["max"], 1.936492, places=5)   # [9,8,7,6]
        self.assertAlmostEqual(loo["per_encounter"][2], 1.742842, places=5)
        lg = J.lag1(rows)
        # pairs (9,8), (7,6), (5,7): r = 2 / sqrt(8 * 2) = 0.5
        self.assertEqual(lg["pairs"], 3)
        self.assertAlmostEqual(lg["r"], 0.5, places=9)

    def test_lag1_needs_adjacency(self):
        rows = [(1, T0, 9.0), (1, T0 + 2 * STEP, 8.0), (2, T0, 7.0), (2, T0 + STEP, 6.0)]
        self.assertEqual(J.lag1(rows)["pairs"], 1)

    def test_per_encounter(self):
        pe = J.per_encounter(rows_from({4: [9, 8], 5: [1]}))
        self.assertAlmostEqual(pe[4]["mean_rank"], 8.5)
        self.assertEqual(pe[5]["n"], 1)
        self.assertAlmostEqual(pe[5]["z"], -4.0 / math.sqrt(80 / 12.0), places=9)


class TestGrade(unittest.TestCase):
    def test_not_gradable_small(self):
        rows = rows_from({1: [9, 8], 2: [7, 6], 3: [5, 7]})
        self.assertEqual(J.analyse(rows)["H2c_R"], "NOT GRADABLE")

    def test_robust_and_not_robust(self):
        cr = {"n": 3000, "E": 20, "z_cl": 3.2}
        self.assertEqual(J.grade(cr, {"min": 3.05}), "ROBUST")
        self.assertEqual(J.grade(cr, {"min": 2.9}), "NOT ROBUST")
        self.assertEqual(J.grade({"n": 3000, "E": 20, "z_cl": 2.99}, {"min": 3.5}), "NOT ROBUST")
        self.assertEqual(J.grade({"n": 3000, "E": 9, "z_cl": 5.0}, {"min": 5.0}), "NOT GRADABLE")
        self.assertEqual(J.grade({"n": 1999, "E": 20, "z_cl": 5.0}, {"min": 5.0}), "NOT GRADABLE")
        self.assertEqual(J.grade(cr, {"min": float("nan")}), "NOT GRADABLE")

    def test_small_thresholds_override(self):
        rows = rows_from({1: [9, 9, 9], 2: [9, 8, 9], 3: [9, 9, 8]})
        res = J.analyse(rows, min_n=1, min_e=2)
        self.assertEqual(res["H2c_R"], "ROBUST")


class TestCsv(unittest.TestCase):
    def test_load_filters(self):
        txt = ("enc,D,start,r_rs,rank_L,L_T,rank_L2,L2_T\n"
               "22,135,2024-12-24T00:00:00Z,10,5,1.0,7,1.1\n"
               "22,135,2024-12-24T02:15:00Z,10,5,1.0,,\n"
               "22,70,2024-12-24T00:00:00Z,10,5,1.0,3,1.1\n"
               "23,135,2025-03-22T00:00:00Z,10,5,1.0,nan,1.1\n"
               "23,135,2025-03-22T02:15:00Z,10,5,1.0,2,1.1\n")
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "w.csv")
            with open(p, "w") as fh:
                fh.write(txt)
            rows, skipped = J.load_rows(p)
        self.assertEqual(len(rows), 2)
        self.assertEqual(skipped, 2)
        self.assertEqual(rows[0][0], 22)
        self.assertEqual(rows[1][2], 2.0)


if __name__ == "__main__":
    unittest.main()
