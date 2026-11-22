#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. Written by AI agents working for MV2 of Munim, Inc.
"""psp5jack 0.1.0 - encounter-clustered robustness check for H2c (PREREG amendment 0.1.2).

Written 2026-11-22, BEFORE any output of psp5min-v1, psp5alf-v1 or psp5curv-v1
existed. It reads psp5curv's windows.csv (one row per tiled window) and asks
whether H2c's rank Z survives two things the plain formula ignores:

1. Windows are not independent. Z = (mean rank - 5) / sqrt(6.667 / n) treats
   ~3,000 adjacent 135-minute windows as independent draws. Windows from the
   same encounter share solar-wind streams, so ranks can be correlated within
   an encounter. We replace the iid variance by the cluster-robust (CR1)
   variance with one cluster per encounter:
       V = E/(E-1) * sum_e (sum_i (r_ei - rbar))^2 / n^2
       Z_cl = (rbar - 5) / sqrt(V),  phi = V / (6.667 / n)
2. One encounter can drive the result. We recompute the iid Z with each
   encounter left out and report the smallest and largest.
Reported only: per-encounter mean ranks and Z, and the lag-1 correlation of
ranks between time-adjacent windows of the same encounter.

Grade H2c-R (fixed in amendment 0.1.2):
  ROBUST      if Z_cl >= 3.0 and min leave-one-encounter-out Z >= 3.0
  NOT ROBUST  otherwise
  NOT GRADABLE if fewer than 2,000 windows or fewer than 10 encounters.
Pure Python (standard library only).
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from datetime import datetime

VERSION = "0.1.0"
NB = 9
D_GRADED = 135
MIN_N = 2000
MIN_E = 10
Z_ROBUST = 3.0


def parse_time(s: str) -> float:
    s = s.strip()
    if s.endswith("Z"):
        s = s[:-1] + "+00:00"
    return datetime.fromisoformat(s).timestamp()


def _num(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def load_rows(path: str, d_min: int = D_GRADED, col: str = "rank_L2"):
    """(enc, t_start, rank) for windows of length d_min with a finite rank."""
    out = []
    skipped = 0
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            d = _num(r.get("D"))
            if d is None or int(round(d)) != d_min:
                continue
            rk = _num(r.get(col))
            enc = _num(r.get("enc"))
            st = r.get("start", "")
            if rk is None or enc is None or not st:
                skipped += 1
                continue
            out.append((int(enc), parse_time(st), rk))
    return out, skipped


def z_iid(ranks, nb: int = NB) -> float:
    n = len(ranks)
    if n == 0:
        return float("nan")
    mu = (nb + 1) / 2.0
    var = (nb * nb - 1) / 12.0
    return (sum(ranks) / n - mu) / math.sqrt(var / n)


def by_encounter(rows):
    g = {}
    for enc, t, rk in rows:
        g.setdefault(enc, []).append((t, rk))
    for enc in g:
        g[enc].sort()
    return g


def cluster_robust(rows, nb: int = NB):
    """CR1 cluster-robust Z of the mean rank, clusters = encounters."""
    g = by_encounter(rows)
    n = len(rows)
    E = len(g)
    mu = (nb + 1) / 2.0
    var_iid = (nb * nb - 1) / 12.0 / n if n else float("nan")
    if n == 0 or E < 2:
        return {"n": n, "E": E, "rbar": float("nan"), "z_cl": float("nan"),
                "phi": float("nan"), "var_cl": float("nan")}
    rbar = sum(rk for _, _, rk in rows) / n
    ssum = 0.0
    for enc, lst in g.items():
        s = sum(rk - rbar for _, rk in lst)
        ssum += s * s
    V = E / (E - 1.0) * ssum / (n * n)
    z = (rbar - mu) / math.sqrt(V) if V > 0 else float("nan")
    return {"n": n, "E": E, "rbar": rbar, "z_cl": z, "phi": V / var_iid, "var_cl": V}


def leave_one_out(rows, nb: int = NB):
    g = by_encounter(rows)
    res = {}
    for enc in sorted(g):
        rest = [rk for e, _, rk in rows if e != enc]
        res[enc] = z_iid(rest, nb)
    if not res:
        return {"min": float("nan"), "max": float("nan"), "argmin": None, "argmax": None,
                "per_encounter": {}}
    amin = min(res, key=lambda e: res[e])
    amax = max(res, key=lambda e: res[e])
    return {"min": res[amin], "max": res[amax], "argmin": amin, "argmax": amax,
            "per_encounter": res}


def lag1(rows, d_min: int = D_GRADED, tol_s: float = 1.0):
    """Pearson correlation of ranks of time-adjacent windows (gap = d_min) within encounters."""
    g = by_encounter(rows)
    xs, ys = [], []
    for lst in g.values():
        for (t0, r0), (t1, r1) in zip(lst, lst[1:]):
            if abs((t1 - t0) - d_min * 60.0) <= tol_s:
                xs.append(r0)
                ys.append(r1)
    m = len(xs)
    if m < 3:
        return {"pairs": m, "r": float("nan")}
    mx, my = sum(xs) / m, sum(ys) / m
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sxx = sum((a - mx) ** 2 for a in xs)
    syy = sum((b - my) ** 2 for b in ys)
    r = sxy / math.sqrt(sxx * syy) if sxx > 0 and syy > 0 else float("nan")
    return {"pairs": m, "r": r}


def per_encounter(rows, nb: int = NB):
    g = by_encounter(rows)
    return {enc: {"n": len(lst), "mean_rank": sum(rk for _, rk in lst) / len(lst),
                  "z": z_iid([rk for _, rk in lst], nb)} for enc, lst in sorted(g.items())}


def grade(cr, loo, min_n: int = MIN_N, min_e: int = MIN_E, zr: float = Z_ROBUST) -> str:
    if cr["n"] < min_n or cr["E"] < min_e:
        return "NOT GRADABLE"
    z, lm = cr["z_cl"], loo["min"]
    if not (math.isfinite(z) and math.isfinite(lm)):
        return "NOT GRADABLE"
    return "ROBUST" if (z >= zr and lm >= zr) else "NOT ROBUST"


def analyse(rows, nb: int = NB, min_n: int = MIN_N, min_e: int = MIN_E):
    cr = cluster_robust(rows, nb)
    loo = leave_one_out(rows, nb)
    return {"version": VERSION, "n": cr["n"], "E": cr["E"], "rbar": cr["rbar"],
            "z_iid": z_iid([rk for _, _, rk in rows], nb), "z_cl": cr["z_cl"],
            "phi": cr["phi"], "loo_min": loo["min"], "loo_min_enc": loo["argmin"],
            "loo_max": loo["max"], "loo_max_enc": loo["argmax"],
            "lag1": lag1(rows), "per_encounter": per_encounter(rows, nb),
            "loo_per_encounter": loo["per_encounter"],
            "H2c_R": grade(cr, loo, min_n, min_e)}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--windows", required=True, help="psp5curv windows.csv")
    ap.add_argument("--out", default="")
    a = ap.parse_args(argv)
    rows, skipped = load_rows(a.windows)
    res = analyse(rows)
    res["skipped_rows_without_rank"] = skipped
    # also the subsets psp5curv reports, information only
    sub = {}
    for name, keep in (("A_E1_E24", lambda e: e <= 24), ("B_E25_E27", lambda e: e > 24)):
        rr = [x for x in rows if keep(x[0])]
        c = cluster_robust(rr)
        sub[name] = {"n": c["n"], "E": c["E"], "z_iid": z_iid([x[2] for x in rr]),
                     "z_cl": c["z_cl"], "phi": c["phi"]}
    res["subsets"] = sub
    txt = json.dumps(res, indent=1, default=str)
    if a.out:
        with open(a.out, "w") as fh:
            fh.write(txt)
    print(txt)
    print(f"PJACK: n={res['n']} E={res['E']} rbar={res['rbar']:.4f} Z_iid={res['z_iid']:.3f} "
          f"Z_cl={res['z_cl']:.3f} phi={res['phi']:.3f} LOO min={res['loo_min']:.3f} "
          f"(E{res['loo_min_enc']}) max={res['loo_max']:.3f} (E{res['loo_max_enc']}) "
          f"lag1 r={res['lag1']['r']:.3f} pairs={res['lag1']['pairs']} | H2c-R {res['H2c_R']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
