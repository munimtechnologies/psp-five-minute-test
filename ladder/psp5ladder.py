#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. Written by AI agents working for MV2 of Munim, Inc.
"""psp5ladder 0.1.0 - how much do the published "5.8 / 6.2 sigma" depend on
choices the paper does not state? (PREREG amendment 0.1.3.)

Claim: Huang, Velli et al., ApJL 999 L4 (2026), arXiv:2511.10906. Two
magnetic wave trains, E22 2024-12-24 09:50-12:05 UTC (3.2 mHz, 5.8 sigma)
and E24 2025-06-19 08:15-09:25 UTC (3.1 mHz, 6.2 sigma), where
sigma = (peak of P*f - mean) / std over the frequency range "containing the
central 50% of the total fluctuation energy". The paper does not state the
taper, detrending, the upper frequency limit of the "total" energy, whether
the peak bin sits inside the baseline, or the std convention.

For the two published intervals only, this script computes the statistic
under a FIXED list of readings (V0-V6, below) and with each interval edge
moved by -10, -5, +5 and +10 minutes. Rules and predictions are fixed in
PREREG-amendment-0.1.3.md before any data were read.

It imports psp5min 0.1.0 (p5min.py at commit a025a506) unchanged for data
access, binning and the validity rule. V0 is psp5min's own reading.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
_SIB = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "psp5min")
if os.path.isdir(_SIB):
    sys.path.append(_SIB)
import p5min  # noqa: E402

VERSION = "0.1.0"
P5MIN_SHA = "a025a506d2092e982e82e87c2c27e4f2a477fa98"
DT = p5min.DT
FPK_HZ = (2.5e-3, 4.0e-3)       # peak search range, as psp5min
R0_F_MHZ = (3.00, 3.35)         # psp5min R0 frequency rule
R0_TOL = 1.5                    # psp5min R0 sigma rule
MARGIN_S = 1800                 # data fetched around each interval
SHIFTS_MIN = (-10, -5, 5, 10)   # edge moves
F_MID_E22_MHZ = 1e3 / 105.0     # paper: 1/f_mid = 105 s for E22

# name -> options. Order matters (RL2 selection rule).
VARIANTS = (
    ("V0", {"detrend": False, "taper": "none", "measure": "P", "fmax": None,
            "ddof": 0, "excl": False},
     "psp5min reading: mean removed, no taper, energy = sum of P over all bins to "
     "Nyquist, 25-75 % first reach, peak bin inside the baseline, population std"),
    ("V1", {"detrend": False, "taper": "none", "measure": "P", "fmax": None,
            "ddof": 1, "excl": False},
     "as V0 with sample std (ddof 1)"),
    ("V2", {"detrend": False, "taper": "none", "measure": "P", "fmax": None,
            "ddof": 0, "excl": True},
     "as V0 with the peak bin and its two neighbours left out of mean and std"),
    ("V3", {"detrend": True, "taper": "none", "measure": "P", "fmax": None,
            "ddof": 0, "excl": False},
     "as V0 with a linear trend removed from each component"),
    ("V4", {"detrend": False, "taper": "hann", "measure": "P", "fmax": None,
            "ddof": 0, "excl": False},
     "as V0 with a Hann taper"),
    ("V5", {"detrend": False, "taper": "none", "measure": "P", "fmax": 0.1,
            "ddof": 0, "excl": False},
     "as V0 with the total energy counted only up to 0.1 Hz"),
    ("V6", {"detrend": False, "taper": "none", "measure": "Pf", "fmax": None,
            "ddof": 0, "excl": False},
     "cumulative of P*f over linear bins (a less defensible reading; reported, "
     "not used in RL1)"),
)
VNAMES = tuple(v[0] for v in VARIANTS)
VOPTS = {v[0]: v[1] for v in VARIANTS}
DEFENSIBLE = ("V0", "V1", "V2", "V3", "V4", "V5")
RL1_MIN_FINITE = 5
RL3_MIN_FINITE = 6


# ---------------------------------------------------------------- spectra
def prep(seg):
    """psp5min's validity rule and linear gap fill; None if the window fails."""
    seg = np.asarray(seg)
    if seg.ndim != 2 or seg.shape[0] < 64:
        return None
    nan = np.isnan(seg).any(axis=1)
    if nan.mean() > p5min.MAX_NAN_FRAC or p5min.longest_run(nan) > p5min.MAX_NAN_RUN:
        return None
    x = np.asarray(seg, np.float64).copy()
    if nan.any():
        idx = np.arange(len(x))
        good = ~nan
        for c in range(x.shape[1]):
            x[nan, c] = np.interp(idx[nan], idx[good], x[good, c])
    return x


def detrend_linear(x):
    """Remove a least-squares straight line from each column."""
    x = np.asarray(x, np.float64)
    n = x.shape[0]
    t = np.arange(n, dtype=np.float64)
    t = t - t.mean()
    out = x - x.mean(axis=0)
    slope = (t[:, None] * out).sum(axis=0) / (t * t).sum()
    return out - t[:, None] * slope[None, :]


def spec(x, taper: str = "none"):
    """Trace PSD, bins k = 1 .. N//2 - 1. With taper 'none' this is p5min.spectrum."""
    x = np.asarray(x, np.float64)
    n = x.shape[0]
    x = x - x.mean(axis=0)
    if taper == "hann":
        w = np.hanning(n)
        X = np.fft.rfft(x * w[:, None], axis=0)
        norm = float((w * w).sum())
    else:
        X = np.fft.rfft(x, axis=0)
        norm = float(n)
    P = (np.abs(X) ** 2).sum(axis=1) * (2.0 * DT / norm)
    k = np.arange(1, n // 2)
    return k / (n * DT), P[1:n // 2]


def baseline(P, f, y, measure: str = "P", fmax=None, ddof: int = 0, excl=None,
             min_bins: int = p5min.MIN_BASE_BINS):
    """Bins from the first 25 % to the first 75 % of cumulative energy.

    measure 'P' sums P over linear bins (energy); 'Pf' sums P*f.
    fmax limits the bins counted in the total. excl = (i0, i1) drops
    those 0-based bins from the mean and std (not from the edges)."""
    m = len(P)
    lim = m if fmax is None else int(np.searchsorted(f, fmax, side="right"))
    if lim < 2:
        return None
    w = np.asarray(P[:lim] if measure == "P" else y[:lim], np.float64)
    c = np.cumsum(w)
    if not c[-1] > 0:
        return None
    c = c / c[-1]
    lo = int(np.searchsorted(c, 0.25))
    hi = min(int(np.searchsorted(c, 0.75)), lim - 1)
    med = min(int(np.searchsorted(c, 0.5)), lim - 1)
    idx = np.arange(lo, hi + 1)
    nb_full = len(idx)
    if excl is not None:
        idx = idx[(idx < excl[0]) | (idx > excl[1])]
    rec = {"lo": lo, "hi": hi, "med": med, "nb": int(len(idx)), "nb_full": nb_full}
    if len(idx) < min_bins or (ddof >= len(idx)):
        rec["ok"] = False
        return rec
    base = np.asarray(y, np.float64)[idx]
    sd = float(base.std(ddof=ddof))
    rec.update({"ok": sd > 0, "mu": float(base.mean()), "sd": sd})
    return rec


def score(seg, vname: str):
    """The paper's statistic at the largest P*f bin in 2.5-4.0 mHz under one reading."""
    o = VOPTS[vname]
    x = prep(seg)
    if x is None:
        return {"variant": vname, "valid": 0}
    if o["detrend"]:
        x = detrend_linear(x)
    f, P = spec(x, o["taper"])
    y = P * f
    sel = np.flatnonzero((f >= FPK_HZ[0]) & (f <= FPK_HZ[1]))
    rec = {"variant": vname, "valid": 1, "n": int(x.shape[0])}
    if len(sel) == 0:
        rec["valid"] = 0
        return rec
    i = int(sel[np.argmax(y[sel])])
    rec["fpk_mhz"] = float(f[i]) * 1e3
    bl = baseline(P, f, y, o["measure"], o["fmax"], o["ddof"],
                  (i - 1, i + 1) if o["excl"] else None)
    rec["base_ok"] = int(bool(bl and bl.get("ok")))
    if bl:
        lo_f, hi_f = float(f[bl["lo"]]) * 1e3, float(f[bl["hi"]]) * 1e3
        rec.update({"base_lo_mhz": lo_f, "base_hi_mhz": hi_f,
                    "mid_arith_mhz": 0.5 * (lo_f + hi_f),
                    "mid_geo_mhz": math.sqrt(lo_f * hi_f) if lo_f > 0 else float("nan"),
                    "median_energy_mhz": float(f[bl["med"]]) * 1e3,
                    "nb": bl["nb"], "nb_full": bl["nb_full"]})
    if rec["base_ok"]:
        s = (float(y[i]) - bl["mu"]) / bl["sd"]
        inside = bool(bl["lo"] <= i <= bl["hi"])
        rec["S_fpk"] = s
        rec["peak_inside_base"] = int(inside)
        rec["samuelson_ceiling"] = (math.sqrt(bl["nb"] - 1)
                                    if inside and not o["excl"] and o["ddof"] == 0
                                    else float("nan"))
    else:
        rec["S_fpk"] = float("nan")
    return rec


# ---------------------------------------------------------------- grading
def _fin(v) -> bool:
    return isinstance(v, (int, float)) and math.isfinite(float(v))


def r0_pass(rec, published_sigma: float) -> bool:
    s = rec.get("S_fpk", float("nan"))
    fp = rec.get("fpk_mhz", float("nan"))
    return (_fin(s) and _fin(fp) and R0_F_MHZ[0] <= fp <= R0_F_MHZ[1]
            and abs(s - published_sigma) <= R0_TOL)


def grade(res):
    """res = {"events": {"E22": {"published_sigma", "valid_V0", "variants": {V: rec},
    "edges": [rec, ...]}, "E24": {...}}}"""
    g = {}
    ev = res["events"]
    names = sorted(ev)
    ok0 = (len(names) == 2 and all(
        ev[e]["variants"].get("V0", {}).get("base_ok") == 1 for e in names))
    g["RL0"] = "PASS" if ok0 else "FAIL"
    # RL1 reading spread over the defensible readings
    spreads, gradable = {}, ok0
    for e in names:
        vals = [ev[e]["variants"][v].get("S_fpk") for v in DEFENSIBLE
                if v in ev[e]["variants"]]
        vals = [float(v) for v in vals if _fin(v)]
        if len(vals) < RL1_MIN_FINITE:
            gradable = False
            spreads[e] = None
        else:
            spreads[e] = {"min": min(vals), "max": max(vals),
                          "spread": max(vals) - min(vals), "n": len(vals)}
    g["RL1_detail"] = spreads
    if not gradable:
        g["RL1"] = "NOT GRADABLE"
    elif any(spreads[e]["spread"] >= 2.0 for e in names):
        g["RL1"] = "FRAGILE"
    elif all(spreads[e]["spread"] < 1.0 for e in names):
        g["RL1"] = "STABLE"
    else:
        g["RL1"] = "MODERATE"
    # RL2 reproduction and the fixed selection rule
    passing = []
    worst = {}
    for v in VNAMES:
        recs = [ev[e]["variants"].get(v, {}) for e in names]
        if len(names) == 2 and all(r0_pass(ev[e]["variants"].get(v, {}),
                                           ev[e]["published_sigma"]) for e in names):
            passing.append(v)
        ss = [r.get("S_fpk") for r in recs]
        if len(names) == 2 and all(_fin(s) for s in ss):
            worst[v] = max(abs(float(ev[e]["variants"][v]["S_fpk"]) - ev[e]["published_sigma"])
                           for e in names)
    g["RL2_passing"] = passing
    g["RL2_worst_abs_dS"] = worst
    if not ok0:
        g["RL2"] = "NOT GRADABLE"
        g["RL2_selected"] = None
    elif passing:
        g["RL2"] = "REPRODUCED" if "V0" in passing else "REPRODUCED-OTHER-READING"
        g["RL2_selected"] = passing[0]
    else:
        g["RL2"] = "NOT REPRODUCED"
        g["RL2_selected"] = (min(worst, key=lambda v: (worst[v], VNAMES.index(v)))
                             if worst else None)
    # RL3 interval edges (V0)
    drops, gradable3 = {}, ok0
    for e in names:
        s0 = ev[e]["variants"].get("V0", {}).get("S_fpk")
        vals = [r.get("S_fpk") for r in ev[e].get("edges", [])]
        vals = [float(v) for v in vals if _fin(v)]
        if not _fin(s0) or len(vals) < RL3_MIN_FINITE:
            gradable3 = False
            drops[e] = None
        else:
            drops[e] = {"S_V0": float(s0), "min": min(vals), "max": max(vals),
                        "drop": float(s0) - min(vals), "n": len(vals)}
    g["RL3_detail"] = drops
    if not gradable3:
        g["RL3"] = "NOT GRADABLE"
    elif any(drops[e]["drop"] >= 2.0 for e in names):
        g["RL3"] = "EDGE-SENSITIVE"
    elif all(drops[e]["drop"] < 1.0 for e in names):
        g["RL3"] = "ROBUST"
    else:
        g["RL3"] = "PARTLY"
    # reported only: which midpoint convention matches the paper's 105 s
    e22 = ev.get("E22", {}).get("variants", {}).get("V0", {})
    mid = {}
    for key in ("mid_arith_mhz", "mid_geo_mhz", "median_energy_mhz"):
        v = e22.get(key)
        mid[key] = v
        mid[key + "_within_25pct"] = (bool(abs(v / F_MID_E22_MHZ - 1) <= 0.25)
                                      if _fin(v) else None)
    mid["paper_mhz"] = F_MID_E22_MHZ
    g["MID_E22_reported"] = mid
    return g


# ---------------------------------------------------------------- data
def load_interval(a: float, b: float, data_dir: str, log):
    g0 = a - MARGIN_S
    g1 = b + MARGIN_S
    nbin = int(round((g1 - g0) / DT))
    try:
        files = p5min.mag_file_list(g0, g1)
        src = "cdas"
    except Exception as e:  # noqa: BLE001
        log(f"LADDER_WARN: CDAS list failed ({e!r}); trying direct URLs")
        files = p5min.direct_file_list(g0, g1)
        src = "direct"
    binner = p5min.Binner(g0, nbin)
    nsamp, nfiles, checked = 0, 0, False
    names = []
    for url, ln in files:
        try:
            p = p5min.download(url, ln, os.path.join(data_dir, "raw"))
            ep, bv = p5min.read_mag(p)
        except Exception as e:  # noqa: BLE001
            log(f"LADDER_WARN: file failed {url.rsplit('/', 1)[-1]}: {e!r}")
            continue
        if len(ep) and not checked:
            d, s = p5min.tt2000_check(int(ep[0]))
            log(f"LADDER_TIMECHECK: ours-minus-cdflib {d} (cdflib '{s}')")
            if d is not None and abs(d) > 0.5:
                p5min._TIME_FIX = -round(d)
                log(f"LADDER_TIMECHECK: applying correction {p5min._TIME_FIX} s")
            checked = True
        nsamp += binner.add(p5min.tt2000_to_unix(ep), bv)
        nfiles += 1
        names.append(url.rsplit("/", 1)[-1])
    B = binner.result()
    return g0, B, {"source": src, "listed": len(files), "files": nfiles,
                   "file_names": names, "samples_in_range": nsamp,
                   "coverage_0p5s": float(np.isfinite(B[:, 0]).mean()) if len(B) else 0.0}


def segment(B, g0: float, a: float, b: float):
    s = int(round((a - g0) / DT))
    n = int(round((b - a) / DT))
    if s < 0 or s + n > B.shape[0]:
        return None
    return B[s:s + n]


def edge_windows(a: float, b: float):
    out = []
    for edge in ("start", "stop"):
        for sh in SHIFTS_MIN:
            if edge == "start":
                out.append((edge, sh, a + sh * 60.0, b))
            else:
                out.append((edge, sh, a, b + sh * 60.0))
    return out


# ---------------------------------------------------------------- output
def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def fmt(v):
    if isinstance(v, (bool, np.bool_)):
        return int(v)
    if isinstance(v, (int, np.integer)):
        return int(v)
    if isinstance(v, (float, np.floating)):
        v = float(v)
        return "nan" if not math.isfinite(v) else f"{v:.6g}"
    return v


def write_rows(path, rows):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: fmt(r.get(k, "")) for k in keys})


def jdefault(v):
    if isinstance(v, np.integer):
        return int(v)
    if isinstance(v, np.floating):
        return float(v)
    if isinstance(v, np.bool_):
        return bool(v)
    return str(v)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    a = ap.parse_args(argv)
    t0 = time.time()
    os.makedirs(a.out_dir, exist_ok=True)

    def log(msg):
        print(msg, flush=True)

    log(f"LADDER_VERSION: {VERSION} numpy {np.__version__} "
        f"sha256 {sha256_file(os.path.abspath(__file__))} "
        f"p5min {p5min.VERSION} sha256 {sha256_file(os.path.abspath(p5min.__file__))}")
    log("LADDER_VARIANTS: " + json.dumps({v[0]: v[2] for v in VARIANTS}))
    res = {"version": VERSION, "p5min_commit": P5MIN_SHA, "events": {}}
    vrows, erows = [], []
    for pub in p5min.PUBLISHED:
        name = f"E{pub['enc']}"
        pa, pb = p5min.to_unix(pub["start"]), p5min.to_unix(pub["stop"])
        g0, B, meta = load_interval(pa, pb, a.data_dir, log)
        log(f"LADDER_DATA_{name}: " + json.dumps(meta, default=jdefault))
        seg = segment(B, g0, pa, pb)
        evr = {"published_sigma": pub["sigma"], "published_f_mhz": pub["f_mhz"],
               "start": pub["start"], "stop": pub["stop"], "data": meta,
               "variants": {}, "edges": []}
        for v in VNAMES:
            rec = score(seg, v) if seg is not None else {"variant": v, "valid": 0}
            rec["r0_pass"] = int(r0_pass(rec, pub["sigma"]))
            evr["variants"][v] = rec
            vrows.append({"event": name, **rec})
            log(f"LADDER_{name}_{v}: " + json.dumps(rec, default=jdefault))
        for edge, sh, wa, wb in edge_windows(pa, pb):
            s2 = segment(B, g0, wa, wb)
            rec = score(s2, "V0") if s2 is not None else {"variant": "V0", "valid": 0}
            rec.update({"edge": edge, "shift_min": sh,
                        "start": p5min.iso(wa), "stop": p5min.iso(wb)})
            evr["edges"].append(rec)
            erows.append({"event": name, **rec})
            log(f"LADDER_{name}_EDGE: " + json.dumps(rec, default=jdefault))
        res["events"][name] = evr
        del B
    res["grades"] = grade(res)
    res["runtime_s"] = round(time.time() - t0, 1)
    with open(os.path.join(a.out_dir, "results.json"), "w") as f:
        json.dump(res, f, indent=1, default=jdefault)
    write_rows(os.path.join(a.out_dir, "variants.csv"), vrows)
    write_rows(os.path.join(a.out_dir, "edges.csv"), erows)
    log("LADDER_GRADES: " + json.dumps(res["grades"], default=jdefault))
    return 0


if __name__ == "__main__":
    sys.exit(main())
