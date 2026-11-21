#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. Written by AI agents working for MV2 of Munim, Inc.
"""psp5alf 0.1.0 - is the Alfvenicity of the published Parker Solar Probe
5-minute wave trains unusual for near-Sun fluctuations?

Companion to psp5min 0.1.0 (same repository, imported unchanged).
Huang, Velli, Panasenco et al. (ApJL 999, 2026; arXiv:2511.10906) support
their two ~3.1-3.2 mHz peaks (E22 2024-12-24 09:50-12:05, E24 2025-06-19
08:15-09:25) with Alfvenicity: cross helicity sigma_c "close to unity" in
time and in frequency, |B| nearly constant, and power concentrated in B_N.
Near the Sun almost all fluctuations are Alfvenic, so this asks whether
those signatures distinguish the published peaks from (a) ordinary
fluctuations at nearby frequencies and (b) chance peaks of the paper's own
statistic at control frequencies. Rules and predictions: PREREG-ALF-v0.1.md.

Data (public, NASA CDAWeb): FIELDS PSP_FLD_L2_MAG_RTN_4_SA_PER_CYC and
SWEAP SPAN-i PSP_SWP_SPI_SF00_L3_MOM (DENS, VEL_RTN_SUN) via HAPI.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
import urllib.parse

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(_HERE), "psp5min"))  # repo layout
sys.path.insert(0, _HERE)                                             # job layout
import p5min as P  # noqa: E402

_TIMECHECK_DONE = False

VERSION = "0.1.0"
SPI_DATASET = "PSP_SWP_SPI_SF00_L3_MOM"
SPI_PARAMS = "QUALITY_FLAG,DENS,VEL_RTN_SUN"   # HAPI order must follow the dataset
ENC_FIRST, ENC_LAST = 17, 27
HALF_SPAN_DAYS_A = 3
R_MAX_RS = 20.0
DT2 = 2.0                  # s, common B/V grid (the paper shows B and V at 2 s)
V_MAX_GAP = 10.0           # s, max spacing of the SPAN-i samples bracketing a grid point
MAX_NAN_FRAC = 0.05
MAX_NAN_RUN2 = 30          # 2 s points = 60 s
ROLL_PTS = 301             # centred 10-minute rolling mean at 2 s (+-150 points)
NEAR_CONTROLS = (1.3, 2.1, 4.8)
S_PEAK = 4.0
MU0 = 4e-7 * math.pi
M_P = 1.67262192e-27
# Alfven speed in km/s for B in nT and n in cm^-3
VA_KMS = 1e-9 / math.sqrt(MU0 * 1e6 * M_P) / 1e3


def va_factor(n_cm3: float) -> float:
    """km/s per nT for proton density n (cm^-3)."""
    return VA_KMS / math.sqrt(n_cm3)


# ------------------------------------------------------------- SPAN-i (HAPI)
def parse_hapi_spi(text: str):
    """HAPI CSV rows: time, QUALITY_FLAG, DENS, VR, VT, VN -> arrays."""
    ts, rows = [], []
    for line in text.splitlines():
        line = line.strip()
        if not line or line[0] in "#{}\"" or line.startswith(","):
            continue
        parts = line.split(",")
        if len(parts) < 6:
            continue
        s = parts[0].rstrip("Z")
        try:
            vals = [float(x) for x in parts[1:6]]
        except ValueError:
            continue
        ts.append(s)
        rows.append(vals)
    if not ts:
        return np.zeros(0), np.zeros((0, 5))
    try:
        tt = np.array(ts, dtype="datetime64[ns]")
    except ValueError:
        return np.zeros(0), np.zeros((0, 5))
    t = (tt - np.datetime64("1970-01-01T00:00:00", "ns")).astype(np.int64) / 1e9
    a = np.asarray(rows, float)
    return t, a


def fetch_spi_day(day_u: float, data_dir: str, log):
    d = os.path.join(data_dir, "spi")
    os.makedirs(d, exist_ok=True)
    tag = P.iso(day_u)[:10].replace("-", "")
    cache = os.path.join(d, f"spi_{tag}.npz")
    if os.path.exists(cache):
        z = np.load(cache)
        return z["t"], z["a"]
    q = {"id": SPI_DATASET, "parameters": SPI_PARAMS,
         "time.min": P.iso(day_u), "time.max": P.iso(day_u + 86400), "format": "csv"}
    try:
        txt = P.http_get(P.HAPI + "?" + urllib.parse.urlencode(q), timeout=300).decode()
    except Exception as e:  # noqa: BLE001
        log(f"PALF_WARN: SPAN-i {tag} failed: {e!r}")
        return np.zeros(0), np.zeros((0, 5))
    t, a = parse_hapi_spi(txt)
    if len(t):
        np.savez(cache, t=t, a=a)
    else:
        log(f"PALF_WARN: SPAN-i {tag} returned no rows")
    return t, a


def interp_to_grid(ts, vals, grid, max_gap: float = V_MAX_GAP):
    """Linear interpolation of (N,k) samples onto grid; NaN where the
    bracketing samples are more than max_gap apart or out of range."""
    ts = np.asarray(ts, float)
    vals = np.asarray(vals, float)
    if vals.ndim == 1:
        vals = vals[:, None]
    ok = np.isfinite(ts) & np.all(np.isfinite(vals), axis=1)
    ts, vals = ts[ok], vals[ok]
    out = np.full((len(grid), vals.shape[1]), np.nan)
    if len(ts) < 2:
        return out
    o = np.argsort(ts, kind="stable")
    ts, vals = ts[o], vals[o]
    keep = np.r_[True, np.diff(ts) > 0]
    ts, vals = ts[keep], vals[keep]
    i = np.searchsorted(ts, grid, side="right")
    inr = (i >= 1) & (i < len(ts))
    j = np.clip(i - 1, 0, len(ts) - 1)
    exact = ts[j] == grid
    gap = np.full(len(grid), np.inf)
    gap[inr] = ts[i[inr]] - ts[i[inr] - 1]
    good = (inr & (gap <= max_gap)) | exact
    for c in range(vals.shape[1]):
        out[good, c] = np.interp(grid[good], ts, vals[:, c])
    return out


def block_mean(x, m: int):
    """Mean of consecutive blocks of m rows, ignoring NaN; NaN if a block
    has no finite value."""
    x = np.asarray(x, float)
    n = (x.shape[0] // m) * m
    y = x[:n].reshape(-1, m, x.shape[1])
    fin = np.isfinite(y)
    cnt = fin.sum(axis=1)
    s = np.where(fin, y, 0.0).sum(axis=1)
    out = np.full(cnt.shape, np.nan)
    np.divide(s, cnt, out=out, where=cnt > 0)
    return out


# ------------------------------------------------------------- statistics
def fill_gaps(seg, max_frac=MAX_NAN_FRAC, max_run=MAX_NAN_RUN2):
    seg = np.asarray(seg, float)
    nan = ~np.all(np.isfinite(seg), axis=1)
    if nan.mean() > max_frac or P.longest_run(nan) > max_run:
        return None
    x = seg.copy()
    if nan.any():
        idx = np.arange(len(x))
        for c in range(x.shape[1]):
            x[nan, c] = np.interp(idx[nan], idx[~nan], x[~nan, c])
    return x


def rolling_detrend(x, w: int = ROLL_PTS):
    """x minus its centred rolling mean over w points (shrinking at edges)."""
    x = np.asarray(x, float)
    h = w // 2
    n = len(x)
    c = np.vstack([np.zeros((1, x.shape[1])), np.cumsum(x, axis=0)])
    i = np.arange(n)
    lo = np.clip(i - h, 0, n)
    hi = np.clip(i + h + 1, 0, n)
    mean = (c[hi] - c[lo]) / (hi - lo)[:, None]
    return x - mean


def sigma_c(dv, db) -> float:
    dv = np.asarray(dv, float)
    db = np.asarray(db, float)
    den = float((dv ** 2).sum() + (db ** 2).sum())
    return 2.0 * float((dv * db).sum()) / den if den > 0 else float("nan")


def band_terms(v, b):
    """FFT of mean-removed v and b (N,3); returns per-bin cross (Re V.B*),
    energy (|V|^2+|B|^2), residual (|V|^2-|B|^2) and per-component |B|^2,
    for bins k = 1 .. N//2 - 1 (index k-1)."""
    v = np.asarray(v, float) - np.mean(v, axis=0)
    b = np.asarray(b, float) - np.mean(b, axis=0)
    n = v.shape[0]
    V = np.fft.rfft(v, axis=0)[1:n // 2]
    B = np.fft.rfft(b, axis=0)[1:n // 2]
    cross = np.real(V * np.conj(B)).sum(axis=1)
    ev = (np.abs(V) ** 2).sum(axis=1)
    eb = (np.abs(B) ** 2).sum(axis=1)
    return cross, ev + eb, ev - eb, np.abs(B) ** 2


def band_values(terms, ks):
    cross, en, res, bc = terms
    idx = np.asarray(ks) - 1
    E = float(en[idx].sum())
    btot = float(bc[idx].sum())
    return {
        "sc": 2.0 * float(cross[idx].sum()) / E if E > 0 else float("nan"),
        "sr": float(res[idx].sum()) / E if E > 0 else float("nan"),
        "bn_share": float(bc[idx, 2].sum()) / btot if btot > 0 else float("nan"),
    }


def frac_below(ref, x) -> float:
    if x is None or not math.isfinite(float(x)):
        return float("nan")
    ref = np.asarray(ref, float)
    ref = ref[np.isfinite(ref)]
    return float((ref < x).mean()) if len(ref) else float("nan")


def alf_window(b2, v2, n2, lay):
    """Alfvenicity of one window on the 2 s grid. Returns dict or None."""
    bx = fill_gaps(b2)
    vx = fill_gaps(v2)
    if bx is None or vx is None:
        return None
    nn = np.asarray(n2, float).ravel()
    nn = nn[np.isfinite(nn) & (nn > 0)]
    if len(nn) < 0.95 * len(n2):
        return None
    nmed = float(np.median(nn))
    bA = bx * va_factor(nmed)
    rec = {"n_med": round(nmed, 2),
           "sc_t": sigma_c(rolling_detrend(vx), rolling_detrend(bA))}
    bmag = np.sqrt((bx ** 2).sum(axis=1))
    rec["dBmag_over_B"] = float(bmag.std() / bmag.mean()) if bmag.mean() > 0 else float("nan")
    terms = band_terms(vx, bA)
    for name, ks in lay["paper"].items():
        bv = band_values(terms, ks)
        rec["sc_" + name] = bv["sc"]
        rec["sr_" + name] = bv["sr"]
        rec["bn_" + name] = bv["bn_share"]
    return rec


# ------------------------------------------------------------- pipeline
def encounter_grids(enc, peri_t, data_dir, log):
    global _TIMECHECK_DONE  # noqa: PLW0603
    g0 = math.floor((peri_t - HALF_SPAN_DAYS_A * 86400) / 60.0) * 60.0
    g1 = g0 + 2 * HALF_SPAN_DAYS_A * 86400
    nbin = int(round((g1 - g0) / P.DT))
    bdir = os.path.join(data_dir, "bins")
    os.makedirs(bdir, exist_ok=True)
    cache = os.path.join(bdir, f"E{enc:02d}_{int(g0)}_pm{HALF_SPAN_DAYS_A}.npy")
    meta = {}
    if os.path.exists(cache):
        B05 = np.load(cache)
        meta["mag_cached"] = True
    else:
        try:
            files = P.mag_file_list(g0, g1)
        except Exception as e:  # noqa: BLE001
            log(f"PALF_WARN: E{enc} CDAS list failed ({e!r}); trying direct URLs")
            files = P.direct_file_list(g0, g1)
        raw = os.path.join(data_dir, "raw")
        binner = P.Binner(g0, nbin)
        nf = 0
        for url, ln in files:
            try:
                p = P.download(url, ln, raw)
                ep, b = P.read_mag(p)
            except Exception as e:  # noqa: BLE001
                log(f"PALF_WARN: E{enc} MAG file failed {url.rsplit('/', 1)[-1]}: {e!r}")
                continue
            if len(ep) and not _TIMECHECK_DONE:
                d, s = P.tt2000_check(int(ep[0]))
                log(f"PALF_TIMECHECK: ours-minus-cdflib {d} (cdflib '{s}')")
                if d is not None and abs(d) > 0.5:
                    P._TIME_FIX = -round(d)
                    log(f"PALF_TIMECHECK: applying correction {P._TIME_FIX} s")
                _TIMECHECK_DONE = True
            binner.add(P.tt2000_to_unix(ep), b)
            nf += 1
            try:
                os.remove(p)
            except OSError:
                pass
        B05 = binner.result()
        if nf:
            np.save(cache, B05)
        meta.update({"mag_files": nf, "mag_listed": len(files)})
    m = int(round(DT2 / P.DT))
    B2 = block_mean(B05, m)
    grid = g0 + (np.arange(B2.shape[0]) + 0.5) * DT2
    ts, As = [], []
    day = math.floor(g0 / 86400) * 86400
    while day < g1:
        t, a = fetch_spi_day(day, data_dir, log)
        ts.append(t)
        As.append(a)
        day += 86400
    t = np.concatenate(ts) if ts else np.zeros(0)
    a = np.concatenate(As) if As else np.zeros((0, 5))
    meta["spi_samples"] = int(len(t))
    if len(t):
        meta["spi_median_dt"] = round(float(np.median(np.diff(np.sort(t)))), 4)
        qf = a[:, 0]
        u, cnt = np.unique(qf[np.isfinite(qf)], return_counts=True)
        top = np.argsort(-cnt)[:5]
        meta["spi_quality_flags_most_common"] = [[int(u[i]), int(cnt[i])] for i in top]
    V2 = interp_to_grid(t, a[:, 2:5], grid) if len(t) else np.full((len(grid), 3), np.nan)
    N2 = interp_to_grid(t, a[:, 1], grid)[:, 0] if len(t) else np.full(len(grid), np.nan)
    meta["mag_cover_0p5s"] = round(float(np.isfinite(B05[:, 0]).mean()), 4)
    meta["v_cover_2s"] = round(float(np.isfinite(V2[:, 0]).mean()), 4)
    return g0, B05, B2, V2, N2, meta


def window_record(enc, D, s, g0, B05, B2, V2, N2, lays, pos_t, pos_r, peri_t):
    n05 = int(round(D * 60 / P.DT))
    m = int(round(DT2 / P.DT))
    s2, n2 = s // m, n05 // m
    ts = g0 + s * P.DT
    te = ts + n05 * P.DT
    rec = {"enc": enc, "D": D, "s": s, "start": P.iso(ts),
           "r_rs": round(float(np.interp((ts + te) / 2, pos_t, pos_r)) * P.RS_PER_AU, 3),
           "hours_from_peri": round(((ts + te) / 2 - peri_t) / 3600.0, 3),
           "pub": P.overlaps_pub(ts, te)}
    rec.update(P.analyse_window(B05[s:s + n05], lays[D]))
    a = alf_window(B2[s2:s2 + n2], V2[s2:s2 + n2], N2[s2:s2 + n2], lays[D])
    rec["alf_ok"] = int(a is not None)
    if a:
        rec.update(a)
    return rec


def eligible(r):
    return (r.get("valid") == 1 and r.get("base_ok") == 1 and r.get("alf_ok") == 1
            and r.get("pub") == 0 and r.get("r_rs", 99) <= R_MAX_RS)


def summarise(rows, pubs):
    out = {"version": VERSION}
    step = int(round(P.STEP_MIN * 60 / P.DT))
    cpub = {}
    for p in pubs:
        cpub[p["enc"]] = abs(p.get("sc_P3.15", float("nan")))
    vals = list(cpub.values())
    cmin = (min(vals) if len(vals) == 2 and all(math.isfinite(v) for v in vals)
            else float("nan"))
    out["c_pub"] = cpub
    out["c_min"] = cmin
    # A1: published |sc| at 3.15 mHz vs near-control |sc| of eligible tiled windows
    a1 = {}
    for p in pubs:
        D = p["D"]
        ref = [abs(r["sc_C" + str(fc)]) for r in rows if r["D"] == D and r["tiled"] == 1
               and eligible(r) for fc in NEAR_CONTROLS]
        ref_all = [abs(r["sc_C" + str(fc)]) for r in rows if r["D"] == D and r["tiled"] == 1
                   and eligible(r) for fc in P.CONTROL_CENTRES]
        ref_p = [abs(r["sc_P3.15"]) for r in rows if r["D"] == D and r["tiled"] == 1
                 and eligible(r)]
        ref_t = [abs(r["sc_t"]) for r in rows if r["D"] == D and r["tiled"] == 1 and eligible(r)]
        a1[f"E{p['enc']}"] = {
            "D": D, "c_pub": cpub[p["enc"]], "n_ref": len(ref),
            "q_near": frac_below(ref, cpub[p["enc"]]),
            "q_all8": frac_below(ref_all, cpub[p["enc"]]),
            "q_paperband_other_windows": frac_below(ref_p, cpub[p["enc"]]),
            "median_ref_near": float(np.median(ref)) if ref else float("nan"),
            "sc_t_pub": p.get("sc_t", float("nan")),
            "q_sc_t": frac_below(ref_t, abs(p.get("sc_t", float("nan")))),
            "bn_pub": p.get("bn_P3.15", float("nan")),
            "q_bn_near": frac_below([r["bn_C" + str(fc)] for r in rows if r["D"] == D
                                     and r["tiled"] == 1 and eligible(r)
                                     for fc in NEAR_CONTROLS], p.get("bn_P3.15", float("nan"))),
            "dBmag_over_B_pub": p.get("dBmag_over_B", float("nan")),
        }
    out["A1"] = a1
    # A2: chance peaks (S >= 4 at a near-control band) - are they as Alfvenic?
    peaks = []
    for D in P.DURATIONS_MIN:
        for fc in NEAR_CONTROLS:
            key = "S_C" + str(fc)
            by_enc = {}
            for r in rows:
                if r["D"] == D and r["sliding"] == 1 and eligible(r) and r[key] >= S_PEAK:
                    by_enc.setdefault(r["enc"], []).append(r)
            for e, hits in sorted(by_enc.items()):
                for c in P.clusters_from_starts([h["s"] for h in hits], step):
                    cs = set(c)
                    best = max((h for h in hits if h["s"] in cs), key=lambda h: h[key])
                    peaks.append({"D": D, "band": "C" + str(fc), "enc": e,
                                  "start": best["start"], "S": round(best[key], 3),
                                  "r_rs": best["r_rs"], "n_windows": len(c),
                                  "abs_sc": abs(best["sc_C" + str(fc)]),
                                  "sc_t": best["sc_t"], "bn": best["bn_C" + str(fc)]})
    sc = np.asarray([p["abs_sc"] for p in peaks], float)
    ok_c = math.isfinite(cmin)
    a2 = {"n_peaks": len(peaks),
          "n_peaks_S_ge_5p8": sum(1 for p in peaks if p["S"] >= P.THRESH_PAPER),
          "frac_ge_cmin": float((sc >= cmin).mean()) if (len(sc) and ok_c) else float("nan"),
          "median_abs_sc": float(np.median(sc)) if len(sc) else float("nan"),
          "frac_ge_cmin_S_ge_5p8": (float(np.mean([p["abs_sc"] >= cmin for p in peaks
                                                   if p["S"] >= P.THRESH_PAPER]))
                                    if ok_c and any(p["S"] >= P.THRESH_PAPER for p in peaks)
                                    else float("nan"))}
    bns = [p["bn"] for p in peaks if math.isfinite(p["bn"])]
    a2["median_bn_share"] = float(np.median(bns)) if bns else float("nan")
    out["A2"] = a2
    out["peaks"] = peaks
    # info: paper-band peaks inside r <= 20 Rs (E17-E27), out of the published windows
    info = {}
    for D in P.DURATIONS_MIN:
        el = [r for r in rows if r["D"] == D and r["sliding"] == 1 and eligible(r)]
        d = {"eligible_sliding": len(el)}
        for bn in ["P3.15"] + ["C" + str(fc) for fc in P.CONTROL_CENTRES]:
            for thr in (P.THRESH_PAPER, S_PEAK):
                k = 0
                by = {}
                for r in el:
                    if r["S_" + bn] >= thr:
                        by.setdefault(r["enc"], []).append(r["s"])
                for e in by:
                    k += len(P.clusters_from_starts(by[e], step))
                d[f"{bn}@{thr}"] = k
        info[f"D{D}"] = d
    out["cluster_counts_r_le_20"] = info
    return out


def grade(res, pubs):
    g = {}
    r1 = [abs(p.get("sc_t", float("nan"))) >= 0.8 for p in pubs]
    g["R1"] = "PASS" if (len(r1) == 2 and all(r1)) else "FAIL"
    qs = [v["q_near"] for v in res["A1"].values()]
    ns = [v["n_ref"] for v in res["A1"].values()]
    if len(qs) != 2 or any(not math.isfinite(q) for q in qs) or min(ns) < 30:
        g["A1"] = "NOT GRADABLE"
    elif all(q < 0.90 for q in qs):
        g["A1"] = "NOT DISTINCTIVE"
    elif all(q >= 0.99 for q in qs):
        g["A1"] = "DISTINCTIVE"
    else:
        g["A1"] = "MIXED"
    a2 = res["A2"]
    if a2["n_peaks"] < 5 or not math.isfinite(a2["frac_ge_cmin"]):
        g["A2"] = "NOT GRADABLE"
    elif a2["frac_ge_cmin"] >= 0.5:
        g["A2"] = "AS ALFVENIC AS PUBLISHED"
    elif a2["frac_ge_cmin"] < 0.2:
        g["A2"] = "LESS ALFVENIC"
    else:
        g["A2"] = "INTERMEDIATE"
    if g["R1"] == "FAIL":
        g["note"] = "R1 FAIL: A1/A2 are graded 'with our reading of the paper's sigma_c'."
    return g


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--budget-s", type=float, default=12000)
    ap.add_argument("--only", default="")
    a = ap.parse_args(argv)
    t_start = time.time()
    os.makedirs(a.out_dir, exist_ok=True)

    def log(msg):
        print(msg, flush=True)

    log(f"PALF_VERSION: {VERSION} p5min {P.VERSION} numpy {np.__version__} "
        f"sha256 {P.sha256_file(os.path.abspath(__file__))} "
        f"p5min_sha256 {P.sha256_file(os.path.abspath(P.__file__))}")
    lays = {d: P.layout(d * 60) for d in P.DURATIONS_MIN}
    pos_t, pos_r = P.fetch_positions()
    peris = P.find_perihelia(pos_t, pos_r)
    allp = [(i + 1, t, r) for i, (t, r) in enumerate(peris)]
    num_ok = all(any(n == e and abs(t - P.to_unix(x)) <= 2 * 86400 for n, t, r in allp)
                 for e, x in P.EXPECTED_PERI.items())
    log(f"PALF_PERIHELIA: n={len(allp)} numbering_ok={num_ok}")
    encs = [(e, t) for e, t, r in allp if ENC_FIRST <= e <= ENC_LAST
            and t + HALF_SPAN_DAYS_A * 86400 <= P.to_unix(P.MAG_STOP)]
    if a.only:
        want = {int(x) for x in a.only.split(",") if x.strip()}
        encs = [x for x in encs if x[0] in want]
    encs.sort(key=lambda x: (0, -x[0]) if x[0] >= 22 else (1, -x[0]))
    rows, pubs, info = [], [], []
    step = int(round(P.STEP_MIN * 60 / P.DT))
    for enc, pt in encs:
        if time.time() - t_start > a.budget_s:
            log(f"PALF_BUDGET: stopping before E{enc}")
            break
        t0 = time.time()
        g0, B05, B2, V2, N2, meta = encounter_grids(enc, pt, a.data_dir, log)
        nbin = B05.shape[0]
        n_enc = 0
        for D in P.DURATIONS_MIN:
            n05 = int(round(D * 60 / P.DT))
            starts = set(range(0, nbin - n05 + 1, step)) | set(range(0, nbin - n05 + 1, n05))
            for s in sorted(starts):
                ts = g0 + s * P.DT
                rmid = float(np.interp(ts + n05 * P.DT / 2, pos_t, pos_r)) * P.RS_PER_AU
                if rmid > R_MAX_RS:
                    continue
                r = window_record(enc, D, s, g0, B05, B2, V2, N2, lays, pos_t, pos_r, pt)
                r["sliding"] = int(s % step == 0)
                r["tiled"] = int(s % n05 == 0)
                rows.append(r)
                n_enc += 1
        for p in P.PUBLISHED:
            if p["enc"] != enc:
                continue
            ta, tb = P.to_unix(p["start"]), P.to_unix(p["stop"])
            D = int(round((tb - ta) / 60))
            s = int(round((ta - g0) / P.DT))
            r = window_record(enc, D, s, g0, B05, B2, V2, N2, lays, pos_t, pos_r, pt)
            r.update({"published_f_mhz": p["f_mhz"], "published_sigma": p["sigma"]})
            pubs.append(r)
            log("PALF_PUBLISHED: " + json.dumps(r, default=P.jdefault))
        meta.update({"enc": enc, "perihelion": P.iso(pt), "windows": n_enc,
                     "alf_ok": sum(1 for r in rows if r["enc"] == enc and r["alf_ok"] == 1),
                     "seconds": round(time.time() - t0, 1)})
        info.append(meta)
        log("PALF_ENC: " + json.dumps(meta, default=P.jdefault))
        del B05, B2, V2, N2
    res = summarise(rows, pubs)
    res["numbering_ok"] = num_ok
    res["encounters"] = info
    res["published"] = pubs
    res["grades"] = grade(res, pubs)
    res["runtime_s"] = round(time.time() - t_start, 1)
    with open(os.path.join(a.out_dir, "results.json"), "w") as f:
        json.dump(res, f, indent=1, default=P.jdefault)
    P.write_rows(os.path.join(a.out_dir, "windows.csv"), rows)
    P.write_rows(os.path.join(a.out_dir, "peaks.csv"), res["peaks"])
    P.write_rows(os.path.join(a.out_dir, "published.csv"), pubs)
    log("PALF_A1: " + json.dumps(res["A1"], default=P.jdefault))
    log("PALF_A2: " + json.dumps(res["A2"], default=P.jdefault))
    log("PALF_COUNTS: " + json.dumps(res["cluster_counts_r_le_20"], default=P.jdefault))
    for p in res["peaks"]:
        log("PALF_PEAK: " + json.dumps(p, default=P.jdefault))
    log("PALF_GRADES: " + json.dumps(res["grades"], default=P.jdefault))
    return 0


if __name__ == "__main__":
    sys.exit(main())
