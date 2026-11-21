#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. Written by AI agents working for MV2 of Munim, Inc.
"""psp5min 0.1.0 - pre-registered calibration of the in-situ 5-minute
oscillation claim from Parker Solar Probe.

Claim tested: Huang, Velli, Panasenco et al., "Evidence of 5 minute
Oscillations from Parker Solar Probe", ApJ Letters 999 (2026),
arXiv:2511.10906. Two magnetic-field wave trains at 9.9 solar radii
(E22 2024-12-24 09:50-12:05 UTC, 3.2 mHz, 5.8 sigma; E24 2025-06-19
08:15-09:25 UTC, 3.1 mHz, 6.2 sigma). Their "sigma" is
(peak of P(f)*f - mean) / std, with mean and std taken over the
frequency range holding the central 50 % of the fluctuation energy.
No false-alarm rate was given; the search covered +-7 days around
perihelion in encounters E1-E24.

This script runs the same statistic over every window of the authors'
search space (E1-E24, +-7 d) and over three perihelia they did not use
(E25-E27), at the paper's frequency and at control frequencies. The
questions, rules and predictions are fixed in PREREG-v0.1.md.

Data: NASA CDAWeb, dataset PSP_FLD_L2_MAG_RTN_4_SA_PER_CYC (FIELDS
fluxgate, RTN) and PSP_HELIO1HR_POSITION. Public, no login.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import sys
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone

import numpy as np

VERSION = "0.1.0"
DT = 0.5                      # s, resampling cadence used by the paper
RS_PER_AU = 215.032
TT2000_UNIX_OFFSET = 946727930.816   # valid after 2017-01-01 (TAI-UTC = 37 s)
UA = "psp5min/0.1 (open research script)"

MAG_DATASET = "PSP_FLD_L2_MAG_RTN_4_SA_PER_CYC"
MAG_VAR = "psp_fld_l2_mag_RTN_4_Sa_per_Cyc"
MAG_EPOCH_FALLBACK = "epoch_mag_RTN_4_Sa_per_Cyc"
CDAS = ("https://cdaweb.gsfc.nasa.gov/WS/cdasr/1/dataviews/sp_phys/datasets/"
        "{ds}/orig_data/{a},{b}")
DIRECT = ("https://cdaweb.gsfc.nasa.gov/sp_phys/data/psp/fields/l2/"
          "mag_rtn_4_per_cycle/{y}/psp_fld_l2_mag_rtn_4_sa_per_cyc_{ymd}_v{vv}.cdf")
HAPI = "https://cdaweb.gsfc.nasa.gov/hapi/data"
POS_DATASET = "PSP_HELIO1HR_POSITION"
POS_START = "2018-08-13T00:00:00Z"
POS_STOP = "2026-06-01T00:00:00Z"
MAG_STOP = "2026-04-30T23:59:59Z"   # HAPI stopDate seen 2026-11-21

HALF_SPAN_DAYS = 7
DURATIONS_MIN = (135, 70)     # the two published interval lengths
STEP_MIN = 10
THRESH_PAPER = 5.8            # the smaller of the two published values
THRESH_LOW = 4.0
MAX_NAN_FRAC = 0.05
MAX_NAN_RUN = 120             # 0.5 s bins = 60 s
MIN_BASE_BINS = 10
PERI_THRESH_AU = 0.25

PAPER_CENTRE = 3.15           # mHz, published peaks at 3.1 and 3.2 mHz
TARGET_CENTRE = 3.2           # mHz, local-ratio target band centre
CONTROL_CENTRES = (1.3, 2.1, 4.8, 6.6, 7.4, 8.2, 9.0, 9.8)
FLANK_CENTRES = (2.1, 4.8)

PUBLISHED = (
    {"enc": 22, "start": "2024-12-24T09:50:00Z", "stop": "2024-12-24T12:05:00Z",
     "f_mhz": 3.2, "sigma": 5.8},
    {"enc": 24, "start": "2025-06-19T08:15:00Z", "stop": "2025-06-19T09:25:00Z",
     "f_mhz": 3.1, "sigma": 6.2},
)
AUTHOR_LAST_ENC = 24
EXPECTED_PERI = {22: "2024-12-24T12:00:00Z", 24: "2025-06-19T03:00:00Z"}

_TIME_FIX = 0.0


# ---------------------------------------------------------------- utilities
def to_unix(s: str) -> float:
    s = s.strip()
    if s.endswith("Z"):
        s = s[:-1] + "+00:00"
    d = datetime.fromisoformat(s)
    if d.tzinfo is None:
        d = d.replace(tzinfo=timezone.utc)
    return d.timestamp()


def iso(u: float) -> str:
    return datetime.fromtimestamp(u, timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def cdas_stamp(u: float) -> str:
    return datetime.fromtimestamp(u, timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def http_get(url: str, accept: str | None = None, tries: int = 4,
             timeout: int = 180) -> bytes:
    last = None
    for i in range(tries):
        try:
            h = {"User-Agent": UA}
            if accept:
                h["Accept"] = accept
            req = urllib.request.Request(url, headers=h)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(5 * (i + 1))
    raise last  # type: ignore[misc]


def sha256_file(p: str) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


# ------------------------------------------------------------- positions
def parse_hapi_csv(text: str):
    ts, vs = [], []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(",")
        if len(parts) < 2:
            continue
        try:
            t = to_unix(parts[0])
            v = float(parts[1])
        except ValueError:
            continue
        if not math.isfinite(v) or v < -1e30 or v <= 0:
            continue
        ts.append(t)
        vs.append(v)
    return np.asarray(ts, float), np.asarray(vs, float)


def fetch_positions(start: str = POS_START, stop: str = POS_STOP):
    a, b = to_unix(start), to_unix(stop)
    out_t, out_v = [], []
    cur = a
    while cur < b:
        nxt = min(b, cur + 366 * 86400)
        q = {"id": POS_DATASET, "parameters": "RAD_AU", "start": iso(cur),
             "stop": iso(nxt), "format": "csv"}
        try:
            txt = http_get(HAPI + "?" + urllib.parse.urlencode(q)).decode()
        except Exception:  # noqa: BLE001  (HAPI 2 parameter names)
            q2 = {"id": POS_DATASET, "parameters": "RAD_AU",
                  "time.min": iso(cur), "time.max": iso(nxt), "format": "csv"}
            txt = http_get(HAPI + "?" + urllib.parse.urlencode(q2)).decode()
        t, v = parse_hapi_csv(txt)
        out_t.append(t)
        out_v.append(v)
        cur = nxt
    t = np.concatenate(out_t)
    v = np.concatenate(out_v)
    o = np.argsort(t, kind="stable")
    t, v = t[o], v[o]
    keep = np.r_[True, np.diff(t) > 0]
    return t[keep], v[keep]


def find_perihelia(t, r_au, thresh: float = PERI_THRESH_AU,
                   merge_days: float = 20.0):
    """One perihelion per contiguous run of r < thresh (argmin of r)."""
    below = r_au < thresh
    peris = []
    n = len(t)
    i = 0
    while i < n:
        if below[i]:
            j = i
            while j < n and below[j]:
                j += 1
            k = i + int(np.argmin(r_au[i:j]))
            peris.append([float(t[k]), float(r_au[k])])
            i = j
        else:
            i += 1
    merged = []
    for p in peris:
        if merged and p[0] - merged[-1][0] < merge_days * 86400:
            if p[1] < merged[-1][1]:
                merged[-1] = p
        else:
            merged.append(p)
    return [(a, b) for a, b in merged]


# ------------------------------------------------------------- magnetometer
def mag_file_list(a: float, b: float):
    url = CDAS.format(ds=MAG_DATASET, a=cdas_stamp(a), b=cdas_stamp(b))
    js = json.loads(http_get(url, accept="application/json").decode())
    fd = js.get("FileDescription", []) if isinstance(js, dict) else []
    out = []
    for d in fd:
        name = d.get("Name")
        if not name:
            continue
        try:
            ln = int(d.get("Length", 0) or 0)
        except (TypeError, ValueError):
            ln = 0
        out.append((name, ln))
    return out


def direct_file_list(a: float, b: float):
    out = []
    day = math.floor(a / 86400) * 86400
    while day < b:
        d = datetime.fromtimestamp(day, timezone.utc)
        for vv in ("02", "03", "01", "04"):
            url = DIRECT.format(y=d.strftime("%Y"), ymd=d.strftime("%Y%m%d"), vv=vv)
            try:
                req = urllib.request.Request(url, method="HEAD",
                                             headers={"User-Agent": UA})
                with urllib.request.urlopen(req, timeout=60) as r:
                    ln = int(r.headers.get("Content-Length", "0") or 0)
                out.append((url, ln))
                break
            except Exception:  # noqa: BLE001
                continue
        day += 86400
    return out


def download(url: str, length: int, dest: str) -> str:
    os.makedirs(dest, exist_ok=True)
    p = os.path.join(dest, url.rsplit("/", 1)[-1])
    if os.path.exists(p) and (length <= 0 or os.path.getsize(p) == length):
        return p
    tmp = p + ".part"
    last = None
    for i in range(4):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=300) as r, open(tmp, "wb") as f:
                shutil.copyfileobj(r, f, 1 << 20)
            if length > 0 and os.path.getsize(tmp) != length:
                raise IOError(f"size {os.path.getsize(tmp)} != {length}")
            os.replace(tmp, p)
            return p
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(10 * (i + 1))
    raise last  # type: ignore[misc]


def read_mag(path: str):
    import cdflib  # noqa: PLC0415
    c = cdflib.CDF(path)
    dep = None
    try:
        att = c.varattsget(MAG_VAR)
        dep = att.get("DEPEND_0")
        if isinstance(dep, (list, tuple, np.ndarray)):
            dep = dep[0]
        if isinstance(dep, bytes):
            dep = dep.decode()
    except Exception:  # noqa: BLE001
        dep = None
    ep = None
    for nm in [x for x in (dep, MAG_EPOCH_FALLBACK) if x]:
        try:
            ep = c.varget(nm)
            if ep is not None:
                break
        except Exception:  # noqa: BLE001
            continue
    try:
        b = c.varget(MAG_VAR)
    except Exception:  # noqa: BLE001
        b = None
    if ep is None or b is None:
        return np.zeros(0, np.int64), np.zeros((0, 3))
    ep = np.asarray(ep).astype(np.int64).ravel()
    b = np.asarray(b, dtype=float).reshape(-1, 3)
    m = min(len(ep), len(b))
    return ep[:m], b[:m]


def tt2000_to_unix(ep):
    return np.asarray(ep, dtype=np.float64) / 1e9 + TT2000_UNIX_OFFSET + _TIME_FIX


def tt2000_check(ep0: int):
    """Compare our TT2000 conversion with cdflib's for one value."""
    try:
        import cdflib  # noqa: PLC0415
        s = cdflib.cdfepoch.encode(int(ep0))
        if isinstance(s, (list, tuple, np.ndarray)):
            s = s[0]
        s = str(s)
        m = re.match(r"(\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d)(?:\.(\d+))?", s)
        if not m:
            return None, s
        u_ref = to_unix(m.group(1) + "+00:00")
        if m.group(2):
            u_ref += float("0." + m.group(2))
    except Exception as e:  # noqa: BLE001
        return None, repr(e)
    u_ours = float(ep0) / 1e9 + TT2000_UNIX_OFFSET
    return u_ours - u_ref, s


class Binner:
    """Mean of each component in 0.5 s bins on a fixed grid."""

    def __init__(self, g0: float, nbin: int):
        self.g0 = g0
        self.nbin = nbin
        self.cnt = np.zeros(nbin, np.int64)
        self.s = np.zeros((3, nbin))

    def add(self, t, b):
        t = np.asarray(t, float)
        b = np.asarray(b, float).reshape(-1, 3)
        ok = (np.isfinite(t) & np.all(np.isfinite(b), axis=1)
              & np.all(np.abs(b) < 1e29, axis=1))
        idx = np.floor((t[ok] - self.g0) / DT).astype(np.int64)
        bb = b[ok]
        inr = (idx >= 0) & (idx < self.nbin)
        idx, bb = idx[inr], bb[inr]
        if len(idx) == 0:
            return 0
        self.cnt += np.bincount(idx, minlength=self.nbin)
        for c in range(3):
            self.s[c] += np.bincount(idx, weights=bb[:, c], minlength=self.nbin)
        return int(len(idx))

    def result(self):
        out = np.full((self.nbin, 3), np.nan, dtype=np.float32)
        m = self.cnt > 0
        for c in range(3):
            out[m, c] = self.s[c][m] / self.cnt[m]
        return out


# ------------------------------------------------------------- spectra
def longest_run(mask) -> int:
    mask = np.asarray(mask, bool)
    if not mask.any():
        return 0
    m = np.concatenate(([0], mask.astype(np.int8), [0]))
    d = np.diff(m)
    starts = np.flatnonzero(d == 1)
    ends = np.flatnonzero(d == -1)
    return int((ends - starts).max())


def spectrum(x):
    """Trace PSD of an (N,3) series, bins k = 1 .. N//2 - 1 (no taper)."""
    x = np.asarray(x, np.float64)
    n = x.shape[0]
    x = x - x.mean(axis=0)
    X = np.fft.rfft(x, axis=0)
    P = (np.abs(X) ** 2).sum(axis=1) * (2.0 * DT / n)
    k = np.arange(1, n // 2)
    return k / (n * DT), P[1:n // 2]


def nearest_bins(fc_mhz: float, dur_s: float, n: int):
    kc = fc_mhz * 1e-3 * dur_s
    lo = max(1, int(math.floor(kc)) - n - 1)
    cand = list(range(lo, int(math.floor(kc)) + n + 2))
    cand.sort(key=lambda k: (abs(k - kc), k))
    return sorted(cand[:n])


def layout(dur_s: float):
    n_p = max(1, int(round(0.2e-3 * dur_s)))
    h = max(1, int(round(0.3e-3 * dur_s)))
    n_l = 2 * h + 1
    paper = {"P3.15": nearest_bins(PAPER_CENTRE, dur_s, n_p)}
    local = {"T3.2": nearest_bins(TARGET_CENTRE, dur_s, n_l)}
    for fc in CONTROL_CENTRES:
        paper[f"C{fc}"] = nearest_bins(fc, dur_s, n_p)
        local[f"C{fc}"] = nearest_bins(fc, dur_s, n_l)
    for group in (paper, local):
        seen = set()
        for ks in group.values():
            if seen & set(ks):
                raise ValueError(f"overlapping bands at duration {dur_s}")
            seen |= set(ks)
    return {"paper": paper, "local": local, "n_paper": n_p, "n_local": n_l}


def paper_baseline(P, y, min_bins: int = MIN_BASE_BINS):
    """mean/std of y = P*f over the bins holding the central 50 % of energy."""
    c = np.cumsum(P)
    if c[-1] <= 0:
        return None
    c = c / c[-1]
    lo = int(np.searchsorted(c, 0.25))
    hi = int(np.searchsorted(c, 0.75))
    hi = min(hi, len(y) - 1)
    if hi - lo + 1 < min_bins:
        return {"lo": lo, "hi": hi, "ok": False}
    base = y[lo:hi + 1]
    sd = float(base.std())
    return {"lo": lo, "hi": hi, "ok": sd > 0, "mu": float(base.mean()), "sd": sd}


def local_ratio(y, k: int) -> float:
    """y at bin k (1-based) over the median of y at 2 <= |j-k| <= w."""
    w = max(6, int(math.ceil(0.3 * k)))
    nmax = len(y)
    js = [j for j in range(k - w, k - 1) if j >= 1] + \
         [j for j in range(k + 2, k + w + 1) if j <= nmax]
    ref = float(np.median(y[np.asarray(js) - 1]))
    return float(y[k - 1] / ref) if ref > 0 else float("nan")


def analyse_window(seg, lay):
    nan = np.isnan(seg).any(axis=1)
    frac = float(nan.mean())
    run = longest_run(nan)
    rec = {"nan_frac": round(frac, 5), "nan_run": run, "valid": 0, "base_ok": 0}
    if frac > MAX_NAN_FRAC or run > MAX_NAN_RUN:
        return rec
    x = np.asarray(seg, np.float64).copy()
    if nan.any():
        idx = np.arange(len(x))
        good = ~nan
        for c in range(3):
            x[nan, c] = np.interp(idx[nan], idx[good], x[good, c])
    f, P = spectrum(x)
    y = P * f
    rec["valid"] = 1
    bl = paper_baseline(P, y)
    ok = bool(bl and bl["ok"])
    rec["base_ok"] = int(ok)
    if bl:
        rec["base_lo_mhz"] = round(float(f[bl["lo"]]) * 1e3, 4)
        rec["base_hi_mhz"] = round(float(f[min(bl["hi"], len(f) - 1)]) * 1e3, 4)
    for name, ks in lay["paper"].items():
        rec["S_" + name] = (max((y[k - 1] - bl["mu"]) / bl["sd"] for k in ks)
                            if ok else float("nan"))
    lvals = []
    for name, ks in lay["local"].items():
        v = max(local_ratio(y, k) for k in ks)
        rec["L_" + name] = v
        lvals.append(v)
    rec["rank_T"] = 1 + sum(1 for v in lvals[1:] if v < lvals[0])
    sel = np.flatnonzero((f >= 2.5e-3) & (f <= 4.0e-3))
    i = int(sel[np.argmax(y[sel])])
    rec["fpk_mhz"] = round(float(f[i]) * 1e3, 4)
    rec["S_fpk"] = (y[i] - bl["mu"]) / bl["sd"] if ok else float("nan")
    return rec


def clusters_from_starts(starts, step: int):
    """Group sorted window starts into runs of consecutive grid steps."""
    out = []
    for s in sorted(starts):
        if out and s - out[-1][-1] == step:
            out[-1].append(s)
        else:
            out.append([s])
    return out


def rank_z(ranks, nb: int = 9):
    ranks = np.asarray(ranks, float)
    n = len(ranks)
    if n == 0:
        return {"n": 0, "mean_rank": float("nan"), "z": float("nan")}
    mu = (nb + 1) / 2.0
    var = (nb * nb - 1) / 12.0
    z = (ranks.mean() - mu) / math.sqrt(var / n)
    return {"n": n, "mean_rank": float(ranks.mean()), "z": float(z)}


def binom_sf(k: int, n: int, p: float) -> float:
    """P(X >= k) for X ~ Binomial(n, p)."""
    if k <= 0:
        return 1.0
    tot = 0.0
    for j in range(k, n + 1):
        tot += math.comb(n, j) * p ** j * (1 - p) ** (n - j)
    return min(1.0, tot)


# ------------------------------------------------------------- pipeline
def encounter_bins(enc, peri_t, data_dir, keep_raw, log):
    g0 = math.floor((peri_t - HALF_SPAN_DAYS * 86400) / 60.0) * 60.0
    g1 = g0 + 2 * HALF_SPAN_DAYS * 86400
    nbin = int(round((g1 - g0) / DT))
    bdir = os.path.join(data_dir, "bins")
    os.makedirs(bdir, exist_ok=True)
    cache = os.path.join(bdir, f"E{enc:02d}_{int(g0)}.npy")
    if os.path.exists(cache):
        return g0, np.load(cache), {"cached": True}
    try:
        files = mag_file_list(g0, g1)
        src = "cdas"
    except Exception as e:  # noqa: BLE001
        log(f"P5MIN_WARN: E{enc} CDAS list failed ({e!r}); trying direct URLs")
        files = direct_file_list(g0, g1)
        src = "direct"
    raw = os.path.join(data_dir, "raw")
    binner = Binner(g0, nbin)
    nsamp, nfiles = 0, 0
    global _TIME_FIX  # noqa: PLW0603
    for url, ln in files:
        try:
            p = download(url, ln, raw)
            ep, b = read_mag(p)
        except Exception as e:  # noqa: BLE001
            log(f"P5MIN_WARN: E{enc} file failed {url.rsplit('/', 1)[-1]}: {e!r}")
            continue
        if len(ep) and enc == 22 and nfiles == 0:
            d, s = tt2000_check(int(ep[0]))
            log(f"P5MIN_TIMECHECK: ours-minus-cdflib {d} (cdflib '{s}')")
            if d is not None and abs(d) > 0.5:
                _TIME_FIX = -round(d)
                log(f"P5MIN_TIMECHECK: applying correction {_TIME_FIX} s")
        nsamp += binner.add(tt2000_to_unix(ep), b)
        nfiles += 1
        if not keep_raw:
            try:
                os.remove(p)
            except OSError:
                pass
    B = binner.result()
    if nfiles:
        np.save(cache, B)
    return g0, B, {"cached": False, "files": nfiles, "listed": len(files),
                   "samples": nsamp, "source": src}


_PUB_IV = [(to_unix(p["start"]), to_unix(p["stop"])) for p in PUBLISHED]


def overlaps_pub(ts: float, te: float):
    for a, b in _PUB_IV:
        if ts < b and te > a:
            return 1
    return 0


def run_encounter(enc, peri_t, B, g0, pos_t, pos_r, lays):
    rows = []
    nbin = B.shape[0]
    step = int(round(STEP_MIN * 60 / DT))
    for dmin in DURATIONS_MIN:
        n = int(round(dmin * 60 / DT))
        lay = lays[dmin]
        starts = set(range(0, nbin - n + 1, step)) | set(range(0, nbin - n + 1, n))
        for s in sorted(starts):
            ts = g0 + s * DT
            te = ts + n * DT
            rec = analyse_window(B[s:s + n], lay)
            rec.update({
                "enc": enc, "D": dmin, "s": s, "start": iso(ts),
                "sliding": int(s % step == 0), "tiled": int(s % n == 0),
                "r_rs": round(float(np.interp(ts + n * DT / 2, pos_t, pos_r)) * RS_PER_AU, 3),
                "pub": overlaps_pub(ts, te),
                "hours_from_peri": round((ts + n * DT / 2 - peri_t) / 3600.0, 3),
            })
            rows.append(rec)
    return rows


def published_windows(enc, B, g0, lays):
    out = []
    for p in PUBLISHED:
        if p["enc"] != enc:
            continue
        a, b = to_unix(p["start"]), to_unix(p["stop"])
        dmin = int(round((b - a) / 60))
        s = int(round((a - g0) / DT))
        n = int(round((b - a) / DT))
        rec = analyse_window(B[s:s + n], lays[dmin])
        rec.update({"enc": enc, "D": dmin, "start": p["start"], "stop": p["stop"],
                    "published_f_mhz": p["f_mhz"], "published_sigma": p["sigma"]})
        out.append(rec)
    return out


def jdefault(v):
    if isinstance(v, np.integer):
        return int(v)
    if isinstance(v, np.floating):
        return float(v)
    if isinstance(v, np.ndarray):
        return v.tolist()
    return str(v)


def fmt(v):
    if isinstance(v, np.integer):
        return int(v)
    if isinstance(v, (float, np.floating)):
        v = float(v)
        if not math.isfinite(v):
            return "nan"
        return f"{v:.6g}"
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


def eligible(r):
    return r.get("valid") == 1 and r.get("base_ok") == 1 and r.get("pub") == 0


def summarise(rows, pubrecs, encs_done, lays, enc_info):
    res = {"version": VERSION, "encounters_done": sorted(encs_done)}
    step = int(round(STEP_MIN * 60 / DT))
    paper_bands = list(lays[DURATIONS_MIN[0]]["paper"].keys())
    clusters = []
    for dmin in DURATIONS_MIN:
        dres = {}
        for samp, lo, hi in (("A_E1_E24", 1, AUTHOR_LAST_ENC),
                             ("B_E25_E27", AUTHOR_LAST_ENC + 1, 99)):
            sl = [r for r in rows if r["D"] == dmin and r["sliding"] == 1
                  and lo <= r["enc"] <= hi]
            el = [r for r in sl if eligible(r)]
            n_eff = len(el) * STEP_MIN / dmin
            by_enc = {}
            for r in el:
                by_enc.setdefault(r["enc"], []).append(r)
            band = {}
            for thr in (THRESH_PAPER, THRESH_LOW):
                for bn in paper_bands:
                    key = "S_" + bn
                    k = 0
                    for e in sorted(by_enc):
                        hits = [r for r in by_enc[e] if r[key] >= thr]
                        cl = clusters_from_starts([r["s"] for r in hits], step)
                        k += len(cl)
                        if thr == THRESH_PAPER and bn in ("P3.15", "C2.1", "C4.8"):
                            for c in cl:
                                cs = set(c)
                                best = max((r for r in hits if r["s"] in cs),
                                           key=lambda r: r[key])
                                clusters.append({"D": dmin, "band": bn, "enc": e,
                                                 "first": iso(to_unix(best["start"]) - (best["s"] - c[0]) * DT),
                                                 "n_windows": len(c), "peak_start": best["start"],
                                                 "peak_score": round(best[key], 3),
                                                 "r_rs": best["r_rs"]})
                    band[f"{bn}@{thr}"] = k
            kF = band[f"C2.1@{THRESH_PAPER}"] + band[f"C4.8@{THRESH_PAPER}"]
            ctrl = [band[f"C{fc}@{THRESH_PAPER}"] for fc in CONTROL_CENTRES]
            dres[samp] = {"sliding_windows": len(sl), "eligible": len(el),
                          "valid_frac": (sum(1 for r in sl if r["valid"] == 1) / len(sl)) if sl else float("nan"),
                          "n_eff_independent": round(n_eff, 1),
                          "clusters": band, "kF_flanks": kF,
                          "kF_per_band": kF / 2.0,
                          "k_controls_mean": float(np.mean(ctrl)),
                          "kT_paper": band[f"P3.15@{THRESH_PAPER}"]}
        # rank test on tiled windows (local ratio statistic)
        rk = {}
        tl = [r for r in rows if r["D"] == dmin and r["tiled"] == 1
              and r["valid"] == 1 and r["pub"] == 0]
        for name, sel in (("all", tl),
                          ("r_le_30", [r for r in tl if r["r_rs"] <= 30]),
                          ("A_E1_E24", [r for r in tl if r["enc"] <= AUTHOR_LAST_ENC]),
                          ("B_E25_E27", [r for r in tl if r["enc"] > AUTHOR_LAST_ENC])):
            rk[name] = rank_z([r["rank_T"] for r in sel])
        # tail test: pooled control 99th percentile of L
        ctrl_vals = np.asarray([r["L_C" + str(fc)] for r in tl for fc in CONTROL_CENTRES], float)
        ctrl_vals = ctrl_vals[np.isfinite(ctrl_vals)]
        if len(ctrl_vals):
            thr = float(np.percentile(ctrl_vals, 99))
            eT = sum(1 for r in tl if r["L_T3.2"] >= thr)
            eC = sum(1 for r in tl for fc in CONTROL_CENTRES if r["L_C" + str(fc)] >= thr)
            tail = {"L99": thr, "exceed_T": eT, "exceed_C_total": eC,
                    "p_one_sided": binom_sf(eT, eT + eC, 1.0 / 9.0)}
        else:
            tail = {}
        pubs = []
        for p in pubrecs:
            if p["D"] != dmin:
                continue
            q = dict(p)
            if len(ctrl_vals) and math.isfinite(p.get("L_T3.2", float("nan"))):
                q["L_T_percentile_among_controls"] = float((ctrl_vals < p["L_T3.2"]).mean() * 100)
            pubs.append(q)
        dres["rank_test"] = rk
        dres["tail_test"] = tail
        dres["published_windows"] = pubs
        res[f"D{dmin}"] = dres
    res["clusters_5.8"] = clusters
    res["encounters"] = enc_info
    return res


def grade(res):
    g = {}
    info = res["encounters"]
    a_ok = sum(1 for e in info if e["enc"] <= 24 and e.get("valid_frac_D135", 0) >= 0.70)
    b_ok = sum(1 for e in info if e["enc"] >= 25 and e.get("valid_frac_D135", 0) >= 0.50)
    num_ok = bool(res.get("numbering_ok"))
    g["P0"] = "PASS" if (num_ok and a_ok >= 20 and b_ok >= 3) else "FAIL"
    g["P0_detail"] = {"numbering_ok": num_ok, "E1_E24_ge70pct": a_ok, "E25_E27_ge50pct": b_ok}
    r0 = []
    for dmin in DURATIONS_MIN:
        for p in res[f"D{dmin}"]["published_windows"]:
            f_ok = 3.0 <= p.get("fpk_mhz", 0) <= 3.35
            s = p.get("S_fpk", float("nan"))
            s_ok = math.isfinite(s) and abs(s - p["published_sigma"]) <= 1.5
            r0.append(f_ok and s_ok)
    g["R0"] = "PASS" if (len(r0) == 2 and all(r0)) else "FAIL"
    kf = {d: res[f"D{d}"]["A_E1_E24"]["kF_flanks"] for d in DURATIONS_MIN}
    if any(v >= 4 for v in kf.values()):
        g["H1"] = "CHANCE-LEVEL"
    elif all(v <= 1 for v in kf.values()):
        g["H1"] = "RARE"
    else:
        g["H1"] = "UNDECIDED"
    z = res[f"D{DURATIONS_MIN[0]}"]["rank_test"]["all"]["z"]
    if not math.isfinite(z):
        g["H2"] = "NOT GRADABLE"
    elif z >= 3.0:
        g["H2"] = "SUPPORT"
    elif z < 2.0:
        g["H2"] = "NO-EXCESS"
    else:
        g["H2"] = "INDETERMINATE"
    if g["P0"] == "FAIL":
        g["note"] = "P0 FAIL: H1 and H2 are NOT GRADABLE (reported indicative only)."
    return g


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--budget-s", type=float, default=15000)
    ap.add_argument("--keep-raw", action="store_true")
    ap.add_argument("--only", default="", help="comma list of encounter numbers")
    a = ap.parse_args(argv)
    t_start = time.time()
    os.makedirs(a.out_dir, exist_ok=True)
    wdir = os.path.join(a.out_dir, "windows")
    os.makedirs(wdir, exist_ok=True)

    def log(msg):
        print(msg, flush=True)

    try:
        import cdflib  # noqa: PLC0415
        cv = cdflib.__version__
    except Exception as e:  # noqa: BLE001
        cv = f"missing ({e!r})"
    log(f"P5MIN_VERSION: {VERSION} numpy {np.__version__} cdflib {cv} "
        f"sha256 {sha256_file(os.path.abspath(__file__))}")
    lays = {d: layout(d * 60) for d in DURATIONS_MIN}
    log("P5MIN_LAYOUT: " + json.dumps({str(d): lays[d] for d in lays}))

    pos_t, pos_r = fetch_positions()
    peris = find_perihelia(pos_t, pos_r)
    mag_stop = to_unix(MAG_STOP)
    allp = [(i + 1, t, r) for i, (t, r) in enumerate(peris)]
    with open(os.path.join(a.out_dir, "perihelia.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["enc", "perihelion_utc", "r_au", "r_rs"])
        for e, t, r in allp:
            w.writerow([e, iso(t), f"{r:.5f}", f"{r * RS_PER_AU:.3f}"])
    num_ok = True
    for e, exp in EXPECTED_PERI.items():
        hit = [t for (n, t, r) in allp if n == e]
        if not hit or abs(hit[0] - to_unix(exp)) > 2 * 86400:
            num_ok = False
    log(f"P5MIN_PERIHELIA: n={len(allp)} numbering_ok={num_ok} "
        + " ".join(f"E{e}={iso(t)[:13]}" for e, t, r in allp))
    encs = [(e, t) for e, t, r in allp if t + HALF_SPAN_DAYS * 86400 <= mag_stop]
    if a.only:
        want = {int(x) for x in a.only.split(",") if x.strip()}
        encs = [x for x in encs if x[0] in want]
    encs.sort(key=lambda x: (0, x[0]) if x[0] >= 22 else (1, -x[0]))

    rows, pubrecs, done, enc_info = [], [], [], []
    for enc, pt in encs:
        if time.time() - t_start > a.budget_s:
            log(f"P5MIN_BUDGET: stopping before E{enc}")
            break
        t0 = time.time()
        g0, B, meta = encounter_bins(enc, pt, a.data_dir, a.keep_raw, log)
        cover = float(np.isfinite(B[:, 0]).mean())
        er = run_encounter(enc, pt, B, g0, pos_t, pos_r, lays)
        pubrecs += published_windows(enc, B, g0, lays)
        write_rows(os.path.join(wdir, f"E{enc:02d}.csv"), er)
        rows += er
        done.append(enc)
        vf = {}
        for d in DURATIONS_MIN:
            sl = [r for r in er if r["D"] == d and r["sliding"] == 1]
            vf[f"valid_frac_D{d}"] = round(sum(r["valid"] for r in sl) / len(sl), 4) if sl else 0.0
        info = {"enc": enc, "perihelion": iso(pt), "coverage_0p5s": round(cover, 4),
                **vf, **meta, "seconds": round(time.time() - t0, 1)}
        enc_info.append(info)
        log("P5MIN_ENC: " + json.dumps(info))
        del B

    res = summarise(rows, pubrecs, done, lays, enc_info)
    res["numbering_ok"] = num_ok
    res["grades"] = grade(res)
    res["runtime_s"] = round(time.time() - t_start, 1)
    with open(os.path.join(a.out_dir, "results.json"), "w") as f:
        json.dump(res, f, indent=1, default=jdefault)
    write_rows(os.path.join(a.out_dir, "clusters_5p8.csv"), res["clusters_5.8"])
    write_rows(os.path.join(a.out_dir, "published_windows.csv"), pubrecs)
    for d in DURATIONS_MIN:
        dr = res[f"D{d}"]
        for samp in ("A_E1_E24", "B_E25_E27"):
            log(f"P5MIN_D{d}_{samp}: " + json.dumps(dr[samp], default=jdefault))
        log(f"P5MIN_D{d}_RANK: " + json.dumps(dr["rank_test"], default=jdefault))
        log(f"P5MIN_D{d}_TAIL: " + json.dumps(dr["tail_test"], default=jdefault))
        log(f"P5MIN_D{d}_PUBLISHED: " + json.dumps(dr["published_windows"], default=jdefault))
    for c in res["clusters_5.8"]:
        log("P5MIN_CLUSTER: " + json.dumps(c, default=jdefault))
    log("P5MIN_GRADES: " + json.dumps(res["grades"], default=jdefault))
    return 0


if __name__ == "__main__":
    sys.exit(main())
