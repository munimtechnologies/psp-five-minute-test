#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. Written by AI agents working for MV2 of Munim, Inc.
"""patchclock 0.1.0 - do Parker Solar Probe's switchback patches have a size
on the Sun, a lifetime, or both?

Two published readings disagree:
  * spatial imprint: patches are fixed structures the size of supergranules
    (5-10 deg of Carrington longitude; Bale et al., Nature 618, 252, 2023;
    Fargette et al. 2021, ApJ, doi:10.3847/1538-4357/ac1112); their duration
    shrinks from ~10 h to ~2 h toward perihelion as the spacecraft sweeps faster;
  * temporal: near corotation, patches still last 4-6 h although the
    spacecraft barely moves over the Sun, so they are bursts with a
    lifetime (Shi et al. 2022, ApJ, doi:10.3847/1538-4357/ac7c11; arXiv:2206.03807).

This script measures both at once. Switchback activity (fraction of 0.5 s
magnetic-field samples deflected more than 60 deg from the Parker-spiral
direction) is formed in 10-minute blocks for encounters E1-E27, and each
pair of blocks up to 12 h apart is placed by its time lag dt and by the
distance |dphi| between the ballistic Carrington footpoints of the two
blocks. Parker's orbit separates the two: at the corotation turn-arounds
dphi stays near zero for many hours, while at the 9.86 Rsun perihelia of
E22-E27 the footpoint sweeps ~120 deg per day. The space-time variogram
gamma(dt, dphi) = c0 + c1 (1 - exp(-dt/tau - |dphi|/theta)) is fitted on a
grid; tau is the patch lifetime and theta the patch angular size.
Uncertainties: encounter-block bootstrap. Calibration: synthetic patches
with known tau and theta, sampled on the real footpoint tracks, must be
recovered (gate G1). Rules: patchclock/PREREG-PATCH-v0.1.md.

Reuses, unchanged and pinned: psp5min/p5min.py 0.1.0 (positions,
perihelion numbering, FIELDS magnetometer binning to 0.5 s) and
alfven/psp5alf.py 0.1.0 (SWEAP SPAN-i moments via HAPI).
Data (public, NASA CDAWeb): PSP_FLD_L2_MAG_RTN_4_SA_PER_CYC,
PSP_SWP_SPI_SF00_L3_MOM (VEL_RTN_SUN), PSP_HELIO1HR_POSITION (RAD_AU,
HGI_LON; HG_LON only as a frame check).
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
import urllib.parse

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in (os.path.join(os.path.dirname(_HERE), "psp5min"),
           os.path.join(os.path.dirname(_HERE), "alfven"), _HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import p5min as P  # noqa: E402

VERSION = "0.1.0"
DT = 0.5                         # s, p5min grid
BLOCK_S = 600                    # 10-minute activity blocks
NB = int(round(BLOCK_S / DT))    # 0.5 s bins per block
HOUR_BLOCKS = 6                  # polarity is decided per hour
NLAG = 72                        # pair lags 1..72 blocks = 10 min .. 12 h
RS_KM = 695700.0
OMEGA_DEG_DAY = 14.1844          # Carrington (sidereal) rotation rate
OMEGA_DEG_S = OMEGA_DEG_DAY / 86400.0
OMEGA_RAD_S = math.radians(OMEGA_DEG_DAY) / 86400.0
V_REF = 350.0                    # km/s, nominal speed for the Parker reference direction
COS60 = 0.5                      # primary deflection threshold (60 deg)
POL_MIN = 0.3                    # |median B.e / median |B|| needed to call an hour's polarity
BLOCK_MIN_FRAC = 0.8             # finite 0.5 s bins needed in a block
R_MAX_RS = 50.0
R_INNER_RS = 25.0                # robustness fit
MIN_ENC_BLOCKS = 200             # usable blocks for an encounter to enter the fit
V_MIN, V_MAX = 50.0, 2000.0
V_WIN_S = 6 * 3600               # ballistic speed: 6-hour median of SPAN-i V_R
V_MIN_N = 30
MAX_POS_GAP_S = 3 * 3600
POS_DATASET = "PSP_HELIO1HR_POSITION"
POS_START = "2018-08-13T00:00:00Z"
POS_STOP = "2026-06-01T00:00:00Z"
DPHI_EDGES = np.array([0, 0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12.5, 15,
                       20, 25, 30, 40, 50, 60, 80, 100, 150, 1e9], float)
NPHI = len(DPHI_EDGES) - 1
TAU_GRID = np.logspace(math.log10(0.25), math.log10(2000.0), 71)      # h
THETA_GRID = np.logspace(math.log10(0.1), math.log10(2000.0), 71)     # deg
CELL_MIN_N = 30
OMEGA_EDGES = (0.0, 5.0, 10.0, 20.0, 40.0, 80.0, 1e9)                 # |dlon_sc/dt|, deg/day
NCLS = len(OMEGA_EDGES) - 1
N_BOOT = 1000
SEED = 20261124
NEW_FIRST_ENC = 22
# injection-recovery scenarios: (name, tau_h, theta_deg, footpoint-error sigma deg);
# inf = absent. S5 (information only) samples a fixed spatial pattern at a
# footpoint that wanders from the ballistic one (OU, sigma 2 deg, 6 h).
SCENARIOS = (("S1_spatial", math.inf, 5.0, 0.0), ("S2_temporal", 5.0, math.inf, 0.0),
             ("S3_mixed", 5.0, 5.0, 0.0), ("S4_long", 36.0, 7.5, 0.0),
             ("S5_spatial_mapnoise", math.inf, 5.0, 2.0))
MAPNOISE_TAU_H = 6.0
INJ_NUGGET = 0.5
INJ_SEEDS = (11, 12)


# ------------------------------------------------------------- positions
def parse_hapi_cols(text: str, ncol: int):
    ts, rows = [], []
    for line in text.splitlines():
        line = line.strip()
        if not line or line[0] in "#{}\"":
            continue
        parts = line.split(",")
        if len(parts) < ncol + 1:
            continue
        try:
            t = P.to_unix(parts[0])
            v = [float(x) for x in parts[1:ncol + 1]]
        except ValueError:
            continue
        if any((not math.isfinite(x)) or x < -1e30 for x in v):
            continue
        ts.append(t)
        rows.append(v)
    return np.asarray(ts, float), np.asarray(rows, float).reshape(-1, ncol)


def fetch_hg(data_dir: str, log):
    """Hourly RAD_AU, HG_LON (rotating) and HGI_LON (inertial) for the mission."""
    cache = os.path.join(data_dir, "patchclock_pos_hgi.npz")
    if os.path.exists(cache):
        z = np.load(cache)
        return z["t"], z["r_au"], z["hg"], z["hgi"]
    a, b = P.to_unix(POS_START), P.to_unix(POS_STOP)
    out_t, out_v = [], []
    cur = a
    while cur < b:
        nxt = min(b, cur + 366 * 86400)
        q = {"id": POS_DATASET, "parameters": "RAD_AU,HG_LON,HGI_LON",
             "time.min": P.iso(cur), "time.max": P.iso(nxt), "format": "csv"}
        try:
            txt = P.http_get(P.HAPI + "?" + urllib.parse.urlencode(q)).decode()
        except Exception:  # noqa: BLE001  (HAPI 3 parameter names)
            q2 = {"id": POS_DATASET, "parameters": "RAD_AU,HG_LON,HGI_LON",
                  "start": P.iso(cur), "stop": P.iso(nxt), "format": "csv"}
            txt = P.http_get(P.HAPI + "?" + urllib.parse.urlencode(q2)).decode()
        t, v = parse_hapi_cols(txt, 3)
        log(f"PCLK_POS: {P.iso(cur)[:10]}..{P.iso(nxt)[:10]} rows={len(t)}")
        out_t.append(t)
        out_v.append(v)
        cur = nxt
    t = np.concatenate(out_t)
    v = np.concatenate(out_v)
    o = np.argsort(t, kind="stable")
    t, v = t[o], v[o]
    keep = np.r_[True, np.diff(t) > 0]
    t, v = t[keep], v[keep]
    np.savez(cache, t=t, r_au=v[:, 0], hg=v[:, 1], hgi=v[:, 2])
    return t, v[:, 0], v[:, 1], v[:, 2]


def _unwrap_rate(lon, t, seg):
    lon_u = np.empty_like(lon)
    rate = np.full(len(t), np.nan)
    for s in np.unique(seg):
        m = seg == s
        lon_u[m] = np.degrees(np.unwrap(np.radians(lon[m])))
        if m.sum() >= 3:
            rate[m] = np.gradient(lon_u[m], t[m]) * 86400.0
    return lon_u, rate


def prep_positions(t, r_au, hgi_lon, hg_lon=None):
    """Longitude in a frame rotating at the Carrington rate, built from the
    inertial HGI_LON (lon_u = unwrapped HGI_LON - Omega t; the constant offset
    from true Carrington longitude cancels in every difference), and its rate
    in deg/day. If HG_LON is given, the median of d(HG)/dt - d(HGI)/dt is
    returned as a frame check (Carrington: -14.18 deg/day)."""
    t = np.asarray(t, float)
    r_au = np.asarray(r_au, float)
    hgi_lon = np.asarray(hgi_lon, float)
    seg = np.concatenate(([0], np.cumsum(np.diff(t) > MAX_POS_GAP_S))).astype(np.int64)
    hgi_u, hgi_rate = _unwrap_rate(hgi_lon, t, seg)
    t_ref = t[0] if len(t) else 0.0
    lon_u = hgi_u - OMEGA_DEG_S * (t - t_ref)
    omega = hgi_rate - OMEGA_DEG_DAY
    frame_check = float("nan")
    if hg_lon is not None and len(t) >= 3:
        _, hg_rate = _unwrap_rate(np.asarray(hg_lon, float), t, seg)
        d = hg_rate - hgi_rate
        d = d[np.isfinite(d)]
        if len(d):
            frame_check = float(np.median(d))
    return {"t": t, "r_rs": r_au * P.RS_PER_AU, "lon_u": lon_u, "omega": omega, "seg": seg,
            "frame_check_deg_day": frame_check}


def interp_gapped(t, x, tq, max_gap_s=MAX_POS_GAP_S):
    t = np.asarray(t, float)
    x = np.asarray(x, float)
    tq = np.asarray(tq, float)
    y = np.interp(tq, t, x)
    j = np.clip(np.searchsorted(t, tq), 1, len(t) - 1)
    gap = t[j] - t[j - 1]
    bad = (tq < t[0]) | (tq > t[-1]) | (gap > max_gap_s)
    y[bad] = np.nan
    return y


def seg_of(pos, tq):
    j = np.clip(np.searchsorted(pos["t"], tq, side="right") - 1, 0, len(pos["t"]) - 1)
    return pos["seg"][j]


# ------------------------------------------------------------- activity
def hour_polarity(proj, bmag, ok, nblk):
    """Per block: sign (+1/-1) of its hour's median B.e_out, 0 if ambiguous."""
    pol = np.zeros(nblk)
    pfrac = np.full(nblk, np.nan)
    for h0 in range(0, nblk, HOUR_BLOCKS):
        sl = slice(h0, min(nblk, h0 + HOUR_BLOCKS))
        m = ok[sl]
        if m.sum() < 0.5 * m.size:
            continue
        mb = float(np.median(bmag[sl][m]))
        if not mb > 0:
            continue
        p = float(np.median(proj[sl][m])) / mb
        pfrac[sl] = p
        if abs(p) >= POL_MIN:
            pol[sl] = 1.0 if p > 0 else -1.0
    return pol, pfrac


def block_activity(B, r_rs):
    """Fractions of 0.5 s samples deflected > 60 deg (f60) and > 90 deg (f90)
    from the polarity-signed Parker-spiral direction, per 10-minute block."""
    nblk = min(B.shape[0] // NB, len(r_rs))
    X = np.asarray(B[: nblk * NB], np.float64).reshape(nblk, NB, 3)
    r = np.asarray(r_rs[:nblk], float)
    tanpsi = OMEGA_RAD_S * r * RS_KM / V_REF
    norm = np.sqrt(1.0 + tanpsi ** 2)
    with np.errstate(invalid="ignore", divide="ignore"):
        bmag = np.sqrt((X ** 2).sum(axis=2))
        proj = (X[:, :, 0] - tanpsi[:, None] * X[:, :, 1]) / norm[:, None]
        c = proj / bmag
    ok = np.isfinite(c) & (bmag > 0)
    pol, pfrac = hour_polarity(proj, bmag, ok, nblk)
    cosd = np.where(ok, c * pol[:, None], 1.0)
    nok = ok.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        f60 = (ok & (cosd < COS60)).sum(axis=1) / nok
        f90 = (ok & (cosd < 0.0)).sum(axis=1) / nok
    valid = (nok >= BLOCK_MIN_FRAC * NB) & (pol != 0) & np.isfinite(r) & (r <= R_MAX_RS)
    f60 = np.where(valid, f60, np.nan)
    f90 = np.where(valid, f90, np.nan)
    return {"f60": f60, "f90": f90, "valid": valid, "pol": pol, "pfrac": pfrac,
            "nok": nok, "nblk": nblk}


def window_median_v(spi_t, spi_vr, tc, win_s=V_WIN_S):
    t = np.asarray(spi_t, float)
    v = np.asarray(spi_vr, float)
    m = np.isfinite(t) & np.isfinite(v) & (v > V_MIN) & (v < V_MAX)
    t, v = t[m], v[m]
    o = np.argsort(t, kind="stable")
    t, v = t[o], v[o]
    tc = np.asarray(tc, float)
    out = np.full(len(tc), np.nan)
    lo = np.searchsorted(t, tc - win_s / 2.0)
    hi = np.searchsorted(t, tc + win_s / 2.0)
    for i in range(len(tc)):
        if hi[i] - lo[i] >= V_MIN_N:
            out[i] = float(np.median(v[lo[i]:hi[i]]))
    return out


def footpoint(lon_sc, r_rs, v_kms):
    """Ballistic Carrington longitude of the source: phi_sc + Omega (r - Rsun) / V."""
    return lon_sc + OMEGA_DEG_S * (np.asarray(r_rs, float) - 1.0) * RS_KM / np.asarray(v_kms, float)


def encounter_blocks(enc, g0, B, pos, spi_t, spi_vr):
    nblk = B.shape[0] // NB
    tc = g0 + (np.arange(nblk) + 0.5) * BLOCK_S
    r_rs = interp_gapped(pos["t"], pos["r_rs"], tc)
    lon = interp_gapped(pos["t"], pos["lon_u"], tc)
    om = interp_gapped(pos["t"], pos["omega"], tc)
    act = block_activity(B, r_rs)
    v6 = window_median_v(spi_t, spi_vr, tc)
    phi_fp = footpoint(lon, r_rs, v6)
    ok_sc = act["valid"] & np.isfinite(lon) & np.isfinite(act["f60"])
    ok_fp = ok_sc & np.isfinite(phi_fp)
    return {"enc": int(enc), "g0": float(g0), "tc": tc, "r_rs": r_rs, "lon": lon, "omega": om,
            "v6": v6, "phi_fp": phi_fp, "seg": seg_of(pos, tc), "f60": act["f60"],
            "f90": act["f90"], "pfrac": act["pfrac"], "ok_sc": ok_sc, "ok_fp": ok_fp}


# ------------------------------------------------------------- variogram
def zscore(a, ok, min_n=MIN_ENC_BLOCKS):
    a = np.asarray(a, float)
    z = np.full(len(a), np.nan)
    ok = np.asarray(ok, bool) & np.isfinite(a)
    if ok.sum() < min_n:
        return z, np.zeros(len(a), bool)
    v = a[ok]
    sd = float(v.std())
    if not sd > 0:
        return z, np.zeros(len(a), bool)
    z[ok] = (v - v.mean()) / sd
    return z, ok


def cell_sums(z, phi, seg, ok):
    """Pair sums by (lag, |dphi| bin): count, sum of 0.5 (z_i - z_j)^2,
    sum of dt (h), sum of |dphi| (deg). Shape (4, NLAG, NPHI)."""
    out = np.zeros((4, NLAG, NPHI))
    n = len(z)
    ok = np.asarray(ok, bool) & np.isfinite(z) & np.isfinite(phi)
    for k in range(1, NLAG + 1):
        if k >= n:
            break
        a, b = slice(0, n - k), slice(k, n)
        m = ok[a] & ok[b] & (seg[a] == seg[b])
        if not m.any():
            continue
        g = 0.5 * (z[a][m] - z[b][m]) ** 2
        dp = np.abs(phi[b][m] - phi[a][m])
        j = np.clip(np.searchsorted(DPHI_EDGES, dp, side="right") - 1, 0, NPHI - 1)
        cnt = np.bincount(j, minlength=NPHI).astype(float)
        out[0, k - 1] += cnt
        out[1, k - 1] += np.bincount(j, weights=g, minlength=NPHI)
        out[2, k - 1] += cnt * (k * BLOCK_S / 3600.0)
        out[3, k - 1] += np.bincount(j, weights=dp, minlength=NPHI)
    return out


def make_design(S):
    """Fixed cells (n >= CELL_MIN_N in the pooled sums S) and the model matrix
    X[g, cell] = 1 - exp(-dt/tau_g - |dphi|/theta_g) over the (tau, theta) grid."""
    n = S[0].ravel()
    keep = n >= CELL_MIN_N
    dt = S[2].ravel()[keep] / n[keep]
    dp = S[3].ravel()[keep] / n[keep]
    tg, thg = np.meshgrid(TAU_GRID, THETA_GRID, indexing="ij")
    tg, thg = tg.ravel(), thg.ravel()
    X = 1.0 - np.exp(-(dt[None, :] / tg[:, None] + dp[None, :] / thg[:, None]))
    return {"keep": keep, "X": X, "X2": X * X, "tau": tg, "theta": thg,
            "dt": dt, "dp": dp, "ncell": int(keep.sum())}


def wls_grid(D, n, gsum):
    """Best (tau, theta) on the grid by weighted least squares of the cell-mean
    semivariance on c0 + c1 X; c1 must be positive."""
    w = np.asarray(n, float)
    W = float(w.sum())
    if W <= 0 or D["ncell"] < 10:
        return None
    y = np.where(w > 0, np.asarray(gsum, float) / np.maximum(w, 1e-300), 0.0)
    wy = w * y
    swy = float(wy.sum())
    Sx = D["X"] @ w
    Sxx = D["X2"] @ w - Sx ** 2 / W
    Sxy = D["X"] @ wy - Sx * swy / W
    with np.errstate(invalid="ignore", divide="ignore"):
        c1 = Sxy / Sxx
        score = np.where((Sxx > 1e-12) & (c1 > 0), Sxy ** 2 / np.maximum(Sxx, 1e-300), -np.inf)
    i = int(np.argmax(score))
    if not np.isfinite(score[i]):
        return None
    c0 = (swy - c1[i] * Sx[i]) / W
    syy = float((wy * y).sum() - swy ** 2 / W)
    nth = len(THETA_GRID)
    it, ith = i // nth, i % nth
    return {"tau_h": float(D["tau"][i]), "theta_deg": float(D["theta"][i]),
            "c0": float(c0), "c1": float(c1[i]), "r2": float(score[i] / syy) if syy > 0 else float("nan"),
            "tau_at_grid_max": int(it == len(TAU_GRID) - 1),
            "theta_at_grid_max": int(ith == nth - 1),
            "tau_at_grid_min": int(it == 0), "theta_at_grid_min": int(ith == 0),
            "pairs": W, "cells": D["ncell"]}


def fit_cells(S):
    D = make_design(S)
    f = wls_grid(D, S[0].ravel()[D["keep"]], S[1].ravel()[D["keep"]])
    return f, D


def bootstrap(cells_list, D, n_boot=N_BOOT, seed=SEED):
    """Encounter-block bootstrap of (tau, theta) with the cells held fixed."""
    keep = D["keep"]
    Ns = np.array([c[0].ravel()[keep] for c in cells_list])
    Gs = np.array([c[1].ravel()[keep] for c in cells_list])
    nE = len(cells_list)
    rng = np.random.default_rng(seed)
    taus, ths = [], []
    for _ in range(n_boot):
        cnt = np.bincount(rng.integers(0, nE, nE), minlength=nE).astype(float)
        f = wls_grid(D, cnt @ Ns, cnt @ Gs)
        if f:
            taus.append(f["tau_h"])
            ths.append(f["theta_deg"])
    if len(taus) < 0.5 * n_boot:
        return {"n_ok": len(taus)}
    lt, lh = np.log(taus), np.log(ths)
    return {"n_ok": len(taus),
            "tau95_h": [float(np.exp(np.percentile(lt, 2.5))), float(np.exp(np.percentile(lt, 97.5)))],
            "theta95_deg": [float(np.exp(np.percentile(lh, 2.5))), float(np.exp(np.percentile(lh, 97.5)))],
            "tau_frac_at_grid_max": float(np.mean(np.asarray(taus) >= TAU_GRID[-1] * 0.999)),
            "theta_frac_at_grid_max": float(np.mean(np.asarray(ths) >= THETA_GRID[-1] * 0.999))}


# ------------------------------------------------------------- model-free curves
def modelfree_sums(z, om, seg, ok):
    """Semivariance by lag for pairs whose two blocks share a |dlon_sc/dt| class.
    Shape (3, NCLS, NLAG): count, sum gamma, sum |omega|."""
    out = np.zeros((3, NCLS, NLAG))
    om = np.asarray(om, float)
    n = len(z)
    okk = np.asarray(ok, bool) & np.isfinite(z) & np.isfinite(om)
    cls = np.searchsorted(np.asarray(OMEGA_EDGES), np.abs(np.where(np.isfinite(om), om, 0.0)),
                          side="right") - 1
    for k in range(1, NLAG + 1):
        if k >= n:
            break
        a, b = slice(0, n - k), slice(k, n)
        m = okk[a] & okk[b] & (seg[a] == seg[b]) & (cls[a] == cls[b])
        if not m.any():
            continue
        c = cls[a][m]
        g = 0.5 * (z[a][m] - z[b][m]) ** 2
        out[0, :, k - 1] += np.bincount(c, minlength=NCLS)[:NCLS]
        out[1, :, k - 1] += np.bincount(c, weights=g, minlength=NCLS)[:NCLS]
        out[2, :, k - 1] += np.bincount(c, weights=np.abs(om[a][m]), minlength=NCLS)[:NCLS]
    return out


def efold_time(gsum, n, min_n=50):
    """Lag (h) where the semivariance has risen 1 - 1/e of the way from its
    10-minute value to its 8-12 h plateau; nan if it never does."""
    g = np.where(np.asarray(n) >= min_n, np.asarray(gsum) / np.maximum(n, 1), np.nan)
    if not np.isfinite(g[0]):
        return float("nan")
    tail = g[47:NLAG]
    if not np.isfinite(tail).any():
        return float("nan")
    plat = float(np.nanmean(tail))
    if not plat > g[0]:
        return float("nan")
    target = g[0] + (1.0 - math.exp(-1.0)) * (plat - g[0])
    prev_k, prev_g = 0, g[0]
    for k in range(1, NLAG):
        if not np.isfinite(g[k]):
            continue
        if g[k] >= target:
            lag0 = (prev_k + 1) * BLOCK_S / 3600.0
            lag1 = (k + 1) * BLOCK_S / 3600.0
            fr = (target - prev_g) / (g[k] - prev_g) if g[k] > prev_g else 1.0
            return float(lag0 + fr * (lag1 - lag0))
        prev_k, prev_g = k, g[k]
    return float("nan")


# ------------------------------------------------------------- injection
def sim_activity(t_h, phi, tau_h, theta_deg, nugget, rng):
    """Gaussian field with covariance exp(-|dt|/tau - |dphi|/theta), sampled
    at (t_h, phi), plus white noise of variance `nugget`."""
    t_h = np.asarray(t_h, float)
    phi = np.asarray(phi, float)
    m = len(t_h)
    r = np.zeros((m, m))
    if math.isfinite(tau_h):
        r += np.abs(t_h[:, None] - t_h[None, :]) / tau_h
    if math.isfinite(theta_deg):
        r += np.abs(phi[:, None] - phi[None, :]) / theta_deg
    C = np.exp(-r)
    del r
    C[np.diag_indices(m)] += 1e-6
    try:
        L = np.linalg.cholesky(C)
        x = L @ rng.standard_normal(m)
    except np.linalg.LinAlgError:
        w, V = np.linalg.eigh(C)
        x = V @ (np.sqrt(np.clip(w, 0.0, None)) * rng.standard_normal(m))
    return x + math.sqrt(nugget) * rng.standard_normal(m)


def within(x, ref, fac=1.5):
    return bool(math.isfinite(x) and ref / fac <= x <= ref * fac)


def scenario_pass(name, f):
    if not f:
        return False
    tau, th = f["tau_h"], f["theta_deg"]
    if name == "S1_spatial":
        return within(th, 5.0) and tau >= 48.0
    if name == "S2_temporal":
        return within(tau, 5.0) and th >= 40.0
    if name == "S3_mixed":
        return within(tau, 5.0) and within(th, 5.0)
    if name == "S4_long":
        return within(th, 7.5) and tau >= 18.0
    if name == "S5_spatial_mapnoise":      # information only, never a gate
        return within(th, 5.0) and tau >= 24.0
    return False


def ou_noise(t_h, sigma, tau_h, rng):
    """Exact Ornstein-Uhlenbeck sample (stationary sd sigma) at increasing times t_h."""
    t_h = np.asarray(t_h, float)
    x = np.zeros(len(t_h))
    if len(t_h) == 0 or sigma <= 0:
        return x
    x[0] = sigma * rng.standard_normal()
    for i in range(1, len(t_h)):
        rho = math.exp(-max(t_h[i] - t_h[i - 1], 0.0) / tau_h)
        x[i] = rho * x[i - 1] + sigma * math.sqrt(1.0 - rho * rho) * rng.standard_normal()
    return x


def run_injection(E, D, log):
    out = {}
    for name, tau, th, mapsig in SCENARIOS:
        out[name] = []
        for sd in INJ_SEEDS:
            rng = np.random.default_rng(SEED + 1000 * sd + len(name))
            tot = np.zeros((4, NLAG, NPHI))
            for e in E:
                ok = e["used_fp"]
                a = np.full(len(ok), np.nan)
                idx = np.flatnonzero(ok)
                if len(idx) == 0:
                    continue
                th_ = e["tc"][idx] / 3600.0
                ph_ = e["phi_fp"][idx] + ou_noise(th_, mapsig, MAPNOISE_TAU_H, rng)
                a[idx] = sim_activity(th_, ph_, tau, th, INJ_NUGGET, rng)
                z, okz = zscore(a, ok)
                tot += cell_sums(z, e["phi_fp"], e["seg"], okz)
            f = wls_grid(D, tot[0].ravel()[D["keep"]], tot[1].ravel()[D["keep"]])
            rec = {"seed": sd, "fit": f, "pass": scenario_pass(name, f)}
            out[name].append(rec)
            log(f"PCLK_INJ: {name} tau_in={tau} theta_in={th} seed={sd} "
                f"tau_fit={f['tau_h'] if f else 'nan'} theta_fit={f['theta_deg'] if f else 'nan'} "
                f"pass={rec['pass']}")
    return out


# ------------------------------------------------------------- grades
def grade_g0(enc_info):
    good = [e for e in enc_info if e.get("used_fp", 0) >= MIN_ENC_BLOCKS]
    n_all = len(good)
    n_new = sum(1 for e in good if e["enc"] >= NEW_FIRST_ENC)
    return ("PASS" if (n_all >= 15 and n_new >= 4) else "FAIL"), n_all, n_new


def grade_g1(inj):
    core = ("S1_spatial", "S2_temporal", "S3_mixed")
    ok = all(inj.get(k) and all(r["pass"] for r in inj[k]) for k in core)
    long_ok = bool(inj.get("S4_long")) and all(r["pass"] for r in inj["S4_long"])
    return ("PASS" if ok else "FAIL"), long_ok


def grade_h1(theta, ci):
    if not (math.isfinite(theta) and ci and len(ci) == 2):
        return "NOT GRADABLE"
    lo, hi = ci
    if 2.0 <= theta <= 15.0 and hi <= 30.0:
        return "SUPERGRANULAR"
    if theta < 2.0 and hi < 4.0:
        return "SMALLER"
    if lo > 15.0:
        return "NONE"
    return "UNRESOLVED"


def grade_h2(tau, ci, long_ok):
    if not (math.isfinite(tau) and ci and len(ci) == 2):
        return "NOT GRADABLE"
    lo, hi = ci
    if tau <= 8.0 and hi <= 16.0:
        return "SHORT"
    if lo >= 24.0:
        return "LONG" if long_ok else "LONG (UNCALIBRATED)"
    if 8.0 < tau < 24.0 and lo >= 4.0 and hi <= 48.0:
        return "MIDDLE"
    return "UNRESOLVED"


def verdict(g0, g1, h1, h2):
    if g0 != "PASS" or g1 != "PASS":
        return "NOT GRADABLE"
    if h1 == "SUPERGRANULAR" and h2 == "LONG":
        return "FIXED SUPERGRANULAR IMPRINT"
    if h1 == "NONE" and h2 == "SHORT":
        return "TRANSIENT BURSTS"
    if h1 in ("SUPERGRANULAR", "SMALLER") and h2 in ("SHORT", "MIDDLE"):
        return "BOTH: SIZED AND SHORT-LIVED"
    return "UNDETERMINED"


# ------------------------------------------------------------- outputs
def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def jd(v):
    return P.jdefault(v)


def write_blocks(path, E):
    cols = ["enc", "t_mid", "r_rs", "lon_sc_deg", "dlon_dt_deg_day", "v6_kms", "phi_fp_deg",
            "seg", "f60", "f90", "polarity_frac", "used_fp", "used_sc"]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for e in E:
            for i in np.flatnonzero(e["used_sc"]):
                w.writerow([e["enc"], P.iso(e["tc"][i]), P.fmt(e["r_rs"][i]), P.fmt(e["lon"][i]),
                            P.fmt(e["omega"][i]), P.fmt(e["v6"][i]), P.fmt(e["phi_fp"][i]),
                            int(e["seg"][i]), P.fmt(e["f60"][i]), P.fmt(e["f90"][i]),
                            P.fmt(e["pfrac"][i]), int(e["used_fp"][i]), int(e["used_sc"][i])])


def write_cells(path, S):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["lag_h", "dphi_lo_deg", "dphi_hi_deg", "pairs", "gamma_mean",
                    "dt_mean_h", "dphi_mean_deg"])
        for k in range(NLAG):
            for j in range(NPHI):
                n = S[0, k, j]
                if n <= 0:
                    continue
                w.writerow([round((k + 1) * BLOCK_S / 3600.0, 4), DPHI_EDGES[j],
                            DPHI_EDGES[j + 1] if DPHI_EDGES[j + 1] < 1e8 else "inf", int(n),
                            P.fmt(S[1, k, j] / n), P.fmt(S[2, k, j] / n), P.fmt(S[3, k, j] / n)])


def sum_cells(E, key, sel=None):
    tot = np.zeros((4, NLAG, NPHI))
    for e in E:
        if sel is None or sel(e):
            tot += e[key]
    return tot


def analyse(E, log, n_boot=N_BOOT, do_injection=True):
    """Fits, bootstrap, robustness, model-free curves, injection, grades."""
    for e in E:
        z, used = zscore(e["f60"], e["ok_fp"])
        e["z60"], e["used_fp"] = z, used
        e["cells_fp"] = cell_sums(z, e["phi_fp"], e["seg"], used)
        z90, u90 = zscore(e["f90"], e["ok_fp"])
        e["cells_fp90"] = cell_sums(z90, e["phi_fp"], e["seg"], u90)
        zs, us = zscore(e["f60"], e["ok_sc"])
        e["used_sc"] = us
        e["cells_sc"] = cell_sums(zs, e["lon"], e["seg"], us)
        inner = used & (e["r_rs"] <= R_INNER_RS)
        e["cells_inner"] = cell_sums(z, e["phi_fp"], e["seg"], inner)
        e["mf"] = modelfree_sums(z, e["omega"], e["seg"], used)
    res = {}
    S = sum_cells(E, "cells_fp")
    f, D = fit_cells(S)
    res["primary"] = f
    res["bootstrap"] = bootstrap([e["cells_fp"] for e in E if e["used_fp"].sum() > 0], D, n_boot) \
        if f else {}
    log("PCLK_FIT: " + json.dumps({"primary": f, "bootstrap": res["bootstrap"]}, default=jd))
    rob = {}
    for name, key, sel in (("f90", "cells_fp90", None), ("spacecraft_longitude", "cells_sc", None),
                           ("r_le_25", "cells_inner", None),
                           ("E1_E21", "cells_fp", lambda e: e["enc"] < NEW_FIRST_ENC),
                           ("E22_E27", "cells_fp", lambda e: e["enc"] >= NEW_FIRST_ENC)):
        fr, _ = fit_cells(sum_cells(E, key, sel))
        rob[name] = fr
        log(f"PCLK_ROBUST: {name} " + json.dumps(fr, default=jd))
    res["robustness"] = rob
    # model-free decorrelation times by spacecraft angular rate
    MF = np.zeros((3, NCLS, NLAG))
    for e in E:
        MF += e["mf"]
    mf_rows = []
    for c in range(NCLS):
        n = MF[0, c]
        om_mean = float(MF[2, c].sum() / n.sum()) if n.sum() > 0 else float("nan")
        te = efold_time(MF[1, c], n)
        pred = float("nan")
        if f and math.isfinite(om_mean):
            rate = 1.0 / f["tau_h"] + (om_mean / 24.0) / f["theta_deg"]
            pred = 1.0 / rate
        row = {"omega_lo": OMEGA_EDGES[c], "omega_hi": OMEGA_EDGES[c + 1] if OMEGA_EDGES[c + 1] < 1e8 else "inf",
               "pairs": int(n.sum()), "mean_abs_omega_deg_day": om_mean, "efold_h": te,
               "fit_pred_efold_h": pred,
               "gamma_by_lag": [float(MF[1, c, k] / MF[0, c, k]) if MF[0, c, k] >= 50 else None
                                for k in range(NLAG)]}
        mf_rows.append(row)
        log("PCLK_MODELFREE: " + json.dumps({k: v for k, v in row.items() if k != "gamma_by_lag"},
                                             default=jd))
    res["modelfree"] = mf_rows
    res["_S"] = S
    if do_injection and f:
        res["injection"] = run_injection(E, D, log)
    else:
        res["injection"] = {}
    return res


def main(argv=None):
    import psp5alf as A  # noqa: PLC0415  (network-facing helpers, unchanged)
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

    log(f"PCLK_VERSION: {VERSION} numpy {np.__version__} sha256 {sha256_file(os.path.abspath(__file__))} "
        f"p5min {P.VERSION} {sha256_file(P.__file__)} psp5alf {A.VERSION} {sha256_file(A.__file__)}")
    pos_t, pos_r = P.fetch_positions()
    peris = P.find_perihelia(pos_t, pos_r)
    allp = [(i + 1, t, r) for i, (t, r) in enumerate(peris)]
    num_ok = True
    for e, exp in P.EXPECTED_PERI.items():
        hit = [t for (nn, t, r) in allp if nn == e]
        if not hit or abs(hit[0] - P.to_unix(exp)) > 2 * 86400:
            num_ok = False
    log(f"PCLK_PERIHELIA: n={len(allp)} numbering_ok={num_ok} "
        + " ".join(f"E{e}={P.iso(t)[:13]}" for e, t, r in allp))
    ht, hr, hg, hgi = fetch_hg(a.data_dir, log)
    pos = prep_positions(ht, hr, hgi, hg)
    log(f"PCLK_HG: rows={len(ht)} segments={int(pos['seg'].max()) + 1 if len(ht) else 0} "
        f"frame_check d(HG)/dt-d(HGI)/dt median={pos['frame_check_deg_day']:.4f} deg/day "
        f"(Carrington -{OMEGA_DEG_DAY})")
    mag_stop = P.to_unix(P.MAG_STOP)
    encs = [(e, t) for e, t, r in allp if t + P.HALF_SPAN_DAYS * 86400 <= mag_stop]
    if a.only:
        want = {int(x) for x in a.only.split(",") if x.strip()}
        encs = [x for x in encs if x[0] in want]
    encs.sort(key=lambda x: (0, 0) if x[0] == 22 else
              ((1, -x[0]) if x[0] >= NEW_FIRST_ENC else (2, -x[0])))
    E, enc_info = [], []
    for enc, pt in encs:
        if time.time() - t_start > a.budget_s:
            log(f"PCLK_BUDGET: stopping before E{enc}")
            break
        t0 = time.time()
        g0, B, meta = P.encounter_bins(enc, pt, a.data_dir, False, log)
        st, sv = [], []
        day = math.floor(g0 / 86400) * 86400
        while day < g0 + 2 * P.HALF_SPAN_DAYS * 86400:
            tt, aa = A.fetch_spi_day(day, a.data_dir, log)
            if len(tt):
                st.append(tt)
                sv.append(aa[:, 2])
            day += 86400
        spi_t = np.concatenate(st) if st else np.zeros(0)
        spi_vr = np.concatenate(sv) if sv else np.zeros(0)
        e = encounter_blocks(enc, g0, B, pos, spi_t, spi_vr)
        del B
        E.append(e)
        info = {"enc": enc, "perihelion": P.iso(pt), "blocks": int(len(e["tc"])),
                "ok_sc": int(e["ok_sc"].sum()), "ok_fp": int(e["ok_fp"].sum()),
                "min_r_rs": float(np.nanmin(e["r_rs"])) if np.isfinite(e["r_rs"]).any() else None,
                "max_abs_dlon_dt": float(np.nanmax(np.abs(e["omega"][e["ok_sc"]]))) if e["ok_sc"].any() else None,
                "corotating_blocks": int((e["ok_sc"] & (np.abs(np.nan_to_num(e["omega"], nan=1e9)) < 5.0)).sum()),
                "spi_samples": int(len(spi_t)), **{k: meta[k] for k in meta},
                "seconds": round(time.time() - t0, 1)}
        enc_info.append(info)
        log("PCLK_ENC: " + json.dumps(info, default=jd))

    res = analyse(E, log) if E else {"primary": None, "bootstrap": {}, "injection": {}}
    for info in enc_info:
        e = next(x for x in E if x["enc"] == info["enc"])
        info["used_fp"] = int(e["used_fp"].sum())
        info["used_sc"] = int(e["used_sc"].sum())
    g0g, n_all, n_new = grade_g0(enc_info)
    g1g, long_ok = grade_g1(res.get("injection", {}))
    f = res.get("primary")
    bs = res.get("bootstrap", {})
    h1 = grade_h1(f["theta_deg"], bs.get("theta95_deg")) if f else "NOT GRADABLE"
    h2 = grade_h2(f["tau_h"], bs.get("tau95_h"), long_ok) if f else "NOT GRADABLE"
    if g0g != "PASS" or g1g != "PASS":
        h1n, h2n = "NOT GRADABLE", "NOT GRADABLE"
    else:
        h1n, h2n = h1, h2
    grades = {"G0": g0g, "G1": g1g, "S4_long_recovered": long_ok, "H1": h1n, "H2": h2n,
              "VERDICT": verdict(g0g, g1g, h1, h2),
              "H1_if_ungated": h1, "H2_if_ungated": h2,
              "G0_detail": {"encounters_ge_min": n_all, "new_E22plus_ge_min": n_new,
                            "numbering_ok": num_ok}}
    S = res.pop("_S", None)
    out = {"version": VERSION, "grades": grades, "encounters": enc_info,
           "frame_check_deg_day": pos["frame_check_deg_day"],
           "runtime_s": round(time.time() - t_start, 1), **res}
    with open(os.path.join(a.out_dir, "results.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=jd)
    if E:
        write_blocks(os.path.join(a.out_dir, "blocks.csv"), E)
    if S is not None:
        write_cells(os.path.join(a.out_dir, "cells.csv"), S)
    log("PCLK_GRADES: " + json.dumps(grades, default=jd))
    return 0


if __name__ == "__main__":
    sys.exit(main())
