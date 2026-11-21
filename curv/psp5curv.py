#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. Written by AI agents working for MV2 of Munim, Inc.
"""psp5curv 0.1.0 - calibration add-on to psp5min (PREREG amendment 0.1.1).

Written after a second review of psp5min 0.1.0 and BEFORE any of its results
existed. It adds three things psp5min cannot do:

1. H2c, a frequency-local statistic with the SAME fit design in every band.
   L2 = y_k / exp(a), where a is the intercept of a least-squares quadratic
   fit of log y over the bin offsets d = +-2..+-6 around bin k (y = P*f).
   psp5min's L divides by a median over 10 to 48 flank bins depending on
   frequency, so target and control bands do not share one null
   distribution and its rank test is not exactly calibrated.
2. Null simulations (frequency domain, no 3 mHz signal) that measure the
   rank-test Z which L and L2 give for the observed number of windows:
   SIM-F flat P*f, SIM-K smooth spectral knees spread log-uniformly over
   1-10 mHz, SIM-K0 a sharp knee exactly at 3.2 mHz (worst case, reported).
3. H4, the Gaussian false-alarm probability of each published peak under
   the paper's own statistic: surrogate spectra drawn around a smooth
   (weighted cubic in log f) fit to the window's own spectrum with
   2.9-3.5 mHz excluded from the fit. Also the Samuelson ceiling of the
   statistic: S <= sqrt(n_b - 1) whenever the peak bin lies inside the
   n_b-bin baseline.
Reported only: the spectral-knee frequency of every tiled 135-min window.

It imports psp5min 0.1.0 (p5min.py at commit a025a506) unchanged. Magnetometer
bins are read from psp5min's cache in the shared data volume (read only); an
encounter whose cache is missing is rebuilt in a private temporary folder.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import sys
import tempfile
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
DURATIONS_MIN = (135, 70)
DFIT = (-6, -5, -4, -3, -2, 2, 3, 4, 5, 6)
KMAX = 260                     # stored bins per window (30 mHz at 135 min = bin 243)
EXCL_MHZ = (2.9, 3.5)          # left out of the smooth fit for the published windows
SMOOTH_FMAX_MHZ = 30.0
SMOOTH_KFLOOR = 8              # fit weight 1/max(k, 8): equal per log-f, capped at low k
FPK_MHZ = (2.5, 4.0)
N_SURR = {135: 20000, 70: 40000}
CHUNK = 500
N_REP = 200
SEED = 20261122
SIM_N_FALLBACK = 3000
MIN_N_H2C = 2000
SIM_GATE = 0.5
KNEE_DELTA = 0.15
KNEE_GRID = np.linspace(math.log10(0.5e-3), math.log10(20e-3), 41)
KNEE_FIT_MHZ = (0.4, 30.0)
KNEE_NEAR_DEX = 0.1
R_BINS = (0, 12, 20, 30, 45, 1e9)
EULER_GAMMA = 0.5772156649015329
LN10 = math.log(10.0)

_XFIT = np.array([[1.0, d, d * d] for d in DFIT])
W_INT = np.linalg.pinv(_XFIT)[0]          # intercept weights of the quadratic fit


# ---------------------------------------------------------------- spectra
def comp_spectrum(x):
    """Per-component PSD, bins k = 1 .. N//2 - 1 (same normalisation as p5min)."""
    x = np.asarray(x, np.float64)
    n = x.shape[0]
    x = x - x.mean(axis=0)
    X = np.fft.rfft(x, axis=0)
    Pc = (np.abs(X) ** 2) * (2.0 * DT / n)
    k = np.arange(1, n // 2)
    return k / (n * DT), Pc[1:n // 2]


def prep(seg):
    """p5min's validity rule and linear gap fill; None if the window fails."""
    seg = np.asarray(seg)
    nan = np.isnan(seg).any(axis=1)
    if nan.mean() > p5min.MAX_NAN_FRAC or p5min.longest_run(nan) > p5min.MAX_NAN_RUN:
        return None
    x = np.asarray(seg, np.float64).copy()
    if nan.any():
        idx = np.arange(len(x))
        good = ~nan
        for c in range(3):
            x[nan, c] = np.interp(idx[nan], idx[good], x[good, c])
    return x


# ------------------------------------------------------- local statistics
def l_flanks(k: int, nmax: int):
    """0-based flank indices used by p5min.local_ratio for bin k (1-based)."""
    w = max(6, int(math.ceil(0.3 * k)))
    js = [j for j in range(k - w, k - 1) if j >= 1] + \
         [j for j in range(k + 2, k + w + 1) if j <= nmax]
    return np.asarray(js, dtype=np.int64) - 1


def band_L(Y, ks, nmax: int):
    """p5min's L (max over the band's bins) for every row of Y (n, K)."""
    vals = []
    for k in ks:
        fl = l_flanks(k, nmax)
        if fl.max() >= Y.shape[1]:
            raise ValueError("Y too short for flanks")
        ref = np.median(Y[:, fl], axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            v = np.where(ref > 0, Y[:, k - 1] / ref, np.nan)
        vals.append(v)
    return np.max(np.stack(vals, axis=1), axis=1)


def l2_usable(ks, K: int) -> bool:
    return all(k + DFIT[0] >= 2 and k + DFIT[-1] <= K for k in ks)


def band_L2(Y, ks):
    """L2 = y_k / exp(intercept of quadratic fit of log y over d = +-2..+-6)."""
    vals = []
    d = np.asarray(DFIT)
    with np.errstate(divide="ignore", invalid="ignore"):
        LY = np.log(Y)
    for k in ks:
        a = LY[:, k + d - 1] @ W_INT
        vals.append(Y[:, k - 1] / np.exp(a))
    return np.max(np.stack(vals, axis=1), axis=1)


def ranks_from(vals):
    """vals (n, nb), column 0 = target. Rank 1 = lowest .. nb = highest."""
    return 1 + np.sum(vals[:, 1:] < vals[:, [0]], axis=1)


def band_table(Y, lay_local, nmax: int):
    """L and L2 values (n, nb) for the usable bands, target first."""
    names = list(lay_local)
    L = np.stack([band_L(Y, lay_local[nm], nmax) for nm in names], axis=1)
    use2 = [nm for nm in names if l2_usable(lay_local[nm], Y.shape[1])]
    L2 = (np.stack([band_L2(Y, lay_local[nm]) for nm in use2], axis=1)
          if names[0] in use2 else None)
    return names, L, use2, L2


# ------------------------------------------------ the paper's statistic, batched
def s_batch(P, f, ks_band, sel_fpk, min_bins: int = p5min.MIN_BASE_BINS):
    """p5min's S for many spectra at once. P (m, nbins) trace PSD."""
    P = np.asarray(P, np.float64)
    m, nb = P.shape
    y = P * f[None, :]
    c = np.cumsum(P, axis=1)
    c = c / c[:, -1:]
    lo = np.argmax(c >= 0.25, axis=1)
    hi = np.argmax(c >= 0.75, axis=1)
    hi = np.minimum(hi, nb - 1)
    nbase = hi - lo + 1
    z = np.zeros((m, 1))
    cy = np.concatenate([z, np.cumsum(y, axis=1)], axis=1)
    cy2 = np.concatenate([z, np.cumsum(y * y, axis=1)], axis=1)
    r = np.arange(m)
    s1 = cy[r, hi + 1] - cy[r, lo]
    s2 = cy2[r, hi + 1] - cy2[r, lo]
    mu = s1 / nbase
    var = np.maximum(s2 / nbase - mu * mu, 0.0)
    sd = np.sqrt(var)
    ok = (nbase >= min_bins) & (sd > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        sb = np.max(np.stack([(y[:, k - 1] - mu) / sd for k in ks_band], axis=1), axis=1)
        ipk = sel_fpk[np.argmax(y[:, sel_fpk], axis=1)]
        sf = (y[r, ipk] - mu) / sd
    sb = np.where(ok, sb, np.nan)
    sf = np.where(ok, sf, np.nan)
    return {"S_band": sb, "S_fpk": sf, "ok": ok, "lo": lo, "hi": hi,
            "nbase": nbase, "ipk": ipk}


def smooth_cubic(f, Pc, excl=EXCL_MHZ, fmax=SMOOTH_FMAX_MHZ):
    """Smooth PSD per component: weighted (1/max(k,8)) least-squares cubic in log10 f
    of log P over f <= fmax with the excl band left out, + Euler gamma
    (E[log P] = log S - gamma for an exponential periodogram ordinate).
    Returns (S (nb,3) with S = 0 above fmax, mask of fitted bins, rms)."""
    fm = f * 1e3
    nb = len(f)
    k = np.arange(1, nb + 1, dtype=float)
    rand = fm <= fmax
    use = rand & ~((fm >= excl[0]) & (fm <= excl[1]))
    u = np.log10(f / 3e-3)
    X = np.vander(u[use], 4, increasing=True)
    wt = 1.0 / np.maximum(k[use], SMOOTH_KFLOOR)
    sw = np.sqrt(wt)
    Xa = np.vander(u[rand], 4, increasing=True)
    S = np.zeros((nb, 3))
    rms = []
    for c in range(3):
        ly = np.log(np.maximum(Pc[use, c], 1e-300))
        coef, *_ = np.linalg.lstsq(X * sw[:, None], ly * sw, rcond=None)
        res = ly - X @ coef
        rms.append(float(np.sqrt(np.sum(wt * res ** 2) / np.sum(wt))))
        S[rand, c] = np.exp(Xa @ coef + EULER_GAMMA)
    return S, rand, rms


def surrogate_test(rng, f, Pc, ks_band, nsurr: int, chunk: int = CHUNK):
    """Gaussian-noise surrogates around the smooth fit; bins above
    SMOOTH_FMAX_MHZ keep the observed values. Returns observed S and the
    surrogate exceedance counts."""
    sel = np.flatnonzero((f >= FPK_MHZ[0] * 1e-3) & (f <= FPK_MHZ[1] * 1e-3))
    obs = s_batch(Pc.sum(axis=1)[None, :], f, ks_band, sel)
    S, rand, rms = smooth_cubic(f, Pc)
    nr = int(rand.sum())
    fixed = Pc.sum(axis=1).copy()
    sb_all, sf_all = [], []
    done = 0
    while done < nsurr:
        m = min(chunk, nsurr - done)
        E = rng.exponential(size=(m, nr, 3))
        P = np.empty((m, len(f)))
        P[:, :] = fixed[None, :]
        P[:, :nr] = (E * S[None, :nr, :]).sum(axis=2)
        r = s_batch(P, f, ks_band, sel)
        sb_all.append(r["S_band"])
        sf_all.append(r["S_fpk"])
        done += m
    sb = np.concatenate(sb_all)
    sf = np.concatenate(sf_all)
    out = {"S_band_obs": float(obs["S_band"][0]), "S_fpk_obs": float(obs["S_fpk"][0]),
           "fpk_mhz": float(f[obs["ipk"][0]] * 1e3),
           "n_base": int(obs["nbase"][0]),
           "base_lo_mhz": float(f[obs["lo"][0]] * 1e3),
           "base_hi_mhz": float(f[obs["hi"][0]] * 1e3),
           "peak_in_base": bool(obs["lo"][0] <= obs["ipk"][0] <= obs["hi"][0]),
           "samuelson_ceiling": float(math.sqrt(max(int(obs["nbase"][0]) - 1, 0))),
           "smooth_rms_ln": [round(v, 4) for v in rms]}
    for key, arr, o in (("band", sb, out["S_band_obs"]), ("fpk", sf, out["S_fpk_obs"])):
        okm = np.isfinite(arr)
        n_ok = int(okm.sum())
        hits = int(np.sum(arr[okm] >= o)) if math.isfinite(o) else -1
        out[f"n_ok_{key}"] = n_ok
        out[f"hits_{key}"] = hits
        out[f"p_{key}"] = ((hits + 1) / (n_ok + 1)) if (hits >= 0 and n_ok) else float("nan")
        out[f"z_{key}"] = gauss_isf(out[f"p_{key}"])
        if n_ok:
            out[f"q999_{key}"] = float(np.percentile(arr[okm], 99.9))
            out[f"q99_{key}"] = float(np.percentile(arr[okm], 99.0))
    return out


def gauss_isf(p: float) -> float:
    """One-sided Gaussian z with P(Z >= z) = p (bisection on erfc)."""
    if not (isinstance(p, float) and math.isfinite(p)) or p <= 0 or p >= 1:
        return float("nan")
    lo, hi = -10.0, 40.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if 0.5 * math.erfc(mid / math.sqrt(2.0)) > p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ------------------------------------------------------------ simulations
def sim_Y(rng, n: int, f, kind: str):
    fk = f[:KMAX]
    if kind == "F":
        S = fk ** -1.0
        w = np.array([1.0, 1.0, 1.0])
    elif kind == "K":
        fb = 10 ** rng.uniform(-3.0, -2.0, size=(n, 1))
        S = fk[None, :] ** -1.0 * (1.0 + (fk[None, :] / fb) ** 2) ** -0.3
        w = np.array([0.2, 1.0, 1.0])
    elif kind == "K0":
        S = np.where(fk <= 3.2e-3, fk ** -1.0, (3.2e-3 ** -1.0) * (fk / 3.2e-3) ** -1.6)
        w = np.array([0.2, 1.0, 1.0])
    else:
        raise ValueError(kind)
    E = rng.exponential(size=(3, n, KMAX))
    P = np.tensordot(w, E, axes=1) * S
    return P * fk


def run_sims(rng, n: int, f, lay_local, nmax: int, kinds=("F", "K", "K0"), nrep: int = N_REP):
    out = {}
    for kind in kinds:
        zL, zL2 = [], []
        for _ in range(nrep):
            Y = sim_Y(rng, n, f, kind)
            names, L, use2, L2 = band_table(Y, lay_local, nmax)
            zL.append(p5min.rank_z(ranks_from(L), L.shape[1])["z"])
            if L2 is not None:
                zL2.append(p5min.rank_z(ranks_from(L2), L2.shape[1])["z"])
        zL, zL2 = np.asarray(zL), np.asarray(zL2)
        out[kind] = {"n_windows": n, "n_rep": nrep,
                     "Z_L_mean": float(zL.mean()), "Z_L_sd": float(zL.std(ddof=1)),
                     "Z_L2_mean": float(zL2.mean()) if len(zL2) else float("nan"),
                     "Z_L2_sd": float(zL2.std(ddof=1)) if len(zL2) > 1 else float("nan")}
    return out


# ------------------------------------------------------------ spectral knee
def knee_fit(f, P):
    """Smooth broken power law per row of P (n, K) over KNEE_FIT_MHZ, grid in
    the knee. Model: log10 P = a + b1 u' + (b2 - b1) D log10(1 + 10^(u'/D)),
    u' = log10 f - u_b. Returns arrays fknee_mhz, b1, b2."""
    K = P.shape[1]
    fk = f[:K]
    sel = np.flatnonzero((fk >= KNEE_FIT_MHZ[0] * 1e-3) & (fk <= KNEE_FIT_MHZ[1] * 1e-3))
    u = np.log10(fk[sel])
    sw = np.sqrt(1.0 / (sel + 1.0))
    with np.errstate(divide="ignore"):
        Yw = (np.log10(P[:, sel]) * sw[None, :]).T      # (m, n)
    n = P.shape[0]
    best = np.full(n, np.inf)
    ub_best = np.full(n, np.nan)
    b1_best = np.full(n, np.nan)
    b2_best = np.full(n, np.nan)
    for ub in KNEE_GRID:
        up = u - ub
        sp = KNEE_DELTA * np.logaddexp(0.0, up / KNEE_DELTA * LN10) / LN10
        X = np.stack([np.ones_like(u), up, sp], axis=1) * sw[:, None]
        coef, *_ = np.linalg.lstsq(X, Yw, rcond=None)
        sse = np.sum((Yw - X @ coef) ** 2, axis=0)
        better = sse < best
        best[better] = sse[better]
        ub_best[better] = ub
        b1_best[better] = coef[1, better]
        b2_best[better] = coef[1, better] + coef[2, better]
    return 10 ** ub_best * 1e3, b1_best, b2_best


# ------------------------------------------------------------ data access
def load_bins(enc, peri_t, data_dir, priv_dir, log):
    g0 = math.floor((peri_t - p5min.HALF_SPAN_DAYS * 86400) / 60.0) * 60.0
    nbin = int(round(2 * p5min.HALF_SPAN_DAYS * 86400 / DT))
    cache = os.path.join(data_dir, "bins", f"E{enc:02d}_{int(g0)}.npy")
    if os.path.exists(cache):
        try:
            B = np.load(cache)
            if B.shape == (nbin, 3):
                return g0, np.asarray(B), "shared-cache"
            log(f"PCURV_WARN: E{enc} cache shape {B.shape}; rebuilding privately")
        except Exception as e:  # noqa: BLE001
            log(f"PCURV_WARN: E{enc} cache unreadable ({e!r}); rebuilding privately")
    g0b, B, meta = p5min.encounter_bins(enc, peri_t, priv_dir, False, log)
    try:
        shutil.rmtree(os.path.join(priv_dir, "bins"), ignore_errors=True)
    except Exception:  # noqa: BLE001
        pass
    return g0b, B, "private-download:" + json.dumps(meta)


def r_at(t, pos_t, pos_r):
    return float(np.interp(t, pos_t, pos_r)) * p5min.RS_PER_AU


# ------------------------------------------------------------ summaries
def rank_tables(rows, ranksL, ranksL2, nbL, nbL2):
    sets = {"all": np.ones(len(rows), bool),
            "r_le_30": np.array([r["r_rs"] <= 30 for r in rows], bool),
            "A_E1_E24": np.array([r["enc"] <= p5min.AUTHOR_LAST_ENC for r in rows], bool),
            "B_E25_E27": np.array([r["enc"] > p5min.AUTHOR_LAST_ENC for r in rows], bool)}
    out = {"L": {}, "L2": {}}
    for nm, m in sets.items():
        out["L"][nm] = p5min.rank_z(ranksL[m], nbL) if len(rows) else p5min.rank_z([], nbL)
        out["L2"][nm] = (p5min.rank_z(ranksL2[m], nbL2) if (ranksL2 is not None and len(rows))
                         else p5min.rank_z([], 9))
    return out


def knee_summary(rows, fkn, b1, b2, ranksL2):
    out = {"n": int(np.isfinite(fkn).sum())}
    if not out["n"]:
        return out
    r = np.array([x["r_rs"] for x in rows])
    lf = np.log10(fkn / 3.2)
    near = np.abs(lf) <= KNEE_NEAR_DEX
    edge = (np.isclose(fkn, 10 ** KNEE_GRID[0] * 1e3)) | (np.isclose(fkn, 10 ** KNEE_GRID[-1] * 1e3))
    out["frac_knee_2p6_3p9"] = float(np.mean((fkn >= 2.6) & (fkn <= 3.9)))
    out["frac_knee_within_0p1dex_of_3p2"] = float(near.mean())
    out["frac_at_grid_edge"] = float(edge.mean())
    by_r = []
    for a, b in zip(R_BINS[:-1], R_BINS[1:]):
        m = (r > a) & (r <= b) & np.isfinite(fkn)
        if m.sum():
            by_r.append({"r_lo": a, "r_hi": b if b < 1e8 else None, "n": int(m.sum()),
                         "median_fknee_mhz": float(np.median(fkn[m])),
                         "median_slope_lo": float(np.median(b1[m])),
                         "median_slope_hi": float(np.median(b2[m]))})
    out["by_r"] = by_r
    if ranksL2 is not None:
        for nm, m in (("near", near), ("far", ~near)):
            out[f"mean_rank_L2_knee_{nm}"] = float(ranksL2[m].mean()) if m.sum() else float("nan")
            out[f"n_knee_{nm}"] = int(m.sum())
    return out


def grade(res):
    g = {}
    sims = res.get("sims", {})

    def simok(kind):
        s = sims.get(kind)
        return bool(s) and math.isfinite(s.get("Z_L2_mean", float("nan"))) \
            and abs(s["Z_L2_mean"]) <= SIM_GATE

    gate = simok("F") and simok("K")
    g["SIMGATE"] = "PASS" if gate else "FAIL"
    rk = res.get("D135", {}).get("rank", {}).get("L2", {}).get("all", {})
    z, n = rk.get("z", float("nan")), rk.get("n", 0)
    if not gate or not (isinstance(z, float) and math.isfinite(z)) or n < MIN_N_H2C:
        g["H2c"] = "NOT GRADABLE"
    elif z >= 3.0:
        g["H2c"] = "SUPPORT"
    elif z < 2.0:
        g["H2c"] = "NO-EXCESS"
    else:
        g["H2c"] = "INDETERMINATE"
    pubs = res.get("published", [])
    ps = [p.get("p_fpk", float("nan")) for p in pubs]
    nok = [p.get("n_ok_fpk", 0) for p in pubs]
    if len(ps) != 2 or not all(isinstance(p, float) and math.isfinite(p) for p in ps) \
            or min(nok) < 10000:
        g["H4"] = "NOT GRADABLE"
    elif any(p > 1e-2 for p in ps):
        g["H4"] = "ORDINARY"
    elif all(p < 1e-4 for p in ps):
        g["H4"] = "EXTREME"
    else:
        g["H4"] = "MODERATE"
    return g


def fmt(v):
    if isinstance(v, (bool, np.bool_)):
        return int(v)
    if isinstance(v, np.integer):
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
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: fmt(r.get(k, "")) for k in keys})


def jdefault(v):
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating,)):
        return float(v)
    if isinstance(v, np.bool_):
        return bool(v)
    if isinstance(v, np.ndarray):
        return v.tolist()
    return str(v)


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


# ------------------------------------------------------------ main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--budget-s", type=float, default=13000)
    ap.add_argument("--n-rep", type=int, default=N_REP)
    ap.add_argument("--only", default="")
    a = ap.parse_args(argv)
    t_start = time.time()
    os.makedirs(a.out_dir, exist_ok=True)
    priv = tempfile.mkdtemp(prefix="psp5curv-")
    rng = np.random.default_rng(SEED)

    def log(msg):
        print(msg, flush=True)

    log(f"PCURV_VERSION: {VERSION} p5min {p5min.VERSION} numpy {np.__version__} "
        f"self_sha256 {sha256_file(os.path.abspath(__file__))} "
        f"p5min_sha256 {sha256_file(os.path.abspath(p5min.__file__))}")
    log("PCURV_W_INT: " + json.dumps([round(float(v), 6) for v in W_INT]))
    lays = {d: p5min.layout(d * 60) for d in DURATIONS_MIN}

    pos_t, pos_r = p5min.fetch_positions()
    peris = p5min.find_perihelia(pos_t, pos_r)
    allp = [(i + 1, t, r) for i, (t, r) in enumerate(peris)]
    mag_stop = p5min.to_unix(p5min.MAG_STOP)
    num_ok = all(any(n == e and abs(t - p5min.to_unix(x)) <= 2 * 86400 for n, t, r in allp)
                 for e, x in p5min.EXPECTED_PERI.items())
    log(f"PCURV_PERIHELIA: n={len(allp)} numbering_ok={num_ok}")
    encs = [(e, t) for e, t, r in allp if t + p5min.HALF_SPAN_DAYS * 86400 <= mag_stop]
    if a.only:
        want = {int(x) for x in a.only.split(",") if x.strip()}
        encs = [x for x in encs if x[0] in want]
    encs.sort(key=lambda x: (0, x[0]) if x[0] >= 22 else (1, -x[0]))

    rows = {d: [] for d in DURATIONS_MIN}
    Ys = {d: [] for d in DURATIONS_MIN}
    fgrid, nmaxs = {}, {}
    published, enc_info = [], []
    for enc, pt in encs:
        if time.time() - t_start > a.budget_s:
            log(f"PCURV_BUDGET: stopping before E{enc}")
            break
        t0 = time.time()
        try:
            g0, B, src = load_bins(enc, pt, a.data_dir, priv, log)
        except Exception as e:  # noqa: BLE001
            log(f"PCURV_WARN: E{enc} bins failed: {e!r}")
            continue
        nwin = {}
        for dmin in DURATIONS_MIN:
            n = int(round(dmin * 60 / DT))
            cnt = 0
            for s in range(0, B.shape[0] - n + 1, n):
                ts = g0 + s * DT
                te = ts + n * DT
                if p5min.overlaps_pub(ts, te):
                    continue
                x = prep(B[s:s + n])
                if x is None:
                    continue
                f, Pc = comp_spectrum(x)
                y = Pc.sum(axis=1) * f
                fgrid[dmin] = f
                nmaxs[dmin] = len(f)
                Ys[dmin].append(y[:KMAX].astype(np.float64))
                rows[dmin].append({"enc": enc, "D": dmin, "start": p5min.iso(ts),
                                   "r_rs": round(r_at(ts + n * DT / 2, pos_t, pos_r), 3)})
                cnt += 1
            nwin[f"D{dmin}"] = cnt
        for p in p5min.PUBLISHED:
            if p["enc"] != enc:
                continue
            ta, tb = p5min.to_unix(p["start"]), p5min.to_unix(p["stop"])
            dmin = int(round((tb - ta) / 60))
            s = int(round((ta - g0) / DT))
            n = int(round((tb - ta) / DT))
            rec = {"enc": enc, "D": dmin, "start": p["start"], "stop": p["stop"],
                   "published_sigma": p["sigma"], "published_f_mhz": p["f_mhz"]}
            x = prep(B[s:s + n])
            if x is None:
                rec["valid"] = 0
            else:
                rec["valid"] = 1
                f, Pc = comp_spectrum(x)
                ref = p5min.analyse_window(B[s:s + n], lays[dmin])
                rec["p5min_S_fpk"] = ref.get("S_fpk", float("nan"))
                rec["p5min_S_P3.15"] = ref.get("S_P3.15", float("nan"))
                rec.update(surrogate_test(rng, f, Pc, lays[dmin]["paper"]["P3.15"],
                                          N_SURR[dmin]))
                y = (Pc.sum(axis=1) * f)[None, :KMAX]
                names, L, use2, L2 = band_table(y, lays[dmin]["local"], len(f))
                rec["L_T3.2"] = float(L[0, 0])
                rec["rank_L"] = int(ranks_from(L)[0])
                if L2 is not None:
                    rec["L2_T3.2"] = float(L2[0, 0])
                    rec["rank_L2"] = int(ranks_from(L2)[0])
            published.append(rec)
            log("PCURV_PUBLISHED: " + json.dumps(rec, default=jdefault))
        info = {"enc": enc, "source": src, **nwin, "seconds": round(time.time() - t0, 1)}
        enc_info.append(info)
        log("PCURV_ENC: " + json.dumps(info))
        del B

    res = {"version": VERSION, "p5min_sha": P5MIN_SHA, "numbering_ok": num_ok,
           "encounters": enc_info, "published": published}
    win_rows = []
    rk135 = None
    for dmin in DURATIONS_MIN:
        dres = {"n_windows": len(rows[dmin])}
        if rows[dmin]:
            Y = np.stack(Ys[dmin])
            names, L, use2, L2 = band_table(Y, lays[dmin]["local"], nmaxs[dmin])
            rL = ranks_from(L)
            rL2 = ranks_from(L2) if L2 is not None else None
            dres["bands_L"] = names
            dres["bands_L2"] = use2
            dres["rank"] = rank_tables(rows[dmin], rL, rL2, L.shape[1],
                                       L2.shape[1] if L2 is not None else 9)
            for i, r in enumerate(rows[dmin]):
                r["rank_L"] = int(rL[i])
                r["L_T"] = float(L[i, 0])
                if rL2 is not None:
                    r["rank_L2"] = int(rL2[i])
                    r["L2_T"] = float(L2[i, 0])
            if dmin == 135:
                rk135 = rL2
                fk, b1, b2 = knee_fit(fgrid[dmin], Y / fgrid[dmin][None, :KMAX])
                for i, r in enumerate(rows[dmin]):
                    r["fknee_mhz"] = float(fk[i])
                    r["slope_lo"] = float(b1[i])
                    r["slope_hi"] = float(b2[i])
                dres["knee"] = knee_summary(rows[dmin], fk, b1, b2, rk135)
            del Y
        res[f"D{dmin}"] = dres
        win_rows += rows[dmin]
        log(f"PCURV_D{dmin}: " + json.dumps(dres, default=jdefault))

    f135 = fgrid.get(135)
    if f135 is None:
        f135 = np.arange(1, 8100) / 8100.0
        nmax135 = 8099
    else:
        nmax135 = nmaxs[135]
    nsim = len(rows[135]) if len(rows[135]) >= 100 else SIM_N_FALLBACK
    res["sims"] = run_sims(rng, nsim, f135, lays[135]["local"], nmax135, nrep=a.n_rep)
    log("PCURV_SIMS: " + json.dumps(res["sims"], default=jdefault))
    res["grades"] = grade(res)
    res["runtime_s"] = round(time.time() - t_start, 1)
    with open(os.path.join(a.out_dir, "results.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=jdefault)
    write_rows(os.path.join(a.out_dir, "windows.csv"), win_rows)
    write_rows(os.path.join(a.out_dir, "published.csv"), published)
    log("PCURV_GRADES: " + json.dumps(res["grades"], default=jdefault))
    shutil.rmtree(priv, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
