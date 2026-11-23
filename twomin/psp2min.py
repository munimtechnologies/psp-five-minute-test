#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. Written by AI agents working for MV2 of Munim, Inc.
"""psp2min 0.1.0 - out-of-sample test of the "2-minute dominance" of
Parker Solar Probe magnetic fluctuations near the Alfven surface.

Claim tested: Huang, Velli, Shi et al., "Dominance of 2-Minute
Oscillations near the Alfven Surface", ApJL 977 L12 (2024),
arXiv:2405.15967. Using 4-hour windows from encounters E1-E19, they
found that the period 1/f_mid at the middle of the band holding the
central 50 % of the magnetic fluctuation energy falls as the solar-wind
advection time t_adv = (R - Rsun)/V_r falls, bottoming out at about
2-3 minutes in the youngest wind, and argue the true floor is "likely
closer to 2 minutes": the launch period of Alfven waves from the solar
atmosphere. This script measures the same quantity in encounters
E20-E27 (2024-06 to 2026-03; E22-E27 at 9.86 Rsun, the deepest and
youngest wind ever sampled), which the paper never used, and repeats it
on E1-E19 as a reproduction check. Rules: twomin/PREREG-2MIN-v0.1.md.

Reuses, unchanged and pinned: psp5min/p5min.py 0.1.0 (positions,
perihelion numbering, FIELDS magnetometer binning to 0.5 s) and
alfven/psp5alf.py 0.1.0 (SWEAP SPAN-i moments via HAPI).
Data (public, NASA CDAWeb): PSP_FLD_L2_MAG_RTN_4_SA_PER_CYC,
PSP_SWP_SPI_SF00_L3_MOM (VEL_RTN_SUN), PSP_HELIO1HR_POSITION.
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

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in (os.path.join(os.path.dirname(_HERE), "psp5min"),
           os.path.join(os.path.dirname(_HERE), "alfven"), _HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import p5min as P  # noqa: E402

VERSION = "0.1.0"
DT = 0.5                       # s, p5min grid
WIN_S = 4 * 3600               # 4-hour windows, as in the paper
RS_KM = 695700.0
# Lowest trusted frequency: where 50 % of a Morlet (w0 = 6) wavelet's
# points fall outside the cone of influence for a window of length T:
# scale s = T / (4 sqrt 2), Fourier period 1.0330 s.
MORLET_FOURIER = 4 * math.pi / (6 + math.sqrt(2 + 36))
F_HI = 1.0                     # Hz, Nyquist of the 0.5 s grid
FPEAK_MAX = 0.05               # Hz, upper end of the f_peak search
MAX_NAN_FRAC = 0.05
MAX_NAN_RUN = 120              # 0.5 s bins = 60 s
V_MIN, V_MAX = 50.0, 2000.0    # km/s, accepted SPAN-i radial speeds
V_COVER_MIN = 0.6              # fraction of 1-minute slots with a V_r sample
DTADV_LO, DTADV_HI = 0.005, 0.2
CLEAN_MAX_REMOVED = 0.15       # if the cut removes more in E1-E19, grade without it
TADV_YOUNG = 10.0              # h
TADV_SPLIT = 7.0               # h
TADV_EDGES = (0, 6, 7, 8, 10, 12, 15, 20, 25, 30, 40, 60, 1e9)
AUTHOR_LAST_ENC = 19
N_BOOT = 10000
SEED = 20261123
UA = "psp2min/0.1 (open research script)"


def f_lo_for(win_s: float) -> float:
    return 1.0 / (MORLET_FOURIER * win_s / (4 * math.sqrt(2)))


# ------------------------------------------------------------- spectra
def detrend_taper(x):
    """Remove a least-squares line from each column, apply a Hann taper."""
    x = np.asarray(x, np.float64)
    n = x.shape[0]
    t = np.arange(n, dtype=np.float64)
    tc = t - t.mean()
    out = x - x.mean(axis=0)
    slope = (tc[:, None] * out).sum(axis=0) / (tc * tc).sum()
    out = out - tc[:, None] * slope[None, :]
    w = np.hanning(n)
    return out * w[:, None], float((w * w).mean())


def trace_psd(x):
    """Hann-tapered, detrended trace PSD of an (N,3) series at DT."""
    xw, wn = detrend_taper(x)
    n = xw.shape[0]
    X = np.fft.rfft(xw, axis=0)
    Pw = (np.abs(X) ** 2).sum(axis=1) * (2.0 * DT / (n * wn))
    f = np.arange(len(Pw)) / (n * DT)
    return f[1:], Pw[1:]


def fmid_stats(f, Pw, f_lo, f_hi=F_HI):
    """f25, f75 of the cumulative energy between f_lo and f_hi; f_mid
    geometric (primary) and arithmetic; f_peak of factor-2-smoothed P*f."""
    sel = (f >= f_lo) & (f <= f_hi)
    fs, ps = f[sel], Pw[sel]
    if len(fs) < 20 or not np.all(np.isfinite(ps)) or ps.sum() <= 0:
        return None
    E = np.cumsum(ps)
    tot = E[-1]
    i25 = int(np.searchsorted(E, 0.25 * tot))
    i75 = int(np.searchsorted(E, 0.75 * tot))
    f25, f75 = float(fs[i25]), float(fs[min(i75, len(fs) - 1)])
    y = ps * fs
    cs = np.concatenate(([0.0], np.cumsum(y)))
    lo = np.searchsorted(fs, fs / math.sqrt(2.0), side="left")
    hi = np.searchsorted(fs, fs * math.sqrt(2.0), side="right")
    ys = (cs[hi] - cs[lo]) / np.maximum(hi - lo, 1)
    pk = fs <= FPEAK_MAX
    ipk = int(np.argmax(np.where(pk, ys, -np.inf)))
    return {"f25": f25, "f75": f75,
            "fmid_geo": math.sqrt(f25 * f75), "fmid_ari": 0.5 * (f25 + f75),
            "fpeak": float(fs[ipk]), "fpeak_at_floor": int(ipk == 0)}


def fill_gaps(seg):
    nan = np.isnan(seg).any(axis=1)
    frac = float(nan.mean())
    run = P.longest_run(nan)
    if frac > MAX_NAN_FRAC or run > MAX_NAN_RUN:
        return None, frac, run
    x = np.asarray(seg, np.float64).copy()
    if nan.any():
        idx = np.arange(len(x))
        good = ~nan
        for c in range(3):
            x[nan, c] = np.interp(idx[nan], idx[good], x[good, c])
    return x, frac, run


# ------------------------------------------------------------- plasma
def window_plasma(ts, vr, t0, t1, pos_t, pos_r):
    """Median V_r, t_adv (h), the 10-minute-median spread ratio of t_adv,
    and the fraction of 1-minute slots holding a sample."""
    ts = np.asarray(ts, float)
    vr = np.asarray(vr, float)
    m = (ts >= t0) & (ts < t1) & np.isfinite(vr) & (vr > V_MIN) & (vr < V_MAX)
    out = {"v_n": int(m.sum()), "v_cover": 0.0, "vr": float("nan"),
           "tadv_h": float("nan"), "dtadv_ratio": float("nan")}
    if m.sum() < 10:
        return out
    tt, vv = ts[m], vr[m]
    nslot = int(round((t1 - t0) / 60.0))
    slots = np.unique(np.floor((tt - t0) / 60.0).astype(np.int64))
    out["v_cover"] = float(len(slots)) / nslot
    vmed = float(np.median(vv))
    rc = float(np.interp(0.5 * (t0 + t1), pos_t, pos_r)) * P.RS_PER_AU
    out["vr"] = vmed
    out["r_rs"] = rc
    out["tadv_h"] = (rc - 1.0) * RS_KM / vmed / 3600.0
    blk = np.floor((tt - t0) / 600.0).astype(np.int64)
    meds = []
    for b in np.unique(blk):
        vb = vv[blk == b]
        if len(vb) >= 5:
            tb = t0 + (b + 0.5) * 600.0
            rb = float(np.interp(tb, pos_t, pos_r)) * P.RS_PER_AU
            meds.append((rb - 1.0) * RS_KM / float(np.median(vb)) / 3600.0)
    if len(meds) >= 6:
        meds = np.asarray(meds)
        out["dtadv_ratio"] = float(meds.std() / meds.mean())
    return out


# ------------------------------------------------------------- statistics
def rankdata(a):
    a = np.asarray(a, float)
    o = np.argsort(a, kind="mergesort")
    r = np.empty(len(a), float)
    r[o] = np.arange(1, len(a) + 1)
    s = a[o]
    i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and s[j + 1] == s[i]:
            j += 1
        if j > i:
            r[o[i:j + 1]] = 0.5 * (i + j) + 1
        i = j + 1
    return r


def spearman(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 5:
        return float("nan")
    rx, ry = rankdata(x[ok]), rankdata(y[ok])
    rx -= rx.mean()
    ry -= ry.mean()
    d = math.sqrt((rx * rx).sum() * (ry * ry).sum())
    return float((rx * ry).sum() / d) if d > 0 else float("nan")


def gmean(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v) & (v > 0)]
    return float(np.exp(np.log(v).mean())) if len(v) else float("nan")


def block_boot_gmean(groups, n_boot=N_BOOT, seed=SEED):
    """Encounter-block bootstrap of the pooled geometric mean.
    groups: {enc: array of positive values}. Returns (lo95, hi95)."""
    keys = sorted(k for k, v in groups.items() if len(v))
    if len(keys) < 2:
        return float("nan"), float("nan")
    sums = np.array([np.log(np.asarray(groups[k], float)).sum() for k in keys])
    cnts = np.array([len(groups[k]) for k in keys], float)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(keys), size=(n_boot, len(keys)))
    gm = np.exp(sums[idx].sum(axis=1) / cnts[idx].sum(axis=1))
    return float(np.percentile(gm, 2.5)), float(np.percentile(gm, 97.5))


def trend_table(rows, key="pmid_s"):
    out = []
    for lo, hi in zip(TADV_EDGES[:-1], TADV_EDGES[1:]):
        v = [r[key] for r in rows if lo <= r["tadv_h"] < hi]
        lv = np.log(np.asarray(v, float)) if v else np.zeros(0)
        out.append({"tadv_lo": lo, "tadv_hi": hi if hi < 1e8 else "inf", "n": len(v),
                    "gmean_s": gmean(v), "log_sd": float(lv.std()) if len(v) > 1 else float("nan")})
    return out


# ------------------------------------------------------------- grading
def grade_r1(rho, g_young, n_young):
    if n_young < 20 or not math.isfinite(rho):
        return "NOT GRADABLE"
    return "PASS" if (rho > 0.3 and 90.0 <= g_young <= 240.0) else "FAIL"


def grade_h1(g, n, n_enc):
    if n < 30 or n_enc < 3 or not math.isfinite(g):
        return "NOT GRADABLE"
    if g < 90.0:
        return "SHORTER"
    if g > 180.0:
        return "LONGER"
    return "CONSISTENT"


def grade_h2(g_lo, n_lo, g_hi, n_hi):
    if n_lo < 20 or n_hi < 20 or not (math.isfinite(g_lo) and math.isfinite(g_hi)):
        return "NOT GRADABLE", float("nan")
    ratio = g_lo / g_hi
    if ratio >= 0.85:
        return "FLOOR", ratio
    if ratio <= 0.70:
        return "DECLINING", ratio
    return "INTERMEDIATE", ratio


# ------------------------------------------------------------- pipeline
def analyse_encounter(enc, peri_t, B, g0, spi_t, spi_vr, pos_t, pos_r):
    n = int(round(WIN_S / DT))
    f_lo = f_lo_for(WIN_S)
    rows = []
    for s in range(0, B.shape[0] - n + 1, n):
        t0 = g0 + s * DT
        t1 = t0 + WIN_S
        rec = {"enc": enc, "start": P.iso(t0),
               "hours_from_peri": round((0.5 * (t0 + t1) - peri_t) / 3600.0, 3)}
        x, frac, run = fill_gaps(B[s:s + n])
        rec["nan_frac"] = round(frac, 5)
        rec["mag_ok"] = int(x is not None)
        pl = window_plasma(spi_t, spi_vr, t0, t1, pos_t, pos_r)
        rec.update({k: pl[k] for k in pl})
        if "r_rs" not in rec:
            rec["r_rs"] = float(np.interp(0.5 * (t0 + t1), pos_t, pos_r)) * P.RS_PER_AU
        rec["v_ok"] = int(pl["v_cover"] >= V_COVER_MIN and math.isfinite(pl["tadv_h"]))
        d = pl["dtadv_ratio"]
        rec["clean_ok"] = int(math.isfinite(d) and DTADV_LO < d < DTADV_HI)
        if x is not None:
            f, Pw = trace_psd(x)
            st = fmid_stats(f, Pw, f_lo)
            if st:
                rec.update({"f25_mhz": st["f25"] * 1e3, "f75_mhz": st["f75"] * 1e3,
                            "pmid_s": 1.0 / st["fmid_geo"], "pmid_ari_s": 1.0 / st["fmid_ari"],
                            "ppeak_s": 1.0 / st["fpeak"], "fpeak_at_floor": st["fpeak_at_floor"]})
        rec["usable"] = int(rec["mag_ok"] and rec["v_ok"] and "pmid_s" in rec)
        rows.append(rec)
    return rows


def kolmogorov_check(seed=5):
    """1/f_mid of synthetic f^-5/3 noise in a 4-hour window (paper: ~13 min)."""
    n = int(round(WIN_S / DT))
    rng = np.random.default_rng(seed)
    f = np.fft.rfftfreq(n, d=DT)
    amp = np.zeros_like(f)
    amp[1:] = f[1:] ** (-5.0 / 6.0)
    x = np.zeros((n, 3))
    for c in range(3):
        x[:, c] = np.fft.irfft(amp * np.exp(1j * rng.uniform(0, 2 * np.pi, len(f))), n=n)
    ff, Pw = trace_psd(x)
    return 1.0 / fmid_stats(ff, Pw, f_lo_for(WIN_S))["fmid_geo"]


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


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
            w.writerow({k: P.fmt(r.get(k, "")) for k in keys})


def summarise(rows, encs_done):
    res = {"version": VERSION, "encounters_done": sorted(encs_done)}
    a_all = [r for r in rows if r["enc"] <= AUTHOR_LAST_ENC and r["usable"]]
    a_clean = [r for r in a_all if r["clean_ok"]]
    removed = 1.0 - len(a_clean) / len(a_all) if a_all else float("nan")
    use_cut = bool(a_all) and removed <= CLEAN_MAX_REMOVED
    res["clean_cut"] = {"E1_E19_usable": len(a_all), "after_cut": len(a_clean),
                        "removed_frac": removed, "cut_used_for_grades": use_cut,
                        "paper_removed_frac": 982 / 26623}

    def sample(lo, hi):
        s = [r for r in rows if lo <= r["enc"] <= hi and r["usable"]]
        return [r for r in s if r["clean_ok"]] if use_cut else s

    A = sample(1, AUTHOR_LAST_ENC)
    B = sample(AUTHOR_LAST_ENC + 1, 99)
    res["n_usable"] = {"A_E1_E19": len(A), "B_E20_E27": len(B)}
    res["trend_A"] = trend_table(A)
    res["trend_B"] = trend_table(B)
    res["trend_B_arith"] = trend_table(B, "pmid_ari_s")
    # R1 reproduction on the paper's own encounters
    rhoA = spearman([r["tadv_h"] for r in A], [r["pmid_s"] for r in A])
    yA = [r["pmid_s"] for r in A if r["tadv_h"] <= TADV_YOUNG]
    res["R1"] = {"spearman_tadv_pmid": rhoA, "n_young": len(yA), "gmean_young_s": gmean(yA),
                 "min_tadv_h": min((r["tadv_h"] for r in A), default=float("nan"))}
    # H1 out-of-sample level in the youngest wind
    yB = {}
    for r in B:
        if r["tadv_h"] <= TADV_YOUNG:
            yB.setdefault(r["enc"], []).append(r["pmid_s"])
    vals = [v for vs in yB.values() for v in vs]
    lo95, hi95 = block_boot_gmean(yB)
    res["H1"] = {"n_young": len(vals), "n_encounters": len(yB), "gmean_s": gmean(vals),
                 "boot95_s": [lo95, hi95],
                 "per_encounter": {str(k): {"n": len(v), "gmean_s": gmean(v)} for k, v in sorted(yB.items())},
                 "gmean_arith_def_s": gmean([r["pmid_ari_s"] for r in B if r["tadv_h"] <= TADV_YOUNG]),
                 "spearman_tadv_pmid_B": spearman([r["tadv_h"] for r in B], [r["pmid_s"] for r in B]),
                 "min_tadv_h": min((r["tadv_h"] for r in B), default=float("nan"))}
    # H2 floor or continued decline (all encounters pooled)
    AB = A + B
    lo = [r["pmid_s"] for r in AB if r["tadv_h"] <= TADV_SPLIT]
    hi = [r["pmid_s"] for r in AB if TADV_SPLIT < r["tadv_h"] <= TADV_YOUNG]
    g2, ratio = grade_h2(gmean(lo), len(lo), gmean(hi), len(hi))
    res["H2"] = {"n_le7": len(lo), "gmean_le7_s": gmean(lo), "n_7_10": len(hi),
                 "gmean_7_10_s": gmean(hi), "ratio": ratio,
                 "B_only_ratio": (gmean([r["pmid_s"] for r in B if r["tadv_h"] <= TADV_SPLIT])
                                  / gmean([r["pmid_s"] for r in B if TADV_SPLIT < r["tadv_h"] <= TADV_YOUNG]))
                 if B else float("nan")}
    res["fpeak_at_floor_frac_B_young"] = (float(np.mean([r["fpeak_at_floor"] for r in B if r["tadv_h"] <= TADV_YOUNG]))
                                          if vals else float("nan"))
    res["gmean_ppeak_B_young_s"] = gmean([r["ppeak_s"] for r in B if r["tadv_h"] <= TADV_YOUNG])
    res["grades"] = {"R1": grade_r1(rhoA, res["R1"]["gmean_young_s"], len(yA)),
                     "H1": grade_h1(res["H1"]["gmean_s"], len(vals), len(yB)),
                     "H2": g2}
    return res


def grade_p0(enc_info, numbering_ok):
    a = sum(1 for e in enc_info if e["enc"] <= AUTHOR_LAST_ENC and e["usable_frac"] >= 0.5)
    b = sum(1 for e in enc_info if e["enc"] > AUTHOR_LAST_ENC and e["usable_frac"] >= 0.5)
    return ("PASS" if (numbering_ok and a >= 15 and b >= 6) else "FAIL"), a, b


def main(argv=None):
    import psp5alf as A  # noqa: PLC0415  (network-facing helpers, unchanged)
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--budget-s", type=float, default=13000)
    ap.add_argument("--only", default="")
    a = ap.parse_args(argv)
    t_start = time.time()
    os.makedirs(a.out_dir, exist_ok=True)

    def log(msg):
        print(msg, flush=True)

    log(f"P2MIN_VERSION: {VERSION} numpy {np.__version__} sha256 {sha256_file(os.path.abspath(__file__))} "
        f"p5min {P.VERSION} {sha256_file(P.__file__)} psp5alf {A.VERSION} {sha256_file(A.__file__)}")
    k = kolmogorov_check()
    log(f"P2MIN_K0: f^-5/3 synthetic 1/f_mid = {k:.1f} s ({k / 60:.2f} min); paper's 4-h saturation ~13 min; "
        f"f_lo = {f_lo_for(WIN_S) * 1e3:.4f} mHz")
    pos_t, pos_r = P.fetch_positions()
    peris = P.find_perihelia(pos_t, pos_r)
    allp = [(i + 1, t, r) for i, (t, r) in enumerate(peris)]
    num_ok = True
    for e, exp in P.EXPECTED_PERI.items():
        hit = [t for (nn, t, r) in allp if nn == e]
        if not hit or abs(hit[0] - P.to_unix(exp)) > 2 * 86400:
            num_ok = False
    log(f"P2MIN_PERIHELIA: n={len(allp)} numbering_ok={num_ok} "
        + " ".join(f"E{e}={P.iso(t)[:13]}" for e, t, r in allp))
    mag_stop = P.to_unix(P.MAG_STOP)
    encs = [(e, t) for e, t, r in allp if t + P.HALF_SPAN_DAYS * 86400 <= mag_stop]
    if a.only:
        want = {int(x) for x in a.only.split(",") if x.strip()}
        encs = [x for x in encs if x[0] in want]
    # E22 first (p5min's TT2000 cross-check runs on E22's first file), then the other
    # out-of-sample encounters newest first, then E19 down to E1
    encs.sort(key=lambda x: (0, 0) if x[0] == 22 else
              ((1, -x[0]) if x[0] > AUTHOR_LAST_ENC else (2, -x[0])))

    rows, done, enc_info = [], [], []
    for enc, pt in encs:
        if time.time() - t_start > a.budget_s:
            log(f"P2MIN_BUDGET: stopping before E{enc}")
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
        er = analyse_encounter(enc, pt, B, g0, spi_t, spi_vr, pos_t, pos_r)
        rows += er
        done.append(enc)
        nu = sum(r["usable"] for r in er)
        info = {"enc": enc, "perihelion": P.iso(pt), "windows": len(er),
                "mag_ok": sum(r["mag_ok"] for r in er), "v_ok": sum(r["v_ok"] for r in er),
                "usable": nu, "usable_frac": round(nu / len(er), 4) if er else 0.0,
                "clean_ok": sum(r["clean_ok"] for r in er if r["usable"]),
                "spi_samples": int(len(spi_t)), **{k2: meta[k2] for k2 in meta},
                "seconds": round(time.time() - t0, 1)}
        enc_info.append(info)
        log("P2MIN_ENC: " + json.dumps(info))
        del B

    res = summarise(rows, done)
    p0, na, nb = grade_p0(enc_info, num_ok)
    res["grades"]["P0"] = p0
    res["P0_detail"] = {"numbering_ok": num_ok, "E1_E19_ge50pct": na, "E20_E27_ge50pct": nb}
    if p0 == "FAIL":
        res["grades"]["note"] = "P0 FAIL: H1 and H2 are NOT GRADABLE (numbers indicative only)."
    if res["grades"]["R1"] != "PASS":
        res["grades"]["note_R1"] = "R1 not PASS: H1/H2 are reported 'with our reading of the paper's metric'."
    res["K0_kolmogorov_pmid_s"] = k
    res["encounters"] = enc_info
    res["runtime_s"] = round(time.time() - t_start, 1)
    with open(os.path.join(a.out_dir, "results.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=P.jdefault)
    write_rows(os.path.join(a.out_dir, "windows.csv"), rows)
    for key in ("clean_cut", "n_usable", "R1", "H1", "H2", "P0_detail"):
        log(f"P2MIN_{key.upper()}: " + json.dumps(res[key], default=P.jdefault))
    for key in ("trend_A", "trend_B", "trend_B_arith"):
        for row in res[key]:
            log(f"P2MIN_{key.upper()}: " + json.dumps(row, default=P.jdefault))
    log(f"P2MIN_PEAK: fpeak_at_floor_frac_B_young={res['fpeak_at_floor_frac_B_young']} "
        f"gmean_ppeak_B_young_s={res['gmean_ppeak_B_young_s']}")
    log("P2MIN_GRADES: " + json.dumps(res["grades"], default=P.jdefault))
    return 0


if __name__ == "__main__":
    sys.exit(main())
