#!/usr/bin/env python3
"""Non-compartmental analysis (NCA) of concentration-time data for BE studies.

Derived from K-Dense-AI/scientific-agent-skills skills/pkpd-modeling/scripts/nca.py
(MIT). Changes against that file, each checked against the guideline text
(see README.md, "NCA module"):

1. lambda_z window: among windows whose adjusted r2 is within 1e-4 of the best,
   the one with the MOST points is chosen (PKNCA pk.calc.half.life rule).
   The original kept the shorter window on near-ties.
2. A blank / NA cell is a missed sample: it is dropped and listed as a protocol
   deviation. Only explicit BLQ tokens (or values < --lloq) are set to zero
   (ICH M13A 2.2.2.2). The original turned blanks into zeros.
3. The trapezoidal rule is a required, reported choice; the default is the
   linear trapezoidal rule, the example ICH M13A 2.2.2.2 gives. The original
   called lin-up/log-down "the usual regulatory choice".
4. Added: study-level 80/20 coverage rule, AUC(0-72h), steady-state CtauSS /
   CminSS / CavSS / fluctuation / swing (ICH M13A glossary), pre-dose > 5% of
   Cmax flag, period-specific baseline correction with negatives set to zero,
   GCC emesis rule (IR: at or before 2 x median tmax), GCC Annex 1 table and a
   long table that feeds bioequivalence.py.

Scope: extravascular single-dose or steady-state profiles, time in hours.

    python nca.py -i conc.csv --format json
    python nca.py -i conc.csv --profile sfda --release ir
    python nca.py -i ss.csv --tau 12
    python nca.py -i conc.csv --truncate-72 --auc-method linup-logdown
"""
from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass, field

import numpy as np

from _common import InputError, Report, add_format_argument, main_wrapper, parse_float, read_table, require_columns

BLQ_TOKENS = {"blq", "bql", "<lloq", "lloq", "bloq", "nq", "<loq", "bllq", "nd"}
MISSING_TOKENS = {"", ".", "na", "n/a", "nan", "null", "none", "-", "missing", "ns", "nr"}
TREATMENT_ALIASES = {"TEST": "T", "REF": "R", "REFERENCE": "R"}
ADJ_R2_TOLERANCE = 1e-4          # PKNCA adj.r.squared.factor / WinNonlin
COVERAGE_MIN_PCT = 80.0          # ICH M13A 2.1.8, 2.2.2.2; GCC 3.1.4, 3.1.8
COVERAGE_MAX_SHARE = 0.20
PREDOSE_MAX_FRACTION = 0.05      # ICH M13A 2.2.3.3; GCC 3.1.8
TRUNCATION_H = 72.0              # ICH M13A 2.1.8 / 2.2.2.2; GCC 3.1.5
DEFAULT_TOLERANCE_H = 10 / 60    # ICH M13A 2.1.8: last sample within 10 min of nominal tau


@dataclass
class Sample:
    time: float
    conc: float
    blq: bool


@dataclass
class Profile:
    subject: str
    period: str = ""
    treatment: str = ""
    sequence: str = ""
    stage: str = ""
    group: str = ""
    dose: float | None = None
    emesis: float | None = None
    samples: list[Sample] = field(default_factory=list)
    missing: list[float] = field(default_factory=list)


@dataclass
class LambdaZ:
    lam: float | None = None
    intercept: float | None = None
    r2: float | None = None
    r2_adj: float | None = None
    n_points: int = 0
    t_first: float | None = None
    t_last: float | None = None
    reason: str = ""

    @property
    def half_life(self):
        return math.log(2.0) / self.lam if self.lam and self.lam > 0 else None


# ----------------------------------------------------------------- parsing


def parse_conc(raw, lloq):
    """Return (value, is_blq); value None means a missed sample."""
    text = (raw or "").strip()
    low = text.lower()
    if low in MISSING_TOKENS:
        return None, False
    if low in BLQ_TOKENS or low.startswith("<"):
        return 0.0, True
    value = parse_float(text, "concentration")
    if value < 0:
        raise InputError(f"negative concentration {text}; baseline correction is done by --baseline, not in the input")
    if value == 0 or (lloq is not None and value < lloq):
        return 0.0, True
    return value, False


def load_profiles(a):
    rows = read_table(a.input)
    require_columns(rows, [a.subject_column, a.time_column, a.conc_column], str(a.input))
    cols = rows[0].keys()
    has = {k: getattr(a, k + "_column") in cols for k in ("period", "treatment", "sequence", "stage", "group", "dose", "emesis")}
    buckets: dict[tuple, Profile] = {}
    for n, row in enumerate(rows, start=2):
        subject = row[a.subject_column].strip()
        if not subject:
            raise InputError(f"row {n}: empty subject")
        period = row[a.period_column].strip() if has["period"] else ""
        treatment = row[a.treatment_column].strip().upper() if has["treatment"] else ""
        if has["period"] and not period:
            raise InputError(f"row {n}: empty {a.period_column} for subject {subject}")
        if has["treatment"] and not treatment:
            raise InputError(f"row {n}: empty {a.treatment_column} for subject {subject}")
        key = (subject, period or treatment)
        p = buckets.get(key)
        if p is None:
            p = buckets[key] = Profile(subject, period, treatment)
        elif treatment and p.treatment != treatment:
            raise InputError(f"subject {subject} period {period}: more than one treatment in one profile")
        if has["sequence"] and row[a.sequence_column].strip():
            p.sequence = row[a.sequence_column].strip().upper()
        for k in ("stage", "group"):   # two-stage design / multi-group study: carried to the BE input
            if has[k] and row[getattr(a, k + "_column")].strip():
                setattr(p, k, row[getattr(a, k + "_column")].strip())
        if has["dose"] and row[a.dose_column].strip():
            p.dose = parse_float(row[a.dose_column], f"{a.dose_column} (row {n})")
        if has["emesis"] and row[a.emesis_column].strip():
            p.emesis = parse_float(row[a.emesis_column], f"{a.emesis_column} (row {n})")
        time_text = row[a.time_column].strip()
        if time_text.lower() in MISSING_TOKENS:
            raise InputError(f"row {n}: missing sampling time; ICH M13A needs the actual time of every sample")
        time = parse_float(time_text, f"{a.time_column} (row {n})")
        value, blq = parse_conc(row[a.conc_column], a.lloq)
        if not blq and "blq" in cols and row["blq"].strip().lower() in {"1", "y", "yes", "true"}:
            value, blq = 0.0, True
        if value is None:
            p.missing.append(time)
            continue
        if any(s.time == time for s in p.samples):
            raise InputError(f"subject {subject}{' period ' + period if period else ''}: two samples at time {time:g}")
        p.samples.append(Sample(time, value, blq))
    profiles = list(buckets.values())
    for p in profiles:
        p.samples.sort(key=lambda s: s.time)
    for k in ("stage", "group"):
        if has[k]:
            seen: dict[str, set] = {}
            for p in profiles:
                seen.setdefault(p.subject, set()).add(getattr(p, k))
            bad = sorted(s for s, v in seen.items() if len(v) != 1 or "" in v)
            if bad:
                raise InputError(f"{k}: every subject needs one {k} value in all its rows (check subject(s) {', '.join(bad[:6])})")
    if has["period"] and has["treatment"]:
        derive_sequences(profiles)
    return profiles


def derive_sequences(profiles):
    """Fill a missing sequence from the treatment order by period (e.g. TR, RTR)."""
    by_subject: dict[str, list[Profile]] = {}
    for p in profiles:
        by_subject.setdefault(p.subject, []).append(p)

    def order(p):
        try:
            return (0, float(p.period), "")
        except ValueError:
            return (1, 0.0, p.period)

    derived = {}
    for sid, group in by_subject.items():
        if all(p.sequence for p in group):
            continue
        derived[sid] = "".join(TREATMENT_ALIASES.get(p.treatment, p.treatment) for p in sorted(group, key=order))
    # A subject with missing periods gets the one full sequence of the study that matches its observed periods.
    full = max((len(v) for v in by_subject.values()), default=0)
    complete = {seq for sid, seq in derived.items() if len(by_subject[sid]) == full}
    for sid, seq in derived.items():
        group = by_subject[sid]
        if len(group) < full and all(p.period.isdigit() for p in group):
            fits = [c for c in complete if all(c[int(p.period) - 1:int(p.period)] == TREATMENT_ALIASES.get(p.treatment, p.treatment)
                                               for p in group)]
            if len(fits) == 1:
                seq = fits[0]
        for p in group:
            p.sequence = p.sequence or seq


# ------------------------------------------------------------- AUC rules


def segment(t0, t1, c0, c1, method):
    dt = t1 - t0
    if dt <= 0:
        return 0.0
    if method == "linup-logdown" and c0 > 0 and c1 > 0 and c1 < c0:
        return (c0 - c1) * dt / math.log(c0 / c1)
    return (c0 + c1) / 2.0 * dt


def interpolate(t, c, target, method):
    i = int(np.searchsorted(t, target))
    if i < len(t) and t[i] == target:
        return float(c[i])
    t0, t1, c0, c1 = t[i - 1], t[i], c[i - 1], c[i]
    if method == "linup-logdown" and c0 > 0 and c1 > 0 and c1 < c0:
        return float(c0 * math.exp(math.log(c1 / c0) * (target - t0) / (t1 - t0)))
    return float(c0 + (c1 - c0) * (target - t0) / (t1 - t0))


def auc_between(t, c, start, end, method):
    """AUC from start to end; both must lie inside the sampled range."""
    knots = [start, *[x for x in t if start < x < end], end]
    values = [interpolate(t, c, x, method) for x in knots]
    return sum(segment(knots[i - 1], knots[i], values[i - 1], values[i], method) for i in range(1, len(knots)))


# --------------------------------------------------------------- lambda_z


def fit_loglinear(t, c):
    y = np.log(c)
    tbar, ybar = t.mean(), y.mean()
    sxx = float(np.sum((t - tbar) ** 2))
    slope = float(np.sum((t - tbar) * (y - ybar)) / sxx)
    intercept = float(ybar - slope * tbar)
    ss_res = float(np.sum((y - intercept - slope * t) ** 2))
    ss_tot = float(np.sum((y - ybar) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return slope, intercept, r2


def estimate_lambda_z(t, c, quantifiable, tmax, min_points=3):
    """Best adjusted r2 over windows ending at Clast; near-ties go to more points.

    Follows PKNCA pk.calc.half.life (R/half.life.R): eligible points are
    quantifiable (BLQ omitted, ICH M13A 2.2.2.2) and after tmax; the best
    adjusted r2 is taken over all windows, and among windows with
    lambda_z > 0 and adjusted r2 > best - 1e-4 the one with most points wins.
    """
    idx = np.flatnonzero(quantifiable & (t > tmax))
    if len(idx) < min_points:
        return LambdaZ(reason=f"only {len(idx)} quantifiable points after tmax; at least {min_points} needed")
    fits = []
    for start in range(len(idx) - min_points, -1, -1):
        w = idx[start:]
        n = len(w)
        slope, intercept, r2 = fit_loglinear(t[w], c[w])
        if not math.isfinite(r2):
            continue
        r2_adj = 1.0 - (1.0 - r2) * (n - 1) / (n - 2)
        fits.append(LambdaZ(-slope, intercept, r2, r2_adj, n, float(t[w][0]), float(t[w][-1]), "best adjusted r2"))
    if not fits:
        return LambdaZ(reason="terminal concentrations are constant; no log-linear fit")
    best = max(f.r2_adj for f in fits)
    ok = [f for f in fits if f.lam > 0 and f.r2_adj > best - ADJ_R2_TOLERANCE]
    if not ok:
        return LambdaZ(reason="no window with a negative terminal slope within 1e-4 of the best adjusted r2")
    return max(ok, key=lambda f: f.n_points)


# ------------------------------------------------------------------- NCA


def label(p):
    return f"subject {p.subject}" + (f" period {p.period}" if p.period else f" {p.treatment}" if p.treatment else "")


def analyse(p, a, findings, notes_seen):
    row = {"subject": p.subject, "period": p.period, "treatment": p.treatment, "sequence": p.sequence,
           **{k: getattr(p, k) for k in ("stage", "group") if getattr(p, k)},
           "n_samples": len(p.samples), "n_blq": sum(s.blq for s in p.samples), "n_missing": len(p.missing)}
    where = label(p)
    t_all = np.array([s.time for s in p.samples], float)
    c_all = np.array([s.conc for s in p.samples], float)
    blq_all = np.array([s.blq for s in p.samples], bool)
    pre = t_all <= 0
    row["predose"] = float(c_all[pre].max()) if pre.any() else None

    def post_dose(c_src, blq_src):
        t, c, blq = t_all[t_all >= 0], c_src[t_all >= 0], blq_src[t_all >= 0]
        if len(t) and t[0] == 0:
            return t, c, blq
        if a.tau is None:
            notes_seen.add("single dose: no sample at time 0, so C(0) = 0 was inserted")
            return np.r_[0.0, t], np.r_[0.0, c], np.r_[True, blq]
        if not pre.any():
            raise InputError(f"{where}: steady-state profile needs a pre-dose or time-0 sample for C(0)")
        k = int(np.flatnonzero(pre)[-1])
        notes_seen.add("steady state: the last pre-dose sample was used as C(0)")
        return np.r_[0.0, t], np.r_[c_src[k], c], np.r_[blq_src[k], blq]

    if a.baseline == "predose-mean":
        if not pre.any():
            raise InputError(f"{where}: baseline correction needs at least one pre-dose sample (time <= 0)")
        base = float(c_all[pre].mean())
        row["baseline"] = base
        tu, cu, bu = post_dose(c_all, blq_all)
        qu = np.flatnonzero((~bu) & (cu > 0))
        if len(qu):
            row["cmax_uncorrected"] = float(cu[qu].max())
            row["auc_0_t_uncorrected"] = auc_between(tu, cu, 0.0, float(tu[qu[-1]]), a.auc_method)
        c_all = np.where(pre, 0.0, np.maximum(c_all - base, 0.0))   # ICH M13A 3.1: negative -> zero
        blq_all = blq_all | (c_all <= 0)

    t, c, blq = post_dose(c_all, blq_all)
    quant = (~blq) & (c > 0)
    if not quant.any():
        row["status"] = "excluded: no quantifiable concentration"
        findings.append(f"{where}: no quantifiable concentration; no parameters (GCC 3.1.8 exception 1 if this is the reference product)")
        return row, LambdaZ(reason="all BLQ")

    peak = int(np.argmax(np.where(quant, c, -np.inf)))
    cmax, tmax = float(c[peak]), float(t[peak])
    last = int(np.flatnonzero(quant)[-1])
    clast, tlast = float(c[last]), float(t[last])
    row.update(cmax=cmax, tmax=tmax, clast=clast, tlast=tlast)
    if p.missing:
        findings.append(f"{where}: missed sample(s) at {', '.join(f'{x:g}' for x in p.missing)} h dropped, not set to zero; report as protocol deviation")
    mid_blq = np.flatnonzero(blq[:last] & (t[:last] > 0) & (np.arange(last) > np.flatnonzero(quant)[0]))
    if len(mid_blq):
        findings.append(f"{where}: BLQ between quantifiable samples at {', '.join(f'{t[i]:g}' for i in mid_blq)} h set to zero (ICH M13A); check the bioanalytical record")
    first_post = np.flatnonzero(t > 0)
    if a.tau is None and len(first_post) and peak == int(first_post[0]):
        findings.append(f"{where}: Cmax is the first post-dose sample; the true peak may precede it (GCC 3.1.4)")
    if last == peak:
        findings.append(f"{where}: Cmax is the last quantifiable sample; the terminal phase is not characterised")

    auc = [0.0]
    for i in range(1, len(t)):
        auc.append(auc[-1] + segment(t[i - 1], t[i], c[i - 1], c[i], a.auc_method))
    row["auc_0_t"] = auc_last = float(auc[last])

    lz = estimate_lambda_z(t, c, quant, tmax, a.lambda_z_points)
    if lz.lam:
        row.update(lambda_z=lz.lam, t_half=lz.half_life, lz_n=lz.n_points, lz_start=lz.t_first, lz_end=lz.t_last,
                   r2_adj=lz.r2_adj)
    if lz.lam and a.tau is None:     # AUC(0-inf) and coverage are single-dose parameters
        auc_inf = auc_last + clast / lz.lam
        row["auc_0_inf"] = auc_inf
        row["auc_0_inf_pred"] = auc_last + math.exp(lz.intercept - lz.lam * tlast) / lz.lam
        row["auc_pct_extrap"] = 100.0 * (auc_inf - auc_last) / auc_inf
        row["auc_coverage_pct"] = 100.0 * auc_last / auc_inf
        if p.dose:
            row["cl_f"] = p.dose / auc_inf
            row["vz_f"] = p.dose / auc_inf / lz.lam
    elif not lz.lam and a.tau is None and not a.truncate_72:
        findings.append(f"{where}: kel not estimable ({lz.reason}); AUC(0-inf) and coverage not available")

    if a.truncate_72:
        end = float(t[-1])
        row["c72_quantifiable"] = tlast >= TRUNCATION_H - a.tolerance
        if a.profile == "sfda" and end >= TRUNCATION_H - a.tolerance and not row["c72_quantifiable"]:
            findings.append(f"{where}: concentration at 72 h is not quantifiable (last quantifiable sample at {tlast:g} h); "
                            "GCC 3.1.5 then requires AUC(0-inf) and the residual area")
        if end >= TRUNCATION_H:
            row["auc_0_72"] = auc_between(t, c, 0.0, TRUNCATION_H, a.auc_method)
        elif end >= TRUNCATION_H - a.tolerance:
            row["auc_0_72"] = float(auc[-1])
            notes_seen.add(f"AUC(0-72h) taken to the last sample when it lies within {a.tolerance * 60:g} min before 72 h")
        else:
            findings.append(f"{where}: last sample at {end:g} h, before 72 h; AUC(0-72h) not computed")

    if a.tau is not None:
        ss_tau(row, t, c, quant, a, where, findings, notes_seen)

    for start, end in a.partial_auc:
        key = f"auc_{start:g}_{end:g}"
        row[key] = auc_between(t, c, start, end, a.auc_method) if t[-1] >= end else None
        if row[key] is None:
            findings.append(f"{where}: partial AUC {start:g}-{end:g} h extends past the last sample")
    return row, lz


def ss_tau(row, t, c, quant, a, where, findings, notes_seen):
    tau, tol = a.tau, a.tolerance
    near = np.flatnonzero(np.abs(t - tau) <= tol)
    if len(near):
        k = int(near[np.argmin(np.abs(t[near] - tau))])
        row["ctau_ss"], row["ctau_time"] = float(c[k]), float(t[k])
    else:
        findings.append(f"{where}: no sample within {tol * 60:g} min of tau = {tau:g} h; CtauSS not reported (program rule; ICH M13A 2.1.8 recommends the last sample within 10 min of the nominal tau)")
    if t[-1] >= tau:
        auc_tau = auc_between(t, c, 0.0, tau, a.auc_method)
    elif t[-1] >= tau - tol:
        auc_tau = auc_between(t, c, 0.0, float(t[-1]), a.auc_method)
        notes_seen.add(f"AUC(0-tau) taken to the last sample when it lies within {tol * 60:g} min before tau")
    else:
        findings.append(f"{where}: last sample at {t[-1]:g} h, before tau = {tau:g} h; AUC(0-tau) not computed")
        return
    window = (t >= 0) & (t <= tau + tol)
    cw, tw = c[window], t[window]
    k = int(np.argmax(cw))
    cmin = float(cw.min())
    cav = auc_tau / tau
    row.update(auc_tau=auc_tau, cmax_ss=float(cw[k]), tmax_ss=float(tw[k]), cmin_ss=cmin, cav_ss=cav,
               fluctuation_pct=100.0 * (float(cw[k]) - cmin) / cav if cav > 0 else None,
               swing=(float(cw[k]) - cmin) / cmin if cmin > 0 else None)
    if cmin <= 0:
        findings.append(f"{where}: CminSS is zero or BLQ; swing is not defined")


# --------------------------------------------------------------- rules


def apply_emesis(rows, profiles, a, findings, notes):
    flagged = [(r, p) for r, p in zip(rows, profiles) if p.emesis is not None]
    if not flagged:
        return
    if a.profile != "sfda":
        for r, p in flagged:
            findings.append(f"{label(p)}: emesis at {p.emesis:g} h; ICH M13A Q&A 2.9 lists emesis within 2 x the expected median tmax as a deviation that may justify exclusion if pre-specified in the protocol (the GCC rule is applied with --profile sfda)")
        return
    if a.release == "mr":
        for r, p in flagged:
            if a.dosing_interval is None or p.emesis <= a.dosing_interval:
                r["status"] = f"excluded: emesis at {p.emesis:g} h (modified release)"
        notes.append("GCC 3.1.8: modified-release subjects with emesis during the labelled dosing interval are excluded"
                     + ("" if a.dosing_interval else "; no --dosing-interval given, so any recorded emesis excludes"))
        return
    observed_used = False
    for r, p in flagged:
        med = a.median_tmax.get(p.treatment, a.median_tmax.get("*")) if a.median_tmax else None
        if med is None:
            others = [x.get("tmax") for x, q in zip(rows, profiles) if q.emesis is None and q.treatment == p.treatment and x.get("tmax") is not None]
            if not others:
                findings.append(f"{label(p)}: emesis recorded but no median tmax available for {p.treatment or 'this product'}; give --median-tmax")
                continue
            med, observed_used = float(np.median(others)), True
        r["emesis_limit"] = 2 * med
        if p.emesis <= 2 * med:
            r["status"] = f"excluded: emesis at {p.emesis:g} h <= 2 x median tmax ({2 * med:g} h)"
    notes.append("GCC 3.1.8 (immediate release): a period is excluded if emesis occurs at or before 2 x median tmax")
    if observed_used:
        findings.append("emesis rule used the median tmax observed in this study (per product, periods without emesis): "
                        "exclusion decisions should be made before bioanalysis (GCC 3.1.8) and ICH M13A Q&A 2.9 refers to the "
                        "expected median tmax; give the protocol value with --median-tmax")


def coverage_rule(rows, a, report):
    if a.tau is not None or a.truncate_72:
        report.note("80/20 coverage rule not applied: " + ("steady-state study" if a.tau is not None else
                    "AUC(0-72h) is the primary AUC (ICH M13A 2.2.2.2, GCC 3.1.8)"))
        return
    included = [r for r in rows if r.get("status") == "included"]
    cov = [r["auc_coverage_pct"] for r in included if r.get("auc_coverage_pct") is not None]
    below = sum(x < COVERAGE_MIN_PCT for x in cov)
    report.scalar("coverage_evaluable", len(cov))
    report.scalar("coverage_not_evaluable", len(included) - len(cov))
    report.scalar("coverage_below_80", below)
    share = below / len(cov) if cov else None
    report.scalar("coverage_below_80_pct", 100 * share if share is not None else None)
    if share is not None and share > COVERAGE_MAX_SHARE:
        report.finding(f"AUC(0-t) covers less than 80% of AUC(0-inf) in {below} of {len(cov)} observations "
                       f"({100 * share:.1f}% > 20%); the validity of the study may need to be discussed (ICH M13A 2.2.2.2, GCC 3.1.8)")


SUMMARY_KEYS = ["cmax", "tmax", "auc_0_t", "auc_0_inf", "auc_0_72", "auc_coverage_pct", "lambda_z", "t_half",
                "auc_tau", "cmax_ss", "tmax_ss", "cmin_ss", "ctau_ss", "cav_ss", "fluctuation_pct", "swing", "cl_f", "vz_f"]


def describe(values):
    arr = np.asarray(values, float)
    n = len(arr)
    out = {"n": n, "mean": float(arr.mean()), "sd": float(arr.std(ddof=1)) if n > 1 else None,
           "median": float(np.median(arr)), "min": float(arr.min()), "max": float(arr.max())}
    out["cv_pct"] = 100 * out["sd"] / out["mean"] if out["sd"] is not None and out["mean"] else None
    if np.all(arr > 0):
        logs = np.log(arr)
        out["geo_mean"] = float(math.exp(logs.mean()))
        out["geo_cv_pct"] = 100 * math.sqrt(math.expm1(float(logs.var(ddof=1)))) if n > 1 else None
    return out


def summary(rows, keys, left_out=()):
    out = []
    groups = sorted({r["treatment"] for r in rows})
    for g in groups:
        for key in keys:
            vals = [r[key] for r in rows if r["treatment"] == g and r.get("status") == "included" and r["subject"] not in left_out
                    and isinstance(r.get(key), (int, float))]
            if vals:
                out.append({"treatment": g or "-", "parameter": key, **describe(vals)})
    return out


def gcc_annex(stats, test, ref):
    lookup = {(s["treatment"], s["parameter"]): s for s in stats}
    rows = []
    for name, key in (("AUC(0-t)", "auc_0_t"), ("AUC(0-inf)", "auc_0_inf"), ("AUC(0-72h)", "auc_0_72"),
                      ("Cmax", "cmax"), ("Tmax", "tmax"), ("T1/2", "t_half")):
        tt, rr = lookup.get((test, key)), lookup.get((ref, key))
        if tt and rr:
            rows.append({"parameter": name, "test_mean": tt["mean"], "test_cv_pct": tt["cv_pct"],
                         "reference_mean": rr["mean"], "reference_cv_pct": rr["cv_pct"]})
    return rows


def clean(v):
    return None if isinstance(v, float) and not math.isfinite(v) else v


def population(rows):
    """Who is in the BE analysis and in the summary statistics.

    Returns (design info, {subject: reason left out of the BE input}, {subject: reason left out of the
    summary statistics}), or (None, {}, {}).
    2x2: a subject with an excluded or missing period is left out of both (GCC 3.1.8: in a 2-period trial
    the subject is removed). Replicate designs: only the excluded periods leave the BE input, and subjects
    with missing periods stay in (bioequivalence.py analyses them with the fixed-effects model, EMA Method A;
    ICH M13A Q&A 2.9); drop-outs are still left out of the summary statistics (GCC 3.1.8).
    Designs with more than two treatments keep the other periods.
    """
    if not rows or not all(r["treatment"] for r in rows):
        return None, {}, {}
    by_subject: dict[str, list[dict]] = {}
    for r in rows:
        by_subject.setdefault(r["subject"], []).append(r)
    parallel = all(len(v) == 1 for v in by_subject.values())
    if not parallel and not all(r["period"] for r in rows):
        return None, {}, {}
    labels = {TREATMENT_ALIASES.get(r["treatment"], r["treatment"]) for r in rows}
    info = {"parallel": parallel, "multi": len(labels) > 2, "periods": max(len(v) for v in by_subject.values())}
    left_out: dict[str, str] = {}
    summary_out: dict[str, str] = {}
    if not parallel and not info["multi"]:
        for subject, group in by_subject.items():
            short = len(group) < info["periods"]
            if info["periods"] == 2 and any(r["status"] != "included" for r in group):
                left_out[subject] = summary_out[subject] = "a period was excluded"
            elif short:
                summary_out[subject] = f"{len(group)} of {info['periods']} periods in the data"
                if info["periods"] == 2:
                    left_out[subject] = summary_out[subject]
    return info, left_out, summary_out


def be_labels(rows, a):
    """Treatment codes for bioequivalence.py: Test/Ref aliases -> T/R, and with two products the
    --test-label / --reference-label codes -> T/R. Designs with more than two products keep their codes."""
    labels = sorted({r["treatment"] for r in rows if r["treatment"]})
    out = {x: TREATMENT_ALIASES.get(x, x) for x in labels}
    t_lab, r_lab = a.test_label.strip().upper(), a.reference_label.strip().upper()
    if len(labels) == 2 and t_lab != r_lab and {t_lab, r_lab} <= set(labels):
        out[t_lab], out[r_lab] = "T", "R"
    return out


def _natural(text):
    return [(0, int(x), "") if x.isdigit() else (1, 0, x) for x in re.split(r"(\d+)", text) if x]


def be_periods(rows):
    """Integer periods for bioequivalence.py: 1, 2.0 -> 1, 2; other codes (P1, P2, ...) numbered in natural order."""
    labels = sorted({r["period"] for r in rows if r["period"]}, key=_natural)

    def as_int(x):
        try:
            f = float(x)
        except ValueError:
            return None
        return int(f) if math.isfinite(f) and f.is_integer() and f >= 1 else None

    ints = {x: as_int(x) for x in labels}
    if all(v is not None for v in ints.values()) and len(set(ints.values())) == len(ints):
        return {x: str(v) for x, v in ints.items()}, any(x != str(v) for x, v in ints.items())
    return {x: str(i + 1) for i, x in enumerate(labels)}, True


def be_sequences(rows, codes, periods, full):
    """Sequence of every subject in T/R codes, in period order. A sequence column already in these codes
    is kept; otherwise (e.g. 1/2, AB/BA, or no column) it is derived from the treatments. A subject with
    missing periods gets the one full sequence of the study that fits its observed periods.
    Returns ({subject: sequence}, {subject: reason it cannot be determined})."""
    by_subject: dict[str, list[dict]] = {}
    for r in rows:
        by_subject.setdefault(r["subject"], []).append(r)
    letters = set(codes.values())
    pos = lambda r: int(periods[r["period"]])
    seq, unknown, derived = {}, {}, {}
    for s, g in by_subject.items():
        given = {r["sequence"] for r in g}
        given = given.pop() if len(given) == 1 else ""
        if given and len(given) == full and set(given) <= letters:
            seq[s] = given
            continue
        derived[s] = "".join(codes[r["treatment"]] for r in sorted(g, key=pos))
    complete = ({v for v in derived.values() if len(v) == full} |
                {v for s, v in seq.items() if len(by_subject[s]) == full})
    for s, d in derived.items():
        if len(d) == full:
            seq[s] = d
            continue
        fits = [c for c in complete if all(c[pos(r) - 1:pos(r)] == codes[r["treatment"]] for r in by_subject[s])]
        if len(fits) == 1:
            seq[s] = fits[0]
        else:
            unknown[s] = "sequence cannot be determined from the observed periods; add a sequence column"
    return seq, unknown


def be_input(rows, a, report, notes):
    """Long table for bioequivalence.py: one row per subject-period (or per subject in a parallel study).

    2x2: a subject with an excluded or missing period is left out (GCC 3.1.8: in a 2-period trial the
    subject is removed). Replicate designs keep the other periods of such subjects; bioequivalence.py
    analyses them with the fixed-effects model (EMA Method A; ICH M13A Q&A 2.9). Designs with more than
    two treatments keep the other periods: each comparison uses the subjects that have both products.
    Treatment codes become T/R, periods integers and sequences T/R strings; stage and group columns
    (two-stage design, multi-group study) are carried over.
    """
    info, left_out, _ = population(rows)
    if info is None:
        return
    parallel, multi = info["parallel"], info["multi"]
    primary = ("auc_tau", "cmax_ss") if a.tau is not None else ("cmax", "auc_0_t", "auc_0_72", "auc_0_inf")
    keys = [k for k in primary if any(k in r for r in rows)] + [f"auc_{s:g}_{e:g}" for s, e in a.partial_auc]
    codes = be_labels(rows, a)
    periods, renumbered = be_periods(rows) if not parallel else ({}, False)
    left_out = dict(left_out)
    seqs = {}
    two = set(codes.values()) == {"T", "R"}
    if not parallel and not multi and two:
        seqs, unknown = be_sequences(rows, codes, periods, info["periods"])
        for s, why in unknown.items():
            left_out.setdefault(s, why)
    extra = [k for k in ("stage", "group") if any(r.get(k) for r in rows)]

    table = []
    for r in rows:
        if r["subject"] in left_out or r["status"] != "included":
            continue
        trt = codes[r["treatment"]]
        seq = trt if parallel else seqs.get(r["subject"], r["sequence"])
        table.append({"subject": r["subject"], "sequence": seq, "period": periods.get(r["period"], r["period"] or "1"),
                      "treatment": trt, **{k: r.get(k, "") for k in extra},
                      **{k: clean(r.get(k)) for k in keys}, "predose": clean(r.get("predose"))})
    report.table("be input", table)
    mapped = {k: v for k, v in codes.items() if k != v}
    if mapped:
        notes.append("BE input: treatment codes " + ", ".join(f"{k} -> {v}" for k, v in sorted(mapped.items())))
    if not multi and not two:
        notes.append("BE input: the treatment codes are not T and R; choose the test and reference labels "
                     "(--test-label / --reference-label) before sending the parameters to the BE analysis")
    if renumbered:
        notes.append("BE input: periods numbered in order: " +
                     ", ".join(f"{k} -> {v}" for k, v in sorted(periods.items(), key=lambda kv: int(kv[1]))))
    if "stage" in extra:
        notes.append("BE input keeps the stage column (two-stage design): the BE analysis fits periods within stage "
                     "with the fixed-effects model; set the pre-specified confidence level (e.g. 94.12%)")
    if "group" in extra:
        notes.append("BE input keeps the group column (multi-group study, ICH M13A 2.2.3.5): the BE analysis adds the group terms")
    if left_out:
        notes.append("BE input leaves out " + "; ".join(f"subject {s}: {why}" for s, why in sorted(left_out.items())))
    if multi and any(r["status"] != "included" for r in rows):
        notes.append("BE input keeps the other periods of subjects with an excluded period: with more than two products, "
                     "each comparison uses only the subjects that have both products")
    if info["periods"] > 2 and not multi and not parallel and (
            any(r["status"] != "included" for r in rows) or
            any(len([x for x in rows if x["subject"] == s]) < info["periods"] for s in {r["subject"] for r in rows})):
        notes.append("BE input (replicate design) keeps the available periods of subjects with an excluded or missing period; "
                     "analyse it with the fixed-effects model (EMA Method A), which uses all available observations")
    periodwise = not parallel and (info["periods"] > 2 or multi)
    for k in keys:
        gaps = sorted({x["subject"] for x in table if x.get(k) is None})
        if gaps:
            notes.append(f"BE input: {k} is not available for subject(s) {', '.join(gaps)}; "
                         + ("only those periods are left out when it is analysed" if periodwise
                            else "they are left out when it is analysed"))


# --------------------------------------------------------------------- CLI


def _range(text):
    try:
        lo, hi = (float(x) for x in text.split("-", 1))
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"expected a range like 0-2, got {text!r}") from exc
    if hi <= lo:
        raise argparse.ArgumentTypeError(f"range end must exceed start, got {text!r}")
    return lo, hi


def _median_tmax(text):
    out = {}
    for part in text.split(","):
        key, _, val = part.rpartition("=")
        try:
            x = float(val)
        except ValueError as exc:
            raise argparse.ArgumentTypeError(f"expected hours, e.g. 1.5 or T=1.5,R=2, got {text!r}") from exc
        if not math.isfinite(x) or x <= 0:
            raise argparse.ArgumentTypeError(f"median tmax must be a positive number of hours, got {text!r}")
        out[key.strip().upper() or "*"] = x
    return out


def build_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", required=True)
    for name, default in (("subject", "subject"), ("time", "time"), ("conc", "conc"), ("period", "period"),
                          ("treatment", "treatment"), ("sequence", "sequence"), ("stage", "stage"), ("group", "group"),
                          ("dose", "dose"), ("emesis", "emesis_time")):
        p.add_argument(f"--{name}-column", default=default)
    p.add_argument("--profile", choices=("none", "sfda"), default="none", help="sfda: GCC DS-G-010 V3.1 rules on top of ICH M13A")
    p.add_argument("--auc-method", choices=("linear", "linup-logdown"), default="linear",
                   help="pre-specify in the protocol; default linear trapezoidal (ICH M13A 2.2.2.2 example)")
    p.add_argument("--lloq", type=float, help="values below it are BLQ (set to zero, omitted from kel)")
    p.add_argument("--tau", type=float, help="dosing interval (h): steady-state parameters")
    p.add_argument("--tolerance", type=float, default=DEFAULT_TOLERANCE_H, help="h; sample-time tolerance for tau and 72 h (default 10 min)")
    p.add_argument("--truncate-72", action="store_true", help="AUC(0-72h) is the primary AUC (long half-life)")
    p.add_argument("--baseline", choices=("none", "predose-mean"), default="none",
                   help="endogenous: subtract the mean pre-dose concentration of the same period")
    p.add_argument("--release", choices=("ir", "mr"), default="ir", help="for the GCC emesis rule")
    p.add_argument("--dosing-interval", type=float, help="labelled dosing interval (h), MR emesis rule")
    p.add_argument("--median-tmax", type=_median_tmax, default={},
                   help="expected median tmax (h) from the protocol for the emesis rule: one value, or per product e.g. T=1.5,R=2")
    p.add_argument("--partial-auc", type=_range, action="append", default=[], help="e.g. 0-2 (early exposure)")
    p.add_argument("--lambda-z-points", type=int, default=3, help="minimum points in the terminal fit (ICH M13A: three or more)")
    p.add_argument("--test-label", default="T")
    p.add_argument("--reference-label", default="R")
    add_format_argument(p)
    return p


def run(argv=None):
    a = build_parser().parse_args(argv)
    for name in ("subject", "time", "conc", "period", "treatment", "sequence", "stage", "group", "dose", "emesis"):
        setattr(a, name + "_column", getattr(a, name + "_column").strip().lower())   # read_table lower-cases headers
    for name in ("tau", "lloq", "tolerance", "dosing_interval"):
        value = getattr(a, name)
        if value is not None and not math.isfinite(value):
            raise InputError(f"--{name.replace('_', '-')} must be a finite number")
    if a.lambda_z_points < 3:
        raise InputError("--lambda-z-points must be at least 3 (ICH M13A 2.1.8)")
    if a.tau is not None and a.tau <= 0:
        raise InputError("--tau must be positive")
    if a.tau is not None and a.truncate_72:
        raise InputError("--truncate-72 applies to single-dose studies, not with --tau")
    profiles = load_profiles(a)
    report = Report()
    findings: list[str] = []
    notes_seen: set[str] = set()
    rows, diag = [], []
    for p in profiles:
        row, lz = analyse(p, a, findings, notes_seen)
        row.setdefault("status", "included")
        rows.append(row)
        diag.append({"subject": p.subject, "period": p.period, "treatment": p.treatment, "lambda_z": lz.lam,
                     "t_half": lz.half_life, "n_points": lz.n_points, "start": lz.t_first, "end": lz.t_last,
                     "r2": lz.r2, "r2_adj": lz.r2_adj, "basis": lz.reason})

    notes: list[str] = []
    apply_emesis(rows, profiles, a, findings, notes)
    single = a.tau is None and a.baseline == "none"
    for r, p in zip(rows, profiles):
        if single and r.get("predose") and r.get("cmax") and r["predose"] > PREDOSE_MAX_FRACTION * r["cmax"]:
            r["predose_pct_cmax"] = 100 * r["predose"] / r["cmax"]
            findings.append(f"{label(p)}: pre-dose concentration is {r['predose_pct_cmax']:.1f}% of Cmax (> 5%); "
                            "exclude this period from the primary analysis (ICH M13A 2.2.3.3, GCC 3.1.8)")

    report.scalar("profile", a.profile)
    report.scalar("auc_method", a.auc_method)
    report.scalar("blq_rule", "zero; omitted from kel")
    report.scalar("n_profiles", len(rows))
    report.scalar("n_included", sum(r["status"] == "included" for r in rows))
    if a.tau is not None:
        report.scalar("tau", a.tau)
    coverage_rule(rows, a, report)

    keys = [k for k in SUMMARY_KEYS if any(k in r for r in rows)] + [f"auc_{s:g}_{e:g}" for s, e in a.partial_auc]
    _, _, summary_out = population(rows)
    stats = summary(rows, keys, summary_out)
    if summary_out:
        notes.append("summary statistics leave out subject(s) " + ", ".join(sorted(summary_out)) +
                     ": drop-outs and excluded subjects are listed individually but not summarised (GCC 3.1.8)")
    report.table("per-profile parameters", [{k: clean(v) for k, v in r.items()} for r in rows])
    report.table("lambda_z diagnostics", [{k: clean(v) for k, v in d.items()} for d in diag])
    report.table("summary statistics", [{k: clean(v) for k, v in s.items()} for s in stats])
    missing = [{"subject": p.subject, "period": p.period, "treatment": p.treatment, "time": x} for p in profiles for x in p.missing]
    if missing:
        report.table("missed samples", missing)
    excluded = [{"subject": r["subject"], "period": r["period"], "treatment": r["treatment"], "reason": r["status"][10:]}
                for r in rows if r["status"] != "included"]
    if excluded:
        report.table("exclusions", excluded)
    if a.profile == "sfda" and a.tau is None:
        codes = be_labels(rows, a) if all(r["treatment"] for r in rows) else {}
        raw = {v: k for k, v in codes.items()}
        annex = gcc_annex(stats, raw.get("T"), raw.get("R")) if {"T", "R"} <= set(raw) else []
        if annex:
            report.table("gcc annex 1", [{k: clean(v) for k, v in x.items()} for x in annex])
            notes.append("GCC Annex 1 section 5: the ratio of geometric means and the confidence interval come from the BE analysis")
        else:
            notes.append("GCC Annex 1 table not produced: it needs a test and a reference product (treatment codes T/R, "
                         "Test/Reference, or --test-label / --reference-label)")

    be_input(rows, a, report, notes)

    report.note(f"trapezoidal rule: {'linear' if a.auc_method == 'linear' else 'linear-up / log-down'}; pre-specify it in the protocol and report it (ICH M13A 2.2.2.2)")
    report.note("BLQ values are zero in AUC and omitted from kel (ICH M13A 2.2.2.2); blank cells are missed samples, dropped")
    report.note(f"kel: best adjusted r2 over windows of >= {a.lambda_z_points} points after tmax ending at Clast; "
                "within 1e-4 of the best, the window with more points is used. Report the number of points (ICH M13A, GCC 3.1.8)")
    report.note("parameters use the actual sampling times in the input; times are taken as hours after dosing")
    if a.baseline == "predose-mean":
        report.note("baseline: mean pre-dose concentration of each period subtracted, negatives set to zero (ICH M13A 3.1; GCC 3.1.5 prefers subtracting the mean pre-dose concentration). "
                    "ICH M13A 3.1 asks for PK and statistical analyses of the uncorrected data as well: run the analysis again without --baseline")
    for n in sorted(notes_seen):
        report.note(n)
    for n in notes:
        report.note(n)
    for f in findings:
        report.finding(f)
    return report.emit(a.format)


if __name__ == "__main__":
    raise SystemExit(main_wrapper(run))
