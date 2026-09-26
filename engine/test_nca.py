"""Tests for nca.py: closed-form profiles, an independent reference
implementation (scipy.stats.linregress + scipy.integrate), the four fixes
against the upstream nca.py, and the ICH M13A / GCC rules."""
import csv
import io
import json
import math
import os
import tempfile
import unittest
from contextlib import redirect_stdout, redirect_stderr

import numpy as np
from scipy.integrate import trapezoid
from scipy.stats import linregress

import bioequivalence as be
import nca

D, V, KA, K = 100.0, 10.0, 1.2, 0.1
TIMES = [0, .25, .5, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24, 36, 48, 72]


def conc(t):
    return D * KA / (V * (KA - K)) * (math.exp(-K * t) - math.exp(-KA * t))


def auc_exact(t):
    """Closed-form AUC(0-t) of the one-compartment oral profile."""
    a = D * KA / (V * (KA - K))
    return a * ((1 - math.exp(-K * t)) / K - (1 - math.exp(-KA * t)) / KA)


def run(rows, *argv, cols=("subject", "time", "conc")):
    fd, path = tempfile.mkstemp(suffix=".csv")
    os.close(fd)
    try:
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(cols)
            w.writerows(rows)
        out, err = io.StringIO(), io.StringIO()
        with redirect_stdout(out), redirect_stderr(err):
            code = nca.main_wrapper(nca.run, ["-i", path, *argv, "--format", "json"])
        return code, (json.loads(out.getvalue()) if out.getvalue() else None), err.getvalue()
    finally:
        os.remove(path)


def table(payload, title):
    return next(t["rows"] for t in payload["tables"] if t["title"] == title)


def first(payload):
    return table(payload, "per-profile parameters")[0]


# ------------------------------------------------ independent reference


def ref_lambda_z(t, c, tmax):
    """Brute force over every terminal window with scipy.stats.linregress; PKNCA selection
    (R/half.life.R): best adj r2 over all windows, then lambda_z > 0 and adj r2 > best - 1e-4, most points."""
    pts = [(ti, ci) for ti, ci in zip(t, c) if ti > tmax and ci > 0]
    fits = []
    for s in range(len(pts) - 2):
        w = pts[s:]
        fit = linregress([p[0] for p in w], [math.log(p[1]) for p in w])
        n = len(w)
        adj = 1 - (1 - fit.rvalue ** 2) * (n - 1) / (n - 2)
        fits.append((adj, n, -fit.slope))
    if not fits:
        return None
    top = max(b[0] for b in fits)
    ok = [b for b in fits if b[2] > 0 and b[0] > top - 1e-4]
    return max(ok, key=lambda b: b[1]) if ok else None


def ref_nca(t, c, blq, method):
    t, c = np.asarray(t, float), np.where(blq, 0.0, c)
    q = np.flatnonzero(~np.asarray(blq) & (c > 0))
    last, peak = q[-1], q[np.argmax(c[q])]
    tt, cc = t[: last + 1], c[: last + 1]
    if method == "linear":
        auc = trapezoid(cc, tt)
    else:
        auc = 0.0
        for i in range(1, len(tt)):
            c0, c1, dt = cc[i - 1], cc[i], tt[i] - tt[i - 1]
            auc += dt * (c0 - c1) / (math.log(c0) - math.log(c1)) if 0 < c1 < c0 else dt * (c0 + c1) / 2
    lz = ref_lambda_z([x for x, b in zip(t, blq) if not b], [x for x, b in zip(c, blq) if not b], t[peak])
    return {"cmax": c[peak], "tmax": t[peak], "auc_0_t": auc, "lambda_z": lz[2] if lz else None,
            "lz_n": lz[1] if lz else None, "auc_0_inf": auc + c[last] / lz[2] if lz else None}


class ClosedForm(unittest.TestCase):
    def test_analytical_profile(self):
        code, p, _ = run([["S1", x, conc(x)] for x in TIMES], "--auc-method", "linup-logdown")
        r = first(p)
        # the tie rule takes 3-72 h (adj r2 0.99997, within 1e-4 of the 3-point fit), where the absorption
        # term still contributes a little: lambda_z is 0.2% low, exactly as PKNCA/WinNonlin would report
        self.assertEqual(r["lz_n"], 10)
        self.assertLess(abs(r["lambda_z"] / K - 1), 0.005)
        self.assertLess(abs(r["auc_0_inf"] / (D / (V * K)) - 1), 0.01)
        self.assertLess(abs(r["auc_0_t"] / auc_exact(72) - 1), 0.01)
        self.assertEqual(r["tmax"], 2)

    def test_linear_rule_is_default_and_reported(self):
        code, p, _ = run([["S1", x, conc(x)] for x in TIMES])
        self.assertEqual(p["scalars"]["auc_method"], "linear")
        self.assertTrue(any("trapezoidal rule: linear" in n for n in p["notes"]))
        # linear trapezoids over-estimate a convex decline
        self.assertGreater(first(p)["auc_0_t"], auc_exact(72))

    def test_auc_0_72_interpolated(self):
        times = TIMES[:-1] + [96]
        for method in ("linear", "linup-logdown"):
            code, p, _ = run([["S1", x, conc(x)] for x in times], "--truncate-72", "--auc-method", method)
            r = first(p)
            c48, c96 = conc(48), conc(96)
            c72 = c48 + (c96 - c48) * 0.5 if method == "linear" else c48 * (c96 / c48) ** 0.5
            tail = (c48 + c72) / 2 * 24 if method == "linear" else (c48 - c72) * 24 / math.log(c48 / c72)
            code2, p2, _ = run([["S1", x, conc(x)] for x in TIMES[:-1]], "--auc-method", method)
            self.assertAlmostEqual(r["auc_0_72"], first(p2)["auc_0_t"] + tail, places=9)
            self.assertTrue(any("not applied" in n for n in p["notes"]))


class Reference(unittest.TestCase):
    def test_random_profiles_match_independent_implementation(self):
        rng = np.random.default_rng(7)
        checked = with_lz = 0
        for trial in range(300):
            ka, k = rng.uniform(.3, 3), rng.uniform(.02, .4)
            if abs(ka - k) < .05:
                continue
            ts = sorted(set(np.round(np.r_[0, rng.uniform(.1, 72, 13)], 3)))
            cs = [0.0 if x == 0 else 50 * ka / (ka - k) * (math.exp(-k * x) - math.exp(-ka * x)) * math.exp(rng.normal(0, .15)) for x in ts]
            lloq = sorted(cs)[len(cs) // 4]
            blq = [x == 0 or ci < lloq for x, ci in zip(ts, cs)]
            if sum(not b for b in blq) < 5:
                continue
            method = "linear" if trial % 2 else "linup-logdown"
            rows = [["S", x, "BLQ" if b else repr(ci)] for x, ci, b in zip(ts, cs, blq)]
            code, p, err = run(rows, "--auc-method", method)
            r, ref = first(p), ref_nca(ts, cs, blq, method)
            checked += 1
            with_lz += ref["lambda_z"] is not None
            for key in ("cmax", "tmax", "auc_0_t"):
                self.assertAlmostEqual(r[key], ref[key], places=9, msg=(trial, key))
            if ref["lambda_z"] is None:
                self.assertNotIn("lambda_z", r)
            else:
                self.assertEqual(r["lz_n"], ref["lz_n"], msg=trial)
                self.assertAlmostEqual(r["lambda_z"], ref["lambda_z"], places=10, msg=trial)
                self.assertAlmostEqual(r["auc_0_inf"], ref["auc_0_inf"], places=8, msg=trial)
        print(f" [{checked} random profiles, {with_lz} with lambda_z]", end=" ")
        self.assertGreater(with_lz, 200)


class UpstreamFixes(unittest.TestCase):
    def test_lambda_z_near_tie_prefers_more_points(self):
        tail_t = [4, 6, 8, 10]
        tail = [80 * math.exp(-.2 * (x - 2)) for x in tail_t]
        tail[0] *= 1.0025
        code, p, _ = run([["S", 0, 0], ["S", 1, 50], ["S", 2, 80]] + [["S", x, c] for x, c in zip(tail_t, tail)])
        diag = table(p, "lambda_z diagnostics")[0]
        self.assertEqual(diag["n_points"], 4)
        self.assertEqual(ref_lambda_z(tail_t, tail, 2)[1], 4)

    def test_blank_is_a_missed_sample_not_zero(self):
        prof = [["S", x, conc(x)] for x in TIMES]
        prof[5][2] = ""
        code, p, _ = run(prof)
        code2, p2, _ = run([r for r in prof if r[2] != ""])
        self.assertAlmostEqual(first(p)["auc_0_t"], first(p2)["auc_0_t"], places=10)
        self.assertEqual(first(p)["n_missing"], 1)
        self.assertEqual(first(p)["n_blq"], 1)          # only the measured 0 at time 0
        self.assertEqual(table(p, "missed samples")[0]["time"], 2)
        self.assertTrue(any("protocol deviation" in f for f in p["findings"]))

    def test_blq_zero_in_auc_and_omitted_from_kel(self):
        prof = [["S", x, conc(x)] for x in TIMES]
        prof[9][2] = "BLQ"                               # 8 h, between quantifiable values
        code, p, _ = run(prof, "--auc-method", "linear")
        r = first(p)
        t = np.array(TIMES, float)
        c = np.array([0 if i == 9 else conc(x) for i, x in enumerate(TIMES)])
        self.assertAlmostEqual(r["auc_0_t"], trapezoid(c, t), places=9)
        self.assertLessEqual(r["lz_start"], 6)          # window spans the BLQ time
        self.assertEqual(r["lz_n"], sum(1 for x in TIMES if x > 2) - 1)
        self.assertTrue(any("BLQ between quantifiable" in f for f in p["findings"]))

    def test_lloq_threshold(self):
        code, p, _ = run([["S", x, conc(x)] for x in TIMES], "--lloq", "0.1")
        r = first(p)
        self.assertEqual(r["tlast"], max(x for x in TIMES if conc(x) >= .1))


class Rules(unittest.TestCase):
    def profiles(self, n, tail_cut=None, predose=None):
        rows = []
        for s in range(n):
            for x in TIMES:
                if tail_cut and s < tail_cut[0] and x > tail_cut[1]:
                    continue
                v = predose if (x == 0 and predose and s == 0) else conc(x)
                rows.append([f"S{s}", x, v])
        return rows

    def test_coverage_80_20_rule(self):
        # truncating after 6 h leaves AUC(0-t) well below 80% of AUC(0-inf)
        code, p, _ = run(self.profiles(10, (2, 6)))
        self.assertEqual(p["scalars"]["coverage_below_80"], 2)
        self.assertFalse(any("80%" in f for f in p["findings"]))          # 20% is not > 20%
        code, p, _ = run(self.profiles(10, (3, 6)))
        self.assertTrue(any("less than 80%" in f and "3 of 10" in f for f in p["findings"]))

    def test_predose_rule(self):
        code, p, _ = run(self.profiles(2, predose=0.06 * conc(2)))
        self.assertAlmostEqual(first(p)["predose_pct_cmax"], 6.0, places=6)
        code, p, _ = run(self.profiles(2, predose=0.04 * conc(2)))
        self.assertNotIn("predose_pct_cmax", first(p))

    def test_steady_state(self):
        tau = 12
        ts = [0, .5, 1, 2, 3, 4, 6, 8, 12.1]
        ss = lambda x: conc(x) + 3.0
        code, p, _ = run([["S", x, ss(x)] for x in ts], "--tau", str(tau))
        r = first(p)
        cs = [ss(x) for x in ts]
        c12 = cs[-2] + (cs[-1] - cs[-2]) * (12 - 8) / (12.1 - 8)
        auc = trapezoid(cs[:-1] + [c12], ts[:-1] + [12])
        self.assertAlmostEqual(r["auc_tau"], auc, places=9)
        self.assertAlmostEqual(r["ctau_ss"], cs[-1])
        self.assertEqual(r["ctau_time"], 12.1)
        self.assertAlmostEqual(r["cmin_ss"], min(cs))
        self.assertAlmostEqual(r["cav_ss"], auc / tau)
        self.assertAlmostEqual(r["fluctuation_pct"], 100 * (max(cs) - min(cs)) / (auc / tau))
        self.assertAlmostEqual(r["swing"], (max(cs) - min(cs)) / min(cs))
        for single_dose_only in ("auc_0_inf", "auc_coverage_pct", "auc_pct_extrap"):
            self.assertNotIn(single_dose_only, r)
        code, p, _ = run([["S", x, ss(x)] for x in ts[:-1] + [12.3]], "--tau", "12")
        self.assertNotIn("ctau_ss", first(p))
        self.assertTrue(any("within 10 min" in f for f in p["findings"]))

    def test_baseline_correction(self):
        base = [["S", -1, 2.0], ["S", -0.5, 4.0]] + [["S", x, 3.0 + conc(x) if x < 24 else 2.0] for x in TIMES[1:]]
        code, p, _ = run(base, "--baseline", "predose-mean")
        r = first(p)
        self.assertAlmostEqual(r["baseline"], 3.0)
        self.assertAlmostEqual(r["cmax"], conc(2))
        self.assertEqual(r["tlast"], 16)                   # 2.0 - 3.0 < 0 -> zero from 24 h
        self.assertNotIn("predose_pct_cmax", r)
        self.assertAlmostEqual(r["cmax_uncorrected"], 3.0 + conc(2))

    def test_gcc_emesis_rule(self):
        rows = []
        for s in range(6):
            for x in TIMES:
                rows.append([f"S{s}", 1, "T", x, conc(x), "1.5" if s == 0 else ("5" if s == 1 else "")])
        cols = ("subject", "period", "treatment", "time", "conc", "emesis_time")
        code, p, _ = run(rows, "--profile", "sfda", cols=cols)
        rs = table(p, "per-profile parameters")
        self.assertTrue(rs[0]["status"].startswith("excluded: emesis at 1.5"))   # 1.5 <= 2 x 2 h
        self.assertEqual(rs[1]["status"], "included")                           # 5 > 4
        code, p, _ = run(rows, cols=cols)
        self.assertTrue(all(r["status"] == "included" for r in table(p, "per-profile parameters")))
        self.assertTrue(any("pre-specified in the protocol" in f for f in p["findings"]))

    def test_gcc_72h_quantifiable(self):
        rows = [["S", x, conc(x) if x < 72 else "BLQ"] for x in TIMES]
        code, p, _ = run(rows, "--truncate-72", "--profile", "sfda")
        self.assertTrue(any("72 h is not quantifiable" in f for f in p["findings"]))

    def test_input_guards(self):
        self.assertEqual(run([["S", 0, 0], ["S", 1, -2]])[0], 2)
        self.assertEqual(run([["S", "", 1]])[0], 2)
        self.assertEqual(run([["S", 1, 1], ["S", 1, 2]])[0], 2)
        self.assertEqual(run([["S", x, conc(x)] for x in TIMES], "--lambda-z-points", "2")[0], 2)

    def test_json_has_no_nan(self):
        code, p, _ = run([["S", 0, 0], ["S", 1, "BLQ"], ["S", 2, 1]])
        json.dumps(p, allow_nan=False)


class BeHandoff(unittest.TestCase):
    def test_be_input_feeds_bioequivalence(self):
        rng = np.random.default_rng(3)
        rows = []
        for s in range(1, 19):
            seq = "TR" if s % 2 else "RT"
            for per, trt in enumerate(seq, 1):
                f = math.exp(rng.normal(0, .1)) * (0.97 if trt == "T" else 1)
                rows += [[f"S{s}", per, trt, x, f * conc(x)] for x in TIMES]
        code, p, _ = run(rows, "--profile", "sfda", cols=("subject", "period", "treatment", "time", "conc"))
        bi = table(p, "be input")
        self.assertEqual(len(bi), 36)
        self.assertEqual({r["sequence"] for r in bi}, {"TR", "RT"})
        fd, path = tempfile.mkstemp(suffix=".csv")
        os.close(fd)
        try:
            with open(path, "w", newline="") as fh:
                w = csv.DictWriter(fh, list(bi[0]))
                w.writeheader()
                w.writerows(bi)
            out = io.StringIO()
            with redirect_stdout(out), redirect_stderr(io.StringIO()):
                code = be.main_wrapper(be.run, ["-i", path, "--value-column", "cmax", "--metric", "cmax",
                                                "--predose-column", "predose", "--profile", "sfda", "--format", "json"])
            res = json.loads(out.getvalue())
        finally:
            os.remove(path)
        self.assertEqual(res["scalars"]["n_subjects"], 18)
        cmax = {(r["subject"], r["treatment"]): r["cmax"] for r in bi}
        diffs = [math.log(cmax[(f"S{s}", "T")]) - math.log(cmax[(f"S{s}", "R")]) for s in range(1, 19)]
        self.assertAlmostEqual(math.exp(np.mean(diffs)) * 100, res["scalars"]["gmr_pct"], places=6)

    def test_steady_state_be_input_offers_ss_parameters_only(self):
        ts = [0, .5, 1, 2, 4, 8, 12]
        rows = [[s, per, trt, x, conc(x) + 3] for s in ("A", "B") for per, trt in ((1, "T"), (2, "R")) for x in ts]
        code, p, _ = run(rows, "--tau", "12", cols=("subject", "period", "treatment", "time", "conc"))
        self.assertEqual(next(t["columns"] for t in p["tables"] if t["title"] == "be input"),
                         ["subject", "sequence", "period", "treatment", "auc_tau", "cmax_ss", "predose"])

    def test_excluded_period_drops_subject_from_be_input(self):
        rows = []
        for s in range(4):
            for per, trt in enumerate("TR", 1):
                rows += [[f"S{s}", per, trt, x, conc(x), "1" if (s == 0 and per == 1) else ""] for x in TIMES]
        code, p, _ = run(rows, "--profile", "sfda", cols=("subject", "period", "treatment", "time", "conc", "emesis_time"))
        self.assertNotIn("S0", {r["subject"] for r in table(p, "be input")})
        self.assertTrue(any("subject S0: a period was excluded" in n for n in p["notes"]))


def be_run(rows, *argv):
    """Write a 'be input' table to CSV and run bioequivalence.py on it; return (exit code, payload, stderr)."""
    fd, path = tempfile.mkstemp(suffix=".csv")
    os.close(fd)
    try:
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        out, err = io.StringIO(), io.StringIO()
        with redirect_stdout(out), redirect_stderr(err):
            code = be.main_wrapper(be.run, ["-i", path, "--format", "json", *argv])
        return code, (json.loads(out.getvalue()) if out.getvalue() else None), err.getvalue()
    finally:
        os.remove(path)


def crossover(n, seqs=("TR", "RT"), labels=None, **tweak):
    """Rows (subject, period, treatment, time, conc) for n subjects; tweak[(subject, period)] = callable(rows)."""
    rng = np.random.default_rng(11)
    rows = []
    for s in range(n):
        seq = seqs[s % len(seqs)]
        for per, trt in enumerate(seq, 1):
            f = math.exp(rng.normal(0, .1))
            lab = labels.get(trt, trt) if labels else trt
            rows += [[f"S{s:02d}", per, lab, x, f * conc(x)] for x in TIMES[:-1]]
    return rows


COLS = ("subject", "period", "treatment", "time", "conc")


class ReviewFixes(unittest.TestCase):
    """Regression tests for the findings of the 2026-09-26 review."""

    def test_profile_quantifiable_only_at_time_zero_does_not_crash(self):
        code, p, err = run([["A", 0, 5], ["A", 1, ""], ["A", 2, ""]])
        self.assertNotIn("Traceback", err)
        self.assertEqual(first(p)["cmax"], 5)

    def test_blank_period_or_treatment_is_an_input_error(self):
        rows = [["A", 1, "T", 0, 0], ["A", "", "T", 1, 10], ["A", 1, "T", 2, 8]]
        self.assertEqual(run(rows, cols=COLS)[0], 2)
        rows = [["A", 1, "T", 0, 0], ["A", 1, "", 1, 10]]
        self.assertEqual(run(rows, cols=COLS)[0], 2)

    def test_test_ref_labels_give_canonical_sequences(self):
        code, p, _ = run(crossover(12, labels={"T": "Test", "R": "Ref"}), cols=COLS)
        bi = table(p, "be input")
        self.assertEqual({r["sequence"] for r in bi}, {"TR", "RT"})
        code, res, err = be_run(bi, "--value-column", "cmax", "--metric", "cmax")
        self.assertIn(code, (0, 1), err)
        self.assertEqual(res["scalars"]["n_subjects"], 12)

    def test_subject_missing_a_period_leaves_the_be_input(self):
        rows = [r for r in crossover(12) if not (r[0] == "S03" and r[1] == 2)]
        code, p, _ = run(rows, cols=COLS)
        self.assertNotIn("S03", {r["subject"] for r in table(p, "be input")})
        self.assertTrue(any("subject S03: 1 of 2 periods" in n for n in p["notes"]))

    def test_column_options_are_case_insensitive(self):
        rows = crossover(2)
        code, p, _ = run(rows, "--conc-column", "Conc", "--period-column", "Period", cols=("subject", "Period", "treatment", "time", "Conc"))
        self.assertEqual(code in (0, 1), True)
        self.assertEqual({r["period"] for r in table(p, "per-profile parameters")}, {"1", "2"})

    def test_non_finite_numbers_are_rejected(self):
        self.assertEqual(run([["A", 0, "inf"], ["A", 1, 5]])[0], 2)
        self.assertEqual(run([["A", 0, 0], ["A", 1, 5], ["A", 2, 3], ["A", 4, 1]], "--tau", "nan")[0], 2)

    def test_nd_and_measured_zero_are_blq(self):
        prof = [["S", x, conc(x)] for x in TIMES]
        prof[9][2] = "ND"
        code, p, _ = run(prof)
        self.assertEqual(first(p)["n_blq"], 2)          # time 0 (measured 0) and 8 h (ND)
        prof[9][2] = 0
        code, p, _ = run(prof)
        self.assertEqual(first(p)["n_blq"], 2)
        self.assertTrue(any("BLQ between quantifiable" in f for f in p["findings"]))

    def test_replicate_predose_drops_the_subject_so_the_be_analysis_runs(self):
        rows = crossover(12, seqs=("TRR", "RTR", "RRT"))
        for r in rows:
            if r[0] == "S01" and r[1] == 3 and r[3] == 0:
                r[4] = 0.08 * conc(2)
        code, p, _ = run(rows, cols=COLS)
        bi = table(p, "be input")
        self.assertNotIn("S01", {r["subject"] for r in bi})
        code, res, err = be_run(bi, "--design", "partial", "--value-column", "auc_0_t", "--metric", "auc",
                                "--predose-column", "predose", "--cmax-column", "cmax")
        self.assertIn(code, (0, 1), err)

    def test_multi_treatment_keeps_subject_with_excluded_unrelated_period(self):
        rows = [r + [""] for r in crossover(18, seqs=("TAB", "ABT", "BTA", "TBA", "ATB", "BAT"), labels={"A": "R1", "B": "R2"})]
        for r in rows:
            if r[0] == "S01" and r[2] == "R2":
                r[5] = "0.5"
        code, p, _ = run(rows, "--profile", "sfda", cols=COLS + ("emesis_time",))
        bi = table(p, "be input")
        self.assertEqual({(r["subject"], r["treatment"]) for r in bi if r["subject"] == "S01"}, {("S01", "T"), ("S01", "R1")})
        code, res, err = be_run(bi, "--design", "multi", "--test-label", "T", "--reference-label", "R1",
                                "--value-column", "cmax", "--metric", "cmax", "--profile", "sfda")
        self.assertIn(code, (0, 1), err)
        self.assertEqual(res["scalars"]["n_subjects"], 18)

    def test_parallel_without_period_column(self):
        rows = [[f"S{s:02d}", "T" if s % 2 else "R", x, conc(x) * (1 + s / 100)] for s in range(24) for x in TIMES[:-1]]
        code, p, _ = run(rows, cols=("subject", "treatment", "time", "conc"))
        bi = table(p, "be input")
        self.assertEqual(len(bi), 24)
        code, res, err = be_run(bi, "--design", "parallel", "--value-column", "cmax", "--metric", "cmax")
        self.assertIn(code, (0, 1), err)

    def test_missing_metric_is_named(self):
        rows = crossover(12)
        rows = [r for r in rows if not (r[0] == "S02" and r[1] == 1 and r[3] > 3)]   # no terminal phase
        code, p, _ = run(rows, cols=COLS)
        self.assertTrue(any("auc_0_inf is not available for subject(s) S02" in n for n in p["notes"]))
        code, res, err = be_run(table(p, "be input"), "--value-column", "auc_0_inf", "--metric", "aucinf")
        self.assertEqual(code, 2)
        self.assertIn("subject S02 period 1", err)

    def test_auc_0_inf_listed_after_primary_parameters_and_partial_auc_included(self):
        code, p, _ = run(crossover(2), "--partial-auc", "0-2", cols=COLS)
        cols = next(t["columns"] for t in p["tables"] if t["title"] == "be input")
        self.assertEqual(cols, ["subject", "sequence", "period", "treatment", "cmax", "auc_0_t", "auc_0_inf", "auc_0_2", "predose"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
