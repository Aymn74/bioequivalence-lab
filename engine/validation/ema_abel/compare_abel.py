"""Compare bioequivalence.py ABEL (EMA Method A) with replicateBE and with EMA's published results.

Run Rscript run_replicatebe.R . first. Writes abel_comparison.json, prints a table, exits 1 on any
numeric difference beyond the tolerance in a dataset the engine accepts, or on a mismatch with the
EMA-published values (rounded as published).
"""
import contextlib
import csv
import io
import json
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent))
import bioequivalence as be  # noqa: E402

TOL = 1e-8
DESIGNS = {frozenset({"TRR", "RTR", "RRT"}): "partial", frozenset({"TRTR", "RTRT"}): "full"}
# EMA Q&A reference datasets, results reported by EMA (SAS Proc GLM, Method A), as quoted in replicateBE 1.1.3
EMA = {"rds01": {"CVwR": 47.0, "PE": 115.66, "CL.lo": 107.11, "CL.hi": 124.89},
       "rds02": {"CVwR": 11.2, "PE": 102.26, "CL.lo": 97.32, "CL.hi": 107.46}}


def inspect(path):
    rows = list(csv.DictReader(open(path, encoding="utf-8")))
    seqs = frozenset(r["sequence"] for r in rows)
    per_subject = defaultdict(int)
    for r in rows:
        if r["value"] != "":
            per_subject[r["subject"]] += 1
    periods = len({r["period"] for r in rows})
    complete = all(v == periods for v in per_subject.values()) and all(r["value"] != "" for r in rows)
    return DESIGNS.get(seqs), complete, "|".join(sorted(seqs))


def engine(path, design):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = be.main_wrapper(be.run, ["-i", str(path), "--design", design, "--scaling", "abel", "--analysis", "ema",
                                        "--metric", "cmax", "--abel-justified", "--format", "json"])
    if code == 2 or not out.getvalue():
        return code, None, err.getvalue().strip().splitlines()[-1] if err.getvalue().strip() else ""
    res = json.loads(out.getvalue())
    abel = next(r for t in res["tables"] if t["title"] == "scaled criteria" for r in t["rows"] if r["criterion"].startswith("EMA"))
    return code, {"CVwR": 100 * abel["cvwr"], "EL.lo": 100 * abel["low"], "EL.hi": 100 * abel["high"], "PE": 100 * abel["gmr"],
                  "CL.lo": 100 * abel["ci_low"], "CL.hi": 100 * abel["ci_high"], "BE": "pass" if abel["met"] else "fail"}, ""


def main():
    ref = {r["set"]: r for r in csv.DictReader(open(HERE / "replicatebe_method_a.csv", encoding="utf-8"))}
    report, problems = [], []
    for name, r in ref.items():
        design, complete, seqs = inspect(HERE / f"{name}.csv")
        entry = {"set": name, "sequences": seqs, "n": int(r["n"]), "complete": complete, "engine_design": design}
        if design is None:
            code, vals, msg = engine(HERE / f"{name}.csv", "full")
            entry.update(status="rejected by the engine (unsupported design)", exit_code=code, message=msg)
            if code != 2:
                problems.append(f"{name}: expected a rejection (exit 2), got exit {code}")
            report.append(entry)
            continue
        code, vals, msg = engine(HERE / f"{name}.csv", design)
        if vals is None:
            problems.append(f"{name}: engine failed: {msg}")
            continue
        diffs = {k: abs(vals[k] - float(r[k])) / abs(float(r[k])) for k in ("CVwR", "EL.lo", "EL.hi", "PE", "CL.lo", "CL.hi")}
        worst = max(diffs.values())
        entry.update(status="compared", engine=vals, replicatebe={k: float(r[k]) for k in diffs}, max_rel_diff=worst,
                     decision_engine=vals["BE"], decision_replicatebe=r["BE"])
        if worst > TOL:
            problems.append(f"{name}: max relative difference {worst:.2e}")
        if vals["BE"] != r["BE"]:
            problems.append(f"{name}: decision {vals['BE']} vs replicateBE {r['BE']}")
        base = name.replace("_complete", "")
        if base in EMA and name == base:
            pub = EMA[base]
            ok = (round(vals["CVwR"], 1) == pub["CVwR"] and all(round(vals[k], 2) == pub[k] for k in ("PE", "CL.lo", "CL.hi")))
            entry["ema_published"] = pub
            entry["ema_match"] = ok
            if not ok:
                problems.append(f"{name}: differs from the EMA-published result {pub}")
        report.append(entry)
    json.dump(report, open(HERE / "abel_comparison.json", "w"), indent=1)
    for e in report:
        if e["status"] == "compared":
            ema = "" if "ema_match" not in e else ("  EMA published: match" if e["ema_match"] else "  EMA published: MISMATCH")
            ema += "" if e["complete"] else "  (incomplete)"
            print(f"{e['set']:16s} {e['sequences']:22s} n={e['n']:3d}  CVwR {e['engine']['CVwR']:7.3f}%  "
                  f"PE {e['engine']['PE']:7.3f}%  CI {e['engine']['CL.lo']:.3f}-{e['engine']['CL.hi']:.3f}%  "
                  f"limits {e['engine']['EL.lo']:.2f}-{e['engine']['EL.hi']:.2f}%  {e['engine']['BE']:4s}  max rel diff {e['max_rel_diff']:.1e}{ema}")
        else:
            print(f"{e['set']:16s} {e['sequences']:22s} n={e['n']:3d}  {e['status']}")
    for p in problems:
        print("PROBLEM:", p)
    n = sum(e["status"] == "compared" for e in report)
    inc = sum(e["status"] == "compared" and not e["complete"] for e in report)
    print(f"{n} datasets compared ({inc} incomplete), {len(report) - n} rejected as expected (unsupported designs)" if not problems else f"FAILED: {len(problems)} problem(s)")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
