"""Compare nca.py with PKNCA output (run make_pknca_data.py, then Rscript run_pknca.R, then this).

Writes pknca_comparison.json, prints one line per parameter, and exits 1 when anything differs
from the documented expectation: a profile present on one side only, a count of compared values
that is not the full set, or a mismatch that is not one of the explained cases below.
"""
import contextlib
import csv
import io
import json
import math
import os
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import nca  # noqa: E402

TOL = 1e-9
SINGLE = {("auclast", "inf"): "auc_0_t", ("cmax", "inf"): "cmax", ("tmax", "inf"): "tmax", ("tlast", "inf"): "tlast",
          ("clast.obs", "inf"): "clast", ("lambda.z", "inf"): "lambda_z", ("half.life", "inf"): "t_half",
          ("lambda.z.n.points", "inf"): "lz_n", ("adj.r.squared", "inf"): "r2_adj", ("aucinf.obs", "inf"): "auc_0_inf",
          ("aucinf.pred", "inf"): "auc_0_inf_pred", ("aucint.all", "2"): "auc_0_2", ("aucint.all", "24"): "auc_0_24"}
# Convention check only: aucint.last sets every concentration after tlast to zero, while nca.py (like
# aucint.all) keeps a later BLQ sample as zero at its own time (ICH M13A 2.2.2.2) and interpolates to it.
LAST = {("aucint.last", "2"): "auc_0_2", ("aucint.last", "24"): "auc_0_24"}
STEADY = {("cmax", "12"): "cmax_ss", ("tmax", "12"): "tmax_ss", ("cmin", "12"): "cmin_ss", ("aucint.last", "12"): "auc_tau",
          ("ctrough", "12"): "ctau_ss", ("cav.int.last", "12"): "cav_ss", ("swing", "12"): "swing",
          ("fluctuation", "12"): "fluctuation_pct"}

# Explained differences (see README.md): {block: {parameter: set of profile ids}}. Anything else fails.
LATE_TROUGH = {f"S{k:02d}" for k in range(50, 60)}                 # trough sample at 12.05 h, outside PKNCA's [0, 12]
EXPECTED_COUNTS = {
    ("single_linear", "auc_0_24"): 15, ("single_linlog", "auc_0_24"): 8,
    ("convention_aucint_last_linear", "auc_0_24"): 28, ("convention_aucint_last_linlog", "auc_0_24"): 28,
}
LATE_TROUGH_KEYS = {"ctau_ss", "cmin_ss", "swing", "fluctuation_pct"}


def run_nca(src, *argv):
    rows = list(csv.DictReader(open(src, encoding="utf-8")))
    fd, path = tempfile.mkstemp(suffix=".csv")
    os.close(fd)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["subject", "time", "conc"])
        w.writerows([r["id"], r["time"], "BLQ" if r["blq"] == "1" else r["conc"]] for r in rows)
    out = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(io.StringIO()):
        nca.main_wrapper(nca.run, ["-i", path, "--format", "json", *argv])
    os.remove(path)
    table = next(t for t in json.loads(out.getvalue())["tables"] if t["title"] == "per-profile parameters")
    return {r["subject"]: r for r in table["rows"]}


def read_pknca(name, mapping):
    ref = defaultdict(dict)
    for r in csv.DictReader(open(HERE / name, encoding="utf-8")):
        ref[r["id"]]  # every profile PKNCA reports, even without mapped parameters
        key = mapping.get((r["PPTESTCD"], r["end"].lower()))
        if key:
            ref[r["id"]][key] = None if r["PPORRES"] in ("", "NA", "Inf") else float(r["PPORRES"])
    if mapping is STEADY:   # PKNCA swing is a percentage; ICH M13A defines swing as a ratio
        for vals in ref.values():
            if vals.get("swing") is not None:
                vals["swing"] /= 100.0
            cmax, cmin, cav = vals.get("cmax_ss"), vals.get("cmin_ss"), vals.get("cav_ss")
            vals["fluctuation_pct"] = 100 * (cmax - cmin) / cav if None not in (cmax, cmin, cav) and cav > 0 else None
    return ref


def compare(ours, ref, keys):
    report = {}
    for key in keys:
        n, both_none, worst, bad = 0, 0, 0.0, []
        for pid, vals in ref.items():
            a, b = ours[pid].get(key), vals.get(key)
            if a is None and b is None:
                both_none += 1
                continue
            if a is None or b is None:
                bad.append({"id": pid, "nca_py": a, "pknca": b})
                continue
            rel = abs(a - b) / max(abs(b), 1e-12)
            if rel > TOL:
                bad.append({"id": pid, "nca_py": a, "pknca": b, "rel_diff": rel})
            else:
                n += 1
                worst = max(worst, rel)
        report[key] = {"compared": n, "both_missing": both_none, "max_rel_diff_matching": worst, "mismatches": bad}
    return report


def check(result, profiles):
    problems = []
    for block, params in result.items():
        total = profiles["ss" if block.startswith("ss") else "single"]
        for key, r in params.items():
            if r["compared"] + r["both_missing"] + len(r["mismatches"]) != total:
                problems.append(f"{block} {key}: {r['compared']} + {r['both_missing']} + {len(r['mismatches'])} != {total} profiles")
            ids = {m["id"] for m in r["mismatches"]}
            expected = EXPECTED_COUNTS.get((block, key))
            if block.startswith("ss") and key in LATE_TROUGH_KEYS:
                if not ids <= LATE_TROUGH:
                    problems.append(f"{block} {key}: unexpected mismatches {sorted(ids - LATE_TROUGH)}")
            elif expected is not None:
                if len(ids) != expected:
                    problems.append(f"{block} {key}: {len(ids)} mismatches, expected {expected}")
            elif ids:
                problems.append(f"{block} {key}: unexpected mismatches {sorted(ids)[:5]}")
    return problems


def main():
    result = {}
    for tag, method in (("linear", "linear"), ("linlog", "linup-logdown")):
        ours = run_nca(HERE / "single_dose.csv", "--auc-method", method, "--partial-auc", "0-2", "--partial-auc", "0-24")
        ref = read_pknca(f"pknca_single_{tag}.csv", SINGLE)
        if set(ours) != set(ref):
            print("profile sets differ (single):", sorted(set(ours) ^ set(ref))[:10])
            return 1
        result[f"single_{tag}"] = compare(ours, ref, SINGLE.values())
        result[f"convention_aucint_last_{tag}"] = compare(ours, read_pknca(f"pknca_single_{tag}.csv", LAST), LAST.values())
        ours = run_nca(HERE / "steady_state.csv", "--auc-method", method, "--tau", "12")
        ref = read_pknca(f"pknca_ss_{tag}.csv", STEADY)
        if set(ours) != set(ref):
            print("profile sets differ (steady state):", sorted(set(ours) ^ set(ref))[:10])
            return 1
        result[f"ss_{tag}"] = compare(ours, ref, STEADY.values())
    json.dump(result, open(HERE / "pknca_comparison.json", "w"), indent=1)
    for block, params in result.items():
        for key, r in params.items():
            print(f"{block:30s} {key:15s} n={r['compared']:4d} max_rel={r['max_rel_diff_matching']:.1e} "
                  f"both_missing={r['both_missing']:3d} mismatches={len(r['mismatches'])}")
    problems = check(result, {"single": 348, "ss": 60})
    for p in problems:
        print("UNEXPECTED:", p)
    print("OK: all differences are the documented ones" if not problems else f"FAILED: {len(problems)} unexpected difference(s)")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
