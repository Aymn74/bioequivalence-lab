"""Compare nca.py with PKNCA output (run make_pknca_data.py, then Rscript run_pknca.R, then this).

Writes pknca_comparison.json and prints one line per parameter: values compared,
largest relative difference, and every mismatch above the tolerance.
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
STEADY = {("cmax", "12"): "cmax_ss", ("tmax", "12"): "tmax_ss", ("cmin", "12"): "cmin_ss", ("aucint.last", "12"): "auc_tau"}


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
        key = mapping.get((r["PPTESTCD"], r["end"].lower()))
        if key:
            ref[r["id"]][key] = None if r["PPORRES"] in ("", "NA") else float(r["PPORRES"])
    return ref


def compare(ours, ref, keys):
    report = {}
    for key in keys:
        n, worst, bad = 0, 0.0, []
        for pid, vals in ref.items():
            a, b = ours[pid].get(key), vals.get(key)
            if a is None and b is None:
                continue
            if a is None or b is None:
                bad.append({"id": pid, "nca_py": a, "pknca": b})
                continue
            n += 1
            rel = abs(a - b) / max(abs(b), 1e-12)
            worst = max(worst, rel)
            if rel > TOL:
                bad.append({"id": pid, "nca_py": a, "pknca": b, "rel_diff": rel})
        report[key] = {"compared": n, "max_rel_diff": worst, "mismatches": bad}
    return report


def main():
    result = {}
    for tag, method in (("linear", "linear"), ("linlog", "linup-logdown")):
        ours = run_nca(HERE / "single_dose.csv", "--auc-method", method, "--partial-auc", "0-2", "--partial-auc", "0-24")
        result[f"single_{tag}"] = compare(ours, read_pknca(f"pknca_single_{tag}.csv", SINGLE), SINGLE.values())
        result[f"convention_aucint_last_{tag}"] = compare(ours, read_pknca(f"pknca_single_{tag}.csv", LAST), LAST.values())
        ours = run_nca(HERE / "steady_state.csv", "--auc-method", method, "--tau", "12")
        result[f"ss_{tag}"] = compare(ours, read_pknca(f"pknca_ss_{tag}.csv", STEADY), STEADY.values())
    json.dump(result, open(HERE / "pknca_comparison.json", "w"), indent=1)
    for block, params in result.items():
        for key, r in params.items():
            print(f"{block:14s} {key:15s} n={r['compared']:4d} max_rel={r['max_rel_diff']:.1e} mismatches={len(r['mismatches'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
