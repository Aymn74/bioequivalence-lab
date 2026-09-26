"""Build the datasets for the nca.py vs PKNCA comparison (see README.md in this folder).

Writes, next to this file:
  single_dose.csv   example_nca.csv profiles + 300 random one-compartment profiles with BLQ values
  steady_state.csv  60 random steady-state profiles (tau = 12 h; trough at 12 h, or at 12.05 h for S50-S59)
Each file has id, time, conc (BLQ written as 0) and blq (1/0); missed samples are omitted.
"""
import csv
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
rng = np.random.default_rng(20260926)

rows = []
with open(HERE.parent / "example_nca.csv", encoding="utf-8") as fh:
    lines = [ln for ln in fh if not ln.startswith("#")]
for r in csv.DictReader(lines):
    txt = r["conc"].strip()
    if not txt:
        continue                                   # missed sample
    blq = txt.upper() == "BLQ"
    rows.append([f"EX-{r['subject']}-P{r['period']}", float(r["time"]), 0.0 if blq else float(txt), int(blq)])

for k in range(300):
    ka, ke = rng.uniform(.3, 3), rng.uniform(.02, .4)
    if abs(ka - ke) < .05:
        ka += .1
    ts = sorted(set(np.round(np.r_[0, rng.uniform(.1, 72, 13)], 3)))
    cs = [0.0 if x == 0 else 50 * ka / (ka - ke) * (math.exp(-ke * x) - math.exp(-ka * x)) * math.exp(rng.normal(0, .15)) for x in ts]
    lloq = sorted(cs)[len(cs) // 4]
    for x, c in zip(ts, cs):
        blq = x == 0 or c < lloq
        rows.append([f"R{k:03d}", float(x), 0.0 if blq else round(c, 6), int(blq)])

with open(HERE / "single_dose.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["id", "time", "conc", "blq"])
    w.writerows(rows)

ss = []
for k in range(60):
    ka, ke, base = rng.uniform(.5, 3), rng.uniform(.05, .3), rng.uniform(1, 5)
    # S00-S49: trough sample exactly at tau = 12 h; S50-S59: trough sample 3 min after tau (12.05 h), which nca.py uses (within 10 min) and PKNCA's [0, 12] interval does not
    ts = [0, .5, 1, 1.5, 2, 3, 4, 6, 8, 10, 12 if k < 50 else 12.05]
    for x in ts:
        # true steady state of a one-compartment oral model (superposition of all previous doses), so C(0) = C(tau)
        c = base * 20 * ka / (ka - ke) * (math.exp(-ke * x) / (1 - math.exp(-ke * 12)) - math.exp(-ka * x) / (1 - math.exp(-ka * 12)))
        ss.append([f"S{k:02d}", x, round(c * math.exp(rng.normal(0, .08)), 6), 0])
with open(HERE / "steady_state.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["id", "time", "conc", "blq"])
    w.writerows(ss)
print(len(rows), "single-dose rows;", len(ss), "steady-state rows")
