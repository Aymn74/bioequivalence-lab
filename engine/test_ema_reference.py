"""ABEL (EMA Method A) against the EMA Q&A reference datasets (see validation/ema_abel/README.md)."""
import contextlib
import io
import json
import os
import unittest

import bioequivalence as be

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "validation", "ema_abel")


def abel(name, design):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = be.main_wrapper(be.run, ["-i", os.path.join(HERE, name), "--design", design, "--scaling", "abel",
                                        "--analysis", "ema", "--metric", "cmax", "--abel-justified", "--format", "json"])
    if not out.getvalue():
        return code, None, err.getvalue()
    res = json.loads(out.getvalue())
    row = next(r for t in res["tables"] if t["title"] == "scaled criteria" for r in t["rows"] if r["criterion"].startswith("EMA"))
    return code, row, err.getvalue()


class EmaReference(unittest.TestCase):
    def test_dataset_ii_partial_replicate_matches_ema(self):
        # EMA Annex III: CVwR 11.2%, PE 102.26%, 90% CI 97.32-107.46% (SAS Proc GLM, Method A)
        code, r, _ = abel("rds02.csv", "partial")
        self.assertEqual(round(100 * r["cvwr"], 1), 11.2)
        self.assertEqual(round(100 * r["gmr"], 2), 102.26)
        self.assertEqual((round(100 * r["ci_low"], 2), round(100 * r["ci_high"], 2)), (97.32, 107.46))
        self.assertFalse(r["widened"])
        self.assertTrue(r["met"])
        # replicateBE 1.1.3 method.A, unrounded
        self.assertAlmostEqual(100 * r["cvwr"], 11.1707611776742, places=9)
        self.assertAlmostEqual(100 * r["ci_low"], 97.3155468707817, places=9)
        self.assertAlmostEqual(100 * r["ci_high"], 107.464919793407, places=9)

    def test_dataset_i_is_incomplete_and_rejected(self):
        # EMA Annex II has subjects with missing periods; the engine analyses complete designs only
        code, r, err = abel("rds01.csv", "full")
        self.assertEqual(code, 2)
        self.assertIsNone(r)


if __name__ == "__main__":
    unittest.main(verbosity=2)
