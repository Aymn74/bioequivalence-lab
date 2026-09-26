"""Regression tests for the pre-launch review (2026-09-26): NCA -> BE hand-off across designs,
subjects without both products, potency correction with scaled criteria, ABEL point-estimate rounding."""
import math
import unittest

import bioequivalence as b
from test_bioequivalence import fixture
from test_nca import COLS, TIMES, be_run, conc, crossover, run, table as nca_table
from test_sfda import run_cli, table


def with_column(rows, name, value_of):
    """Append a column to NCA rows (subject is rows[i][0])."""
    return [r + [value_of(r[0])] for r in rows]


class HandOff(unittest.TestCase):
    def test_stage_column_reaches_the_be_input_and_the_model(self):
        rows = with_column(crossover(24), "stage", lambda s: "1" if int(s[1:]) < 12 else "2")
        code, p, _ = run(rows, cols=COLS + ("stage",))
        bi = nca_table(p, "be input")
        self.assertEqual({r["stage"] for r in bi}, {"1", "2"})
        self.assertTrue(any("keeps the stage column" in n for n in p["notes"]))
        code, res, err = be_run(bi, "--value-column", "cmax", "--metric", "cmax", "--analysis", "ema", "--anova", "--ci-level", "94.12")
        self.assertIn(code, (0, 1), err)
        self.assertIn("period(stage)", {r["source"] for r in table(res, "anova")})
        direct = b.ema_analysis([dict(subject=r["subject"], treatment=r["treatment"], sequence=r["sequence"], period=r["period"],
                                      stage=r["stage"], value=r["cmax"], logvalue=math.log(r["cmax"])) for r in bi], "2x2", 94.12)
        self.assertAlmostEqual(res["scalars"]["ci90_low_pct"], 100 * direct.ci_low, places=10)

    def test_group_column_reaches_the_be_input_and_the_interaction_table(self):
        rows = with_column(crossover(24), "group", lambda s: "A" if int(s[1:]) % 3 else "B")
        code, p, _ = run(rows, cols=COLS + ("group",))
        bi = nca_table(p, "be input")
        code, res, err = be_run(bi, "--value-column", "cmax", "--metric", "cmax", "--analysis", "ema")
        self.assertIn(code, (0, 1), err)
        self.assertEqual(len(table(res, "group interaction")), 1)

    def test_stage_must_be_constant_within_a_subject(self):
        rows = [r + ["1" if r[1] == 1 else "2"] for r in crossover(4)]
        code, _, err = run(rows, cols=COLS + ("stage",))
        self.assertEqual(code, 2)
        self.assertIn("one stage value", err)

    def test_period_codes_sequence_codes_and_treatment_labels_are_normalised(self):
        rows = []
        for r in crossover(12, labels={"R": "C"}):
            rows.append([r[0], f"P{r[1]}", r[2], r[3], r[4], "1" if int(r[0][1:]) % 2 == 0 else "2"])
        cols = ("subject", "period", "treatment", "time", "conc", "sequence")
        code, p, _ = run(rows, "--reference-label", "C", cols=cols)
        bi = nca_table(p, "be input")
        self.assertEqual({r["period"] for r in bi}, {"1", "2"})
        self.assertEqual({r["treatment"] for r in bi}, {"T", "R"})
        self.assertEqual({r["sequence"] for r in bi}, {"TR", "RT"})
        code, res, err = be_run(bi, "--value-column", "cmax", "--metric", "cmax")
        self.assertIn(code, (0, 1), err)
        self.assertEqual(res["scalars"]["n_subjects"], 12)
        code, p, _ = run(rows, cols=cols)                          # C not declared as the reference
        self.assertTrue(any("not T and R" in n for n in p["notes"]))

    def test_partial_replicate_without_sequence_column_fits_drop_outs(self):
        rows = [r for r in crossover(12, seqs=("TRR", "RTR", "RRT")) if not (r[0] == "S01" and r[1] > 1)]   # RTR, only period 1 (R)
        code, p, _ = run(rows, cols=COLS)
        self.assertNotIn("S01", {r["subject"] for r in nca_table(p, "be input")})
        self.assertTrue(any("S01: sequence cannot be determined" in n for n in p["notes"]))
        rows = [r for r in crossover(12, seqs=("TRR", "RTR", "RRT")) if not (r[0] == "S00" and r[1] == 3)]   # TRR: T first, fits
        code, p, _ = run(rows, cols=COLS)
        self.assertEqual({r["sequence"] for r in nca_table(p, "be input") if r["subject"] == "S00"}, {"TRR"})

    def test_gcc_annex_with_test_reference_words(self):
        code, p, _ = run(crossover(12, labels={"T": "Test", "R": "Reference"}), "--profile", "sfda", cols=COLS)
        self.assertTrue(nca_table(p, "gcc annex 1"))

    def test_williams_drop_out_is_left_out_of_that_comparison_only(self):
        seqs = ("TRX", "RXT", "XTR", "TXR", "RTX", "XRT")
        rows = [r for r in crossover(18, seqs=seqs) if not (r[0] == "S00" and r[1] == 2)]   # S00 (TRX) misses R
        code, p, _ = run(rows, cols=COLS)
        bi = nca_table(p, "be input")
        argv = ("--design", "multi", "--value-column", "cmax", "--metric", "cmax")
        code, res, err = be_run(bi, *argv, "--test-label", "T", "--reference-label", "R")
        self.assertIn(code, (0, 1), err)
        self.assertEqual(res["scalars"]["n_subjects"], 17)
        self.assertEqual([r["subject"] for r in table(res, "subjects without both products")], ["S00"])
        code, res, err = be_run(bi, *argv, "--test-label", "T", "--reference-label", "X")
        self.assertEqual(res["scalars"]["n_subjects"], 18)


class OneProductSubjects(unittest.TestCase):
    def test_2x2_fixed_effects_leaves_out_subjects_with_one_period(self):
        rows = fixture('2x2', [12, 12])
        lone = [r for r in rows if not (r['subject'] in ('TR-0', 'RT-3') and r['period'] == '2')]
        complete = [r for r in rows if r['subject'] not in ('TR-0', 'RT-3')]
        _, a, _ = run_cli(lone, '--analysis', 'ema', '--anova')
        _, c, _ = run_cli(complete, '--analysis', 'ema', '--anova')
        self.assertEqual(a['scalars']['n_subjects'], 22)
        self.assertAlmostEqual(a['scalars']['ci90_low_pct'], c['scalars']['ci90_low_pct'], places=12)
        for x, y in zip(table(a, 'anova'), table(c, 'anova')):           # sequence and subject(sequence) too
            self.assertEqual((x['source'], x['df']), (y['source'], y['df']))
            self.assertAlmostEqual(x['ss'], y['ss'], places=10)
        self.assertEqual({r['subject'] for r in table(a, 'subjects without both products')}, {'TR-0', 'RT-3'})


class PotencyWithScaledCriteria(unittest.TestCase):
    def test_uncorrected_data_get_the_same_scaled_criterion(self):
        rows = fixture('partial', [8, 8, 8])
        argv = ('--design', 'partial', '--metric', 'cmax', '--scaling', 'abel', '--abel-justified')
        _, raw, _ = run_cli(rows, *argv)
        _, cor, _ = run_cli(rows, *argv, '--potency-test', '94', '--potency-reference', '101', '--potency-correction')
        pc = table(cor, 'potency correction')
        self.assertEqual([(r['data'], r['criterion']) for r in pc],
                         [('uncorrected', 'average BE (fixed limits)'), ('uncorrected', 'EMA ABEL (fixed-effects model)'),
                          ('potency-corrected', 'average BE (fixed limits)'), ('potency-corrected', 'EMA ABEL (fixed-effects model)')])
        abel_raw = table(raw, 'scaled criteria')[0]
        self.assertAlmostEqual(pc[1]['ci_low_pct'], 100 * abel_raw['ci_low'], places=10)
        self.assertAlmostEqual(pc[1]['limit_low_pct'], 100 * abel_raw['low'], places=10)
        self.assertEqual(pc[1]['met'], abel_raw['met'])


class AbelPointEstimate(unittest.TestCase):
    def test_gmr_constraint_compared_at_two_decimals_and_outlier_note(self):
        rows = fixture('full', [12, 12], cv=.6)
        argv = ('--design', 'full', '--metric', 'cmax', '--scaling', 'abel', '--abel-justified')
        _, p, _ = run_cli(rows, *argv)
        k = 1.250049 / table(p, 'scaled criteria')[0]['gmr']             # GMR 125.0049% -> 125.00% at two decimals
        shifted = [dict(r, value=r['value'] * k) if r['treatment'] == 'T' else r for r in rows]
        _, q, _ = run_cli(shifted, *argv)
        row = table(q, 'scaled criteria')[0]
        self.assertGreater(row['gmr'], 1.25)
        self.assertTrue(row['passes_gmr_constraint'])
        self.assertTrue(any('not the result of outliers' in n for n in q['notes']))


class ReviewRound2(unittest.TestCase):
    """Independent review of the pre-launch fixes (2026-09-26, second round)."""

    def test_parallel_with_stage_is_refused_not_ignored(self):
        rows = [[f"S{s:02d}", "T" if s % 2 else "R", x, conc(x) * (1 + s / 100), "1" if s < 12 else "2"] for s in range(24) for x in TIMES[:-1]]
        code, p, _ = run(rows, cols=("subject", "treatment", "time", "conc", "stage"))
        self.assertTrue(any("parallel BE model has no stage or group term" in n for n in p["notes"]))
        code, res, err = be_run(nca_table(p, "be input"), "--design", "parallel", "--value-column", "cmax", "--metric", "cmax")
        self.assertEqual(code, 2)
        self.assertIn("not supported for a parallel design", err)

    def test_2x2_coding_error_is_reported_not_excluded(self):
        rows = fixture('2x2', [12, 12])
        for r in rows:
            if r['subject'] == 'TR-1' and r['period'] == '2':
                r['treatment'] = 'T'
        code, _, err = run_cli(rows, '--analysis', 'ema')
        self.assertEqual(code, 2)
        self.assertIn('treatment conflicts with sequence/period', err)

    def test_swapped_test_and_reference_labels(self):
        _, p0, _ = run(crossover(12), cols=COLS)
        _, p1, _ = run(crossover(12), "--test-label", "R", "--reference-label", "T", cols=COLS)
        argv = ("--value-column", "cmax", "--metric", "cmax")
        _, a, _ = be_run(nca_table(p0, "be input"), *argv)
        code, b_, err = be_run(nca_table(p1, "be input"), *argv)
        self.assertIn(code, (0, 1), err)
        self.assertAlmostEqual(a["scalars"]["gmr_pct"] * b_["scalars"]["gmr_pct"], 1e4, places=6)
        _, p2, _ = run(crossover(12, labels={"T": "Test", "R": "Ref"}), "--test-label", "Ref", "--reference-label", "Test", cols=COLS)
        code, c, err = be_run(nca_table(p2, "be input"), *argv)
        self.assertIn(code, (0, 1), err)
        self.assertAlmostEqual(c["scalars"]["gmr_pct"], b_["scalars"]["gmr_pct"], places=9)

    def test_gcc_annex_for_multi_treatment_with_named_labels(self):
        seqs = ("TRX", "RXT", "XTR", "TXR", "RTX", "XRT")
        rows = crossover(18, seqs=seqs, labels={"T": "A", "R": "B", "X": "C"})
        code, p, _ = run(rows, "--profile", "sfda", "--test-label", "A", "--reference-label", "B", cols=COLS)
        self.assertTrue(nca_table(p, "gcc annex 1"))

    def test_stage_conflict_within_a_profile_and_blank_stage_column(self):
        rows = [r + ["1" if r[3] < 4 else "2"] for r in crossover(4)]
        code, _, err = run(rows, cols=COLS + ("stage",))
        self.assertEqual(code, 2)
        self.assertIn("two stage values in one profile", err)
        code, p, _ = run([r + [""] for r in crossover(4)], cols=COLS + ("stage",))
        self.assertNotIn("stage", nca_table(p, "be input")[0])

    def test_sequence_whose_subjects_all_dropped_out_is_fitted(self):
        rows = [r for r in crossover(12, seqs=("TRR", "RTR", "RRT")) if not (int(r[0][1:]) % 3 == 2 and r[1] == 3)]
        code, p, _ = run(rows, cols=COLS)
        bi = nca_table(p, "be input")
        self.assertEqual({r["sequence"] for r in bi if int(r["subject"][1:]) % 3 == 2}, {"RRT"})
        code, res, err = be_run(bi, "--design", "partial", "--value-column", "cmax", "--metric", "cmax", "--analysis", "ema")
        self.assertIn(code, (0, 1), err)

if __name__ == '__main__':
    unittest.main(verbosity=2)
