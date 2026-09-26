"""Potency (assay content) correction: ICH M13A 2.2.2.3, GCC DS-G-010 3.1.2 / 3.1.8."""
import math
import unittest

from test_bioequivalence import fixture
from test_sfda import run_cli, table

PT, PR = 94.0, 101.0          # 7 percentage points apart


def analyse(rows, *argv, extra=()):
    code, p, err = run_cli(rows, *argv, extra=extra)
    return code, p, err


class PotencyCorrection(unittest.TestCase):
    def setUp(self):
        self.rows = fixture('2x2', [12, 12])

    def test_correction_shifts_the_estimate_by_the_potency_ratio_only(self):
        _, raw, _ = analyse(self.rows, '--metric', 'auc')
        code, cor, _ = analyse(self.rows, '--metric', 'auc', '--potency-test', str(PT), '--potency-reference', str(PR), '--potency-correction')
        shift = math.log(PR / PT)
        self.assertAlmostEqual(cor['scalars']['estimate'], raw['scalars']['estimate'] + shift, places=12)
        self.assertAlmostEqual(cor['scalars']['se'], raw['scalars']['se'], places=14)          # CI width unchanged
        for k in ('gmr_pct', 'ci90_low_pct', 'ci90_high_pct'):
            self.assertAlmostEqual(cor['scalars'][k], raw['scalars'][k] * PR / PT, places=9)
        both = table(cor, 'potency correction')
        self.assertEqual([r['data'] for r in both], ['uncorrected', 'potency-corrected'])
        self.assertAlmostEqual(both[0]['gmr_pct'], raw['scalars']['gmr_pct'], places=10)
        self.assertTrue(cor['scalars']['potency_corrected'])

    def test_fixed_effects_and_parallel_models_shift_the_same_way(self):
        par = [dict(subject=f'P{i}', treatment='T' if i % 2 else 'R', value=math.exp(4 + 0.3 * math.sin(i)), sequence='', period='1')
               for i in range(28)]
        for rows, argv in ((self.rows, ('--analysis', 'ema')), (par, ('--design', 'parallel'))):
            _, raw, _ = analyse(rows, '--metric', 'auc', *argv)
            _, cor, _ = analyse(rows, '--metric', 'auc', *argv, '--potency-test', str(PT), '--potency-reference', str(PR), '--potency-correction')
            self.assertAlmostEqual(cor['scalars']['gmr_pct'], raw['scalars']['gmr_pct'] * PR / PT, places=9)

    def test_abel_uses_corrected_estimate_and_unchanged_cvwr(self):
        rows = fixture('partial', [8, 8, 8])
        argv = ('--design', 'partial', '--metric', 'cmax', '--scaling', 'abel', '--abel-justified')
        _, raw, _ = analyse(rows, *argv)
        _, cor, _ = analyse(rows, *argv, '--potency-test', str(PT), '--potency-reference', str(PR), '--potency-correction')
        a, b = table(raw, 'scaled criteria')[0], table(cor, 'scaled criteria')[0]
        self.assertAlmostEqual(a['cvwr'], b['cvwr'], places=12)
        self.assertAlmostEqual(b['gmr'], a['gmr'] * PR / PT, places=10)

    def test_within_five_points_no_correction(self):
        code, p, _ = analyse(self.rows, '--metric', 'auc', '--potency-test', '98', '--potency-reference', '102', '--potency-correction')
        self.assertFalse(p['scalars']['potency_corrected'])
        self.assertTrue(any('within 5%' in n and 'ignored' in n for n in p['notes']))
        self.assertFalse(any(t['title'] == 'potency correction' for t in p['tables']))

    def test_difference_above_five_points_without_correction_is_a_finding(self):
        code, p, _ = analyse(self.rows, '--metric', 'auc', '--potency-test', str(PT), '--potency-reference', str(PR))
        self.assertFalse(p['scalars']['potency_corrected'])
        self.assertTrue(any('differ by 7.00 percentage points' in f for f in p['findings']))

    def test_predose_rule_uses_measured_values(self):
        rows = fixture('2x2', [12, 12])
        for r in rows: r['predose'] = 0.0
        rows[0]['predose'] = 0.051 * rows[0]['value']        # just above 5% of the measured Cmax
        code, p, _ = analyse(rows, '--metric', 'cmax', '--potency-test', '80', '--potency-reference', '100', '--potency-correction', extra=('predose',))
        self.assertAlmostEqual(table(p, 'pre-dose exclusions')[0]['ratio_pct'], 5.1, places=6)

    def test_input_guards(self):
        self.assertEqual(analyse(self.rows, '--potency-test', '95')[0], 2)
        self.assertEqual(analyse(self.rows, '--potency-correction')[0], 2)
        self.assertEqual(analyse(self.rows, '--potency-test', '0', '--potency-reference', '100')[0], 2)


if __name__ == '__main__':
    unittest.main(verbosity=2)
