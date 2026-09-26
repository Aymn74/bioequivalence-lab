"""Tests for the GCC/SFDA profile, fixed-effects ANOVA table, pre-dose rule,
adjusted confidence level and two-stage (stage) term."""
import csv
import io
import json
import math
import os
import tempfile
import unittest
from contextlib import redirect_stdout, redirect_stderr
import numpy as np
import bioequivalence as b
from test_bioequivalence import fixture


def run_cli(rows, *argv, extra=()):
    """Write rows to a temporary CSV, run the CLI with --format json, return (exit, payload, stderr)."""
    cols = ['subject', 'sequence', 'period', 'treatment', 'value', *extra]
    fd, path = tempfile.mkstemp(suffix='.csv'); os.close(fd)
    try:
        with open(path, 'w', newline='', encoding='utf-8') as fh:
            w = csv.DictWriter(fh, cols, extrasaction='ignore'); w.writeheader(); w.writerows(rows)
        out, err = io.StringIO(), io.StringIO()
        with redirect_stdout(out), redirect_stderr(err):
            code = b.main_wrapper(b.run, ['-i', path, *argv, '--format', 'json'])
        return code, (json.loads(out.getvalue()) if out.getvalue() else None), err.getvalue()
    finally:
        os.remove(path)


def table(payload, title):
    return next(t['rows'] for t in payload['tables'] if t['title'] == title)


class AnovaTests(unittest.TestCase):
    def test_2x2_unbalanced_textbook_partition(self):
        rows = fixture('2x2', [14, 9]); at = {r['source']: r for r in b.anova_table(rows, '2x2')}
        subj = {}
        for r in rows: subj.setdefault((r['subject'], r['sequence']), []).append(r['logvalue'])
        totals = {k: sum(v) for k, v in subj.items()}
        by_seq = {}
        for (s, q), v in totals.items(): by_seq.setdefault(q, []).append(v)
        ss_subj = sum((v - np.mean(vs)) ** 2 for vs in by_seq.values() for v in vs) / 2
        grand = np.mean(list(totals.values()))
        ss_seq = sum(len(vs) * (np.mean(vs) - grand) ** 2 for vs in by_seq.values()) / 2
        self.assertAlmostEqual(at['subject(sequence)']['ss'], ss_subj, 10)
        self.assertAlmostEqual(at['sequence']['ss'], ss_seq, 10)
        self.assertEqual(at['subject(sequence)']['df'], 21); self.assertEqual(at['residual']['df'], 21)
        # Balanced 2x2: all terms orthogonal, so the sums of squares add up to the total.
        bal = {r['source']: r for r in b.anova_table(fixture('2x2', [12, 12]), '2x2')}
        total = sum(bal[k]['ss'] for k in ('sequence', 'subject(sequence)', 'period', 'formulation', 'residual'))
        self.assertAlmostEqual(total, bal['total (corrected)']['ss'], 8)

    def test_formulation_line_matches_fixed_model(self):
        for design, counts in (('2x2', [11, 8]), ('partial', [7, 9, 6]), ('full', [10, 7])):
            rows = fixture(design, counts); at = {r['source']: r for r in b.anova_table(rows, design)}
            e, se, df, mse = b.fixed_fit(rows)
            self.assertAlmostEqual(at['residual']['ms'], mse, 12)
            self.assertEqual(at['residual']['df'], df)
            self.assertAlmostEqual(at['formulation']['f'], (e / se) ** 2, 8)

    def test_parallel_one_way(self):
        rng = np.random.default_rng(3)
        rows = [dict(subject=f'P{i}', treatment=tr, logvalue=rng.normal(4, .3), period='', sequence='') for i, tr in enumerate('TR' * 10)]
        at = {r['source']: r for r in b.anova_table(rows, 'parallel')}
        res = b.parallel_design(rows)
        self.assertEqual(at['residual']['df'], 18)
        self.assertAlmostEqual(at['formulation']['f'], (res.estimate / res.se) ** 2, 8)


class SfdaProfileTests(unittest.TestCase):
    def test_profile_uses_anova_rounding_and_table(self):
        rows = fixture('partial', [8, 8, 8])
        code, p, _ = run_cli(rows, '--design', 'partial', '--metric', 'cmax', '--profile', 'sfda')
        self.assertIn(code, (0, 1))
        self.assertEqual(p['scalars']['profile'], 'sfda')
        self.assertTrue(p['scalars']['method'].startswith('fixed subject'))
        self.assertIn('2 decimal', p['scalars']['ci_rounding'])
        self.assertEqual([r['source'] for r in table(p, 'anova')][:4], ['sequence', 'subject(sequence)', 'period', 'formulation'])
        ema = b.ema_analysis(rows, 'partial')
        self.assertAlmostEqual(p['scalars']['estimate'], ema.estimate, 12)

    def test_profile_rejects_non_sfda_methods(self):
        rows = fixture('partial', [8, 8, 8])
        for extra in (['--scaling', 'rsabe'], ['--scaling', 'both', '--abel-justified'], ['--analysis', 'contrast']):
            code, _, err = run_cli(rows, '--design', 'partial', '--metric', 'cmax', '--profile', 'sfda', *extra)
            self.assertEqual(code, 2); self.assertIn('SFDA profile', err)

    def test_abel_allowed_for_cmax_only(self):
        rows = fixture('full', [10, 10], cv=.6)
        code, p, _ = run_cli(rows, '--design', 'full', '--metric', 'cmax', '--profile', 'sfda', '--scaling', 'abel', '--abel-justified')
        self.assertIn(code, (0, 1)); self.assertEqual(len(table(p, 'scaled criteria')), 1)
        code, _, err = run_cli(rows, '--design', 'full', '--metric', 'auc', '--profile', 'sfda', '--scaling', 'abel', '--abel-justified')
        self.assertEqual(code, 2); self.assertIn('cmax', err)

    def test_minimum_18_subjects(self):
        small = fixture('2x2', [8, 8], cv=.05)
        code, p, _ = run_cli(small, '--design', '2x2', '--profile', 'sfda')
        self.assertEqual(code, 1); self.assertTrue(any('at least 18' in f for f in p['findings']))
        code, p, _ = run_cli(small, '--design', '2x2')  # no profile: no such rule
        self.assertEqual(code, 0)
        code, p, _ = run_cli(fixture('2x2', [9, 9], cv=.05), '--design', '2x2', '--profile', 'sfda')
        self.assertFalse(any('at least 18' in f for f in p['findings']))

    def test_predose_exclusion_2x2(self):
        rows = fixture('2x2', [12, 12])
        for r in rows: r['predose'] = 0.0
        bad = next(r for r in rows if r['period'] == '2'); bad['predose'] = 0.06 * bad['value']
        code, p, _ = run_cli(rows, '--design', '2x2', '--profile', 'sfda', extra=('predose',))
        self.assertEqual(p['scalars']['n_subjects'], 23)
        ex = table(p, 'pre-dose exclusions'); self.assertEqual(ex[0]['subject'], bad['subject'])
        code, p, _ = run_cli(rows, '--design', '2x2', '--profile', 'sfda', '--endogenous', '--baseline-corrected', extra=('predose',))
        self.assertEqual(p['scalars']['n_subjects'], 24)
        bad['predose'] = 0.05 * bad['value']  # exactly 5% is not "greater than"
        code, p, _ = run_cli(rows, '--design', '2x2', '--profile', 'sfda', extra=('predose',))
        self.assertEqual(p['scalars']['n_subjects'], 24)

    def test_predose_uses_cmax_column_for_auc(self):
        rows = fixture('2x2', [10, 10])
        for r in rows: r['predose'] = 1.0; r['cmax'] = 100.0
        rows[3]['cmax'] = 10.0  # 1.0 > 5% of 10
        code, p, _ = run_cli(rows, '--design', '2x2', '--metric', 'auc', '--profile', 'sfda', '--cmax-column', 'cmax', extra=('predose', 'cmax'))
        self.assertEqual(p['scalars']['n_subjects'], 19)

    def test_predose_in_replicate_removes_only_that_period(self):
        rows = fixture('partial', [7, 7, 7])
        for r in rows: r['predose'] = 0.0
        rows[4]['predose'] = rows[4]['value']
        code, p, err = run_cli(rows, '--design', 'partial', '--profile', 'sfda', extra=('predose',))
        self.assertIn(code, (0, 1), err)
        self.assertEqual(len(table(p, 'pre-dose exclusions')), 1)
        self.assertEqual(p['scalars']['n_subjects'], 21)                  # the subject keeps its other periods
        self.assertEqual(table(p, 'incomplete subjects')[0]['subject'], rows[4]['subject'])
        # the complete-data contrast analysis cannot use the remaining periods: explicit error
        code, _, err = run_cli(rows, '--design', 'partial', '--analysis', 'contrast', extra=('predose',))
        self.assertEqual(code, 2); self.assertIn('incomplete replicate', err)

    def test_power_minimum_and_anova_df(self):
        n, p = b.sample_size(.10, .97, .8, '2x2', analysis='ema', min_n=18)
        self.assertEqual(n, 18); self.assertGreater(p, .8)
        self.assertEqual(b.sample_size(.10, .97, .8, '2x2')[0] < 18, True)
        self.assertAlmostEqual(b.tost_power(24, .3, .95, 'partial', analysis='ema'), 0.7250, 3)
        with self.assertRaises(b.InputError):
            b.run(['--power', '--cv', '.2', '--n', '12', '--profile', 'sfda'])

    def test_power_cli_with_profile(self):
        out = io.StringIO()
        with redirect_stdout(out), redirect_stderr(io.StringIO()):
            b.run(['--power', '--cv', '.1', '--gmr', '.97', '--profile', 'sfda', '--format', 'json'])
        s = json.loads(out.getvalue())['scalars']
        self.assertEqual((s['n_total'], s['analysis'], s['profile']), (18, 'ema', 'sfda'))
        out = io.StringIO()
        with redirect_stdout(out), redirect_stderr(io.StringIO()):
            b.run(['--power', '--design', 'partial', '--cv', '.3', '--gmr', '.95', '--n', '24', '--profile', 'sfda', '--format', 'json'])
        self.assertAlmostEqual(json.loads(out.getvalue())['scalars']['power'], 0.7250, 3)


class CiLevelAndStageTests(unittest.TestCase):
    def test_ci_level_widens_interval(self):
        rows = fixture('2x2', [12, 12])
        r90 = b.ema_analysis(rows, '2x2'); r94 = b.ema_analysis(rows, '2x2', 94.12)
        self.assertLess(r94.ci_low, r90.ci_low); self.assertGreater(r94.ci_high, r90.ci_high)
        from scipy.stats import t
        self.assertAlmostEqual(math.log(r94.ci_high) - r94.estimate, t.ppf(1 - .0588 / 2, r94.df) * r94.se, 12)

    def test_stage_nests_periods(self):
        rows = fixture('2x2', [12, 12]); rng = np.random.default_rng(5)
        for r in rows:
            r['stage'] = '2' if r['subject'].endswith(('7', '8', '9', '10', '11')) else '1'
            if r['stage'] == '2' and r['period'] == '2': r['logvalue'] += .3; r['value'] = math.exp(r['logvalue'])
        e, se, df, _ = b.fixed_fit(rows)
        subjects = sorted({r['subject'] for r in rows}); keys = sorted({(r['stage'], r['period']) for r in rows})
        x = np.array([[float(r['subject'] == s) for s in subjects] + [float((r['stage'], r['period']) == k) for k in keys[1:]] + [float(r['treatment'] == 'T')] for r in rows])
        y = np.array([r['logvalue'] for r in rows]); beta = np.linalg.lstsq(x, y, rcond=None)[0]
        self.assertAlmostEqual(e, beta[-1], 12); self.assertEqual(df, len(y) - np.linalg.matrix_rank(x))
        code, p, _ = run_cli(rows, '--design', '2x2', '--profile', 'sfda', '--ci-level', '94.12', extra=('stage',))
        self.assertIn('period(stage)', [r['source'] for r in table(p, 'anova')])
        self.assertEqual(p['scalars']['ci_level'], 94.12)
        code, _, err = run_cli(rows, '--design', '2x2', extra=('stage',))
        self.assertEqual(code, 2); self.assertIn('stage or group column', err)

    def test_stage_terms_in_anova(self):
        rows = fixture('2x2', [12, 12])
        for r in rows: r['stage'] = '2' if int(r['subject'].split('-')[1]) >= 8 else '1'
        at = {r['source']: r for r in b.anova_table(rows, '2x2')}
        for term in ('stage', 'sequence', 'sequence×stage', 'subject(sequence×stage)', 'period(stage)', 'formulation'):
            self.assertIn(term, at)
        self.assertEqual((at['stage']['df'], at['sequence']['df'], at['sequence×stage']['df'], at['period(stage)']['df']), (1, 1, 1, 2))
        self.assertEqual(at['subject(sequence×stage)']['df'], 24 - 4)
        # stage/sequence terms are tested against subject(sequence x stage)
        self.assertAlmostEqual(at['stage']['f'], at['stage']['ms'] / at['subject(sequence×stage)']['ms'], 12)
        e, se, df, _ = b.fixed_fit(rows)
        self.assertAlmostEqual(at['formulation']['f'], (e / se) ** 2, 8)
        between = sum(at[k]['ss'] for k in ('stage', 'sequence', 'sequence×stage', 'subject(sequence×stage)'))
        within = at['period(stage)']['ss'] + at['formulation']['ss'] + at['residual']['ss']
        self.assertAlmostEqual(between + within, at['total (corrected)']['ss'], 8)  # balanced stages: orthogonal


class ExclusionFlagTests(unittest.TestCase):
    def test_low_reference_auc_is_flagged_not_excluded(self):
        rows = fixture('2x2', [10, 10])
        low = next(r for r in rows if r['treatment'] == 'R'); low['value'] = 1e-3; low['logvalue'] = math.log(1e-3)
        code, p, err = run_cli(rows, '--design', '2x2', '--metric', 'auc', '--profile', 'sfda')
        flagged = table(p, 'low AUC')
        self.assertEqual([f['subject'] for f in flagged], [low['subject']])
        self.assertEqual(p['scalars']['n_subjects'], 20)  # not excluded
        code, p, _ = run_cli(rows, '--design', '2x2', '--metric', 'cmax', '--profile', 'sfda')
        self.assertFalse(any(t['title'] == 'low AUC' for t in p['tables']))

    def test_endogenous_and_predose_notes(self):
        rows = fixture('2x2', [10, 10])
        for r in rows: r['predose'] = 0.0
        _, p, _ = run_cli(rows, '--design', '2x2', '--profile', 'sfda', extra=('predose',))
        self.assertTrue(any('single-dose' in n for n in p['notes']))
        _, p, _ = run_cli(rows, '--design', '2x2', '--profile', 'sfda', '--endogenous', '--baseline-corrected', extra=('predose',))
        self.assertTrue(any('already baseline-corrected' in n for n in p['notes']))


class IchM13aTests(unittest.TestCase):
    """Default (ICH M13A / EMA / FDA) framework: minimum 12 evaluable subjects,
    pre-dose > 5% of Cmax exclusion (single dose), low-AUC flag for T and R."""
    def test_minimum_12_crossover(self):
        code, p, _ = run_cli(fixture('2x2', [5, 5], cv=.05), '--design', '2x2')
        self.assertEqual(code, 1); self.assertTrue(any('at least 12' in f for f in p['findings']))
        code, p, _ = run_cli(fixture('2x2', [6, 6], cv=.05), '--design', '2x2')
        self.assertFalse(any('at least 12' in f for f in p['findings']))

    def test_minimum_12_per_arm_parallel(self):
        rng = np.random.default_rng(1)
        rows = [dict(subject=f'P{i}', treatment=tr, value=math.exp(rng.normal(4, .05))) for i, tr in enumerate('T' * 11 + 'R' * 14)]
        code, p, _ = run_cli(rows, '--design', 'parallel')
        self.assertEqual(code, 1); self.assertTrue(any('per arm' in f for f in p['findings']))

    def test_predose_rule_and_multiple_dose(self):
        rows = fixture('2x2', [12, 12])
        for r in rows: r['predose'] = 0.0
        bad = next(r for r in rows if r['period'] == '2'); bad['predose'] = 0.06 * bad['value']
        _, p, _ = run_cli(rows, '--design', '2x2', extra=('predose',))
        self.assertEqual(p['scalars']['n_subjects'], 23)
        self.assertTrue(any('ICH M13A' in n for n in p['notes']))
        _, p, _ = run_cli(rows, '--design', '2x2', '--multiple-dose', extra=('predose',))
        self.assertEqual(p['scalars']['n_subjects'], 24)
        _, p, _ = run_cli(rows, '--design', '2x2', '--profile', 'sfda', '--multiple-dose', extra=('predose',))
        self.assertEqual(p['scalars']['n_subjects'], 24)

    def test_low_auc_covers_test_under_ich_reference_only_under_sfda(self):
        rows = fixture('2x2', [10, 10])
        low = next(r for r in rows if r['treatment'] == 'T'); low['value'] = 1e-3; low['logvalue'] = math.log(1e-3)
        _, p, _ = run_cli(rows, '--design', '2x2', '--metric', 'auc')
        self.assertEqual([(f['subject'], f['treatment']) for f in table(p, 'low AUC')], [(low['subject'], 'T')])
        _, p, _ = run_cli(rows, '--design', '2x2', '--metric', 'auc', '--profile', 'sfda')
        self.assertFalse(any(t['title'] == 'low AUC' for t in p['tables']))

    def test_power_floor_12(self):
        out = io.StringIO()
        with redirect_stdout(out), redirect_stderr(io.StringIO()):
            b.run(['--power', '--cv', '.10', '--gmr', '.97', '--format', 'json'])
        self.assertEqual(json.loads(out.getvalue())['scalars']['n_total'], 12)
        with self.assertRaises(b.InputError): b.run(['--power', '--cv', '.2', '--n', '10'])
        out = io.StringIO()
        with redirect_stdout(out), redirect_stderr(io.StringIO()):
            b.run(['--power', '--design', 'parallel', '--cv', '.05', '--gmr', '.97', '--format', 'json'])
        self.assertEqual(json.loads(out.getvalue())['scalars']['n_total'], 24)


def endogenous_rows(added_t=35.0, added_r=50.0, baseline=100.0, n=12, seed=9):
    """2x2 study: every observation = endogenous baseline + dose-related increase."""
    rng = np.random.default_rng(seed); rows = []
    for seq in ('TR', 'RT'):
        for j in range(n):
            base = baseline * math.exp(rng.normal(0, .05))
            for period, tr in enumerate(seq, 1):
                add = (added_t if tr == 'T' else added_r) * math.exp(rng.normal(0, .08))
                rows.append(dict(subject=f'{seq}-{j}', sequence=seq, period=str(period), treatment=tr,
                                 value=base + add, baseline=base, baseline_auc=base * 24))
    return rows


class BaselineCorrectionTests(unittest.TestCase):
    def test_required_for_endogenous(self):
        rows = endogenous_rows()
        code, _, err = run_cli(rows, '--design', '2x2', '--metric', 'cmax', '--endogenous')
        self.assertEqual(code, 2); self.assertIn('baseline correction is required', err)
        # explicit column + confirmation is contradictory
        code, _, err = run_cli(rows, '--design', '2x2', '--metric', 'cmax', '--endogenous', '--baseline-corrected', '--baseline-column', 'baseline', extra=('baseline',))
        self.assertEqual(code, 2); self.assertIn('choose one', err)
        # confirmation while the file merely contains a default-named baseline column: runs, warns, does not subtract
        code, p, _ = run_cli(rows, '--design', '2x2', '--metric', 'cmax', '--endogenous', '--baseline-corrected', extra=('baseline',))
        self.assertEqual(code, 0); self.assertAlmostEqual(p['scalars']['gmr_pct'], 90, delta=2)
        self.assertTrue(any(n.startswith('WARNING') and 'NOT subtracted' in n for n in p['notes']))
        code, _, err = run_cli(rows, '--design', '2x2', '--metric', 'cmax', '--baseline-corrected')
        self.assertEqual(code, 2)

    def test_uncorrected_hides_difference_corrected_reveals_it(self):
        rows = endogenous_rows()
        # Without correction the ratio is ~(100+35)/(100+50) = 90%: looks bioequivalent.
        code, p, _ = run_cli(rows, '--design', '2x2', '--metric', 'cmax', '--endogenous', '--baseline-corrected')
        self.assertEqual(code, 0); self.assertAlmostEqual(p['scalars']['gmr_pct'], 90, delta=2)
        # With subtraction the true ratio 35/50 = 70% appears and BE fails.
        code, p, _ = run_cli(rows, '--design', '2x2', '--metric', 'cmax', '--endogenous', extra=('baseline',))
        self.assertEqual(code, 1); self.assertAlmostEqual(p['scalars']['gmr_pct'], 70, delta=3)
        self.assertIn('pre-dose concentration', p['scalars']['baseline_correction'])
        corr = table(p, 'baseline correction')
        self.assertEqual(len(corr), 48)
        r0, s0 = corr[0], rows[0]
        self.assertAlmostEqual(r0['corrected'], s0['value'] - s0['baseline'], 9)

    def test_corrected_values_enter_the_model(self):
        rows = endogenous_rows(added_t=40, added_r=40)
        code, p, _ = run_cli(rows, '--design', '2x2', '--metric', 'cmax', '--endogenous', '--analysis', 'ema', extra=('baseline',))
        manual = [dict(r, value=r['value'] - r['baseline'], logvalue=math.log(r['value'] - r['baseline'])) for r in rows]
        self.assertAlmostEqual(p['scalars']['estimate'], b.ema_analysis(manual, '2x2').estimate, 12)

    def test_auc_correction_options(self):
        rows = endogenous_rows()
        for r in rows: r['value'] = r['value'] * 24  # AUC over 24 h
        code, p, _ = run_cli(rows, '--design', '2x2', '--metric', 'auc', '--endogenous', extra=('baseline_auc',))
        self.assertIn('baseline) AUC', p['scalars']['baseline_correction'])
        code2, p2, _ = run_cli(rows, '--design', '2x2', '--metric', 'auc', '--endogenous', '--auc-hours', '24', extra=('baseline',))
        self.assertAlmostEqual(p['scalars']['estimate'], p2['scalars']['estimate'], 10)
        code, _, err = run_cli(rows, '--design', '2x2', '--metric', 'auc', '--endogenous', extra=('baseline',))
        self.assertEqual(code, 2); self.assertIn('--auc-hours', err)

    def test_nonpositive_after_correction_is_explicit_error(self):
        rows = endogenous_rows()
        rows[5]['baseline'] = rows[5]['value'] + 1
        code, _, err = run_cli(rows, '--design', '2x2', '--metric', 'cmax', '--endogenous', extra=('baseline',))
        self.assertEqual(code, 2); self.assertIn('zero or negative', err); self.assertIn(rows[5]['subject'], err)

    def test_sfda_profile_with_endogenous_skips_predose_and_corrects(self):
        rows = endogenous_rows()
        for r in rows: r['predose'] = r['baseline']  # endogenous pre-dose levels are naturally high
        code, p, _ = run_cli(rows, '--design', '2x2', '--metric', 'cmax', '--profile', 'sfda', '--endogenous', extra=('baseline', 'predose'))
        self.assertEqual(p['scalars']['n_subjects'], 24)
        self.assertFalse(any(t['title'] == 'pre-dose exclusions' for t in p['tables']))
        self.assertTrue(any(t['title'] == 'baseline correction' for t in p['tables']))


if __name__ == '__main__':
    unittest.main(verbosity=2)
