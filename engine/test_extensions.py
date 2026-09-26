"""FDA NTI method, multiple-treatment (e.g. Williams) designs and multi-group studies.

Each result is checked against an independent computation: the FDA NTI steps are
re-implemented literally from the SAS example in Appendix F of FDA's "Statistical
Approaches to Establishing Bioequivalence" (May 2026); the fixed-effects fits are
compared with ordinary least squares in NumPy."""
import math
import unittest
import numpy as np
from scipy.stats import chi2, f as f_dist, t
import bioequivalence as b
from test_sfda import run_cli, table


def full_replicate(n=(12, 12), gmr=.98, swr=.08, swt=.08, seed=4):
    rng = np.random.default_rng(seed); rows = []
    for seq, k in zip(('TRTR', 'RTRT'), n):
        for j in range(k):
            s = rng.normal(4, .3)
            for per, tr in enumerate(seq, 1):
                sd = swt if tr == 'T' else swr
                y = s + .02 * per + (math.log(gmr) if tr == 'T' else 0) + rng.normal(0, sd)
                rows.append(dict(subject=f'{seq}-{j}', sequence=seq, period=str(per), treatment=tr, logvalue=y, value=math.exp(y)))
    return rows


def sas_appendix_f(rows):
    """Literal re-implementation of the Appendix F SAS steps (ilat, dlat, tlat)."""
    subj = {}
    for r in rows: subj.setdefault(r['subject'], {'seq': r['sequence'], 'T': [], 'R': []})[r['treatment']].append((int(r['period']), r['logvalue']))
    by = {}
    for v in subj.values():
        T = [x for _, x in sorted(v['T'])]; R = [x for _, x in sorted(v['R'])]
        by.setdefault(v['seq'], []).append((.5 * (T[0] + T[1] - R[0] - R[1]), R[0] - R[1], T[0] - T[1]))
    seqs = sorted(by); m = len(seqs); N = sum(len(by[s]) for s in seqs); df = N - m
    def seq_model(idx):
        means = [np.mean([x[idx] for x in by[s]]) for s in seqs]
        mse = sum(sum((x[idx] - mu) ** 2 for x in by[s]) for s, mu in zip(seqs, means)) / df
        return np.mean(means), math.sqrt(mse * sum(1 / len(by[s]) for s in seqs) / m ** 2), mse
    est, se, _ = seq_model(0)
    lower, upper = est - t.ppf(.95, df) * se, est + t.ppf(.95, df) * se
    x = est ** 2 - se ** 2; boundx = max(abs(lower), abs(upper)) ** 2
    s2wr = seq_model(1)[2] / 2; s2wt = seq_model(2)[2] / 2
    theta = (math.log(1 / 0.9) / 0.1) ** 2
    y = -theta * s2wr; boundy = y * df / chi2.ppf(.95, df)
    crit = (x + y) + math.sqrt((boundx - x) ** 2 + (boundy - y) ** 2)
    ratio = math.sqrt(s2wt / s2wr)
    return dict(estimate=est, se=se, crit=crit, s2wr=s2wr, s2wt=s2wt, ratio_hi=ratio / math.sqrt(f_dist.ppf(.05, df, df)), ratio_lo=ratio / math.sqrt(f_dist.ppf(.95, df, df)))


class FdaNtiTests(unittest.TestCase):
    def test_theta_constant(self):
        self.assertAlmostEqual(b.THETA_FDA_NTI, 1.1100, 3)

    def test_matches_appendix_f_steps(self):
        for seed, n in ((4, (12, 12)), (9, (14, 10)), (21, (9, 13))):
            rows = full_replicate(n, seed=seed); ref = sas_appendix_f(rows); got = b.fda_nti(rows, 'full')
            self.assertAlmostEqual(got['estimate'], ref['estimate'], 12)
            self.assertAlmostEqual(got['se'], ref['se'], 12)
            self.assertAlmostEqual(got['criterion_95_upper_bound'], ref['crit'], 12)
            self.assertAlmostEqual(got['s2wr'], ref['s2wr'], 12); self.assertAlmostEqual(got['s2wt'], ref['s2wt'], 12)
            self.assertAlmostEqual(got['swt_swr_ci90_high'], ref['ratio_hi'], 12)
            self.assertAlmostEqual(got['swt_swr_ci90_low'], ref['ratio_lo'], 12)

    def test_three_conditions(self):
        good = b.fda_nti(full_replicate(gmr=.99, swr=.08, swt=.08), 'full')
        self.assertTrue(good['met'] and good['passes_scaled_criterion'] and good['passes_unscaled_abe'] and good['passes_variability_ratio'])
        # test product far more variable than the reference: the sWT/sWR condition fails
        noisy = b.fda_nti(full_replicate(gmr=.99, swr=.06, swt=.30, seed=5), 'full')
        self.assertFalse(noisy['passes_variability_ratio']); self.assertFalse(noisy['met'])
        # 8% difference with low variability: fixed 80-125 passes, the tight scaled criterion fails
        shifted = b.fda_nti(full_replicate(gmr=.92, swr=.05, swt=.05, seed=6), 'full')
        self.assertTrue(shifted['passes_unscaled_abe']); self.assertFalse(shifted['passes_scaled_criterion']); self.assertFalse(shifted['met'])

    def test_cli_and_guards(self):
        rows = full_replicate()
        code, p, _ = run_cli(rows, '--design', 'full', '--metric', 'auc', '--scaling', 'fda-nti')
        self.assertEqual(code, 0); self.assertEqual(table(p, 'scaled criteria')[0]['criterion'], 'FDA NTI (reference-scaled, full replicate)')
        code, _, err = run_cli(rows, '--design', 'full', '--scaling', 'fda-nti', '--nti')
        self.assertEqual(code, 2)
        code, _, err = run_cli(rows, '--design', 'full', '--scaling', 'fda-nti', '--profile', 'sfda')
        self.assertEqual(code, 2); self.assertIn('not an SFDA method', err)
        code, _, err = run_cli(fixture_partial(), '--design', 'partial', '--scaling', 'fda-nti')
        self.assertEqual(code, 2); self.assertIn('full replicate', err)


def fixture_partial():
    from test_bioequivalence import fixture
    return fixture('partial', [8, 8, 8])


WILLIAMS = ('ABC', 'BCA', 'CAB', 'ACB', 'BAC', 'CBA')

def williams(n=4, seed=3, effects=dict(A=0, B=.03, C=-.02)):
    rng = np.random.default_rng(seed); rows = []
    for q in WILLIAMS:
        for j in range(n):
            s = rng.normal(4, .3)
            for per, tr in enumerate(q, 1):
                y = s + .04 * per + effects[tr] + rng.normal(0, .15)
                rows.append(dict(subject=f'{q}-{j}', sequence=q, period=str(per), treatment=tr, value=math.exp(y)))
    return rows


def ols_formulation(rows, test, ref, period_key=lambda r: r['period']):
    rows = [r for r in rows if r['treatment'] in (test, ref)]
    subs = sorted({r['subject'] for r in rows}); keys = sorted({period_key(r) for r in rows})
    x = np.array([[float(r['subject'] == s) for s in subs] + [float(period_key(r) == k) for k in keys[1:]] + [float(r['treatment'] == test)] for r in rows])
    y = np.array([math.log(r['value']) for r in rows]); beta = np.linalg.lstsq(x, y, rcond=None)[0]
    rk = np.linalg.matrix_rank(x); df = len(y) - rk; mse = np.sum((y - x @ beta) ** 2) / df
    se = math.sqrt(mse * np.linalg.pinv(x.T @ x)[-1, -1])
    return beta[-1], se, df


class MultiTreatmentTests(unittest.TestCase):
    def test_williams_pairs_match_ols(self):
        rows = williams()
        for test, ref in (('A', 'B'), ('A', 'C'), ('C', 'B')):
            code, p, _ = run_cli(rows, '--design', 'multi', '--test-label', test, '--reference-label', ref)
            est, se, df = ols_formulation(rows, test, ref)
            self.assertAlmostEqual(p['scalars']['estimate'], est, 10)
            self.assertAlmostEqual(p['scalars']['se'], se, 10)
            self.assertEqual(p['scalars']['df'], df)
            self.assertEqual(p['scalars']['n_subjects'], 24)

    def test_anova_and_sfda(self):
        rows = williams()
        code, p, _ = run_cli(rows, '--design', 'multi', '--test-label', 'A', '--reference-label', 'B', '--profile', 'sfda')
        at = {r['source']: r for r in table(p, 'anova')}
        self.assertEqual(at['sequence']['df'], 5); self.assertEqual(at['formulation']['df'], 1)
        e, se, _ = ols_formulation(rows, 'A', 'B')
        self.assertAlmostEqual(at['formulation']['f'], (e / se) ** 2, 8)
        self.assertTrue(any('comparison A (test) vs B' in n for n in p['notes']))

    def test_guards(self):
        rows = williams()
        code, _, err = run_cli(rows, '--design', 'multi', '--test-label', 'A', '--reference-label', 'Z')
        self.assertEqual(code, 2); self.assertIn('"Z" not found', err)
        code, _, err = run_cli(rows, '--design', 'multi', '--test-label', 'A', '--reference-label', 'B', '--analysis', 'contrast')
        self.assertEqual(code, 2)
        code, _, err = run_cli(rows, '--design', 'multi', '--test-label', 'A', '--reference-label', 'B', '--scaling', 'abel', '--abel-justified', '--metric', 'cmax')
        self.assertEqual(code, 2)

    def test_predose_drops_period_then_incomplete_subject(self):
        rows = williams()
        for r in rows: r['predose'] = 0.0
        bad = next(r for r in rows if r['treatment'] == 'B' and r['period'] != '1'); bad['predose'] = bad['value']
        code, p, _ = run_cli(rows, '--design', 'multi', '--test-label', 'A', '--reference-label', 'B', extra=('predose',))
        self.assertEqual(p['scalars']['n_subjects'], 23)
        # a flagged period of the excluded treatment C does not affect the A-B comparison
        for r in rows: r['predose'] = 0.0
        other = next(r for r in rows if r['treatment'] == 'C' and r['period'] != '1'); other['predose'] = other['value']
        code, p, _ = run_cli(rows, '--design', 'multi', '--test-label', 'A', '--reference-label', 'B', extra=('predose',))
        self.assertEqual(p['scalars']['n_subjects'], 24)


class MultiGroupTests(unittest.TestCase):
    def rows(self):
        from test_bioequivalence import fixture
        rows = fixture('2x2', [15, 15])
        for r in rows:
            j = int(r['subject'].split('-')[1]); r['group'] = 'G1' if j < 6 else 'G2' if j < 11 else 'G3'
            if r['group'] == 'G2' and r['period'] == '2': r['logvalue'] += .2; r['value'] = math.exp(r['logvalue'])
        return rows

    def test_model_and_terms(self):
        rows = self.rows()
        code, p, _ = run_cli(rows, '--design', '2x2', '--analysis', 'ema', '--anova', extra=('group',))
        est, se, df = ols_formulation(rows, 'T', 'R', period_key=lambda r: (r['group'], r['period']))
        self.assertAlmostEqual(p['scalars']['estimate'], est, 10); self.assertAlmostEqual(p['scalars']['se'], se, 10)
        at = {r['source']: r for r in table(p, 'anova')}
        for term in ('group', 'sequence', 'sequence×group', 'subject(sequence×group)', 'period(group)', 'formulation'):
            self.assertIn(term, at)
        self.assertEqual((at['group']['df'], at['period(group)']['df']), (2, 3))

    def test_supportive_interaction(self):
        rows = self.rows()
        code, p, _ = run_cli(rows, '--design', '2x2', '--analysis', 'ema', extra=('group',))
        gi = table(p, 'group interaction')[0]
        # independent: compare residual SS with and without group x formulation columns
        subs = sorted({r['subject'] for r in rows}); keys = sorted({(r['group'], r['period']) for r in rows}); groups = sorted({r['group'] for r in rows})
        def X(inter):
            out = []
            for r in rows:
                t_ = float(r['treatment'] == 'T')
                out.append([float(r['subject'] == s) for s in subs] + [float((r['group'], r['period']) == k) for k in keys[1:]] + [t_] + ([float(r['group'] == g) * t_ for g in groups[1:]] if inter else []))
            return np.array(out)
        y = np.array([r['logvalue'] for r in rows])
        def rss(x): beta = np.linalg.lstsq(x, y, rcond=None)[0]; return np.sum((y - x @ beta) ** 2), np.linalg.matrix_rank(x)
        (r0, k0), (r1, k1) = rss(X(False)), rss(X(True))
        F = ((r0 - r1) / (k1 - k0)) / (r1 / (len(y) - k1))
        self.assertEqual(gi['df'], 2); self.assertAlmostEqual(gi['f'], F, 8)

    def test_group_requires_fixed_model(self):
        code, _, err = run_cli(self.rows(), '--design', '2x2', extra=('group',))
        self.assertEqual(code, 2); self.assertIn('fixed-effects model', err)


if __name__ == '__main__':
    unittest.main(verbosity=2)
