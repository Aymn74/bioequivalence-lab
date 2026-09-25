import copy
import io
import math
import unittest
from contextlib import redirect_stdout, redirect_stderr
import numpy as np
from scipy.stats import chi2,t
import bioequivalence as b


def fixture(design='partial',counts=None,seed=72819,cv=.4):
    rng=np.random.default_rng(seed); rows=[]
    sequences=b.DESIGNS[design]; counts=counts or [12]*len(sequences)
    for seq,n in zip(sequences,counts):
        for j in range(n):
            subj=f'{seq}-{j}'; intercept=rng.normal(4,.5)
            for period,tr in enumerate(seq,1):
                y=intercept+.13*period+(math.log(.95) if tr=='T' else 0)+rng.normal(0,math.sqrt(math.log1p(cv**2)))
                rows.append(dict(subject=subj,sequence=seq,period=str(period),treatment=tr,logvalue=y,value=math.exp(y)))
    return rows

class Tests(unittest.TestCase):
    def test_2x2_against_independent_full_ols(self):
        rows=fixture('2x2',[14,9]); result=b.contrast_analysis(rows)
        subjects=sorted({r['subject'] for r in rows})
        x=np.array([[float(r['subject']==s) for s in subjects]+[float(r['period']=='2'),float(r['treatment']=='T')] for r in rows])
        y=np.array([r['logvalue'] for r in rows]); beta=np.linalg.lstsq(x,y,rcond=None)[0]
        df=len(y)-np.linalg.matrix_rank(x); mse=np.sum((y-x@beta)**2)/df
        self.assertAlmostEqual(result.estimate,beta[-1],12)
        self.assertAlmostEqual(result.se,math.sqrt(mse*np.linalg.inv(x.T@x)[-1,-1]),12)
        self.assertEqual(result.df,df)
    def test_three_sequence_against_glm(self):
        rows=fixture(counts=[7,11,9]); result=b.contrast_analysis(rows,'partial')
        _,groups=b.validate(rows,'partial'); seqs=b.DESIGNS['partial']
        x=np.array([[float(g[0]['sequence']==s) for s in seqs] for g in groups.values()])
        y=np.array([np.mean([r['logvalue'] for r in g if r['treatment']=='T'])-np.mean([r['logvalue'] for r in g if r['treatment']=='R']) for g in groups.values()])
        beta=np.linalg.lstsq(x,y,rcond=None)[0]; c=np.ones(3)/3
        mse=np.sum((y-x@beta)**2)/(len(y)-3)
        self.assertAlmostEqual(result.estimate,c@beta,12)
        self.assertAlmostEqual(result.se,math.sqrt(mse*(c@np.linalg.inv(x.T@x)@c)),12)
    def test_reference_period_invariance(self):
        for design in ('partial','full'):
            for method in ('contrast','ema'):
                rows=fixture(design); rv=b.reference_variability(rows,method,design)
                changed=copy.deepcopy(rows)
                for r in changed:r['logvalue']+=int(r['period'])**2*.7
                other=b.reference_variability(changed,method,design)
                self.assertAlmostEqual(rv.s2wr,other.s2wr,12)
                self.assertEqual(rv.df,other.df)
    def test_reference_df(self):
        rows=fixture('partial')
        self.assertEqual(b.reference_variability(rows,'contrast').df,33)
        self.assertEqual(b.reference_variability(rows,'ema').df,34)
    def test_rsabe_fda_reference_algebra(self):
        # Literal independent translation of FDA 2011 SAS, pages 5-6.
        for delta,se,df,s2,dfw in [(.1,.05,21,.13,21),(-.12,.07,33,.2,33),(0,.04,15,.1,15)]:
            x=delta*delta-se*se; low=delta-t.ppf(.95,df)*se; high=delta+t.ppf(.95,df)*se
            boundx=max(abs(low),abs(high))**2
            theta=(np.log(1.25)/.25)**2; y=-theta*s2; boundy=y*dfw/chi2.ppf(.95,dfw)
            expected=x+y+np.sqrt((boundx-x)**2+(boundy-y)**2)
            actual=b.rsabe_bound(delta,se,df,b.ReferenceVariability(s2,dfw,dfw+3))
            self.assertAlmostEqual(actual['criterion_95_upper_bound'],expected,14)
    def test_malformed_data(self):
        base=fixture()
        mutations=[]
        for key,val in [('sequence',''),('sequence','?'),('period','9'),('period','1.5'),('treatment','R'),('subject',''),('logvalue',float('nan'))]:
            bad=copy.deepcopy(base);bad[0][key]=val;mutations.append(bad)
        mutations.extend([base[:-1],base+[base[0]]])
        for bad in mutations:
            with self.assertRaises(b.InputError):b.contrast_analysis(bad,'partial')
    def test_power_known_and_designs(self):
        self.assertAlmostEqual(b.tost_power(24,.3,.95,'partial'),.708053175348,10)
        self.assertGreater(b.tost_power(24,.3,.95,'full'),b.tost_power(24,.3,.95,'partial'))
        self.assertEqual(b.power_parameters(24,'partial')[0],21)
        self.assertEqual(b.power_parameters(24,'full')[0],22)
        self.assertEqual(b.power_parameters(24,'partial','ema')[0],45)
        for n,d in [(25,'partial'),(25,'2x2'),(24,'replicate')]:
            with self.assertRaises(b.InputError):b.tost_power(n,.3,.95,d)
    def test_sample_size_minimal(self):
        for d in ('2x2','partial','full'):
            n,p=b.sample_size(.3,.95,.8,d);m=len(b.DESIGNS[d])
            self.assertGreaterEqual(p,.8)
            self.assertLess(b.tost_power(n-m,.3,.95,d),.8)
            self.assertEqual(n%m,0)
    def test_rounding_and_cap(self):
        self.assertEqual(str(b.rounded_pct(.79996)),'80.00')
        rv=b.ReferenceVariability(math.log1p(.8**2),20,23)
        low,high,_=b.abel_limits(rv)
        self.assertAlmostEqual(low,.6984,4);self.assertAlmostEqual(high,1.4319,4)
    def test_fixed_model_against_lstsq(self):
        for d in ('partial','full'):
            rows=fixture(d); e,se,df,mse=b.fixed_fit(rows)
            self.assertGreater(se,0);self.assertEqual(df,(len({r['subject'] for r in rows})-1)*(len(b.DESIGNS[d][0])-1)-1)
            shifted=copy.deepcopy(rows)
            for r in shifted:r['logvalue']+=.7 if r['treatment']=='T' else 0
            e2,se2,_,_=b.fixed_fit(shifted)
            self.assertAlmostEqual(e2-e,.7,12);self.assertAlmostEqual(se,se2,12)
    def test_parallel_welch(self):
        from scipy.stats import ttest_ind
        rows=[dict(subject=str(i),treatment='T' if i<5 else 'R',logvalue=y) for i,y in enumerate([1,2,3,4,7,2,2.1,2.2,2.4,2.5,2.6])]
        result=b.parallel_design(rows,True)
        ref=ttest_ind([r['logvalue'] for r in rows[:5]],[r['logvalue'] for r in rows[5:]],equal_var=False)
        self.assertAlmostEqual(result.estimate/result.se,ref.statistic,12);self.assertAlmostEqual(result.df,ref.df,12)
    def test_cli_guards(self):
        for argv in [['--power','--cv','.3','--design','replicate'],['--power','--cv','.3','--scaling','rsabe','--design','partial'],['--power','--cv','.3','--nti','--limits','.7,1.4']]:
            with self.assertRaises(b.InputError):b.run(argv)

if __name__=='__main__':unittest.main(verbosity=2)
