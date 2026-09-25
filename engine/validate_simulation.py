"""Seeded raw-log-data simulation; verifies the actual analysis functions.
No attempt to validate regulatory acceptance or RSABE type-I error globally.
"""
import json, math, platform
import numpy as np
import scipy
from scipy.stats import t
import bioequivalence as b

SEED=20260925
REPETITIONS=4000
rng=np.random.default_rng(SEED)
results=[]
for design,n in [('2x2',24),('partial',24),('full',24)]:
    seqs=b.DESIGNS[design]; sigma=math.sqrt(math.log1p(.3**2)); delta=math.log(.95)
    covered=passed=0
    for iteration in range(REPETITIONS):
        rows=[]
        for seq in seqs:
            for j in range(n//len(seqs)):
                intercept=rng.normal(4,.5)
                for period,tr in enumerate(seq,1):
                    y=intercept+.08*period+(delta if tr=='T' else 0)+rng.normal(0,sigma)
                    rows.append(dict(subject=f'{seq}-{j}',sequence=seq,period=str(period),treatment=tr,logvalue=y))
        result=b.contrast_analysis(rows,design)
        covered+=result.ci_low<=.95<=result.ci_high
        passed+=b.passes_ci(result)
    theoretical=b.tost_power(n,.3,.95,design)
    empirical=passed/REPETITIONS; coverage=covered/REPETITIONS
    mcse=math.sqrt(theoretical*(1-theoretical)/REPETITIONS)
    assert abs(empirical-theoretical)<4*mcse
    assert abs(coverage-.9)<4*math.sqrt(.9*.1/REPETITIONS)
    results.append(dict(design=design,n=n,theoretical_power=theoretical,empirical_power=empirical,power_mc_standard_error=mcse,ci90_coverage=coverage))
output=dict(seed=SEED,repetitions_per_design=REPETITIONS,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,assumptions=dict(true_gmr=.95,cv_t=.3,cv_r=.3,between_subject_log_sd=.5,period_effect='.08 * period',interaction=0,missing_data=False,allocation='balanced',analysis='equal-sequence subject contrasts'),results=results)
print(json.dumps(output,indent=2))
