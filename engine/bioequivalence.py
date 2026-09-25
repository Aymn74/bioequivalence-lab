#!/usr/bin/env python3
"""Auditable BE calculations for complete, canonical crossover designs.
See README.md for model assumptions, regulatory scope and intentional API changes.
"""
from __future__ import annotations
import argparse
import math
from collections import defaultdict
from dataclasses import dataclass, asdict
from decimal import Decimal, ROUND_HALF_UP
import numpy as np
from scipy.integrate import quad
from scipy.linalg import qr as pivoted_qr
from scipy.stats import chi2, norm, t
from _common import InputError, Report, read_table, require_columns, add_format_argument, main_wrapper

DESIGNS = {'2x2': ('TR', 'RT'), 'partial': ('TRR', 'RTR', 'RRT'), 'full': ('TRTR', 'RTRT')}
THETA_FDA = (math.log(1.25)/0.25)**2

@dataclass
class AverageBE:
    estimate: float
    se: float
    df: float
    n_subjects: int
    method: str
    @property
    def gmr(self): return math.exp(self.estimate)
    @property
    def ci_low(self): return math.exp(self.estimate-t.ppf(.95,self.df)*self.se)
    @property
    def ci_high(self): return math.exp(self.estimate+t.ppf(.95,self.df)*self.se)

@dataclass
class ReferenceVariability:
    s2wr: float
    df: int
    n_subjects: int
    @property
    def swr(self): return math.sqrt(self.s2wr)
    @property
    def cvwr(self): return math.sqrt(math.expm1(self.s2wr))

def positive(value, label):
    try: x = float(value)
    except (ValueError, TypeError): raise InputError(f'{label}: expected a finite positive number')
    if not math.isfinite(x) or x <= 0: raise InputError(f'{label}: expected a finite positive number')
    return x

def validate(records, design):
    """No silent exclusion, imputation, sequence guessing or duplicate averaging."""
    groups = defaultdict(list)
    for r in records:
        if not r.get('subject'): raise InputError('empty subject ID')
        if r.get('treatment') not in ('T','R'): raise InputError('treatment must be T or R')
        if not math.isfinite(r['logvalue']): raise InputError('non-finite log value')
        groups[r['subject']].append(r)
    if design == 'parallel':
        if any(len(g)!=1 for g in groups.values()): raise InputError('parallel: one observation per unique subject required')
        if any(sum(r['treatment']==tr for r in records)<2 for tr in ('T','R')): raise InputError('parallel: at least two subjects per arm')
        return design, groups
    seqs = {r.get('sequence','') for r in records}
    if design == 'replicate':
        matches = [d for d in ('partial','full') if seqs == set(DESIGNS[d])]
        if not matches: raise InputError('replicate requires all TRR/RTR/RRT or TRTR/RTRT sequences')
        design = matches[0]
    if design not in DESIGNS or seqs != set(DESIGNS[design]):
        raise InputError(f'{design}: all canonical sequences required; missing/unknown sequence found')
    for sid,g in groups.items():
        seq = g[0]['sequence']
        if any(r['sequence']!=seq for r in g): raise InputError(f'{sid}: inconsistent sequence')
        if len(g)!=len(seq): raise InputError(f'{sid}: incomplete or duplicate records')
        try: periods = [int(str(r['period'])) for r in g]
        except (ValueError,KeyError): raise InputError(f'{sid}: integer period required')
        if sorted(periods)!=list(range(1,len(seq)+1)): raise InputError(f'{sid}: missing/duplicate/invalid period')
        for r,p in zip(g,periods):
            if seq[p-1]!=r['treatment']: raise InputError(f'{sid}: treatment conflicts with sequence/period')
        g.sort(key=lambda r:int(r['period']))
    counts = [sum(g[0]['sequence']==s for g in groups.values()) for s in DESIGNS[design]]
    if min(counts)<2: raise InputError('at least two complete subjects in every sequence required')
    return design, groups

def contrast_analysis(records, design='2x2'):
    design,groups = validate(records,design)
    ds = defaultdict(list)
    for g in groups.values():
        diff = np.mean([r['logvalue'] for r in g if r['treatment']=='T'])-np.mean([r['logvalue'] for r in g if r['treatment']=='R'])
        ds[g[0]['sequence']].append(diff)
    m=len(ds); df=len(groups)-m
    mse=sum(sum((np.asarray(v)-np.mean(v))**2) for v in ds.values())/df
    se=math.sqrt(mse*sum(1/len(v) for v in ds.values())/m**2)
    return AverageBE(float(np.mean([np.mean(v) for v in ds.values()])),se,df,len(groups),'equal-sequence subject contrasts; common contrast variance')

# Compatibility name, with stricter validation and corrected m-squared denominator.
crossover_2x2 = contrast_analysis

def fixed_fit(records, treatment=True):
    """OLS subject + period (+ treatment); sequence is absorbed by subject effects.
    Reduced QR avoids forming the inverse of X'X. Complete designs only.
    """
    subjects=sorted({r['subject'] for r in records}); periods=sorted({int(r['period']) for r in records})
    x=np.array([[1.]+[float(r['subject']==s) for s in subjects[1:]]+[float(int(r['period'])==p) for p in periods[1:]]+([float(r['treatment']=='T')] if treatment else []) for r in records])
    y=np.array([r['logvalue'] for r in records])
    rank=np.linalg.matrix_rank(x)
    if rank!=x.shape[1]:
        if treatment: raise InputError('rank-deficient treatment model')
        # Reference-only full replication has disconnected period/subject blocks.
        # Remove aliased nuisance columns; retain their complete column space.
        _,_,pivot=pivoted_qr(x,mode='economic',pivoting=True)
        x=x[:,pivot[:rank]]
    df=len(y)-x.shape[1]
    if df<=0: raise InputError('no residual degrees of freedom')
    q,r=np.linalg.qr(x,mode='reduced'); beta=np.linalg.solve(r,q.T@y)
    mse=float(np.sum((y-x@beta)**2)/df)
    z=np.linalg.solve(r.T,np.eye(len(beta))[:,-1])
    return float(beta[-1]), math.sqrt(mse*float(z@z)), df, mse

def ema_analysis(records, design):
    _,groups=validate(records,design)
    e,se,df,_=fixed_fit(records)
    return AverageBE(e,se,df,len(groups),'fixed subject + period + treatment; common residual variance (EMA Method A structure)')

def reference_variability(records, method='contrast', design='replicate'):
    design,groups=validate(records,design)
    if design not in ('partial','full'): raise InputError('replicated reference required')
    if method=='ema':
        _,_,df,mse=fixed_fit([r for r in records if r['treatment']=='R'],False)
        return ReferenceVariability(mse,df,len(groups))
    ds=defaultdict(list)
    for g in groups.values():
        vals=[r['logvalue'] for r in g if r['treatment']=='R']
        ds[g[0]['sequence']].append(vals[0]-vals[1])
    df=len(groups)-len(ds)
    s2=sum(sum((np.asarray(v)-np.mean(v))**2) for v in ds.values())/(2*df)
    return ReferenceVariability(float(s2),df,len(groups))

def rsabe_bound(estimate,se,df_point,rv):
    if not all(math.isfinite(x) for x in (estimate,se,df_point,rv.s2wr,rv.df)) or se<0 or rv.s2wr<0 or min(df_point,rv.df)<=0:
        raise InputError('invalid RSABE inputs')
    x=estimate**2-se**2
    bx=(abs(estimate)+t.ppf(.95,df_point)*se)**2
    y=-THETA_FDA*rv.s2wr
    by=y*rv.df/chi2.ppf(.95,rv.df)
    bound=float(x+y+math.hypot(bx-x,by-y))
    return {'criterion_point_estimate':x+y,'criterion_95_upper_bound':bound,'passes_scaled_criterion':bool(bound<=0)}

def abel_limits(rv):
    if rv.cvwr<=.30: return .8,1.25,False
    s=min(rv.swr,math.sqrt(math.log1p(.5**2)))
    return math.exp(-.760*s),math.exp(.760*s),True

def rounded_pct(x): return Decimal(str(100*x)).quantize(Decimal('.01'),rounding=ROUND_HALF_UP)

def passes_ci(result,limits=(.8,1.25),rounding=False):
    if rounding: return rounded_pct(result.ci_low)>=rounded_pct(limits[0]) and rounded_pct(result.ci_high)<=rounded_pct(limits[1])
    return bool(result.ci_low>=limits[0] and result.ci_high<=limits[1])

def parallel_design(records,welch=False):
    _,groups=validate(records,'parallel')
    a,b=[np.array([r['logvalue'] for r in records if r['treatment']==tr]) for tr in ('T','R')]
    n1,n2=len(a),len(b); v1,v2=a.var(ddof=1),b.var(ddof=1)
    if welch:
        v=v1/n1+v2/n2
        if v<=0: raise InputError('Welch requires positive estimated variance')
        df=v**2/((v1/n1)**2/(n1-1)+(v2/n2)**2/(n2-1))
    else:
        df=n1+n2-2; v=((n1-1)*v1+(n2-1)*v2)/df*(1/n1+1/n2)
    return AverageBE(float(a.mean()-b.mean()),math.sqrt(v),float(df),len(groups),'parallel Welch' if welch else 'parallel pooled variance')

def power_parameters(n,design,analysis='contrast'):
    if design=='replicate': raise InputError('power: specify partial or full explicitly')
    m=3 if design=='partial' else 2
    if design not in (*DESIGNS,'parallel') or not isinstance(n,(int,np.integer)) or n<2*m or n%m:
        raise InputError('power requires an integer N, at least two subjects per sequence/arm and balanced allocation')
    if analysis not in ('contrast','ema'): raise InputError('unknown analysis')
    df=n-m
    if analysis=='ema' and design in ('partial','full'):
        p=len(DESIGNS[design][0]); df=(n-1)*(p-1)-1
    factor=math.sqrt({'parallel':4,'2x2':2,'partial':1.5,'full':1}[design]/n)
    return df,factor

def tost_power(n_total,cv,gmr,design='2x2',limits=(.8,1.25),analysis='contrast'):
    sigma=math.sqrt(math.log1p(positive(cv,'CV')**2)); delta=math.log(positive(gmr,'GMR'))
    lo,hi=map(math.log,check_limits(limits)); df,factor=power_parameters(n_total,design,analysis)
    crit=t.ppf(.95,df); sd=sigma*factor
    max_q=((hi-lo)/(2*crit*sd))**2*df
    def integrand(q):
        se=sd*math.sqrt(q/df)
        return max(0.,norm.cdf((hi-crit*se-delta)/sd)-norm.cdf((lo+crit*se-delta)/sd))*chi2.pdf(q,df)
    value,error=quad(integrand,0,max_q,epsabs=1e-9,epsrel=1e-8,limit=200)
    if error>1e-6: raise InputError('power integration did not converge sufficiently')
    return float(np.clip(value,0,1))

def check_limits(limits):
    try: lo,hi=limits; lo=positive(lo,'lower limit'); hi=positive(hi,'upper limit')
    except (ValueError,TypeError): raise InputError('limits must be two positive numbers')
    if lo>=hi: raise InputError('lower limit must be below upper limit')
    return lo,hi

def sample_size(cv,gmr,target,design='2x2',limits=(.8,1.25),analysis='contrast'):
    if not math.isfinite(target) or not 0<target<1: raise InputError('target power must be between 0 and 1')
    if design=='replicate': raise InputError('specify partial or full for power')
    lo,hi=check_limits(limits)
    if not lo<gmr<hi: raise InputError('sample-size solver requires a GMR strictly inside the limits')
    m=3 if design=='partial' else 2
    for n in range(2*m,5001,m):
        p=tost_power(n,cv,gmr,design,limits,analysis)
        if p>=target:return n,p
    raise InputError('target power not reached by N=5000')

def load_records(args):
    rows=read_table(args.input)
    cols=[args.subject_column,args.treatment_column,args.value_column]
    require_columns(rows,cols,str(args.input))
    if args.design!='parallel': require_columns(rows,['sequence','period'],str(args.input))
    out=[]
    for row in rows:
        tr=row[args.treatment_column].upper()
        tr={'TEST':'T','REF':'R','REFERENCE':'R'}.get(tr,tr)
        val=positive(row[args.value_column],'value')
        out.append(dict(subject=row[args.subject_column].strip(),treatment=tr,value=val,logvalue=math.log(val),sequence=row.get('sequence','').upper(),period=row.get('period','')))
    return out

def build_parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('-i','--input'); p.add_argument('--design',choices=(*DESIGNS,'replicate','parallel'),default='2x2')
    p.add_argument('--analysis',choices=('contrast','ema'),default='contrast')
    p.add_argument('--metric',default=''); p.add_argument('--scaling',choices=('none','abel','rsabe','both'),default='none')
    p.add_argument('--limits',default='0.80,1.25'); p.add_argument('--nti',action='store_true',help='fixed 90.00-111.11%% limits only; NOT a full NTI methodology')
    p.add_argument('--ema-rounding',action='store_true',help='compare CI endpoints as percentages rounded to 2 decimals')
    p.add_argument('--abel-justified',action='store_true',help='acknowledge clinical justification and prospective protocol specification')
    p.add_argument('--welch',action='store_true')
    for c in ('subject','treatment','value'): p.add_argument('--'+c+'-column',default=c)
    p.add_argument('--power',action='store_true'); p.add_argument('--cv',type=float); p.add_argument('--gmr',type=float,default=.95)
    p.add_argument('--target-power',type=float,default=.8); p.add_argument('--n',type=int)
    add_format_argument(p); return p

def run(argv=None):
    a=build_parser().parse_args(argv); report=Report()
    limits=check_limits((.9,1.1111) if a.nti else a.limits.split(','))
    if a.nti and a.limits!='0.80,1.25': raise InputError('--nti and custom --limits cannot be combined')
    if a.scaling!='none' and (a.nti or limits!=(.8,1.25)): raise InputError('scaled HVD criteria cannot be combined with NTI/custom limits')
    if a.scaling!='none' and a.design not in ('partial','full','replicate'): raise InputError('scaling requires a replicate design')
    if a.welch and a.design!='parallel': raise InputError('--welch is for parallel analysis only')
    if a.power:
        if a.scaling!='none' or a.welch or a.ema_rounding: raise InputError('power supports unrounded fixed-limit ABE with equal variances only; no ABEL/RSABE/Welch power')
        if a.cv is None: raise InputError('--power requires --cv')
        n,p=(a.n,tost_power(a.n,a.cv,a.gmr,a.design,limits,a.analysis)) if a.n is not None else sample_size(a.cv,a.gmr,a.target_power,a.design,limits,a.analysis)
        report.scalar('n_total',n); report.scalar('power',p); report.scalar('design',a.design); report.scalar('analysis',a.analysis)
        report.scalar('assumed_cv',a.cv); report.scalar('assumed_gmr',a.gmr); report.scalar('limits',limits)
        report.note('Numerical integration of fixed-limit TOST power; balanced complete design, normal log data, equal T/R within-subject variance, no subject-by-treatment interaction. Parallel CV is between-subject. No dropout allowance.')
        return report.emit(a.format)
    if not a.input: raise InputError('provide --input or --power')
    records=load_records(a); design,_=validate(records,a.design)
    result=parallel_design(records,a.welch) if design=='parallel' else (ema_analysis(records,design) if a.analysis=='ema' else contrast_analysis(records,design))
    report.scalar('design',design); report.scalar('metric',a.metric or a.value_column)
    for key,val in asdict(result).items(): report.scalar(key,val)
    report.scalar('gmr_pct',100*result.gmr); report.scalar('ci90_low_pct',100*result.ci_low); report.scalar('ci90_high_pct',100*result.ci_high)
    passed=passes_ci(result,limits,a.ema_rounding)
    report.scalar('average_be_met',passed); report.scalar('limits',limits)
    report.scalar('ci_rounding','percent, 2 decimal places, half-up' if a.ema_rounding else 'none')
    if a.nti: report.note('NTI flag changes fixed limits only. It does NOT implement FDA NTI scaling/variance comparison or all EMA endpoint-specific requirements.')
    rows=[]
    if a.scaling in ('abel','both'):
        if a.metric.lower()!='cmax': raise InputError('ABEL is restricted to --metric cmax')
        if not a.abel_justified: raise InputError('ABEL requires --abel-justified: clinical justification and prospective protocol specification must exist')
        ema=ema_analysis(records,design); rv=reference_variability(records,'ema',design); low,high,widened=abel_limits(rv)
        met=passes_ci(ema,(low,high),True) and .8<=ema.gmr<=1.25
        rows.append(dict(criterion='EMA ABEL (fixed-effects model)',met=bool(met),gmr=ema.gmr,ci_low=ema.ci_low,ci_high=ema.ci_high,df=ema.df,cvwr=rv.cvwr,s2wr=rv.s2wr,df_wr=rv.df,low=low,high=high,widened=widened))
    if a.scaling in ('rsabe','both'):
        contrast=contrast_analysis(records,design); rv=reference_variability(records,'contrast',design)
        if rv.swr<.294: raise InputError('sWR < 0.294: FDA unscaled heterogeneous mixed-model fallback is not implemented; no FDA decision produced')
        bound=rsabe_bound(contrast.estimate,contrast.se,contrast.df,rv)
        met=bound['passes_scaled_criterion'] and .8<=contrast.gmr<=1.25
        rows.append(dict(criterion='FDA HVD RSABE (complete-data contrast model)',met=bool(met),estimate=contrast.estimate,se=contrast.se,df=contrast.df,s2wr=rv.s2wr,df_wr=rv.df,cvwr=rv.cvwr,**bound))
    if rows:
        report.table('scaled criteria',rows)
        for row in rows:
            if not row['met']:report.finding(row['criterion']+': criterion not met')
    elif not passed: report.finding('Fixed-limit average BE criterion not met')
    report.note('Research implementation, not regulatory certification. Complete canonical designs only; see README for assumptions, unsupported cases and references.')
    return report.emit(a.format)

if __name__=='__main__': raise SystemExit(main_wrapper(run))
