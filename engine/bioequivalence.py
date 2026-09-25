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
from scipy.stats import chi2, f as f_dist, norm, t
from _common import InputError, Report, read_table, require_columns, add_format_argument, main_wrapper

DESIGNS = {'2x2': ('TR', 'RT'), 'partial': ('TRR', 'RTR', 'RRT'), 'full': ('TRTR', 'RTRT')}
THETA_FDA = (math.log(1.25)/0.25)**2
SFDA_MIN_SUBJECTS = 18

@dataclass
class AverageBE:
    estimate: float
    se: float
    df: float
    n_subjects: int
    method: str
    ci_level: float = 90.0
    @property
    def gmr(self): return math.exp(self.estimate)
    @property
    def _q(self): return t.ppf(1-(1-self.ci_level/100)/2,self.df)
    @property
    def ci_low(self): return math.exp(self.estimate-self._q*self.se)
    @property
    def ci_high(self): return math.exp(self.estimate+self._q*self.se)

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
        if len({r.get('stage','') for r in g})!=1: raise InputError(f'{sid}: inconsistent stage')
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

def period_key(r):
    """Period effect key; in a two-stage study periods are nested within stage."""
    return (str(r.get('stage','')), int(r['period']))

def fixed_fit(records, treatment=True):
    """OLS subject + period (+ treatment); sequence is absorbed by subject effects.
    With a stage column, periods are fitted within stage and the stage main
    effect is absorbed by the (stage-nested) subject effects.
    Reduced QR avoids forming the inverse of X'X. Complete designs only.
    """
    subjects=sorted({r['subject'] for r in records}); keys=sorted({period_key(r) for r in records})
    # One reference period per stage: the stage main effect lies in the subject space.
    first={}
    for s,p in keys: first.setdefault(s,p)
    periods=[k for k in keys if k[1]!=first[k[0]]]
    x=np.array([[1.]+[float(r['subject']==s) for s in subjects[1:]]+[float(period_key(r)==p) for p in periods]+([float(r['treatment']=='T')] if treatment else []) for r in records])
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

def ema_analysis(records, design, ci_level=90.0):
    _,groups=validate(records,design)
    e,se,df,_=fixed_fit(records)
    return AverageBE(e,se,df,len(groups),'fixed subject + period + treatment; common residual variance (EMA Method A structure)',ci_level)

def _rss(cols, y):
    x=np.column_stack(cols); beta,_,rank,_=np.linalg.lstsq(x,y,rcond=None)
    return float(np.sum((y-x@beta)**2)), int(rank)

def _dummies(keys, values):
    return [np.array([float(v==k) for v in values]) for k in keys[1:]]

def anova_table(records, design):
    """Fixed-effects ANOVA of log values: sequence, subject(sequence), period, formulation.

    Sequence and subject(sequence) partition the between-subject sum of squares
    (sequence is tested against subject(sequence)); period and formulation are
    adjusted for all other terms and tested against the residual. For complete
    crossover designs the formulation line equals the Type I, II and III sums of
    squares. With a stage column (two-stage design) the between-subject stratum is
    stage, sequence, sequence x stage and subject(sequence x stage), and periods
    are nested within stage; the formulation estimate equals that of fixed_fit.
    """
    design,groups=validate(records,design)
    y=np.array([r['logvalue'] for r in records]); n=len(y); one=np.ones(n)
    frm=np.array([float(r['treatment']=='T') for r in records])
    rows=[]
    def row(source,ss,df,err=None):
        rows.append(dict(source=source,df=int(df),ss=float(ss),ms=float(ss/df) if df else float('nan'),f=None,p=None,_err=err))
    if design=='parallel':
        rss_full,rk=_rss([one,frm],y); rss0,_=_rss([one],y)
        row('formulation',rss0-rss_full,1,'residual'); row('residual',rss_full,n-rk)
    else:
        col=lambda key:[key(r) for r in records]
        subs=sorted(set(col(lambda r:r['subject']))); U=_dummies(subs,col(lambda r:r['subject']))
        # Within-subject period effects: one reference period per stage (the stage
        # main effect is between-subject and is modelled explicitly below).
        keys=sorted(set(col(period_key))); first={}
        for s,p in keys: first.setdefault(s,p)
        P=[np.array([float(period_key(r)==k) for r in records]) for k in keys if k[1]!=first[k[0]]]
        full,rk_full=_rss([one,*U,*P,frm],y)
        no_f,rk_nf=_rss([one,*U,*P],y); no_p,rk_np=_rss([one,*U,frm],y)
        staged=len(first)>1
        # Between-subject stratum, sequential: [stage], sequence, [sequence x stage], subject(...).
        between=[('stage',lambda r:r.get('stage',''))] if staged else []
        between+=[('sequence',lambda r:r['sequence'])]
        if staged: between+=[('sequence×stage',lambda r:(r['sequence'],r.get('stage','')))]
        # Between-subject terms use subject-level contrasts only (within-subject
        # terms excluded): with stages, within-stage period dummies are not
        # orthogonal to the stage effect.
        prev,rk_prev=_rss([one],y); fixed=[]
        for name,key in between:
            vals=col(key); fixed=fixed+_dummies(sorted(set(vals),key=str),vals)
            cur,rk_cur=_rss([one,*fixed],y)
            row(name,prev-cur,rk_cur-rk_prev,'subject')
            prev,rk_prev=cur,rk_cur
        subj_only,rk_u=_rss([one,*U],y)
        subject_term='subject(sequence×stage)' if staged else 'subject(sequence)'
        row(subject_term,prev-subj_only,rk_u-rk_prev,'residual')
        for r in rows:
            if r['_err']=='subject': r['_err']=subject_term
        row('period(stage)' if staged else 'period',no_p-full,rk_full-rk_np,'residual')
        row('formulation',no_f-full,rk_full-rk_nf,'residual')
        row('residual',full,n-rk_full)
    ms={r['source']:r['ms'] for r in rows}; dfs={r['source']:r['df'] for r in rows}
    for r in rows:
        e=r.pop('_err')
        if e and r['df']>0 and ms[e]>0:
            r['f']=r['ms']/ms[e]; r['p']=float(f_dist.sf(r['f'],r['df'],dfs[e]))
    rows.append(dict(source='total (corrected)',df=n-1,ss=float(np.sum((y-y.mean())**2)),ms=None,f=None,p=None))
    return rows

def low_reference_auc(records):
    """GCC/SFDA guideline exception: a reference AUC below 5% of the reference
    geometric mean (computed without that subject). Flagged, never excluded
    automatically: exclusion is exceptional and must be pre-specified."""
    ref=[r for r in records if r['treatment']=='R']; out=[]
    for r in ref:
        others=[x['logvalue'] for x in ref if x['subject']!=r['subject']]
        if not others: continue
        gm=math.exp(float(np.mean(others)))
        if r['value']<0.05*gm: out.append(dict(subject=r['subject'],period=r['period'],reference_auc=r['value'],reference_gm=gm,pct_of_gm=100*r['value']/gm))
    return out

def predose_screen(records, design):
    """GCC/SFDA BE guideline: exclude a subject-period whose pre-dose concentration
    exceeds 5% of that period's Cmax. Returns (kept_records, excluded_rows)."""
    flagged=[r for r in records if r.get('predose') is not None and r['predose']>0.05*r['cmax_ref']]
    if not flagged: return records,[]
    info=[dict(subject=r['subject'],period=r['period'],predose=r['predose'],cmax=r['cmax_ref'],ratio_pct=100*r['predose']/r['cmax_ref']) for r in flagged]
    if design not in ('2x2','parallel'):
        raise InputError('pre-dose > 5% of Cmax in '+', '.join(f"{i['subject']} period {i['period']}" for i in info)+
                         ': the guideline excludes only that period, which leaves an incomplete replicate design; incomplete designs are not supported')
    drop={r['subject'] for r in flagged}
    return [r for r in records if r['subject'] not in drop],info

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

def sample_size(cv,gmr,target,design='2x2',limits=(.8,1.25),analysis='contrast',min_n=0):
    if not math.isfinite(target) or not 0<target<1: raise InputError('target power must be between 0 and 1')
    if design=='replicate': raise InputError('specify partial or full for power')
    lo,hi=check_limits(limits)
    if not lo<gmr<hi: raise InputError('sample-size solver requires a GMR strictly inside the limits')
    m=3 if design=='partial' else 2
    for n in range(max(2*m,-(-min_n//m)*m),5001,m):
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
        rec=dict(subject=row[args.subject_column].strip(),treatment=tr,value=val,logvalue=math.log(val),sequence=row.get('sequence','').upper(),period=row.get('period',''))
        if row.get('stage','')!='': rec['stage']=row['stage'].strip()
        if getattr(args,'predose_check',False) and row.get(args.predose_column,'')!='':
            try: pre=float(row[args.predose_column])
            except ValueError: raise InputError('pre-dose concentration must be a number')
            if not math.isfinite(pre) or pre<0: raise InputError('pre-dose concentration must be a finite non-negative number')
            ref=val if args.cmax_column in ('',args.value_column) else positive(row.get(args.cmax_column),'cmax (pre-dose check)')
            rec['predose']=pre; rec['cmax_ref']=ref
        out.append(rec)
    return out

def build_parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('-i','--input'); p.add_argument('--design',choices=(*DESIGNS,'replicate','parallel'),default='2x2')
    p.add_argument('--analysis',choices=('contrast','ema'),default=None,help='default: contrast; ema (fixed-effects ANOVA) under --profile sfda')
    p.add_argument('--profile',choices=('none','sfda'),default='none',help='sfda: GCC/SFDA BE guideline DS-G-010 V3.1 rules (see README)')
    p.add_argument('--ci-level',type=float,default=90.0,help='confidence level in percent, e.g. 94.12 for a pre-specified two-stage design')
    p.add_argument('--anova',action='store_true',help='add the fixed-effects ANOVA table (always on with --profile sfda)')
    p.add_argument('--predose-column',default='predose'); p.add_argument('--cmax-column',default='',help='Cmax column for the pre-dose check when the analysed value is not Cmax')
    p.add_argument('--endogenous',action='store_true',help='endogenous substance: skip the pre-dose exclusion rule')
    p.add_argument('--metric',default=''); p.add_argument('--scaling',choices=('none','abel','rsabe','both'),default='none')
    p.add_argument('--limits',default='0.80,1.25'); p.add_argument('--nti',action='store_true',help='fixed 90.00-111.11%% limits only; NOT a full NTI methodology')
    p.add_argument('--ema-rounding',action='store_true',help='compare CI endpoints as percentages rounded to 2 decimals')
    p.add_argument('--abel-justified',action='store_true',help='acknowledge clinical justification and prospective protocol specification')
    p.add_argument('--welch',action='store_true')
    for c in ('subject','treatment','value'): p.add_argument('--'+c+'-column',default=c)
    p.add_argument('--power',action='store_true'); p.add_argument('--cv',type=float); p.add_argument('--gmr',type=float,default=.95)
    p.add_argument('--target-power',type=float,default=.8); p.add_argument('--n',type=int)
    add_format_argument(p); return p

def apply_profile(a):
    """Resolve defaults and enforce the GCC/SFDA guideline (DS-G-010 V3.1) under --profile sfda."""
    if not 50<a.ci_level<100: raise InputError('--ci-level must be between 50 and 100')
    if a.profile!='sfda':
        a.analysis=a.analysis or 'contrast'; a.predose_check=False; return
    if a.analysis=='contrast': raise InputError('SFDA profile: the guideline requires a fixed-effects ANOVA (sequence, subject within sequence, period, formulation); use --analysis ema or omit it')
    if a.scaling in ('rsabe','both'): raise InputError('SFDA profile: reference-scaled ABE (RSABE) is not an SFDA method; the guideline provides ABEL (widened limits for Cmax) only')
    if a.welch: raise InputError('SFDA profile: the guideline requires the ANOVA model; Welch is not provided for')
    a.analysis='ema'; a.ema_rounding=True; a.anova=True; a.predose_check=not a.endogenous

def run(argv=None):
    a=build_parser().parse_args(argv); report=Report(); apply_profile(a)
    limits=check_limits((.9,1.1111) if a.nti else a.limits.split(','))
    if a.nti and a.limits!='0.80,1.25': raise InputError('--nti and custom --limits cannot be combined')
    if a.scaling!='none' and (a.nti or limits!=(.8,1.25)): raise InputError('scaled HVD criteria cannot be combined with NTI/custom limits')
    if a.scaling!='none' and a.design not in ('partial','full','replicate'): raise InputError('scaling requires a replicate design')
    if a.welch and a.design!='parallel': raise InputError('--welch is for parallel analysis only')
    if a.power:
        # The SFDA profile turns CI rounding on for analyses; prospective power ignores rounding (see README).
        if a.scaling!='none' or a.welch or (a.ema_rounding and a.profile!='sfda'): raise InputError('power supports unrounded fixed-limit ABE with equal variances only; no ABEL/RSABE/Welch power')
        if a.cv is None: raise InputError('--power requires --cv')
        min_n=SFDA_MIN_SUBJECTS if a.profile=='sfda' else 0
        if a.n is not None and a.n<min_n: raise InputError(f'SFDA profile: at least {min_n} evaluable subjects are required')
        n,p=(a.n,tost_power(a.n,a.cv,a.gmr,a.design,limits,a.analysis)) if a.n is not None else sample_size(a.cv,a.gmr,a.target_power,a.design,limits,a.analysis,min_n)
        report.scalar('n_total',n); report.scalar('power',p); report.scalar('design',a.design); report.scalar('analysis',a.analysis)
        report.scalar('assumed_cv',a.cv); report.scalar('assumed_gmr',a.gmr); report.scalar('limits',limits); report.scalar('profile',a.profile)
        if a.profile=='sfda': report.note(f'SFDA profile: fixed-effects ANOVA degrees of freedom; minimum {SFDA_MIN_SUBJECTS} evaluable subjects (24 generally recommended). Plan for drop-outs separately.')
        report.note('Numerical integration of fixed-limit TOST power; balanced complete design, normal log data, equal T/R within-subject variance, no subject-by-treatment interaction. Parallel CV is between-subject. No dropout allowance.')
        return report.emit(a.format)
    if not a.input: raise InputError('provide --input or --power')
    records=load_records(a); design,_=validate(records,a.design)
    if any(r.get('stage') for r in records) and a.analysis!='ema' and design!='parallel':
        raise InputError('a stage column requires the fixed-effects model (--analysis ema or --profile sfda)')
    excluded=[]
    if a.predose_check:
        records,excluded=predose_screen(records,design)
        if excluded: design,_=validate(records,design)
    if design=='parallel': result=parallel_design(records,a.welch)
    else: result=ema_analysis(records,design,a.ci_level) if a.analysis=='ema' else contrast_analysis(records,design)
    result.ci_level=a.ci_level
    report.scalar('design',design); report.scalar('metric',a.metric or a.value_column)
    for key,val in asdict(result).items(): report.scalar(key,val)
    report.scalar('gmr_pct',100*result.gmr); report.scalar('ci90_low_pct',100*result.ci_low); report.scalar('ci90_high_pct',100*result.ci_high)
    report.scalar('profile',a.profile)
    passed=passes_ci(result,limits,a.ema_rounding)
    report.scalar('average_be_met',passed); report.scalar('limits',limits)
    report.scalar('ci_rounding','percent, 2 decimal places, half-up' if a.ema_rounding else 'none')
    if a.nti: report.note('NTI flag changes fixed limits only. It does NOT implement FDA NTI scaling/variance comparison or all EMA endpoint-specific requirements.')
    if a.ci_level!=90: report.note(f'{a.ci_level:g}% confidence interval (the ci90_* fields hold this interval). Adjusted levels are for a pre-specified two-stage design.')
    if excluded:
        report.table('pre-dose exclusions',excluded)
        report.note(f'Pre-dose concentration > 5% of Cmax: {len(excluded)} subject(s) removed from the analysis (GCC/SFDA BE guideline, carry-over).')
    if a.predose_check and any('predose' in r for r in records):
        report.note('Pre-dose rule applied as for single-dose studies. Do not supply a pre-dose column for steady-state (multiple-dose) studies, where pre-dose concentrations are expected.')
    if a.anova:
        at=anova_table(records,design); report.table('anova',at)
        res=next(r for r in at if r['source']=='residual')
        if design!='parallel': report.scalar('cv_intra_pct',100*math.sqrt(math.expm1(res['ms'])))
    if a.profile=='sfda':
        if 'auc' in (a.metric or '').lower():
            low=low_reference_auc(records)
            if low:
                report.table('low reference AUC',low)
                report.note('Reference AUC < 5% of the reference geometric mean: listed only. The guideline allows exclusion in exceptional cases, pre-specified in the protocol; this program does not exclude them.')
        if a.endogenous: report.note('Endogenous substance: values must already be baseline-corrected by the pre-specified method (subtraction of the mean pre-dose concentration or of the pre-dose AUC is preferred); this program does not perform baseline correction.')
        if result.n_subjects<SFDA_MIN_SUBJECTS: report.finding(f'SFDA: {result.n_subjects} evaluable subjects; the guideline requires at least {SFDA_MIN_SUBJECTS}')
        report.note('SFDA profile (GCC BE guideline DS-G-010 V3.1): fixed-effects ANOVA, CI bounds rounded to two decimals, ABEL for Cmax only, no RSABE, minimum 18 evaluable subjects, pre-dose > 5% of Cmax exclusion'+(' (skipped: endogenous)' if a.endogenous else '')+'.')
    rows=[]
    if a.scaling in ('abel','both'):
        if a.metric.lower()!='cmax': raise InputError('ABEL is restricted to --metric cmax')
        if not a.abel_justified: raise InputError('ABEL requires --abel-justified: clinical justification and prospective protocol specification must exist')
        ema=ema_analysis(records,design,a.ci_level); rv=reference_variability(records,'ema',design); low,high,widened=abel_limits(rv)
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
