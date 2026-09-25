# Bioequivalence — corrected research implementation

This package corrects the reviewed K-Dense script and replaces its statistical core with explicit, separate models. It is not regulatory certification, a full validation package, or a substitute for the applicable current product-specific guidance/statistical analysis plan.

## Quick start

Python 3.10+ with NumPy and SciPy. No network connection is used by the program.

```bash
python -m pip install -r requirements.txt
python -m unittest -v
python validate_simulation.py > simulation_results.json
python bioequivalence.py -i example_partial.csv --design partial --metric cmax --scaling both --abel-justified --format json
python bioequivalence.py --power --design partial --cv .30 --gmr .95 --n 24 --format json
python bioequivalence.py --power --design full --cv .30 --gmr .95 --target-power .80
python bioequivalence.py -i study.csv --design 2x2 --format json
```

`--abel-justified` is an acknowledgement by the analyst that clinical justification and prospective specification exist. The software cannot verify these documents. The example CSV is simulated, NOT a clinical study.

## Changes

- RSABE: subtract SE squared from the squared contrast; use chi-square 0.95, not 0.05; report the inputs and upper bound. Uses the FDA 2011 progesterone example as an explicit computational reference, not a claim that the historical guidance applies to every product today.
- Subject contrast estimator: variance denominator is the square of the number of sequences (9 for partial; 4 for two-sequence designs).
- Separate `partial` (TRR/RTR/RRT) and `full` (TRTR/RTRT) designs. `replicate` remains a data-analysis alias which infers only these designs. Power requires the explicit design.
- Complete-case input is required: no silent dropping, imputing sequence labels or averaging duplicates. Reject missing labels, nonfinite/nonpositive values, conflicting sequence assignments, missing/duplicate periods, and treatment/sequence/period mismatches. All canonical sequences and at least two complete subjects per sequence are required. Parallel input requires unique subject IDs and one record per subject.
- FDA reference variance: ordered R differences, centered within sequence, with residual degrees of freedom N minus number of sequences. Period effects are removed through sequence-specific contrasts.
- EMA reference variance: separate reference-only fixed subject + period model. Degrees of freedom come from matrix rank; they are NOT forced to N minus number of sequences. In full replication, aliased nuisance columns are removed without changing their column space.
- EMA ABEL: separate fixed-effects subject + period + treatment model (sequence is absorbed by subject), reference-only variance, Cmax restriction, 30% trigger, 50% CV cap, fixed 80–125% point-estimate constraint. CI comparisons use percentages rounded to two decimals (half-up). No expansion for AUC.
- `--nti` changes fixed limits only. It does not implement FDA NTI scaling, variance-ratio testing, or endpoint-specific EMA NTI policy. Scaling plus NTI/custom limits is rejected.
- Parallel pooled variance remains available; optional `--welch` implements unequal-variance inference.
- Power integrates over the chi-square distribution with numerical error checking. Fixed-limit ABE only: ABEL/RSABE power is not implemented and requests are rejected. Balanced allocation is enforced, with sample-size increments of 3 for partial designs and 2 otherwise. No dropout adjustment.

## Models and limits

Default `--analysis contrast`: equally weighted sequence means of individual T-minus-R means; common variance of subject contrasts; residual df N-m. Matches ordinary complete 2x2 analysis and the complete-data contrast structure in the FDA HVD example. Supported full-design contrast fitting assumes common residual contrast variance; it does not implement arbitrary mixed covariance structures or incomplete subjects.

`--analysis ema`: fixed subject, period and treatment model, common residual variance. ABEL always uses this model, regardless of the primary analysis flag. RSABE always uses the contrast model. With `--scaling both`, the two rows intentionally have different estimators/variance estimates; neither is derived by merely changing the other's acceptance limits.

For FDA HVD RSABE, the historical reference threshold is sWR >= 0.294. Below that threshold the program rejects a FDA decision because the unscaled heterogeneous mixed-model fallback is not implemented. You may perform an explicitly selected fixed-limit analysis separately, but it must not be described as that FDA fallback.

Power assumes independent normal log errors, equal T/R within-subject variances, balanced complete allocation and no subject-by-treatment interaction. Between-subject CV is required for a parallel design. Contrast power uses df N-m and variance factors 2/N, 1.5/N and 1/N for 2x2, partial and full. EMA fixed-model power uses df (N-1)(P-1)-1, with the same balanced variance factors. The supplied Monte Carlo checks exercise contrast power only. CI rounding is not incorporated in prospective power.

Incomplete studies, other sequence sets, carryover models, FDA NTI methods, RSABE/ABEL power, Welch power and regulatory end-to-end validation are outside this implementation. Failures are explicit. API signatures and JSON fields changed; this is not an unconditional drop-in replacement for existing Python callers. Table/TSV/JSON output and exit codes use the original `_common.py` helper: 0 pass, 1 criterion not met, 2 invalid/unsupported input. For scaling, exit status follows requested scaled criteria rather than descriptive ordinary ABE.

## Verification

`test_bioequivalence.py` covers independent full-data OLS comparison for unbalanced 2x2, independent sequence-GLM contrast comparison, literal FDA reference algebra, invariance of reference variance to period effects, degrees of freedom, invalid data, balanced allocation/sample-size minimality, rounding/cap, Welch comparison to SciPy and CLI guards.

`validate_simulation.py` generates raw log observations and calls the actual contrast analysis functions for 4,000 studies per design (12,000 total), seed 20260925. The included JSON records assumptions, versions, power, CI coverage and Monte Carlo error. Tests are internal verification, not independent SAS/R validation. No real clinical dataset was used, no full EMA published-dataset validation was performed, and no broad RSABE type-I-error study was performed.

## Provenance

Original source fetched 2026-09-25 from:
https://github.com/K-Dense-AI/scientific-agent-skills/blob/main/skills/pkpd-modeling/scripts/bioequivalence.py

`original_bioequivalence.py` preserves the fetched original; `_common.py` is its upstream dependency. `provenance.json` contains SHA-256 hashes. This package does not claim the fetched main branch is commit 49c6e97; that earlier review identifier was not independently resolved. `changes.patch` provides a unified diff. See upstream license.

Computational references:
- FDA 2011 progesterone example, pages 5–8: https://www.accessdata.fda.gov/drugsatfda_docs/psg/Progesterone_caps_19781_RC02-11.pdf
- EMA bioequivalence guideline, sections 4.1.8–4.1.10: https://www.ema.europa.eu/en/documents/scientific-guideline/guideline-investigation-bioequivalence-rev1_en.pdf
- EMA Method A discussion: https://www.ema.europa.eu/en/human-regulatory-overview/research-development/scientific-guidelines/clinical-pharmacology-pharmacokinetics-guidelines/clinical-pharmacology-pharmacokinetics-questions-answers

These references establish the implemented model/formula choices; they do not establish regulatory acceptability of this software for a particular submission.
