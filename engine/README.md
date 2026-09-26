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

`--abel-justified` is an acknowledgement by the analyst that clinical justification and prospective specification exist, and that CVwR is a reliable estimate, not the result of outliers (EMA 4.1.10; GCC 3.1.10; the program does not test for outliers). The software cannot verify these documents. The example CSV is simulated, NOT a clinical study.

## Changes

- RSABE: subtract SE squared from the squared contrast; use chi-square 0.95, not 0.05; report the inputs and upper bound. Uses the FDA 2011 progesterone example as an explicit computational reference, not a claim that the historical guidance applies to every product today.
- Subject contrast estimator: variance denominator is the square of the number of sequences (9 for partial; 4 for two-sequence designs).
- Separate `partial` (TRR/RTR/RRT) and `full` (TRTR/RTRT) designs. `replicate` remains a data-analysis alias which infers only these designs. Power requires the explicit design.
- Complete-case input is required (later relaxed for the fixed-effects model, see *Additions* below): no silent dropping, imputing sequence labels or averaging duplicates. Reject missing labels, nonfinite/nonpositive values, conflicting sequence assignments, missing/duplicate periods, and treatment/sequence/period mismatches. All canonical sequences and at least two complete subjects per sequence are required. Parallel input requires unique subject IDs and one record per subject.
- FDA reference variance: ordered R differences, centered within sequence, with residual degrees of freedom N minus number of sequences. Period effects are removed through sequence-specific contrasts.
- EMA reference variance: separate reference-only fixed subject + period model. Degrees of freedom come from matrix rank; they are NOT forced to N minus number of sequences. In full replication, aliased nuisance columns are removed without changing their column space.
- EMA ABEL: separate fixed-effects subject + period + treatment model (sequence is absorbed by subject), reference-only variance, Cmax restriction, 30% trigger, 50% CV cap, fixed 80–125% point-estimate constraint. CI comparisons use percentages rounded to two decimals (half-up); since version 1.0.0 the point estimate is compared the same way (80.00–125.00%). No expansion for AUC.
- `--nti` changes fixed limits only. It does not implement FDA NTI scaling, variance-ratio testing, or endpoint-specific EMA NTI policy. Scaling plus NTI/custom limits is rejected.
- Parallel pooled variance remains available; optional `--welch` implements unequal-variance inference.
- Power integrates over the chi-square distribution with numerical error checking. Fixed-limit ABE only: ABEL/RSABE power is not implemented and requests are rejected. Balanced allocation is enforced, with sample-size increments of 3 for partial designs and 2 otherwise. No dropout adjustment.

## Models and limits

Default `--analysis contrast`: equally weighted sequence means of individual T-minus-R means; common variance of subject contrasts; residual df N-m. Matches ordinary complete 2x2 analysis and the complete-data contrast structure in the FDA HVD example. Supported full-design contrast fitting assumes common residual contrast variance; it does not implement arbitrary mixed covariance structures or incomplete subjects.

`--analysis ema`: fixed subject, period and treatment model, common residual variance. It accepts subjects with missing periods (absent rows, or a blank / `NA` value) and uses all available observations, as EMA Method A does; such subjects are listed (`incomplete subjects`, `missing values`), `n_subjects` counts subjects with both products, and in replicate designs a pre-dose > 5% period is removed alone. The contrast analysis, RSABE and the FDA NTI method still require complete subjects and say so. ABEL always uses this model, regardless of the primary analysis flag. RSABE always uses the contrast model. With `--scaling both`, the two rows intentionally have different estimators/variance estimates; neither is derived by merely changing the other's acceptance limits.

For FDA HVD RSABE, the historical reference threshold is sWR >= 0.294. Below that threshold the program rejects a FDA decision because the unscaled heterogeneous mixed-model fallback is not implemented. You may perform an explicitly selected fixed-limit analysis separately, but it must not be described as that FDA fallback.

Power assumes independent normal log errors, equal T/R within-subject variances, balanced complete allocation and no subject-by-treatment interaction. Between-subject CV is required for a parallel design. Contrast power uses df N-m and variance factors 2/N, 1.5/N and 1/N for 2x2, partial and full. EMA fixed-model power uses df (N-1)(P-1)-1, with the same balanced variance factors. The supplied Monte Carlo checks exercise contrast power only. CI rounding is not incorporated in prospective power.

Incomplete subjects are supported only by the fixed-effects model (`--analysis ema`, ABEL, `--profile sfda`; see below). Other sequence sets, carryover models, EMA Method B, RSABE/ABEL power, Welch power and regulatory end-to-end validation are outside this implementation. Failures are explicit. API signatures and JSON fields changed; this is not an unconditional drop-in replacement for existing Python callers. Table/TSV/JSON output and exit codes use the original `_common.py` helper: 0 pass, 1 criterion not met, 2 invalid/unsupported input. For scaling, exit status follows requested scaled criteria rather than descriptive ordinary ABE.

## Verification

`test_bioequivalence.py` covers independent full-data OLS comparison for unbalanced 2x2, independent sequence-GLM contrast comparison, literal FDA reference algebra, invariance of reference variance to period effects, degrees of freedom, invalid data, balanced allocation/sample-size minimality, rounding/cap, Welch comparison to SciPy and CLI guards.

`validate_simulation.py` generates raw log observations and calls the actual contrast analysis functions for 4,000 studies per design (12,000 total), seed 20260925. The included JSON records assumptions, versions, power, CI coverage and Monte Carlo error. Tests are internal verification, not independent SAS/R validation. No real clinical dataset was used, no full EMA published-dataset validation was performed, and no broad RSABE type-I-error study was performed.

## Provenance

Original source fetched 2026-09-25 from:
https://github.com/K-Dense-AI/scientific-agent-skills/blob/main/skills/pkpd-modeling/scripts/bioequivalence.py

`original_bioequivalence.py` and `original_nca.py` preserve the fetched originals; `_common.py` is their upstream dependency (changed once, in commit 7f42c0b: `parse_float` rejects non-finite numbers). `provenance.json` records the SHA-256 of each upstream file as fetched. This package does not claim the fetched main branch is commit 49c6e97; that earlier review identifier was not independently resolved. `changes.patch` is the unified diff from `original_bioequivalence.py` to the current `bioequivalence.py` (regenerate with `git diff --no-index original_bioequivalence.py bioequivalence.py`). See upstream license.

Computational references:
- FDA 2011 progesterone example, pages 5–8: https://www.accessdata.fda.gov/drugsatfda_docs/psg/Progesterone_caps_19781_RC02-11.pdf
- EMA bioequivalence guideline, sections 4.1.8–4.1.10: https://www.ema.europa.eu/en/documents/scientific-guideline/guideline-investigation-bioequivalence-rev1_en.pdf
- EMA Method A discussion: https://www.ema.europa.eu/en/human-regulatory-overview/research-development/scientific-guidelines/clinical-pharmacology-pharmacokinetics-guidelines/clinical-pharmacology-pharmacokinetics-questions-answers

These references establish the implemented model/formula choices; they do not establish regulatory acceptability of this software for a particular submission.

## Additions in the bioequivalence-lab repository

This copy of the engine extends the corrected package above (the `provenance.json` hashes refer to the package as received):

- `--profile sfda`: GCC/SFDA bioequivalence guideline DS-G-010 V3.1 rules — fixed-effects ANOVA, CI rounding to two decimals, ABEL for Cmax only, no RSABE, minimum 18 evaluable subjects, pre-dose > 5% of Cmax exclusion, ANOVA table output.
- `--anova`: fixed-effects ANOVA table (sequence, subject(sequence), period, formulation) for any analysis.
- `--ci-level`: adjusted confidence level (e.g. 94.12 for a pre-specified two-stage design); an optional `stage` column nests periods within stage.
- `--predose-column`, `--cmax-column`, `--endogenous`: pre-dose screening inputs.
- `sample_size(..., min_n=)`: minimum total N (18 under the SFDA profile).

- `--potency-test`, `--potency-reference`, `--potency-correction`: batch assay content (% of label claim); correction (value x 100 / content) only when declared pre-specified and the batches differ by more than 5 points (ICH M13A 2.2.2.3; GCC 3.1.8). The decision uses the corrected data; the `potency correction` table gives, for the uncorrected and the corrected data, the same criteria (fixed limits and any ABEL, RSABE or FDA NTI row) with the limits used, as ICH M13A 2.2.2.3 asks. Tests: `test_potency.py`, `test_launch_review.py`.
- Subjects without data for both products are left out, and listed in `subjects without both products`, in multi-treatment comparisons and in 2×2 fixed-effects analyses (ICH M13A 2.2.3.2; GCC 3.1.8). Replicate designs keep them (EMA Method A: they inform the period effects and CVwR).

ABEL is checked against the EMA reference datasets and replicateBE in [`validation/ema_abel/`](validation/ema_abel/README.md) (`test_ema_reference.py`).

Default behaviour without these options is unchanged (reference output in `example_results.json` reproduced to floating-point precision). New tests: `test_sfda.py`.

## NCA module (`nca.py`)

`nca.py` derives concentration-time parameters for BE studies and feeds `bioequivalence.py`. It is derived from upstream `skills/pkpd-modeling/scripts/nca.py`; `original_nca.py` is the file as fetched on 2026-09-26 (SHA-256 `a97bbc4a…23b11`, byte-identical to upstream commit `4c90a52`).

Fixes against the upstream file, each reproduced before the change:

1. **kel window:** upstream kept the *shorter* window when adjusted r² values were within 1e-4. `nca.py` follows PKNCA `pk.calc.half.life` (checked in its source, `R/half.life.R`, 2026-09-26, and in the installed 0.12.1; results match PKNCA, see below): best adjusted r² over all terminal windows after tmax, then among windows with kel > 0 and adjusted r² > best − 1e-4 the one with the most points. The same rule is commonly attributed to Phoenix WinNonlin; that was not checked here.
2. **Missed samples:** upstream turned a blank cell into a BLQ zero (AUC 5.8% low in the test case). A blank / `NA` is now a missed sample: dropped, listed in `missed samples`, and flagged as a protocol deviation. Only BLQ tokens or values below `--lloq` become zero (ICH M13A 2.2.2.2).
3. **Trapezoidal rule:** upstream called linear-up/log-down "the usual regulatory choice". ICH M13A 2.2.2.2 gives the linear rule as its example and asks for the method to be reported. Default is `linear`; `--auc-method linup-logdown` is available; the method is always reported.
4. **Missing rules added:** study-level 80/20 coverage rule; AUC(0-72h) (`--truncate-72`); steady-state CtauSS, CminSS, CavSS, fluctuation and swing per the ICH M13A glossary (`--tau`). CtauSS is reported only when a sample lies within 10 min of tau: a program rule built on ICH M13A 2.1.8, which recommends the last sample within 10 min of the nominal tau.

Rules implemented, with sources checked in the guideline text:

| Rule | Source |
|---|---|
| BLQ = 0 in PK calculations; BLQ omitted from kel and t½ | ICH M13A 2.2.2.2 |
| Actual sampling times | ICH M13A 2.2.2.2; GCC 3.1.5 |
| kel from ≥ 3 terminal points; number of points reported | ICH M13A 2.1.8, 2.2.2.2; GCC 3.1.4, 3.1.8 |
| AUC(0-t)/AUC(0-inf) < 80% in > 20% of observations → validity may need discussion | ICH M13A 2.2.2.2; GCC 3.1.8 |
| AUC(0-72h) replaces AUC(0-t); AUC(0-inf), kel, t½ not required | ICH M13A 2.2.2.2 (GCC 3.1.5 adds: only when C(72 h) is quantifiable — flagged under `--profile sfda`) |
| Pre-dose > 5% of Cmax → exclude the period (single dose) | ICH M13A 2.2.3.3; GCC 3.1.8 |
| Baseline correction per period, negatives set to zero; PK and statistical analyses of corrected **and** uncorrected data | ICH M13A 3.1 (mean or median, or time-matched, pre-specified); GCC 3.1.5 prefers subtracting the mean pre-dose concentration or pre-dose AUC. `nca.py` subtracts the mean and reports uncorrected Cmax and AUC(0-t); the full uncorrected analysis is a second run without `--baseline` |
| Emesis at or before 2 × median tmax (IR) / during the dosing interval (MR) → exclude | GCC 3.1.8 states the rule (`--profile sfda`). ICH M13A Q&A 2.9 lists emesis within 2 × the *expected* median tmax as an example deviation that may justify removal if pre-specified; without `--profile sfda` emesis is flagged only |
| Summary statistics: n, geometric mean, CV, median, arithmetic mean, SD, min, max | ICH M13A 2.2.2.2; GCC 3.1.8 |
| Annex 1 section 5 table (arithmetic mean and CV% for T and R) | GCC DS-G-010 V3.1 Annex 1 |

Emesis rule: `--median-tmax` takes the expected median tmax from the protocol (one value, or per product, e.g. `T=1.5,R=2`), as ICH M13A Q&A 2.9 ("expected median tmax") and GCC 3.1.8 (exclusion decided before bioanalysis) imply. Without it, the median observed in the study (per product, periods without emesis) is used and a finding says so. The pre-dose value is the highest sample at time ≤ 0.

Summary statistics cover the analysis population: subjects left out of the BE input (drop-outs with missing periods, subjects with an excluded period in 2×2 and replicate designs) stay in the individual listings but not in the summary (GCC 3.1.8: drop-out data "should not be included in the summary statistics"). The 80/20 coverage count also reports the profiles without an estimable kel.

AUC(0-72h) under `--profile sfda`: C(72 h) counts as quantifiable only when the last quantifiable sample is at or after 72 h (10-minute tolerance, a program convention); otherwise GCC 3.1.5 requires AUC(0-∞) and the residual area, and a finding says so. ICH M13A 2.1.8.2 limits truncation to immediate-release drugs with a half-life of 24 h or longer; that condition is the user's to check.

`be input` table: one row per included subject-period (subject, sequence, period, treatment, parameters, pre-dose), ready for `bioequivalence.py`. Columns: Cmax, AUC(0-t), AUC(0-72h) when computed, AUC(0-∞) (an *additional* parameter in ICH M13A 2.2.2.2, offered for supportive analysis), and any `--partial-auc` (early exposure, primary where applicable); at steady state AUC(0-τ) and Cmax,ss.

- **2×2:** a subject with an excluded or missing period is left out (GCC 3.1.8; ICH M13A Q&A 2.9: in a 2-way crossover the subject leaves the analysis). A pre-dose > 5% period stays in; `bioequivalence.py` applies and reports the pre-dose rule.
- **Replicate designs:** only the excluded periods leave; subjects with an excluded or missing period stay in and are analysed by the fixed-effects model (EMA Method A), which the UI selects. A replicate subject with pre-dose > 5% of Cmax keeps its other periods.
- **More than two products** (e.g. Williams): only the excluded periods are removed; each comparison keeps the subjects that have both products.
- **Parallel studies** (one profile per subject, no period column) are supported.
- A parameter that is not available for a subject (e.g. AUC(0-∞) without an estimable kel) is named in a note. When that parameter is sent, the UI leaves out the subject in 2×2 and parallel studies and only that period in replicate and multi-treatment designs; `bioequivalence.py` names subject and period if an empty value reaches it.
- **Codes:** treatment codes become T/R (`Test`/`Ref`/`Reference`, or two other codes named with `--test-label` / `--reference-label`, chosen in the NCA tab); periods become integers (`2.0` → 2, `P1`/`P2` → 1/2 in natural order); a `sequence` column that is not in T/R codes (e.g. 1/2, AB/BA) is replaced by the sequence derived from the treatment order. A subject with missing periods gets the one full sequence that fits its observed periods, or is left out with a note when none or several fit. Every mapping is stated in a note.
- **`stage` and `group` columns** (two-stage design, multi-group study) are read, must be constant within a subject, and are carried into the `be input`; the UI then selects the fixed-effects model (and 94.12% for two stages).
- **GCC Annex 1 table:** built for T/R, Test/Reference or the chosen labels; a note says when it cannot be built.

Verification (`test_nca.py`, 37 tests): 297 random profiles with BLQ values (296 with an estimable kel) match a second implementation (`scipy.stats.linregress` over every terminal window, `scipy.integrate.trapezoid`) — AUC to 9 decimals, lambda_z to 10, AUC(0-∞) to 8. That reference re-implements the same kel selection rule, so for that rule the independent evidence is the PKNCA comparison below; closed-form one-compartment profiles (AUC(0-inf), AUC(0-72h), lambda_z); each upstream fix; steady-state parameters against hand calculation; baseline, emesis, coverage and 72 h rules; hand-off into `bioequivalence.py` (GMR equals the geometric mean of within-subject ratios); and 17 regression tests for the findings of the review below. **Browser engine against local Python** (`validation/browser_parity.mjs`): 20 command lines (5 NCA runs including the 348 single-dose and 60 steady-state validation profiles, 15 BE runs covering every design, framework, scaled criterion, potency correction and power) run through the page's Pyodide worker and through the local engines; 24,146 values compared. Exit codes, text, decisions and table structure are identical; NCA values agree to 2.6 × 10⁻¹⁶ and their SHA-256 fingerprint at 10 significant digits (23,507 values) is equal; the largest BE difference is 4.6 × 10⁻⁹ (relative), in the RSABE upper bound near zero, where Pyodide's SciPy and the local SciPy differ in the last digits. The script exits 1 on any difference above 1 × 10⁻⁸ or a different NCA fingerprint; the result of the last run is in `browser_parity.json`.

**Review of 2026-09-26.** Five independent reviews (engine, regulatory text, UI, NCA→BE integration, tests and docs) found, and this version fixes: a crash when the only quantifiable sample is at time 0; a blank period cell silently splitting a profile (now an input error); `Test`/`Ref` labels producing non-canonical sequences; subjects with a missing period left in the BE input; column options being case-sensitive (`--period-column Period` was silently ignored); `inf`/`nan` accepted and producing invalid JSON; `ND` treated as missing and a measured 0 not counted as BLQ; replicate subjects with pre-dose > 5% blocking the BE analysis; Williams designs losing subjects; parallel studies without a period column; unnamed empty values; the wrong citation "ICH M13A 2.1.5" for baseline correction (it is 3.1); and the statement that ICH M13A has no emesis rule. A second batch added the protocol median tmax for the emesis rule, the analysis population for summary statistics, the corrected 72 h quantifiability check, a count of profiles without kel in the coverage rule, a PKNCA comparison script that fails on any undocumented difference, and true steady-state validation data with CtauSS, CavSS, swing and fluctuation compared.

**Pre-launch review of 2026-09-26.** Six reviews (engines, regulatory citations, UI, NCA → BE integration across all designs, validation and documentation, public release) led to version 1.0.0; the changes are listed in the repository [README](../README.md#changes-in-100-pre-launch-review) and tested in `test_launch_review.py` (16 tests: 10 for the first round and 6 for the independent second review, each failing on the version before its fix).

Not implemented: intravenous routes, manual kel windows, mean curves on nominal times, urinary data. **Comparison with PKNCA 0.12.1** ([`validation/`](validation/README.md)): on 348 single-dose and 60 steady-state profiles, with both trapezoidal rules, Cmax, tmax, AUC(0-t), kel, t½, kel points, adjusted r², AUC(0-∞), AUC(0-τ), Cav,ss and (with the trough sample at τ) Cτ,ss, Cmin,ss, swing and fluctuation match to 4 × 10⁻¹⁴. The only differences are partial AUCs ending after Clast, all explained: PKNCA returns NA when the end time falls between two BLQ samples, and integrates the post-Clast segment log-linearly even under the linear rule. PKNCA's default BLQ handling drops BLQ values between quantifiable ones, unlike ICH M13A; with that default AUC(0-t) changed in 33 of 36 affected profiles, by up to 21.1%. No comparison with Phoenix WinNonlin has been made.
