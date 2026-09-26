# ABEL against the EMA reference datasets · التحقق من ABEL على بيانات EMA المرجعية

Independent check of the EMA ABEL analysis in `engine/bioequivalence.py` (`--scaling abel --analysis ema`), run on 2026-09-26.

## Sources

- **EMA reference datasets.** The EMA published two datasets with its questions and answers on ABEL: dataset I (Annex II, full replicate TRTR/RTRT, 77 subjects, incomplete) and dataset II (Annex III, partial replicate TRR/RTR/RRT, 24 subjects, complete), with results computed in SAS Proc GLM (Method A). They are distributed, with the published results, as `rds01` and `rds02` in the R package **replicateBE 1.1.3** (Schütz H, CRAN, GPL ≥ 3, doi:10.32614/CRAN.package.replicateBE; D. Labes and M. Tomashevskiy are contributors), which also implements Method A. `rds01.csv` and `rds02.csv` in this folder are these EMA datasets (Annex II and Annex III of EMA/582648/2016, 21 September 2016, as cited in the replicateBE documentation), exported from replicateBE by `run_replicatebe.R`; they are EMA data, reproduced with the source acknowledged, not replicateBE code. `replicatebe_method_a.csv` holds replicateBE results (program output).
- **Other datasets.** replicateBE ships 28 further reference datasets (public-domain, edited or simulated). They are used here only through replicateBE's own Method A results; they are not copied into this repository (`run_replicatebe.R` regenerates them). The published EMA values were taken from the replicateBE documentation, not re-read from the EMA annexes.

## Reproduce

```bash
Rscript run_replicatebe.R .     # needs install.packages("replicateBE"); writes rdsNN.csv and replicatebe_method_a.csv
python compare_abel.py          # writes abel_comparison.json; exit 0 only if everything matches
```

`engine/test_ema_reference.py` keeps the dataset I and II results as unit tests.

## Results

| Dataset | Design | Result |
|---|---|---|
| **EMA dataset II (rds02)** | TRR/RTR/RRT, 24 subjects | **Matches the EMA-published result** — CVwR 11.2%, PE 102.26%, 90% CI 97.32–107.46%, no widening — and replicateBE to 6 × 10⁻¹⁵ |
| **EMA dataset I (rds01)** | TRTR/RTRT, 77 subjects, 8 with missing periods | **Matches the EMA-published result** — CVwR 47.0% (46.96%), PE 115.66%, 90% CI 107.11–124.89%, widened limits 71.23–140.40% — and replicateBE to 4 × 10⁻¹⁵ |
| 29 more datasets and complete-subject subsets (partial and full replicate; 12 incomplete, including missing values coded `NA`) | 7 to 360 subjects | CVwR, widened limits, PE and 90% CI match replicateBE to 7.5 × 10⁻¹⁵; the pass/fail decision matches in every case (including CVwR above the 50% cap, e.g. 231%, and the 80.00–125.00% point-estimate constraint) |
| 15 datasets with other designs (TRT/RTR, TRRT/RTTR, TRR/RTT, TRR/RTR, 4-sequence designs, …) | — | Rejected explicitly (exit 2), as documented |

## What this shows, and what it does not

- The ABEL computation (reference-only ANOVA for CVwR, widened limits with k = 0.760 and the 50% cap, fixed-effects Method A point estimate and CI, decision rules) agrees with replicateBE to machine precision on every complete dataset in a supported design, and with the published EMA result for dataset II.
- Incomplete data: the fixed-effects analysis (EMA Method A) uses all available observations, like SAS Proc GLM and replicateBE; EMA dataset I and 12 other incomplete datasets reproduce the reference results. The contrast analysis, RSABE and the FDA NTI method still require complete subjects.
- Method B (mixed model) and FDA RSABE are not covered by this check.

---

**بالعربية:** قورن تحليل ABEL (طريقة EMA A) في `bioequivalence.py` ببيانات EMA المرجعية وبحزمة replicateBE. المجموعتان المرجعيتان لـ EMA تطابقان النتائج التي نشرتها EMA تمامًا: الأولى (تصميم كامل، 77 مشاركًا، ثمانية منهم تنقصهم فترات) والثانية (تصميم جزئي، 24 مشاركًا). و31 مجموعة مدعومة التصميم، 13 منها ببيانات ناقصة، تطابق replicateBE حتى 7.5 × 10⁻¹⁵ مع القرار نفسه في كل الحالات. البيانات الناقصة مدعومة في نموذج التأثيرات الثابتة وABEL وإطار SFDA، أما تحليل فروق الأفراد وRSABE وطريقة FDA للأدوية ضيقة المؤشر فتتطلب مشاركين مكتملين.
