# ABEL against the EMA reference datasets · التحقق من ABEL على بيانات EMA المرجعية

Independent check of the EMA ABEL analysis in `engine/bioequivalence.py` (`--scaling abel --analysis ema`), run on 2026-09-26.

## Sources

- **EMA reference datasets.** The EMA published two datasets with its questions and answers on ABEL: dataset I (Annex II, full replicate TRTR/RTRT, 77 subjects, incomplete) and dataset II (Annex III, partial replicate TRR/RTR/RRT, 24 subjects, complete), with results computed in SAS Proc GLM (Method A). They are distributed, with the published results, as `rds01` and `rds02` in the R package **replicateBE 1.1.3** (Labes & Schütz, CRAN, GPL-3), which also implements Method A.
- **Other datasets.** replicateBE ships 28 further reference datasets (public-domain, edited or simulated). They are used here only through replicateBE's own Method A results; they are not copied into this repository (`run_replicatebe.R` regenerates them). The published EMA values were taken from the replicateBE documentation, not re-read from the EMA annexes.

## Reproduce

```bash
Rscript run_replicatebe.R .     # needs install.packages("replicateBE"); writes rdsNN.csv and replicatebe_method_a.csv
python compare_abel.py          # writes abel_comparison.json; exit 0 only if everything matches
```

`engine/test_ema_reference.py` keeps the dataset II result as a unit test.

## Results

| Dataset | Design | Result |
|---|---|---|
| **EMA dataset II (rds02)** | TRR/RTR/RRT, 24 subjects | **Matches the EMA-published result** — CVwR 11.2%, PE 102.26%, 90% CI 97.32–107.46%, no widening — and replicateBE to 6 × 10⁻¹⁵ |
| **EMA dataset I (rds01)** | TRTR/RTRT, 77 subjects, 8 with missing periods | **Rejected**: the engine analyses complete designs only. On the 69 complete subjects it matches replicateBE to 3 × 10⁻¹⁵ (CVwR 47.57%, PE 115.46%, CI 106.49–125.19%, limits 70.94–140.96%), but that is not the case EMA published (CVwR 47.0%, PE 115.66%, CI 107.11–124.89%) |
| 16 more datasets or complete-subject subsets (partial and full replicate) | 7 to 360 subjects | CVwR, widened limits, PE and 90% CI match replicateBE to 7.5 × 10⁻¹⁵; the pass/fail decision matches in every case (including CVwR above the 50% cap, e.g. 231%, and the 80.00–125.00% point-estimate constraint) |
| 28 datasets with other designs (TRT/RTR, TRRT/RTTR, TRR/RTT, 4-sequence designs, …) or missing data | — | Rejected explicitly (exit 2), as documented |

## What this shows, and what it does not

- The ABEL computation (reference-only ANOVA for CVwR, widened limits with k = 0.760 and the 50% cap, fixed-effects Method A point estimate and CI, decision rules) agrees with replicateBE to machine precision on every complete dataset in a supported design, and with the published EMA result for dataset II.
- **Gap:** the engine cannot analyse EMA dataset I as published, because incomplete replicate data are not supported. Method A is an ordinary fixed-effects model and handles missing periods naturally, so supporting incomplete data for the EMA analysis is the natural next step; until then, dataset I is verified only on its complete subjects.
- Method B (mixed model) and FDA RSABE are not covered by this check.

---

**بالعربية:** قورن تحليل ABEL (طريقة EMA A) في `bioequivalence.py` ببيانات EMA المرجعية وبحزمة replicateBE. المجموعة الثانية لـ EMA (تصميم جزئي، 24 مشاركًا) تطابق النتيجة التي نشرتها EMA تمامًا، وتطابق replicateBE حتى 6 × 10⁻¹⁵. و16 مجموعة أخرى مدعومة التصميم تطابق replicateBE حتى 7.5 × 10⁻¹⁵ مع القرار نفسه في كل الحالات. أما المجموعة الأولى لـ EMA (77 مشاركًا، ثمانية منهم تنقصهم فترات) فيرفضها المحرك لأنه لا يقبل البيانات الناقصة؛ تطابق على المشاركين المكتملين فقط، وهذا ليس ما نشرته EMA. دعم البيانات الناقصة في تحليل EMA هو الخطوة التالية الطبيعية.
