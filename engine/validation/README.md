# nca.py against PKNCA · مقارنة nca.py مع PKNCA

Independent check of `engine/nca.py` against **PKNCA 0.12.1** (R 4.6.1, CRAN), run on 2026-09-26.

## Reproduce

```bash
python make_pknca_data.py          # writes single_dose.csv, steady_state.csv
Rscript run_pknca.R .              # writes pknca_*.csv (needs install.packages("PKNCA"))
python compare_pknca.py            # writes pknca_comparison.json, prints the summary
```

## Data

- `single_dose.csv`: 348 profiles — the 48 profiles of `engine/example_nca.csv` (missed sample removed) and 300 random one-compartment profiles with actual times and about a quarter of values BLQ (written as 0).
- `steady_state.csv`: 60 random steady-state profiles, tau = 12 h, last sample within ±0.1 h of 12 h.

## PKNCA settings

The installed PKNCA 0.12.1 defaults were read before the run: `adj.r.squared.factor = 1e-4`, `allow.tmax.in.half.life = FALSE`, `min.hl.points = 3`, `first.tmax = TRUE`. The body of `pk.calc.half.life` in the installed version was checked: best adjusted r² over all windows, then windows with lambda.z > 0 and adjusted r² > best − 1e-4, then the most points; zero concentrations are left out of the fit. `nca.py` implements the same rule.

One setting was changed: PKNCA's default `conc.blq` drops BLQ values between quantifiable ones (`middle = "drop"`). ICH M13A 2.2.2.2 treats BLQ as zero, so the comparison uses `middle = "keep"`. Both trapezoidal rules were run (`linear`, `lin up/log down`).

## Results (tolerance: relative difference 1e-9)

| Parameters | Profiles | Largest relative difference | Mismatches |
|---|---|---|---|
| Cmax, tmax, tlast, Clast, AUC(0-t), kel, t½, kel points, adjusted r², AUC(0-∞) obs and pred, AUC(0-2) | 348 × 2 rules | 4 × 10⁻¹⁴ | 0 |
| Steady state: Cmax,ss, tmax,ss, Cmin,ss, AUC(0-τ) | 60 × 2 rules | 2 × 10⁻¹⁴ | 0 |
| AUC(0-24) against `aucint.all` | 348 × 2 rules | 5 × 10⁻¹⁵ where both give a value | linear 15, lin-up/log-down 8 |

Every AUC(0-24) mismatch has 24 h after Clast, and all are explained exactly:

- **8 profiles (both rules): PKNCA returns NA.** 24 h lies between two BLQ samples after Clast. The area after Clast is zero there, and `nca.py` returns AUC(0-t); PKNCA cannot interpolate between two post-Clast BLQ values.
- **7 profiles (linear rule only): different integration after Clast.** PKNCA integrates the segment from Clast to the interpolated 24 h value with the log trapezoid even when the linear rule is selected (`interval_method <- "log"` for times after tlast in `pk.calc.aucint`). `nca.py` applies the pre-specified rule to every segment (ICH M13A 2.2.2.2 asks for the method to be stated). Replacing that one segment with the log formula reproduces PKNCA to 0 at 12 decimals in all 7.
- Against `aucint.last` (all concentrations after Clast set to zero), 28 profiles differ, as expected: `nca.py` keeps a later BLQ sample as zero at its own time, like `aucint.all`.

These partial-AUC cases do not arise for AUC(0-t), AUC(0-∞) or AUC(0-τ), which match exactly.

## Effect of PKNCA's default BLQ rule

36 of the 348 profiles have a BLQ value between quantifiable ones. With PKNCA's default (`middle = "drop"`), AUC(0-t) changes in 33 of them, by up to 21.1%. A study analysed with PKNCA defaults would therefore not follow ICH M13A 2.2.2.2 for these profiles unless `conc.blq` is set.

## Not covered

AUC(0-72h) (the random profiles end before 72 h), baseline correction, the emesis and coverage rules (rule logic, tested in `test_nca.py`), and Phoenix WinNonlin.

---

**بالعربية:** قورن `nca.py` بحزمة PKNCA 0.12.1 على 348 منحنى جرعة مفردة و60 منحنى حالة مستقرة، بطريقتي شبه المنحرف. تطابقت كل المعايير الأساسية حتى 4 × 10⁻¹⁴: Cmax وtmax وAUC(0-t) وkel وt½ وعدد نقاط kel وR² المعدّل وAUC(0-∞) وAUC(0-τ) وCmin,ss. الفروق الوحيدة في المساحة الجزئية AUC(0-24) حين تقع 24 ساعة بعد آخر تركيز مقاس، وكلها مفسَّرة: PKNCA يعيد NA حين تقع 24 ساعة بين عينتين تحت حد القياس، ويحسب المقطع بعد آخر تركيز مقاس بالطريقة اللوغاريتمية حتى لو اختيرت الخطية، بينما nca.py يلتزم بالطريقة المحددة مسبقًا. تنبيه: إعداد PKNCA الافتراضي يحذف القيمة تحت حد القياس الواقعة بين قيمتين مقاستين بدل عدّها صفرًا كما في ICH M13A، وهذا غيّر AUC(0-t) في 33 منحنى بفرق يصل إلى 21.1%.
