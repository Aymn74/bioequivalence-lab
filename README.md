# Bioequivalence Lab · مختبر التكافؤ الحيوي

[![Live demo](https://img.shields.io/badge/demo-live-0f766e?style=flat-square&logo=vercel&logoColor=white)](https://bioequivalence-lab.vercel.app)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue?style=flat-square)](LICENSE)
![Python](https://img.shields.io/badge/python-3.11-3776AB?style=flat-square&logo=python&logoColor=white)
![Pyodide](https://img.shields.io/badge/runs%20in%20browser-Pyodide%200.25-654FF0?style=flat-square&logo=webassembly&logoColor=white)
![NumPy](https://img.shields.io/badge/NumPy-013243?style=flat-square&logo=numpy&logoColor=white)
![SciPy](https://img.shields.io/badge/SciPy-8CAAE6?style=flat-square&logo=scipy&logoColor=white)
![Languages](https://img.shields.io/badge/UI-%D8%A7%D9%84%D8%B9%D8%B1%D8%A8%D9%8A%D8%A9%20%7C%20English-b7791f?style=flat-square)
![Privacy](https://img.shields.io/badge/data-stays%20in%20browser-177245?style=flat-square)
![Status](https://img.shields.io/badge/status-research%20use%20only-orange?style=flat-square)

A bilingual (Arabic / English) web interface for average bioequivalence analysis — ABE, EMA ABEL, FDA RSABE, and power / sample size.
The page runs the Python engine (`engine/bioequivalence.py`) **unmodified in the browser** via [Pyodide](https://pyodide.org). No data leaves the user's device.

> ⚠️ **Research implementation, not regulatory certification.** Complete canonical designs only. It does not replace a statistical analysis plan or product-specific guidance. See [`engine/README.md`](engine/README.md) for model assumptions, unsupported cases and references.

**Live app:** https://bioequivalence-lab.vercel.app

## Screenshots

| Arabic interface — example study, all criteria met | Highly variable drug — EMA ABEL + FDA RSABE |
|---|---|
| ![Arabic results](docs/screenshots/02-results-ar.png) | ![ABEL and RSABE](docs/screenshots/03-hvd-abel-rsabe-en.png) |
| **Same study without scaling — ABE not met (exit 1)** | **Unsupported case rejected explicitly (exit 2)** |
| ![Not met](docs/screenshots/04-not-met-ar.png) | ![Rejected input](docs/screenshots/05-rejected-input-ar.png) |

**Power & sample size** — partial replicate, CV 30%, GMR 0.95 → N = 30 (81.0%)

![Power curve](docs/screenshots/06-power-en.png)

<details>
<summary>Full page and mobile views</summary>

![Full page, Arabic](docs/screenshots/01-overview-ar.png)

<img src="docs/screenshots/07-mobile-ar.png" width="320" alt="Mobile, Arabic">

</details>

All screenshots use simulated data from `engine/example_partial.csv`, `test_data/` and `test_cases/`.

## Features

- **Designs:** 2×2 (TR/RT), partial replicate (TRR/RTR/RRT), full replicate (TRTR/RTRT), auto-detected replicate, parallel.
- **Models:** equal-sequence subject contrasts (default) or EMA fixed effects (subject + period + treatment); Welch for parallel.
- **Limits:** 80.00–125.00%, NTI 90.00–111.11% (fixed limits only), or custom; optional half-up CI rounding.
- **Highly variable drugs:** EMA ABEL (Cmax, capped at CV 50%, point-estimate constraint) and FDA RSABE (Howe bound, sWR ≥ 0.294).
- **Power & sample size:** fixed-limit TOST power by numerical integration, with a power-vs-N curve.
- **Output:** verdict with exit code (0 pass · 1 not met · 2 invalid/unsupported), CI chart, scaled-criteria cards, JSON / CSV download.
- Arabic (RTL) and English UI, light and dark themes; engine messages are translated with the English original preserved.

## SFDA / GCC profile

Selecting **SFDA / GCC** in the UI (CLI: `--profile sfda`) applies the GCC bioequivalence guideline **DS-G-010 V3.1** (adopted by the SFDA), section by section:

| Guideline rule | What the profile does |
|---|---|
| ANOVA with fixed effects: sequence, subject within sequence, period, formulation | Forces the fixed-effects model; rejects `--analysis contrast` and Welch |
| CI bounds ≥ 80.00% and ≤ 125.00% **after rounding to two decimals** | Rounding always on |
| Widening (ABEL, k = 0.760, max 69.84–143.19%) for **Cmax only**; AUC stays 80.00–125.00% | ABEL allowed for Cmax only, with prospective justification |
| No reference-scaled (RSABE) approach | RSABE / "both" rejected |
| Minimum **18 evaluable subjects** | Fewer than 18 → finding, exit code 1; power solver starts at 18 |
| Pre-dose concentration **> 5% of Cmax** → exclude that subject-period (not for endogenous substances) | Optional `predose` column; 2×2 / parallel: subject removed and listed; replicate: explicit error (incomplete design unsupported); `--endogenous` skips the rule |
| Two-stage design: adjusted CI (e.g. 94.12%) and a stage term | `--ci-level 94.12`; an optional `stage` column nests periods within stage |
| ANOVA tables must be submitted | Full ANOVA table (SS, df, MS, F, p) and intra-subject CV in the output |

![SFDA results](docs/screenshots/08-sfda-results-ar.png)

Not implemented: endogenous baseline correction, dissolution similarity f2, product-specific guidance lookup, and multiple-reference comparisons (analyse each T/R pair as a separate file, as the guideline requires).

**Verification of the ANOVA:** estimate, 90% CI and formulation F-test match `statsmodels` OLS to ~1e-15 on balanced and unbalanced 2×2, partial and full replicate data; the subject(sequence) sum of squares matches the textbook between-subject formula and an explicit nested regression (statsmodels' Type I output is incorrect for this singular parameterisation). Tests: `engine/test_sfda.py` (13 tests) plus the original 12.

**Power note:** for a partial replicate with N = 24, CV 30%, GMR 0.95, power is 70.8% with the contrast model (df = N − 3) and 72.5% with the fixed-effects ANOVA model the SFDA profile uses (df = (N − 1)(P − 1) − 1). Both are correct for their model.

## Run

**Online:** deploy the repository root as a static site (Vercel, GitHub Pages, Cloudflare Pages). No build step is needed — `index.html` is self-contained.

**Locally (Windows):** double-click `run.bat`, or:

```bash
python -m http.server 8765
```

Then open http://localhost:8765. The first visit downloads Pyodide + NumPy + SciPy from jsDelivr (~20 MB, then cached).

## Data format

```csv
subject,sequence,period,treatment,value
S01,TR,1,T,933.70
S01,TR,2,R,1079.54
```

Column names are case-insensitive. `subject`, `treatment` and `value` can be remapped in the UI; `sequence` and `period` must use these names. Treatment labels `T`/`R`, `Test`/`Ref`/`Reference` are accepted. Values must be positive; the engine never silently drops, imputes or averages records.

Sample files: [`test_data/`](test_data) (two simulated studies) and [`test_cases/`](test_cases) (edge cases; files prefixed `ERR_` must be rejected). All data are simulated — no real clinical study.

## Repository layout

| Path | Contents |
|---|---|
| `index.html` | Built, self-contained app (deploy this) |
| `src/template.html`, `src/build.js` | UI source; `node src/build.js` regenerates `index.html` after editing the template or engine |
| `engine/` | Corrected Python engine, tests, simulation validation, provenance and the original upstream script |
| `test_data/`, `test_cases/` | Simulated CSV files for trying the app |

## Engine verification

```bash
cd engine
python -m pip install -r requirements.txt
python -m unittest -v
python validate_simulation.py
```

The web UI was checked against the engine's reference output (`engine/example_results.json`) and against an independent 2×2 computation (identical to 4 decimals).

---

## بالعربية

واجهة ويب ثنائية اللغة لتحليل التكافؤ الحيوي: ABE وEMA ABEL وFDA RSABE، وحساب القوة وحجم العينة. تشغّل الواجهة ملف `engine/bioequivalence.py` كما هو داخل المتصفح عبر Pyodide، ولا تُرسل البيانات إلى أي خادم.

> ⚠️ **تنفيذ بحثي وليس اعتمادًا تنظيميًا.** يقبل التصاميم الكاملة القياسية فقط، ولا يغني عن خطة التحليل الإحصائي أو الإرشادات الخاصة بالمنتج. راجع [`engine/اقرأني.md`](engine/اقرأني.md).

**الرابط المباشر:** https://bioequivalence-lab.vercel.app — لقطات الشاشة في الأعلى، وجميعها ببيانات مُحاكاة.

**إطار SFDA / الخليجي:** اختيار «SFDA / الخليجي» في الواجهة (أو `--profile sfda`) يطبّق دليل التكافؤ الحيوي الخليجي DS-G-010 V3.1: تحليل ANOVA بتأثيرات ثابتة، وتقريب حدود فترة الثقة لمنزلتين، وتوسيع الحدود (ABEL) لـ Cmax فقط، وعدم استخدام RSABE، وحد أدنى 18 مشاركًا قابلًا للتقييم، واستبعاد الفترة التي يتجاوز فيها تركيز ما قبل الجرعة 5% من Cmax (عمود `predose` اختياري)، وفترة ثقة معدّلة للتصميم على مرحلتين (مثل 94.12%) مع عمود `stage`، وإخراج جدول ANOVA كاملًا. غير منفّذ: تصحيح خط الأساس، ومعامل f2، وقاعدة بيانات الأدلة الخاصة بالمنتجات، والمقارنة مع أكثر من مرجع.

**التشغيل:** انشر المجلد كموقع ثابت (Vercel أو GitHub Pages)، أو شغّل `run.bat` محليًا وافتح http://localhost:8765. أول زيارة تحمّل Python وNumPy وSciPy (نحو 20 ميغابايت).

**تنسيق البيانات:** الأعمدة `subject, sequence, period, treatment, value`. القيم موجبة، ولا حذف أو تعويض صامت للبيانات. ملفات للتجربة في `test_data` و`test_cases` (كلها بيانات مُحاكاة).

## License

MIT — see [LICENSE](LICENSE). The original engine is © 2025 K-Dense Inc. ([scientific-agent-skills](https://github.com/K-Dense-AI/scientific-agent-skills)), used under the MIT License.
