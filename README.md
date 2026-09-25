# Bioequivalence Lab · مختبر التكافؤ الحيوي

A bilingual (Arabic / English) web interface for average bioequivalence analysis — ABE, EMA ABEL, FDA RSABE, and power / sample size.
The page runs the Python engine (`engine/bioequivalence.py`) **unmodified in the browser** via [Pyodide](https://pyodide.org). No data leaves the user's device.

> ⚠️ **Research implementation, not regulatory certification.** Complete canonical designs only. It does not replace a statistical analysis plan or product-specific guidance. See [`engine/README.md`](engine/README.md) for model assumptions, unsupported cases and references.

## Features

- **Designs:** 2×2 (TR/RT), partial replicate (TRR/RTR/RRT), full replicate (TRTR/RTRT), auto-detected replicate, parallel.
- **Models:** equal-sequence subject contrasts (default) or EMA fixed effects (subject + period + treatment); Welch for parallel.
- **Limits:** 80.00–125.00%, NTI 90.00–111.11% (fixed limits only), or custom; optional half-up CI rounding.
- **Highly variable drugs:** EMA ABEL (Cmax, capped at CV 50%, point-estimate constraint) and FDA RSABE (Howe bound, sWR ≥ 0.294).
- **Power & sample size:** fixed-limit TOST power by numerical integration, with a power-vs-N curve.
- **Output:** verdict with exit code (0 pass · 1 not met · 2 invalid/unsupported), CI chart, scaled-criteria cards, JSON / CSV download.
- Arabic (RTL) and English UI, light and dark themes; engine messages are translated with the English original preserved.

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

**التشغيل:** انشر المجلد كموقع ثابت (Vercel أو GitHub Pages)، أو شغّل `run.bat` محليًا وافتح http://localhost:8765. أول زيارة تحمّل Python وNumPy وSciPy (نحو 20 ميغابايت).

**تنسيق البيانات:** الأعمدة `subject, sequence, period, treatment, value`. القيم موجبة، ولا حذف أو تعويض صامت للبيانات. ملفات للتجربة في `test_data` و`test_cases` (كلها بيانات مُحاكاة).

## License

MIT — see [LICENSE](LICENSE). The original engine is © 2025 K-Dense Inc. ([scientific-agent-skills](https://github.com/K-Dense-AI/scientific-agent-skills)), used under the MIT License.
