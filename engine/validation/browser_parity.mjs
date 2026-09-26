// Browser engine vs local Python: runs the same command lines through the page's Pyodide worker and
// through the local Python engines, then compares every value of the JSON outputs.
//
// Usage: serve the repository root (python -m http.server 8766) or use the deployed site, then
//   node engine/validation/browser_parity.mjs [base-url] [python] [path-to-chrome]
// Defaults: http://localhost:8766, "python", Chrome at its usual Windows path. Needs Node >= 22.
// Writes engine/validation/browser_parity.json and exits 1 on any difference beyond the tolerance.
//
// Numbers are compared as relative differences (tolerance 1e-8: the browser uses Pyodide's own NumPy/SciPy,
// the local run whatever is installed, and distribution quantiles such as chi2.ppf differ in the last digits
// between SciPy versions; the largest difference seen, 4.6e-9, is the RSABE upper bound near zero).
// Strings, booleans, nulls, table and key structure and exit codes must be identical, so every decision
// is. The fingerprint is the SHA-256 of every NCA number printed to 10 significant digits, in document
// order, computed separately for each side (NCA uses no distribution quantiles, so it must be equal).
import { spawn, spawnSync } from 'node:child_process';
import { mkdtempSync, readFileSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';
import { createHash } from 'node:crypto';

const BASE = process.argv[2] || 'http://localhost:8766';
const PYTHON = process.argv[3] || 'python';
const CHROME = process.argv[4] || 'C:/Program Files/Google/Chrome/Application/chrome.exe';
const HERE = dirname(fileURLToPath(import.meta.url));
const ENGINE = join(HERE, '..');
const ROOT = join(ENGINE, '..');
const TOL = 1e-8;
const PORT = 9334;
const sleep = ms => new Promise(r => setTimeout(r, ms));

// [engine, input file (repo-relative, or null), arguments]
const CASES = [
  ['nca', 'engine/example_nca.csv', ['--profile', 'sfda']],
  ['nca', 'engine/example_nca.csv', ['--auc-method', 'linup-logdown', '--partial-auc', '0-2', '--lloq', '0.5']],
  ['nca', 'engine/validation/single_dose.csv', ['--subject-column', 'id']],
  ['nca', 'engine/validation/single_dose.csv', ['--subject-column', 'id', '--auc-method', 'linup-logdown', '--truncate-72']],
  ['nca', 'engine/validation/steady_state.csv', ['--subject-column', 'id', '--tau', '12']],
  ['be', 'engine/example_partial.csv', ['--design', 'partial', '--metric', 'cmax', '--scaling', 'both', '--abel-justified']],
  ['be', 'engine/example_partial.csv', ['--design', 'partial', '--metric', 'cmax', '--profile', 'sfda', '--scaling', 'abel', '--abel-justified']],
  ['be', 'engine/example_partial.csv', ['--design', 'partial', '--metric', 'cmax', '--scaling', 'both', '--abel-justified',
    '--potency-test', '94', '--potency-reference', '101', '--potency-correction']],
  ['be', 'test_data/study_2x2_AUC.csv', ['--metric', 'auc', '--anova']],
  ['be', 'test_data/study_2x2_AUC.csv', ['--metric', 'auc', '--analysis', 'ema', '--anova', '--ema-rounding']],
  ['be', 'test_data/study_endogenous_2x2_Cmax_baseline.csv', ['--metric', 'cmax', '--endogenous', '--profile', 'sfda']],
  ['be', 'test_data/study_full_replicate_Cmax.csv', ['--design', 'full', '--metric', 'cmax', '--scaling', 'both', '--abel-justified']],
  ['be', 'test_data/study_multigroup_2x2_Cmax.csv', ['--metric', 'cmax', '--analysis', 'ema', '--anova']],
  ['be', 'test_data/study_nti_full_replicate_AUC.csv', ['--design', 'full', '--metric', 'auc', '--scaling', 'fda-nti']],
  ['be', 'test_data/study_sfda_2x2_Cmax_predose.csv', ['--metric', 'cmax', '--profile', 'sfda']],
  ['be', 'test_data/study_williams_T_vs_GCC_US.csv', ['--design', 'multi', '--test-label', 'T', '--reference-label', 'R_GCC', '--metric', 'cmax']],
  ['be', 'test_data/study_williams_T_vs_GCC_US.csv', ['--design', 'multi', '--test-label', 'T', '--reference-label', 'R_US', '--metric', 'cmax', '--profile', 'sfda']],
  ['be', 'test_cases/04_parallel.csv', ['--design', 'parallel', '--welch']],
  ['be', null, ['--power', '--design', 'partial', '--cv', '0.3', '--gmr', '0.95']],
  ['be', null, ['--power', '--design', 'full', '--cv', '0.35', '--gmr', '0.95', '--n', '24', '--profile', 'sfda']],
];

// ---------------------------------------------------------------- comparison
function walk(a, b, path, acc) {
  const ta = a === null ? 'null' : Array.isArray(a) ? 'array' : typeof a;
  const tb = b === null ? 'null' : Array.isArray(b) ? 'array' : typeof b;
  if (ta !== tb) { acc.mismatch.push(`${path}: ${ta} vs ${tb}`); return; }
  if (ta === 'number') {
    acc.n++; acc.fa.push(a.toPrecision(10)); acc.fb.push(b.toPrecision(10));
    const d = a === b ? 0 : Math.abs(a - b) / Math.max(Math.abs(a), Math.abs(b), 1e-300);
    if (d > acc.max) { acc.max = d; acc.where = path; }
    if (d > TOL) acc.mismatch.push(`${path}: ${a} vs ${b}`);
  } else if (ta === 'array') {
    if (a.length !== b.length) { acc.mismatch.push(`${path}: length ${a.length} vs ${b.length}`); return; }
    a.forEach((x, i) => walk(x, b[i], `${path}[${i}]`, acc));
  } else if (ta === 'object') {
    const ka = Object.keys(a), kb = Object.keys(b);
    if (ka.join('|') !== kb.join('|')) { acc.mismatch.push(`${path}: keys differ`); return; }
    ka.forEach(k => walk(a[k], b[k], `${path}.${k}`, acc));
  } else if (a !== b) acc.mismatch.push(`${path}: ${JSON.stringify(a)} vs ${JSON.stringify(b)}`);
}
const sha = xs => createHash('sha256').update(xs.join('\n')).digest('hex');

// ---------------------------------------------------------------- local Python
function runPython(engine, file, args) {
  const script = engine === 'nca' ? 'nca.py' : 'bioequivalence.py';
  const argv = [script, ...(file ? ['-i', join(ROOT, file)] : []), ...args, '--format', 'json'];
  const r = spawnSync(PYTHON, argv, { cwd: ENGINE, encoding: 'utf8', maxBuffer: 1 << 28 });
  if (r.error) throw r.error;
  return { code: r.status, out: r.stdout.trim() ? JSON.parse(r.stdout) : null, err: r.stderr };
}

// ---------------------------------------------------------------- browser (headless Chrome over DevTools)
const chrome = spawn(CHROME, ['--headless=new', `--remote-debugging-port=${PORT}`, `--user-data-dir=${mkdtempSync(join(tmpdir(), 'cdp-'))}`,
  '--no-first-run', 'about:blank'], { stdio: 'ignore' });
let ws, seq = 0; const pending = {};
async function connect() {
  let url;
  for (let i = 0; i < 100 && !url; i++) {
    try {
      url = (await (await fetch(`http://127.0.0.1:${PORT}/json`)).json()).find(t => t.type === 'page')?.webSocketDebuggerUrl;
    } catch { /* not up yet */ }
    if (!url) await sleep(200);
  }
  if (!url) throw new Error('headless Chrome did not start');
  ws = new WebSocket(url);
  await new Promise((res, rej) => { ws.onopen = res; ws.onerror = rej; });
  ws.onmessage = e => { const m = JSON.parse(e.data); if (m.id && pending[m.id]) { pending[m.id](m); delete pending[m.id]; } };
}
const send = (method, params = {}) => new Promise((res, rej) => {
  const id = ++seq; pending[id] = m => m.error ? rej(new Error(`${method}: ${m.error.message}`)) : res(m.result);
  ws.send(JSON.stringify({ id, method, params }));
});
async function evaluate(expr) {
  const r = await send('Runtime.evaluate', { expression: `(async()=>{${expr}})()`, awaitPromise: true, returnByValue: true });
  if (r.exceptionDetails) throw new Error(r.exceptionDetails.exception?.description || r.exceptionDetails.text);
  return r.result.value;
}
async function runBrowser(engine, file, args) {
  const text = file ? readFileSync(join(ROOT, file), 'utf8') : '';
  const ext = file && /\.(tsv|tab)$/i.test(file) ? '.tsv' : '.csv';
  const r = await evaluate(`
    const text=${JSON.stringify(text)},path='/be/parity${ext}';
    if(text) await call('write',{path,text});
    return await call(${JSON.stringify(engine === 'nca' ? 'nca' : 'cli')}, [...(text?['-i',path]:[]), ...${JSON.stringify(args)}, '--format','json']);`);
  return { code: r.code, out: r.stdout.trim() ? JSON.parse(r.stdout) : null, err: r.stderr };
}

// ---------------------------------------------------------------- main
let failed = false;
const results = [];
try {
  await connect();
  await send('Page.navigate', { url: `${BASE}/?parity=${Date.now()}` });
  const t0 = Date.now();
  while (!(await evaluate(`return typeof ready!=='undefined'&&ready===true`).catch(() => false))) {
    if (Date.now() - t0 > 180000) throw new Error('engine did not load in 3 min');
    await sleep(500);
  }
  let allA = [], allB = [], total = 0, ncaN = 0;
  for (const [engine, file, args] of CASES) {
    const py = runPython(engine, file, args), br = await runBrowser(engine, file, args);
    const acc = { n: 0, max: 0, where: '', mismatch: [], fa: [], fb: [] };
    if (py.code !== br.code) acc.mismatch.push(`exit code ${py.code} (python) vs ${br.code} (browser)`);
    walk(py.out, br.out, '$', acc);
    const ok = acc.mismatch.length === 0;
    failed ||= !ok;
    total += acc.n;
    if (engine === 'nca') { ncaN += acc.n; allA = allA.concat(acc.fa); allB = allB.concat(acc.fb); }
    const label = `${engine} ${file || '(no input)'} ${args.join(' ')}`;
    results.push({ case: label, exit_code: py.code, values: acc.n, max_rel_diff: acc.max, at: acc.where, ok, mismatches: acc.mismatch.slice(0, 10) });
    console.log(`${ok ? 'OK  ' : 'DIFF'} values=${String(acc.n).padStart(5)} max_rel=${acc.max.toExponential(1)} exit=${py.code}  ${label}`);
    if (!ok) acc.mismatch.slice(0, 5).forEach(m => console.log('     ', m));
  }
  const summary = { base: BASE, date: new Date().toISOString(), tolerance: TOL, cases: CASES.length, values: total,
    max_rel_diff: Math.max(...results.map(r => r.max_rel_diff)),
    nca_values: ncaN, nca_fingerprint_python: sha(allA), nca_fingerprint_browser: sha(allB), nca_fingerprints_equal: sha(allA) === sha(allB), results };
  failed ||= !summary.nca_fingerprints_equal;
  writeFileSync(join(HERE, 'browser_parity.json'), JSON.stringify(summary, null, 1) + '\n');
  console.log(`\n${CASES.length} cases, ${total} values, max relative difference ${summary.max_rel_diff.toExponential(1)}; ` +
    `NCA fingerprint over ${ncaN} values (10 significant digits) ${summary.nca_fingerprints_equal ? 'equal' : 'DIFFERS'}: ${summary.nca_fingerprint_python.slice(0, 16)}…`);
  console.log(failed ? 'FAILED: differences beyond the tolerance' : 'OK: the browser engine reproduces the local Python output');
} catch (e) {
  console.error(e); failed = true;
} finally {
  try { ws && ws.close(); } catch { /* ignore */ }
  chrome.kill();
}
process.exit(failed ? 1 : 0);
