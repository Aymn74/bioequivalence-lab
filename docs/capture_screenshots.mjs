// Regenerates docs/screenshots/*.png from the running app (simulated data only).
// Usage: serve the repository root (python -m http.server 8766), then
//   node docs/capture_screenshots.mjs [http://localhost:8766] [path-to-chrome]
// Drives headless Chrome over the DevTools protocol; needs Node >= 22 (built-in WebSocket).
import { spawn } from 'node:child_process';
import { mkdtempSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';

const BASE = process.argv[2] || 'http://localhost:8766';
const CHROME = process.argv[3] || 'C:/Program Files/Google/Chrome/Application/chrome.exe';
const OUT = join(dirname(fileURLToPath(import.meta.url)), 'screenshots');
const PORT = 9333;
const sleep = ms => new Promise(r => setTimeout(r, ms));

const chrome = spawn(CHROME, ['--headless=new', `--remote-debugging-port=${PORT}`, `--user-data-dir=${mkdtempSync(join(tmpdir(), 'cdp-'))}`,
  '--no-first-run', '--hide-scrollbars', '--font-render-hinting=none', 'about:blank'], { stdio: 'ignore' });

let ws, seq = 0; const pending = {};
async function connect() {
  for (let i = 0; i < 50; i++) {
    try {
      const tabs = await (await fetch(`http://127.0.0.1:${PORT}/json`)).json();
      const page = tabs.find(t => t.type === 'page');
      if (page) { ws = new WebSocket(page.webSocketDebuggerUrl); break; }
    } catch { /* not up yet */ }
    await sleep(200);
  }
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

// page helpers injected before every scenario
const HELPERS = `
window.W=async(c,ms=60000)=>{const t=Date.now();while(!c()){if(Date.now()-t>ms)throw new Error('timeout');await new Promise(r=>setTimeout(r,100));}};
window.runBE=async()=>{S.last=null;$('#runBtn').click();await W(()=>S.last&&!$('#runBtn').classList.contains('busy'));await new Promise(r=>setTimeout(r,300));};
window.runNCA=async()=>{$('#ncaRun').click();await new Promise(r=>setTimeout(r,100));await W(()=>!N.running);await new Promise(r=>setTimeout(r,300));};
window.tab=n=>document.querySelector('.tab[data-tab="'+n+'"]').click();
window.loadFile=async(path)=>{const t=await (await fetch(path)).text();setData(t,path.split('/').pop());};
`;

async function open({ lang = 'ar', theme = 'light', width = 1320, height = 900, mobile = false }) {
  await send('Emulation.setDeviceMetricsOverride', { width, height, deviceScaleFactor: 2, mobile });
  await send('Page.navigate', { url: `${BASE}/?shot=${Date.now()}` });
  await sleep(300);
  await evaluate(`localStorage.setItem('be_lang','${lang}');localStorage.setItem('be_theme','${theme}');`);
  await send('Page.reload', { ignoreCache: false });
  await sleep(500);
  await evaluate(`await new Promise(r=>{const f=()=>typeof $==='function'&&document.querySelector('#engDot.ok')?r():setTimeout(f,200);f();});` + HELPERS
    + `document.querySelector('#boot')&&(document.querySelector('#boot').style.display='none');`);
}

async function shot(name, selector) {
  await evaluate(`const t=document.querySelector('#toast');t.classList.remove('on');t.style.display='none';`);
  await sleep(400);
  let clip;
  if (selector) {
    clip = await evaluate(`const e=document.querySelector('${selector}');e.scrollIntoView();const r=e.getBoundingClientRect();
      return {x:r.left+scrollX,y:r.top+scrollY,width:r.width,height:r.height,scale:1};`);
  } else {
    const m = await send('Page.getLayoutMetrics');
    clip = { x: 0, y: 0, width: m.cssContentSize.width, height: m.cssContentSize.height, scale: 1 };
  }
  await evaluate('scrollTo(0,0);');
  const { data } = await send('Page.captureScreenshot', { format: 'png', clip, captureBeyondViewport: true });
  writeFileSync(join(OUT, name), Buffer.from(data, 'base64'));
  console.log('saved', name, Math.round(clip.width * 2), 'x', Math.round(clip.height * 2));
}

const EXAMPLE_BE = `tab('analysis');$('#loadExample').click();await new Promise(r=>setTimeout(r,200));`;

try {
  await connect();
  await send('Page.enable'); await send('Runtime.enable');

  // 01 full page, 02 result card: partial-replicate example, ABEL + RSABE
  await open({});
  await evaluate(EXAMPLE_BE + 'await runBE();');
  await shot('01-overview-ar.png');
  await shot('02-results-ar.png', '#resCard');

  // 03 same study, English
  await open({ lang: 'en' });
  await evaluate(EXAMPLE_BE + 'await runBE();');
  await shot('03-hvd-abel-rsabe-en.png', '#resCard');

  // 04 same study without scaling: ABE not met
  await open({});
  await evaluate(EXAMPLE_BE + `document.querySelector('#scalingSeg button[data-v="none"]').click();syncForm();await runBE();`);
  await shot('04-not-met-ar.png', '#resCard');

  // 05 unsupported case rejected (dark theme): RSABE with low reference variability
  await open({ theme: 'dark' });
  await evaluate(`tab('analysis');await loadFile('test_cases/03_partial_low_var.csv');
    document.querySelector('#scalingSeg button[data-v="rsabe"]').click();syncForm();await runBE();`);
  await shot('05-rejected-input-ar.png', '#resCard');

  // 06 power and sample size, English
  await open({ lang: 'en' });
  await evaluate(`tab('power');document.querySelector('#pDesignSeg button[data-v="partial"]').click();
    $('#pCV').value='30';$('#pGMR').value='0.95';$('#pTarget').value='80';$('#pRun').click();
    await W(()=>S.plast&&S.plast.curve&&!$('#pRun').classList.contains('busy'));`);
  await shot('06-power-en.png', '#pane-power');

  // 07 mobile, Arabic, NCA tab
  await open({ width: 390, height: 844, mobile: true });
  await evaluate(`tab('nca');document.querySelector('#toast').style.display='none';`);
  const { data } = await send('Page.captureScreenshot', { format: 'png' });
  writeFileSync(join(OUT, '07-mobile-ar.png'), Buffer.from(data, 'base64')); console.log('saved 07-mobile-ar.png');

  // 08 SFDA profile: 2x2 with a pre-dose exclusion
  await open({});
  await evaluate(`tab('analysis');await loadFile('test_data/study_sfda_2x2_Cmax_predose.csv');
    document.querySelector('#profileSeg button[data-v="sfda"]').click();$('#metric').value='cmax';syncForm();await runBE();`);
  await shot('08-sfda-results-ar.png', '#resCard');

  // 09 NCA results (SFDA profile, simulated 2x2 example), 10 hand-off to the BE analysis (AUC(0-t))
  await open({});
  await evaluate(`tab('nca');$('#ncaExample').click();$('#ncaMedT').value='1.5';await runNCA();`);
  await shot('09-nca-results-ar.png', '#ncaResCard');
  await evaluate(`$('#ncaBeMetric').value='auc_0_t';$('#ncaBeGo').click();await new Promise(r=>setTimeout(r,300));await runBE();`);
  await shot('10-nca-to-be-ar.png', '#pane-analysis');
} catch (e) {
  console.error(e); process.exitCode = 1;
} finally {
  try { ws?.close(); } catch {}
  chrome.kill();
}
