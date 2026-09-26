// Builds ../index.html from src/template.html, embedding the Python engines and example data.
const fs=require('fs'),p=require('path');const root=p.join(__dirname,'..'),eng=p.join(root,'engine');
// JSON string safe inside <script>: '<' becomes <, so '</script' in an embedded file cannot end the element
const esc=s=>JSON.stringify(s).replace(/</g,'\\u003c');
// LF line endings, so the build is byte-identical whatever the checkout's core.autocrlf
const read=f=>fs.readFileSync(p.join(eng,f),'utf8').replace(/\r\n/g,'\n');
// MIT licence notice carried by every embedded engine file (the page can be saved and shared on its own)
const NOTICE='# Bioequivalence Lab - https://bioequivalence-lab.vercel.app\n'+
  '# Copyright (c) 2025 K-Dense Inc. (original pkpd-modeling scripts, https://github.com/K-Dense-AI/scientific-agent-skills)\n'+
  '# Copyright (c) 2026 Ayman Mohammed ALQasem (Aymn74) (corrected engines)\n'+
  '# MIT License: permission is granted free of charge to use, copy, modify, merge, publish, distribute, sublicense and/or\n'+
  '# sell copies, provided this notice is included. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND. Full text: LICENSE\n';
const py=f=>()=>esc(NOTICE+read(f));
const data=f=>()=>esc(read(f));
let h=fs.readFileSync(p.join(__dirname,'template.html'),'utf8').replace(/\r\n/g,'\n');
h=h.replace('__BIOEQUIVALENCE_PY__',py('bioequivalence.py'))
 .replace('__NCA_PY__',py('nca.py'))
 .replace('__COMMON_PY__',py('_common.py'))
 .replace('__EXAMPLE_CSV__',data('example_partial.csv'))
 .replace('__EXAMPLE_NCA_CSV__',data('example_nca.csv'));
fs.writeFileSync(p.join(root,'index.html'),h);console.log('built index.html',h.length);
