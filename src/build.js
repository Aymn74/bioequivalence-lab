// Builds ../index.html from src/template.html, embedding the Python engine and example data.
const fs=require('fs'),p=require('path');const root=p.join(__dirname,'..'),eng=p.join(root,'engine');
const esc=s=>JSON.stringify(s).replace(/</g,'\u003c');
let h=fs.readFileSync(p.join(__dirname,'template.html'),'utf8');
h=h.replace('__BIOEQUIVALENCE_PY__',()=>esc(fs.readFileSync(p.join(eng,'bioequivalence.py'),'utf8')))
 .replace('__COMMON_PY__',()=>esc(fs.readFileSync(p.join(eng,'_common.py'),'utf8')))
 .replace('__EXAMPLE_CSV__',()=>esc(fs.readFileSync(p.join(eng,'example_partial.csv'),'utf8')));
fs.writeFileSync(p.join(root,'index.html'),h);console.log('built index.html',h.length);
