// Builds ../index.html from src/template.html, embedding the Python engines and example data.
const fs=require('fs'),p=require('path');const root=p.join(__dirname,'..'),eng=p.join(root,'engine');
const esc=s=>JSON.stringify(s).replace(/</g,'<');
const embed=f=>()=>esc(fs.readFileSync(p.join(eng,f),'utf8'));
let h=fs.readFileSync(p.join(__dirname,'template.html'),'utf8');
h=h.replace('__BIOEQUIVALENCE_PY__',embed('bioequivalence.py'))
 .replace('__NCA_PY__',embed('nca.py'))
 .replace('__COMMON_PY__',embed('_common.py'))
 .replace('__EXAMPLE_CSV__',embed('example_partial.csv'))
 .replace('__EXAMPLE_NCA_CSV__',embed('example_nca.csv'));
fs.writeFileSync(p.join(root,'index.html'),h);console.log('built index.html',h.length);
