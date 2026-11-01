const fs=require('node:fs'),vm=require('node:vm');
const html=fs.readFileSync('index.html','utf8');
let count=0;
for(const script of html.matchAll(/<script(?:\s[^>]*)?>([\s\S]*?)<\/script>/g)){new vm.Script(script[1]);count++;}
for(const dir of ['assets/js'])for(const file of fs.readdirSync(dir).filter(f=>f.endsWith('.js')))new vm.Script(fs.readFileSync(dir+'/'+file,'utf8'));
new vm.Script(fs.readFileSync('service-worker.js','utf8'));
for(const file of ['package.json','manifest.json'])JSON.parse(fs.readFileSync(file,'utf8'));
for(const script of html.matchAll(/<script src="([^"?]+)/g))if(!fs.existsSync(script[1]))throw Error('Missing script '+script[1]);
if(/https?:\/\/[^"\s]+\.js/.test(html))throw Error('Unexpected runtime remote script');
console.log(count+' inline scripts and local modules parsed; local assets present.');
