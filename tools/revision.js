const fs=require('node:fs'),crypto=require('node:crypto');
const file='service-worker.js',source=fs.readFileSync(file,'utf8');
const paths=[...new Set(['index.html',...Array.from(source.matchAll(/new URL\('\.\/([^']+)'/g),m=>m[1])])].sort();
const hash=crypto.createHash('sha256');for(const p of paths)hash.update(p+'\0').update(fs.readFileSync(p));
const revision='v4-'+hash.digest('hex').slice(0,16);
const updated=source.replace(/const CACHE_NAME = CACHE_PREFIX \+ '[^']+';/,"const CACHE_NAME = CACHE_PREFIX + '"+revision+"';");
if(process.argv.includes('--check')){if(updated!==source)throw Error('App shell changed: run node tools/revision.js');}
else fs.writeFileSync(file,updated);
console.log('App shell revision:',revision,'('+paths.length+' resources)');
