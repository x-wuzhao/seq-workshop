const assert=require('node:assert/strict');
const fs=require('node:fs');
const {groups,routes,resolve}=require('../assets/js/module-config');
const html=fs.readFileSync('index.html','utf8');
const ids=[...html.matchAll(/\bid="([^"]+)"/g)].map(m=>m[1]);
assert.equal(new Set(ids).size,ids.length,'all DOM ids must remain unique');
assert.equal(groups.length,6);
const owned=new Set();
for(const r of routes){
  assert(groups.some(g=>g.id===r.group),'route has a parent module');
  assert.equal(resolve('#/'+r.id),r);
  for(const id of r.panels){assert(ids.includes(id),'route panel exists: '+id);assert(!owned.has(id),'panel belongs to one route');owned.add(id);assert.equal(resolve('#'+id),r);}
}
for(const g of groups)assert.equal(resolve(g.defaultRoute).group,g.id);
assert.equal(resolve('#unknown').id,'workspace');
assert.equal(resolve('#nav-primer').id,'design/primers');
assert.equal(resolve('#learn').id,'learn/glossary');
assert(routes.filter(r=>r.group==='learn').every(r=>!r.requiresSequence),'learning works without analyzed DNA');
assert.doesNotMatch(html,/<script>\s*[^<]/,'application logic has editable source modules');
console.log('Module routes, legacy links, panel ownership and independent learning passed.');
