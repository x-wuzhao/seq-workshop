const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict');
const code=fs.readFileSync('service-worker.js','utf8');
const events={},store=new Map(),deleted=[];
const unrelated='seq-workshop-shell-v3';store.set(unrelated,new Map());
function cache(name){if(!store.has(name))store.set(name,new Map());const entries=store.get(name);return {async addAll(urls){for(const u of urls)entries.set(u,new Response('cached '+u));},async match(request){const value=entries.get(typeof request==='string'?request:request.url);return value?.clone();},async put(key,value){entries.set(typeof key==='string'?key:key.url,value);}};}
const ctx={URL,Response,console,fetch:async()=>{throw Error('offline');},caches:{async open(name){return cache(name);},async keys(){return [...store.keys()];},async delete(name){deleted.push(name);return store.delete(name);}},self:{registration:{scope:'https://example.test/next/'},location:{origin:'https://example.test'},clients:{async claim(){}},async skipWaiting(){},addEventListener(name,fn){events[name]=fn;}}};
vm.createContext(ctx);vm.runInContext(code,ctx);
async function dispatch(name,event={}){let work;events[name]({...event,waitUntil(p){work=p;},respondWith(p){work=p;}});return await work;}
(async()=>{
  await dispatch('install');await dispatch('activate');
  assert.equal(store.has(unrelated),true,'another app edition cache must survive activation');assert.deepEqual(deleted,[]);
  const response=await dispatch('fetch',{request:{method:'GET',mode:'navigate',url:'https://example.test/next/'}});
  assert.equal(response.status,200);assert.match(await response.text(),/cached https:\/\/example.test\/next/);
  const core=await dispatch('fetch',{request:{method:'GET',mode:'cors',url:'https://example.test/next/assets/js/sequence-core.js'}});
  assert.equal(core.status,200);
  let intercepted=false;events.fetch({request:{method:'POST',url:'https://example.test/next/'},respondWith(){intercepted=true;}});assert.equal(intercepted,false);
  events.fetch({request:{method:'GET',url:'https://another.test/'},respondWith(){intercepted=true;}});assert.equal(intercepted,false);
  console.log('Service Worker offline navigation, local core cache and edition isolation passed.');
})().catch(error=>{console.error(error);process.exitCode=1;});
