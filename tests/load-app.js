const fs=require('node:fs'),vm=require('node:vm');
module.exports=function loadApp(){
  function element(){return {value:'',innerHTML:'',textContent:'',hidden:false,disabled:false,style:{},children:[],parentNode:{replaceChild(){}},className:'',classList:{add(){},remove(){},toggle(){},contains(){return false;}},setAttribute(){},getAttribute(){return null;},addEventListener(){},appendChild(e){this.children.push(e);},removeChild(){},querySelector(){return element();},querySelectorAll(){return [];},cloneNode(){return element();},scrollIntoView(){},closest(){return null;},focus(){},click(){},getBoundingClientRect(){return {top:0,bottom:1};}};}
  const elements={};
  const document={documentElement:element(),body:element(),getElementById(id){return elements[id]||(elements[id]=element());},querySelector(){return element();},querySelectorAll(){return [];},createElement(){return element();},addEventListener(){}};
  const storage={getItem(){return null;},setItem(){},removeItem(){}};
  const window={localStorage:storage,addEventListener(){},scrollY:0,scrollTo(){}};
  const ctx={SeqCore:require('../assets/js/sequence-core'),SeqWorkspaceFormat:require('../assets/js/workspace-format'),document,window,console,navigator:{serviceWorker:null},location:{protocol:'http:',hostname:'localhost'},localStorage:storage,setTimeout(){return 0;},clearTimeout(){},setInterval(){return 0;},clearInterval(){},requestAnimationFrame(){return 0;},cancelAnimationFrame(){},Blob:function(){},URL:{createObjectURL(){return '';},revokeObjectURL(){}},alert(){},matchMedia(){return {matches:false};}};
  ctx.globalThis=ctx;vm.createContext(ctx);
  const html=fs.readFileSync('index.html','utf8');
  for(const match of html.matchAll(/<script(?:\s[^>]*)?>([\s\S]*?)<\/script>/g))vm.runInContext(match[1],ctx);
  return ctx;
};
