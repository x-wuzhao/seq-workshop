/* One registry drives routes, navigation labels and panel ownership. */
(function(root,factory){
  var value=factory();
  if(typeof module==='object'&&module.exports)module.exports=value;
  else root.SeqModules=value;
})(typeof globalThis!=='undefined'?globalThis:this,function(){
  'use strict';
  var groups=[
    {id:'workspace',label:'工作区',defaultRoute:'workspace',description:'载入 DNA 序列，在地图中浏览和标注，再进入相应工具。'},
    {id:'sequence',label:'序列分析',defaultRoute:'sequence/transcription',description:'从互补链到蛋白质，分步查看序列的转录、翻译和组成。'},
    {id:'design',label:'实验设计',defaultRoute:'design/digest',description:'在同一份模板上检查限制性酶切，筛选 PCR 引物。'},
    {id:'compare',label:'序列比对',defaultRoute:'compare',description:'以当前序列为参考，比较另一条序列的相同与差异。'},
    {id:'explore',label:'结构与突变',defaultRoute:'explore/mutation',description:'观察突变如何改变翻译结果，结合三维模型理解结构。'},
    {id:'learn',label:'学习中心',defaultRoute:'learn/dogma',description:'用分步演示、练习与概念速查，理解序列背后的生物学。'}
  ];
  var routes=[
    {id:'workspace',group:'workspace',label:'工作区',panels:['nav-map'],requiresSequence:false},
    {id:'sequence/transcription',group:'sequence',label:'互补与转录',panels:['nav-complement']},
    {id:'sequence/translation',group:'sequence',label:'六框翻译',panels:['nav-frames']},
    {id:'sequence/protein',group:'sequence',label:'蛋白性质',panels:['nav-protein']},
    {id:'sequence/codons',group:'sequence',label:'密码子频率',panels:['nav-codon']},
    {id:'sequence/optimize',group:'sequence',label:'密码子优化',panels:['nav-optimize']},
    {id:'sequence/charts',group:'sequence',label:'组成图表',panels:['nav-charts','nav-aa']},
    {id:'design/digest',group:'design',label:'限制性酶切',panels:['nav-enzyme']},
    {id:'design/primers',group:'design',label:'PCR 引物',panels:['nav-primer']},
    {id:'compare',group:'compare',label:'序列比对',panels:['nav-align']},
    {id:'explore/mutation',group:'explore',label:'突变模拟',panels:['nav-mutation']},
    {id:'explore/structure',group:'explore',label:'三维结构',panels:['nav-structure']},
    {id:'learn/dogma',group:'learn',label:'中心法则',panels:['nav-dogma'],requiresSequence:false},
    {id:'learn/quiz',group:'learn',label:'知识练习',panels:['panel-quiz'],requiresSequence:false},
    {id:'learn/glossary',group:'learn',label:'概念速查',panels:['learn'],requiresSequence:false}
  ];
  routes.forEach(function(r){if(r.requiresSequence===undefined)r.requiresSequence=true;});
  function resolve(value){
    var id=String(value||'').replace(/^#\/?/,'').replace(/\/$/,'');
    var exact=routes.find(function(r){return r.id===id;});
    if(exact)return exact;
    var legacy=routes.find(function(r){return r.panels.indexOf(id)>=0;});
    if(legacy)return legacy;
    var group=groups.find(function(g){return g.id===id;});
    return routes.find(function(r){return r.id===(group?group.defaultRoute:'workspace');});
  }
  return {groups:groups,routes:routes,resolve:resolve};
});
