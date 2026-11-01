/* Task views share the existing workspace; switching never recreates tool DOM. */
(function(){
  'use strict';
  var registry=SeqModules,current=registry.resolve(location.hash),pendingRoute=null;
  var $=function(id){return document.getElementById(id);};
  var panelIds=registry.routes.reduce(function(ids,route){return ids.concat(route.panels);},[]);
  function hasSequence(){
    return !!window._currentSeq && !$('results').hidden && SeqCore.parseDNA($('seq-input').value).seq===window._currentSeq;
  }
  function context(){
    var valid=hasSequence(),input=$('seq-input').value,seq=valid?window._currentSeq:'';
    var header=input.match(/^\s*>[^\r\n]+/),name='当前 DNA 序列';
    if(header)name=header[0].replace(/^\s*>/,'').trim();
    else if(valid){
      var sample=SAMPLE_SEQUENCES.find(function(s){return (s.seq||SAMPLE_DNA)===seq;});
      if(sample)name=sample.name;
    }
    $('contextName').textContent=valid?name:(input.trim()?'输入待分析':'尚未载入序列');
    $('contextMeta').textContent=valid?seq.length.toLocaleString()+' bp · GC '+SeqCore.composition(seq).gcPct.toFixed(1)+'% · 所有模块共用':(input.trim()?'完成分析后，工具与结果会使用这份新序列':'先加载一个示例，或输入自己的 DNA 序列');
    $('contextDot').classList.toggle('ready',valid);
    $('editSequence').textContent=valid?'更换序列':'输入序列';
    $('backupWorkspace').disabled=!valid;
    $('downloadAnalysisReport').disabled=!valid;
    $('sequenceEditor').querySelector('summary').firstChild.textContent=valid?'输入与更换序列 ':'输入 DNA 序列 ';
  }
  function sync(options){
    options=options||{};
    var valid=hasSequence(),group=registry.groups.find(function(g){return g.id===current.group;});
    panelIds.forEach(function(id){$(id).hidden=!(current.panels.indexOf(id)>=0 && (!current.requiresSequence||valid) && (id!=='nav-map'||valid));});
    $('inputSection').hidden=current.group!=='workspace';
    $('sequenceSummary').hidden=current.group!=='workspace'||!valid;
    $('moduleEmpty').hidden=!current.requiresSequence||valid;
    $('moduleTitle').textContent=group.label;
    $('moduleEyebrow').textContent='序列工坊 / '+(current.group==='workspace'?'本地 DNA 工作台':group.label);
    $('moduleDescription').textContent=group.description;
    document.title=(current.group==='workspace'?group.label:current.label)+' · 序列工坊 Next';
    document.body.dataset.module=current.group;
    document.querySelectorAll('#moduleNav a, #siteLinks a').forEach(function(a){
      var selected=registry.resolve(a.hash).group===current.group;
      a.classList.toggle('active',selected);
      if(selected)a.setAttribute('aria-current','page');else a.removeAttribute('aria-current');
    });
    var siblings=registry.routes.filter(function(r){return r.group===current.group;});
    $('moduleTabs').hidden=siblings.length<2;
    // Links are native navigation, retaining standard keyboard and history behavior.
    if($('moduleTabs').dataset.group!==current.group){
      $('moduleTabs').replaceChildren();
      siblings.forEach(function(r){var a=document.createElement('a');a.href='#/'+r.id;a.textContent=r.label;$('moduleTabs').appendChild(a);});
      $('moduleTabs').dataset.group=current.group;
    }
    $('moduleTabs').querySelectorAll('a').forEach(function(a){
      var selected=registry.resolve(a.hash).id===current.id;
      a.classList.toggle('active',selected);
      if(selected)a.setAttribute('aria-current','page');else a.removeAttribute('aria-current');
    });
    if(current.id!=='learn/dogma')stopDogmaAuto();
    context();
    if(current.id==='learn/dogma'){
      $('dogmaSequenceNote').textContent=valid?'演示使用当前工作区的序列。':'尚未分析序列；演示使用内置教学示例。';
      renderDogmaStep();
    }
    requestAnimationFrame(function(){
      if(current.id==='workspace'&&hasSequence())renderSeqMap();
      if(current.id==='sequence/charts'&&hasSequence()&&lastChartArgs)doChartsRender(lastChartArgs.seq,lastChartArgs.avgGC,lastChartArgs.aaCounts);
      if(current.id==='explore/structure'&&hasSequence()){
        renderStructure3D();
        window.dispatchEvent(new Event('resize'));
      }
    });
    if(options.focus){$('moduleTitle').focus({preventScroll:true});window.scrollTo({top:0,behavior:'instant'});}
  }
  function navigate(value,options){
    var target=registry.resolve(value),hash='#/'+target.id;
    if(options&&options.replace){history.replaceState(null,'',hash);current=target;sync(options);}
    else if(location.hash!==hash)location.hash=hash;
    else{current=target;sync({focus:true});}
  }
  function openInput(){
    if(current.id!=='workspace')pendingRoute=current.id;
    navigate('workspace');
    $('sequenceEditor').open=true;
    requestAnimationFrame(function(){$('seq-input').focus({preventScroll:true});$('inputSection').scrollIntoView({block:'start',behavior:'smooth'});});
  }
  function afterAnalyze(){
    if(!hasSequence()){sync();return;}
    $('sequenceEditor').open=false;
    if(pendingRoute){var route=pendingRoute;pendingRoute=null;navigate(route);}
    else sync();
  }
  window.SeqModuleRouter={navigate:navigate,sync:sync,afterAnalyze:afterAnalyze,hasSequence:hasSequence,get current(){return current.id;}};
  var note=document.createElement('p');note.id='dogmaSequenceNote';note.className='hint dogma-source';
  $('nav-dogma').querySelector('h2').after(note);
  $('editSequence').addEventListener('click',openInput);
  $('emptyInput').addEventListener('click',openInput);
  $('emptySample').addEventListener('click',function(){
    $('seq-input').value=SAMPLE_DNA;
    $('seq-input').dispatchEvent(new Event('input',{bubbles:true}));
    analyze();
  });
  $('emptyRestore').addEventListener('click',function(){$('workspaceFile').click();});
  $('seq-input').addEventListener('input',function(){sync();});
  $('btnClear').addEventListener('click',function(){pendingRoute=null;$('sequenceEditor').open=true;sync();});
  $('workspaceFile').addEventListener('change',function(){setTimeout(sync,0);});
  window.addEventListener('hashchange',function(){navigate(location.hash,{replace:true,focus:true});});
  // Normalize both legacy anchor links and unknown hashes to a stable local route.
  document.addEventListener('click',function(e){
    if(e.target.closest&&e.target.closest('.skip-link')){e.preventDefault();$('moduleTitle').focus();return;}
    var anchor=e.target.closest&&e.target.closest('a[href^="#nav-"],a[href="#panel-quiz"],a[href="#learn"]');
    if(anchor){e.preventDefault();navigate(anchor.getAttribute('href'));}
  });
  navigate(current.id,{replace:true});
})();
