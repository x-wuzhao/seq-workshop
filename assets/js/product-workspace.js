/* Portable backups are data only. Computed results are always invalidated/recomputed. */
(function(){
  'use strict';
  function status(text){document.getElementById('workspaceStatus').textContent=text;}
  function download(name,value){
    var blob=new Blob([JSON.stringify(value,null,2)],{type:'application/json;charset=utf-8'});
    var url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=name;a.click();
    setTimeout(function(){URL.revokeObjectURL(url);},1000);
  }
  function currentWorkspace(){
    if(!window._currentSeq)throw Error('请先完成一次有效的序列分析。');
    if(SeqCore.parseDNA(document.getElementById('seq-input').value).seq!==window._currentSeq)throw Error('输入已改变，请重新分析后再备份。');
    return SeqWorkspaceFormat.create(window._currentSeq,mapFeatures,{codonTable:activeCodonTable,topology:seqMapState.circular?'circular':'linear'});
  }
  function importWorkspace(value){
    // Validate everything before touching the current workspace.
    var clean=SeqWorkspaceFormat.validate(value);
    document.getElementById('seq-input').value=clean.sequence;
    activeCodonTable=clean.settings.codonTable;
    document.getElementById('codonTableSel').value=activeCodonTable;
    window._pendingFeatures=clean.features;
    // Force pending annotations to replace those of an identical sequence, too.
    window._lastAnalyzedSeq=null;
    analyze();
    seqMapState.circular=clean.settings.topology==='circular';
    document.getElementById('digestTopo').value=clean.settings.topology;
    document.getElementById('mapLinear').classList.toggle('active',!seqMapState.circular);
    document.getElementById('mapCircular').classList.toggle('active',seqMapState.circular);
    renderSeqMap();renderMapFeatures();
    if(lastSeq&&lastEnz)renderDigest(lastSeq,lastEnz);
    saveWorkspaceNow();
    status('已恢复 '+clean.sequence.length+' bp、'+clean.features.length+' 条注释；分析结果已按当前模型重新计算。');
  }
  document.getElementById('seq-input').addEventListener('input',function(){
    if(window._currentSeq&&SeqCore.parseDNA(this.value).seq!==window._currentSeq){
      document.getElementById('results').hidden=true;
      status('输入已改变，请重新分析。上一次分析结果已隐藏，避免混用。');
    }
  });
  document.getElementById('backupWorkspace').addEventListener('click',function(){
    try{download('sequence-workshop-next-workspace.json',currentWorkspace());status('工作区备份已导出。文件包含原始序列与注释，请保存在可信位置。');}
    catch(error){status(error.message);}
  });
  function analysisRecord(){
    var workspace=currentWorkspace(),seq=workspace.sequence,reverse=SeqCore.reverseComplement(seq),frames=[];
    ['+','-'].forEach(function(sign){for(var f=0;f<3;f++){
      var candidates=findORFs(sign==='+'?seq:reverse,f),best=null;
      candidates.forEach(function(orf){if(!best||orf.len>best.len)best=orf;});
      frames.push({frame:sign+(f+1),completeORFCount:candidates.length,longest:best?{start:sign==='+'?best.start+1:seq.length-best.end+1,end:sign==='+'?best.end:seq.length-best.start,coordinates:'1-based inclusive, stop codon included',aminoAcids:best.prot,stopCodon:best.stopCodon}:null});
    }});
    return {format:'seq-workshop-analysis-record',schemaVersion:1,productVersion:SeqCore.version,createdAt:new Date().toISOString(),workspace:workspace,composition:SeqCore.composition(seq),readingFrames:frames,primerCandidates:window._lastPrimers||[],models:{geneticCode:activeCodonTable==='standard'?1:2,ORF:'complete start-to-stop; alternative CDS start translated as M; nested ORFs retained',meltingTemperature:{parameters:'SantaLucia 1998 / DNA_NN3',Na_mM:50,totalOligo_nM:250,saltCorrection:'entropy',includesMg:false},codonPreference:'rank-derived teaching score; not standard CAI',structures:'idealized B-DNA or bundled experimental PDB structures; not structure prediction'},sources:['https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi','https://github.com/biopython/biopython/blob/biopython-186/Bio/SeqUtils/MeltingTemp.py'],limitations:['Learning and preliminary analysis only; not an experimental protocol.','Reference checks cover specified fixtures, not all possible inputs.']};
  }
  document.getElementById('downloadAnalysisReport').addEventListener('click',function(){
    try{download('sequence-workshop-next-analysis.json',analysisRecord());status('计算记录已导出：包含序列、参数、模型版本和结果，可用于复核。');}catch(error){status(error.message);}
  });
  function networkStatus(){
    var label=navigator.onLine?'浏览器报告在线':'浏览器报告离线';
    document.getElementById('connectivityStatus').textContent=label+(navigator.serviceWorker&&navigator.serviceWorker.controller?' · 离线应用壳已接管当前页面':' · 离线缓存需首次联网访问完成');
  }
  window.addEventListener('online',networkStatus);window.addEventListener('offline',networkStatus);
  if(navigator.serviceWorker)navigator.serviceWorker.addEventListener('controllerchange',networkStatus);
  networkStatus();
  var file=document.getElementById('workspaceFile');
  document.getElementById('importWorkspace').addEventListener('click',function(){file.click();});
  file.addEventListener('change',async function(){
    var selected=this.files&&this.files[0];this.value='';if(!selected)return;
    if(selected.size>2*1024*1024){status('备份超过 2 MB，未导入；当前工作区保持原样。');return;}
    try{importWorkspace(JSON.parse(await selected.text()));}catch(error){status('未导入：'+error.message);}
  });
  window.SeqProduct={currentWorkspace:currentWorkspace,importWorkspace:importWorkspace,analysisRecord:analysisRecord};
})();
