(function(root,factory){var api=factory();if(typeof module==='object'&&module.exports)module.exports=api;else root.SeqWorkspaceFormat=api;})(typeof globalThis!=='undefined'?globalThis:this,function(){
  'use strict';
  var FORMAT='seq-workshop-workspace',VERSION=1;
  function validate(value){
    if(!value||value.format!==FORMAT||value.schemaVersion!==VERSION)throw Error('不是支持的序列工坊工作区备份（schemaVersion=1）。');
    var seq=value.sequence;
    if(typeof seq!=='string'||seq.length<3||seq.length>50000||!/^[ATCGN]+$/.test(seq))throw Error('备份序列无效：须为 3–50,000 bp 的 A/T/C/G/N。');
    var features=value.features||[];
    if(!Array.isArray(features)||features.length>1000)throw Error('注释数量超出限制（最多 1,000 条）。');
    features=features.map(function(f){
      if(!f||!Number.isInteger(f.start)||!Number.isInteger(f.end)||f.start<1||f.end<f.start||f.end>seq.length||!['+','-'].includes(f.strand))throw Error('备份含越界或方向无效的注释。');
      if(typeof f.label!=='string'||f.label.length>200||(f.note&&typeof f.note!=='string')||(f.note||'').length>2000)throw Error('注释名称或备注超出限制。');
      return {start:f.start,end:f.end,strand:f.strand,label:f.label,note:f.note||'',key:'misc_feature'};
    });
    var settings=value.settings||{};
    if(settings.codonTable&&!['standard','vertebrate_mito'].includes(settings.codonTable))throw Error('备份的遗传密码表不受支持。');
    if(settings.topology&&!['linear','circular'].includes(settings.topology))throw Error('备份的序列拓扑不受支持。');
    return {format:FORMAT,schemaVersion:VERSION,productVersion:'4.0.0',sequence:seq,features:features,settings:{codonTable:settings.codonTable||'standard',topology:settings.topology||'linear'},savedAt:typeof value.savedAt==='string'?value.savedAt:''};
  }
  function create(sequence,features,settings){return validate({format:FORMAT,schemaVersion:VERSION,sequence:sequence,features:features,settings:settings,savedAt:new Date().toISOString()});}
  return {validate:validate,create:create};
});
