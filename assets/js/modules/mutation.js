function applyMutation(seq, type, pos, base){
  pos = Number(pos) - 1; // 转为0-based，拒绝小数位置
  if(!Number.isInteger(pos) || pos<0 || pos>=seq.length) return {error:'位置超出序列范围（1-' + seq.length + '）'};
  var mutated = seq;
  var detail = '';
  if(type==='point'){
    var old = seq[pos];
    if(old===base) return {error:'该位置已经是 ' + base + '，无需替换'};
    mutated = seq.substring(0,pos) + base + seq.substring(pos+1);
    detail = '第 '+(pos+1)+' 位：' + old + ' → ' + base;
  } else if(type==='insert'){
    mutated = seq.substring(0,pos) + base + seq.substring(pos);
    detail = '在第 '+(pos+1)+' 位前插入 ' + base;
  } else if(type==='delete'){
    var del = seq[pos];
    mutated = seq.substring(0,pos) + seq.substring(pos+1);
    detail = '删除第 '+(pos+1)+' 位的 ' + del;
  }
  return {seq:mutated, detail:detail, type:type, pos:pos, base:base};
}

function classifyMutation(origSeq, mutSeq, type, pos, frame){
  var table=typeof codonMap==='function'?codonMap():CODON;
  var starts=typeof CODON_TABLES!=='undefined'?CODON_TABLES[activeCodonTable].starts:['ATG'];
  return SeqCore.mutationEffect(origSeq,mutSeq,type,pos,frame||'+1',table,starts);
}

function renderMutationResult(){
  var origSeq = window._currentSeq || '';
  if(!origSeq){document.getElementById('mutResult').innerHTML='<div class="msg info">请先在上方输入序列并点击"开始分析"</div>';return;}
  var type = document.getElementById('mutType').value;
  var pos = document.getElementById('mutPos').value;
  var baseEl=document.getElementById('mutBase');
  var base=baseEl?baseEl.value:'A';
  var result = applyMutation(origSeq, type, pos, base);
  if(result.error){document.getElementById('mutResult').innerHTML='<div class="msg err">'+result.error+'</div>';return;}
  var frame=document.getElementById('mutationFrame').value;
  var cls = classifyMutation(origSeq,result.seq,type,result.pos,frame);
  var p1=cls.originalProtein,p2=cls.mutantProtein;

  // 高亮突变位置
  var mutPos = result.pos;
  var origDisplay = '';
  for(var i=0;i<origSeq.length;i++){
    var ch = origSeq[i];
    var cls2 = 'b'+ch;
    if(i===mutPos && type==='point') cls2 += '" style="background:var(--danger-bg);padding:1px 3px;border-radius:3px;font-weight:700';
    if(i===mutPos && type==='delete') cls2 += '" style="background:#FBE0E0;padding:1px 3px;border-radius:3px;text-decoration:line-through';
    origDisplay += '<span class="'+cls2+'">'+ch+'</span>';
  }
  var mutDisplay = '';
  for(var i=0;i<result.seq.length;i++){
    var ch = result.seq[i];
    var cls2 = 'b'+ch;
    if(i===mutPos && type==='point') cls2 += '" style="background:#E0F0E0;padding:1px 3px;border-radius:3px;font-weight:700';
    if(i===mutPos && type==='insert') cls2 += '" style="background:#E0F0E0;padding:1px 3px;border-radius:3px;font-weight:700';
    mutDisplay += '<span class="'+cls2+'">'+ch+'</span>';
  }
  window._seqMapMutation = { pos: result.pos, type: type, label: cls.label, detail: result.detail };
  if (typeof renderSeqMap === 'function') renderSeqMap();
  if (typeof scheduleWorkspaceSave === 'function') scheduleWorkspaceSave();
  var html = '<div class="mut-explain"><strong>突变详情：</strong>'+result.detail+' <span class="mut-badge '+cls.type+'">'+cls.label+'</span></div>';
  html += '<div class="mut-result-grid">';
  html += '<div class="mut-compare"><h4>原始 DNA 序列</h4><div class="seq">'+origDisplay+'</div><h4 style="margin-top:12px">原始蛋白（'+frame+' 原始 ORF，'+p1.length+'aa）</h4><div class="seq mono" style="font-size:11px">'+p1+'</div></div>';
  html += '<div class="mut-compare"><h4>突变后 DNA 序列</h4><div class="seq">'+mutDisplay+'</div><h4 style="margin-top:12px">突变后蛋白（相同起点，'+p2.length+'aa）</h4><div class="seq mono" style="font-size:11px">'+(cls.type==='start-loss'?'起始位点改变，不作确定的产物预测':(p2||'所选读码框无可比较的完整 ORF'))+'</div></div>';
  html += '</div>';
  html += '<div class="mut-explain" style="margin-top:14px"><strong>'+cls.label+'说明：</strong>'+cls.explain+'</div>';
  document.getElementById('mutResult').innerHTML = html;
}

function initMutationSimulator(){
  document.getElementById('mutApply').onclick = renderMutationResult;
  document.getElementById('mutationFrame').onchange=function(){document.getElementById('mutResult').textContent='读码框已改变，请重新应用突变。';};
  document.getElementById('mutReset').onclick = function(){
    document.getElementById('mutResult').innerHTML = '<div class="msg info">已恢复原序列，点击"应用突变"查看效果</div>';
    window._seqMapMutation = null;
    if (typeof renderSeqMap === 'function') renderSeqMap();
  };
  document.getElementById('mutType').onchange = function(){
    var v = this.value;
    var label = document.getElementById('mutBaseLabel');
    if(v==='point') label.innerHTML='替换为：<select id="mutBase"><option value="A">A 腺嘌呤</option><option value="T">T 胸腺嘧啶</option><option value="G">G 鸟嘌呤</option><option value="C">C 胞嘧啶</option></select>';
    else if(v==='insert') label.innerHTML='插入：<select id="mutBase"><option value="A">A 腺嘌呤</option><option value="T">T 胸腺嘧啶</option><option value="G">G 鸟嘌呤</option><option value="C">C 胞嘧啶</option></select>';
    else label.innerHTML='（删除无需选择碱基）';
  };
}

/* ============================================================
   五、引物设计助手
   ============================================================ */
