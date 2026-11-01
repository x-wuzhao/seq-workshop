function optimizeCodons(seq, species){
  var bias = SPECIES_CODON_BIAS[species];
  var map=codonMap();
  var mrna = transcribe(seq);
  var optimized = '';
  var changes = 0;
  for(var i=0;i+3<=mrna.length;i+=3){
    var codon = mrna.substring(i,i+3);
    var dnaCodon = codon.replace(/U/g,'T');
    var aa = map[dnaCodon];
    if(!aa || !bias[aa]){optimized+=dnaCodon;continue;}
    var preferred = aa==='*'?dnaCodon:bias[aa][0]; // 取该物种最高频密码子
    if(preferred !== dnaCodon) changes++;
    optimized += preferred;
  }
  // 处理末尾不足3个碱基的部分
  if(mrna.length % 3 !== 0) optimized += mrna.substring(Math.floor(mrna.length/3)*3);
  // 转回DNA（U→T）
  var optDNA = optimized.replace(/U/g,'T');
  return {seq:optDNA, changes:changes, totalCodons:Math.floor(mrna.length/3)};
}

function calcCodonFreq(seq){
  var freq = {};
  var total = 0;
  for(var i=0;i+3<=seq.length;i+=3){
    var c = seq.substring(i,i+3);
    freq[c] = (freq[c]||0)+1;
    total++;
  }
  return {freq:freq, total:total};
}

function renderOptimizeResult(){
  var seq = window._currentSeq || '';
  if(!seq){document.getElementById('optResult').innerHTML='<div class="msg info">请先在上方输入序列并点击"开始分析"</div>';return;}
  var species = document.getElementById('optSpecies').value;
  var speciesName = {ecoli:'大肠杆菌',yeast:'酿酒酵母',human:'人'}[species];
  var result = optimizeCodons(seq, species);

  var bias=SPECIES_CODON_BIAS[species];
  var map=codonMap();
  if(seq.length%3!==0 || CODON_TABLES[activeCodonTable].starts.indexOf(seq.slice(0,3))<0 || translateFrame(seq.slice(0,-3),0).indexOf('*')>=0){
    document.getElementById('optResult').textContent='请提供从起始密码子开始、长度为 3 的倍数、无内部终止密码子的 +1 编码序列。';return;
  }
  var caiBefore=SeqCore.preferenceScore(seq,bias,map),caiAfter=SeqCore.preferenceScore(result.seq,bias,map);
  function scoreText(value){return value===null?'未知':value.toFixed(3);}


  var proteinSame=translateFrame(seq,0)===translateFrame(result.seq,0);

  var html = '<div class="opt-compare">';
  html += '<div class="opt-box"><h4>原始序列（'+seq.length+'bp）</h4><div class="seq">'+seqHTML(seq)+'</div></div>';
  html += '<div class="opt-box"><h4>优化后序列（'+result.seq.length+'bp，'+speciesName+'偏好）</h4><div class="seq">'+seqHTML(result.seq)+'</div></div>';
  html += '</div>';
  html += '<div class="opt-stats">';
  html += '<div class="opt-stat"><div class="k">密码子改变</div><div class="v">'+result.changes+'/'+result.totalCodons+'</div></div>';
  html += '<div class="opt-stat"><div class="k">替换前排名分数</div><div class="v">'+scoreText(caiBefore)+'</div></div>';
  html += '<div class="opt-stat"><div class="k">替换后排名分数</div><div class="v" style="color:var(--teal)">'+scoreText(caiAfter)+'</div></div>';
  html += '<div class="opt-stat"><div class="k">蛋白序列</div><div class="v" style="color:'+(proteinSame?'#2A7A2A':'#9A2A2A')+'">'+(proteinSame?'一致 ✓':'改变 ✗')+'</div></div>';
  html += '</div>';
  if(proteinSame){
    html += '<div class="mut-explain" style="margin-top:14px"><strong>优化说明：</strong>本工具在不改变蛋白质序列的前提下，按'+speciesName+'的简化密码子偏好进行替换。排名分数是各有义密码子 1/排名 的几何平均值（不是标准 CAI），仅用于本表内教学对比。频率表未经过表达系统校准，不预测实验表达量；含未知密码子时不评分。</div>';
  } else {
    html += '<div class="mut-explain" style="margin-top:14px;background:#FBE0E0;color:#9A2A2A"><strong>注意：</strong>优化后蛋白质序列发生了改变！这可能是因为序列中存在终止密码子或序列长度不是3的倍数。请检查输入序列是否为完整的编码区。</div>';
  }
  document.getElementById('optResult').innerHTML = html;
  window._optimizedSeq = result.seq;
}

function initCodonOptimizer(){
  document.getElementById('optApply').onclick = renderOptimizeResult;
  document.getElementById('optCopy').onclick = function(){
    if(!window._optimizedSeq){alert('请先点击"优化密码子"');return;}
    var ta = document.createElement('textarea');
    ta.value = window._optimizedSeq;
    document.body.appendChild(ta);
    ta.select();
    document.execCommand('copy');
    document.body.removeChild(ta);
    this.textContent = '已复制 ✓';
    var btn = this;
    setTimeout(function(){btn.textContent='复制优化后序列';},2000);
  };
}

/* ============================================================
   七、导出功能
   ============================================================ */
