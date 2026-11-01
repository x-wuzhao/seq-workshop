function calcPrimerTm(seq){
  // 简化版 Wallace 公式 + 最近邻修正
  var c = baseCounts(seq);
  var len = seq.length;
  if(len<=13) return (c.A+c.T)*2 + (c.G+c.C)*4;
  // 超过13bp用更精确的公式
  return 64.9 + 41*(c.G+c.C-16.4)/len;
}

function checkDimer(p1, p2){
  // 简化二聚体检测：检查互补碱基配对
  var maxMatch = 0;
  for(var offset=-p1.length+1; offset<p2.length; offset++){
    var match = 0;
    for(var i=0;i<p1.length;i++){
      var j = i + offset;
      if(j>=0 && j<p2.length){
        var a=p1[i], b=p2[p2.length-1-j]; // p2反向
        if((a==='A'&&b==='T')||(a==='T'&&b==='A')||(a==='G'&&b==='C')||(a==='C'&&b==='G')) match++;
      }
    }
    if(match>maxMatch) maxMatch=match;
  }
  return maxMatch;
}

function checkHairpin(seq){
  // 简化发夹结构检测
  var maxStem = 0;
  for(var stemLen=3; stemLen<=8; stemLen++){
    for(var loopLen=3; loopLen<=10; loopLen++){
      for(var start=0; start+stemLen*2+loopLen<=seq.length; start++){
        var s1=seq.substring(start,start+stemLen);
        var s2=seq.substring(start+stemLen+loopLen,start+stemLen*2+loopLen);
        var match=0;
        for(var i=0;i<stemLen;i++){
          var a=s1[i],b=s2[stemLen-1-i];
          if((a==='A'&&b==='T')||(a==='T'&&b==='A')||(a==='G'&&b==='C')||(a==='C'&&b==='G')) match++;
        }
        if(match===stemLen && stemLen>maxStem) maxStem=stemLen;
      }
    }
  }
  return maxStem;
}

function designPrimers(){
  var seq = window._currentSeq || '';
  if(!seq){document.getElementById('primerResult').innerHTML='<div class="msg info">请先在上方输入序列并点击"开始分析"</div>';return;}
  var len = parseInt(document.getElementById('primerLen').value);
  var targetTm = parseFloat(document.getElementById('primerTm').value);
  var productLen = parseInt(document.getElementById('primerProduct').value);
  var offset = parseInt(document.getElementById('primerOffset').value);

  if(seq.length < productLen + len*2 + offset){
    document.getElementById('primerResult').innerHTML='<div class="msg err">序列长度不足，无法设计产物长度为 '+productLen+'bp 的引物</div>';
    return;
  }

  // 正向引物：从offset开始，找Tm最接近targetTm的len-bp序列
  var fStart = offset;
  var fSeq = seq.substring(fStart, fStart+len);
  var fTm = calcPrimerTm(fSeq);
  // 微调起始位置使Tm更接近目标
  var bestFDiff = Math.abs(fTm-targetTm), bestFStart = fStart;
  for(var s=Math.max(0,offset-10); s<=offset+10 && s+len<=seq.length; s++){
    var t = calcPrimerTm(seq.substring(s,s+len));
    if(Math.abs(t-targetTm)<bestFDiff){bestFDiff=Math.abs(t-targetTm);bestFStart=s;}
  }
  fSeq = seq.substring(bestFStart, bestFStart+len);
  fTm = calcPrimerTm(fSeq);

  // 反向引物：在fStart+productLen附近，取反向互补
  var rEnd = bestFStart + productLen;
  var rStart = rEnd - len;
  var rSeqRaw = seq.substring(rStart, rEnd);
  var rSeq = reverseComplement(rSeqRaw);
  var rTm = calcPrimerTm(rSeq);
  // 微调
  var bestRDiff = Math.abs(rTm-targetTm), bestREnd = rEnd;
  for(var e=Math.max(len,rEnd-10); e<=rEnd+10 && e<=seq.length; e++){
    var rs = reverseComplement(seq.substring(e-len,e));
    var t = calcPrimerTm(rs);
    if(Math.abs(t-targetTm)<bestRDiff){bestRDiff=Math.abs(t-targetTm);bestREnd=e;}
  }
  rSeqRaw = seq.substring(bestREnd-len, bestREnd);
  rSeq = reverseComplement(rSeqRaw);
  rTm = calcPrimerTm(rSeq);

  var actualProduct = bestREnd - bestFStart;
  var fc = baseCounts(fSeq), rc = baseCounts(rSeq);
  var fGC = ((fc.G+fc.C)/len*100).toFixed(1);
  var rGC = ((rc.G+rc.C)/len*100).toFixed(1);
  var dimerScore = checkDimer(fSeq, rSeq);
  var fHairpin = checkHairpin(fSeq), rHairpin = checkHairpin(rSeq);

  var warnings = [];
  if(Math.abs(fTm-targetTm)>3) warnings.push('正向引物Tm（'+fTm.toFixed(1)+'°C）偏离目标超过3°C');
  if(Math.abs(rTm-targetTm)>3) warnings.push('反向引物Tm（'+rTm.toFixed(1)+'°C）偏离目标超过3°C');
  if(Math.abs(fTm-rTm)>5) warnings.push('正反向引物Tm差异超过5°C（'+fTm.toFixed(1)+' vs '+rTm.toFixed(1)+'°C）');
  if(dimerScore>=5) warnings.push('引物二聚体风险较高（互补配对'+dimerScore+'bp）');
  if(fHairpin>=4) warnings.push('正向引物可能形成发夹结构（茎长'+fHairpin+'bp）');
  if(rHairpin>=4) warnings.push('反向引物可能形成发夹结构（茎长'+rHairpin+'bp）');
  if(fGC<40||fGC>60) warnings.push('正向引物GC含量（'+fGC+'%）不在40-60%理想范围');
  if(rGC<40||rGC>60) warnings.push('反向引物GC含量（'+rGC+'%）不在40-60%理想范围');

  var html = '<div class="primer-card"><div class="primer-name">正向引物 (Forward)</div>';
  html += '<div class="primer-seq">5\'-'+fSeq+'-3\'</div>';
  html += '<div class="primer-stats"><span>长度：<strong>'+len+' bp</strong></span><span>Tm：<strong>'+fTm.toFixed(1)+'°C</strong></span><span>GC：<strong>'+fGC+'%</strong></span><span>起始位置：<strong>'+(bestFStart+1)+'</strong></span></div></div>';
  html += '<div class="primer-card"><div class="primer-name">反向引物 (Reverse)</div>';
  html += '<div class="primer-seq">5\'-'+rSeq+'-3\'</div>';
  html += '<div class="primer-stats"><span>长度：<strong>'+len+' bp</strong></span><span>Tm：<strong>'+rTm.toFixed(1)+'°C</strong></span><span>GC：<strong>'+rGC+'%</strong></span><span>结合位置：<strong>'+(bestREnd-len+1)+'-'+bestREnd+'</strong></span></div></div>';
  html += '<div class="primer-card"><div class="primer-name">产物信息</div>';
  html += '<div class="primer-stats"><span>产物长度：<strong>'+actualProduct+' bp</strong></span><span>Tm差异：<strong>'+Math.abs(fTm-rTm).toFixed(1)+'°C</strong></span><span>二聚体评分：<strong>'+dimerScore+' bp</strong></span></div></div>';
  if(warnings.length>0){
    html += '<div class="primer-card" style="border-color:#F0C9C9;background:#FDF5F5"><div class="primer-name primer-warn">⚠ 注意事项</div><ul style="margin:8px 0 0 18px;font-size:13px;color:var(--ink-soft)">';
    for(var i=0;i<warnings.length;i++) html += '<li>'+warnings[i]+'</li>';
    html += '</ul></div>';
  } else {
    html += '<div class="primer-card" style="border-color:#C4DEDC;background:#F0F8F7"><div class="primer-name" style="color:#2A7A2A">✓ 候选引物通过基础筛选</div><p style="font-size:13px;color:var(--ink-soft);margin-top:6px">各项简化指标处于目标范围内；实际 PCR 前仍需用专业工具检查特异性并进行实验验证。</p></div>';
  }
  document.getElementById('primerResult').innerHTML = html;
}

function initPrimerDesign(){
  document.getElementById('primerDesign').onclick = designPrimers;
}

/* ============================================================
   六、密码子优化器
   ============================================================ */
