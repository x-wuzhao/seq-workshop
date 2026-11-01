function resetSequenceDependentPanels(){
  ['mutResult','primerResult','alignResult','optResult'].forEach(function(id){
    var el=document.getElementById(id);
    if(el)el.innerHTML='';
  });
  var alignInput=document.getElementById('alignInput2');
  if(alignInput)alignInput.value='';
  window._optimizedSeq='';
  window._lastPrimers=[];
  window._seqMapMutation=null;
  if(typeof seqMapState!=='undefined')seqMapState.sel=null;
  lastSeq='';
  lastEnz=null;
  lastChartArgs=null;
  if(charts.gc){charts.gc.dispose();charts.gc=null;}
  if(charts.aa){charts.aa.dispose();charts.aa=null;}
  stopDogmaAuto();
  dogmaStep=0;
  if(typeof s3d!=='undefined'&&s3d.viewer){
    s3d.viewer.setAutoRotate(false);
    s3d.viewer.resetView();
    s3d.mode='seq';
    s3d.posTouched=false;
    var spinBtn=document.getElementById('s3dSpin');
    if(spinBtn){spinBtn.setAttribute('aria-pressed','false');spinBtn.textContent='⟳ 自动旋转';}
    var winSlider=document.getElementById('structurePosition');
    if(winSlider){winSlider.max=1;winSlider.value=1;}
  }
}

function analyze(){
  var raw=document.getElementById('seq-input').value;
  var parsed=enhancedParseInput(raw);
  if(parsed.error){
    window._currentSeq='';
    resetSequenceDependentPanels();
    setMsg('err',parsed.error);
    document.getElementById('results').hidden=true;
    if(window.SeqModuleRouter)window.SeqModuleRouter.sync();
    return;
  }
  if(parsed.warning)setMsg('info',parsed.warning);else clearMsg();
  var seq=parsed.seq;
  resetSequenceDependentPanels();
  window._currentSeq=seq;
  if (seq !== window._lastAnalyzedSeq || window._pendingFeatures) {
    mapFeatures = (window._pendingFeatures || []).slice();
    window._lastAnalyzedSeq = seq;
  }
  window._pendingFeatures = null;
  renderResults(seq);
  initNewFeatures();
  renderDogmaStep();
  if(window.SeqModuleRouter)window.SeqModuleRouter.afterAnalyze();
}

document.getElementById('btnAnalyze').addEventListener('click',analyze);
document.getElementById('btnSample').addEventListener('click',function(e){
  e.stopPropagation();
  document.getElementById('sampleMenu').classList.toggle('show');
});
document.getElementById('btnClear').addEventListener('click',function(){
  document.getElementById('seq-input').value='';
  window._currentSeq='';
  clearWorkspace();
  mapFeatures=[];
  window._lastAnalyzedSeq=null;
  resetSequenceDependentPanels();
  clearMsg();
  document.getElementById('results').hidden=true;
  renderStructure3D();
});

// 在禁用存储或隐私模式下也保持主流程可用
function safeStorageGet(key){try{return window.localStorage.getItem(key);}catch(e){return null;}}
function safeStorageSet(key,value){try{window.localStorage.setItem(key,value);}catch(e){} }

// 深色模式切换（localStorage 持久化）
(function(){
  var root=document.documentElement;
  var saved=safeStorageGet('seq-theme');
  if(saved==='dark')root.setAttribute('data-theme','dark');
  document.getElementById('btnTheme').addEventListener('click',function(){
    var isDark=root.getAttribute('data-theme')==='dark';
    if(isDark){root.removeAttribute('data-theme');safeStorageSet('seq-theme','light');}
    else{root.setAttribute('data-theme','dark');safeStorageSet('seq-theme','dark');}
    if(!document.getElementById('results').hidden && typeof echarts!=='undefined' && lastChartArgs){
      doChartsRender(lastChartArgs.seq,lastChartArgs.avgGC,lastChartArgs.aaCounts);
    }
  });
})();

// 移动端抽屉导航
(function(){
  var nav=document.getElementById('top'),btn=document.getElementById('btnMenu');
  if(!nav||!btn)return;
  function close(){nav.classList.remove('menu-open');btn.setAttribute('aria-expanded','false');btn.setAttribute('aria-label','打开导航菜单');}
  btn.addEventListener('click',function(){
    var open=nav.classList.toggle('menu-open');
    btn.setAttribute('aria-expanded',String(open));
    btn.setAttribute('aria-label',open?'关闭导航菜单':'打开导航菜单');
  });
  document.querySelectorAll('#siteLinks a').forEach(function(a){
    a.addEventListener('click',function(e){
      close();
    });
  });
  document.addEventListener('keydown',function(e){if(e.key==='Escape')close();});
})();

// 回到顶部按钮
(function(){
  var btn=document.getElementById('toTop');
  window.addEventListener('scroll',function(){
    if(window.scrollY>420)btn.classList.add('show');else btn.classList.remove('show');
  });
  btn.addEventListener('click',function(){window.scrollTo({top:0,behavior:'smooth'});});
})();

// 酶切表格筛选
document.getElementById('enzFilter').addEventListener('change',function(){
  if(lastEnz)renderEnzTable(lastSeq,lastEnz,this.value);
});
// 凝胶酶选择
document.getElementById('gelEnzSelect').addEventListener('change',function(){
  if(lastEnz)renderGel(lastSeq,lastEnz,this.value);
});

// 比对序列着色
function colorDNA(s){
  var h='';
  for(var i=0;i<s.length;i++){
    var ch=s.charAt(i);
    if(ch==='-')h+='<span class="gapch">-</span>';
    else h+='<span class="b'+ch+'">'+ch+'</span>';
  }
  return h;
}
function colorMid(s){
  var h='';
  for(var i=0;i<s.length;i++){
    var ch=s.charAt(i);
    if(ch==='|')h+='<span class="match">|</span>';
    else if(ch==='·')h+='<span class="mismatch">·</span>';
    else h+=ch;
  }
  return h;
}

// 双序列全局比对
document.getElementById('btnAlign').addEventListener('click',function(){
  var out=document.getElementById('alignResult');
  var p1=parseInput(document.getElementById('seq-input').value);
  var p2=parseInput(document.getElementById('alignInput2').value);
  if(p1.error||p2.error){out.textContent=p1.error||p2.error;return;}
  if(!p1.seq||!p2.seq){
    out.innerHTML='<p class="empty-note">需要两条序列才能比对：请确认上方主序列与第二条序列都已填写。</p>';
    return;
  }
  if(p1.seq.length>1200||p2.seq.length>1200){
    out.innerHTML='<p class="empty-note">为保证浏览器流畅，本工具比对长度上限为 1200 bp，当前序列超出。</p>';
    return;
  }
  var sc=+document.getElementById('pMatch').value;
  var mm=+document.getElementById('pMismatch').value;
  var gp=+document.getElementById('pGap').value;
  var r=nwAlign(p1.seq,p2.seq,sc,mm,gp);
  var identity=(r.match/r.len*100);
  var gaps=(r.a.match(/-/g)||[]).length;
  var h='<div class="align-stats">'+
    '<span>比对得分 <b>'+r.score+'</b></span>'+
    '<span>一致位点 <b>'+r.match+'</b></span>'+
    '<span>空位数 <b>'+gaps+'</b></span>'+
    '<span>序列相似度 <b>'+identity.toFixed(1)+'%</b></span></div>';
  h+='<div class="align-view">'+
    '<div>'+colorDNA(r.a)+'</div>'+
    '<div>'+colorMid(r.mid)+'</div>'+
    '<div>'+colorDNA(r.b)+'</div></div>';
  out.innerHTML=h;
});

/* ============================================================
   二、新增功能：数据与工具函数
   ============================================================ */

// 更多示例序列
