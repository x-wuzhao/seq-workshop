/* ===== 三维结构查看器：WebGL（three.js）+ 真实 PDB 坐标 ===== */
var s3d = { viewer:null, ready:false, mode:'seq', pdb:'1BNA', bare:false, win:24, style:'cartoon', colorMode:'ss', posTouched:false };
var S3D_LIB = ['1BNA','1EHZ','1L2Y','1GFL','2HHB'];
var S3D_SS_NAME = { H:'α 螺旋', E:'β 折叠', C:'无规卷曲' };

function s3dIsDark(){ return document.documentElement.getAttribute('data-theme')==='dark'; }

function s3dShowFallback(msg){
  var fb=document.getElementById('s3dFallback');
  if(fb){fb.hidden=false;fb.innerHTML=msg;}
  var em=document.getElementById('s3dEmpty'); if(em)em.hidden=true;
}

function s3dChip(text){
  var el=document.getElementById('s3dChip');
  if(el)el.textContent=text;
}

function s3dLegend(kind,chains){
  var el=document.getElementById('s3dLegend');
  if(!el)return;
  if(kind==='dna'||kind==='base'){
    el.innerHTML='<span><i style="background:#2E8B4F"></i>A</span><span><i style="background:#D64545"></i>T</span>'+
      '<span><i style="background:#D98A16"></i>G</span><span><i style="background:#2F7DC1"></i>C</span>'+
      '<span><i style="background:#0E7C7B"></i>骨架</span>';
  }else if(kind==='chain'){
    var colors=['#0E7C7B','#B9781B','#2F7DC1','#2E8B4F','#D64545'];
    var list=chains&&chains.length?chains:['A','B','C','D'];
    var html='';
    for(var i=0;i<Math.min(list.length,5);i++){
      html+='<span><i style="background:'+colors[i]+'"></i>链 '+list[i]+'</span>';
    }
    el.innerHTML=html;
  }else{
    el.innerHTML='<span><i style="background:#0E7C7B"></i>α 螺旋</span><span><i style="background:#B9781B"></i>β 折叠</span>'+
      '<span><i style="background:#7C8B89"></i>无规卷曲</span>';
  }
}

function s3dStructureInfo(id){
  var d=window.STRUCTURES&&window.STRUCTURES[id];
  var box=document.getElementById('s3dStructInfo');
  if(!d||!box)return;
  var res=d.resolution?('分辨率 '+d.resolution+' Å'):'溶液核磁（多模型）';
  box.innerHTML='<strong>'+d.title+'</strong><br>'+
    d.method+' · '+res+' · '+d.n_res+' 个残基 · '+d.n_atoms+' 个重原子<br>'+
    d.note+'<br><a href="'+d.url+'" target="_blank" rel="noopener">在 RCSB PDB 查看原始记录 ↗</a>';
}

function s3dBuildLib(){
  var box=document.getElementById('s3dLib');
  if(!box)return;
  var cur=window._currentSeq||'';
  box.innerHTML='';
  for(var i=0;i<S3D_LIB.length;i++){
    (function(id){
      var d=window.STRUCTURES&&window.STRUCTURES[id];
      if(!d)return;
      var btn=document.createElement('button');
      btn.type='button';
      btn.className='s3d-libcard'+(s3d.mode==='lib'&&s3d.pdb===id?' active':'');
      var meta=d.method+(d.resolution?(' · '+d.resolution+' Å'):'')+' · '+d.n_res+' 残基';
      var match=(id==='1BNA'&&cur.indexOf('CGCGAATTCGCG')>=0)?
        '<span class="s3d-match">与当前序列一致</span>':'';
      btn.innerHTML='<b>'+d.title+'</b><span>'+d.subtitle+'</span><em>'+meta+'</em>'+match;
      btn.addEventListener('click',function(){ openS3dStructure(id); });
      box.appendChild(btn);
    })(S3D_LIB[i]);
  }
}

function openS3dStructure(id){
  var d=window.STRUCTURES&&window.STRUCTURES[id];
  if(!d||!s3d.ready)return;
  s3d.mode='lib'; s3d.pdb=id;
  s3d.style=(d.kind==='protein')?'cartoon':(d.mode==='full'?'stick':'cartoon');
  s3d.colorMode=(d.kind==='protein'&&d.chains&&d.chains.length>1)?'chain':'ss';
  renderStructure3D();
  s3dStructureInfo(id);
  s3dBuildLib();
}

function s3dSwitchTab(name){
  var tabs=document.querySelectorAll('.s3d-tab');
  for(var i=0;i<tabs.length;i++){
    var on=tabs[i].getAttribute('data-tab')===name;
    tabs[i].classList.toggle('active',on);
    tabs[i].setAttribute('aria-selected',on?'true':'false');
  }
  var panes=document.querySelectorAll('.s3d-pane');
  for(var j=0;j<panes.length;j++)panes[j].hidden=panes[j].getAttribute('data-pane')!==name;
  if(name==='seq'){ s3d.mode='seq'; renderStructure3D(); }
  else if(!(s3d.mode==='lib'&&s3d.pdb)){ openS3dStructure(s3d.pdb||'1BNA'); }
  else { renderStructure3D(); }
}

function s3dInit(){
  if(s3d.ready)return;
  var host=document.getElementById('s3dHost');
  if(!host)return;
  if(!window.THREE||!window.Struct3D){
    s3dShowFallback('三维模块未能加载：请确认 assets/vendor/three.min.js、assets/structures/structures.js、assets/js/structure3d.js 与页面一起部署。');
    return;
  }
  var viewer=window.Struct3D.create(host,{dark:s3dIsDark()});
  if(!viewer){
    s3dShowFallback('当前浏览器未启用 WebGL，无法显示三维模型。请换用最新版 Chrome / Edge / Safari，或开启浏览器硬件加速后重试。');
    return;
  }
  s3d.viewer=viewer;
  s3d.ready=true;

  viewer.onPick=function(data,hit){
    if(!data)return;
    if(data.kind==='bp'){
      var idx=data.index+1;
      var slider=document.getElementById('structurePosition');
      if(slider)slider.value=idx;
      s3d.posTouched=true;
      var pv=document.getElementById('structurePosValue');
      if(pv)pv.textContent=idx+' / '+((window._currentSeq||'').length);
      var nH=(data.base==='G'||data.base==='C')?3:(data.base==='A'||data.base==='T'?2:null);
      var info=document.getElementById('structureInfo');
      if(info)info.innerHTML='<strong>第 '+idx+' 位：</strong>'+data.base+'–'+data.comp+' 碱基对，'+
        (nH===null?'氢键数量未知':nH+' 条氢键'+(nH===3?'（G≡C，三键更稳定）':'（A=T，两键）'))+'。';
    }else if(data.kind==='residue'){
      var r=data.residues[hit.instanceId];
      if(!r)return;
      var d=window.STRUCTURES[s3d.pdb]||{};
      var ss=d.ssmap?(d.ssmap[r[2]+':'+r[3]]||''):'';
      var box=document.getElementById('s3dStructInfo');
      if(box)box.innerHTML='<strong>'+r[0]+'</strong>（'+r[1]+'）· 链 '+r[2]+' · 残基 '+r[3]+
        (ss?(' · '+(S3D_SS_NAME[ss]||ss)):'')+'<br><span style="color:var(--ink-soft);font-size:12px">'+d.title+'</span>';
    }
  };

  var tabs=document.querySelectorAll('.s3d-tab');
  for(var i=0;i<tabs.length;i++){
    tabs[i].addEventListener('click',function(){ s3dSwitchTab(this.getAttribute('data-tab')); });
  }
  var spin=document.getElementById('s3dSpin');
  if(spin)spin.addEventListener('click',function(){
    var on=this.getAttribute('aria-pressed')!=='true';
    this.setAttribute('aria-pressed',on?'true':'false');
    this.textContent=on?'⏸ 停止旋转':'⟳ 自动旋转';
    s3d.viewer.setAutoRotate(on);
  });
  var reset=document.getElementById('s3dReset');
  if(reset)reset.addEventListener('click',function(){ s3d.viewer.resetView(); });

  var slider=document.getElementById('structurePosition');
  if(slider)slider.addEventListener('input',function(){ s3d.posTouched=true; renderStructure3D(); });
  var win=document.getElementById('s3dWin');
  if(win)win.addEventListener('input',function(){
    s3d.win=Number(this.value)||24;
    var out=document.getElementById('s3dWinValue');
    if(out)out.textContent=s3d.win+' bp';
    renderStructure3D();
  });
  var bare=document.getElementById('s3dBare');
  if(bare)bare.addEventListener('click',function(){
    s3d.bare=!s3d.bare;
    this.setAttribute('aria-pressed',s3d.bare?'true':'false');
    this.textContent=s3d.bare?'显示碱基对':'隐藏碱基对（观察大小沟）';
    renderStructure3D();
  });
  var libStyle=document.getElementById('s3dLibStyle');
  if(libStyle)libStyle.addEventListener('click',function(){
    var d=window.STRUCTURES[s3d.pdb]||{};
    if(d.mode==='ca'&&s3d.style==='cartoon')s3d.style='stick';
    else if(d.mode==='ca'&&s3d.style==='stick')s3d.style='cartoon';
    else s3d.style=(s3d.style==='stick')?'sphere':'stick';
    renderStructure3D(); s3dUpdateStyleLabel();
  });
  var libColor=document.getElementById('s3dLibColor');
  if(libColor)libColor.addEventListener('click',function(){
    s3d.colorMode=(s3d.colorMode==='chain')?'ss':'chain';
    renderStructure3D(); s3dUpdateStyleLabel();
  });

  if(document.getElementById('seqInputShow')){
    document.getElementById('seqInputShow').addEventListener('click',function(e){
      var target=e.target.closest('[data-seq-index]');
      if(!target)return;
      var slider2=document.getElementById('structurePosition');
      if(slider2){
        slider2.value=Number(target.getAttribute('data-seq-index'))+1;
        s3d.mode='seq';
        s3dSwitchTab('seq');
      }
      if(window.SeqModuleRouter)window.SeqModuleRouter.navigate('explore/structure');
    });
  }

  // 主题切换时同步三维背景
  var themeBtn=document.getElementById('btnTheme');
  if(themeBtn)themeBtn.addEventListener('click',function(){
    if(s3d.viewer)setTimeout(function(){ s3d.viewer.setTheme(s3dIsDark()); },0);
  });

  s3dBuildLib();
  renderStructure3D();
}

function s3dUpdateStyleLabel(){
  var a=document.getElementById('s3dLibStyle');
  var d=window.STRUCTURES[s3d.pdb]||{};
  if(a){
    if(d.mode==='ca')a.textContent='显示：'+(s3d.style==='cartoon'?'卡通 / 主链连线':'');
    else a.textContent='显示：'+(s3d.style==='stick'?'球棍':'球模型');
  }
  var b=document.getElementById('s3dLibColor');
  if(b)b.textContent='着色：'+(s3d.colorMode==='chain'?'按链':'按二级结构');
}

function renderStructure3D(){
  if(document.getElementById('nav-structure').hidden)return;
  if(!s3d.ready)s3dInit();
  if(!s3d.ready)return;
  var seq=window._currentSeq||'';
  var empty=document.getElementById('s3dEmpty');
  if(s3d.mode==='seq'&&!seq){
    s3d.viewer.clear();
    if(empty)empty.hidden=false;
    var sl=document.getElementById('structurePosition');
    if(sl){sl.disabled=true;sl.max=1;sl.value=1;}
    var pv=document.getElementById('structurePosValue'); if(pv)pv.textContent='—';
    var info=document.getElementById('structureInfo'); if(info)info.textContent='请先输入 DNA 序列并开始分析。';
    s3dChip('B-DNA 双螺旋');
    return;
  }
  if(empty)empty.hidden=true;

  if(s3d.mode==='seq'){
    var slider=document.getElementById('structurePosition');
    var mid=Math.ceil(seq.length/2);
    var center=s3d.posTouched?Math.min(seq.length,Math.max(1,Number(slider&&slider.value)||mid)):mid;
    if(slider){slider.disabled=false;slider.max=seq.length;slider.value=center;}
    var posValue=document.getElementById('structurePosValue');
    if(posValue)posValue.textContent=center+' / '+seq.length;
    s3d.viewer.showDNA(seq,center,s3d.win,s3d.bare);
    s3dChip('B-DNA 双螺旋 · 由你的序列生成');
    s3dLegend('dna');
    var info2=document.getElementById('structureInfo');
    if(info2){
      var b=seq.charAt(center-1),pair=reverseComplement(b);
      var nH=(b==='G'||b==='C')?3:(b==='A'||b==='T'?2:null);
      info2.innerHTML='<strong>第 '+center+' 位：</strong>'+b+'–'+pair+' 碱基对，'+(nH===null?'氢键数量未知':nH+' 条氢键')+'。'+
        '每个碱基对上升 3.38 Å、旋转 36°，约 10.5 bp 转一圈。';
    }
  }else{
    var d=window.STRUCTURES&&window.STRUCTURES[s3d.pdb];
    if(!d)return;
    s3d.viewer.showStructure(s3d.pdb,s3d.style,s3d.colorMode);
    s3dChip(d.title);
    s3dLegend(d.kind==='protein'?(s3d.colorMode==='chain'?'chain':'ss'):'dna', d.chains);
    s3dStructureInfo(s3d.pdb);
    s3dUpdateStyleLabel();
  }
}

// 兼容旧调用点
function initStructureViewer(){ if(!document.getElementById('nav-structure').hidden)s3dInit(); }

// 酶切线性图谱（SVG，悬停显示酶名与位置）
