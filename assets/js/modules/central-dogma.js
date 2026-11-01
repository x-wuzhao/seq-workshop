var dogmaStep = 0;
var dogmaAutoTimer = null;
var dogmaSteps = [
  {title:'第一步：DNA 双链', desc:'DNA 是由两条反向平行的多核苷酸链组成的双螺旋结构。两条链通过碱基互补配对（A-T、G-C）连接，A 和 T 之间有2个氢键，G 和 C 之间有3个氢键。', visualFn:function(seq){
    var rc = reverseComplement(seq);
    var top = seqHTML(seq.substring(0,40));
    var bot = seqHTML(complement(seq.substring(0,40)));
    return '<div style="text-align:center"><div style="font-size:11px;color:var(--ink-soft);margin-bottom:4px">编码链 5\'→3\'</div>'+top+'<div style="font-size:11px;color:var(--ink-soft);margin:8px 0 4px">模板链 3\'→5\'</div>'+bot+'</div>';
  }},
  {title:'第二步：DNA 复制', desc:'DNA 复制是半保留复制：双链解开后，每条母链作为模板，按照碱基互补原则合成一条新的子链。最终一个 DNA 分子变成两个完全相同的 DNA 分子，每个分子都含一条母链和一条新合成的子链。', visualFn:function(seq){
    var rc = reverseComplement(seq);
    return '<div style="text-align:center"><div style="color:var(--ink-soft);font-size:11px;margin-bottom:4px">亲代 DNA（母链）</div>'+seqHTML(seq.substring(0,30))+'<div style="margin:8px 0;color:var(--teal);font-size:12px;font-weight:700">↓ 半保留复制 ↓</div><div style="display:flex;justify-content:space-around;gap:20px"><div><div style="font-size:11px;color:var(--ink-soft);margin-bottom:4px">子代 DNA 1（母链+新链）</div>'+seqHTML(seq.substring(0,30))+'</div><div><div style="font-size:11px;color:var(--ink-soft);margin-bottom:4px">子代 DNA 2（新链+母链）</div>'+seqHTML(seq.substring(0,30))+'</div></div></div>';
  }},
  {title:'第三步：DNA 解旋', desc:'在转录开始前，DNA 双链在 RNA 聚合酶的作用下局部解开，形成转录泡。模板链（反义链）将作为合成 mRNA 的模板，编码链（正义链）不参与转录。', visualFn:function(seq){
    var rc = reverseComplement(seq);
    return '<div style="text-align:center"><div style="color:var(--ink-soft);font-size:11px;margin-bottom:4px">编码链（不参与转录）</div>'+seqHTML(seq.substring(0,40))+'<div style="margin:6px 0;color:var(--red);font-size:12px;font-weight:700">↓ RNA 聚合酶结合，双链解开 ↓</div><div style="color:var(--ink-soft);font-size:11px;margin-bottom:4px">模板链（转录模板）</div>'+seqHTML(complement(seq.substring(0,40)))+'</div>';
  }},
  {title:'第四步：转录成 mRNA', desc:'RNA 聚合酶沿着模板链 3\'→5\' 方向移动，按照碱基互补原则合成 mRNA。注意：mRNA 中用尿嘧啶（U）替代了胸腺嘧啶（T）。最终 mRNA 的序列与 DNA 编码链几乎相同，只是 T→U。', visualFn:function(seq){
    var mrna = transcribe(seq);
    return '<div style="text-align:center"><div style="color:var(--ink-soft);font-size:11px;margin-bottom:4px">DNA 编码链</div>'+seqHTML(seq.substring(0,40))+'<div style="margin:6px 0;color:var(--teal);font-size:12px;font-weight:700">↓ 转录（T→U）↓</div><div style="color:var(--ink-soft);font-size:11px;margin-bottom:4px">mRNA</div>'+seqHTML(mrna.substring(0,40),true)+'</div>';
  }},
  {title:'第五步：mRNA 出核', desc:'在真核生物中，转录完成后，mRNA 经过加工（5\'加帽、3\'加 polyA 尾、剪接去除内含子），从细胞核进入细胞质，与核糖体结合，准备开始翻译。', visualFn:function(seq){
    var mrna = transcribe(seq);
    return '<div style="text-align:center;padding:16px 0"><div style="display:inline-block;background:var(--teal-soft);border:2px solid var(--teal);border-radius:12px;padding:20px 30px;margin-bottom:14px"><div style="display:flex;align-items:center;gap:12px"><div style="width:40px;height:40px;border-radius:50%;background:var(--teal);display:flex;align-items:center;justify-content:center;color:#fff;font-weight:700;font-size:18px">N</div><div style="color:var(--ink-soft);font-size:20px">→</div><div style="width:40px;height:40px;border-radius:50%;background:var(--amber);display:flex;align-items:center;justify-content:center;color:#fff;font-weight:700;font-size:18px">C</div></div><div style="margin-top:8px;font-size:12px;color:var(--teal-dark);font-weight:600">细胞核 → 细胞质 → 核糖体</div></div><div style="font-family:JetBrains Mono,monospace;font-size:12px;color:var(--ink-soft)">5\'-m7G cap  '+mrna.substring(0,20)+'...  AAAAAAAA-3\'</div></div>';
  }},
  {title:'第六步：翻译起始', desc:'核糖体识别 mRNA 上的起始密码子 AUG（编码甲硫氨酸 M），翻译开始。tRNA 携带对应的氨基酸进入核糖体，通过反密码子与 mRNA 上的密码子配对。', visualFn:function(seq){
    var mrna = transcribe(seq);
    var codons = [];
    for(var i=0;i<Math.min(9,Math.floor(mrna.length/3));i++){
      var c = mrna.substring(i*3,i*3+3);
      var aa = CODON[c.replace(/U/g,'T')] || '?';
      codons.push('<span style="display:inline-block;margin:2px;padding:3px 8px;background:var(--card);border:1px solid var(--line);border-radius:5px"><span class="mono" style="color:var(--bU)">'+c+'</span> → <strong>'+aa+'</strong></span>');
    }
    return '<div style="text-align:center"><div style="color:var(--amber);font-weight:700;margin-bottom:10px">起始密码子 AUG → 甲硫氨酸（M）</div>'+codons.join('')+'</div>';
  }},
  {title:'第七步：肽链延伸', desc:'核糖体沿着 mRNA 5\'→3\' 方向移动，每读取一个密码子，对应的 tRNA 携带氨基酸进入，氨基酸之间通过肽键连接，肽链不断延伸。直到遇到终止密码子（UAA、UAG、UGA），翻译结束。', visualFn:function(seq){
    var mrna = transcribe(seq);
    var prot = '';
    for(var i=0;i+3<=mrna.length;i+=3){
      var aa = CODON[mrna.substring(i,i+3).replace(/U/g,'T')];
      if(aa==='*'){prot += '<span style="color:var(--red);font-weight:700">*</span>';break;}
      prot += aa || '?';
    }
    return '<div style="text-align:center"><div style="color:var(--ink-soft);font-size:11px;margin-bottom:6px">mRNA（每3个碱基=1个密码子）</div><div class="mono" style="font-size:12px;letter-spacing:1px">'+mrna.substring(0,60)+'...</div><div style="margin:12px 0;color:var(--teal);font-weight:700">↓ 翻译 ↓</div><div style="color:var(--ink-soft);font-size:11px;margin-bottom:6px">多肽链（氨基酸序列）</div><div class="mono" style="font-size:14px;letter-spacing:2px;color:var(--teal-dark);font-weight:600">'+prot.substring(0,40)+'...</div></div>';
  }},
  {title:'中心法则完成！', desc:'遗传信息从 DNA → DNA（复制）→ mRNA（转录）→ 蛋白质（翻译）的完整流动过程结束了。这就是分子生物学的中心法则。序列工坊让你可以亲手拨动这个过程的每一步——从粘贴一段 DNA 开始，亲眼看着它变成蛋白质。', visualFn:function(seq){
    return '<div style="text-align:center;padding:24px 0"><div style="display:inline-flex;align-items:center;gap:10px;margin-bottom:18px"><div style="width:56px;height:56px;border-radius:12px;background:linear-gradient(135deg,var(--teal),var(--teal-dark));display:flex;align-items:center;justify-content:center;color:#fff;font-weight:900;font-size:22px">DNA</div><div style="color:var(--teal);font-size:24px;font-weight:700">→</div><div style="width:56px;height:56px;border-radius:12px;background:linear-gradient(135deg,var(--amber),#d4941a);display:flex;align-items:center;justify-content:center;color:#fff;font-weight:900;font-size:18px">RNA</div><div style="color:var(--teal);font-size:24px;font-weight:700">→</div><div style="width:56px;height:56px;border-radius:12px;background:linear-gradient(135deg,var(--red),#a02828);display:flex;align-items:center;justify-content:center;color:#fff;font-weight:900;font-size:16px">Pro</div></div><div style="font-size:18px;font-weight:700;color:var(--teal-dark);margin-bottom:6px">DNA → mRNA → 蛋白质</div><div style="font-size:14px;color:var(--ink-soft)">中心法则，看得见、摸得着</div></div>';
  }}
];

function renderDogmaDots(){
  var box = document.getElementById('dogmaDots');
  if(!box) return;
  box.innerHTML = '';
  for(var i=0;i<dogmaSteps.length;i++){
    (function(idx){
      var d = document.createElement('button');
      d.type='button';
      d.className = 'dogma-dot' + (idx===dogmaStep?' active':'') + (idx<dogmaStep?' done':'');
      d.textContent = idx+1;
      d.title = dogmaSteps[idx].title;
      d.setAttribute('aria-label',dogmaSteps[idx].title);
      d.setAttribute('aria-current',idx===dogmaStep?'step':'false');
      d.onclick = function(){
        stopDogmaAuto();
        dogmaStep = idx;
        renderDogmaStep();
      };
      box.appendChild(d);
    })(i);
  }
}

function renderDogmaStep(){
  var step = dogmaSteps[dogmaStep];
  var seq = (window.SeqModuleRouter && !window.SeqModuleRouter.hasSequence()) ? SAMPLE_DNA : (window._currentSeq || SAMPLE_DNA);
  var content = document.getElementById('dogmaContent');
  var visual = document.getElementById('dogmaVisual');
  content.classList.add('dogma-fade');
  visual.classList.add('dogma-fade');
  setTimeout(function(){
    content.innerHTML = '<strong style="font-size:16px;color:var(--teal-dark)">'+step.title+'</strong><p style="margin-top:8px">'+step.desc+'</p>';
    visual.innerHTML = step.visualFn(seq);
    content.classList.remove('dogma-fade');
    visual.classList.remove('dogma-fade');
  }, 200);
  document.getElementById('dogmaCounter').textContent = '第 '+(dogmaStep+1)+' / '+dogmaSteps.length+' 步';
  var timeline=document.getElementById('dogmaTimeline');
  var timelineValue=document.getElementById('dogmaTimelineValue');
  if(timeline)timeline.value=dogmaStep;
  if(timelineValue)timelineValue.textContent=(dogmaStep+1)+' / '+dogmaSteps.length;
  document.getElementById('dogmaPrev').disabled = dogmaStep===0;
  document.getElementById('dogmaNext').textContent = dogmaStep===dogmaSteps.length-1?'完成':'下一步';
  renderDogmaDots();
}

function stopDogmaAuto(){
  if(dogmaAutoTimer){
    clearInterval(dogmaAutoTimer);
    dogmaAutoTimer = null;
    var btn = document.getElementById('dogmaAuto');
    if(btn){ btn.textContent = '▶ 自动播放'; btn.classList.remove('dogma-auto-playing'); }
  }
}

function startDogmaAuto(){
  if(dogmaStep >= dogmaSteps.length-1){ dogmaStep = 0; renderDogmaStep(); }
  var btn = document.getElementById('dogmaAuto');
  if(btn){ btn.textContent = '⏸ 暂停'; btn.classList.add('dogma-auto-playing'); }
  dogmaAutoTimer = setInterval(function(){
    if(dogmaStep < dogmaSteps.length-1){
      dogmaStep++;
      renderDogmaStep();
    } else {
      stopDogmaAuto();
    }
  }, 2200);
}

function initDogmaAnimation(){
  document.getElementById('dogmaNext').onclick = function(){
    stopDogmaAuto();
    if(dogmaStep<dogmaSteps.length-1){dogmaStep++;renderDogmaStep();}
  };
  document.getElementById('dogmaPrev').onclick = function(){
    stopDogmaAuto();
    if(dogmaStep>0){dogmaStep--;renderDogmaStep();}
  };
  document.getElementById('dogmaReset').onclick = function(){stopDogmaAuto();dogmaStep=0;renderDogmaStep();};
  document.getElementById('dogmaAuto').onclick = function(){
    if(dogmaAutoTimer){ stopDogmaAuto(); } else { startDogmaAuto(); }
  };
  document.getElementById('dogmaTimeline').addEventListener('input',function(){
    stopDogmaAuto();
    dogmaStep=Number(this.value)||0;
    renderDogmaStep();
  });
  document.addEventListener('keydown', function(e){
    var panel = document.getElementById('nav-dogma');
    if(!panel || panel.hidden || /INPUT|TEXTAREA|SELECT/.test(e.target.tagName) || e.target.isContentEditable) return;
    var rect = panel.getBoundingClientRect();
    var visible = rect.top < window.innerHeight && rect.bottom > 0;
    if(!visible) return;
    if(e.key === 'ArrowRight' || e.key === 'ArrowDown'){
      e.preventDefault();
      stopDogmaAuto();
      if(dogmaStep<dogmaSteps.length-1){dogmaStep++;renderDogmaStep();}
    } else if(e.key === 'ArrowLeft' || e.key === 'ArrowUp'){
      e.preventDefault();
      stopDogmaAuto();
      if(dogmaStep>0){dogmaStep--;renderDogmaStep();}
    } else if(e.key === ' '){
      e.preventDefault();
      if(dogmaAutoTimer){ stopDogmaAuto(); } else { startDogmaAuto(); }
    }
  });
  renderDogmaStep();
}

/* ============================================================
   四、突变模拟器
   ============================================================ */
