function buildSampleMenu(){
  var menu = document.getElementById('sampleMenu');
  if(!menu) return;
  if(menu.getAttribute('data-bound')==='1') return;
  menu.setAttribute('data-bound','1');
  menu.innerHTML = '';
  for(var i=0;i<SAMPLE_SEQUENCES.length;i++){
    (function(idx){
      var btn = document.createElement('button');
      btn.innerHTML = SAMPLE_SEQUENCES[idx].name + '<small>'+SAMPLE_SEQUENCES[idx].desc+'</small>';
      btn.onclick = function(){
        var seq = SAMPLE_SEQUENCES[idx].seq;
        if(!seq) seq = buildSampleDNA();
        document.getElementById('seq-input').value = seq;
        document.getElementById('sampleMenu').classList.remove('show');
        analyze();
      };
      menu.appendChild(btn);
    })(i);
  }
  // 点击外部关闭
  document.addEventListener('click', function(e){
    if(!e.target.closest('.sample-dropdown')){
      document.getElementById('sampleMenu').classList.remove('show');
    }
  });
}

/* ============================================================
   九、引导教程
   ============================================================ */
var tutorialStep = 0;
var tutorialReturnFocus=null;
function showTutorial(){
  if(safeStorageGet('seqworkshop_tutorial_done')) return;
  tutorialStep = 0;
  renderTutorial();
  tutorialReturnFocus=document.activeElement;
  document.getElementById('tutorialOverlay').classList.add('show');
  document.getElementById('tutorialNext').focus();
}
function renderTutorial(){
  var step = TUTORIAL_STEPS[tutorialStep];
  document.getElementById('tutorialTitle').textContent = step.title;
  document.getElementById('tutorialText').textContent = step.text;
  var prog = document.getElementById('tutorialProgress');
  prog.innerHTML = '';
  for(var i=0;i<TUTORIAL_STEPS.length;i++){
    var d = document.createElement('div');
    d.className = 'dot' + (i<=tutorialStep?' active':'');
    prog.appendChild(d);
  }
  document.getElementById('tutorialNext').textContent = tutorialStep===TUTORIAL_STEPS.length-1?'开始使用':'下一步';
}
function initTutorial(){
  document.getElementById('tutorialNext').onclick = function(){
    if(tutorialStep<TUTORIAL_STEPS.length-1){tutorialStep++;renderTutorial();}
    else{closeTutorial();}
  };
  document.getElementById('tutorialSkip').onclick = closeTutorial;
  document.getElementById('tutorialOverlay').addEventListener('keydown',function(event){
    if(event.key==='Escape'){event.preventDefault();closeTutorial();return;}
    if(event.key!=='Tab')return;
    var first=document.getElementById('tutorialSkip'),last=document.getElementById('tutorialNext');
    if(event.shiftKey && document.activeElement===first){event.preventDefault();last.focus();}
    else if(!event.shiftKey && document.activeElement===last){event.preventDefault();first.focus();}
  });
}
function closeTutorial(){
  document.getElementById('tutorialOverlay').classList.remove('show');
  safeStorageSet('seqworkshop_tutorial_done','1');
  if(tutorialReturnFocus && typeof tutorialReturnFocus.focus==='function')tutorialReturnFocus.focus();
}

/* ============================================================
   十一、增强输入校验
   ============================================================ */
function enhancedParseInput(raw){return SeqCore.parseDNA(raw);}

/* ============================================================
   十二、初始化所有新功能（在 analyze() 中调用）
   ============================================================ */
var newFeaturesInitialized=false;
function initNewFeatures(){
  if(newFeaturesInitialized)return;
  newFeaturesInitialized=true;
  initDogmaAnimation();
  initMutationSimulator();
  initPrimerDesign();
  initCodonOptimizer();
  initStructureViewer();
  buildSampleMenu();
  initTutorial();
  renderQuiz();
  addExportButtons();
  // 示例按钮下拉菜单切换已在页面加载时通过 addEventListener 绑定，此处不重复绑定
}

/* ============================================================
   页面加载初始化（不依赖序列分析的功能提前就绪）
   ============================================================ */
(function(){
  // 所有增强功能只初始化一次，避免首次分析时重复绑定事件。
  if(typeof initNewFeatures==='function') initNewFeatures();
  // 首次使用自动显示引导（延迟1秒，等页面渲染完成）
  setTimeout(function(){
    if(typeof showTutorial==='function') showTutorial();
  }, 1000);
})();

// GitHub Pages / localhost 下注册 PWA；直接双击 file:// 时跳过，避免误报错
(function(){
  if(!('serviceWorker' in navigator))return;
  if(location.protocol!=='https:' && location.hostname!=='localhost')return;
  window.addEventListener('load',function(){
    navigator.serviceWorker.register('./service-worker.js',{scope:'./'}).catch(function(err){
      console.warn('Service worker registration failed:',err);
    });
  });
})();


/* ============================================================
   十三、v3 生信增强（对标 SnapGene / Benchling / Primer3 / BLAST 的基础能力）
   ------------------------------------------------------------
   1) 遗传密码表变体（标准 / 脊椎动物线粒体）
   2) 最近邻热力学：Tm、自身二聚体 ΔG、发夹 ΔG（SantaLucia 1998）
   3) 引物设计 v2：候选枚举 + 罚分排序 + 模板唯一性
   4) 限制性酶切 v2：线性/环状、单/双酶切、黏性末端、片段表
   5) GenBank 导出（可被 SnapGene / Benchling / UGENE 直接读取）
   6) 局部比对（Smith–Waterman）+ 点阵图
   ============================================================ */

/* ---------------- 1. 遗传密码表变体 ---------------- */
